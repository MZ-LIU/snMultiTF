#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Identify candidate driver TFs from the baseline GRN using:
1. Edge-weighted GRN:
   edge_weight = percentile_rank(log1p(n_motif_positive_peaks))
                 * abs(expression_assoc), with zero motif score for n_peaks <= 0

2. TF ranking metrics:
   - NetworkSupport = sum(abs(edge_weight))
   - ActDyn from R mgcv GAM:
       rank_scaled(sd(fitted activity over a uniform pseudotime grid))
       * rank_scaled(deviance explained by pseudotime)
   - ProgramCouple from R mgcv GAM:
       for TFs with at least 10 target genes, keep the top-K target genes by
       edge weight; normalize the complete RNA count matrix to 1e4 counts per
       cell followed by log1p; z-score each selected target gene across cells;
       fit Gaussian GAMs for TF activity and target program score against
       pseudotime, predict both models on the same 100-point grid, then calculate:
       abs(cor(predicted TF activity, predicted target program score))

3. Rank integration:
   RobustRankAggreg::aggregateRanks from R.

Expected companion scripts in the repo-level tf_driver folder:
    run_mgcv_tf_driver.R
    run_rra_tf_driver.R

Run:
    python workflow_labeled/7.2_TF recognition.py
"""

from __future__ import annotations

import json
import re
import sys
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import sparse


# ---------------------------------------------------------------------------
# User configuration
# ---------------------------------------------------------------------------
# R helper scripts are shared from the repo-level tf_driver folder. The input
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_DIR = SCRIPT_DIR.parent
WORKFLOW_DIR = REPO_DIR / 'workflow_unlabeled_quality control'

NETWORK_DIR = WORKFLOW_DIR / "results_sc" / "sc_multiome" / "network_analysis"
FULL_TRAJECTORY_INPUT_DIR = NETWORK_DIR / "full_trajectory_inputs"
GRN_TABLE_PATH: Path | None = None



OUTPUT_DIR = NETWORK_DIR / "tf_driver_ranking"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

PSEUDOTIME_KEY = "pseudotime"
RNA_MODALITY = "RNA"
TF_ACTIVITY_MODALITY = "chromvar_activity"
EXPRESSION_LAYER = None  # None means use RNA.X.
RNA_COUNTS_LAYER = "counts"
RNA_NORMALIZE_TARGET_SUM = 1e4

# ProgramCouple / GAM settings
MIN_TARGETS_FOR_PROGRAM = 10
TOP_K_TARGETS_FOR_PROGRAM = 100
SPLINE_K = 5
EPSILON = 1e-8
GAM_GRID_SIZE = 100

# Driver class thresholds
HIGH_CONFIDENCE_FDR = 0.10
TOP_N_DRIVER = 10
TOP_PERCENTILE_SUPPORT = 0.25
MIN_PROGRAM_COUPLE = 0.10

# R integration
R_SCRIPT_BIN = "Rscript"
R_HELPER_DIR = REPO_DIR / "tf_driver"
R_MGCV_SCRIPT = R_HELPER_DIR / "run_mgcv_tf_driver.R"
R_RRA_SCRIPT = R_HELPER_DIR / "run_rra_tf_driver.R"


# ---------------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------------
def find_full_trajectory_input_dir(input_dir: Path) -> Path:
    manifest_path = input_dir / "switchtfi_inputs_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(
            "Missing full-trajectory input manifest: "
            f"{manifest_path}. Please generate the full_trajectory_inputs export first."
        )
    return input_dir


def read_manifest(input_dir: Path) -> dict:
    manifest_path = input_dir / "switchtfi_inputs_manifest.json"
    if not manifest_path.exists():
        return {}
    return json.loads(manifest_path.read_text(encoding="utf-8"))


def resolve_path(raw: str | None, base_dir: Path, fallback: Path) -> Path:
    if raw:
        p = Path(raw)
        if p.exists():
            return p
        if not p.is_absolute() and (base_dir / p).exists():
            return base_dir / p
    return fallback


def zscore(values: Iterable[float]) -> np.ndarray:
    arr = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)
    mu = np.nanmean(arr)
    sd = np.nanstd(arr)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(arr, dtype=float)
    out = (arr - mu) / sd
    out[~np.isfinite(out)] = 0.0
    return out


def bh_fdr(pvalues: Iterable[float]) -> np.ndarray:
    p = pd.to_numeric(pd.Series(pvalues), errors="coerce").to_numpy(dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return out

    pv = np.clip(p[valid], 0.0, 1.0)
    order = np.argsort(pv)
    ranked = pv[order]
    n = len(ranked)

    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0.0, 1.0)

    tmp = np.empty_like(adj)
    tmp[order] = adj
    out[valid] = tmp
    return out


def normalize_tf_name(name: str) -> str:
    return re.sub(r"\(var\.\d+\)$", "", str(name)).upper()


def as_dense(x) -> np.ndarray:
    if sparse.issparse(x):
        return x.toarray()
    return np.asarray(x)


def matrix_subset_to_frame(adata, genes: list[str], layer: str | None = None) -> pd.DataFrame:
    genes = [g for g in genes if g in adata.var_names]
    if not genes:
        return pd.DataFrame(index=adata.obs_names.astype(str))

    x = adata[:, genes].layers[layer] if layer else adata[:, genes].X
    return pd.DataFrame(as_dense(x), index=adata.obs_names.astype(str), columns=genes)


def check_r_scripts() -> None:
    missing = []
    if not R_MGCV_SCRIPT.exists():
        missing.append(str(R_MGCV_SCRIPT))
    if not R_RRA_SCRIPT.exists():
        missing.append(str(R_RRA_SCRIPT))

    if missing:
        raise FileNotFoundError(
            "Missing required R helper script(s):\n"
            + "\n".join(missing)
            + f"\nExpected R helper directory: {R_HELPER_DIR}"
        )


def cleanup_stale_intermediate_files(output_dir: Path) -> None:
    """Remove intermediate files left by older versions of this workflow."""
    stale_names = [
        "tmp_actdyn_mgcv_input.csv",
        "tmp_program_couple_mgcv_input.csv",
        "rank_network_support_for_rra.csv",
        "rank_actdyn_for_rra.csv",
        "rank_program_for_rra.csv",
        "rra_universe.csv",
        "rra_from_RobustRankAggreg.csv",
    ]
    for name in stale_names:
        path = output_dir / name
        if path.exists():
            path.unlink()


# ---------------------------------------------------------------------------
# Baseline weighted GRN
# ---------------------------------------------------------------------------
def load_baseline_grn(input_dir: Path, manifest: dict) -> pd.DataFrame:
    paths = manifest.get("paths", {})
    if GRN_TABLE_PATH is not None:
        grn_path = Path(GRN_TABLE_PATH)
    else:
        grn_path = resolve_path(
            paths.get("baseline_grn_full_csv"),
            input_dir,
            input_dir / "grn_full_baseline.csv",
        )

    if not grn_path.exists():
        raise FileNotFoundError(
            "Missing baseline GRN for the selected SwitchTFI input directory: "
            f"{grn_path}. Please provide GRN_TABLE_PATH explicitly if you intend "
            "to use a different GRN table."
        )

    grn = pd.read_csv(grn_path)
    grn = grn.rename(columns={"TF": "tf", "target": "gene", "Gene": "gene"})

    required = {"tf", "gene", "correlation"}
    missing = required.difference(grn.columns)
    if missing:
        raise ValueError(f"Baseline GRN missing required columns: {sorted(missing)}")

    if "n_peaks" not in grn.columns:
        grn["n_peaks"] = np.nan

    grn["tf"] = grn["tf"].astype(str).map(normalize_tf_name)
    grn["gene"] = grn["gene"].astype(str)
    return grn


def add_edge_weights(grn: pd.DataFrame) -> pd.DataFrame:
    """
    Build non-negative edge weights from motif-positive peak support and
    expression association strength.

    motif_peak_count = n_peaks
    motif_score = 0 for n_peaks <= 0
    motif_score = percentile_rank(log1p(motif_peak_count)) for n_peaks > 0
    expr_score = abs(expression_assoc)
    edge_weight = motif_score * expr_score
    """
    out = grn.copy()

    out["motif_peak_count"] = pd.to_numeric(out["n_peaks"], errors="coerce")
    if out["motif_peak_count"].isna().all():
        raise ValueError(
            "Cannot compute motif_peak_count: n_peaks is missing. "
            "The current 5.2 baseline should contain this column."
        )

    out["motif_peak_count"] = out["motif_peak_count"].fillna(0.0)
    out["expression_assoc"] = pd.to_numeric(out["correlation"], errors="coerce").fillna(0.0)

    # Motif/peak support: keep unsupported edges at zero, rank supported edges continuously.
    out["motif_score_log1p"] = np.log1p(out["motif_peak_count"])
    out["motif_score"] = 0.0
    positive_motif = out["motif_peak_count"] > 0
    if positive_motif.any():
        out.loc[positive_motif, "motif_score"] = (
            out.loc[positive_motif, "motif_score_log1p"]
            .rank(method="average", pct=True)
            .astype(float)
        )

    # Expression association strength; sign is not used for driver ranking edge strength.
    out["expr_score"] = np.abs(out["expression_assoc"])

    out["edge_weight"] = out["motif_score"] * out["expr_score"]
    out["edge_weight_abs"] = np.abs(out["edge_weight"])

    out["TF"] = out["tf"]
    out["target"] = out["gene"]

    zero_fraction = float((out["edge_weight"] == 0).mean())
    out.attrs["edge_weight_zero_fraction"] = zero_fraction
    return out


def rank_network_support(weighted_grn: pd.DataFrame) -> pd.DataFrame:
    rank = (
        weighted_grn.groupby("tf", as_index=False)
        .agg(
            NetworkSupport=("edge_weight_abs", "sum"),
            weighted_outdegree=("edge_weight", "sum"),
            weighted_outdegree_abs=("edge_weight_abs", "sum"),
            n_targets=("gene", "nunique"),
            n_edges=("gene", "size"),
            mean_edge_weight=("edge_weight", "mean"),
            mean_abs_edge_weight=("edge_weight_abs", "mean"),
        )
        .sort_values(["NetworkSupport", "weighted_outdegree_abs"], ascending=False)
        .reset_index(drop=True)
    )

    rank.insert(0, "NetworkSupport_rank", np.arange(1, len(rank) + 1))
    rank = rank.rename(columns={"tf": "TF"})
    return rank


# ---------------------------------------------------------------------------
# Omics loading
# ---------------------------------------------------------------------------
def load_multiome(input_dir: Path, manifest: dict):
    try:
        import muon as mu
    except Exception as exc:
        raise ImportError(
            "Reading trajectory_multiome.h5mu requires the muon package in the "
            "same Python environment used for 5.2_network.py."
        ) from exc

    paths = manifest.get("paths", {})
    h5mu_path = resolve_path(
        paths.get("trajectory_multiome_h5mu"),
        input_dir,
        input_dir / "trajectory_multiome.h5mu",
    )
    if not h5mu_path.exists():
        raise FileNotFoundError(f"Missing trajectory h5mu: {h5mu_path}")

    mdata = mu.read_h5mu(h5mu_path)
    if RNA_MODALITY not in mdata.mod:
        raise ValueError(f"h5mu missing modality '{RNA_MODALITY}'. Available: {list(mdata.mod)}")
    if TF_ACTIVITY_MODALITY not in mdata.mod:
        raise ValueError(
            f"h5mu missing modality '{TF_ACTIVITY_MODALITY}'. Available: {list(mdata.mod)}"
        )

    rna = mdata[RNA_MODALITY].copy()
    activity = mdata[TF_ACTIVITY_MODALITY].copy()

    common = rna.obs_names.astype(str).intersection(activity.obs_names.astype(str))
    if len(common) == 0:
        raise ValueError("RNA and TF activity modalities have no overlapping cells.")

    rna = rna[common, :].copy()
    activity = activity[common, :].copy()

    if PSEUDOTIME_KEY not in rna.obs.columns:
        if PSEUDOTIME_KEY in mdata.obs.columns:
            rna.obs[PSEUDOTIME_KEY] = mdata.obs.loc[rna.obs_names, PSEUDOTIME_KEY].values
        else:
            raise ValueError(f"Missing pseudotime column '{PSEUDOTIME_KEY}' in h5mu.")

    pseudotime = pd.to_numeric(rna.obs[PSEUDOTIME_KEY], errors="coerce")
    keep = pseudotime.notna().to_numpy()

    rna = rna[keep, :].copy()
    activity = activity[keep, :].copy()
    pseudotime = pd.to_numeric(rna.obs[PSEUDOTIME_KEY], errors="coerce").to_numpy(dtype=float)

    act_x = as_dense(activity.X)
    act_df = pd.DataFrame(
        act_x,
        index=activity.obs_names.astype(str),
        columns=[normalize_tf_name(v) for v in activity.var_names.astype(str)],
    )

    # Collapse duplicated TF names after normalization.
    act_df = act_df.T.groupby(level=0).mean().T

    return rna, act_df, pseudotime


def prepare_rna_for_target_program(rna):
    """Create log-normalized RNA.X from full-gene raw counts.

    Cell-level library-size normalization must be performed before selecting a
    TF-specific target-gene subset. Gene-wise z-scoring remains downstream in
    ``compute_program_couple_r_mgcv`` so that each selected target contributes
    on a comparable scale to the unweighted target-program mean.
    """
    try:
        import scanpy as sc
    except Exception as exc:
        raise ImportError(
            "Preparing target-program expression requires the scanpy package."
        ) from exc

    if RNA_COUNTS_LAYER not in rna.layers:
        raise ValueError(
            f"RNA modality is missing raw-count layer '{RNA_COUNTS_LAYER}'. "
            "ProgramCouple requires full-gene raw counts for cell-level "
            "library-size normalization."
        )

    out = rna.copy()
    out.X = out.layers[RNA_COUNTS_LAYER].copy()

    # The exported trajectory h5mu may retain stale log1p metadata even though
    # RNA.X was explicitly replaced by raw counts. Remove it before rebuilding
    # the normalized matrix to keep AnnData metadata consistent with RNA.X.
    out.uns.pop("log1p", None)

    sc.pp.normalize_total(out, target_sum=RNA_NORMALIZE_TARGET_SUM)
    sc.pp.log1p(out)
    out.uns["tf_driver_program_expression_preprocessing"] = {
        "source": f"RNA.layers['{RNA_COUNTS_LAYER}']",
        "normalize_total_target_sum": float(RNA_NORMALIZE_TARGET_SUM),
        "log1p": True,
        "target_gene_zscore": True,
        "target_program_aggregation": "unweighted_mean",
    }
    return out


# ---------------------------------------------------------------------------
# R mgcv wrappers
# ---------------------------------------------------------------------------
def compute_actdyn_r_mgcv(
    activity_df: pd.DataFrame,
    pseudotime: np.ndarray,
    tfs: list[str],
    work_dir: Path,
) -> pd.DataFrame:
    rows = []

    for tf in tfs:
        if tf not in activity_df.columns:
            continue

        y = pd.to_numeric(activity_df[tf], errors="coerce").to_numpy(dtype=float)

        rows.append(
            pd.DataFrame(
                {
                    "TF": tf,
                    "cell": activity_df.index.astype(str),
                    "pseudotime": pseudotime,
                    "activity": y,
                }
            )
        )

    if not rows:
        return pd.DataFrame()

    actdyn_input = pd.concat(rows, ignore_index=True)

    input_csv = work_dir / "tmp_actdyn_mgcv_input.csv"
    output_csv = work_dir / "rank_actdyn_mgcv.csv"
    actdyn_input.to_csv(input_csv, index=False)

    subprocess.run(
        [
            R_SCRIPT_BIN,
            str(R_MGCV_SCRIPT),
            "actdyn",
            str(input_csv),
            str(output_csv),
            str(SPLINE_K),
            str(EPSILON),
            str(GAM_GRID_SIZE),
        ],
        check=True,
    )

    out = pd.read_csv(output_csv)
    out["TF"] = out["TF"].astype(str).map(normalize_tf_name)
    return out


def compute_program_couple_r_mgcv(
    rna,
    activity_df: pd.DataFrame,
    pseudotime: np.ndarray,
    activity_aware_grn: pd.DataFrame,
    tfs: list[str],
    work_dir: Path,
) -> pd.DataFrame:
    rows = []
    gene_set = set(map(str, rna.var_names))

    for tf in tfs:
        if tf not in activity_df.columns:
            continue

        tf_edges = activity_aware_grn.loc[activity_aware_grn["tf"] == tf].copy()
        if tf_edges.empty:
            targets = []
        else:
            tf_edges["gene"] = tf_edges["gene"].astype(str)
            tf_edges = tf_edges[tf_edges["gene"].isin(gene_set)].copy()

            sort_cols = [
                c
                for c in ["edge_weight_abs", "expr_score", "motif_peak_count", "gene"]
                if c in tf_edges.columns
            ]
            ascending = [False] * (len(sort_cols) - 1) + [True] if sort_cols else True
            if sort_cols:
                tf_edges = tf_edges.sort_values(sort_cols, ascending=ascending)

            targets = (
                tf_edges.drop_duplicates("gene")
                .head(TOP_K_TARGETS_FOR_PROGRAM)["gene"]
                .tolist()
            )
        n_targets = len(targets)

        if n_targets < MIN_TARGETS_FOR_PROGRAM:
            # Include one placeholder row for this TF so R can return too_few_targets.
            rows.append(
                pd.DataFrame(
                    {
                        "TF": [tf],
                        "cell": ["placeholder"],
                        "pseudotime": [np.nan],
                        "activity": [np.nan],
                        "program_score": [np.nan],
                        "n_targets": [n_targets],
                    }
                )
            )
            continue

        expr = matrix_subset_to_frame(rna, targets, EXPRESSION_LAYER)
        expr_z = expr.apply(zscore, axis=0)
        program_score = expr_z.mean(axis=1).to_numpy(dtype=float)

        rows.append(
            pd.DataFrame(
                {
                    "TF": tf,
                    "cell": expr.index.astype(str),
                    "pseudotime": pseudotime,
                    "activity": activity_df.loc[expr.index, tf].to_numpy(dtype=float),
                    "program_score": program_score,
                    "n_targets": n_targets,
                }
            )
        )

    if not rows:
        return pd.DataFrame()

    program_input = pd.concat(rows, ignore_index=True)

    input_csv = work_dir / "tmp_program_couple_mgcv_input.csv"
    output_csv = work_dir / "rank_program_couple_mgcv.csv"
    program_input.to_csv(input_csv, index=False)

    subprocess.run(
        [
            R_SCRIPT_BIN,
            str(R_MGCV_SCRIPT),
            "program",
            str(input_csv),
            str(output_csv),
            str(SPLINE_K),
            str(EPSILON),
            str(GAM_GRID_SIZE),
        ],
        check=True,
    )

    out = pd.read_csv(output_csv)
    out["TF"] = out["TF"].astype(str).map(normalize_tf_name)
    out["target_program_weighted"] = False
    out["program_couple_shared_k"] = True
    out["program_couple_grid_size"] = GAM_GRID_SIZE
    out["top_k_targets_for_program"] = TOP_K_TARGETS_FOR_PROGRAM
    out["target_program_expression_source"] = f"RNA.layers['{RNA_COUNTS_LAYER}']"
    out["target_program_expression_preprocessing"] = (
        f"normalize_total(target_sum={RNA_NORMALIZE_TARGET_SUM:g}) -> log1p -> "
        "target_gene_zscore -> unweighted_mean"
    )
    return out


# ---------------------------------------------------------------------------
# R RobustRankAggreg wrapper
# ---------------------------------------------------------------------------
def run_rra_r_package(
    weighted_rank: pd.DataFrame,
    actdyn_rank: pd.DataFrame,
    program_rank: pd.DataFrame,
    universe: list[str],
    work_dir: Path,
) -> pd.DataFrame:
    network_csv = work_dir / "rank_network_support_for_rra.csv"
    actdyn_csv = work_dir / "rank_actdyn_for_rra.csv"
    program_csv = work_dir / "rank_program_for_rra.csv"
    universe_csv = work_dir / "rra_universe.csv"
    output_csv = work_dir / "rra_from_RobustRankAggreg.csv"

    weighted_rank[["NetworkSupport_rank", "TF"]].rename(
        columns={"NetworkSupport_rank": "rank"}
    ).to_csv(network_csv, index=False)

    actdyn_rank[["ActDyn_rank", "TF"]].rename(
        columns={"ActDyn_rank": "rank"}
    ).to_csv(actdyn_csv, index=False)

    program_cols = ["ProgramCouple_rank", "TF", "ProgramCouple"]
    program_rank[program_cols].rename(
        columns={"ProgramCouple_rank": "rank"}
    ).to_csv(program_csv, index=False)

    pd.DataFrame({"TF": universe}).to_csv(universe_csv, index=False)

    subprocess.run(
        [
            R_SCRIPT_BIN,
            str(R_RRA_SCRIPT),
            str(network_csv),
            str(actdyn_csv),
            str(program_csv),
            str(universe_csv),
            str(output_csv),
        ],
        check=True,
    )

    out = pd.read_csv(output_csv)
    out["TF"] = out["TF"].astype(str).map(normalize_tf_name)
    return out


# ---------------------------------------------------------------------------
# Final table and driver classes
# ---------------------------------------------------------------------------
def merge_metric_tables(
    rra: pd.DataFrame,
    weighted_rank: pd.DataFrame,
    actdyn_rank: pd.DataFrame,
    program_rank: pd.DataFrame,
) -> pd.DataFrame:
    merged = rra.copy()

    keep_weight = [
        "TF",
        "NetworkSupport_rank",
        "NetworkSupport",
        "weighted_outdegree",
        "weighted_outdegree_abs",
        "n_targets",
        "n_edges",
        "mean_edge_weight",
        "mean_abs_edge_weight",
    ]

    keep_act = [
        "TF",
        "ActDyn_rank",
        "gam_pvalue",
        "FDR_GAM",
        "Amp",
        "Amp_rank_scaled",
        "deviance_explained",
        "DevExpl_rank_scaled",
        "ActDyn",
        "activity_std",
    ]

    keep_prog = [
        "TF",
        "ProgramCouple_rank",
        "ProgramCouple",
        "program_couple_pvalue",
        "target_program_n_targets",
        "program_couple_spline_k",
        "program_couple_grid_size",
        "program_couple_shared_k",
        "target_program_weighted",
        "target_program_expression_source",
        "target_program_expression_preprocessing",
        "status",
    ]

    merged = merged.merge(weighted_rank[keep_weight], on="TF", how="left")
    merged = merged.merge(actdyn_rank[keep_act], on="TF", how="left")
    merged = merged.merge(program_rank[keep_prog], on="TF", how="left")
    return merged


def assign_driver_classes(final_table: pd.DataFrame) -> pd.DataFrame:
    out = final_table.copy()

    n_tfs = len(out)
    top25_cutoff = int(np.ceil(n_tfs * TOP_PERCENTILE_SUPPORT))

    for col in ["NetworkSupport_rank", "ActDyn_rank", "ProgramCouple_rank"]:
        if col not in out.columns:
            out[col] = np.nan

    out["n_top25_metrics"] = (
        (out["NetworkSupport_rank"] <= top25_cutoff).fillna(False).astype(int)
        + (out["ActDyn_rank"] <= top25_cutoff).fillna(False).astype(int)
        + (out["ProgramCouple_rank"] <= top25_cutoff).fillna(False).astype(int)
    )

    out["driver_class"] = "other"

    top_candidate_mask = (
        (out["RRA_FDR"] < HIGH_CONFIDENCE_FDR)
        & (out["n_top25_metrics"] >= 2)
        & (out["ProgramCouple"] >= MIN_PROGRAM_COUPLE)
        & (out["target_program_n_targets"] >= MIN_TARGETS_FOR_PROGRAM)
    )

    high_priority_mask = (
        (out["rank"] <= TOP_N_DRIVER)
        & (out["n_top25_metrics"] >= 2)
        & (out["ProgramCouple"] >= MIN_PROGRAM_COUPLE)
        & (out["target_program_n_targets"] >= MIN_TARGETS_FOR_PROGRAM)
    )

    hub_like_mask = (
        (out["NetworkSupport_rank"] <= TOP_N_DRIVER)
        & (
            (out["ActDyn_rank"] > top25_cutoff)
            | (out["ProgramCouple_rank"] > top25_cutoff)
            | (out["ProgramCouple"] < MIN_PROGRAM_COUPLE)
            | out["ProgramCouple"].isna()
        )
    )

    out.loc[high_priority_mask, "driver_class"] = "high_priority_dynamic_coupled"
    out.loc[top_candidate_mask, "driver_class"] = "top_candidate"
    out.loc[hub_like_mask, "driver_class"] = "hub_like_network_driven"

    out["is_high_confidence_driver"] = out["driver_class"].isin(
        ["top_candidate", "high_priority_dynamic_coupled"]
    )

    return out


# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------
def main() -> None:
    print("=" * 80)
    print("7.2 TF driver identification from baseline GRN")
    print("Using R mgcv for GAM and R RobustRankAggreg for rank aggregation")
    print("=" * 80)

    check_r_scripts()
    cleanup_stale_intermediate_files(OUTPUT_DIR)

    input_dir = find_full_trajectory_input_dir(FULL_TRAJECTORY_INPUT_DIR)
    manifest = read_manifest(input_dir)
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {OUTPUT_DIR}")

    grn = load_baseline_grn(input_dir, manifest)
    weighted_grn = add_edge_weights(grn)

    zero_fraction = weighted_grn.attrs.get("edge_weight_zero_fraction", np.nan)
    print(f"Fraction of zero edge weights after motif clipping: {zero_fraction:.3f}")
    if np.isfinite(zero_fraction) and zero_fraction > 0.70:
        print(
            "WARNING: More than 70% edge weights are zero after motif clipping. "
            "Consider using percentile/rank scaling for motif_score."
        )

    weighted_grn_path = OUTPUT_DIR / "activity_aware_transition_grn_weighted.csv"
    weighted_grn.to_csv(weighted_grn_path, index=False)
    print(f"Weighted baseline GRN: {len(weighted_grn)} edges -> {weighted_grn_path}")

    weighted_rank = rank_network_support(weighted_grn)
    weighted_rank_path = OUTPUT_DIR / "rank_network_support.csv"
    weighted_rank.to_csv(weighted_rank_path, index=False)
    print(f"NetworkSupport ranking: {len(weighted_rank)} TFs -> {weighted_rank_path}")

    rna, activity_df, pseudotime = load_multiome(input_dir, manifest)
    rna_program = prepare_rna_for_target_program(rna)
    candidate_tfs = sorted(set(weighted_grn["tf"]).intersection(activity_df.columns))
    if not candidate_tfs:
        raise ValueError("No overlap between GRN TFs and chromVAR activity TFs.")

    print(f"Candidate TFs with activity: {len(candidate_tfs)}")

    with tempfile.TemporaryDirectory(prefix="tf_driver_7_2_") as tmp_dir:
        tmp_dir = Path(tmp_dir)
        actdyn_rank = compute_actdyn_r_mgcv(
            activity_df,
            pseudotime,
            candidate_tfs,
            work_dir=tmp_dir,
        )
        actdyn_path = OUTPUT_DIR / "rank_actdyn_mgcv.csv"
        actdyn_rank.to_csv(actdyn_path, index=False)
        print(f"ActDyn ranking from R mgcv -> {actdyn_path}")

        program_rank = compute_program_couple_r_mgcv(
            rna=rna_program,
            activity_df=activity_df,
            pseudotime=pseudotime,
            activity_aware_grn=weighted_grn,
            tfs=candidate_tfs,
            work_dir=tmp_dir,
        )
        program_path = OUTPUT_DIR / "rank_program_couple_mgcv.csv"
        program_rank.to_csv(program_path, index=False)
        print(f"ProgramCouple ranking from R mgcv -> {program_path}")

        universe = sorted(
            set(weighted_rank["TF"].astype(str).map(normalize_tf_name))
            | set(actdyn_rank["TF"].astype(str).map(normalize_tf_name))
            | set(program_rank["TF"].astype(str).map(normalize_tf_name))
        )

        rra = run_rra_r_package(
            weighted_rank=weighted_rank,
            actdyn_rank=actdyn_rank,
            program_rank=program_rank,
            universe=universe,
            work_dir=tmp_dir,
        )

    final_table = merge_metric_tables(rra, weighted_rank, actdyn_rank, program_rank)
    final_table = assign_driver_classes(final_table)

    final_path = OUTPUT_DIR / "tf_driver_rra_ranking.csv"
    final_table.to_csv(final_path, index=False)

    high_conf = final_table[final_table["is_high_confidence_driver"]].copy()
    high_conf_path = OUTPUT_DIR / "high_confidence_driver_tfs.csv"
    high_conf.to_csv(high_conf_path, index=False)

    summary = {
        "input_dir": str(input_dir),
        "output_dir": str(OUTPUT_DIR),
        "baseline_edges": int(len(grn)),
        "weighted_edges": int(len(weighted_grn)),
        "candidate_tfs": int(len(candidate_tfs)),
        "min_targets_for_program": MIN_TARGETS_FOR_PROGRAM,
        "top_k_targets_for_program": TOP_K_TARGETS_FOR_PROGRAM,
        "spline_k": SPLINE_K,
        "program_couple_grid_size": GAM_GRID_SIZE,
        "program_expression_source": f"RNA.layers['{RNA_COUNTS_LAYER}']",
        "program_expression_preprocessing": (
            f"full-gene normalize_total(target_sum={RNA_NORMALIZE_TARGET_SUM:g}); "
            "log1p; target-gene z-score across cells; unweighted target mean"
        ),
        "edge_weight_definition": (
            "percentile_rank(log1p(n_peaks)) among n_peaks > 0, else 0; "
            "multiplied by abs(correlation)"
        ),
        "edge_weight_zero_fraction": zero_fraction,
        "network_support_definition": "sum(abs(edge_weight)) per TF",
        "actdyn_definition": (
            "R mgcv GAM: gam(activity ~ s(pseudotime, k), family=gaussian, method='REML'); "
            "ActDyn = rank_scaled(sd(fitted activity over uniform pseudotime grid)) "
            "* rank_scaled(deviance explained by pseudotime)"
        ),
        "program_couple_definition": (
            "For TFs with at least 10 target genes, keep the top-K targets by "
            "edge_weight and fit Gaussian GAMs for both TF activity and unweighted "
            "target program score against pseudotime; "
            "predict both models on the same 100-point pseudotime grid; "
            "ProgramCouple = abs(Pearson correlation of the two predicted curves)"
        ),
        "rra_definition": "R RobustRankAggreg::aggregateRanks with N = candidate TF universe size",
        "high_confidence_fdr_threshold": HIGH_CONFIDENCE_FDR,
        "top_n_driver": TOP_N_DRIVER,
        "top_percentile_support": TOP_PERCENTILE_SUPPORT,
        "min_program_couple": MIN_PROGRAM_COUPLE,
        "high_confidence_tfs": int(len(high_conf)),
    }

    (OUTPUT_DIR / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    print("\nTop integrated TFs:")
    cols = [
        "rank",
        "TF",
        "RRA_pvalue",
        "RRA_FDR",
        "driver_class",
        "is_high_confidence_driver",
        "NetworkSupport_rank",
        "ActDyn_rank",
        "ProgramCouple_rank",
        "ProgramCouple",
        "target_program_n_targets",
    ]
    cols = [c for c in cols if c in final_table.columns]
    print(final_table[cols].head(20).to_string(index=False))

    print(f"\nFinal ranking -> {final_path}")
    print(f"High-confidence TFs -> {high_conf_path}")
    print("=" * 80)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        raise
