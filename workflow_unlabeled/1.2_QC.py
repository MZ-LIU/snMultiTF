"""
1_data quality control.py
================
Perform simultaneous quality control on label-free multi-omics data (RNA + ATAC) to generate data that can be directly used in subsequent processes.
`mdata_qc.h5mu`。

Pre-steps: First run 1_Quality Control Index Calculation.py to view the QC chart, and adjust the filter parameters below accordingly.

process:
  1. Load 10X multiome data (h5 format + fragments file)
  2. Use quality_control_separate to perform synchronous quality control on RNA/ATAC
  3. Generate filtered QC visualization
  4. Save MuData after quality control
"""

import gc
import os
import sys
import argparse
import json
from pathlib import Path

import muon as mu

# Compatible with package import paths when running scripts directly
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_DIR, ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import python_scmega as pymega


def parse_args():
    parser = argparse.ArgumentParser(description='Data quality control')
    parser.add_argument("--config", type=str, help='parameters.json path')
    parser.add_argument("--run_dir", type=str, help='Run output directory')
    parser.add_argument("--result_dir", type=str, help='Results directory')
    return parser.parse_args()


def resolve_paths(args):
    if args.config and args.run_dir:
        with open(args.config, "r", encoding="utf-8") as f:
            params = json.load(f)
        data_path = params.get("h5mu_path", "")
        run_dir = Path(args.run_dir)
        output_dir = run_dir
        output_dir.mkdir(parents=True, exist_ok=True)
        qc = params.get("qc", {})
        return data_path, None, output_dir, qc
    else:
        data_path = os.path.normpath(os.path.join(CURRENT_DIR, "..", "..", "..", "data", "Pancreas", "GSM5979673_week2.filtered_feature_bc_matrix.h5"))
        fragments_path = os.path.normpath(os.path.join(CURRENT_DIR, "..", "..", "..", "data", "Pancreas", "GSM5979672_week2.atac_fragments.tsv.gz"))
        output_dir = Path(CURRENT_DIR) / "results_sc" / "sc_multiome"
        output_dir.mkdir(parents=True, exist_ok=True)
        qc = {}
        return data_path, fragments_path, output_dir, qc


def build_qc_params(qc_from_config: dict) -> tuple:
    """Construct RNA and ATAC quality control parameters, giving priority to the parameters passed in by the online tool."""

    rna_params = {
        "min_nFeature_RNA": qc_from_config.get("min_genes", 2000),
        "max_nFeature_RNA": qc_from_config.get("max_genes", 6000),
        "min_nCount_RNA": qc_from_config.get("min_counts", 200),
        "max_nCount_RNA": qc_from_config.get("max_counts", 25000),
        "max_percent_mt": qc_from_config.get("max_mito_percent", 30.0),
        "min_cells_per_gene": qc_from_config.get("min_cells_per_gene", 10),
        "mt_pattern": "^MT-",
        "exclude_gene_patterns": ["^MT-", "^RP"],
    }
    atac_params = {
        "min_cells_per_peak": 10,
        "add_gene_annotation": False,
        "genome": "hg38",
    }
    return rna_params, atac_params


# =========================
# =========================
def load_mudata(data_path: str, fragments_path: str) -> mu.MuData:
    """Load 10X multiome h5 data and record fragments path."""
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"The data file {data_path} does not exist, please check whether the path is correct.")
    if not os.path.exists(fragments_path):
        raise FileNotFoundError(f"The fragment file {fragments_path} does not exist, please check whether the path is correct.")

    mdata = mu.read_10x_h5(data_path)
    print(f"Included modalities: {mdata.mod.keys()}")

    for mod in mdata.mod:
        if mdata[mod].var_names.duplicated().any():
            n_dup = mdata[mod].var_names.duplicated().sum()
            print(f"  {mod}: {n_dup} duplicate variable names were found and are being removed...")
            mdata[mod].var_names_make_unique()

    mdata["atac"].uns["fragments"] = {
        "path": str(Path(fragments_path).absolute()),
        "genome": "hg38",
    }
    print(f"  Fragments path has been recorded: {mdata['atac'].uns['fragments']['path']}")

    print(f"Raw data overview:")
    print(f"  RNA data shape: {mdata['rna'].shape}")
    print(f"  ATAC data shape: {mdata['atac'].shape}")
    return mdata


# =========================
# =========================
def main() -> None:
    args = parse_args()
    data_path, fragments_path, output_dir, qc_config = resolve_paths(args)
    RNA_QC_PARAMS, ATAC_QC_PARAMS = build_qc_params(qc_config)

    print("=" * 70)
    print('Label-free multi-omics data simultaneous quality control (RNA + ATAC)')
    print("=" * 70)
    print(f"Data path: {data_path}")
    print(f"Output directory: {output_dir}")

    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    print('\n[1/4] Loading 10X multiome data...')
    mdata = load_mudata(data_path, fragments_path)
    rna = mdata["rna"]
    atac = mdata["atac"]

    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    print('\n[2/4] Perform synchronous quality control...')
    print(f"  Before filtration: RNA {rna.shape}, ATAC {atac.shape}")

    gc.collect()

    mdata = pymega.quality_control_separate(
        rna,
        atac,
        min_nFeature_RNA=RNA_QC_PARAMS["min_nFeature_RNA"],
        max_nFeature_RNA=RNA_QC_PARAMS["max_nFeature_RNA"],
        min_nCount_RNA=RNA_QC_PARAMS["min_nCount_RNA"],
        max_nCount_RNA=RNA_QC_PARAMS["max_nCount_RNA"],
        max_percent_mt=RNA_QC_PARAMS["max_percent_mt"],
        min_cells_per_gene=RNA_QC_PARAMS["min_cells_per_gene"],
        mt_pattern=RNA_QC_PARAMS["mt_pattern"],
        exclude_gene_patterns=RNA_QC_PARAMS["exclude_gene_patterns"],
        add_gene_annotation=ATAC_QC_PARAMS["add_gene_annotation"],
        genome=ATAC_QC_PARAMS["genome"],
    )

    print(f"After \n quality control: RNA {mdata['rna'].shape}, ATAC {mdata['atac'].shape}")

    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    print('\n[3/4] Generate filtered QC visualization charts...')
    try:
        plots_dir = output_dir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)
        pymega.plot_qc_metrics(
            mdata["rna"],
            save_path=str(plots_dir / "QC_after_filtering.png"),
        )
        print('  QC visualization chart saved')
    except Exception as e:
        print(f"  QC visualization failed: {e}")
        print('  Continue with subsequent analysis...')

    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    print('\n[4/4] Save MuData after quality control...')
    out_path = output_dir / "mdata_qc.h5mu"

    # Fix non-string data in varmap (compatibility issue of muon writing to h5mu)
    for key in mdata.varmap:
        mdata.varmap[key] = mdata.varmap[key].astype(str)

    mdata.write_h5mu(out_path)
    print(f"  Saved: {out_path}")

    print('\nThe process is complete.')


if __name__ == "__main__":
    main()
