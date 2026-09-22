"""
1_data quality control.py
================
Read multi-omics data, perform quality control and generate data that can be used directly in subsequent processes
`mdata_qc.h5mu`。
"""

import os
import sys
from datetime import datetime
from pathlib import Path

import anndata as ad
import muon as mu
import scanpy as sc

# Compatible with package import paths when running scripts directly
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_DIR, ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from python_scmega.data_processing.quality_control import (
    filter_atac_data,
)


# =========================
# =========================
DATA_PATH = r"E:/data/10x_pbmc_celltype/pbmc_multiome.h5mu"
OUTPUT_DIR = Path('E:\日志\snMultiTF\workflow_labeled\results\pbmc_multiome')
os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"Results saving directory: {OUTPUT_DIR}\n")


# =========================
# =========================
RNA_FILTER_PARAMS = {
    "min_cells_per_gene": 10,
    "exclude_gene_patterns": [r"^MT-", r"^RP"],
}

ATAC_QC_PARAMS = {
    "min_cells": 1,
    "add_gene_annotation": False,
    "genome": "hg38",
}


def _load_mudata(path: str) -> mu.MuData:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Input data not found: {path}")
    ext = os.path.splitext(path)[1].lower()
    if ext == ".h5mu":
        mdata = mu.read_h5mu(path)
    elif ext == ".h5ad":
        adata = ad.read_h5ad(path)
        mdata = mu.MuData({"rna": adata})
    else:
        raise ValueError(f"Unsupported input format: {ext}")

    if "rna" not in mdata.mod:
        raise ValueError('The input data is missing RNA modes.')
    if "atac" not in mdata.mod:
        raise ValueError('The input data is missing ATAC mode.')
    return mdata


def filter_features(
    adata,
    min_cells_per_gene: int = 0,
    exclude_gene_patterns: list = None,
    add_qc_flag: bool = False  # Optional: whether to retain QC marks
):
    """
    Filtering RNA cells and genes (corrected full version)
    Process: Copy data → Cell QC filtering → Gene filtering → Gene pattern filtering → Return clean data
    """
    # The first step: you must first copy the original data and define filtered!
    filtered = adata.copy()

    n_cells_before = filtered.n_obs
    n_genes_before = filtered.n_vars

    # ============================================================
    # ============================================================
    if 'pass_rnaQC' in filtered.obs.columns:
        n_before_cells = filtered.n_obs
        filtered = filtered[filtered.obs['pass_rnaQC'] == True].copy()
        n_after_cells = filtered.n_obs
        if n_before_cells > n_after_cells:
            print(f"[Cell QC filter] Delete cells that failed QC: {n_before_cells} → {n_after_cells} (-{n_before_cells - n_after_cells})")

    # ============================================================
    # ============================================================
    n_genes_current = filtered.n_vars

    if min_cells_per_gene > 0:
        sc.pp.filter_genes(filtered, min_cells=min_cells_per_gene)
        print(f"[Gene filter] min_cells={min_cells_per_gene}: {n_genes_current} → {filtered.n_vars}")
        n_genes_current = filtered.n_vars

    if exclude_gene_patterns:
        for pattern in exclude_gene_patterns:
            mask = ~filtered.var_names.str.contains(pattern, case=False, regex=True)
            filtered = filtered[:, mask]

        print(f"[Gene Exclusion] Remove {exclude_gene_patterns}: {n_genes_current} → {filtered.n_vars}")
        n_genes_current = filtered.n_vars

    # ============================================================
    # ============================================================
    print(f"\n=== Filtering completed ===")
    print(f"Cell number: {n_cells_before} → {filtered.n_obs}")
    print(f"Number of genes: {n_genes_before} → {filtered.n_vars}")

    return filtered


def main() -> None:
    print("=" * 70)
    print('RNA/ATAC quality control')
    print("=" * 70)

    print('\n[1/4] Reading multi-omics data...')
    mdata = _load_mudata(DATA_PATH)
    print(f"  RNA: {mdata['rna'].shape}, ATAC: {mdata['atac'].shape}")

    rna = mdata["rna"].copy()
    atac = mdata["atac"].copy()

    print('\n[2/4] RNA gene filtering (cells are not filtered according to QC indicators)...')
    rna_qc = filter_features(
        rna,
        min_cells_per_gene=RNA_FILTER_PARAMS["min_cells_per_gene"],
        exclude_gene_patterns=RNA_FILTER_PARAMS["exclude_gene_patterns"],
    )
    print(f"  After RNA filtration: {rna_qc.shape}")

    print('\n[3/4] Synchronize ATAC by RNA-preserving cells and perform ATAC peaks filtering...')
    shared_cells = rna_qc.obs_names.intersection(atac.obs_names)
    if len(shared_cells) == 0:
        raise ValueError('After RNA filtration, there are no cells in common with ATAC and cannot be synchronized.')
    atac_synced = atac[shared_cells, :].copy()
    rna_qc = rna_qc[shared_cells, :].copy()
    print(f"  Number of cells after synchronization: {len(shared_cells)}")

    atac_qc = filter_atac_data(
        atac_synced,
        min_cells=ATAC_QC_PARAMS["min_cells"],
        add_gene_annotation=ATAC_QC_PARAMS["add_gene_annotation"],
        genome=ATAC_QC_PARAMS["genome"],
    )
    atac_qc = atac_qc[rna_qc.obs_names, :].copy()
    print(f"  After ATAC filtering: {atac_qc.shape}")

    if "counts" not in rna_qc.layers and "counts" in rna.layers:
        rna_qc.layers["counts"] = rna[rna_qc.obs_names, rna_qc.var_names].layers["counts"].copy()
    if "counts" not in atac_qc.layers and "counts" in atac.layers:
        atac_qc.layers["counts"] = atac[atac_qc.obs_names, atac_qc.var_names].layers["counts"].copy()

    print('\n[4/4] Save mdata_qc.h5mu...')
    out_path = OUTPUT_DIR / "mdata_qc.h5mu"

    mdata_out = mu.MuData({"rna": rna_qc, "atac": atac_qc})
    if "fragments" in atac.uns:
        mdata_out["atac"].uns["fragments"] = atac.uns["fragments"]
    mdata_out.write_h5mu(out_path)

    print(f"  Saved: {out_path}")
    print('\nThe process is complete.')


if __name__ == "__main__":
    main()
