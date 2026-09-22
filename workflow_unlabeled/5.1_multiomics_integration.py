"""
5.1_omics data integration.py
===================
Multi-omics data normalization, MOJITOO-CCA integration, and integrated UMAP plots were performed independently.

This script saves the complete intermediate data required for downstream network construction:
  1. mdata_integrated.h5mu
     - RNA/ATAC raw counts are stored in layers['counts']
     - RNA PCA、ATAC LSI
     - mdata.obsm['X_integrated']
     - mdata.obsm['X_umap']
     - mdata.obs['cell_type']
"""

import os
import sys
from pathlib import Path
import warnings

warnings.filterwarnings("ignore")

import anndata as ad
import colorcet as cc
import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import muon as mu
import pandas as pd
import scanpy as sc


# ============================================================
# ============================================================
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, ".."))

if project_root not in sys.path:
    sys.path.insert(0, project_root)
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

import python_scmega as pymega


# ============================================================
# ============================================================
RESULT_BASE_DIR = str(Path(current_dir) / "results_sc" / "sc_multiome")

# os.makedirs(RESULT_DIR,)
MUDATA_PATH = str(Path(RESULT_BASE_DIR) / "mdata_with_celltype.h5mu")

OUTPUT_BASE_DIR = str(Path(RESULT_BASE_DIR) / "network_analysis")

# The downstream 5.2_network.py should read this complete integration file and perform trajectory analysis and network construction separately.
INTEGRATED_MDATA_PATH = str(Path(RESULT_BASE_DIR) / "mdata_integrated.h5mu")

SHOW_PLOTS = False

# ============================================================
# ============================================================
def make_glasbey_palette(categories):
    k = len(categories)
    glasbey_colors = cc.glasbey[:k]
    return [mcolors.to_hex(color) for color in glasbey_colors]


def show_plot_if_enabled():
    if SHOW_PLOTS:
        plt.show(block=False)
        plt.pause(5)
        plt.close()


def _require_counts_layer(mdata, mod):
    if mod not in mdata.mod:
        raise ValueError(f"{mod} modal is missing from MuData.")
    if "counts" not in mdata[mod].layers:
        raise ValueError(
            f"{mod}.layers['counts'] does not exist. 5.1 requires renormalization and integration starting from the original counts."
        )


def _make_counts_adata(mdata, mod):
    _require_counts_layer(mdata, mod)
    src = mdata[mod]
    return ad.AnnData(
        X=src.layers["counts"].copy(),
        obs=src.obs.copy(),
        var=src.var.copy(),
        uns=src.uns.copy(),
        dtype=src.layers["counts"].dtype,
    )


def _sync_cell_type_to_global_obs(mdata):
    if "cell_type" not in mdata["rna"].obs.columns:
        raise ValueError('cell_type is missing in RNA obs, please run 4.2_celltype_writeback.py first.')

    mdata.obs["cell_type"] = mdata["rna"].obs["cell_type"].astype(str).copy()
    cats = sorted(mdata.obs["cell_type"].astype(str).unique())
    mdata.obs["cell_type"] = pd.Categorical(mdata.obs["cell_type"], categories=cats)
    mdata.uns["cell_type_colors"] = make_glasbey_palette(cats)

    for mod in ["rna", "atac"]:
        if mod in mdata.mod:
            common = mdata[mod].obs_names.intersection(mdata.obs_names)
            mdata[mod].obs["cell_type"] = "Unknown"
            mdata[mod].obs.loc[common, "cell_type"] = mdata.obs.loc[common, "cell_type"].astype(str).values
            mdata[mod].obs["cell_type"] = pd.Categorical(
                mdata[mod].obs["cell_type"].astype(str),
                categories=cats,
            )
            mdata[mod].uns["cell_type_colors"] = mdata.uns["cell_type_colors"]

    return cats


def _plot_integrated_umap(mdata, output_dir):
    cats = sorted(mdata.obs["cell_type"].astype(str).unique())
    mdata.obs["cell_type"] = pd.Categorical(mdata.obs["cell_type"].astype(str), categories=cats)
    cmap = plt.cm.get_cmap("tab20", len(cats))
    type_to_color = {cell_type: cmap(i) for i, cell_type in enumerate(cats)}
    mdata.uns["cell_type_colors"] = [
        mcolors.to_hex(type_to_color[cell_type]) for cell_type in cats
    ]

    X_umap = mdata.obsm["X_umap"]
    cell_types = mdata.obs["cell_type"].astype(str).to_numpy()

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    ax.set_title("Integrated UMAP (celltype)", fontsize=13, fontweight="bold")
    for cell_type in cats:
        mask = cell_types == cell_type
        ax.scatter(
            X_umap[mask, 0],
            X_umap[mask, 1],
            c=[type_to_color[cell_type]],
            s=3,
            alpha=0.6,
            label=cell_type,
        )
    ax.set_xlabel("UMAP-1")
    ax.set_ylabel("UMAP-2")
    ax.legend(
        loc="upper right",
        fontsize=6,
        markerscale=2,
        framealpha=0.7,
        ncol=2,
    )
    fig.tight_layout()
    fig.savefig(output_dir / "04_Integrated_UMAP.png", dpi=150, bbox_inches="tight")
    show_plot_if_enabled()
    plt.close(fig)


# ============================================================
# ============================================================
def main():
    print("=" * 70)
    print('5.1 Multi-omics data integration and UMAP')
    print("=" * 70)

    output_dir = Path(RESULT_BASE_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Results saving directory: {output_dir}")

    if not os.path.exists(MUDATA_PATH):
        raise FileNotFoundError(f"MuData file does not exist: {MUDATA_PATH}")

    print('\nLoad MuData with cell_type annotation...')
    file_ext = os.path.splitext(MUDATA_PATH)[1].lower()
    if file_ext == ".h5mu":
        source_mdata = mu.read_h5mu(MUDATA_PATH)
    elif file_ext == ".h5ad":
        adata = ad.read_h5ad(MUDATA_PATH)
        source_mdata = mu.MuData({"rna": adata})
    else:
        raise ValueError(f"Unsupported file format: {file_ext}")
    print(f"Data loading completed: {source_mdata}")

    print('\nRecreate RNA/ATAC objects based on original counts...')
    rna_data = _make_counts_adata(source_mdata, "rna")
    atac_data = _make_counts_adata(source_mdata, "atac")

    if "cell_type" not in atac_data.obs.columns:
        atac_data.obs["cell_type"] = "Unknown"
        common_cells = atac_data.obs_names.intersection(rna_data.obs_names)
        if len(common_cells) > 0:
            atac_data.obs.loc[common_cells, "cell_type"] = (
                rna_data.obs.loc[common_cells, "cell_type"].astype(str).values
            )

    multiome_for_integration = mu.MuData({"rna": rna_data, "atac": atac_data})
    print(f"RNA counts: {rna_data.shape}")
    print(f"ATAC counts: {atac_data.shape}")

    print('\nNormalization and dimensionality reduction for multi-omics integration...')
    multiome_for_integration = pymega.normalize_multiome_muon(
        multiome_for_integration,
        rna_params={
            "n_top_genes": 3000,
            "n_pcs": 50,
            "run_umap": False,
            "umap_dims": 30,
            "random_state": 42,
        },
        atac_params={
            "n_components": 50,
            "run_umap": False,
            "umap_dims": 30,
        },
        inplace=True,
    )

    print('\nPerforming MOJITOO-CCA multi-omics integration...')
    mdata = pymega.coembed_data(
        rna_data=multiome_for_integration["rna"],
        atac_data=multiome_for_integration["atac"],
        rna_reduction="pca",
        atac_reduction="lsi",
        verbose=True,
    )

    cats = _sync_cell_type_to_global_obs(mdata)
    print('\nIntegrated result verification:')
    print(f"  Total cell number: {mdata.n_obs}")
    print(f"  Mode number: {len(mdata.mod)} ({', '.join(mdata.mod.keys())})")
    print(f"  RNA: {mdata['rna'].n_obs} cells x {mdata['rna'].n_vars} genes")
    print(f"  ATAC: {mdata['atac'].n_obs} cells x {mdata['atac'].n_vars} peaks")
    print(f"  Number of cell types: {len(cats)}")
    print(f"  Integration method: {mdata.uns.get('coembedding', {}).get('method', 'unknown')}")

    print('\nCompute integrated UMAP based on X_integrated...')
    sc.pp.neighbors(mdata, use_rep="X_integrated", n_neighbors=30, metric="cosine")
    sc.tl.umap(mdata, min_dist=0.3, spread=1.0, random_state=42)
    print(f"UMAP completed: {mdata.obsm['X_umap'].shape}")
    _plot_integrated_umap(mdata, output_dir)
    print('Integrated UMAP Saved: 04_Integrated_UMAP.png')

    integrated_path = Path(INTEGRATED_MDATA_PATH)
    integrated_path.parent.mkdir(parents=True, exist_ok=True)
    mdata.write_h5mu(integrated_path)
    print(f"\n fully integrates MuData Saved: {integrated_path}")

    summary_lines = [
        '# 5.1 Multi-omics integration results',
        "",
        f"- Input data: {MUDATA_PATH}",
        f"- Output integration MuData: {integrated_path}",
        f"- Total cell number: {mdata.n_obs}",
        f"- RNA: {mdata['rna'].n_obs} cells x {mdata['rna'].n_vars} genes",
        f"- ATAC: {mdata['atac'].n_obs} cells x {mdata['atac'].n_vars} peaks",
        f"- Completion time: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
    ]
    with open(output_dir / "5.1_integration_summary.md", "w", encoding="utf-8") as f:
        f.write("\n".join(summary_lines))

    print('\n5.1 The process is completed. Downstream 5.2_network.py can read mdata_integrated.h5mu and then perform trajectory analysis and network construction.')


if __name__ == "__main__":
    main()
