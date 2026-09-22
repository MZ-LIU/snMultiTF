"""
5.2_network.py
================
Starting from the multi-omics integration results saved in 5.1, execute the trajectory inference and gene regulatory network construction process

Function description:
  1. Read mdata_integrated.h5mu
  2. Trajectory inference analysis
  3. Transcription factor and gene screening (chromVAR)
  4. TF-gene correlation analysis
  5. Regulation network construction and visualization

How to run:
  python 5.2_network.py
"""

import os
import json
import re
import sys
import warnings
warnings.filterwarnings('ignore')

import gzip
import traceback
from pathlib import Path
from datetime import datetime

import pandas as pd

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import scanpy as sc
import muon as mu
import anndata as ad

import colorcet as cc
import matplotlib.colors as mcolors

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
MUDATA_PATH = os.path.join(current_dir, "results_sc", "sc_multiome", "mdata_integrated.h5mu")
OUTPUT_BASE_DIR = os.path.join(current_dir, "results_sc", "sc_multiome", "network_analysis")
os.makedirs(OUTPUT_BASE_DIR , exist_ok=True)

# SwitchTFI input/output folder name under OUTPUT_BASE_DIR.
# Set to None to use the automatic name: switchTFI_inputs_<start>_to_<end>.
# SWITCHTFI_OUTPUT_DIR_NAME = "switchTFI_inputs_CD4T"
SWITCHTFI_OUTPUT_DIR_NAME = None
FULL_TRAJECTORY_OUTPUT_DIR_NAME = "full_trajectory_inputs"
start = "Enterochromaffin cell"
end = "Pancreatic beta cell"
TRAJECTORY_PATH = [start, end]  # Supporting cells->tumor blood vessel-related cells
TRAJECTORY_BACKEND = "archr"  # Set to "archr" to call R ArchR::addTrajectory()."python"
ARCHR_PROJECT_DIR = None       # Optional; if None, ArchRProject is built from FRAGMENTS_PATH.
ARCHR_R_EXECUTABLE = r"E:\program files\R\R-4.5.0\bin\x64\Rscript.exe"
ARCHR_THREADS = 1

CHROMVAR_CONFIG = {
    'motif_database': 'JASPAR2020',  # scMEGA PBMC vignette: library(JASPAR2020)
    'species': None,                 # scMEGA uses tax_group='vertebrates', not species='Homo sapiens'
    'tax_group': 'vertebrates',
    'motif_collection': 'CORE',
    'genome_fasta': 'hg38',
}

FRAGMENTS_PATH = r"E:\data\Pancreas\GSM5979672_week2.atac_fragments.tsv.gz"  

SHOW_PLOTS = False  # Set to True to automatically close after 5 seconds of display, and False to only save without displaying.


# ============================================================
# ============================================================
def make_glasbey_palette(categories):
    """
    Use Glasbey Color Matching - Highly Distinguished Color Matching Designed for a Large Range of Categories
    """
    k = len(categories)
    glasbey_colors = cc.glasbey[:k]
    return [mcolors.to_hex(color) for color in glasbey_colors]


def show_plot_if_enabled():
    """Determine whether to display images based on configuration"""
    if SHOW_PLOTS:
        plt.show(block=False)
        plt.pause(5)
        plt.close()


def sanitize_path_part(value):
    """Convert a cell type/state name into a path-safe token."""
    token = re.sub(r"[^\w]+", "_", str(value).strip(), flags=re.UNICODE)
    return token.strip("_") or "unknown"


def make_switchtfi_output_dir_name(start_state, end_state):
    """Folder format: switchTFI_inputs_<start>_to_<end>."""
    return (
        "switchTFI_inputs_"
        f"{sanitize_path_part(start_state)}_to_{sanitize_path_part(end_state)}"
    )


def save_full_trajectory_comparison_inputs(
    data,
    grn_df,
    output_dir,
    trajectory_path,
    trajectory_key="Trajectory",
    state_key="cell_type",
    label_col_name="switch_label",
):
    """Save all trajectory cells for methods that require continuous pseudotime.

    Unlike ``save_switchtfi_inputs_from_baseline()``, this export retains every
    trajectory state, including intermediate states when present.
    """
    from python_scmega.network_inference.grn_inference import (
        _aggregate_chromvar_activity_by_tf,
        _apply_motif_names_from_atac,
    )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    required_modalities = ("rna", "atac", "chromvar")
    missing_modalities = [name for name in required_modalities if name not in data.mod]
    if missing_modalities:
        raise ValueError(
            'Full trajectory export is missing required modalities: ' + ", ".join(missing_modalities)
        )

    rna_source = data["rna"]
    if state_key not in rna_source.obs.columns:
        raise ValueError(f"RNA.obs['{state_key}'] is missing from full trajectory export.")
    if trajectory_key not in rna_source.obs.columns:
        raise ValueError(f"RNA.obs['{trajectory_key}'] is missing from full trajectory export.")
    if "pseudotime" not in rna_source.obs.columns:
        raise ValueError("Full trajectory export is missing RNA.obs['pseudotime'].")
    if "counts" not in rna_source.layers:
        raise ValueError("Full trajectory export requires RNA.layers['counts'] raw counts.")

    pseudotime = pd.to_numeric(rna_source.obs["pseudotime"], errors="coerce")
    if pseudotime.isna().any():
        raise ValueError(
            f"{int(pseudotime.isna().sum())} cells in the full trajectory object are missing pseudotime."
        )

    observed_states = rna_source.obs[state_key].astype(str)
    missing_states = [state for state in trajectory_path if state not in set(observed_states)]
    if missing_states:
        raise ValueError(
            'Complete track object is missing track state: ' + ", ".join(missing_states)
        )

    selected_obs = rna_source.obs_names
    if not selected_obs.is_unique:
        raise ValueError('Full track RNA obs_names must be unique.')

    rna_out = rna_source.copy()
    rna_out.X = rna_source.layers["counts"].copy()
    rna_out.obs[label_col_name] = pd.Categorical(
        observed_states,
        categories=list(trajectory_path),
        ordered=True,
    )
    rna_out.obs["trajectory_state"] = rna_out.obs[label_col_name].copy()

    missing_atac_obs = selected_obs.difference(data["atac"].obs_names)
    if len(missing_atac_obs) > 0:
        raise ValueError(
            'ATAC mode lacks complete trajectory cells: ' + ", ".join(map(str, missing_atac_obs[:5]))
        )
    atac_out = data["atac"][selected_obs].copy()
    if "counts" not in atac_out.layers:
        raise ValueError("Full trajectory export requires ATAC.layers['counts'] raw counts.")
    atac_out.X = atac_out.layers["counts"].copy()
    atac_out.obs = rna_out.obs.copy()

    missing_chromvar_obs = selected_obs.difference(data["chromvar"].obs_names)
    if len(missing_chromvar_obs) > 0:
        raise ValueError(
            'The chromVAR modality lacks complete trajectory cells: '
            + ", ".join(map(str, missing_chromvar_obs[:5]))
        )
    chromvar_source = data["chromvar"][selected_obs].copy()
    chromvar_source = _apply_motif_names_from_atac(chromvar_source, data)
    chromvar_out = _aggregate_chromvar_activity_by_tf(chromvar_source)
    chromvar_out.obs = pd.DataFrame(index=chromvar_out.obs_names.copy())

    full_trajectory = mu.MuData(
        {
            "RNA": rna_out,
            "ATAC": atac_out,
            "chromvar_activity": chromvar_out,
        }
    )
    for obs_key in (trajectory_key, "pseudotime", label_col_name, "trajectory_state"):
        full_trajectory.obs[obs_key] = rna_out.obs[obs_key].to_numpy()
    full_trajectory.uns["trajectory_analysis"] = {
        trajectory_key: {
            "group_by": state_key,
            "trajectory_groups": list(trajectory_path),
        }
    }
    state_counts = {
        str(state): int(count)
        for state, count in observed_states.value_counts().items()
    }
    full_trajectory.uns["full_trajectory_export"] = {
        "scope": "all_trajectory_cells",
        "intended_methods": ["DrivAER", "RegEnrich", "tf_driver_7_2"],
        "state_key": state_key,
        "trajectory_states": list(trajectory_path),
        "state_counts": state_counts,
    }

    h5mu_path = output_dir / "trajectory_multiome.h5mu"
    grn_path = output_dir / "grn_full_baseline.csv"
    manifest_path = output_dir / "switchtfi_inputs_manifest.json"

    full_trajectory.write_h5mu(h5mu_path)
    grn_df.to_csv(grn_path, index=False)

    manifest = {
        "input_scope": "full_trajectory",
        "intended_methods": ["DrivAER", "RegEnrich", "tf_driver_7_2"],
        "trajectory_key": trajectory_key,
        "trajectory_states": list(trajectory_path),
        "state_key_used": state_key,
        "label_col_name": label_col_name,
        "cell_count": int(full_trajectory.n_obs),
        "state_counts": state_counts,
        "require_multimodal": True,
        "h5mu_keys": {
            "rna_modality": "RNA",
            "atac_modality": "ATAC",
            "tf_activity_modality": "chromvar_activity",
            "rna_matrix": "mod['RNA'].X",
            "atac_matrix": "mod['ATAC'].X",
            "tf_activity_matrix": "mod['chromvar_activity'].X",
            "trajectory_key": trajectory_key,
            "label_key": label_col_name,
            "pseudotime_key": "pseudotime",
        },
        "paths": {
            "trajectory_multiome_h5mu": str(h5mu_path.resolve()),
            "baseline_grn_full_csv": str(grn_path.resolve()),
        },
    }
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    return {
        "output_dir": output_dir,
        "h5mu_path": h5mu_path,
        "grn_path": grn_path,
        "manifest_path": manifest_path,
        "cell_count": int(full_trajectory.n_obs),
        "state_counts": state_counts,
    }


def _iter_fragment_barcodes(fragments_path, max_lines=None):
    """Yield barcodes from a 10x fragments.tsv(.gz) file."""
    path = Path(fragments_path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="ascii", errors="ignore") as handle:
        for line_no, line in enumerate(handle, start=1):
            if max_lines is not None and line_no > max_lines:
                break
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4:
                yield parts[3]


def _diagnose_fragment_barcode_overlap(atac_obs_names, fragments_path, max_lines=500000):
    """Compare ATAC obs names with fragment barcodes, including 10x '-1' suffix variants."""
    atac_cells = set(map(str, atac_obs_names))
    fragment_cells = set(_iter_fragment_barcodes(fragments_path, max_lines=max_lines))
    fragment_cells_no_suffix = {cell.removesuffix("-1") for cell in fragment_cells}
    atac_cells_with_suffix = {
        cell if cell.endswith("-1") else f"{cell}-1"
        for cell in atac_cells
    }

    direct_matches = atac_cells & fragment_cells
    no_suffix_matches = atac_cells & fragment_cells_no_suffix
    with_suffix_matches = atac_cells_with_suffix & fragment_cells

    return {
        "atac_cells": len(atac_cells),
        "fragment_cells_sample": len(fragment_cells),
        "direct_matches": len(direct_matches),
        "no_suffix_matches": len(no_suffix_matches),
        "with_suffix_matches": len(with_suffix_matches),
        "direct_rate": len(direct_matches) / len(atac_cells) if atac_cells else 0.0,
        "no_suffix_rate": len(no_suffix_matches) / len(atac_cells) if atac_cells else 0.0,
        "with_suffix_rate": len(with_suffix_matches) / len(atac_cells) if atac_cells else 0.0,
    }


def _prepare_chromvar_input_for_fragments(mdata_trajectory, fragments_path):
    """
    Signac validates fragment barcodes against ATAC cell names exactly.
    The 10x h5 matrix used here stores cells without '-1', while fragments.tsv.gz
    stores them with '-1'. Use a temporary copy with ATAC obs names suffixed only
    for AddMotifs/RunChromVAR, then copy results back to the original object.
    """
    if fragments_path is None:
        return mdata_trajectory, False

    diagnostic = _diagnose_fragment_barcode_overlap(
        mdata_trajectory["atac"].obs_names,
        fragments_path,
    )
    print('\nFragments barcode compatibility check:')
    print(f"  Fragment cells(sample): {diagnostic['fragment_cells_sample']}")
    print(f"  ATAC cells: {diagnostic['atac_cells']}")
    print(f"  Direct matching cells: {diagnostic['direct_matches']} ({diagnostic['direct_rate']:.1%})")
    print(f"  After removing fragments '-1', match: {diagnostic['no_suffix_matches']} ({diagnostic['no_suffix_rate']:.1%})")
    print(f"  Add '-1' to ATAC and match: {diagnostic['with_suffix_matches']} ({diagnostic['with_suffix_rate']:.1%})")

    if diagnostic["direct_matches"] > 0:
        return mdata_trajectory, False

    if diagnostic["with_suffix_rate"] < 0.1:
        print('  Warning: fragments still barely match ATAC cells, please confirm if FRAGMENTS_PATH belongs to the same batch of data.')
        return mdata_trajectory, False

    print("  10x '-1' suffix difference detected: Create temporary compatible copy for AddMotifs/RunChromVAR.")
    chromvar_input = mdata_trajectory.copy()
    original_atac_names = chromvar_input["atac"].obs_names.astype(str)
    chromvar_input["atac"].obs_names = [
        name if name.endswith("-1") else f"{name}-1"
        for name in original_atac_names
    ]
    chromvar_input["atac"].obs["original_cell_id"] = list(original_atac_names)
    return chromvar_input, True


# ============================================================
# ============================================================
def main():
    print("="*70)
    print('Multi-omics integration and gene regulatory network analysis')
    print("="*70)

    # timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    # output_dir = Path(OUTPUT_BASE_DIR) / timestamp
    # output_dir.mkdir(parents=True, exist_ok=True)
    output_dir = Path(OUTPUT_BASE_DIR)
    print(f" Results saving directory: {output_dir}\n")

    # ============================================================
    # ============================================================
    print("\n" + "="*70)
    print('5. Read multi-omics integration results')
    print("="*70)

    print('\nLoad MuData data generated by 5.1...')
    if not os.path.exists(MUDATA_PATH):
        raise FileNotFoundError(f"MuData file does not exist: {MUDATA_PATH}")

    file_ext = os.path.splitext(MUDATA_PATH)[1].lower()
    if file_ext == '.h5mu':
        mdata = mu.read_h5mu(MUDATA_PATH)
    elif file_ext == '.h5ad':
        adata = ad.read_h5ad(MUDATA_PATH)
        mdata = mu.MuData({"rna": adata})
    else:
        raise ValueError(f"Unsupported file format: {file_ext}")

    print(f" Data loading completed: {mdata}")
    required_mods = ["rna", "atac"]
    missing_mods = [mod for mod in required_mods if mod not in mdata.mod]
    if missing_mods:
        raise ValueError(f"Integrated data is missing required modalities: {missing_mods}")
    if "cell_type" not in mdata.obs.columns:
        raise ValueError('cell_type is missing in integrated data mdata.obs, please check the 5.1 output.')
    if "X_umap" not in mdata.obsm:
        raise ValueError("The integrated data is missing mdata.obsm['X_umap'], please run 5.1_omics data integration.py first.")
    if "X_integrated" not in mdata.obsm:
        raise ValueError("Integrated data is missing mdata.obsm['X_integrated'], please run 5.1_omics data integration.py first.")
    for mod in required_mods:
        if "counts" not in mdata[mod].layers:
            raise ValueError(f"{mod}.layers['counts'] does not exist, and subsequent SwitchTFI export requires the original counts.")

    cats = sorted(mdata.obs["cell_type"].astype(str).unique())
    mdata.obs["cell_type"] = pd.Categorical(mdata.obs["cell_type"].astype(str), categories=cats)
    mdata.uns["cell_type_colors"] = make_glasbey_palette(cats)

    print('\nIntegrated data validation:')
    print(f"  Total cell number: {mdata.n_obs}")
    print(f"  Mode number: {len(mdata.mod)} ({', '.join(mdata.mod.keys())})")
    print(f"  RNA: {mdata['rna'].n_obs} cells × {mdata['rna'].n_vars} genes")
    print(f"  ATAC: {mdata['atac'].n_obs} cells × {mdata['atac'].n_vars} peaks")
    print(f"  X_integrated: {mdata.obsm['X_integrated'].shape}")
    print(f"  X_umap: {mdata.obsm['X_umap'].shape}")
    print(f"  Number of cell types: {len(cats)}")

    # ============================================================
    # ============================================================
    print("\n" + "="*70)
    print('6. Trajectory inference analysis')
    print("="*70)

    trajectory_path = TRAJECTORY_PATH
    print(f"Track path: {' → '.join(trajectory_path)}")

    if 'cell_type' in mdata.obs.columns:
        available_celltypes = mdata.obs['cell_type'].unique()
        print(f"Available cell types: {list(available_celltypes)}")

        missing_types = [ct for ct in trajectory_path if ct not in available_celltypes]
        if missing_types:
            print(f"Warning: The following track cell type was not found: {missing_types}")
    else:
        raise ValueError("The 'cell_type' column was not found in the integrated data. Please complete the cell type annotation first.")

    print('  Trajectory inference using integrated data')
    mdata = pymega.add_trajectory(
        mdata,
        trajectory=trajectory_path,
        group_by='cell_type',
        reduction="X_umap",
        dims=[0, 1],
        pre_filter_quantile=0.9,
        post_filter_quantile=0.9,
        use_all=False,
        dof=250,
        spar=1.0,
        name="Trajectory",
        seed=42,
        backend=TRAJECTORY_BACKEND,
        archr_project_dir=ARCHR_PROJECT_DIR,
        archr_output_dir=output_dir / "archr_trajectory",
        archr_r_executable=ARCHR_R_EXECUTABLE,
        archr_fragments_path=FRAGMENTS_PATH,
        archr_genome=CHROMVAR_CONFIG.get("genome_fasta", "hg38"),
        archr_threads=ARCHR_THREADS
    )

    trajectory_values = mdata.obs['Trajectory']
    trajectory_cells_mask = ~pd.isna(trajectory_values)
    n_trajectory_cells = trajectory_cells_mask.sum()

    print(f"Trajectory inference is completed, number of trajectory cells: {n_trajectory_cells} / {mdata.n_obs}")
    if n_trajectory_cells > 0:
        valid_trajectory_values = trajectory_values[trajectory_cells_mask]
        print(f"Trajectory pseudo-time range: {valid_trajectory_values.min():.1f} - {valid_trajectory_values.max():.1f}")
    else:
        raise ValueError('Trajectory inference failed: no cells were assigned to the trajectory')

    trajectory_numeric = pd.to_numeric(mdata.obs["Trajectory"], errors="coerce")
    pseudotime = pd.Series(index=mdata.obs_names, dtype="float64")
    traj_min = float(trajectory_numeric[trajectory_cells_mask].min())
    traj_max = float(trajectory_numeric[trajectory_cells_mask].max())
    if traj_max > traj_min:
        pseudotime.loc[trajectory_cells_mask] = (
            trajectory_numeric.loc[trajectory_cells_mask] - traj_min
        ) / (traj_max - traj_min)
    else:
        pseudotime.loc[trajectory_cells_mask] = 0.0
    mdata.obs["pseudotime"] = pseudotime
    for mod in ("rna", "atac"):
        mdata[mod].obs["Trajectory"] = mdata.obs.loc[mdata[mod].obs_names, "Trajectory"].to_numpy()
        mdata[mod].obs["pseudotime"] = mdata.obs.loc[mdata[mod].obs_names, "pseudotime"].to_numpy()

    mdata_trajectory = mdata[trajectory_cells_mask].copy()
    mdata_trajectory.uns.setdefault("trajectory_analysis", {})
    mdata_trajectory.uns["trajectory_analysis"]["Trajectory"] = {
        "group_by": "cell_type",
        "trajectory_groups": list(trajectory_path),
    }
    print(f"Track cell subset: {mdata_trajectory.n_obs} cells")
    print(f"  - RNA: {mdata_trajectory['rna'].shape}")
    print(f"  - ATAC: {mdata_trajectory['atac'].shape}")

    # Visualized trajectories: Use temporary AnnData for drawing to prevent MuData string labels from being treated as color values by matplotlib.
    trajectory_plot_obs = mdata.obs.loc[trajectory_cells_mask, ["cell_type", "Trajectory"]].copy()
    trajectory_plot_adata = ad.AnnData(obs=trajectory_plot_obs)
    trajectory_plot_adata.obsm["X_umap"] = mdata.obsm["X_umap"][trajectory_cells_mask.to_numpy(), :].copy()
    trajectory_cats = sorted(trajectory_plot_adata.obs["cell_type"].astype(str).unique())
    trajectory_plot_adata.obs["cell_type"] = pd.Categorical(
        trajectory_plot_adata.obs["cell_type"].astype(str),
        categories=trajectory_cats,
    )
    trajectory_plot_adata.uns["cell_type_colors"] = make_glasbey_palette(trajectory_cats)

    _, axes = plt.subplots(1, 2, figsize=(20, 8))
    sc.pl.umap(
        trajectory_plot_adata,
        color='cell_type',
        ax=axes[0],
        show=False,
        frameon=True,
        legend_loc='on data',
        legend_fontsize=10,
        sort_order=False,
    )
    axes[0].set_title("Trajectory Cells-(cell_type)", fontsize=14)
    axes[0].set_xlabel("UMAP 1", fontsize=12)
    axes[0].set_ylabel("UMAP 2", fontsize=12)

    if 'Trajectory' in trajectory_plot_adata.obs.columns:
        sc.pl.umap(
            trajectory_plot_adata,
            color='Trajectory',
            ax=axes[1],
            show=False,
            frameon=True,
            cmap='viridis',
        )
        axes[1].set_title("Trajectory Cells - Pseudotime", fontsize=14)
        axes[1].set_xlabel("UMAP 1", fontsize=12)
        axes[1].set_ylabel("UMAP 2", fontsize=12)
    else:
        axes[1].text(0.5, 0.5, 'Trajectory information\nnot available',
                    ha='center', va='center', transform=axes[1].transAxes,
                    fontsize=12, color='red')
        axes[1].set_title("Trajectory - Information Missing", fontsize=14)

    plt.tight_layout()
    plt.savefig(output_dir / "05_trajectory_analysis.png", dpi=300, bbox_inches='tight')
    show_plot_if_enabled()

    print('Trajectory inference statistics:')
    print(f"  Track path: {' → '.join(trajectory_path)}")
    print(f"  Number of track cells: {n_trajectory_cells} / {mdata.n_obs}")
    print(f"  Track coverage: {n_trajectory_cells/mdata.n_obs*100:.1f}%")
    if n_trajectory_cells > 0:
        print(f"  Pseudo time range: {valid_trajectory_values.min():.1f} - {valid_trajectory_values.max():.1f}")

    print("="*50)

    # ============================================================
    # ============================================================
    print("\n" + "="*70)
    print('7. Transcription factor and gene screening (using rpy2 to call R’s chromVAR)')
    print("="*70)

    print("\n[chromVAR] Start using rpy2 to call R's chromVAR workflow...")
    print(
        f"Configuration: motif_db={CHROMVAR_CONFIG['motif_database']}, "
        f"tax_group={CHROMVAR_CONFIG['tax_group']}"
    )

    print('\nStep 1: Use rpy2 to call R to load the JASPAR motif database...')
    print(f"  JASPAR version: {CHROMVAR_CONFIG['motif_database']}")
    print(f"  tax_group: {CHROMVAR_CONFIG['tax_group']}")
    print(f"  Collection: {CHROMVAR_CONFIG['motif_collection']}")

    try:
        motif_db_r = pymega.load_jaspar_motifs_r(
            collection=CHROMVAR_CONFIG['motif_collection'],
            species=CHROMVAR_CONFIG['species'],
            tax_group=CHROMVAR_CONFIG['tax_group'],
            version=CHROMVAR_CONFIG['motif_database']
        )
        print(f"  Successfully loaded the JASPAR motif database through rpy2")
    except Exception as e:
        print(f"  Failed to load motif database: {e}")
        print('\nPlease make sure you have the following R packages installed:')
        print("  install.packages('BiocManager')")
        print(f"  BiocManager::install(c('{CHROMVAR_CONFIG['motif_database']}', 'TFBSTools'))")
        raise

    print('\n=== Debug information ===')
    print(f"mdata.uns content: {mdata.uns.keys() if hasattr(mdata, 'uns') else 'N/A'}")
    print(f"mdata_trajectory.uns content: {mdata_trajectory.uns.keys() if hasattr(mdata_trajectory, 'uns') else 'N/A'}")
    print(f"mdata_trajectory['atac'].uns keys: {list(mdata_trajectory['atac'].uns.keys())}")

    fragments_path = FRAGMENTS_PATH
    if fragments_path is None and 'fragments' in mdata_trajectory['atac'].uns:
        fragments_path = mdata_trajectory['atac'].uns['fragments']['path']
    print(f"Externally passed fragments_path: {fragments_path}")
    chromvar_input, using_fragment_compatible_copy = _prepare_chromvar_input_for_fragments(
        mdata_trajectory,
        fragments_path,
    )
    print("=================\n")

    print('\nStep 2: Use rpy2 to call R’s AddMotifs to scan the motif binding sites in the peak...')
    print('  This step will:')
    print('    1. Call R’s AddMotifs function')
    print('    2. Extract the DNA sequence of each peak region in R')
    print('    3. Use motif PWM scan sequence')
    print('    4. Generate motif matching matrix (peak × TF)')

    genome_fasta_path = CHROMVAR_CONFIG['genome_fasta']
    if 'hg38' in genome_fasta_path.lower():
        genome_build = 'hg38'
    elif 'hg19' in genome_fasta_path.lower():
        genome_build = 'hg19'
    else:
        print(f"  Warning: Unable to infer genome version from path {genome_fasta_path}, using default hg38")
        genome_build = 'hg38'

    print(f"  Genome used: {genome_build}")

    try:
        pymega.add_motifs_r(
            chromvar_input,
            motifs=motif_db_r,
            genome=genome_build,
            assay='ATAC',
            fragments_path=fragments_path
        )

        if using_fragment_compatible_copy:
            mdata_trajectory['atac'].uns['motifs'] = chromvar_input['atac'].uns['motifs']

        motif_matching_matrix = mdata_trajectory['atac'].uns['motifs']['motif_matching']
        print(f"  AddMotifs completed: {motif_matching_matrix.shape[0]} peaks × {motif_matching_matrix.shape[1]} motifs")
        print(f"  On average, each peak has {motif_matching_matrix.sum(axis=1).mean():.1f} motif binding sites.")

    except Exception as e:
        print(f"  AddMotifs failed: {e}")
        print('\nPlease make sure you have the following R packages installed:')
        print("  BiocManager::install(c('Signac', 'Seurat', 'BSgenome.Hsapiens.UCSC.hg38'))")
        raise

    print("\nStep 4: Use rpy2 to call R's RunChromVAR to calculate the motif activity score...")
    print('  This step will:')
    print('    1. Call R’s RunChromVAR function')
    print('    2. Using motif matching matrices in R')
    print('    3. Use GC content to match background peaks')
    print('    4. Calculate the activity score of each motif in each cell')

    try:
        pymega.run_chromvar_r(
            chromvar_input,
            genome=genome_build,
            assay='ATAC'
        )

        if using_fragment_compatible_copy:
            if 'motifs' in chromvar_input['atac'].uns:
                mdata_trajectory['atac'].uns['motifs'] = chromvar_input['atac'].uns['motifs']
            if 'chromvar' in chromvar_input.mod:
                mdata_trajectory.mod['chromvar'] = chromvar_input['chromvar']

        if 'chromvar' in mdata_trajectory.mod:
            chromvar_adata = mdata_trajectory['chromvar']
            print(f"  RunChromVAR completed: {chromvar_adata.n_vars} motifs")
            print(f"  Activity matrix dimension: {chromvar_adata.X.shape}")
            print(f"  The results are stored in: mdata_trajectory['chromvar']")
        else:
            raise ValueError('The chromVAR modality was not created. Please check run_chromvar_r().')

        print(f"\n chromvar has been synchronized to the trajectory object")
        print(f"     - Cell number: {chromvar_adata.n_obs}")
        print(f"     - TF number: {chromvar_adata.n_vars}")

    except Exception as e:
        print(f"   RunChromVAR failed: {e}")
        print('\nPlease make sure you have the following R packages installed:')
        print("  BiocManager::install(c('Signac', 'chromVAR'))")
        raise

    print('\nScreening important transcription factors...')
    tf_selection_result = pymega.select_tfs(
        mdata_trajectory,
        correlation_threshold=0.5,
        return_dict=True,
        return_heatmap=True
    )

    tfs_df = tf_selection_result['tfs']
    selected_tfs = tfs_df['tfs'].tolist()
    print(f"Selected {len(selected_tfs)} important transcription factors")
    tfs_df.to_csv(output_dir / "selected_tfs.csv", index=False)
    tf_cor_df = tf_selection_result.get('tf_cor_df')
    if tf_cor_df is not None:
        tf_cor_df.to_csv(output_dir / "tf_cor_df_all.csv", index=False)
        print(f"Saved full TF correlation table: {len(tf_cor_df)} TFs")

    print('\nScreen target genes...')
    gene_selection_result = pymega.select_genes(
        mdata_trajectory,
        var_cutoff_gene=0.9,
        distance_cutoff=2000,
        cor_cutoff=0.0,
        fdr_cutoff=1e-4,
        return_heatmap=True
    )

    peak_gene_df = gene_selection_result['p2g']
    selected_genes = gene_selection_result['selected_genes']
    print(f"   Selected {len(selected_genes)} target genes")

    selected_genes_df = pd.DataFrame({'gene': selected_genes})
    selected_genes_df.to_csv(output_dir / "selected_genes.csv", index=False)

    peak_gene_df.to_csv(output_dir / "peak_gene_links.csv", index=False)
    print(f"Peak-gene connections have been saved, containing {len(peak_gene_df)} connections")

    if 'heatmap' in tf_selection_result and tf_selection_result['heatmap'] is not None:
        tf_selection_result['heatmap'].savefig(
            output_dir / "06_TF_selection_heatmap.png",
            dpi=300, bbox_inches='tight'
        )
        plt.close(tf_selection_result['heatmap'])
        print('   TF selection heatmap has been saved: 06_TF_selection_heatmap.png')
    else:
        print('  Note: heatmap is not generated (needs to be implemented in select_tfs)')

    if 'heatmap' in gene_selection_result and gene_selection_result['heatmap'] is not None:
        gene_selection_result['heatmap'].savefig(
            output_dir / "07_Gene_selection_heatmap.png",
            dpi=300, bbox_inches='tight'
        )
        plt.close(gene_selection_result['heatmap'])
        print('   Gene selection heatmap saved: 07_Gene_selection_heatmap.png')
    else:
        print('  Note: heatmap is not generated (needs to be implemented in select_genes)')

    # ============================================================
    # ============================================================
    print("\n" + "="*70)
    print('9. TF-gene correlation analysis (R-compatible workflow)')
    print("="*70)

    print("\nStep 1: Use the motif activity calculated by chromVAR as the TF activity (tf.assay='chromvar')")

    if 'chromvar' not in mdata_trajectory.mod:
        raise ValueError(
            'chromVAR data does not exist!\n'
            "Please make sure run_chromvar_r() completed successfully and created the mdata['chromvar'] modality"
        )

    print(f"   chromVAR active matrix: {chromvar_adata.X.shape}")
    print(f"   Number of available TFs: {chromvar_adata.n_vars}")

    available_tfs = set(chromvar_adata.var_names)
    selected_tf_set = set(selected_tfs)
    missing_tfs = selected_tf_set - available_tfs
    if missing_tfs:
        print(f"   Warning: {len(missing_tfs)} TFs are missing from chromVAR results")

    print('\nStep 2: Calculate TF-gene correlation (along the trajectory)...')
    print('  This step will:')
    print('    1. Obtain the change pattern of TF activity along the trajectory')
    print('    2. Obtain the change pattern of gene expression along the trajectory')
    print('    3. Calculate the Pearson correlation coefficient between the two')

    tf_gene_result = pymega.get_tf_gene_correlation(
        mdata_trajectory,
        tf_list=selected_tfs,
        gene_list=selected_genes,
        tf_assay="chromvar",
        gene_assay="RNA",
        trajectory_name="Trajectory",
        group_every=1,
        smooth_window=7,
        method='pearson',
        return_heatmap=True,
    )

    tf_gene_df = tf_gene_result['correlation']

    print(f"\n calculation completed: {len(tf_gene_df)} TF-gene correlations")
    print(f"   Number of TF included: {tf_gene_df['tf'].nunique()}")
    print(f"   Number of genes included: {tf_gene_df['gene'].nunique()}")

    tf_gene_df.to_csv(output_dir / "tf_gene_correlation.csv", index=False)
    print(f"   TF-gene correlation saved")

    if tf_gene_result['heatmap'] is not None:
        tf_gene_result['heatmap'].savefig(
            output_dir / "08_GRN_correlation_heatmap.png",
            dpi=300, bbox_inches='tight'
        )
        plt.close(tf_gene_result['heatmap'])
        print('   TF-gene correlation heatmap saved')

    # ============================================================
    # ============================================================
    print("\n" + "="*70)
    print('10. Construct a gene regulatory network (R-compatible GetGRN)')
    print("="*70)

    print('\nStep 1: Prepare three input data for GRN construction...')
    print("="*70)

    print('\n[1] Prepare motif.matching...')
    print("  Python: mdata_trajectory.atac.uns['motifs']['motif_matching']")

    motif_matching = mdata_trajectory['atac'].uns['motifs']['motif_matching'].copy()
    print(f"   Original matrix: {motif_matching.shape[0]} peaks × {motif_matching.shape[1]} motifs")

    print('\n[2] Prepare tf.gene.cor...')

    if 'tf_gene_df' not in locals() or tf_gene_df is None or tf_gene_df.empty:
        raise ValueError('tf_gene_df is undefined or empty! Please check the get_tf_gene_correlation step.')

    print(f"   TF-gene correlation: {len(tf_gene_df)} pair ({tf_gene_df['tf'].nunique()} TFs × {tf_gene_df['gene'].nunique()} genes)")

    print('\n[3] Prepare df.p2g...')

    if 'peak_gene_df' not in locals() or peak_gene_df is None or peak_gene_df.empty:
        raise ValueError('peak_gene_df is undefined or empty! Please check the select_genes step.')

    print(f"   Peak-gene connection: {len(peak_gene_df)} pair ({peak_gene_df['peak'].nunique()} peaks × {peak_gene_df['gene'].nunique()} genes)")

    print('\n[4] Filtering motif.matching...')

    unique_tfs_in_cor = tf_gene_df['tf'].unique()
    available_tfs = [tf for tf in unique_tfs_in_cor if tf in motif_matching.columns]
    if len(available_tfs) == 0:
        raise ValueError(
            f"No TF match! \n"
            f"  TF in tf.gene.cor: {list(unique_tfs_in_cor[:5])}\n"
            f"  TF in motif.matching: {list(motif_matching.columns[:5])}"
        )
    motif_matching = motif_matching[available_tfs]
    print(f"   TF filter: {len(available_tfs)}/{len(unique_tfs_in_cor)} TFs reserved")

    unique_peaks_in_p2g = peak_gene_df['peak'].unique()
    available_peaks = [p for p in unique_peaks_in_p2g if p in motif_matching.index]
    if len(available_peaks) == 0:
        raise ValueError(
            f"No peak match! \n"
            f"  peak in df.p2g: {list(unique_peaks_in_p2g[:5])}\n"
            f"  peak in motif.matching: {list(motif_matching.index[:5])}"
        )
    motif_matching = motif_matching.loc[available_peaks]
    print(f"   Peak filtering: {len(available_peaks)}/{len(unique_peaks_in_p2g)} peaks retained")

    print(f"\n final motif.matching: {motif_matching.shape[0]} peaks × {motif_matching.shape[1]} TFs")
    motif_matching.to_csv(output_dir / "motif_matching.csv")
    print(f"   The motif matching matrix has been saved: motif_matching.csv")
    print("="*70)

    print('\nStep 2: Construct a gene regulatory network (GetGRN)...')
    # switchtfi_dir_name = make_switchtfi_output_dir_name(
    #     trajectory_path[0],
    #     trajectory_path[-1],
    # )
    switchtfi_dir_name = (
        SWITCHTFI_OUTPUT_DIR_NAME
        if SWITCHTFI_OUTPUT_DIR_NAME
        else make_switchtfi_output_dir_name(trajectory_path[0], trajectory_path[-1])
    )
    print(f"  SwitchTFI input directory: {switchtfi_dir_name}")

    switchtfi_output_dir = output_dir / switchtfi_dir_name

    df_grn = pymega.infer_grn(
        motif_matching=motif_matching,
        tf_gene_cor=tf_gene_df,
        peak_gene_links=peak_gene_df,
        correlation_threshold=None,
        return_dict=False,
        save_switchtfi_inputs=False,
    )

    if df_grn is None:
        raise ValueError('infer_grn() did not return a GRN result, df_grn is None. Please check the motif_matching, tf_gene_df and peak_gene_df inputs.')

    if df_grn.empty:
        print('  Warning: Failed to build GRN, possibly because:')
        print('    1. Insufficient Motif matching data')
        print('    2. TF-gene correlation threshold is too high')
        print('    3. Insufficient Peak-gene connection data')
        print('    Please check input data and adjust parameters')
    else:
        print(f"\n GRN is constructed, including {len(df_grn)} control edges")

    n_before_corr_filter = len(df_grn)
    n_removed_corr_missing = 0
    if not df_grn.empty:
        if "correlation" not in df_grn.columns:
            raise ValueError('The GRN result is missing the correlation column, and the correlation missing edges cannot be filtered.')

        df_grn = df_grn[df_grn["correlation"].notna()].reset_index(drop=True)
        n_removed_corr_missing = n_before_corr_filter - len(df_grn)
        print(
            '  Correlation missing edge filtering: '
            f"{n_before_corr_filter} -> {len(df_grn)} "
            f"(Remove article {n_removed_corr_missing})"
        )

        if df_grn.empty:
            raise ValueError(
                'GRN is empty after filtering correlation missing edges. Please check if the TF-gene correlation table matches the motif/peak support.'
            )

    full_trajectory_output_dir = output_dir / FULL_TRAJECTORY_OUTPUT_DIR_NAME
    full_trajectory_export = save_full_trajectory_comparison_inputs(
        data=mdata_trajectory,
        grn_df=df_grn,
        output_dir=full_trajectory_output_dir,
        trajectory_path=trajectory_path,
        trajectory_key="Trajectory",
        state_key="cell_type",
        label_col_name="switch_label",
    )
    print(
        '  Complete trajectory comparison input saved: '
        f"{full_trajectory_export['h5mu_path']} "
        f"({full_trajectory_export['cell_count']} cells)"
    )
    print(f"  Complete trajectory status count: {full_trajectory_export['state_counts']}")

    pymega.save_switchtfi_inputs_from_baseline(
        data=mdata_trajectory,
        grn_df=df_grn,
        output_dir=switchtfi_output_dir,
        trajectory_key='Trajectory',
        start_state=trajectory_path[0],
        end_state=trajectory_path[-1],
        label_col_name='switch_label',
        save_multimodal_h5ad=False,
        save_multimodal_h5mu=True,
        multiome_h5mu_filename='trajectory_multiome.h5mu',
        switchtfi_network_filename='grn_switchTFI.csv',
        save_supporting_tables=False,
        full_grn_filename='grn_full_baseline.csv',
        require_multimodal=True,
    )

    required_switchtfi_files = [
        switchtfi_output_dir / "trajectory_multiome.h5mu",
        switchtfi_output_dir / "grn_switchTFI.csv",
        switchtfi_output_dir / "grn_full_baseline.csv",
        switchtfi_output_dir / "switchtfi_inputs_manifest.json",
    ]
    missing_switchtfi_files = [
        str(path) for path in required_switchtfi_files if not path.exists()
    ]
    if missing_switchtfi_files:
        raise FileNotFoundError(
            "SwitchTFI input files were not saved:\n"
            + "\n".join(f"  - {path}" for path in missing_switchtfi_files)
        )

    df_grn.to_csv(output_dir / 'grn_table.csv', index=False)
    print(f"  The complete GRN has been saved: grn_table.csv ({len(df_grn)} edges)")

    if not df_grn.empty:
        print("\n" + "="*70)
        print('GRN network statistics')
        print("="*70)
        print(f"  Total number of control edges: {len(df_grn)}")
        print(f"  Number of TFs involved: {df_grn['tf'].nunique()}")
        print(f"  Number of target genes involved: {df_grn['gene'].nunique()}")
        print(f"  Average number of target genes per TF: {len(df_grn) / df_grn['tf'].nunique():.1f}")
        print(f"  Average correlation: {df_grn['correlation'].mean():.3f}")
        print(f"  Correlation range: [{df_grn['correlation'].min():.3f}, {df_grn['correlation'].max():.3f}]")

        if 'n_peaks' in df_grn.columns:
            print(f"  Average number of supported peaks per TF-gene pair: {df_grn['n_peaks'].mean():.1f}")

        print('\n  Top 10 regulatory relationships (ordered by relevance):')
        top_cols = ['tf', 'gene', 'correlation']
        if 'n_peaks' in df_grn.columns:
            top_cols.append('n_peaks')

        top_regulations = df_grn.nlargest(10, 'correlation')[top_cols]
        for _, row in top_regulations.iterrows():
            if 'n_peaks' in row:
                print(f"    {row['tf']:15s} → {row['gene']:15s}  (r={row['correlation']:.3f}, peaks={int(row['n_peaks'])})")
            else:
                print(f"    {row['tf']:15s} → {row['gene']:15s}  (r={row['correlation']:.3f})")

        print("="*70)

    # ============================================================
    # ============================================================
    print("\n" + "="*70)
    print('Prepare data for visualization')
    print("="*70)

    if not tfs_df.empty:
        tfs_df_sorted = tfs_df.sort_values('time_point').reset_index(drop=True)
        tfs_timepoint_dict = dict(zip(tfs_df_sorted['tfs'], tfs_df_sorted['time_point']))
        print(f"   TF time point data has been prepared: {len(tfs_timepoint_dict)} TFs")
    else:
        tfs_timepoint_dict = {}
        print('   Warning: TF time point data is not available')

    if not df_grn.empty:
        # grn_for_plot = df_grn[
        #     df_grn['correlation'].notna() &
        #     (df_grn['correlation'] > 0.5)
        # ][['tf', 'gene', 'correlation']].copy()
        # grn_for_plot['weights'] = grn_for_plot['correlation']
        grn_for_plot = df_grn[
            df_grn['correlation'].notna()
        ][['tf', 'gene', 'correlation']].copy()
        grn_for_plot['weights'] = grn_for_plot['correlation']

        grn_for_plot.to_csv(output_dir / "filtered_grn_for_plotting.csv", index=False)
        print(f"   GRN for visualization has been saved: filtered_grn_for_plotting.csv, correlation non-empty edges: {len(grn_for_plot)} strips")
    else:
        grn_for_plot = pd.DataFrame(columns=['tf', 'gene', 'correlation', 'weights'])
        print('   GRN is empty, unable to prepare visualization data')

    print("="*70)

    # ============================================================
    # ============================================================
    print('11. Regulatory network visualization')
    print("="*70)

    if grn_for_plot.empty:
        print('   GRN is empty and the visualization step cannot be performed')
    else:
        print(f"  Start visualizing {len(grn_for_plot)} regulatory edges...")

        try:
            pymega.plot_grn_network_circular(
                grn_network=grn_for_plot,
                tfs_timepoint=tfs_timepoint_dict,
                tf_node_size=500,
                gene_node_size=300,
                inner_radius=0.3,
                outer_radius=0.7,
                tf_radius_spread=0.05,
                gene_radius_spread=0.1,
                show_tf_labels=True,
                edge_alpha=0.15,
                edge_width=0.3,
                figsize=(15, 15),
                title='Gene Regulatory Network',
                seed=42,
                save_path=None
            )

            plt.savefig(output_dir / "09_Final_GRN_Plot_Circular.png",
                    dpi=300, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
            show_plot_if_enabled()
            print('   Circular network diagram has been saved: 09_Final_GRN_Plot_Circular.png')

        except Exception as e:
            print(f"   Network visualization failed: {e}")
            traceback.print_exc()

        print("="*70)

    # ============================================================
    # ============================================================
    print("\n" + "="*70)
    print('12. Generate analysis summary report')
    print("="*70)

    summary_lines = [
        '# PyMEGA gene regulatory network analysis report',
        "",
        '## Data Overview',
        f"- RNA data: {mdata['rna'].shape[0]} cells × {mdata['rna'].shape[1]} genes",
        f"- ATAC data: {mdata['atac'].shape[0]} cells × {mdata['atac'].n_vars} peaks",
        f"- Number of track cells: {mdata_trajectory['rna'].shape[0]}",
        ""
    ]

    if 'cell_type' in mdata['rna'].obs.columns:
        summary_lines.extend([
            '## Cell type annotation',
            f"- Number of identified cell types: {mdata['rna'].obs['cell_type'].nunique()} species",
            ""
        ])

    summary_lines.extend([
        '## Trajectory analysis',
        f"- Track path: {' → '.join(trajectory_path)}",
        f"- Number of track cells: {n_trajectory_cells}",
        "",
        '## Feature filter results',
        f"- Number of transcription factors screened: {len(selected_tfs)}",
        f"- Number of genes screened: {len(selected_genes)}",
        f"- Peak-gene connection number: {len(peak_gene_df)}",
        ""
    ])

    if not df_grn.empty:
        summary_lines.extend([
            '## Gene regulatory network',
            f"- Number of edges to control after construction: {n_before_corr_filter}",
            f"- correlation number of missing edges removed: {n_removed_corr_missing}",
            f"- Final number of control edges: {len(df_grn)}",
            f"- Number of TFs involved: {df_grn['tf'].nunique()}",
            f"- Number of target genes involved: {df_grn['gene'].nunique()}",
            ""
        ])

    summary_lines.extend([
        f"## Analysis completion time",
        f"{pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}"
    ])
    summary_report = "\n".join(summary_lines)

    with open(output_dir / "Analysis_Summary_Report.md", 'w', encoding='utf-8') as f:
        f.write(summary_report)
    print('\n Summary report saved: Analysis_Summary_Report.md')

    print("="*70)
    print(f"The whole process of \n is completed and the results are saved to: {output_dir}")
    print("="*70)


if __name__ == "__main__":
    main()
