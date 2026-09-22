"""
The whole process_GPU.py
================
Multi-K value consistent clustering main process - GPU accelerated version (Standalone Script)

Function description:
  1. Read h5ad/h5 data
  2. Preprocessing + PCA (CPU)
  3. Execute WARI shared iteration (sampling → PCA → four basic clustering → round-by-round caching)
  4. Supports recovery from shared_iteration_cache and only compensates for missing rounds.
  5. Replay round-by-round cache and build WARI weighted consistency matrix
  6. Summarize gene scores, save results, and perform final clustering

How to run:
  python whole process_GPU.py
"""

import os
import sys
import gc
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
import muon as mu
import anndata as ad

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend to avoid tkinter thread problems and pop-up interrupt scripts
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")  # Block all warnings

# ============================================================
# ============================================================
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, ".."))

if project_root not in sys.path:
    sys.path.insert(0, project_root)
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

# import python_scmega as pymega
# from python_scmega.data_processing.quality_control import filter_rna_cells_features
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]      
CLUSTER_DIR = ROOT / "ensemble_clustering"               

for p in (ROOT, CLUSTER_DIR):
    p = str(p)
    if p not in sys.path:
        sys.path.insert(0, p)
        
from ensemble_clustering.wari.function_wari import (
    get_gpu_info, clear_gpu_memory, GPU_AVAILABLE,
    run_wari_integration, load_consensus_triu_values
)

# ============================================================
# ============================================================
DATA_PATH = str(Path(current_dir) / "results_sc" / "sc_multiome" / "mdata_qc.h5mu")
BASE_DIR  = str(Path(current_dir) / "results_sc" / "sc_multiome")
os.makedirs(BASE_DIR , exist_ok=True)
N_ITER          = 1000  # Number of iterations (it is recommended to set it to 100-1000 for formal running)
SAMPLING_RATIO  = 0.8  # Sample 80% of cells each time
N_TOP_GENES     = 3000  # Number of hypervariable genes
N_PCS           = 30  # PCA principal components
N_PCS_FOR_SCORE = 30  # Number of PCs used for gene scoring
N_JOBS          = 3  # Number of parallel jobs (-1 uses all cores)
K_VALUES = [8]  # The K value range to be tested, such as [8, 9, 10, 11, 12]
CLUSTER_PARAMS = {
    "n_clusters" : 14,  # Number of pruning iterations (subsequent loops will be overwritten by K_VALUE)
    # "resolution" : 0.8, # Control the Seurat algorithm, which has been modified to automatically search based on the incoming k value
    "n_neighbors": 20,  # ← New
}

USE_SHARED_CACHE = True
RESUME = False
RESUME_FROM_DIR = None
USE_GPU = True  # Whether to use GPU acceleration
USE_MEMMAP = False


def main():
    print("="*70)
    print("="*70)

    # timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    # output_dir = Path(BASE_DIR) / timestamp
    # os.makedirs(output_dir, exist_ok=True)
    output_dir = Path(BASE_DIR) / 'Clustering_results'
    os.makedirs(output_dir, exist_ok=True)
    print(f"Results saving directory: {output_dir}\n")

    print("="*60)
    print('GPU environment check')
    print("="*60)
    if GPU_AVAILABLE and USE_GPU:
        print('GPU is available, PyTorch will be used to accelerate matrix operations')
        get_gpu_info()
    else:
        print('Tip: Will run using CPU (GPU is unavailable or disabled)')
    print("="*60 + "\n")

    # ============================================================
    # ============================================================
    if not os.path.exists(DATA_PATH):
        raise FileNotFoundError(f"Data file does not exist: {DATA_PATH}")

    # mdata = mu.read_h5mu(DATA_PATH)
    file_ext = os.path.splitext(DATA_PATH)[1].lower()
    if file_ext == '.h5':
        print(f".h5 file detected, read using mu.read_10x_h5()...")
        mdata = mu.read_10x_h5(DATA_PATH)
    elif file_ext == '.h5mu':
        print(f".h5mu file detected, read using mu.read_h5mu()...")
        mdata = mu.read_h5mu(DATA_PATH)
    elif file_ext == '.h5ad':
        print(f".h5ad file detected, read using ad.read_h5ad()...")
        adata = ad.read_h5ad(DATA_PATH)
        # Wrap AnnData as MuData to keep subsequent code compatible
        mdata = mu.MuData({"rna": adata})
    else:
        raise ValueError(f"Unsupported file formats: {file_ext}, only .h5, .h5mu, .h5ad supported")
    print(f"Data reading completed: {mdata}")

    print(f"RNA: {mdata['rna'].shape}")
    if 'atac' in mdata.mod:
        print(f"ATAC: {mdata['atac'].shape}")

    # ============================================================
    # ============================================================
    # No filtering is repeated here, nor is the counts layer copied; function_wari is in the counts layer
    # If not present, the current RNA.X (the raw count retained after QC) will be used directly.
    adata_qc = mdata['rna'].copy()
    print(f" Directly use the 1.2 quality control results: {adata_qc.n_obs} cells, {adata_qc.n_vars} genes")

    del mdata
    gc.collect()
    print(' Raw data memory released')

    # ============================================================
    # ============================================================
    k_data_plot = {}

    for k_idx, k_value in enumerate(K_VALUES, 1):
        print(f"\n [K value cycle] No. {k_idx}/{len(K_VALUES)} - Processing K={k_value}")
        k_dir = os.path.join(output_dir, f"k_{k_value}")
        os.makedirs(k_dir, exist_ok=True)

        resume_k_dir = None
        if RESUME_FROM_DIR is not None:
            resume_k_dir = os.path.join(RESUME_FROM_DIR, f"k_{k_value}")

        # The final result is still written directly to k_dir, keeping the existing downstream path unchanged.
        run_wari_integration(
            adata=adata_qc,
            base_output_dir=k_dir,
            n_clusters=k_value,
            n_iter=N_ITER,
            true_labels_path=None,
            n_jobs=N_JOBS,
            use_gpu=USE_GPU and GPU_AVAILABLE,
            ratio=SAMPLING_RATIO,
            n_pcs=N_PCS,
            n_pcs_use_for_score=N_PCS_FOR_SCORE,
            n_top_genes=N_TOP_GENES,
            resolution=CLUSTER_PARAMS.get("resolution", 0.8),
            n_neighbors=CLUSTER_PARAMS.get("n_neighbors", 20),
            use_memmap=USE_MEMMAP,
            use_shared_cache=USE_SHARED_CACHE,
            resume=RESUME,
            resume_from_dir=resume_k_dir,
        )

        # Read K-selection metrics from the upper triangular consistency matrix of WARI results without reloading the full n×n matrix.
        consensus_path = os.path.join(k_dir, "consensus_matrix_triu.npy")
        if not os.path.isfile(consensus_path):
            raise FileNotFoundError(f"The WARI consistency matrix of K={k_value} does not exist: {consensus_path}")
        triu_vals = load_consensus_triu_values(consensus_path)
        cdf = np.arange(1, len(triu_vals) + 1) / len(triu_vals)
        auc = np.trapz(y=cdf, x=triu_vals)
        k_data_plot[k_value] = {
            "triu_vals": triu_vals,
            "cdf": cdf,
            "auc": auc,
        }

        if GPU_AVAILABLE: clear_gpu_memory()

    print(f"\n{'='*50}\nThe whole process is completed and the results are saved to: {output_dir}\n{'='*50}")

if __name__ == "__main__":
    main()
