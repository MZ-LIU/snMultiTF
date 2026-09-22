"""
1.1_K value determined.py
================
Clustering label-free multi-omics data using the ACM_WEP clustering method,
Determine the optimal K value with assistance from beilv curves and UMAP.

process:
  1. Read the h5mu data after quality control and extract the RNA pattern
  2. Preprocessing + PCA
  3. Clustering PCA embeddings with ACM_WEP
  4. Save beilv curve
  5. Calculate UMAP coordinates and plot clustering results
"""

import os
import sys
import argparse
import json
from pathlib import Path
import numpy as np
import scanpy as sc
import muon as mu
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
CLUSTER_DIR = ROOT / 'ensemble clustering'
for p in (ROOT, CLUSTER_DIR):
    p = str(p)
    if p not in sys.path:
        sys.path.insert(0, p)

from ensemble_clustering.jullei._01_05_01_05_ACM_WEP import cluster_demo_pro


def parse_args():
    parser = argparse.ArgumentParser(description='K value determined')
    parser.add_argument("--config", type=str, help='parameters.json path')
    parser.add_argument("--run_dir", type=str, help='Run output directory')
    parser.add_argument("--result_dir", type=str, help='Results directory')
    return parser.parse_args()


def resolve_paths(args):
    if args.config and args.run_dir:
        with open(args.config, "r", encoding="utf-8") as f:
            params = json.load(f)
        run_dir = Path(args.run_dir)
        data_path = str(run_dir / "mdata_qc.h5mu")
        output_dir = run_dir / "cell_k"
        output_dir.mkdir(parents=True, exist_ok=True)
        return data_path, output_dir
    else:
        data_path = str(Path(__file__).resolve().parent / "results_sc" / "sc_multiome" / "mdata_qc.h5mu")
        output_dir = Path(__file__).resolve().parent / "results_sc" / "sc_multiome" / "sc_multiome" / "K"
        output_dir.mkdir(parents=True, exist_ok=True)
        return data_path, output_dir


N_HVG         = 3000
N_PCS         = 30
CLUSTER_COUNT = 20  # ACM_WEP cluster number
RANDOM_STATE  = 42

# ============================================================
# ============================================================
args = parse_args()
DATA_PATH, OUTPUT_DIR = resolve_paths(args)

print("=" * 60)
print(f"[1/5] Read data: {DATA_PATH}")
file_ext = os.path.splitext(DATA_PATH)[1].lower()
if file_ext == '.h5mu':
    mdata = mu.read_h5mu(DATA_PATH)
elif file_ext == '.h5ad':
    adata = ad.read_h5ad(DATA_PATH)
    mdata = mu.MuData({"rna": adata})
elif file_ext == '.h5':
    mdata = mu.read_10x_h5(DATA_PATH)
else:
    raise ValueError(f"Unsupported file format: {file_ext}")

adata = mdata['rna'].copy()
print(f"  Data dimension: {adata.n_obs} cells x {adata.n_vars} genes (quality controlled)")

# ============================================================
# ============================================================
print('\n[2/5] Preprocessing + PCA...')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG)
adata_hvg = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=N_PCS, random_state=RANDOM_STATE)
X_pca = adata_hvg.obsm['X_pca']
print(f"  PCA completed: {X_pca.shape}")

# ============================================================
# ============================================================
print(f"\n[3/5] ACM_WEP clustering (number of clusters={CLUSTER_COUNT})...")
labels_raw, beilv = cluster_demo_pro(X_pca, CLUSTER_COUNT)
labels = np.array(labels_raw, dtype=int)
print(f"  Clustering completed, total {len(np.unique(labels))} clusters")

plt.figure(figsize=(10, 6))
plt.plot(range(1, len(beilv) + 1), beilv, marker='o', markersize=4)
plt.xlabel('Iteration')
plt.ylabel('Beilv')
plt.title('Beilv Curve')
plt.grid(True, alpha=0.3)
plt.tight_layout()
beilv_path = os.path.join(OUTPUT_DIR, f"beilv_curve_K{CLUSTER_COUNT}.png")
plt.savefig(beilv_path, dpi=150)
plt.close()
print(f"  beilv curve graph has been saved to: {beilv_path}")

# ============================================================
# ============================================================
print('\n[4/5] Calculate UMAP...')
sc.pp.neighbors(adata_hvg, n_neighbors=15, random_state=RANDOM_STATE)
sc.tl.umap(adata_hvg, min_dist=0.1, random_state=RANDOM_STATE)
X_umap = adata_hvg.obsm['X_umap']

# ============================================================
# ============================================================
print('\n[5/5] Draw UMAP...')
unique_clusters = sorted(np.unique(labels))
n_clusters = len(unique_clusters)
cmap = plt.cm.get_cmap('tab20', n_clusters)
cluster_colors = {cid: cmap(i) for i, cid in enumerate(unique_clusters)}

fig, ax = plt.subplots(figsize=(10, 8))
ax.set_title(
    f"ACM_WEP Clustering (K={CLUSTER_COUNT})",
    fontsize=13, fontweight='bold'
)
for cid in unique_clusters:
    mask = labels == cid
    ax.scatter(
        X_umap[mask, 0], X_umap[mask, 1],
        c=[cluster_colors[cid]],
        s=3, alpha=0.6,
        label=f"Cluster {cid} ({mask.sum()} cells)",
    )
ax.set_xlabel("UMAP-1")
ax.set_ylabel("UMAP-2")
ax.legend(loc='upper right', fontsize=7, markerscale=2,
          framealpha=0.7, ncol=2)

plt.tight_layout()
save_path = os.path.join(OUTPUT_DIR, f"UMAP_ACMWEP_K{CLUSTER_COUNT}.png")
plt.savefig(save_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"  UMAP saved to: {save_path}")

print('\n[Done] All processes are completed.')
print(f"  Please check the beilv curve and UMAP to evaluate whether K={CLUSTER_COUNT} is appropriate.")
