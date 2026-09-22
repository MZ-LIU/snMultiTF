"""
Use the cluster results obtained by clustering to perform difference analysis. The difference analysis method is a two-sided rank sum test.
The result is taken as top10

process:
  1. Directly call scipy.stats.mannwhitneyu(x, y, alternative='two-sided') to obtain the U statistic and p value;
  2. Use scipy.stats.false_discovery_control (BH FDR) for multiple correction;
  3. Screen genes with FDR < 0.05 and log2FC > 0.25;
  4. Sort by p-value in ascending order (main order), log2FC in descending order (order), and take top10.
"""

# ============================================================
# ============================================================
import os
import pandas as pd
import scanpy as sc
import muon as mu
import numpy as np
import sys
import scipy.sparse as sp
import anndata as ad
from typing import Optional
                                              
import scipy.stats
# 

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, ".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

# ============================================================
# ============================================================
OUTPUT_DIR     = os.path.join(current_dir, "results_sc", "sc_multiome", "Clustering_results", "k_8")
STABILITY_PATH = os.path.join(OUTPUT_DIR, "gene_stability_result.csv")
LABEL_PATH     = os.path.join(OUTPUT_DIR, "ACMWEP", "final_labels.csv")
DATA_PATH      = os.path.join(current_dir, "results_sc", "sc_multiome", "mdata_qc.h5mu")

os.makedirs(OUTPUT_DIR, exist_ok=True)

markers_save_path = os.path.join(OUTPUT_DIR, 'cluster_marker_genes_10.csv')

print('Extract genes from stability results...')
df_stability = pd.read_csv(STABILITY_PATH)
upto_genes = (
    df_stability[df_stability['mean_score_when_detected'] > 0.2256].sort_values('mean_score_when_detected', ascending=False)['gene'].tolist()
)
print(f" {len(upto_genes)} qualified stability genes were screened out.")

print('Loading and quality control data...')
if DATA_PATH.endswith(".h5mu"):
    mdata = mu.read_h5mu(DATA_PATH)
else:
    mdata = mu.read_10x_h5(DATA_PATH)

adata_raw = mdata['rna'].copy()
adata_raw.var_names_make_unique()
sc.pp.normalize_total(adata_raw, target_sum=1e4)
sc.pp.log1p(adata_raw)

print('Loading cell labels and filtering cells...')
try:
    df_labels = pd.read_csv(LABEL_PATH)
    if {'cell_name', 'final_label'}.issubset(df_labels.columns):
        df_labels = df_labels[['cell_name', 'final_label']].rename(columns={'final_label': 'cluster'})
    elif {'cell_name', 'cluster'}.issubset(df_labels.columns):
        df_labels = df_labels[['cell_name', 'cluster']]
    else:
        df_labels = pd.read_csv(LABEL_PATH, header=None, names=['cell_name', 'cluster'])
    df_labels['cell_name'] = df_labels['cell_name'].astype(str).str.strip()
    df_labels.set_index('cell_name', inplace=True)
except Exception as e:
    print(f"Failed to read label file: {e}")
    df_labels = pd.read_csv(LABEL_PATH)
    if 'cell_name' in df_labels.columns:
        df_labels.set_index('cell_name', inplace=True)
    else:
        df_labels.index = df_labels.iloc[:, 0].astype(str).str.strip()
        df_labels.columns = ['cluster']

adata_raw.obs['final_cluster'] = (
    adata_raw.obs_names.astype(str).str.strip().map(df_labels.iloc[:, 0])
)

matched_count = adata_raw.obs['final_cluster'].notna().sum()
print(f"   Tag matching completed: {matched_count} / {adata_raw.n_obs} cells have matched tags")
if matched_count == 0:
    print(f"  Example adata cell name: {adata_raw.obs_names[:3].tolist()}")
    print(f"  Example Tag Cell Name: {df_labels.index[:3].tolist()}")

adata_final = adata_raw[
    adata_raw.obs['final_cluster'].notna() &
    (adata_raw.obs['final_cluster'].astype(float).astype(int) != -1)
].copy()

if adata_final.n_obs > 0:
    adata_final.obs['final_cluster'] = (
        adata_final.obs['final_cluster'].astype(float).astype(int).astype(str)
    )
else:
    print('   Fatal error: No cells remaining after filtering, please check tag content or -1 filtering logic.')


def cluster_sort_key(x):
    s = str(x)
    try:
        return (0, int(s))
    except ValueError:
        return (1, s)


# ============================================================
# ============================================================


def two_sided_wilcoxon_one_vs_rest(
    X_dense: np.ndarray,
    gene_names: list,
    labels: np.ndarray,
    target_cluster: str,
) -> pd.DataFrame:
    """
    Do a one-vs-rest two-sided Mann-Whitney U test on target_cluster.

    Call scipy.stats.mannwhitneyu(x, y, alternative='two-sided') directly,
    And use scipy.stats.false_discovery_control to do BH FDR correction.

    return
    ----
    DataFrame contains columns: gene, log2FC, pval, pval_adj
    """
    mask_tgt  = labels == target_cluster
    mask_rest = ~mask_tgt

    X_tgt  = X_dense[mask_tgt, :]
    X_rest = X_dense[mask_rest, :]

    n_genes = X_dense.shape[1]
    pvals = np.ones(n_genes)

    for j in range(n_genes):
        tgt_j  = X_tgt[:, j]
        rest_j = X_rest[:, j]

        if tgt_j.max() == tgt_j.min() == rest_j.max() == rest_j.min():
            continue

        _u_stat, p = scipy.stats.mannwhitneyu(tgt_j, rest_j, alternative='two-sided')
        pvals[j] = p

    fdrs = scipy.stats.false_discovery_control(pvals, method='bh')

    log2fc = (X_tgt.mean(axis=0) - X_rest.mean(axis=0)) / np.log(2)

    df = pd.DataFrame({
        "gene"    : gene_names,
        "log2FC"  : log2fc,
        "pval"    : pvals,
        "pval_adj": fdrs,
    })

    df = df.sort_values("pval_adj", ascending=True).reset_index(drop=True)
    return df



# 1. Only keep candidate genes
candidate_genes = [g for g in upto_genes if g in adata_final.var_names]
adata_candidate = adata_final[:, candidate_genes].copy()
print(f"Candidate gene pool: {len(candidate_genes)}")

MIN_CELLS = 10
cluster_counts = adata_candidate.obs['final_cluster'].value_counts()
valid_clusters  = cluster_counts[cluster_counts >= MIN_CELLS].index.tolist()
adata_candidate = adata_candidate[
    adata_candidate.obs['final_cluster'].isin(valid_clusters)
].copy()
removed = set(cluster_counts.index) - set(valid_clusters)
if removed:
    print(f"  Cluster: {removed} with cell number < {MIN_CELLS} has been removed")
print(f"  There are {adata_candidate.obs['final_cluster'].nunique()} clusters remaining, "
      f"{adata_candidate.n_obs} cells")

# 3. Extract the dense matrix (once to avoid repeated conversions in loops)
if sp.issparse(adata_candidate.X):
    X_dense = adata_candidate.X.toarray()
elif isinstance(adata_candidate.X, np.matrix):
    X_dense = np.asarray(adata_candidate.X)
else:
    X_dense = adata_candidate.X.copy()

labels_arr = adata_candidate.obs['final_cluster'].values.astype(str)
gene_names  = adata_candidate.var_names.tolist()

FDR_THRESHOLD    = 0.05
LOG2FC_THRESHOLD = 0.25
TOP_MARKER_GENES = 10
all_clusters    = sorted(adata_candidate.obs['final_cluster'].unique(), key=cluster_sort_key)

cluster_markers = {}
marker_details  = []

for cluster in all_clusters:
    df_de = two_sided_wilcoxon_one_vs_rest(
        X_dense, gene_names, labels_arr, cluster
    )
    # Step 2: log2FC filtering (only up-regulated genes with log2FC > 0.25 are retained)
    markers = (
        df_de[
            (df_de['pval_adj'] < FDR_THRESHOLD) &
            (df_de['log2FC'] > LOG2FC_THRESHOLD)
        ]
        .sort_values(['pval', 'log2FC'], ascending=[True, False])
        .head(TOP_MARKER_GENES)
    )
    cluster_markers[cluster] = markers['gene'].tolist()
    for _, row in markers.iterrows():
        marker_details.append({
            'cluster' : cluster,
            'gene'    : row['gene'],
            'log2FC'  : row['log2FC'],
            'pval'    : row['pval'],
            'pval_adj': row['pval_adj'],
        })
    print(f"  Cluster {cluster}: {len(markers)} specific marker genes")

final_genes = sorted(set(g for genes in cluster_markers.values() for g in genes))
print(f"\n final specific genes: {len(final_genes)} (from {len(candidate_genes)} candidate genes)")

df_markers = pd.DataFrame(marker_details)
df_markers.to_csv(markers_save_path, index=False)
print(f"Marker gene details saved: {markers_save_path}")

upto_gene_names = final_genes
print('\n Bilateral Wilcoxon differential gene analysis completed!')
