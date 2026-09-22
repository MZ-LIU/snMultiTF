"""W-ARI-only GPU accelerated consistent clustering.

Each round runs Kmeans, Seurat/Leiden, ACM-WEP and
Hierarchical four basic clustering, and use pairwise ARI of four results on the full amount of data
The average is used as the method weight. The file only implements this one result aggregation strategy.

Round-by-round caching, breakpoint recovery, gene scoring, consistency matrix saving, and final clustering output are preserved.
"""
import os
import sys
import random
import numpy as np
import pandas as pd
import scanpy as sc
from tqdm import tqdm
import scipy.sparse as sp
import torch
import anndata as ad
from typing import Optional
from scipy import sparse
import gc
from pathlib import Path

CURRENT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = CURRENT_DIR.parent
if str(CURRENT_DIR) not in sys.path:
    sys.path.insert(0, str(CURRENT_DIR))
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))

try:
    GPU_AVAILABLE = torch.cuda.is_available()
except Exception:
    GPU_AVAILABLE = False


def set_iteration_seed(seed: int):
    """Set per-iteration seeds inside joblib workers for reproducibility."""
    random.seed(seed)
    np.random.seed(seed)
    try:
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed(seed)
            torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.benchmark = False
        torch.backends.cudnn.deterministic = True
    except Exception:
        pass

import matplotlib.pyplot as plt
import os as _os
import colorsys
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from sklearn.metrics import silhouette_score, adjusted_rand_score, normalized_mutual_info_score, calinski_harabasz_score, davies_bouldin_score

# sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), '..', '..'))
# from python_scmega.data_processing.normalization import normalize_rna
def _is_raw_counts(X, sample_size: int = 1000) -> bool:
    """
    Check if matrix may be original counts
    """
    if X is None:
        return False
    if sparse.issparse(X):
        if X.data.size == 0:
            return True
        if X.dtype in [np.int32, np.int64, np.uint32, np.uint64]:
            return True
        data_sample = X.data[:min(sample_size, len(X.data))]
        return np.allclose(data_sample, np.round(data_sample))
    else:
        flat = np.asarray(X).ravel()
        if flat.size == 0:
            return True
        sample = flat[:min(sample_size, len(flat))]
        return np.allclose(sample, np.round(sample))
def normalize_rna(adata: ad.AnnData,
                  target_sum: Optional[float] = None,
                  log_transform: bool = True,
                  scale_data: bool = True,
                  n_top_genes: Optional[int] = None,
                  run_pca: bool = True,
                  n_pcs: Optional[int] = None,
                  run_umap: bool = True,
                  umap_dims: Optional[int] = None,
                  inplace: bool = False) -> ad.AnnData:
    """
    RNA normalization using muon/scanpy - matching Seurat behavior

    Data storage structure:
    - .X: log normalized data
    - .layers['counts']: original counts (retained permanently)
    - .layers['scaled']: scaled matrix (only HVG column is valid)
    - .obsm['X_pca']: PCA dimensionality reduction results
    - .obsm['X_umap']: UMAP dimensionality reduction result
    """
    # Step one: Make sure layers['counts'] exists (only on first call)
    if 'counts' not in adata.layers:
        if not _is_raw_counts(adata.X):
            raise ValueError(
                "First-time normalization requires raw counts in .X, "
                "but detected non-integer values. Please ensure input is raw count matrix."
            )
        adata.layers['counts'] = adata.X.copy()

    target = adata if inplace else adata.copy()

    if target_sum is None:
        target_sum = 1e4  # Default RNA scale factor
    if n_top_genes is None:
        n_top_genes = 3000  # Default number of highly variable genes
    if n_pcs is None:
        n_pcs = 50  # Default number of PCA components
    if umap_dims is None:
        umap_dims = 30

    target.X = target.layers['counts'].copy()

    # Step 1: NormalizeData
    sc.pp.normalize_total(target, target_sum=target_sum)
    if log_transform:
        sc.pp.log1p(target)

    # Step 2: FindVariableFeatures
    sc.pp.highly_variable_genes(
        target,
        n_top_genes=n_top_genes,
        subset=False,
        flavor='seurat'
        # flavor='seurat'
    )

    if scale_data:
        
        # vars_to_regress = [v for v in ['nCount_RNA', 'percent.mt'] if v in target.obs.columns]

        hvg_mask = target.var['highly_variable'].values
        target_hvg = target[:, hvg_mask].copy()

        # if vars_to_regress:
        # # Run regression only on HVG subset, speed up by more than 10 times
        #     sc.pp.regress_out(target_hvg, vars_to_regress)

        sc.pp.scale(target_hvg, max_value=10)
        
        # Stored in .uns (to avoid AnnData’s layers dimension limitation)
        target.uns['scaled_hvg'] = {
            'data': target_hvg.X,
            'genes': target.var_names[hvg_mask].tolist()
        }
        
        print(f"   Scaled {hvg_mask.sum()} HVGs (stored in .uns['scaled_hvg'])")

    if run_pca:
        if not scale_data or 'scaled_hvg' not in target.uns:
            raise ValueError("PCA requires scaled data. Set scale_data=True.")
        
        target_for_pca = target[:, target.var['highly_variable']].copy()
        target_for_pca.X = target.uns['scaled_hvg']['data']
        
        sc.tl.pca(target_for_pca, n_comps=n_pcs, svd_solver='arpack')
        target.obsm['X_pca'] = target_for_pca.obsm['X_pca']
        target.varm['PCs'] = np.zeros((target.n_vars, n_pcs), dtype=np.float32)
        target.varm['PCs'][target.var['highly_variable'].values] = target_for_pca.varm['PCs']
        target.uns['pca'] = target_for_pca.uns['pca']
        print(f"   PCA completed: {n_pcs} components")

    if run_umap and run_pca:
        sc.pp.neighbors(target, n_pcs=umap_dims, use_rep='X_pca', n_neighbors=30, metric='cosine')  
        # sc.pp.neighbors(target, n_pcs=umap_dims, use_rep='X_pca', n_neighbors=30, metric='euclidean')  
        sc.tl.umap(target, min_dist=0.3, spread=1.0, random_state=42)

    return target

from ensemble_clustering.jullei._01_03_01_03_Kmeans import k_meanscsp
from ensemble_clustering.jullei._01_04_01_04_seurat import csp_seurat_from_array
from jullei._01_05_01_05_ACM_WEP import  cluster_demo_pro 
from jullei._01_06_01_06_hierarchical import cengcicsp
try:
    from .leiden_resolution_search import (
        find_leiden_resolution_for_k,
    )
except ImportError:
    from ensemble_clustering.wari.leiden_resolution_search import (
        find_leiden_resolution_for_k,
    )


# ============================================================
# ============================================================

def get_gpu_info():
    """Get GPU information"""
    if not GPU_AVAILABLE:
        return None
    try:
        device_name = torch.cuda.get_device_name(0)
        total_memory = torch.cuda.get_device_properties(0).total_memory
        print(f"GPU device: {device_name}")
        print(f"GPU memory: {total_memory / 1024**3:.2f} GB")
        return device_name
    except Exception as e:
        print(f"Failed to obtain GPU information: {e}")
        return None


def to_gpu(arr):
    """Transfer array to GPU"""
    if not GPU_AVAILABLE:
        return arr
    if isinstance(arr, torch.Tensor):
        return arr.cuda()
    if isinstance(arr, np.ndarray):
        return torch.from_numpy(arr).cuda()
    return torch.tensor(arr, device='cuda')


def to_cpu(arr):
    """Transfer the array back to the CPU"""
    if isinstance(arr, torch.Tensor):
        return arr.detach().cpu().numpy()
    return arr


def clear_gpu_memory():
    """Clean GPU memory"""
    if GPU_AVAILABLE:
        torch.cuda.empty_cache()


# ============================================================
# ============================================================

def random_sampling(adata, ratio: float = 0.8, random_state=None):
    n_cells = adata.n_obs
    n_sample = int(n_cells * ratio)

    rng = np.random.default_rng(random_state)
    sampled_idx = rng.choice(n_cells, size=n_sample, replace=False)
    sampled_idx = np.sort(sampled_idx)

    sampled_adata = adata[sampled_idx].copy()
    return sampled_adata, sampled_idx


# ============================================================
# ============================================================

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

def prepare_sampled_data(sampled_adata,
                         n_top_genes: int = 3000,
                         n_pcs: int = 50,
                         verbose: bool = False):
    adata = sampled_adata

    if verbose:
        print(f"    Standardization + PCA (HVG={n_top_genes}, PCs={n_pcs})...")

    try:
        adata = normalize_rna(
            adata,
            target_sum=1e4,
            log_transform=True,
            scale_data=True,
            n_top_genes=n_top_genes,
            run_pca=True,
            n_pcs=n_pcs,
            run_umap=False,
            inplace=True
        )
        if verbose:
            print(f"  The preprocessing is completed and {np.sum(adata.var['highly_variable'])} HVGs are retained.")
        return adata
    except Exception as e:
        print(f"  Preprocessing failed: {e}")
        raise


def extract_pca_result(adata):
    if 'X_pca' not in adata.obsm:
        raise ValueError('X_pca not detected, please run normalize_rna(run_pca=True) first')

    pca_result = adata.obsm['X_pca']
    hvg_mask = adata.var['highly_variable'].values
    gene_names = adata.var_names[hvg_mask].tolist()

    if 'PCs' in adata.varm:
        components = adata.varm['PCs'][hvg_mask, :].T
    else:
        components = np.eye(len(hvg_mask), pca_result.shape[1])

    if 'pca' in adata.uns and 'variance_ratio' in adata.uns['pca']:
        variance_ratios = adata.uns['pca']['variance_ratio']
    else:
        raise ValueError("Failed to obtain variance: adata.uns['pca']['variance_ratio'] does not exist")

    return pca_result, components, gene_names, variance_ratios

def extract_gene_scores(components, variance_ratios, gene_names, n_pcs_use: int = None, normalize: bool = True, use_gpu: bool = True):
    if use_gpu and GPU_AVAILABLE:
        comp_torch = to_gpu(components.astype(np.float32))
        var_torch = to_gpu(variance_ratios.astype(np.float32))

        if n_pcs_use is not None:
            comp_torch = comp_torch[:n_pcs_use]
            var_torch = var_torch[:n_pcs_use]

        squared_components = comp_torch ** 2
        weighted_components = squared_components * var_torch.unsqueeze(1)
        scores_torch = weighted_components.sum(dim=0)

        if normalize:
            s_min = scores_torch.min()
            s_max = scores_torch.max()
            if s_max > s_min:
                scores_torch = (scores_torch - s_min) / (s_max - s_min)
            else:
                scores_torch = torch.zeros_like(scores_torch)

        scores = to_cpu(scores_torch)
        del comp_torch, var_torch, scores_torch
    else:
        if n_pcs_use is not None:
            components = components[:n_pcs_use]
            variance_ratios = variance_ratios[:n_pcs_use]

        squared_components = components ** 2
        weighted_components = squared_components * variance_ratios[:, np.newaxis]
        scores = weighted_components.sum(axis=0)

        if normalize:
            score_min = scores.min()
            score_max = scores.max()
            if score_max > score_min:
                scores = (scores - score_min) / (score_max - score_min)
            else:
                scores = np.zeros_like(scores)

    return pd.Series(scores, index=gene_names)


def _set_global_leiden_resolution(counts_matrix, var_names, obs_names,
                                  n_clusters: int, n_top_genes: int, n_pcs: int,
                                  output_dir: str, cluster_params: dict = None,
                                  n_neighbors: int = 20, random_state: int = 42,
                                  verbose: bool = True):
    """
    Search Leiden resolution once on the full dataset, save the search log,
    and return cluster_params with the fixed resolution injected.
    """
    cluster_params = dict(cluster_params or {})
    os.makedirs(output_dir, exist_ok=True)

    print(f"[Leiden] Use full data to search for K={n_clusters} resolution...")
    full_adata = ad.AnnData(
        X=counts_matrix,
        var=pd.DataFrame(index=var_names),
        obs=pd.DataFrame(index=obs_names),
    )
    full_adata.layers['counts'] = full_adata.X.copy()
    full_adata = prepare_sampled_data(
        full_adata,
        n_top_genes=n_top_genes,
        n_pcs=n_pcs,
        verbose=verbose,
    )
    full_pca_result, _, _, _ = extract_pca_result(full_adata)
    full_data_eval = np.ascontiguousarray(full_pca_result)

    best_resolution, actual_k, search_log = find_leiden_resolution_for_k(
        data_eval=full_data_eval,
        target_k=n_clusters,
        n_neighbors=n_neighbors,
        random_state=random_state,
        verbose=verbose,
    )
    if actual_k == -1:
        raise RuntimeError('Leiden resolution search failed: all candidate resolutions failed to cluster successfully')

    cluster_params.update({
        'resolution': best_resolution,
        'leiden_resolution': best_resolution,
        'leiden_resolution_actual_k': actual_k,
        'n_neighbors': n_neighbors,
        'random_state': random_state,
    })

    log_file = os.path.join(output_dir, "leiden_resolution_search_log.csv")
    log_df = pd.DataFrame(search_log)
    log_df.insert(0, "search_scope", "full_data")
    log_df.insert(1, "n_neighbors", n_neighbors)
    log_df["selected_resolution"] = best_resolution
    log_df["selected_actual_k"] = actual_k
    log_df["selected"] = (
        (log_df["resolution"] == best_resolution)
        & (log_df["actual_k"] == actual_k)
    )
    log_df.to_csv(log_file, index=False, encoding="utf-8-sig")

    print(
        f"[Leiden] Full search completed: resolution={best_resolution:.4f}, "
        f"Actual K={actual_k}, target K={n_clusters}"
    )
    print(f"[Leiden] Search log saved: {log_file}")

    del full_adata, full_pca_result, full_data_eval
    gc.collect()
    return cluster_params, best_resolution, actual_k, search_log


# ============================================================
# ============================================================

def init_count_files(n_cells: int, output_dir: str = None, use_memmap: bool = False, use_gpu: bool = True):
    """
    Initialize the co-clustering count matrix.

    Fix: dtype changed from int16 (max 32767) to int32 (max ~2.1 billion),
    Avoid integer overflow when the number of iterations is large or the sampling rate is high.
    """
    if use_gpu and GPU_AVAILABLE:
        # Fix 2: int16 → int32 to prevent underflow at high iteration times
        co_cluster = torch.zeros((n_cells, n_cells), dtype=torch.int32, device='cuda')
        co_sample  = torch.zeros((n_cells, n_cells), dtype=torch.int32, device='cuda')
        print(f"  GPU dense counting matrix initialized ({n_cells} × {n_cells}, int32)")
    else:
        # Fix 2: uint16 → int32
        dtype = np.int32
        if use_memmap and output_dir is not None:
            os.makedirs(output_dir, exist_ok=True)
            co_cluster = np.memmap(
                os.path.join(output_dir, "co_cluster_memmap.npy"),
                dtype=dtype, mode='w+', shape=(n_cells, n_cells)
            )
            co_sample = np.memmap(
                os.path.join(output_dir, "co_sample_memmap.npy"),
                dtype=dtype, mode='w+', shape=(n_cells, n_cells)
            )
        else:
            co_cluster = np.zeros((n_cells, n_cells), dtype=dtype)
            co_sample  = np.zeros((n_cells, n_cells), dtype=dtype)
        print(f"    Dense counting matrix initialized ({n_cells} × {n_cells}, int32)")

    return co_cluster, co_sample


def save_count_matrices(co_cluster, co_sample, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)
    cluster_path = os.path.join(output_dir, "co_cluster.npy")
    sample_path  = os.path.join(output_dir, "co_sample.npy")

    np.save(cluster_path, to_cpu(co_cluster))
    np.save(sample_path,  to_cpu(co_sample))

    print(f"       Dense counting matrix saved: {cluster_path}")
    print(f"       Dense counting matrix saved: {sample_path}")


def load_count_matrices(output_dir: str, use_gpu: bool = True):
    cluster_path = os.path.join(output_dir, "co_cluster.npy")
    sample_path  = os.path.join(output_dir, "co_sample.npy")

    if not os.path.exists(cluster_path) or not os.path.exists(sample_path):
        raise FileNotFoundError('The dense counting matrix file does not exist and cannot be resumed.')

    co_cluster = np.load(cluster_path)
    co_sample  = np.load(sample_path)

    if use_gpu and GPU_AVAILABLE:
        # Fix 2: also use int32 when loading
        co_cluster = to_gpu(co_cluster.astype(np.int32))
        co_sample  = to_gpu(co_sample.astype(np.int32))
        print(f"       GPU dense counting matrix loaded")
    else:
        print(f"       Dense counting matrix loaded")

    return co_cluster, co_sample


def update_count_matrices(co_cluster, co_sample, sampled_idx, labels):
    if GPU_AVAILABLE and isinstance(co_cluster, torch.Tensor):
        labels_torch      = to_gpu(labels)
        sampled_idx_torch = to_gpu(sampled_idx.astype(np.int64))

        valid_mask   = labels_torch != -1
        valid_idx    = sampled_idx_torch[valid_mask]
        valid_labels = labels_torch[valid_mask]

        if len(valid_idx) < 2:
            return

        co_sample[valid_idx.unsqueeze(1), valid_idx.unsqueeze(0)] += 1

        unique_clusters = torch.unique(valid_labels)
        for cluster_id in unique_clusters:
            cluster_cell_idx = valid_idx[valid_labels == cluster_id]
            if len(cluster_cell_idx) > 1:
                co_cluster[cluster_cell_idx.unsqueeze(1), cluster_cell_idx.unsqueeze(0)] += 1

        del labels_torch, sampled_idx_torch
    else:
        valid_mask   = labels != -1
        valid_idx    = sampled_idx[valid_mask]
        valid_labels = labels[valid_mask]

        if len(valid_idx) < 2:
            return

        co_sample[np.ix_(valid_idx, valid_idx)] += 1

        unique_clusters = np.unique(valid_labels)
        for cluster_id in unique_clusters:
            cluster_cell_idx = valid_idx[valid_labels == cluster_id]
            if len(cluster_cell_idx) > 1:
                co_cluster[np.ix_(cluster_cell_idx, cluster_cell_idx)] += 1


# ============================================================
# ============================================================

METHOD_NAMES = ['Kmeans', 'Seurat', 'ACM_WEP', 'Hierarchical']










def _update_wari_total_counts(total_cluster, total_sample, sampled_idx, all_labels, weights):
    """
    Total weighted counts are directly accumulated by W-ARI method weights.

    total_cluster[i,j] = Σ_m w_m * I(labels_m[i] == labels_m[j])
    total_sample[i,j] = Σ_m w_m * I(i,j are both valid in this method)
    """
    for labels, w in zip(all_labels, weights):
        if w == 0:
            continue

        valid_mask = labels != -1
        v_idx = sampled_idx[valid_mask]
        v_labels = labels[valid_mask]

        if len(v_idx) < 2:
            continue

        total_sample[np.ix_(v_idx, v_idx)] += w

        for cid in np.unique(v_labels):
            cidx = v_idx[v_labels == cid]
            if len(cidx) > 1:
                total_cluster[np.ix_(cidx, cidx)] += w


def build_wari_consensus_from_totals(total_cluster, total_sample, n_cells):
    """Construct a consistency matrix from W-ARI total weighted homocluster counts/total weighted co-occurrence counts."""
    print('[3/5] Construct a weighted consistency matrix...')
    consensus = np.divide(
        total_cluster,
        total_sample,
        out=np.zeros((n_cells, n_cells), dtype=np.float32),
        where=total_sample != 0
    )
    np.fill_diagonal(consensus, 1.0)
    print(f"  The weighted consistency matrix is constructed ({n_cells} × {n_cells})")
    return consensus


def _save_wari_iteration_cache(cache_dir, result):
    """Saves the sampling index, four algorithm labels, and gene scores required for one round of W-ARI."""
    iter_idx = result['iter_idx']
    sampled_idx = result['sampled_idx']
    method_labels = result['method_labels']

    if sampled_idx is None or not method_labels:
        return None

    method_names = [name for name in METHOD_NAMES if name in method_labels]
    scores = result.get('scores')
    payload = {
        'sampled_idx': np.asarray(sampled_idx, dtype=np.int32),
        'method_names': np.asarray(method_names, dtype='<U32'),
    }

    for name in method_names:
        payload[f'labels_{name}'] = np.asarray(method_labels[name], dtype=np.int32)

    if scores is not None:
        payload['score_genes'] = np.asarray(scores.index.astype(str), dtype='<U128')
        payload['score_values'] = np.asarray(scores.values, dtype=np.float32)

    path = os.path.join(cache_dir, f"iter_{iter_idx:06d}.npz")
    os.makedirs(cache_dir, exist_ok=True)
    np.savez_compressed(path, **payload)
    return path


def _load_wari_iteration_cache(path):
    """
    Read a round of W-ARI cache.

    Non-W-ARI additional fields in legacy caches are safely ignored, so there are
    shared_iteration_cache can still be used.
    """
    with np.load(path, allow_pickle=False) as data:
        sampled_idx = data['sampled_idx']
        method_names = [str(x) for x in data['method_names']]
        method_labels = {
            name: data[f'labels_{name}']
            for name in method_names
        }
        scores = None
        if 'score_genes' in data.files and 'score_values' in data.files:
            scores = pd.Series(
                data['score_values'].astype(np.float32),
                index=[str(x) for x in data['score_genes']]
            )

    iter_idx = int(os.path.basename(path).split('_')[1].split('.')[0])
    return {
        'iter_idx': iter_idx,
        'sampled_idx': sampled_idx,
        'scores': scores,
        'method_labels': method_labels,
        'all_methods_ok': all(name in method_labels for name in METHOD_NAMES),
        'method_errors': {},
    }










def _precompute_w_ari_global_weights(counts_matrix, var_names, obs_names,
                                     n_clusters, n_top_genes, n_pcs,
                                     resolution=0.8, n_neighbors=20):
    """
    Consistent with single-mode w_ari: run four algorithms on the full data, using pairwise ARI means to normalize to weights.

    resolution uses the fixed value obtained by pre-searching the entire amount of data. Seurat/Leiden does not repeat the search here.
    """
    print('[w_ari] Precompute ARI weights for four algorithms on global data...')
    full_adata = ad.AnnData(
        X=counts_matrix,
        var=pd.DataFrame(index=var_names),
        obs=pd.DataFrame(index=obs_names)
    )
    full_adata.layers['counts'] = full_adata.X.copy()

    try:
        full_adata_prep = prepare_sampled_data(
            full_adata,
            n_top_genes=n_top_genes,
            n_pcs=n_pcs,
            verbose=False
        )
        pca_full, _, _, _ = extract_pca_result(full_adata_prep)
        data_full = np.ascontiguousarray(pca_full)

        gl_kmeans = k_meanscsp(data_full, n_clusters, random_state=42)
        gl_seurat, _, _ = csp_seurat_from_array(
            data_full,
            resolution=float(resolution),
            n_neighbors=n_neighbors,
            method='leiden',
            random_state=42,
        )
        gl_seurat = np.asarray(gl_seurat)
        gl_acm, _ = cluster_demo_pro(data_full, n_clusters)
        gl_hc = cengcicsp(data_full, n_clusters, 'ward')
        gl_labels = [gl_kmeans, gl_seurat, gl_acm, gl_hc]

        ari_matrix = np.zeros((4, 4))
        for ai in range(4):
            for bi in range(4):
                if ai != bi:
                    la, lb = gl_labels[ai], gl_labels[bi]
                    valid = (la != -1) & (lb != -1)
                    if valid.sum() > 0:
                        ari_matrix[ai, bi] = adjusted_rand_score(la[valid], lb[valid])

        raw_w = ari_matrix.sum(axis=1) / 3.0
        raw_w = np.clip(raw_w, 0.0, None)
        total = raw_w.sum()
        if total <= 0:
            raise ValueError('The sum of the non-negative average ARI weights of the four base clusters is 0')
        weights = (raw_w / total).tolist()
        print(f"  [w_ari] Global weight: Kmeans={weights[0]:.4f}, "
              f"Seurat={weights[1]:.4f}, "
              f"ACM_WEP={weights[2]:.4f}, "
              f"Hierarchical={weights[3]:.4f}")
        del full_adata, full_adata_prep, data_full
        return weights
    except Exception as e:
        raise RuntimeError(
            'W-ARI global weight calculation failed; to avoid silently switching to other aggregation strategies, this run was terminated.'
        ) from e


def _run_wari_iteration(iter_idx, counts_matrix, var_names, obs_names,
                        n_pcs, n_pcs_use_for_score,
                        n_top_genes, ratio, n_clusters,
                        resolution, n_neighbors,
                        random_state_base, use_gpu: bool = True):
    """
    Performs one round of W-ARI basic calculations: fixed subsampling, one preprocessing, four basic clusterings, and gene scoring.

    The four basic clusterings together form the input to W-ARI.
    Resolution is a fixed value obtained by pre-searching the entire amount of data, and this value is reused in each round.
    """
    random_state = random_state_base + iter_idx
    set_iteration_seed(random_state)

    try:
        n_cells = counts_matrix.shape[0]
        if n_cells == 0:
            raise ValueError("counts_matrix is empty")

        n_sample = int(n_cells * ratio)
        rng = np.random.default_rng(random_state)
        sampled_idx = np.sort(
            rng.choice(n_cells, size=n_sample, replace=False)
        )

        sampled_adata = ad.AnnData(
            X=counts_matrix[sampled_idx],
            var=pd.DataFrame(index=var_names),
            obs=pd.DataFrame(index=obs_names[sampled_idx])
        )
        sampled_adata.layers['counts'] = sampled_adata.X.copy()

        sampled_adata = prepare_sampled_data(
            sampled_adata,
            n_top_genes=n_top_genes,
            n_pcs=n_pcs,
            verbose=False,
        )
        pca_result, components, gene_names, variance_ratios = extract_pca_result(sampled_adata)
        data_eval = np.ascontiguousarray(pca_result)

        method_labels = {}
        method_errors = {}

        def _save_labels(method_name, labels):
            labels = np.asarray(labels)
            if len(labels) != sampled_adata.n_obs:
                raise ValueError(
                    f"Cluster output mismatch for {method_name}: "
                    f"labels={len(labels)}, cells={sampled_adata.n_obs}"
                )
            method_labels[method_name] = labels

        try:
            _save_labels('Kmeans', k_meanscsp(data_eval, n_clusters, random_state=42))
        except Exception as e:
            method_errors['Kmeans'] = str(e)

        try:
            labels_seurat, _, _ = csp_seurat_from_array(
                data_eval,
                resolution=float(resolution),
                n_neighbors=n_neighbors,
                method='leiden',
                random_state=42,
            )
            _save_labels('Seurat', labels_seurat)
        except Exception as e:
            method_errors['Seurat'] = str(e)

        try:
            labels_acm, _ = cluster_demo_pro(data_eval, n_clusters)
            _save_labels('ACM_WEP', labels_acm)
        except Exception as e:
            method_errors['ACM_WEP'] = str(e)

        try:
            _save_labels('Hierarchical', cengcicsp(data_eval, n_clusters, 'ward'))
        except Exception as e:
            method_errors['Hierarchical'] = str(e)

        scores = extract_gene_scores(
            components,
            variance_ratios,
            gene_names,
            n_pcs_use=n_pcs_use_for_score,
            use_gpu=use_gpu,
        )

        return {
            'iter_idx': iter_idx,
            'sampled_idx': sampled_idx,
            'scores': scores,
            'method_labels': method_labels,
            'all_methods_ok': all(name in method_labels for name in METHOD_NAMES),
            'method_errors': method_errors,
        }

    except Exception as e:
        print(f"\n W-ARI iteration {iter_idx} failed: {e}")
        import traceback
        traceback.print_exc()
        return {
            'iter_idx': iter_idx,
            'sampled_idx': None,
            'scores': None,
            'method_labels': {},
            'all_methods_ok': False,
            'method_errors': {'iteration': str(e)},
        }





def run_wari_integration(adata, base_output_dir, n_clusters,
                       n_iter=100, true_labels_path=None,
                       n_jobs=-1, use_gpu=True,
                       ratio=0.8, n_pcs=30,
                       n_pcs_use_for_score=30, n_top_genes=3000,
                       resolution=0.8, n_neighbors=20,
                       use_memmap=False,
                       use_shared_cache=False,
                       resume=False,
                       resume_from_dir=None):
    """
    W-ARI unique summary entry.

    Each round runs Kmeans, Seurat/Leiden, ACM-WEP and
    Hierarchical four basic clustering, and then use the pairwise ARI obtained from the full amount of data
    Average weight summary consistency matrix. The final result is written directly to base_output_dir.
    """
    from joblib import Parallel, delayed
    import shutil

    n_cells = adata.n_obs

    if 'counts' in adata.layers:
        counts_matrix = adata.layers['counts']
    else:
        counts_matrix = adata.X
    if sp.issparse(counts_matrix):
        counts_matrix = counts_matrix.toarray()

    var_names = np.asarray(adata.var_names)
    obs_names = np.asarray(adata.obs_names)

    _, best_leiden_resolution, leiden_actual_k, _ = _set_global_leiden_resolution(
        counts_matrix=counts_matrix,
        var_names=var_names,
        obs_names=obs_names,
        n_clusters=n_clusters,
        n_top_genes=n_top_genes,
        n_pcs=n_pcs,
        output_dir=base_output_dir,
        cluster_params={
            "n_clusters": n_clusters,
            "n_neighbors": n_neighbors,
            "random_state": 42,
        },
        n_neighbors=n_neighbors,
        random_state=42,
        verbose=True,
    )
    resolution = best_leiden_resolution

    print(f"\n{'='*55}")
    print('  W-ARI Weighted Consistency Summary')
    print(f"  Number of iterations: {n_iter}, basic algorithm: {', '.join(METHOD_NAMES)}")
    print(f"{'='*55}")

    if use_gpu and GPU_AVAILABLE:
        get_gpu_info()

    global_weights = _precompute_w_ari_global_weights(
        counts_matrix=counts_matrix,
        var_names=var_names,
        obs_names=obs_names,
        n_clusters=n_clusters,
        n_top_genes=n_top_genes,
        n_pcs=n_pcs,
        resolution=resolution,
        n_neighbors=n_neighbors,
    )

    shared_results = []
    failure_records = []
    cache_paths = {}

    cache_output_dir = os.path.join(base_output_dir, "shared_iteration_cache")
    cache_source_dir = None
    if use_shared_cache:
        cache_source_dir = (
            os.path.join(resume_from_dir, "shared_iteration_cache")
            if resume and resume_from_dir
            else cache_output_dir
        )
        os.makedirs(cache_output_dir, exist_ok=True)

    missing_iters = list(range(n_iter))
    if use_shared_cache and resume and cache_source_dir and os.path.isdir(cache_source_dir):
        missing_iters = []
        for iter_idx in range(n_iter):
            cache_path = os.path.join(cache_source_dir, f"iter_{iter_idx:06d}.npz")
            if os.path.isfile(cache_path):
                cache_paths[iter_idx] = cache_path
            else:
                missing_iters.append(iter_idx)
        print(
            f"  W-ARI round-by-round cache: {len(cache_paths)}/{n_iter} round found,"
            f"Missing {len(missing_iters)} iterations"
        )
    elif use_shared_cache:
        print(f"  W-ARI round-by-round cache: generated from scratch, save directory {cache_output_dir}")

    parallel_gen = Parallel(n_jobs=n_jobs, return_as="generator")(
        delayed(_run_wari_iteration)(
            iter_idx,
            counts_matrix,
            var_names,
            obs_names,
            n_pcs,
            n_pcs_use_for_score,
            n_top_genes,
            ratio,
            n_clusters,
            resolution,
            n_neighbors,
            42,
            use_gpu=use_gpu,
        )
        for iter_idx in missing_iters
    )

    print('\n[2/5] W-ARI round-by-round basic clustering calculation...')
    for result in tqdm(
        parallel_gen,
        total=len(missing_iters),
        desc='W-ARI calculation round by round',
        unit="iter",
    ):
        iter_idx = result['iter_idx']
        method_errors = result.get('method_errors', {})
        if method_errors:
            failure_records.append({
                'iter_idx': iter_idx,
                'failed_methods': ';'.join(method_errors.keys()),
                'error_detail': ' | '.join(
                    f"{method}: {message}"
                    for method, message in method_errors.items()
                ),
            })

        if result['sampled_idx'] is None or result['scores'] is None:
            continue

        if use_shared_cache:
            cache_path = _save_wari_iteration_cache(cache_output_dir, result)
            if cache_path is not None:
                cache_paths[iter_idx] = cache_path
        else:
            shared_results.append(result)

    available_results = len(cache_paths) if use_shared_cache else len(shared_results)
    print(f"[2/5] W-ARI calculation is completed round by round (can be summarized: {available_results}/{n_iter})")

    def iter_wari_results():
        if use_shared_cache:
            for iter_idx in sorted(cache_paths):
                yield _load_wari_iteration_cache(cache_paths[iter_idx])
        else:
            yield from shared_results

    memmap_dir = os.path.join(base_output_dir, "_memmap") if use_memmap else None
    if memmap_dir:
        os.makedirs(memmap_dir, exist_ok=True)
        total_cluster = np.memmap(
            os.path.join(memmap_dir, "wari_total_cluster.npy"),
            dtype=np.float64,
            mode='w+',
            shape=(n_cells, n_cells),
        )
        total_sample = np.memmap(
            os.path.join(memmap_dir, "wari_total_sample.npy"),
            dtype=np.float64,
            mode='w+',
            shape=(n_cells, n_cells),
        )
        total_cluster[:] = 0.0
        total_sample[:] = 0.0
        print(f"  W-ARI summary matrix uses memmap: {memmap_dir}")
    else:
        total_cluster = np.zeros((n_cells, n_cells), dtype=np.float64)
        total_sample = np.zeros((n_cells, n_cells), dtype=np.float64)

    wari_gene_scores = {}
    wari_results = []

    print(f"\n{'='*55}")
    print(f"  Summary W-ARI weighted results → {base_output_dir}")
    print(f"{'='*55}")

    for result in tqdm(
        iter_wari_results(),
        total=available_results,
        desc='W-ARI weighted summary',
        unit="iter",
    ):
        iter_idx = result['iter_idx']
        sampled_idx = result['sampled_idx']
        method_labels = result['method_labels']
        missing_methods = [
            name for name in METHOD_NAMES
            if name not in method_labels
        ]
        if missing_methods:
            failure_records.append({
                'iter_idx': iter_idx,
                'failed_methods': ';'.join(missing_methods),
                'error_detail': 'W-ARI requires all four base clustering labels',
            })
            continue

        all_labels = [method_labels[name] for name in METHOD_NAMES]
        _update_wari_total_counts(
            total_cluster,
            total_sample,
            sampled_idx,
            all_labels,
            global_weights,
        )

        scores = result.get('scores')
        if scores is not None:
            wari_gene_scores[iter_idx] = scores
        wari_results.append((iter_idx, sampled_idx, all_labels[0], None))

    failure_log_path = os.path.join(base_output_dir, "shared_iteration_failures.csv")
    failure_columns = ['iter_idx', 'failed_methods', 'error_detail']
    if failure_records:
        failure_df = pd.DataFrame(failure_records, columns=failure_columns)
        failure_df = failure_df.drop_duplicates()
    else:
        failure_df = pd.DataFrame(columns=failure_columns)
    failure_df.to_csv(failure_log_path, index=False, encoding='utf-8-sig')

    if not wari_results:
        raise RuntimeError(
            'W-ARI does not have complete iterations available for aggregation; please check shared_iteration_failures.csv.'
        )

    counts_file = os.path.join(base_output_dir, "wari_integration_summary.txt")
    with open(counts_file, 'w', encoding='utf-8') as file:
        file.write("Integration Mode: W-ARI\n")
        file.write("==================================\n")
        file.write(
            f"Leiden resolution: {best_leiden_resolution:.4f} "
            f"(full-data search, actual K={leiden_actual_k}, target K={n_clusters})\n"
        )
        file.write("Base methods: Kmeans, Seurat, ACM_WEP, Hierarchical\n")
        file.write(
            "Global W-ARI weights: "
            + ", ".join(
                f"{name}={weight:.6f}"
                for name, weight in zip(METHOD_NAMES, global_weights)
            )
            + "\n"
        )
        file.write("==================================\n")
        file.write(f"Total Successful Iterations: {len(wari_gene_scores)}\n")
    print(f"  W-ARI run summary saved: {counts_file}")

    aggregate_cell_sampling(
        n_cells,
        n_iter,
        wari_results,
        adata.obs_names,
        base_output_dir,
    )
    aggregate_gene_scores(
        wari_gene_scores,
        n_iter,
        base_output_dir,
        use_gpu=use_gpu,
    )

    consensus = build_wari_consensus_from_totals(
        total_cluster,
        total_sample,
        n_cells,
    )
    del total_cluster, total_sample

    save_consensus_matrix(consensus, adata.obs_names, base_output_dir)
    final_results = final_clustering_from_consensus(
        consensus,
        n_clusters=n_clusters,
        cell_names=adata.obs_names,
        output_dir=base_output_dir,
        true_labels_path=true_labels_path,
    )

    del consensus, wari_gene_scores, wari_results
    if not use_shared_cache:
        del shared_results
    gc.collect()

    if use_gpu and GPU_AVAILABLE:
        clear_gpu_memory()

    if memmap_dir and os.path.isdir(memmap_dir):
        shutil.rmtree(memmap_dir, ignore_errors=True)
        print('  W-ARI memmap temporary files cleaned')

    print(f"\n{'='*55}")
    print('  W-ARI weighted consistency summary completed')
    print(f"  Root directory: {base_output_dir}")
    print(f"{'='*55}")

    return {
        'w_ari': {
            'output_dir': base_output_dir,
            'weights': dict(zip(METHOD_NAMES, global_weights)),
            'final_results': final_results,
        }
    }

def build_consensus_matrix(co_cluster, co_sample, n_cells: int = None, use_gpu: bool = True):
    print('[3/5] Construct a consistency matrix...')

    if n_cells is None:
        n_cells = co_cluster.shape[0]

    if use_gpu and GPU_AVAILABLE and isinstance(co_cluster, torch.Tensor):
        co_cluster_f = co_cluster.float()
        co_sample_f  = co_sample.float()

        consensus_torch = torch.zeros_like(co_cluster_f)
        mask = co_sample_f != 0
        consensus_torch[mask] = co_cluster_f[mask] / co_sample_f[mask]
        consensus_torch.fill_diagonal_(1.0)

        consensus = to_cpu(consensus_torch)
        del co_cluster_f, co_sample_f, consensus_torch
        clear_gpu_memory()
    else:
        co_cluster_arr = np.asarray(co_cluster, dtype=np.float32)
        co_sample_arr  = np.asarray(co_sample,  dtype=np.float32)

        consensus = np.divide(
            co_cluster_arr,
            co_sample_arr,
            out=np.zeros_like(co_cluster_arr),
            where=co_sample_arr != 0
        )
        np.fill_diagonal(consensus, 1.0)

    print(f"  The consistency matrix is constructed ({n_cells} × {n_cells})")
    return consensus


def aggregate_gene_scores(
    gene_score_matrix: dict,
    n_iter: int,
    output_dir: str = "./output",
    use_gpu: bool = True
):
    """
    Summarize the gene scores of each iteration + calculate stability (optimize large-scale gene versions)

    Optimization points:
      1. CSV saves all genes (top is no longer truncated)
      2. Frequency statistics use numpy argpartition (O(n))
      3. Use pandas Series instead of Python dict
    """

    import os
    import numpy as np
    import pandas as pd

    print('[3/5] Summarize gene scores & calculate gene stability (optimized version)')
    os.makedirs(output_dir, exist_ok=True)

    # =========================
    # =========================
    all_genes = sorted(set(
        gene
        for scores in gene_score_matrix.values()
        for gene in scores.index
    ))

    n_genes = len(all_genes)
    n_successful = len(gene_score_matrix)

    print(f"      All round genes: {n_genes} (successful iteration: {n_successful} round)")

    # =========================
    # =========================
    # Use -1.0 as padding value because the normalized score range is [0, 1]
    # This can distinguish between "genes with a score of 0" and "genes that are not selected"
    scores_array = np.full((n_genes, n_successful), -1.0, dtype=np.float32)

    for col_idx, (iter_idx, scores) in enumerate(sorted(gene_score_matrix.items())):
        scores_array[:, col_idx] = scores.reindex(all_genes, fill_value=-1.0).values

    gene_matrix_df = pd.DataFrame(
        scores_array,
        index=all_genes,
        columns=[f'iter_{k}' for k in sorted(gene_score_matrix.keys())]
    )

    del scores_array

    matrix_path = os.path.join(output_dir, "gene_scores_by_iteration.parquet")
    # In order to keep subsequent calculations convenient, you can restore -1 back to 0 when saving, or save the original matrix directly
    # Here we save it directly, but pay attention to handling -1 when calculating statistical values.
    gene_matrix_df.to_parquet(matrix_path, index=True)
    print(f"      Gene × iteration matrix saved: {matrix_path}")

    # =========================
    # 2. Compute statistics (ignoring -1)
    # =========================
    temp_df = gene_matrix_df.replace(-1.0, np.nan)

    max_score = temp_df.max(axis=1).fillna(0)
    min_score = temp_df.min(axis=1).fillna(0)
    median_score = temp_df.median(axis=1).fillna(0)

    # =========================
    # 3. Cumulative score (treat -1 as 0, then sum)
    # =========================
    cumulative_score = gene_matrix_df.clip(lower=0).sum(axis=1)

    # =========================
    # 3.5 Detection frequency (for computing average)
    # =========================
    freq_detected = (gene_matrix_df >= 0).sum(axis=1)
    mean_score_when_detected = (cumulative_score / freq_detected.replace(0, np.nan)).fillna(0)

    # =========================
    # 4. Frequency (optimized)
    # =========================
    top_k = 300
    genes_index = gene_matrix_df.index
    freq_top300 = pd.Series(0, index=genes_index, dtype=np.int32)
    freq_nonzero = pd.Series(0, index=genes_index, dtype=np.int32)

    for col in gene_matrix_df.columns:
        values = gene_matrix_df[col].values
        genes = genes_index.values
        # Top300 frequency
        if len(values) > top_k:
            idx = np.argpartition(-values, top_k)[:top_k]
            top_genes = genes[idx]
        else:
            top_genes = genes
        freq_top300[top_genes] += 1
        # Detection frequency (score >= 0)
        nonzero_genes = genes[values >= 0]
        freq_nonzero[nonzero_genes] += 1

    # =========================
    # 5. Build stability result
    # =========================
    stability_result = pd.DataFrame({
        "gene": genes_index,
        "max_score": max_score,
        "min_score": min_score,
        "median_score": median_score,
        "mean_score_when_detected": mean_score_when_detected,
        "cumulative_score": cumulative_score,
        "freq_top300": freq_top300,
        "top300_rate(%)": (freq_top300 / n_iter * 100).round(2),
        "freq_detected": freq_nonzero,
        "detection_rate(%)": (freq_nonzero / n_iter * 100).round(2),
    }).sort_values(
        ["median_score", "mean_score_when_detected"],
        ascending=[False, False]
    )

    # =========================
    # =========================
    stability_path = os.path.join(output_dir, "gene_stability_result.csv")
    stability_result.to_csv(stability_path, index=False)

    print(f"      Gene stability results have been saved: {stability_path}")

    final_gene_score = median_score.sort_values(ascending=False)

    return gene_matrix_df, final_gene_score, stability_result


def save_consensus_matrix(consensus: np.ndarray, cell_names, output_dir: str = "./output"):
    os.makedirs(output_dir, exist_ok=True)

    n_cells = consensus.shape[0]
    n_pairs = n_cells * (n_cells - 1) // 2
    triu_path = os.path.join(output_dir, "consensus_matrix_triu.npy")
    meta_path = os.path.join(output_dir, "consensus_matrix_meta.csv")

    triu_values = np.lib.format.open_memmap(
        triu_path,
        mode="w+",
        dtype=consensus.dtype,
        shape=(n_pairs,),
    )
    offset = 0
    for row_idx in range(n_cells - 1):
        row_vals = consensus[row_idx, row_idx + 1:]
        next_offset = offset + row_vals.size
        triu_values[offset:next_offset] = row_vals
        offset = next_offset
    triu_values.flush()
    del triu_values

    pd.DataFrame([{
        "storage": "upper_triangle",
        "n_cells": n_cells,
        "n_pairs": n_pairs,
        "dtype": str(consensus.dtype),
    }]).to_csv(meta_path, index=False, encoding="utf-8-sig")

    full_path = os.path.join(output_dir, "consensus_matrix.npy")
    if os.path.exists(full_path):
        os.remove(full_path)

    print(f"      The upper triangle of the consistency matrix has been saved: {triu_path}")


def aggregate_cell_labels(n_cells, n_iter, results, cell_names, output_dir: str = "./output"):
    os.makedirs(output_dir, exist_ok=True)
    print('[5/5] Generate cell × iterative label matrix...')

    label_matrix_np = np.full((n_cells, n_iter), np.nan, dtype=np.float32)

    for iter_idx, sampled_idx, labels, _ in results:
        if sampled_idx is not None and labels is not None:
            label_matrix_np[sampled_idx, iter_idx] = labels

    cell_label_matrix = pd.DataFrame(
        label_matrix_np,
        index=cell_names,
        columns=[f'iter_{i}' for i in range(n_iter)]
    )

    csv_path = os.path.join(output_dir, "cell_labels_by_iteration.parquet")
    cell_label_matrix.to_parquet(csv_path, index=True)
    print(f"      Cell labeling matrix has been saved: {csv_path}")
    return cell_label_matrix


def aggregate_cell_sampling(n_cells, n_iter, results, cell_names, output_dir: str = "./output"):
    """
    Save whether each cell was selected in each iteration.
    Generate a matrix of cells × iterations (1=drawn, 0=not drawn).
    """
    os.makedirs(output_dir, exist_ok=True)
    print('[5/5] Generate cell × iterative sampling matrix...')

    sampling_matrix = np.zeros((n_cells, n_iter), dtype=np.int8)

    for iter_idx, sampled_idx, labels, _ in results:
        if sampled_idx is not None:
            sampling_matrix[sampled_idx, iter_idx] = 1

    sampling_df = pd.DataFrame(
        sampling_matrix,
        index=cell_names,
        columns=[f'iter_{i}' for i in range(n_iter)]
    )

    parquet_path = os.path.join(output_dir, "cell_sampling_by_iteration.parquet")
    sampling_df.to_parquet(parquet_path, index=True)
    print(f"      Cell sampling matrix has been saved: {parquet_path}")

    sample_counts = sampling_matrix.sum(axis=1)
    counts_df = pd.DataFrame({
        'cell_name': cell_names,
        'sampled_count': sample_counts,
        'sampled_rate(%)': (sample_counts / n_iter * 100).round(2)
    })
    counts_path = os.path.join(output_dir, "cell_sampling_counts.csv")
    counts_df.to_csv(counts_path, index=False)
    print(f"      The statistics of the number of cell extractions have been saved: {counts_path}")

    return sampling_df


def evaluate_clustering_metrics(final_result, true_df, algo_name, algo_dir, consensus=None):
    """
    Calculate clustering indicators (ARI, NMI, purity, consistency, etc.), count cluster sizes, and save the results to cluster_counts.txt.

    parameter:
        final_result: pd.DataFrame, contains 'cell_name' and 'final_label'
        true_df: pd.DataFrame, contains 'cell_name' and 'true_label' (evaluation is skipped if None)
        algo_name: str, algorithm name
        algo_dir: str, directory where results are saved
        consensus: np.ndarray, consensus matrix used to calculate internal metrics (optional)
    """
    import os
    import numpy as np
    import pandas as pd
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, fowlkes_mallows_score
    from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
    from scipy.optimize import linear_sum_assignment

    accuracy_info = []
    internal_metrics = []

    if true_df is not None and not true_df.empty:
        merged = pd.merge(final_result, true_df, on='cell_name', how='inner')

        if merged.empty:
            accuracy_info.append("Merge failed: no common cell_names between prediction and ground truth.")
        else:
            ari = adjusted_rand_score(merged['true_label'], merged['final_label'])
            nmi = normalized_mutual_info_score(merged['true_label'], merged['final_label'])
            fmi = fowlkes_mallows_score(merged['true_label'], merged['final_label'])
            accuracy_info.append(f"ARI: {ari:.4f}")
            accuracy_info.append(f"NMI: {nmi:.4f}")
            accuracy_info.append(f"FMI: {fmi:.4f}")
            accuracy_info.append("----------------------------------")

            # Partial match ACC
            true_labels_uniq = np.array(sorted(merged['true_label'].unique()))
            pred_labels_uniq = np.array(sorted(merged['final_label'].unique()))
            L = len(true_labels_uniq)
            K = len(pred_labels_uniq)

            true_idx = {v: i for i, v in enumerate(true_labels_uniq)}
            pred_idx = {v: i for i, v in enumerate(pred_labels_uniq)}
            C = np.zeros((L, K), dtype=int)
            for t, p in zip(merged['true_label'], merged['final_label']):
                C[true_idx[t], pred_idx[p]] += 1

            row, col = linear_sum_assignment(-C)
            matched = C[row, col].sum()
            total = len(merged)
            acc = matched / total

            accuracy_info.append(f"ACC (Partial Hungarian): {acc:.4f}")
            accuracy_info.append(f"Matched pairs: {len(row)} / min({L}, {K}) = {min(L, K)}")
            accuracy_info.append(f"Matched cells: {matched} / {total}")
            accuracy_info.append("Per-cluster matching:")
            for r, c in zip(row, col):
                true_name = true_labels_uniq[r]
                pred_name = pred_labels_uniq[c]
                hit = C[r, c]
                cluster_total = C[:, c].sum()
                accuracy_info.append(
                    f"Cluster C{pred_name} → {true_name}: "
                    f"{hit}/{cluster_total} (Hit={hit/cluster_total:.2%}, "
                    f"Contamination={1 - hit/cluster_total:.2%})"
                )
            if K < L:
                accuracy_info.append(f"  Unmatched true labels: {L - K} (all cells counted as miss)")
            elif K > L:
                accuracy_info.append(f"  Unmatched clusters: {K - L} (excess clusters not assigned)")

            # Intra-Cluster Consistency
            same_pair_probs = []
            true_label_lines = []
            unique_true_labels = merged['true_label'].unique()
            for tl in unique_true_labels:
                sub = merged[merged['true_label'] == tl]
                n_sub = len(sub)
                if n_sub < 2: continue
                pred_counts = sub['final_label'].value_counts()
                same_pred_pairs = sum([c * (c - 1) / 2 for c in pred_counts])
                total_pairs = n_sub * (n_sub - 1) / 2
                prob = same_pred_pairs / total_pairs
                same_pair_probs.append(prob)
                true_label_lines.append(f"True Label '{tl}': {prob:.4f} (cells={n_sub})")

            if same_pair_probs:
                avg_prob = np.mean(same_pair_probs)
                accuracy_info.append("----------------------------------")
                accuracy_info.append("")
                accuracy_info.append(f"Average Intra-Cluster Consistency: {avg_prob:.4f}")
                accuracy_info.append("----------------------------------")
                accuracy_info.extend(true_label_lines)

            # Purity / Contamination
            accuracy_info.append("")
            accuracy_info.append("Predicted Cluster Purity:")
            accuracy_info.append("----------------------------------")

            cluster_purity_lines = []
            weighted_purity_sum = 0
            for cl in sorted(merged['final_label'].unique()):
                sub = merged[merged['final_label'] == cl]
                total = len(sub)
                if total == 0: continue
                true_counts = sub['true_label'].value_counts()
                major_label = true_counts.idxmax()
                major_count = true_counts.max()
                purity = major_count / total
                weighted_purity_sum += purity * total
                cluster_purity_lines.append(f"Cluster C{cl}: Purity={purity:.4f}, Contamination={1-purity:.4f}, Major={major_label} ({major_count}/{total})")

            accuracy_info.append(f"Overall Purity: {weighted_purity_sum / len(merged):.4f}")
            accuracy_info.append(f"Overall Contamination: {1 - (weighted_purity_sum / len(merged)):.4f}")
            accuracy_info.append("----------------------------------")
            accuracy_info.extend(cluster_purity_lines)

    counts = final_result['final_label'].value_counts().sort_index()

    counts_path = os.path.join(algo_dir, "cluster_counts.txt")
    with open(counts_path, 'w', encoding='utf-8') as f:
        f.write(f"Algorithm: {algo_name}\n")

        if internal_metrics:
            f.write("\nInternal Clustering Metrics:\n")
            f.write("----------------------------------\n")
            for line in internal_metrics:
                f.write(line + "\n")
            f.write("\n")

        if accuracy_info:
            f.write("Comparison with Ground Truth:\n")
            for line in accuracy_info:
                f.write(line + "\n")
            f.write("========================\n\n")

        f.write("Predicted Cluster Counts Summary:\n")
        f.write("========================\n")
        for cluster_id, count in counts.items():
            f.write(f"Cluster C{cluster_id}: {count} cells\n")
        f.write("========================\n")
        f.write(f"Total: {counts.sum()} cells\n")

    return accuracy_info


def final_clustering_from_consensus(consensus: np.ndarray, n_clusters: int, cell_names, output_dir: str = "./output", true_labels_path: str = None):
    """
    Final clustering is performed based on the consistency matrix.
    """
    import os
    import pandas as pd
    import numpy as np
    from ensemble_clustering.jullei._01_08_01_08_ACM_WEP_matrix import cluster_demo_pro as acm_final

    true_df = None
    if true_labels_path and os.path.exists(true_labels_path):
        try:
            true_df = pd.read_csv(true_labels_path)
            first_col = true_df.columns[0]
            if first_col != 'cell_name':
                true_df.rename(columns={first_col: 'cell_name'}, inplace=True)
            if 'celltype' in true_df.columns:
                true_df.rename(columns={'celltype': 'true_label'}, inplace=True)
            print(f"  Successfully loaded original tag: {true_labels_path}")
        except Exception as e:
            print(f"  Failed to load original tag: {e}")

    algorithms = {
        'ACMWEP': lambda cons, k: acm_final(cons, k, input_type='consensus')[0]
    }

    results = {}
    for algo_name, algo_func in algorithms.items():
        print(f"[4/5] Run the final clustering algorithm: {algo_name} (K={n_clusters})")
        algo_dir = os.path.join(output_dir, algo_name)
        os.makedirs(algo_dir, exist_ok=True)

        try:
            final_labels = algo_func(consensus, n_clusters)
            final_result = pd.DataFrame({'cell_name': cell_names, 'final_label': final_labels})

            final_result.to_csv(os.path.join(algo_dir, "final_labels.csv"), index=False, header=True)

            evaluate_clustering_metrics(final_result, true_df, algo_name, algo_dir)

            plot_consensus_heatmap_raw(consensus, n_clusters, algo_dir, filename_prefix=algo_name)
            plot_consensus_heatmap(consensus, n_clusters, final_labels, cell_names, algo_dir, filename_prefix=algo_name)

            results[algo_name] = final_result
            print(f"  {algo_name} Results saved to: {algo_dir}")
        except Exception as e:
            print(f"  {algo_name} Clustering failed: {e}")
            import traceback
            traceback.print_exc()

    return results


def load_consensus_triu_values(consensus_path, row_block=512):
    """
    Load upper triangular values (sorted) from a consistency matrix file.
    Supports consensus_matrix_triu.npy (1D) and consensus_matrix.npy (full matrix).
    """
    if os.path.basename(consensus_path) == "consensus_matrix_triu.npy":
        triu_values = np.load(consensus_path, mmap_mode="r")
        sorted_vals = np.sort(triu_values)
        del triu_values
        return sorted_vals

    consensus = np.load(consensus_path, mmap_mode="r")
    n_cells = consensus.shape[0]
    n_pairs = n_cells * (n_cells - 1) // 2
    triu_values = np.empty(n_pairs, dtype=consensus.dtype)
    offset = 0
    for block_start in range(0, max(n_cells - 1, 0), row_block):
        block_end = min(n_cells - 1, block_start + row_block)
        for row_idx in range(block_start, block_end):
            row_vals = consensus[row_idx, row_idx + 1:]
            next_offset = offset + row_vals.size
            triu_values[offset:next_offset] = row_vals
            offset = next_offset
    del consensus
    return np.sort(triu_values)





def plot_consensus_heatmap_raw(consensus, k, output_dir, filename_prefix="consensus"):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    import os

    # heatmap_cmap = LinearSegmentedColormap.from_list('white_red', ['white', '#A50026'])
    from matplotlib.colors import LinearSegmentedColormap
    heatmap_cmap = LinearSegmentedColormap.from_list(
        "white_yellow_red", ["white", "yellow", "red"])
    plt.figure(figsize=(10, 8))
    plt.imshow(consensus, cmap=heatmap_cmap, vmin=0, vmax=1, aspect='auto', interpolation='nearest')
    cbar = plt.colorbar()
    cbar.set_label('Consensus', fontsize=12, fontweight='bold')
    plt.xlabel('Patient Index', fontsize=12)
    plt.ylabel('Patient Index', fontsize=12)
    plt.title(f'Consensus Matrix (K={k})', fontweight='bold', fontsize=14)
    plt.xticks([])
    plt.yticks([])
    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"{filename_prefix}_raw_heatmap_k{k}.png")
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f" Original consistency matrix heat map completed: {save_path}")

def plot_consensus_heatmap(
    consensus, k, labels, cell_names, output_dir,
    filename_prefix="consensus"
):
    """
     The most stable heatmap sorting version (recommended for production use)

     Mean similarity sorting within cluster (stable + fast)
     block clear structure
     10k cells available
     Does not rely on PCA/MST/hierarchical
    """

    import os
    import numpy as np
    import matplotlib.pyplot as plt
    import colorsys
    import matplotlib.patches as mpatches
    from matplotlib.colors import LinearSegmentedColormap

    # =========================
    # =========================
    n_total = consensus.shape[0]
    labels_sample = labels
    consensus_sample = consensus

    # =========================
    # 2. Sorting within the cluster (introducing hierarchical clustering optimization to improve block clarity)
    # =========================
    from scipy.cluster.hierarchy import linkage, leaves_list
    from scipy.spatial.distance import squareform

    unique_labels = sorted(np.unique(labels_sample))
    order_idx = []

    for lb in unique_labels:
        lb_idx = np.where(labels_sample == lb)[0]

        if len(lb_idx) <= 2:
            order_idx.extend(lb_idx.tolist())
            continue

        sub_matrix = consensus_sample[np.ix_(lb_idx, lb_idx)]

        # Convert similarity to distance (1 - similarity), ensuring symmetry
        dist_mat = 1.0 - (sub_matrix + sub_matrix.T) / 2.0
        np.fill_diagonal(dist_mat, 0)

        try:
            # This allows "more similar" cells to be closer together in the heatmap, resulting in a more compact clump
            Z = linkage(squareform(dist_mat, checks=False), method='average')
            sub_order = leaves_list(Z)
            order_idx.extend(lb_idx[sub_order].tolist())
        except:
            # Fault tolerance: if hierarchical clustering fails, fall back to average similarity ranking
            score = sub_matrix.mean(axis=1)
            order_idx.extend(lb_idx[np.argsort(-score)].tolist())

    order_idx = np.array(order_idx)

    consensus_ordered = consensus_sample[np.ix_(order_idx, order_idx)]
    labels_ordered = labels_sample[order_idx]
    n_cells = len(labels_ordered)

    # =========================
    # =========================
    unique_valid = sorted(np.unique(labels[labels != -1]))
    n_clusters = len(unique_valid)

    def generate_rainbow_colors(n, s=0.7, v=0.85):
        return [colorsys.hsv_to_rgb(i / n, s, v) for i in range(n)]

    rainbow_colors = generate_rainbow_colors(n_clusters)
    label_to_color = {
        label: rainbow_colors[i] for i, label in enumerate(unique_valid)
    }
    label_to_color[-1] = (0.827, 0.827, 0.827)

    cluster_colors_rgba = np.array([
        list(label_to_color[lb]) + [1.0]
        for lb in labels_ordered
    ])

    # =========================
    # =========================
    fig = plt.figure(figsize=(14, 13))

    ax_left = fig.add_axes([0.08, 0.12, 0.02, 0.75])
    ax_left.imshow(cluster_colors_rgba.reshape(-1, 1, 4), aspect='auto', interpolation='nearest')
    ax_left.set_xticks([])
    ax_left.set_yticks([])

    ax_bottom = fig.add_axes([0.10, 0.10, 0.75, 0.02])
    ax_bottom.imshow(cluster_colors_rgba.reshape(1, -1, 4), aspect='auto', interpolation='nearest')
    ax_bottom.set_xticks([])
    ax_bottom.set_yticks([])

    ax_heatmap = fig.add_axes([0.1, 0.12, 0.75, 0.75])

    cmap = LinearSegmentedColormap.from_list(
        "white_yellow_red", ["white", "yellow", "red"]
    )

    im = ax_heatmap.imshow(
        consensus_ordered,
        cmap=cmap,
        vmin=0,
        vmax=1,
        aspect='auto',
        interpolation='nearest',
        rasterized=True
    )

    ax_heatmap.set_xticks([])
    ax_heatmap.set_yticks([])
    ax_heatmap.set_xlabel("Cell Index")
    ax_heatmap.set_ylabel("Cell Index")

    # =========================
    # 5. colorbar
    # =========================
    ax_cbar = fig.add_axes([0.87, 0.12, 0.015, 0.20])
    cbar = plt.colorbar(im, cax=ax_cbar)
    cbar.set_label("Consensus")

    # =========================
    # 6. legend
    # =========================
    ax_legend = fig.add_axes([0.87, 0.8, 0.12, 0.10])
    ax_legend.axis("off")

    handles = [
        mpatches.Patch(
            facecolor=label_to_color[lb],
            edgecolor="black",
            label=f"C{lb}"
        )
        for lb in unique_valid
    ]

    if -1 in labels:
        handles.append(
            mpatches.Patch(
                facecolor=label_to_color[-1],
                edgecolor="black",
                label="Noise"
            )
        )

    ax_legend.legend(
        handles=handles,
        loc="center",
        fontsize=9,
        frameon=True,
        ncol=2,
        title="Clusters"
    )

    # =========================
    # 7. save
    # =========================
    # fig.suptitle(
    #     f"Consensus Matrix Heatmap  (K={k})",
    #     fontsize=20,
    #     fontweight='bold',
    # )
    # plt.tight_layout()

    save_path = os.path.join(
        output_dir,
        f"{filename_prefix}_heatmap_k{k}.png"
    )

    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()

    print(f" Heatmap completed ({n_cells} cells): {save_path}")
