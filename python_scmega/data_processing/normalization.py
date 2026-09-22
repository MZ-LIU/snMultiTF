"""
Normalization Module - R-Compatible Implementation with Muon

Use muon/scanpy to handle the standardization of multi-omics data, as close as possible to the R version Seurat/Signac behavior
"""

import numpy as np
import anndata as ad
from typing import Optional, Dict, Any, List, Union
import warnings
from scipy import sparse
from scipy.sparse.linalg import svds
import scanpy as sc

try:
    import muon as mu
    MUON_AVAILABLE = True
except ImportError:
    mu = None
    MUON_AVAILABLE = False
    warnings.warn("muon not available. Please install: pip install muon")


# ============================================================================
# ============================================================================

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


# ============================================================================
# ============================================================================

def normalize_rna(adata: ad.AnnData,
                  target_sum: Optional[float] = None,
                  log_transform: bool = True,
                  scale_data: bool = True,
                  n_top_genes: Optional[int] = None,
                  run_pca: bool = True,
                  n_pcs: Optional[int] = None,
                  run_umap: bool = True,
                  umap_dims: Optional[int] = None,
                  random_state: int = 42,
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
        hvg_mask = target.var['highly_variable'].values
        target_hvg = target[:, hvg_mask].copy()

        sc.pp.scale(target_hvg, max_value=10)
        
        # Stored in .uns (to avoid AnnData’s layers dimension limitation)
        target.uns['scaled_hvg'] = {
            'data': target_hvg.X.copy(),
            'genes': target.var_names[hvg_mask].tolist()
        }
        
        print(f"   Scaled {hvg_mask.sum()} HVGs (stored in .uns['scaled_hvg'])")

    if run_pca:
        if not scale_data or 'scaled_hvg' not in target.uns:
            raise ValueError("PCA requires scaled data. Set scale_data=True.")
        
        target_for_pca = target[:, target.var['highly_variable']].copy()
        target_for_pca.X = target.uns['scaled_hvg']['data']
        
        sc.tl.pca(
            target_for_pca,
            n_comps=n_pcs,
            svd_solver='arpack',
            random_state=random_state,
        )
        target.obsm['X_pca'] = target_for_pca.obsm['X_pca']
        target.varm['PCs'] = np.zeros((target.n_vars, n_pcs), dtype=np.float32)
        target.varm['PCs'][target.var['highly_variable'].values] = target_for_pca.varm['PCs']
        target.uns['pca'] = target_for_pca.uns['pca']
        print(f"   PCA completed: {n_pcs} components")

    if run_umap and run_pca:
        sc.pp.neighbors(target, n_pcs=umap_dims, use_rep='X_pca', n_neighbors=30, metric='euclidean')
        sc.tl.umap(target, min_dist=0.3, spread=1.0, random_state=random_state)

    return target


# ============================================================================
# ============================================================================

def _find_top_features_signac_like(adata: ad.AnnData,
                                   n_top_features: Optional[int] = None,
                                   min_cutoff: str = 'q0'):
    """
    Closer to Signac::FindTopFeatures:
    - Use the "accessibility frequency" (number of cells counted >0) of the original counts (or layers['counts']) as the filtering basis
    - Supports quantity (n_top_features) or quantile threshold (min_cutoff='q0', 'q0.9', etc.)
    """
    X = adata.layers['counts'] if 'counts' in adata.layers else adata.X
    if sparse.issparse(X):
        freq = np.array((X > 0).sum(axis=0)).ravel()
    else:
        freq = (np.asarray(X) > 0).sum(axis=0)

    if n_top_features is not None:
        idx = np.argsort(freq)[::-1][:n_top_features]
    else:
        if isinstance(min_cutoff, str) and min_cutoff.startswith('q'):
            try:
                q = float(min_cutoff[1:]) if len(min_cutoff) > 1 else 0.0
            except Exception:
                q = 0.0
            thr = np.quantile(freq, q)
            idx = np.where(freq >= thr)[0]
        else:
            idx = np.arange(len(freq))

    if 'top_feature' not in adata.var.columns:
        adata.var['top_feature'] = False
    adata.var['top_feature'] = False
    adata.var.iloc[idx, adata.var.columns.get_loc('top_feature')] = True
    adata.var['accessibility_score'] = freq.astype(np.int64)



def _run_lsi_with_seed(adata: ad.AnnData,
                       n_comps: int = 50,
                       random_state: int = 42,
                       scale_embeddings: bool = False) -> None:
    """
    Reproduce the main logic of Signac RunSVD: store the original U matrix as cell embeddings without additional scaling.
    Consistent with the R scMEGA tutorial (Signac does not scale LSI embeddings by default).

    The result is written as:
    - adata.obsm['X_lsi']
    - adata.uns['lsi']
    - adata.varm['LSI']
    """
    X = adata.X
    if X is None:
        raise ValueError("LSI requires a non-empty matrix in adata.X.")

    max_comps = min(X.shape) - 1
    if max_comps < 1:
        raise ValueError(
            f"LSI requires at least 2 cells and 2 features, got shape={X.shape}."
        )
    k = min(int(n_comps), max_comps)

    # Both writing methods explicitly fix random initialization to avoid muon's internal default random starting point.
    try:
        u, s, vt = svds(
            X,
            k=k,
            which='LM',
            solver='arpack',
            rng=np.random.default_rng(random_state),
        )
    except TypeError:
        u, s, vt = svds(
            X,
            k=k,
            which='LM',
            solver='arpack',
            random_state=random_state,
        )

    # svds does not guarantee the return order, and is arranged from largest to smallest singular values.
    order = np.argsort(s)[::-1]
    u = u[:, order]
    s = s[order]
    vt = vt[order, :]

    # Fix the sign direction of each component to avoid equivalent overall positive and negative flips
    max_abs_rows = np.argmax(np.abs(u), axis=0)
    signs = np.sign(u[max_abs_rows, np.arange(u.shape[1])])
    signs[signs == 0] = 1.0
    u *= signs
    vt *= signs[:, None]

    if scale_embeddings:
        means = u.mean(axis=0)
        stds = u.std(axis=0)
        stds[stds == 0] = 1.0
        u = (u - means) / stds

    adata.obsm['X_lsi'] = u
    adata.uns['lsi'] = {
        'stdev': s / np.sqrt(X.shape[0] - 1),
        'params': {
            'n_comps': k,
            'scale_embeddings': scale_embeddings,
            'random_state': int(random_state),
            'solver': 'arpack',
        },
    }
    adata.varm['LSI'] = vt.T


def normalize_atac(adata: ad.AnnData,
                   run_tfidf: bool = True,
                   n_top_features: Optional[int] = None,
                   run_lsi: bool = True,
                   n_components: Optional[int] = None,
                   run_umap: bool = True,
                   umap_dims: Optional[int] = None,
                   random_state: int = 42,
                   inplace: bool = False) -> ad.AnnData:
    """
    ATAC normalization using muon - closer to Signac behavior:
    - RunTFIDF (does not do additional 1e4 scaling)
    - FindTopFeatures (accessibility frequency based on original counts)
    - RunSVD/LSI (on Top features)
    - UMAP (in LSI space, use component 2:30; neighbors use cosine, n_neighbors=30; UMAP parameters are more robust)
    """
    if not MUON_AVAILABLE:
        raise ImportError("muon package is required for ATAC normalization. Install with: pip install muon")

    # Step one: Make sure layers['counts'] exists (only on first call)
    if 'counts' not in adata.layers:
        if not _is_raw_counts(adata.X):
            raise ValueError(
                "First-time normalization requires raw counts in .X for ATAC, "
                "but detected non-integer values. Please ensure input is raw count matrix."
            )
        adata.layers['counts'] = adata.X.copy()

    target = adata if inplace else adata.copy()

    target.X = target.layers['counts'].copy()

    if n_components is None:
        n_components = 50  # Default number of SVD components
    if umap_dims is None:
        umap_dims = 30

    if run_tfidf:
        mu.atac.pp.tfidf(target)

    # Step 2: FindTopFeatures（Signac-like）
    if n_top_features is not None:
        _find_top_features_signac_like(target, n_top_features=n_top_features)
    else:
        # The default q0 is equivalent to full retention; if you want quantile filtering, you can change it to q0.9, etc.
        _find_top_features_signac_like(target, n_top_features=None, min_cutoff='q0')

    if run_lsi:
        if 'top_feature' in target.var.columns:
            lsi_target = target[:, target.var['top_feature']].copy()
        else:
            lsi_target = target

        _run_lsi_with_seed(
            lsi_target,
            n_comps=n_components,
            random_state=random_state,
        )

        target.obsm['X_lsi'] = lsi_target.obsm['X_lsi']
        target.uns['lsi'] = lsi_target.uns['lsi']

    if run_umap and run_lsi:
        lsi_for_umap = target.obsm['X_lsi'][:, 1:umap_dims]
        target.obsm['X_lsi_for_umap'] = lsi_for_umap

        sc.pp.neighbors(
            target,
            use_rep='X_lsi_for_umap',
            n_neighbors=30,
            metric='cosine',
            random_state=random_state,
        )  # 30 is a bit high, consider 15

        sc.tl.umap(
            target,
            min_dist=0.3,
            spread=1.0,
            random_state=random_state,
        )

    return target


# ============================================================================
# ============================================================================

def normalize_multiome_muon(mdata: 'mu.MuData',
                            rna_params: Optional[Dict[str, Any]] = None,
                            atac_params: Optional[Dict[str, Any]] = None,
                            inplace: bool = False) -> 'mu.MuData':
    """
    Use muon to complete the standardization process of multi-omics data
    """
    if not MUON_AVAILABLE:
        raise ImportError("muon package is required. Install with: pip install muon")

    target = mdata if inplace else mdata.copy()

    if rna_params is None:
        rna_params = {}
    if atac_params is None:
        atac_params = {}

    if 'rna' in target.mod:
        normalize_rna(target['rna'], inplace=True, **rna_params)

    if 'atac' in target.mod:
        normalize_atac(target['atac'], inplace=True, **atac_params)

    target.uns['normalized'] = True
    target.uns['normalization_params'] = {
        'rna_params': rna_params,
        'atac_params': atac_params
    }
    return target



# ============================================================================
# ============================================================================

def get_normalization_summary_muon(mdata: 'mu.MuData') -> str:
    """
    Gets the normalized status summary of a MuData object
    
    Args:
        mdata: MuData object
        
    Returns:
        Formatted summary string
    """
    summary_lines = []
    summary_lines.append("\n" + "="*60)
    summary_lines.append("Normalization Summary")
    summary_lines.append("="*60)
    
    if 'rna' in mdata.mod:
        rna = mdata['rna']
        summary_lines.append("\nRNA Modality:")
        summary_lines.append(f"  Cells: {rna.n_obs:,}")
        summary_lines.append(f"  Genes: {rna.n_vars:,}")
        
        if 'counts' in rna.layers:
            summary_lines.append(f"   Raw counts preserved in layers['counts']")
        
        if 'highly_variable' in rna.var.columns:
            n_hvg = rna.var['highly_variable'].sum()
            summary_lines.append(f"   Highly variable genes: {n_hvg:,}")
        
        if 'scaled' in rna.layers:
            summary_lines.append(f"   Scaled data in layers['scaled']")
        
        if 'X_pca' in rna.obsm:
            n_pcs = rna.obsm['X_pca'].shape[1]
            summary_lines.append(f"   PCA: {n_pcs} components")
        
        if 'X_umap' in rna.obsm:
            summary_lines.append(f"  UMAP computed")
    
    if 'atac' in mdata.mod:
        atac = mdata['atac']
        summary_lines.append("\nATAC Modality:")
        summary_lines.append(f"  Cells: {atac.n_obs:,}")
        summary_lines.append(f"  Peaks: {atac.n_vars:,}")
        
        if 'counts' in atac.layers:
            summary_lines.append(f"   Raw counts preserved in layers['counts']")
        
        if 'top_feature' in atac.var.columns:
            n_top = atac.var['top_feature'].sum()
            summary_lines.append(f"   Top accessible features: {n_top:,}")
        
        if 'X_lsi' in atac.obsm:
            n_comps = atac.obsm['X_lsi'].shape[1]
            summary_lines.append(f"   LSI: {n_comps} components")
        
        if 'X_umap' in atac.obsm:
            summary_lines.append(f"   UMAP computed")
    
    summary_lines.append("="*60 + "\n")
    
    return "\n".join(summary_lines)
