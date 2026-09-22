"""
Co-embedding Module

Implements multi-modal data co-embedding algorithms equivalent to Seurat's
CoembedData and CreatePairedObject functions for integrating RNA and ATAC data.

Key functions:
- coembed_data: Create joint embedding of RNA and ATAC data


References:
- Original R code: CoembedData, CreatePairedObject functions
- Seurat multimodal integration workflow
- Stuart et al. (2019) Comprehensive Integration of Single-Cell Data
"""

import numpy as np
import anndata as ad
from typing import Optional,  List
import scanpy as sc


# Try to import MuData, set to None if not available
try:
    import muon as mu
    MUON_AVAILABLE = True
except ImportError:
    mu = None
    MUON_AVAILABLE = False


def coembed_data(rna_data: ad.AnnData,
                atac_data: ad.AnnData,
                gene_activity: Optional[ad.AnnData] = None,
                rna_reduction: str = "pca",
                atac_reduction: str = "lsi",
                rna_dims: Optional[List[int]] = None,
                atac_dims: Optional[List[int]] = None,
                p_value_threshold: float = 0.05,
                verbose: bool = True) -> 'mu.MuData':
    """
    Create joint co-embedding of RNA and ATAC data.
    
    This function replicates R's MOJITOO integration method using CCA:
    ```R
    pbmc <- mojitoo(
        object = pbmc,
        reduction.list = list("RNA_PCA", "lsi"),
        dims.list = list(1:50, 2:50),
        reduction.name = 'MOJITOO',
        assay = "RNA"
    )
    ```
    
    The integration uses Canonical Correlation Analysis (CCA):
    1. Extract specified dimensions from RNA and ATAC reductions
    2. Perform CCA to find canonical variables
    3. Calculate Pearson correlations and p-values for each pair
    4. Apply BH correction for multiple testing
    5. Select significant dimensions (p < threshold)
    6. Fuse canonical variables by summation
    
    Args:
        rna_data: RNA AnnData object with normalized expression data
        atac_data: ATAC AnnData object with chromatin accessibility data
        gene_activity: **OPTIONAL** - Gene activity matrix for alternative integration.
                      For MOJITOO-style integration (default), this is NOT required.
        rna_reduction: Reduction to use for RNA data (default: "pca")
                      Corresponds to R's reduction.list[1]
        atac_reduction: Reduction to use for ATAC data (default: "lsi")
                       Corresponds to R's reduction.list[2]
        rna_dims: Dimensions to use from RNA reduction (default: 0-49, equivalent to R's 1:50)
        atac_dims: Dimensions to use from ATAC reduction (default: 1-49, equivalent to R's 2:50)
                  Note: By default excludes first LSI dimension (technical artifact)
        p_value_threshold: Threshold for selecting significant CCA dimensions (default: 0.05)
                          Corresponds to R's significance testing threshold
        verbose: Whether to print progress messages (default: True)
        
    Returns:
        MuData object with joint embedding stored in obsm['X_integrated'] (global level)
        
    Example:
        ```python
        # MOJITOO-style CCA integration (R-compatible)
        integrated_data = coembed_data(
            rna_data=multiome_data.rna,
            atac_data=multiome_data.atac,
            rna_reduction="pca",        # RNA_PCA in R
            atac_reduction="lsi",       # lsi in R  
            rna_dims=list(range(50)),   # 1:50 in R
            atac_dims=list(range(1,50)),# 2:50 in R (exclude 1st LSI dimension)
            p_value_threshold=0.05      # Significance threshold
        )
        ```
        
    Note:
        This function now implements true MOJITOO CCA integration,
        which does NOT require gene_activity. This matches the R version
        workflow used in scMEGA's network9.25.Rmd.
        
        R equivalent workflow:
        ```R
        # RNA preprocessing
        pbmc <- pbmc %>% RunPCA(npcs=50, reduction.name="RNA_PCA")
        
        # ATAC preprocessing  
        pbmc <- pbmc %>% RunTFIDF() %>% RunSVD()
        
        # MOJITOO CCA integration
        pbmc <- mojitoo(
            object = pbmc,
            reduction.list = list("RNA_PCA", "lsi"),
            dims.list = list(1:50, 2:50)
        )
        # Internally: CCA → correlation testing → BH correction → fusion
        ```
    """
    if verbose:
        print("Creating joint co-embedding of RNA and ATAC data...")
        print("  Method: MOJITOO-style integration (R-compatible)")
    
    # =====================================================
    # R-COMPATIBLE: MOJITOO-style integration (gene_activity is OPTIONAL)
    # =====================================================
    # In R version: mojitoo(object, reduction.list, dims.list, ...)
    # gene_activity is NOT required for MOJITOO integration
    
    if gene_activity is not None:
        if verbose:
            print(f"  ℹ Gene activity matrix provided: {gene_activity.shape}")
            print(f"    Note: For MOJITOO-style integration, gene_activity is optional")
            print(f"    Will proceed with standard MOJITOO integration")
    else:
        if verbose:
            print("   Using MOJITOO-style integration (no gene_activity required)")
    
    # Set default dimensions (R-compatible)
    # R version: dims.list = list(1:50, 2:50)
    # Python 0-indexed: RNA uses 0-49, ATAC uses 1-49 (skip 1st LSI dimension)
    if rna_dims is None:
        rna_dims = list(range(50))  # 0-49, equivalent to R's 1:50
        if verbose:
            print(f"  Using default RNA dimensions: 0-49 (50 dimensions)")
    
    if atac_dims is None:
        atac_dims = list(range(1, 50))  # 1-49, equivalent to R's 2:50 (exclude 1st LSI)
        if verbose:
            print(f"  Using default ATAC dimensions: 1-49 (49 dimensions, excluding 1st LSI)")
    else:
        # User provided atac_dims
        if verbose:
            print(f"  Using custom ATAC dimensions: {atac_dims[0]}-{atac_dims[-1]} ({len(atac_dims)} dimensions)")
    
    # Get embeddings
    rna_embedding_key = f'X_{rna_reduction}'
    atac_embedding_key = f'X_{atac_reduction}'
    
    if rna_embedding_key not in rna_data.obsm:
        if verbose:
            print(f"Computing {rna_reduction} for RNA data...")
        sc.tl.pca(rna_data, n_comps=max(rna_dims) + 1, random_state=42)
        rna_embedding_key = 'X_pca'
    
    if atac_embedding_key not in atac_data.obsm:
        if verbose:
            print(f"Computing {atac_reduction} for ATAC data...")
        sc.tl.pca(atac_data, n_comps=max(atac_dims) + 1, random_state=42)
        atac_embedding_key = 'X_pca'
    
    # Extract embeddings with specified dimensions (R-compatible)
    # R version extracts: reduction.list[[1]][, dims.list[[1]]] and reduction.list[[2]][, dims.list[[2]]]
    rna_embedding = rna_data.obsm[rna_embedding_key][:, rna_dims]
    atac_embedding = atac_data.obsm[atac_embedding_key][:, atac_dims]
    
    if verbose:
        print(f"  Extracted RNA embedding ({rna_reduction}): {rna_embedding.shape}")
        print(f"  Extracted ATAC embedding ({atac_reduction}): {atac_embedding.shape}")
        print(f"  Total dimensions for integration: {rna_embedding.shape[1] + atac_embedding.shape[1]}")
    
    # Find shared cells while preserving RNA cell order for reproducible outputs.
    shared_cells = rna_data.obs_names.intersection(atac_data.obs_names).tolist()
    
    if len(shared_cells) == 0:
        raise ValueError(
            "No shared cells found between RNA and ATAC data. "
            "scMEGA requires paired multiome data with shared cell identifiers."
        )
    
    if verbose:
        print(f"Found {len(shared_cells)} shared cells")
    
    # Get embeddings for shared cells
    rna_shared_indices = [rna_data.obs.index.get_loc(cell) for cell in shared_cells]
    atac_shared_indices = [atac_data.obs.index.get_loc(cell) for cell in shared_cells]
    
    rna_shared_embedding = rna_embedding[rna_shared_indices]
    atac_shared_embedding = atac_embedding[atac_shared_indices]
    
    # Create joint embedding using MOJITOO-style CCA integration
    # R version: MOJITOO uses CCA to combine the two reduction embeddings
    if verbose:
        print("  Creating MOJITOO-style joint embedding with CCA...")
        print(f"    Combining {rna_shared_embedding.shape[1]} RNA dims + {atac_shared_embedding.shape[1]} ATAC dims")
        print(f"    P-value threshold: {p_value_threshold}")
    
    joint_embedding = _create_joint_embedding(
        rna_shared_embedding, atac_shared_embedding,
        p_value_threshold=p_value_threshold,
        verbose=verbose
    )
    
    # Create MuData object
    rna_shared = rna_data[shared_cells].copy()
    atac_shared = atac_data[shared_cells].copy()
    
    # Create MuData and set joint embedding using new data structure
    multiome = mu.MuData({'rna': rna_shared, 'atac': atac_shared})
    
    # Store joint embedding in centralized reductions
    multiome.obsm['X_integrated'] = joint_embedding
    
    # Also store as X_joint for backward compatibility
    # rna_shared.obsm['X_joint'] = joint_embedding
    # atac_shared.obsm['X_joint'] = joint_embedding
    multiome.uns['shared_cells'] = shared_cells
    
    # Store co-embedding parameters (R-compatible)
    if not hasattr(multiome, 'uns') or multiome.uns is None:
        multiome.uns = {}
    multiome.uns['coembedding'] = {
        'method': 'MOJITOO-CCA',  # R-compatible CCA integration method
        'reduction_list': [rna_reduction, atac_reduction],  # Equivalent to R's reduction.list
        'dims_list': {'rna': rna_dims, 'atac': atac_dims},
        'rna_reduction': rna_reduction,
        'atac_reduction': atac_reduction,
        'rna_dims': rna_dims,
        'atac_dims': atac_dims,
        'p_value_threshold': p_value_threshold,  # CCA significance threshold
        'reduction_name': 'joint'  # Equivalent to R's reduction.name parameter
    }
    if verbose:
        print(f"\n MOJITOO-CCA co-embedding completed!")
        print(f"  Integrated cells: {joint_embedding.shape[0]}")
        print(f"  Joint embedding dimensions: {joint_embedding.shape[1]} (after CCA and significance filtering)")
        print(f"  Method: CCA with p-value<{p_value_threshold}")
        print(f"  Stored in: mdata.obsm['X_integrated'] (global level)") 
    
    return multiome



def _create_joint_embedding(rna_embedding: np.ndarray,
                           atac_embedding: np.ndarray,
                           p_value_threshold: float = 0.05,
                           verbose: bool = True) -> np.ndarray:
    """
    Create joint embedding using CCA (Canonical Correlation Analysis).
    
    This implements the MOJITOO integration method:
    1. Perform CCA on RNA and ATAC embeddings
    2. Calculate Pearson correlation coefficients for each canonical variable pair
    3. Perform significance testing and BH correction
    4. Select significant dimensions (p-value < threshold)
    5. Fuse canonical variables by summation
    
    Args:
        rna_embedding: RNA embedding matrix (n_cells × n_rna_dims)
        atac_embedding: ATAC embedding matrix (n_cells × n_atac_dims)
        p_value_threshold: Threshold for selecting significant dimensions (default: 0.05)
        verbose: Whether to print progress information
        
    Returns:
        joint_embedding: Fused embedding matrix (n_cells × n_significant_dims)
    """
    from sklearn.cross_decomposition import CCA
    from scipy.stats import pearsonr
    from statsmodels.stats.multitest import multipletests
    
    if verbose:
        print("  Performing CCA integration (MOJITOO method)...")
        print(f"    Input: RNA {rna_embedding.shape}, ATAC {atac_embedding.shape}")
    
    # Step 1: Determine number of CCA components
    # Number of canonical variables = min(n_features_rna, n_features_atac)
    # Note: CCA with scale=True (default) handles standardization internally,
    # consistent with the R MOJITOO workflow (is_reduction_center=False, is_reduction_scale=False).
    n_components = min(rna_embedding.shape[1], atac_embedding.shape[1])

    if verbose:
        print(f"    CCA components: {n_components} (min of {rna_embedding.shape[1]} and {atac_embedding.shape[1]})")

    # Step 2: Perform CCA
    cca = CCA(n_components=n_components, max_iter=1000, tol=1e-6)

    try:
        # Fit CCA and transform both embeddings to canonical space
        rna_canonical, atac_canonical = cca.fit_transform(rna_embedding, atac_embedding)
        
        if verbose:
            print(f" CCA completed: {rna_canonical.shape}")

    except Exception as e:
        print(f"    CCA failed: {e}")
        print(f"    This may be due to numerical instability in the data.")
        print(f"    Troubleshooting suggestions:")
        print(f"      1. Check for NaN/Inf values in embeddings")
        print(f"      2. Try reducing n_components")
        print(f"      3. Ensure sufficient cell numbers (n_cells >> n_dims)")
        raise RuntimeError(
            f"CCA integration failed. Cannot proceed with simple concatenation "
            f"as it violates the dimension filtering logic. "
            f"Expected output: n_significant_dims (>=4), Got fallback: {rna_embedding.shape[1] + atac_embedding.shape[1]}"
        ) from e
    
    # Step 4: Calculate Pearson correlation and p-values for each canonical variable pair
    correlations = []
    p_values = []
    
    for i in range(n_components):
        # Calculate Pearson correlation between RNA and ATAC canonical variables
        corr, p_val = pearsonr(rna_canonical[:, i], atac_canonical[:, i])
        correlations.append(corr)
        p_values.append(p_val)
    
    correlations = np.array(correlations)
    p_values = np.array(p_values)
    
    if verbose:
        print(f"    Correlation range: [{correlations.min():.3f}, {correlations.max():.3f}]")
        print(f"    Mean correlation: {correlations.mean():.3f}")
    
    # Step 5: BH (Benjamini-Hochberg) correction for multiple testing
    reject, p_values_corrected, _, _ = multipletests(
        p_values, 
        alpha=p_value_threshold, 
        method='fdr_bh'  # Benjamini-Hochberg procedure
    )
    
    # Step 6: Select significant dimensions
    significant_dims = np.where(reject)[0]
    if len(significant_dims) == 0:
    # If no significant dimensions, use top correlated dimensions
        n_top = max(4, n_components // 10)  # Ensure at least 4 dimensions
        significant_dims = np.argsort(np.abs(correlations))[-n_top:]
        if verbose:
            print(f"    No significant dimensions at p<{p_value_threshold}")
            print(f"    Using top {n_top} dimensions by correlation (minimum 4)")
    elif len(significant_dims) < 4:
        # Supplement significant dims with top-correlated dims to reach minimum 4
        if verbose:
            print(f"    Only {len(significant_dims)} significant dimension(s), supplementing to minimum 4")
        top_by_corr = np.argsort(np.abs(correlations))[-4:]
        significant_dims = np.union1d(significant_dims, top_by_corr)[:4]
        if verbose:
            print(f"    Using {len(significant_dims)} dimensions (significant + top-correlated)")
    
    # Step 7: Fuse canonical variables by summation
    # For each significant dimension, combine RNA and ATAC canonical variables
    rna_selected = rna_canonical[:, significant_dims]
    atac_selected = atac_canonical[:, significant_dims]
    
    # Sum the canonical variables (as described in MOJITOO method)
    joint_embedding = rna_selected + atac_selected
    
    if verbose:
        print(f" Joint embedding created: {joint_embedding.shape}")
        print(f" Final dimensions: {joint_embedding.shape[1]} (from {n_components} CCA components)")
    
    return joint_embedding
