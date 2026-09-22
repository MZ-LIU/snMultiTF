"""
Gene Selection Module

Implements gene selection algorithms equivalent to R scMEGA's SelectGenes function.
This module provides peak-to-gene based gene selection for GRN inference.

Key functions:
- select_genes: Main gene selection function (SelectGenes equivalent)
- _apply_variance_cutoff: Apply variance-based filtering on trajectory data

References:
- Original R code: SelectGenes function in scMEGA R/select_tf_gene.R
- ArchR: Granja et al. (2021) ArchR is a scalable software package for integrative analysis
"""

import pandas as pd
import anndata as ad
from typing import Optional, Dict, Any, List, Union
import warnings

from ..trajectory_analysis.pseudotime_analysis import get_trajectory_data
from ..peak_gene_linking.peak_to_gene import link_peaks_to_genes


# Try to import MuData, set to None if not available
try:
    import muon as mu
    MUON_AVAILABLE = True
except ImportError:
    mu = None
    MUON_AVAILABLE = False


def select_genes(data: Union[ad.AnnData, 'mu.MuData'],
                atac_assay: str = "ATAC",
                rna_assay: str = "RNA",
                var_cutoff_gene: float = 0.9,
                trajectory_name: str = "Trajectory",
                distance_cutoff: int = 2000,  # Add this line: Minimum distance threshold when filtering
                group_every: int = 1,
                cor_cutoff: float = 0.0,
                fdr_cutoff: float = 1e-04,
                return_heatmap: bool = False,
                genome: str = "hg38",
                smooth_window: int = 7,
                peak_annotation: Optional[pd.DataFrame] = None,
                gene_annotation: Optional[pd.DataFrame] = None) -> Dict[str, Any]:
    """
    Select genes for GRN inference based on peak-to-gene links.
    
    This function replicates scMEGA's SelectGenes:
    ```R
    gene_results <- SelectGenes(
        object = obj,
        atac.assay = "ATAC",
        rna.assay = "RNA",
        var.cutoff.gene = 0.9,
        trajectory.name = "Trajectory",
        distance.cutoff = 2000,
        groupEvery = 1,
        cor.cutoff = 0,
        fdr.cutoff = 1e-04,
        return.heatmap = TRUE,
        genome = "hg38"
    )
    ```
    
    Workflow:
    1. Get RNA trajectory data using GetTrajectory
    2. Get ATAC trajectory data using GetTrajectory
    3. Select top variable genes using TrajectoryHeatmap (varCutOff = 0.9, top 10%)
    4. Process ATAC data without variance filtering (varCutOff = 0)
    5. Link peaks to genes using PeakToGene
    6. Filter by distance, correlation, and FDR
    
    Args:
        data: MultiomeData or AnnData object with both RNA and ATAC data
        atac_assay: The assay name for chromatin accessibility. Default: "ATAC"
        rna_assay: The assay name for gene expression. Default: "RNA"
        var_cutoff_gene: The cutoff of variation to select genes. Default: 0.9 (top 10%)
        trajectory_name: The trajectory name used for analysis
        max_distance: Maximum search distance for PeakToGene (bp). Default: 250000
        distance_cutoff: Minimum distance to filter peak-gene pairs (bp). Default: 2000
                        Filters out peaks too close to genes (promoter region)
        group_every: The number of sequential percentiles to group together
        cor_cutoff: The cutoff of peak-to-gene correlation. Default: 0
        fdr_cutoff: The cutoff of peak-to-gene p-value. Default: 1e-04
        return_heatmap: Whether or not return the heatmap for visualization
        genome: Reference genome ("hg38", "hg19", "mm10", "mm9")
        smooth_window: Smoothing window size for trajectory. Default: 7
        peak_annotation: Optional peak annotations (chr, start, end)
        gene_annotation: Optional gene annotations (chr, start, end, strand)
        
    Returns:
        Dictionary containing:
        - 'p2g': DataFrame with peak-to-gene links
        - 'heatmap': Optional heatmap data if return_heatmap=True
        - 'selected_genes': List of selected gene names
        - 'selected_peaks': List of selected peak names
        - 'trajectory_rna': RNA trajectory matrix
        - 'trajectory_atac': ATAC trajectory matrix
        
    Example:
        ```python
        # Equivalent to R: SelectGenes(obj, var.cutoff.gene = 0.9)
        gene_results = select_genes(
            multiome,
            var_cutoff_gene=0.9,
            trajectory_name="Trajectory",
            distance_cutoff=2000,
            cor_cutoff=0.0,
            fdr_cutoff=1e-04
        )
        
        # Access results
        p2g_links = gene_results['p2g']
        selected_genes = gene_results['selected_genes']
        ```
    """
    print(f"Selecting genes for GRN inference (R SelectGenes equivalent)...")
    print(f"Parameters: var_cutoff={var_cutoff_gene}, distance>{distance_cutoff}, "
          f"cor>{cor_cutoff}, FDR<{fdr_cutoff}")
    
    if MUON_AVAILABLE and isinstance(data, mu.MuData):
        if 'rna' not in data.mod or 'atac' not in data.mod:
            raise ValueError("MuData must contain both 'rna' and 'atac' modalities")
        rna_data = data['rna']
        atac_data = data['atac']
    else:
        raise TypeError(
            "select_genes requires MuData object with separate RNA and ATAC data.\n"
            "Please provide MuData object: mdata = mu.MuData({'rna': rna_adata, 'atac': atac_adata})\n"
            "Single AnnData objects are not supported for gene selection."
        )
    if MUON_AVAILABLE and isinstance(data, mu.MuData):
        if trajectory_name not in data.obs.columns and trajectory_name not in rna_data.obs.columns:
            raise ValueError(f"Trajectory '{trajectory_name}' not found in MuData. "
                            "Please run add_trajectory first.")
    else:
        if trajectory_name not in rna_data.obs.columns:
            raise ValueError(f"Trajectory '{trajectory_name}' not found. "
                            "Please run add_trajectory first.")
    
    print(f"\n=== Step 1: Get RNA trajectory data ===")
    # Get RNA trajectory data (equivalent to R's GetTrajectory for RNA)
    traj_rna_result = get_trajectory_data(
        rna_data,
        trajectory_name=trajectory_name,
        assay=rna_assay,
        slot="counts",
        group_every=group_every,
        log2_norm=True,
        scale_to=10000,
        smooth_window=smooth_window,
        return_matrix=False
    )
    
    traj_rna = traj_rna_result['smooth_matrix'] if 'smooth_matrix' in traj_rna_result else traj_rna_result['group_matrix']
    print(f"RNA trajectory matrix: {traj_rna.shape} (genes x time bins)")
    
    print(f"\n=== Step 2: Get ATAC trajectory data ===")
    # Get ATAC trajectory data (equivalent to R's GetTrajectory for ATAC)
    traj_atac_result = get_trajectory_data(
        atac_data,
        trajectory_name=trajectory_name,
        assay=atac_assay,
        slot="X",
        group_every=group_every,
        log2_norm=False,   # ATAC LSI data has negative values; log2 not applicable
        scale_to=None,     # ATAC LSI data; depth normalization not applicable
        smooth_window=smooth_window,
        return_matrix=False
    )
    
    traj_atac = traj_atac_result['smooth_matrix'] if 'smooth_matrix' in traj_atac_result else traj_atac_result['group_matrix']
    print(f"ATAC trajectory matrix: {traj_atac.shape} (peaks x time bins)")
    
    print(f"\n=== Step 3: Select top variable genes (TrajectoryHeatmap varCutOff={var_cutoff_gene}) ===")
    # Apply variance cutoff to RNA data (equivalent to TrajectoryHeatmap with varCutOff=0.9)
    # This selects the top 10% (1 - 0.9) most variable genes across the trajectory
    group_mat_rna = _apply_variance_cutoff(
        traj_rna,
        var_cutoff=var_cutoff_gene,
        max_features=None
    )
    print(f"Selected {len(group_mat_rna)} highly variable genes (top {int((1-var_cutoff_gene)*100)}%)")
    
    print(f"\n=== Step 4: Process ATAC data (no variance filtering) ===")
    # Apply variance cutoff to ATAC data (varCutOff=0 means keep all)
    group_mat_atac = _apply_variance_cutoff(
        traj_atac,
        var_cutoff=0,
        max_features=len(traj_atac)
    )
    print(f"Using all {len(group_mat_atac)} peaks")
    
    print(f"\n=== Step 5: Link peaks to genes (PeakToGene) ===")
    # Link peaks to genes based on trajectory matrices (R PeakToGene equivalent)
    df_p2g = link_peaks_to_genes(
        peak_mat=group_mat_atac,
        gene_mat=group_mat_rna,
        genome=genome,
        peak_annotation=peak_annotation,
        gene_annotation=gene_annotation,
        max_distance=250000  # Suggested addition: Pass parameters explicitly
    )
    
    print(f"Initial peak-gene links: {len(df_p2g)}")
    
    print(f"\n=== Step 6: Filter peak-gene links ===")
    # Filter by distance, correlation, and FDR (equivalent to R's subset operations)
    df_p2g_filtered = df_p2g[
        (df_p2g['distance'] > distance_cutoff) &
        (df_p2g['Correlation'] > cor_cutoff) &
        (df_p2g['FDR'] < fdr_cutoff)
    ].copy()
    
    print(f"Filtered peak-gene links: {len(df_p2g_filtered)}")
    print(f"  - After distance filter (>{distance_cutoff}bp): {len(df_p2g[df_p2g['distance'] > distance_cutoff])}")
    print(f"  - After correlation filter (>{cor_cutoff}): {len(df_p2g[df_p2g['Correlation'] > cor_cutoff])}")
    print(f"  - After FDR filter (<{fdr_cutoff}): {len(df_p2g[df_p2g['FDR'] < fdr_cutoff])}")
    
    if len(df_p2g_filtered) == 0:
        warnings.warn("No peak-gene links passed filtering criteria. Consider relaxing thresholds.")
        return {
            'p2g': pd.DataFrame(),
            'selected_genes': [],
            'selected_peaks': [],
            'trajectory_rna': traj_rna,
            'trajectory_atac': traj_atac
        }
    
    # Subset trajectory matrices to filtered peaks and genes
    valid_pairs = df_p2g_filtered[
        df_p2g_filtered['peak'].isin(traj_atac.index) & 
        df_p2g_filtered['gene'].isin(traj_rna.index)
    ]
    if len(valid_pairs) == 0:
        print(f"   There are no valid peak-gene pairs in the trajectory data")
        traj_atac_filtered = pd.DataFrame()
        traj_rna_filtered = pd.DataFrame()
    else:
        traj_atac_filtered = traj_atac.loc[valid_pairs['peak']]
        traj_rna_filtered = traj_rna.loc[valid_pairs['gene']]
    
    print(f"\n genetic screening completed:")
    print(f"  Number of unique genes selected: {df_p2g_filtered['gene'].nunique()}")
    print(f"  The number of unique peaks selected: {df_p2g_filtered['peak'].nunique()}")
    print(f"  Peak-gene link number: {len(df_p2g_filtered)}")
    print(f"  Trajectory data used for visualization: {len(traj_atac_filtered)} pair")
    
    # Prepare results
    results = {
        'p2g': df_p2g_filtered,
        'selected_genes': df_p2g_filtered['gene'].unique().tolist(),
        'selected_peaks': df_p2g_filtered['peak'].unique().tolist(),
        'trajectory_rna': traj_rna_filtered,
        'trajectory_atac': traj_atac_filtered,
        'parameters': {
            'var_cutoff_gene': var_cutoff_gene,
            'trajectory_name': trajectory_name,
            'distance_cutoff': distance_cutoff,
            'cor_cutoff': cor_cutoff,
            'fdr_cutoff': fdr_cutoff,
            'genome': genome
        }
    }
    
    heatmap_fig = None
    # Add heatmap data if requested
    if return_heatmap:
        print("\n  Generating gene selection heatmap...")
    
        try:
            from python_scmega.visualization.heatmaps import gene_heatmap
            
            
            heatmap_fig = gene_heatmap(
                chromatin_accessibility=traj_atac_filtered,
                gene_expression=traj_rna_filtered,
                scale_rows=True,
                limits=(-2, 2),
                label_rows=False,
                label_top_accessibility=10,
                label_top_expression=10,
                name1="Chromatin accessibility",
                name2="Gene expression"
            )
            results['heatmap'] = heatmap_fig
            print(f"   Heatmap generated for {len(df_p2g_filtered)} peak-gene pairs")
            
        # except Exception as e:
        #     print(f"   Heatmap generation failed: {e}")
        #     results['heatmap'] = None
        except Exception as e:
            print(f"   Heatmap generation failed in gene_selection: {e}")
            import traceback
            traceback.print_exc()  # Print complete stack information
            results['heatmap'] = None
    else:
        results['heatmap'] = None
                
             
    
    # Show top links
    if len(df_p2g_filtered) > 0:
        print(f"\nTop 5 peak-gene links:")
        top_links = df_p2g_filtered.nlargest(5, 'Correlation')
        for i, (_, row) in enumerate(top_links.iterrows()):
            print(f"  {i+1}. {row['peak']} -> {row['gene']}: "
                  f"distance={row['distance']:.0f}bp, r={row['Correlation']:.3f}, FDR={row['FDR']:.2e}")
        
    return results
    

def _apply_variance_cutoff(trajectory_matrix: pd.DataFrame,
                           var_cutoff: float = 0.9,
                           max_features: Optional[int] = None) -> pd.DataFrame:
    """
    Apply variance-based filtering on trajectory matrix.
    
    This replicates the R TrajectoryHeatmap variance filtering logic:
    - Calculate variance for each feature (gene/peak) across trajectory
    - Convert to quantiles
    - Select features above variance quantile cutoff
    
    R code equivalent:
    ```R
    varQ <- ArchR:::.getQuantiles(matrixStats::rowVars(mat))
    mat <- mat[order(varQ, decreasing = TRUE), ]
    n <- (1 - varCutOff) * nrow(mat)
    mat <- mat[head(seq_len(nrow(mat)), n), ]
    ```
    
    Args:
        trajectory_matrix: DataFrame with features as rows, time bins as columns
        var_cutoff: Variance quantile cutoff (0.9 means keep top 10% most variable)
        max_features: Maximum number of features to retain
        
    Returns:
        Filtered trajectory matrix
    """
    mat = trajectory_matrix.copy()
    
    # Remove rows with NA values
    mat = mat.dropna()
    
    # Remove rows with zero standard deviation
    row_stds = mat.std(axis=1)
    mat = mat[row_stds != 0]
    
    if len(mat) == 0:
        raise ValueError("No features remaining after filtering NA and zero-variance rows")
    
    # Calculate variance for each row (feature)
    row_vars = mat.var(axis=1)
    
    # Calculate variance quantiles (equivalent to ArchR:::.getQuantiles)
    var_quantiles = row_vars.rank(pct=True)
    
    # Sort by variance quantile (descending)
    mat = mat.loc[var_quantiles.sort_values(ascending=False).index]
    
    # Determine number of features to keep
    if var_cutoff is None and max_features is None:
        n = len(mat)
    elif var_cutoff is None:
        n = max_features
    elif max_features is None:
        n = int((1 - var_cutoff) * len(mat))
    else:
        n = min(int((1 - var_cutoff) * len(mat)), max_features)
    
    n = min(n, len(mat))
    
    # Select top n features
    mat_filtered = mat.iloc[:n, :]
    
    print(f"  Variance filtering: {len(trajectory_matrix)} -> {len(mat_filtered)} features")
    if n < len(trajectory_matrix):
        print(f"  Kept features with variance quantile > {var_cutoff:.2f}")
    
    return mat_filtered



# # Keep backward compatibility - alias for main function
# SelectGenes = select_genes

