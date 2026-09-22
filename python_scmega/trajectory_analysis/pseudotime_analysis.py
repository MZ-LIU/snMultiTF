"""
Pseudotime Analysis Module

Implements pseudotime analysis functions equivalent to Seurat's GetTrajectory
and related functions for analyzing gene expression dynamics along trajectories.

Key functions:
- get_trajectory_data: Extract trajectory data for analysis
- analyze_pseudotime_dynamics: Analyze gene expression dynamics along pseudotime
- fit_trajectory_curves: Fit smooth curves to gene expression along trajectory
- identify_trajectory_markers: Identify genes with significant trajectory dynamics
- plot_trajectory_heatmap: Visualize gene expression along trajectory

References:
- Original R code: GetTrajectory function in scMEGA
- Monocle 3: Cao et al. (2019) The single-cell transcriptional landscape
- tradeSeq: Van den Berge et al. (2020) Trajectory-based differential expression
"""

import numpy as np
import pandas as pd
import anndata as ad
from typing import Optional, Dict, Any, List, Union, Tuple
import warnings
from scipy import sparse  
from .trajectory_inference import center_roll_mean  # Import the utility function

# Try to import MuData, set to None if not available
try:
    import muon as mu
    MUON_AVAILABLE = True
except ImportError:
    mu = None
    MUON_AVAILABLE = False

def _extract_chromvar_as_adata(data: Union[ad.AnnData, 'mu.MuData']) -> ad.AnnData:
    """
    Extracts chromVAR data from MuData, returning the AnnData of the chromVAR modality.
    The current version of chromVAR data is stored as an independent chromvar mode, corresponding to the R version:
    object@assays[["chromvar"]]
    
    Data structure (created by run_chromvar_r):
    - mdata['chromvar']: independent AnnData modal
    - mdata['chromvar'].X: TF activity matrix (cells × TFs)
    - mdata['chromvar'].layers['z_scores']: Z-scores
    - mdata['chromvar'].var_names: TF/motif names
    
    Args:
        data: MuData or AnnData object
    Returns:
        AnnData object of chromVAR modal, containing:
        - X: TF activity matrix
        - var.index: motif/TF name
        - obs: cell metadata (including trajectory information)
    Raises:
        ValueError: if chromVAR modal not found
        TypeError: if input is not a MuData object
        
    Example:
        ```python
        # Get chromVAR data for trajectory analysis
        chromvar_adata = _extract_chromvar_as_adata(mdata)
        # chromvar_adata.X contains TF activity matrix
        # chromvar_adata.var_names contains TF names
        ```
    """

    if not MUON_AVAILABLE or not isinstance(data, mu.MuData):
        raise TypeError(
            "chromVAR data extraction requires MuData object!\n"
            "Current data structure stores chromVAR as an independent modality.\n"
            "Expected: mdata['chromvar'] where mdata is a MuData object.\n"
        )
    
    if 'chromvar' not in data.mod:
        raise ValueError(
            f"chromVAR modality not found in MuData!\n"
            f"Available modalities: {list(data.mod.keys())}\n"
            f"Please run run_chromvar_r() first."
        )
    
    chromvar_adata = data['chromvar']
    
    if chromvar_adata.X is None or chromvar_adata.X.shape[0] == 0:
        raise ValueError(
            "chromVAR modality exists but is empty!\n"
            f"Shape: {chromvar_adata.shape}\n"
            "Please re-run chromVAR analysis."
        )
    if chromvar_adata.n_vars == 0:
        raise ValueError(
            "chromVAR modality has no features (TFs/motifs)!\n"
            "Please check chromVAR analysis results."
        )
    if chromvar_adata.var_names is None or len(chromvar_adata.var_names) == 0:
        warnings.warn(
            "chromVAR modality has no var_names (TF/motif names)!\n"
            "Using numeric indices instead. This may affect TF-gene matching.\n"
            "Expected: mdata['chromvar'].var_names should contain TF/motif names."
        )
    
    # Note: the chromvar modality should share the same obs (cell metadata) as the rna modality
    if 'rna' in data.mod:
        if chromvar_adata.n_obs != data['rna'].n_obs:
            warnings.warn(
                f"chromVAR cell count ({chromvar_adata.n_obs}) differs from "
                f"RNA cell count ({data['rna'].n_obs})!\n"
                "This may indicate data inconsistency."
            )
    print(f"   chromVAR data extracted from mdata['chromvar']")
    print(f"    - Shape: {chromvar_adata.shape} (cells × TFs)")
    print(f"    - TF/motif names: {chromvar_adata.n_vars} features")
    if chromvar_adata.layers:
        print(f"    - Available layers: {', '.join(chromvar_adata.layers.keys())}")
    
    return chromvar_adata
    


def get_trajectory_data(data: Union[ad.AnnData, 'mu.MuData'],
                       trajectory_name: str = "Trajectory",
                       assay: str = "RNA",
                       slot: str = "X",
                       group_every: int = 1,
                       log2_norm: bool = True,
                       scale_to: Optional[float] = 10000,
                       smooth_window: Optional[int] = 11,
                       genes: Optional[List[str]] = None,
                       motif_names: Optional[Dict[str, str]] = None,
                       return_matrix: bool = False) -> Union[Dict[str, Any], pd.DataFrame]:
    """
    Get data along the trajectory using Seurat object as input - R-compatible version.
    
    This function is a faithful Python implementation of the R GetTrajectory function:
    ```R
    trajMM <- GetTrajectory(
        object = obj,
        trajectory.name = "Trajectory",
        assay = "chromvar",      # or "RNA"
        slot = "data",
        groupEvery = 1,
        log2Norm = FALSE,        # FALSE for chromVAR activity
        scaleTo = 10000,
        smoothWindow = 7
    )
    
    # CRITICAL R step: Map row names to motif names for chromVAR data
    if (assay == "chromvar") {
        rownames(trajMM) <- object@assays[["ATAC"]]@motifs@motif.names
    }
    ```
    
    **IMPORTANT for chromVAR/TF analysis**:
    When using chromVAR assay (TF activity), R version replaces row names with motif names.
    This is essential for matching TF activity to TF gene names in correlation analysis.
    
    Args:
        data: MultiomeData or AnnData object with trajectory information
        trajectory_name: The name of trajectory inferred by AddTrajectory
        assay: Which assay is used for data collection ("RNA", "ATAC", "chromvar")
        slot: Name of slot used to collect data ("X", "raw", "counts")
        group_every: The number of sequential percentiles to group together when generating a trajectory.
                    This is similar to smoothing via a non-overlapping sliding window across pseudo-time.
        log2_norm: Whether or not the data should be log2 transformed. This can be set to TRUE
                  if gene expression or chromatin accessibility data is used,
                  otherwise FALSE if data from TF activity is used.
        scale_to: Once the sequential trajectory matrix is created,
                 each column in that matrix will be normalized to a column sum indicated by scaleTo.
        smooth_window: An integer value indicating the smoothing window in size
                      (relative to groupEvery) for the sequential trajectory matrix to better reveal
                      temporal dynamics.
        genes: Specific genes to extract (if None, use all genes)
        motif_names: **Deprecated** - motif names are now automatically extracted from uns
        return_matrix: If True, return only the matrix (smooth_matrix if available, else group_matrix).
                      If False (default), return the complete dictionary with all metadata.
                      Similar to R's TrajectoryHeatmap(..., returnMatrix=TRUE).
        
    Returns:
        If return_matrix=False: Dictionary containing trajectory matrices and metadata (default)
        If return_matrix=True: pd.DataFrame matrix only (smoothed if available, otherwise grouped)
        
    Example:
        ```python
        # Return complete dictionary (default):
        trajectory_result = get_trajectory_data(
            multiome, 
            trajectory_name="Trajectory",
            assay="RNA",
            smooth_window=7,
            log2_norm=True
        )
        # Access: trajectory_result['smooth_matrix']
        
        # Return matrix directly (convenient for most use cases):
        traj_matrix = get_trajectory_data(
            multiome,
            trajectory_name="Trajectory",
            assay="chromvar",
            smooth_window=7,
            log2_norm=False,
            scale_to=None,
            return_matrix=True # Return the matrix directly
        )
        # Can use directly: traj_matrix.shape
        ```
    """
    print(f"Creating Trajectory Group Matrix for '{trajectory_name}'...")
    
    # ========================================================================
    # ========================================================================
    # If it is a chromvar assay, first call the auxiliary function to convert the data, and then continue the standard process
    if assay.lower() in ["chromvar", "tf"]:
        print(f"  Detected chromVAR/TF assay - extracting from mdata['chromvar']...")  # ← Modify here
        data = _extract_chromvar_as_adata(data)
        assay = "RNA"
        print(f"   ChromVAR data extracted from independent modality")  #  # ← Optional: Update here
    
    # ========================================================================
    # ========================================================================
    if MUON_AVAILABLE and isinstance(data, mu.MuData):
        if assay == "RNA" and 'rna' in data.mod:
            adata = data['rna']
        elif assay == "ATAC" and 'atac' in data.mod:
            adata = data['atac']
        else:
            raise ValueError(f"Assay '{assay}' not found in MuData modalities: {list(data.mod.keys())}")
    else:
        adata = data
    
    # Validate assay (for MuData only)
    if MUON_AVAILABLE and isinstance(data, mu.MuData):
        if assay == "RNA" and 'rna' not in data.mod:
            raise ValueError(f"RNA modality not found in MuData")
        elif assay == "ATAC" and 'atac' not in data.mod:
            raise ValueError(f"ATAC modality not found in MuData")
    
    # Check trajectory information
    if trajectory_name not in adata.obs.columns:
        raise ValueError(f"Cannot find trajectory {trajectory_name}!")
    
    # Get trajectory values
    trajectory_values = adata.obs[trajectory_name].copy()
    
    # Remove cells with NaN trajectory values
    trajectory_mask = ~pd.isna(trajectory_values)
    trajectory_values = trajectory_values[trajectory_mask]
    
    if len(trajectory_values) == 0:
        raise ValueError("No cells found in trajectory")
    
    # Validate trajectory values (must be numeric and in 0-100 range)
    if not pd.api.types.is_numeric_dtype(trajectory_values):
        raise ValueError("Trajectory must be numeric. Did you add the trajectory with add_trajectory?")
    
    if not ((trajectory_values >= 0) & (trajectory_values <= 100)).all():
        raise ValueError("Trajectory values must be between 0 and 100. Did you add the trajectory with add_trajectory?")
    
    # trajectory_data = adata[trajectory_mask].copy()
    # print(f"Found {len(trajectory_data)} cells in trajectory")
    # print(f"Trajectory range: {trajectory_values.min():.1f} - {trajectory_values.max():.1f}")
    # Do not copy the entire object, only get the effective cell index (avoid deep copy uns causing memory errors)
    valid_cell_indices = adata.obs.index[trajectory_mask]
    n_trajectory_cells = len(valid_cell_indices)
    print(f"Found {n_trajectory_cells} cells in trajectory")
    print(f"Trajectory range: {trajectory_values.min():.1f} - {trajectory_values.max():.1f}")


    # Create trajectory bins (matching R logic exactly)
    breaks = np.arange(0, 100 + group_every, group_every)  # seq(0, 100, groupEvery)
    
    print(f"Creating {len(breaks)-1} trajectory groups with groupEvery={group_every}")
    
    # Group cells by trajectory bins
    group_list = []
    group_names = []
    
    for i in range(1, len(breaks)):
        # Find cells in this bin: trajectory[, 1] > breaks[x-1] & trajectory[, 1] <= breaks[x]
        bin_mask = (trajectory_values > breaks[i-1]) & (trajectory_values <= breaks[i])
        # bin_cells = trajectory_data.obs.index[bin_mask].tolist()
        bin_cells = valid_cell_indices[bin_mask].tolist()
        if len(bin_cells) > 0:
            group_list.append(bin_cells)
            group_names.append(f"T.{breaks[i-1]}_{breaks[i]}")
        else:
            group_list.append([])
            group_names.append(f"T.{breaks[i-1]}_{breaks[i]}")
    
    print(f"Created {len(group_list)} trajectory groups")
    

    if genes is not None:
        available_genes = [g for g in genes if g in adata.var.index]
        if len(available_genes) == 0:
            raise ValueError("None of the specified genes found in data")
        
        if len(available_genes) < len(genes):
            missing = set(genes) - set(available_genes)
            warnings.warn(f"Genes not found: {missing}")
        
        gene_mask = adata.var.index.isin(available_genes)
        print(f"Using {len(available_genes)} specified genes")
    else:
        gene_mask = np.ones(adata.n_vars, dtype=bool)
        available_genes = adata.var.index.tolist()
        print(f"Using all {adata.n_vars} genes")

    # Get expression data based on slot
    if slot == "X":
        data_use = adata.X
    elif slot == "raw":
        if adata.raw is not None:
            data_use = adata.raw.X
        else:
            print("Raw data not available, using .X")
            data_use = adata.X
    elif slot == "counts":
        if 'counts' in adata.layers:
            data_use = adata.layers['counts']
        else:
            print("Counts layer not available, using .X")
            data_use = adata.X
    else:
        raise ValueError(f"Unknown slot: {slot}")

    # Apply gene filtering
    if genes is not None:
        data_use = data_use[:, gene_mask]

    # Keep sparse matrices in sparse format - do not convert to dense array to save memory
    is_sparse = sparse.issparse(data_use)
        
    # Calculate group means for each trajectory bin
    # group_matrix = np.zeros((trajectory_data.n_vars, len(group_list)))
    gene_names = adata.var.index[gene_mask] if genes is not None else adata.var.index
    group_matrix = np.zeros((len(gene_names), len(group_list)))

    for i, cell_group in enumerate(group_list):
        if len(cell_group) > 0:
            # Get indices of cells in this group
            # group_indices = [trajectory_data.obs.index.get_loc(cell) for cell in cell_group 
            #             if cell in trajectory_data.obs.index]
            group_indices = [adata.obs.index.get_loc(cell) for cell in cell_group 
                        if cell in adata.obs.index]

            if len(group_indices) > 0:
                # Calculate mean expression for this group
                group_data = data_use[group_indices, :]
                
                if is_sparse:
                    # For sparse matrices, use the .mean() method to calculate column means, returning a matrix of shape (1, n_features)
                    group_matrix[:, i] = np.array(group_data.mean(axis=0)).flatten()
                else:
                    group_matrix[:, i] = np.mean(group_data, axis=0)
        # If no cells in group, values remain 0
    
    # # Create group matrix DataFrame
    # group_df = pd.DataFrame(
    #     group_matrix,
    #     index=trajectory_data.var.index,
    #     columns=group_names
    # )
        # Create group matrix DataFrame
    group_df = pd.DataFrame(
        group_matrix,
        index=gene_names,
        columns=group_names
    )
    
    # ========================================================================
    # ========================================================================
    if motif_names is not None:
        warnings.warn(
            "The 'motif_names' parameter is deprecated. "
            "Motif names are now automatically extracted from uns['chromvar']['motif_names']. "
            "This parameter will be ignored.",
            DeprecationWarning
        )
    
    print(f"Group matrix shape: {group_df.shape}")
    
    # Scale normalization (scaleTo)
    if scale_to is not None:
        if (group_df < 0).any().any():
            print("Some values are below 0, this could be the Motif activity matrix in which scaleTo should be set = None.")
            print("Continuing without depth normalization!")
        else:
            # Normalize each column to scaleTo
            col_sums = group_df.sum(axis=0)
            col_sums[col_sums == 0] = 1  # Avoid division by zero
            group_df = group_df.div(col_sums, axis=1) * scale_to
            print(f"Applied scaleTo normalization: {scale_to}")
    
    # Log2 normalization
    if log2_norm:
        if (group_df < 0).any().any():
            print("Some values are below 0, this could be a Motif activity matrix in which log2Norm should be set = FALSE.")
            print("Continuing without log2 normalization!")
        else:
            group_df = np.log2(group_df + 1)
            print("Applied log2 normalization")
    
    # Smoothing
    result_dict = {
        'group_matrix': group_df,
        'group_names': group_names,
        'trajectory_breaks': breaks,
        #'n_trajectory_cells': len(trajectory_data),
        'n_trajectory_cells': n_trajectory_cells,
        'parameters': {
            'trajectory_name': trajectory_name,
            'assay': assay,
            'slot': slot,
            'group_every': group_every,
            'log2_norm': log2_norm,
            'scale_to': scale_to,
            'smooth_window': smooth_window
        }
    }
    
    if smooth_window is not None:
        print("Smoothing...")
        
        # Apply center-aligned rolling mean to each gene (row)
        smooth_matrix = np.zeros_like(group_matrix)
        
        for i in range(group_matrix.shape[0]):
            gene_trajectory = group_df.iloc[i].values
            smooth_matrix[i, :] = center_roll_mean(gene_trajectory, smooth_window)
        
        smooth_df = pd.DataFrame(
            smooth_matrix,
            index=group_df.index,  # Use the index of group_df (when using chromvar, it is already the motif names)
            columns=group_names
        )
        
        result_dict['smooth_matrix'] = smooth_df
        print(f"Applied smoothing with window size: {smooth_window}")
    
    print(f" Trajectory data processing completed")
    
    if return_matrix:
        # Return the smoothing matrix first, if not, return the grouping matrix
        return result_dict.get('smooth_matrix', result_dict['group_matrix'])
    else:
        return result_dict
