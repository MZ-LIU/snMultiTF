"""
TF Selection Module

Implements transcription factor (TF) selection algorithms equivalent to Seurat's 
SelectTFs function, with focus on trajectory-based TF identification and ranking.

Key functions:
- select_tfs: Main TF selection function


References:
- Original R code: SelectTFs function in scMEGA
- SCENIC: Aibar et al. (2017) SCENIC: single-cell regulatory network inference
- ChromVAR: Schep et al. (2017) chromVAR: inferring transcription-factor
"""

import numpy as np
import pandas as pd
import anndata as ad
from typing import Optional, Dict, Any, List, Union
from scipy.stats import pearsonr

from ..trajectory_analysis.pseudotime_analysis import get_trajectory_data

# Try to import MuData, set to None if not available
try:
    import muon as mu
    MUON_AVAILABLE = True
except ImportError:
    mu = None
    MUON_AVAILABLE = False


def select_tfs(data: Union[ad.AnnData, 'mu.MuData'],
               assay: str = "RNA",
               n_tfs: Optional[int] = None,
               tf_database: Optional[List[str]] = None,
               return_dict: bool = True,
               correlation_threshold: float = 0.3,
               p_cutoff: float = 0.01,
               trajectory_name: str = "Trajectory",
               tf_assay: str = "chromvar",
               group_every: int = 1,
               smooth_window: int = 7,
               return_heatmap: bool = False, 
               **kwargs) -> Union[List[str], Dict[str, Any]]:
    """
    Select transcription factors based on TF activity-expression correlation along trajectory (R-compatible).
    
    This function replicates scMEGA's SelectTFs R logic:
    ```R
    # R equivalent workflow (R/select_tf_gene.R:19-96):
    SelectTFs <- function(object, tf.assay = "chromvar", rna.assay = "RNA",
                          trajectory.name = "Trajectory", groupEvery = 1,
                          p.cutoff = 0.01, cor.cutoff = 0.3, ...) {
        
        # Step 1: Get TF activity along trajectory
        trajMM <- GetTrajectory(object, assay = tf.assay, trajectory.name = trajectory.name,
                               groupEvery = groupEvery, slot = "data", smoothWindow = 7, log2Norm = FALSE)
        rownames(trajMM) <- object@assays$ATAC@motifs@motif.names
        
        # Step 2: Get TF expression along trajectory  
        trajRNA <- GetTrajectory(object, assay = rna.assay, trajectory.name = trajectory.name,
                                groupEvery = groupEvery, slot = "data", smoothWindow = 7, log2Norm = TRUE)
        
        # Step 3: Calculate correlation between TF activity and expression
        df.cor <- GetCorrelation(trajMM, trajRNA)
        
        # Step 4: Filter by p-value and correlation threshold
        df.cor <- df.cor[df.cor$adj_p < p.cutoff & df.cor$correlation > cor.cutoff, ]
        
        return(df.cor)
    }
    ```
    
    Key differences from old implementation:
    - **NEW**: Uses trajectory-based correlation (GetTrajectory + GetCorrelation)
    - **NEW**: Calculates p-values and FDR-adjusted p-values
    - **CHANGED**: Default correlation_threshold = 0.3 (was 0.6)
    - **NEW**: Filters by both correlation AND p-value
    
    Args:
        data: MultiomeData or AnnData object with trajectory and chromVAR results
        assay: RNA assay to use ("RNA" for TF expression)
        n_tfs: Optional maximum number of TFs to select after filtering. By
            default, keep all TFs that pass the correlation and adjusted
            p-value thresholds. This argument is retained for backward
            compatibility when an explicit top-N selection is desired.
        tf_database: Custom list of TF genes to consider
        return_dict: If True, return R-compatible dict; if False, return TF list only
        correlation_threshold: Minimum correlation (default 0.3, matching R's cor.cutoff)
        p_cutoff: Maximum adjusted p-value (default 0.01, matching R's p.cutoff)
        trajectory_name: Name of trajectory to use (default "Trajectory")
        tf_assay: Assay for TF activity (default "chromvar", matching R's tf.assay)
        group_every: Trajectory binning parameter (default 1, matching R's groupEvery)
        smooth_window: Smoothing window size (default 7, matching R's smoothWindow)
        **kwargs: Additional parameters
        
    Returns:
        If return_dict=True (default):
            Dictionary with keys:
            - 'selected_tfs': List of selected TF names
            - 'tf_info': DataFrame with TF correlation data
            - 'tf_timepoint': Dict mapping {tf_name: timepoint_value}
        If return_dict=False:
            List of selected TF gene names
        
    Example:
        ```python
        # Standard usage:
        result = select_tfs(
            multiome_data,
            correlation_threshold=0.6
        )
        selected_tfs = result['selected_tfs']
        
        # Simple usage:
        tf_list = select_tfs(multiome_data, return_dict=False)
        ```
    """
    print("Selecting TFs by trajectory-based activity-expression correlation (R-compatible method)...")
    print(f"Parameters: cor.cutoff={correlation_threshold}, p.cutoff={p_cutoff}")
    
    # ========================================================================
    # Step 0: Validate inputs
    # ========================================================================
    if MUON_AVAILABLE and isinstance(data, mu.MuData):
        if 'rna' not in data.mod:
            raise ValueError("RNA modality not found in MuData")
        adata_rna = data['rna']
        adata_atac = data['atac'] if 'atac' in data.mod else None
    else:
        adata_rna = data
        adata_atac = data
    
    # Check if trajectory exists
    if trajectory_name not in adata_rna.obs.columns:
        raise ValueError(
            f"Trajectory '{trajectory_name}' not found. Please run add_trajectory() first.\n"
            f"Available columns: {list(adata_rna.obs.columns)}"
        )
    
    # Check if chromVAR activity data is available
    if MUON_AVAILABLE and isinstance(data, mu.MuData):
        if 'chromvar' not in data.mod:
            raise ValueError(
                "ChromVAR modality not found in MuData!\n"
                f"Available modalities: {list(data.mod.keys())}\n"
                "Please run chromVAR analysis first using pymega.run_chromvar_r()."
            )
        print(f"   Found chromVAR modality: {data['chromvar'].shape}")
    
    # Prefer reading from var_names of chromvar modal (standard location)
    if MUON_AVAILABLE and isinstance(data, mu.MuData) and 'chromvar' in data.mod:
        motif_names = list(data['chromvar'].var_names)
        print(f"   Using motif names from data['chromvar'].var_names")
    elif MUON_AVAILABLE and isinstance(data, mu.MuData) and 'atac' in data.mod:
        if 'motifs' in data['atac'].uns and 'motif_names' in data['atac'].uns['motifs']:
            motif_names = data['atac'].uns['motifs']['motif_names']
            print(f"  Using motif names from data['atac'].uns['motifs']['motif_names']")
        else:
            raise ValueError(
                "Motif names not found in data['atac'].uns['motifs'].\n"
                "Please ensure add_motifs_r() was run successfully."
            )
    else:
        raise ValueError(
            "Cannot access motif names. Expected MuData structure with 'chromvar' or 'atac' modality."
        )

    print(f"Found {len(motif_names)} motifs from chromVAR analysis")
    
    # ========================================================================
    # Step 1: Get TF activity along trajectory (R: trajMM <- GetTrajectory(..., assay=tf.assay))
    # ========================================================================
    print(f"\nStep 1: Getting TF activity along trajectory (assay='{tf_assay}')...")
    print(f"  groupEvery={group_every}, smoothWindow={smooth_window}, log2Norm=FALSE")
    
    # from ..trajectory_analysis.pseudotime_analysis import get_trajectory_data
    
    # Get TF activity trajectory (chromVAR)
    # R equivalent: trajMM <- GetTrajectory(object, assay = "chromvar", log2Norm = FALSE, scaleTo = NULL)
    # FIXED: chromVAR data is stored in obsm['X_chromvar'], not in main X matrix
    try:
        print(f"Step 1: Getting TF activity along trajectory (assay='chromvar')...")
        print(f"  groupEvery={group_every}, smoothWindow={smooth_window}, log2Norm=FALSE")
        
        traj_activity = get_trajectory_data(
            data,  # Directly pass in the original data (MultiomeData or AnnData)
            trajectory_name=trajectory_name,
            assay="chromvar",  # Specify the use of chromvar assay
            slot="X",  # chromVAR data will be placed in X of temporary AnnData
            group_every=group_every,
            log2_norm=False,               # FALSE for chromVAR activity
            scale_to=None,                 # NULL for chromVAR activity  
            smooth_window=smooth_window,
            return_matrix=True
        )
        
        # CRITICAL: Apply motif name mapping (R equivalent: rownames(trajMM) <- motif.names)
        # This ensures row names match TF gene names for correlation
        if hasattr(traj_activity, 'index'):
            # If it's a DataFrame, update the index
            # Map feature indices to motif names
            if len(traj_activity) == len(motif_names):
                traj_activity.index = motif_names
                print(f"   Applied motif names to TF activity matrix")
            else:
                print(f"  Warning: Motif count mismatch: {len(traj_activity)} features vs {len(motif_names)} motifs")
        
        print(f"   TF activity trajectory shape: {traj_activity.shape}")
        
    except Exception as e:
        raise ValueError(f"Failed to get TF activity trajectory: {e}")
    
    # ========================================================================
    # Step 2: Get TF expression along trajectory (R: trajRNA <- GetTrajectory(..., assay=rna.assay))
    # ========================================================================
    print(f"\nStep 2: Getting TF expression along trajectory (assay='RNA')...")
    print(f"  groupEvery={group_every}, smoothWindow={smooth_window}, log2Norm=TRUE")
    
    try:
        traj_expression = get_trajectory_data(
            data,
            trajectory_name=trajectory_name,
            assay="RNA",
            slot="counts",
            group_every=group_every,
            log2_norm=True,  # TRUE for RNA expression
            scale_to=10000,  # Standard normalization
            smooth_window=smooth_window,
            return_matrix=True
        )
        
        print(f"   TF expression trajectory shape: {traj_expression.shape}")
        
    except Exception as e:
        raise ValueError(f"Failed to get TF expression trajectory: {e}")
    
    # ========================================================================
    # Step 3: Calculate correlations (R: df.cor <- GetCorrelation(trajMM, trajRNA))
    # ========================================================================
    print(f"\nStep 3: Calculating TF activity-expression correlations along trajectory...")

    # R equivalent: features.use <- intersect(rownames(mat1), rownames(mat2))
    # Try direct intersection first (assuming the RNA row names are already gene symbols, consistent with R)
    motif_index = pd.Index(map(str, traj_activity.index))
    expr_index = pd.Index(map(str, traj_expression.index))

    available_tfs = list(set(motif_index) & set(expr_index))

    # If there is no intersection and there is a gene symbol column, map the RNA row name to a symbol and try again.
    if len(available_tfs) == 0:
        print(f"  No direct intersection found. Checking for gene symbol column...")
        symbol_col = next(
            (c for c in ['gene_symbols', 'gene_symbol', 'gene_name', 'symbol', 'SYMBOL']
            if c in adata_rna.var.columns),
            None
        )
        if symbol_col is not None:
            print(f"  Found symbol column: '{symbol_col}'. Mapping RNA row names to gene symbols...")
            id_to_symbol = dict(zip(
                adata_rna.var.index.astype(str),
                adata_rna.var[symbol_col].astype(str)
            ))
            traj_expression.index = [id_to_symbol.get(str(x), str(x)) for x in expr_index]
            expr_index = pd.Index(map(str, traj_expression.index))
            available_tfs = list(set(motif_index) & set(expr_index))
            
            if len(available_tfs) > 0:
                print(f"   After mapping to gene symbols, found {len(available_tfs)} TFs")

    if len(available_tfs) == 0:
        sample_motif = list(motif_index[:5])
        sample_expr = list(expr_index[:5])
        hint_cols = [c for c in ['gene_symbols','gene_symbol','gene_name','symbol','SYMBOL'] 
                    if c in adata_rna.var.columns]
        raise ValueError(
            "No TFs found with both activity and expression after name alignment.\n"
            f"Activity (motif) rows: {len(motif_index)} | Expression (RNA) rows: {len(expr_index)}\n"
            f"Sample motif names: {sample_motif}\n"
            f"Sample RNA row names: {sample_expr}\n"
            f"Gene symbol columns in RNA: {hint_cols if hint_cols else 'None found'}\n"
            "Tips: Ensure RNA uses gene symbols or provide a symbol column in adata_rna.var."
        )

    print(f"  Found {len(available_tfs)} TFs with both activity and expression data (after alignment)")

    filtered_tfs = available_tfs
    if len(filtered_tfs) == 0:
        raise ValueError("No TFs passed expression filters")
    
    # Calculate correlations with p-values (R equivalent: GetCorrelation)
    tf_correlations = {}
    skip_stats = {
        'missing_data': 0,
        'nan_inf': 0,
        'zero_variance': 0,
        'correlation_failed': 0
    }

    for tf in filtered_tfs:
        if tf not in traj_activity.index or tf not in traj_expression.index:
            skip_stats['missing_data'] += 1
            continue
        
        activity_traj = traj_activity.loc[tf].values
        expression_traj = traj_expression.loc[tf].values

        if np.any(~np.isfinite(activity_traj)) or np.any(~np.isfinite(expression_traj)):
            skip_stats['nan_inf'] += 1
            continue
        
        if np.var(activity_traj) == 0 or np.var(expression_traj) == 0:
            skip_stats['zero_variance'] += 1
            continue
        
        try:
            correlation, p_value = pearsonr(activity_traj, expression_traj)
            
            if not np.isnan(correlation) and not np.isnan(p_value):
                tf_correlations[tf] = {
                    'correlation': correlation,
                    'p_value': p_value
                }
            else:
                skip_stats['correlation_failed'] += 1
        except Exception as e:
            skip_stats['correlation_failed'] += 1
            continue

    print(f"  Calculated correlations for {len(tf_correlations)} TFs")
    print(f"  Skipped TFs breakdown:")
    print(f"    - Missing in trajectory data: {skip_stats['missing_data']}")
    print(f"    - Contains NaN/inf: {skip_stats['nan_inf']}")
    print(f"    - Zero variance: {skip_stats['zero_variance']}")
    print(f"    - Correlation failed: {skip_stats['correlation_failed']}")
    total_skipped = sum(skip_stats.values())
    print(f"    Total skipped: {total_skipped} / {len(filtered_tfs)}")
    
    # ========================================================================
    # Step 4: FDR correction and filtering (R: p.adjust + filter)
    # ========================================================================
    print(f"\nStep 4: Applying FDR correction and filtering...")
    print(f"  Thresholds: correlation > {correlation_threshold}, adj_p < {p_cutoff}")
    
    if len(tf_correlations) == 0:
        raise ValueError("No TF correlations calculated. Check your data and trajectory.")
    
    # Create DataFrame for easier manipulation
    tf_cor_df = pd.DataFrame([
        {
            'tf': tf_name,
            'correlation': tf_data['correlation'],
            'p_value': tf_data['p_value']
        }
        for tf_name, tf_data in tf_correlations.items()
    ])
    
    # FDR correction (R equivalent: p.adjust(p_value, method = "fdr"))
    try:
        from statsmodels.stats.multitest import multipletests
        _, adj_p_values, _, _ = multipletests(tf_cor_df['p_value'], method='fdr_bh')
        tf_cor_df['adj_p'] = adj_p_values
        print(f"   Applied FDR correction (Benjamini-Hochberg)")
    except ImportError:
        # Fallback: manual Benjamini-Hochberg procedure
        print("  Warning: statsmodels not available, using manual FDR correction")
        sorted_indices = np.argsort(tf_cor_df['p_value'])
        n = len(tf_cor_df)
        adj_p = np.zeros(n)
        
        for i, idx in enumerate(sorted_indices):
            rank = i + 1
            adj_p[idx] = min(1.0, tf_cor_df['p_value'].iloc[idx] * n / rank)
        
        tf_cor_df['adj_p'] = adj_p
        print(f"   Applied manual FDR correction")
    
    # Filter by correlation and p-value (R: df.cor <- df.cor[df.cor$adj_p < p.cutoff & df.cor$correlation > cor.cutoff, ])
    tf_cor_df['pass_adj_p'] = tf_cor_df['adj_p'] < p_cutoff
    tf_cor_df['pass_correlation'] = tf_cor_df['correlation'] > correlation_threshold
    tf_cor_df['pass_filters'] = tf_cor_df['pass_adj_p'] & tf_cor_df['pass_correlation']

    filtered_df = tf_cor_df[
        tf_cor_df['pass_filters']
    ].copy()
    
    print(f"   {len(filtered_df)} TFs passed both filters")
    print(f"    - adj_p < {p_cutoff}: {(tf_cor_df['adj_p'] < p_cutoff).sum()} TFs")
    print(f"    - correlation > {correlation_threshold}: {(tf_cor_df['correlation'] > correlation_threshold).sum()} TFs")
    
    if len(filtered_df) == 0:
        raise ValueError(
            "No TFs passed both selection thresholds: "
            f"correlation > {correlation_threshold} and adj_p < {p_cutoff}."
        )

    # Sort by correlation (descending). By default, retain every TF that
    # passes both filters, matching scMEGA's SelectTFs behavior.
    filtered_df = filtered_df.sort_values('correlation', ascending=False)
    if n_tfs is None:
        selected_tfs = filtered_df['tf'].tolist()
    else:
        if not isinstance(n_tfs, (int, np.integer)) or n_tfs <= 0:
            raise ValueError("n_tfs must be a positive integer or None")
        selected_tfs = filtered_df.head(int(n_tfs))['tf'].tolist()

    tf_cor_df_all = tf_cor_df.copy()
    tf_cor_df_all['selected'] = tf_cor_df_all['tf'].isin(selected_tfs)
    tf_cor_df_all['correlation_rank_all'] = (
        tf_cor_df_all['correlation'].rank(method='min', ascending=False).astype(int)
    )
    tf_cor_df_all['correlation_rank_passed'] = np.nan
    passed_rank = (
        filtered_df['correlation'].rank(method='min', ascending=False).astype(int)
    )
    tf_cor_df_all.loc[filtered_df.index, 'correlation_rank_passed'] = passed_rank
    tf_cor_df_all = tf_cor_df_all.sort_values(
        ['pass_filters', 'correlation'],
        ascending=[False, False]
    ).reset_index(drop=True)
    
    print(f"\n Selected {len(selected_tfs)} TFs")
    
    if len(selected_tfs) > 0:
        print(f"\nTop 5 TFs by correlation:")
        for i, row in filtered_df.head(5).iterrows():
            print(f"  {row.name+1}. {row['tf']}: r={row['correlation']:.3f}, adj_p={row['adj_p']:.2e}")
    
    print(f"\n TF selection completed using trajectory-based correlation method")
    
    # ========================================================================
    # R-COMPATIBLE: Calculate TF timepoints and prepare return structure
    # ========================================================================
    #                                           time_point = seq(1, 100, length.out = nrow(matMM)))
    
    if return_dict:
        print('\n  Calculating TF time points (based on peak activity time, R compatible mode)...')
        
        # R equivalent: apply(mat, 1, which.max)
        peak_positions = {}
        for tf in traj_activity.index:
            tf_trajectory = traj_activity.loc[tf].values
            if np.var(tf_trajectory) > 0:
                # R equivalent: which.max(tf_trajectory)
                peak_pos = np.argmax(tf_trajectory)
                peak_positions[tf] = peak_pos
            else:
                peak_positions[tf] = 0
        
        print(f"  The peak activity positions of {len(peak_positions)} TFs were calculated.")
        
        # R equivalent: idx <- order(apply(mat, 1, which.max))
        sorted_tfs = sorted(peak_positions.items(), key=lambda x: x[1])
        sorted_tf_names = [tf for tf, pos in sorted_tfs]
        
        # R: seq(1, 100, length.out = nrow(matMM))
        n_all_tfs = len(sorted_tf_names)
        all_time_points = np.linspace(1, 100, n_all_tfs)
        
        tf_to_timepoint = {tf: all_time_points[i] for i, tf in enumerate(sorted_tf_names)}
        
        print(f"   {n_all_tfs} TFs were assigned time points based on peak activity location")
        print(f"    Time point range: 1.0 - 100.0")
        
        # R: df_tf_time_point <- df_tf_time_point[df.cor$tfs,]
        tfs_data = []
        for tf in selected_tfs:
            time_point = tf_to_timepoint.get(tf, np.nan)
            
            tf_data = tf_correlations.get(tf, {'correlation': 0, 'p_value': 1})
            
            tf_row = filtered_df[filtered_df['tf'] == tf]
            adj_p = tf_row['adj_p'].values[0] if len(tf_row) > 0 else np.nan
            
            tfs_data.append({
                'tfs': tf,  # R column name: "tfs"
                'correlation': tf_data['correlation'],  # R column name: "correlation"
                'p_value': tf_data['p_value'],  # R column name: "p_value"
                'adj_p': adj_p,  # R column name: "adj_p"
                'time_point': time_point  # R column name: "time_point"
            })
        
        # Create a DataFrame that matches R’s df.cor structure
        tfs_df = pd.DataFrame(tfs_data)
        
        # R: df.cor <- df.cor[order(df.cor$time_point),]
        tfs_df = tfs_df.sort_values('time_point').reset_index(drop=True)
        
        print(f"   Time points extracted for {len(selected_tfs)} selected TFs")
        print(f"    The time point range of the selected TF: {tfs_df['time_point'].min():.1f} - {tfs_df['time_point'].max():.1f}")
        print(f"   Created an R-compatible DataFrame with column name: {list(tfs_df.columns)}")


        heatmap_fig = None
        if return_heatmap:
            print("\n  Generating TF selection heatmap...")
        
            try:
                from python_scmega.visualization.heatmaps import tf_heatmap
                
                selected_tf_names = tfs_df['tfs'].tolist()
                
                if all(tf in traj_activity.index for tf in selected_tf_names[:min(3, len(selected_tf_names))]):
                    trajMM_selected = traj_activity.loc[selected_tf_names]
                    trajRNA_selected = traj_expression.loc[selected_tf_names]
                else:
                    # If indexes do not match, skip heatmap generation
                    print(f"   Index mismatch: traj_activity uses a motif ID but requires a TF name")
                    print(f"   The first 3 indexes of traj_activity: {list(traj_activity.index[:3])}")
                    print(f"   The first 3 required TF names: {selected_tf_names[:3]}")
                    raise ValueError("Index mismatch between trajectory data and selected TFs")
                
                heatmap_fig = tf_heatmap(
                    tf_activity=trajMM_selected,
                    tf_expression=trajRNA_selected,
                    scale_rows=True,
                    limits=(-2, 2),
                    label_top=min(30, len(selected_tfs)),
                    name1="TF activity",
                    name2="TF expression"
                )
                print(f"   Heatmap generated for {len(selected_tfs)} TFs")
            except Exception as e:
                print(f"   Heatmap generation failed: {e}")
                heatmap_fig = None
                
        # R-compatible return structure matching select_tf_gene.R line 87-94
        # R returns: list("tfs" = df.cor, "heatmap" = ht)
        result = {
            'tfs': tfs_df,      # Main result, matching R's "tfs" key
            'tf_cor_df': tf_cor_df_all,
            'heatmap': heatmap_fig      # Placeholder for heatmap (can be implemented later)
        }
        
        print(f"\nReturning R-compatible dictionary:")
        print(f"  - 'tfs': DataFrame with {len(tfs_df)} TFs")
        print(f"  - 'tf_cor_df': DataFrame with {len(tf_cor_df_all)} TFs before filtering")
        print(f"    Columns: {list(tfs_df.columns)}")
        
        return result
    else:
        # Legacy behavior: return only TF list
        return selected_tfs
