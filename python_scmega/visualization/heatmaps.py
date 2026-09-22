"""
Heatmap Visualization Module

Implements comprehensive heatmap visualization functions equivalent to scMEGA's
GRNHeatmap, TrajectoryHeatmap, and CorrelationHeatmap functions, providing
high-quality heatmaps for expression data, correlations, and trajectory analysis.

Key functions:
- tf_gene_correlation_heatmap: TF-Gene correlation heatmaps (GRNHeatmap equivalent)
- tf_heatmap: TF activity and expression heatmaps (CorrelationHeatmap equivalent)
- gene_heatmap: Chromatin accessibility and gene expression heatmaps

References:
- Original R code: GRNHeatmap, TrajectoryHeatmap, CorrelationHeatmap functions in scMEGA
- Seaborn documentation for advanced heatmap visualization
- Matplotlib documentation for custom heatmap styling
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, Dict, Any, List, Union, Tuple
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.cluster import KMeans
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable


def tf_gene_correlation_heatmap(
    tf_gene_cor: pd.DataFrame,
    tf_timepoint: Optional[Dict[str, float]] = None,
    n_clusters: int = 1,
    figsize: Tuple[int, int] = (12, 10),
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    show_row_names: bool = False,
    show_column_names: bool = True,
    cluster_method: str = "ward",
    cmap: str = "RdBu_r",
    **kwargs
) -> plt.Figure:
    """
    Create a TF-gene correlation heatmap (corresponding to the R version of GRNHeatmap)
    
    This function draws a correlation heatmap between TFs and target genes, showing the patterns in the gene regulatory network.
    Correlation pattern. Annotations can be made based on the time points of TF, and hierarchical clustering and k-means clustering are supported.
    
    Args:
        tf_gene_cor: TF-gene correlation data frame, including 'tf', 'gene', 'correlation' columns
        tf_timepoint: TF time point dictionary, the key is the TF name and the value is the time point
        n_clusters: Number of clusters for k-means clustering (corresponding to the row_km parameter of R)
        figsize: graphic size
        title: graphic title
        save_path: save path
        show_row_names: whether to display row names (gene names)
        show_column_names: whether to display column names (TF names)
        cluster_method: clustering method ('ward' corresponds to R's 'ward.D2')
        cmap: Heatmap color scheme
        **kwargs: other parameters
        
    Returns:
        matplotlib.Figure: TF-gene correlation heatmap
        
    Example:
        >>> fig = tf_gene_correlation_heatmap(
        ...     tf_gene_cor=df_cor,
        ...     tf_timepoint={'TF1': 10, 'TF2': 50, 'TF3': 90},
        ...     n_clusters=3
        ... )
    """
    if not all(col in tf_gene_cor.columns for col in ['tf', 'gene', 'correlation']):
        raise ValueError("tf_gene_cor must contain 'tf', 'gene', 'correlation' columns")
    
    mat_cor = tf_gene_cor.pivot_table(
        index='gene',
        columns='tf',
        values='correlation',
        aggfunc='first'  # If there are duplicates, take the first value (aligned R version)
    )
    print(f"  [DEBUG] mat_cor shape: {mat_cor.shape}")
    print(f"  [DEBUG] Number of valid values of mat_cor: {mat_cor.notna().sum().sum()}")
    print(f"  [DEBUG] mat_cor value range: [{mat_cor.min().min():.3f}, {mat_cor.max().max():.3f}]")
    # Handle missing values - fill only when necessary
    if mat_cor.isna().any().any():
        print(f"  [WARNING] {mat_cor.isna().sum().sum()} NaN values detected, padded with 0s")
        mat_cor = mat_cor.fillna(0)
    if mat_cor.empty:
        raise ValueError('Correlation matrix is empty and heatmap cannot be drawn')
    if mat_cor.nunique().sum() <= 1:
        raise ValueError(f"Correlation matrix data anomaly: all values are the same (value={mat_cor.iloc[0,0] if not mat_cor.empty else 'N/A'})")
    
    if mat_cor.empty:
        raise ValueError('Correlation matrix is empty and heatmap cannot be drawn')
    
    # 2. If time points are provided, arrange the TF columns in time point order (corresponding to cluster_columns=FALSE in R)
    col_colors = None
    if tf_timepoint is not None:
        # Only keep TF with time point information
        tf_with_time = [tf for tf in mat_cor.columns if tf in tf_timepoint]
        if len(tf_with_time) > 0:
            sorted_tfs = sorted(tf_with_time, key=lambda x: tf_timepoint[x])
            mat_cor = mat_cor[sorted_tfs]
        
        timepoints = [tf_timepoint[tf] for tf in mat_cor.columns]
        
        colors_list = ['#2166ac', '#4393c3', '#92c5de', '#d1e5f0', 
                       '#fefebe', '#fddbc7', '#f4a582', '#d6604d', '#b2182b']
        cmap_time = LinearSegmentedColormap.from_list('blueYellow', colors_list)
        
        norm = Normalize(vmin=min(timepoints), vmax=max(timepoints))
        sm = ScalarMappable(norm=norm, cmap=cmap_time)
        
        col_colors = pd.Series(
            [sm.to_rgba(tp) for tp in timepoints],
            index=mat_cor.columns
        )
    
    row_cluster = True
    if n_clusters > 1:
        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        row_clusters = kmeans.fit_predict(mat_cor.values)
        
        # ComplexHeatmap will perform hierarchical clustering after k-means, which is simplified to sorting by clusters.
        cluster_order = np.argsort(row_clusters)
        mat_cor = mat_cor.iloc[cluster_order, :]
        row_cluster = False  # Already sorted, no need to cluster again
    
    fig = plt.figure(figsize=figsize)
    
    if col_colors is not None:
        g = sns.clustermap(
            mat_cor,
            method=cluster_method,
            cmap=cmap,
            center=0,
            col_cluster=False,  # R version: cluster_columns = FALSE
            row_cluster=row_cluster,
            linewidths=0.5,
            cbar_kws={'label': 'Correlation'},
            figsize=figsize,
            col_colors=col_colors,
            yticklabels=show_row_names,
            xticklabels=show_column_names,
            **kwargs
        )
    else:
        g = sns.clustermap(
            mat_cor,
            method=cluster_method,
            cmap=cmap,
            center=0,
            col_cluster=False,  # R version default: cluster_columns = FALSE
            row_cluster=row_cluster,
            linewidths=0.5,
            cbar_kws={'label': 'Correlation'},
            figsize=figsize,
            yticklabels=show_row_names,
            xticklabels=show_column_names,
            **kwargs
        )
    
    if title:
        g.fig.suptitle(title, fontsize=14, fontweight='bold', y=0.98)
    else:
        g.fig.suptitle('TF-Gene Correlation Heatmap', fontsize=14, fontweight='bold', y=0.98)
    
    g.ax_heatmap.set_xlabel('Transcription Factors', fontsize=12)
    g.ax_heatmap.set_ylabel('Genes', fontsize=12)
    
    plt.tight_layout()
    
    if save_path:
        g.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return g.fig


def tf_heatmap(
    tf_activity: pd.DataFrame,
    tf_expression: pd.DataFrame,
    scale_rows: bool = True,
    limits: Tuple[float, float] = (-2, 2),
    label_rows: bool = False,
    label_top: int = 50,
    figsize: Tuple[int, int] = (16, 10),
    name1: str = "TF Activity",
    name2: str = "TF Expression",
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    **kwargs
) -> plt.Figure:
    """
    Create a combined heatmap of TF activity and TF expression (corresponding to the R version of CorrelationHeatmap)
    
    This function creates two side-by-side heatmaps showing changes in TF activity and TF expression along pseudo-time.
    Both heatmaps share the same row order, determined by merging the normalized values of the two matrices.
    
    Args:
        tf_activity: TF activity matrix (TF x pseudo time bins)
        tf_expression: TF expression matrix (TF x pseudo-time bins)
        scale_rows: whether to z-score normalize rows
        limits: color range limits
        label_rows: whether to display all row labels
        label_top: Display the labels of the top N highly variable genes
        figsize: graphic size
        name1: the name of the first heat map
        name2: the name of the second heatmap
        title: overall title
        save_path: save path
        **kwargs: other parameters
        
    Returns:
        matplotlib.Figure: combined heatmap
        
    Example:
        >>> fig = tf_heatmap(
        ...     tf_activity=activity_matrix,
        ...     tf_expression=expression_matrix,
        ...     label_top=30
        ... )
    """
    # 1. Make sure both matrices have the same TF
    common_tfs = tf_activity.index.intersection(tf_expression.index)
    if len(common_tfs) == 0:
        raise ValueError('TF activity and TF expression matrix have no TFs in common')
    
    tf_activity = tf_activity.loc[common_tfs]
    tf_expression = tf_expression.loc[common_tfs]
    
    if scale_rows:
        scaler = StandardScaler()
        activity_scaled = pd.DataFrame(
            scaler.fit_transform(tf_activity.T).T,
            index=tf_activity.index,
            columns=tf_activity.columns
        )
        expression_scaled = pd.DataFrame(
            scaler.fit_transform(tf_expression.T).T,
            index=tf_expression.index,
            columns=tf_expression.columns
        )
    else:
        activity_scaled = tf_activity.copy()
        expression_scaled = tf_expression.copy()
    
    activity_scaled = activity_scaled.clip(lower=limits[0], upper=limits[1])
    expression_scaled = expression_scaled.clip(lower=limits[0], upper=limits[1])
    
    # 3. Determine the row order: merge two normalized matrices (corresponding to the core logic of the R version)
    combined_mat = activity_scaled + expression_scaled
    
    row_vars = combined_mat.var(axis=1)
    
    max_positions = combined_mat.idxmax(axis=1)
    col_positions = {col: i for i, col in enumerate(combined_mat.columns)}
    max_pos_numeric = max_positions.map(col_positions)
    
    row_order = max_pos_numeric.sort_values().index.tolist()
    
    activity_ordered = activity_scaled.loc[row_order]
    expression_ordered = expression_scaled.loc[row_order]
    
    if label_top > 0 and not label_rows:
        top_var_tfs = row_vars.nlargest(label_top).index.tolist()
        label_indices = [i for i, tf in enumerate(row_order) if tf in top_var_tfs]
    else:
        label_indices = None
    
    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)
    
    colors_solar = ['#3d155f', '#5a187b', '#801f95', '#a82ba9', '#cf3db9', 
                    '#f25cb4', '#fc82a8', '#ffa598', '#ffc58b', '#ffe085']
    cmap_solar = LinearSegmentedColormap.from_list('solarExtra', colors_solar)
    
    colors_horizon = ['#000436', '#011271', '#021e90', '#042aa1', '#0636b0', 
                      '#0942bd', '#0c4fc8', '#105dd2', '#156bdb', '#1a79e3']
    cmap_horizon = LinearSegmentedColormap.from_list('horizonExtra', colors_horizon)
    
    im1 = axes[0].imshow(
        activity_ordered.values,
        aspect='auto',
        cmap=cmap_solar,
        vmin=limits[0],
        vmax=limits[1],
        interpolation='nearest'
    )
    axes[0].set_title(name1, fontsize=12, fontweight='bold')
    axes[0].set_xlabel('Pseudotime', fontsize=10)
    axes[0].set_ylabel('Transcription Factors', fontsize=10)
    
    n_bins = len(activity_ordered.columns)
    n_ticks = min(10, n_bins)
    tick_positions = np.linspace(0, n_bins-1, n_ticks, dtype=int)
    axes[0].set_xticks(tick_positions)
    axes[0].set_xticklabels([f'{int(i*100//(n_bins-1))}%' for i in tick_positions])
    
    cbar1 = plt.colorbar(im1, ax=axes[0])
    cbar1.set_label(name1, rotation=270, labelpad=15)
    
    im2 = axes[1].imshow(
        expression_ordered.values,
        aspect='auto',
        cmap=cmap_horizon,
        vmin=limits[0],
        vmax=limits[1],
        interpolation='nearest'
    )
    axes[1].set_title(name2, fontsize=12, fontweight='bold')
    axes[1].set_xlabel('Pseudotime', fontsize=10)
    
    axes[1].set_xticks(tick_positions)
    axes[1].set_xticklabels([f'{int(i*100//(n_bins-1))}%' for i in tick_positions])
    
    cbar2 = plt.colorbar(im2, ax=axes[1])
    cbar2.set_label(name2, rotation=270, labelpad=15)
    
    if label_rows:
        axes[0].set_yticks(range(len(row_order)))
        axes[0].set_yticklabels(row_order, fontsize=6)
    elif label_indices is not None and len(label_indices) > 0:
        axes[0].set_yticks(label_indices)
        axes[0].set_yticklabels([row_order[i] for i in label_indices], fontsize=8)
    else:
        axes[0].set_yticks([])
    
    if title:
        fig.suptitle(title, fontsize=14, fontweight='bold', y=0.98)
    else:
        fig.suptitle('TF Activity and Expression Along Trajectory', 
                     fontsize=14, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig


def gene_heatmap(
    chromatin_accessibility: pd.DataFrame,
    gene_expression: pd.DataFrame,
    scale_rows: bool = True,
    limits: Tuple[float, float] = (-2, 2),
    label_rows: bool = False,
    label_top_accessibility: int = 10,
    label_top_expression: int = 10,
    figsize: Tuple[int, int] = (16, 10),
    name1: str = "Chromatin Accessibility",
    name2: str = "Gene Expression",
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    **kwargs
) -> plt.Figure:
    """
    Create a combined heatmap of chromatin accessibility and gene expression (correlationHeatmap in the R version of SelectGenes)
    
     Important: This function assumes that chromatin_accessibility and gene_expression are in one-to-one correspondence!
    That is, the i-th peak corresponds to the i-th gene (from the same peak-to-gene connection).
    Both heatmaps use the same row order to show this correspondence.
    
    Args:
        chromatin_accessibility: chromatin accessibility matrix (peaks x pseudo-time bins)
        gene_expression: gene expression matrix (genes x pseudo-time bins)
        scale_rows: whether to z-score normalize rows
        limits: color range limits
        label_rows: whether to display all row labels
        label_top_accessibility: The number of top labels displayed by the chromatin accessibility heat map
        label_top_expression: The number of top labels displayed in the gene expression heat map
        figsize: graphic size
        name1: the name of the first heat map
        name2: the name of the second heatmap
        title: overall title
        save_path: save path
        **kwargs: other parameters
        
    Returns:
        matplotlib.Figure: combined heatmap
        
    Example:
        >>> # Note: peak and gene must correspond one to one!
        >>> peak_matrix = trajATAC.loc[df_p2g['peak']]
        >>> gene_matrix = trajRNA.loc[df_p2g['gene']]
        >>> fig = gene_heatmap(
        ...     chromatin_accessibility=peak_matrix,
        ...     gene_expression=gene_matrix,
        ...     label_top_accessibility=10,
        ...     label_top_expression=10
        ... )
    """
    # 0. Make sure the number of rows is the same (one-to-one correspondence)
    if len(chromatin_accessibility) != len(gene_expression):
        raise ValueError(
            f"Peak and Gene numbers do not match! Peak: {len(chromatin_accessibility)}, "
            f"Gene: {len(gene_expression)}. "
            'These two matrices should come from peak-to-gene connections and correspond one to one.'
        )
    
    if scale_rows:
        scaler = StandardScaler()
        accessibility_scaled = pd.DataFrame(
            scaler.fit_transform(chromatin_accessibility.T).T,
            index=chromatin_accessibility.index,
            columns=chromatin_accessibility.columns
        )
        expression_scaled = pd.DataFrame(
            scaler.fit_transform(gene_expression.T).T,
            index=gene_expression.index,
            columns=gene_expression.columns
        )
    else:
        accessibility_scaled = chromatin_accessibility.copy()
        expression_scaled = gene_expression.copy()
    
    accessibility_scaled = accessibility_scaled.clip(lower=limits[0], upper=limits[1])
    expression_scaled = expression_scaled.clip(lower=limits[0], upper=limits[1])
    
    # 2. Key: Use the merge matrix to determine a unified row order (corresponding to the R version of CorrelationHeatmap logic)
    combined_mat = accessibility_scaled.values + expression_scaled.values
    combined_df = pd.DataFrame(
        combined_mat,
        index=accessibility_scaled.index,
        columns=accessibility_scaled.columns
    )
    
    max_positions = combined_df.idxmax(axis=1)
    col_positions = {col: i for i, col in enumerate(combined_df.columns)}
    max_pos_numeric = max_positions.map(col_positions)
    
    row_order_idx = max_pos_numeric.sort_values().index
    
    # accessibility_ordered = accessibility_scaled.loc[row_order_idx]
    # expression_ordered = expression_scaled.loc[row_order_idx]
    peaks_in_order = [idx for idx in row_order_idx if idx in accessibility_scaled.index]
    accessibility_ordered = accessibility_scaled.loc[peaks_in_order]
    genes_in_order = [idx for idx in row_order_idx if idx in expression_scaled.index]
    expression_ordered = expression_scaled.loc[genes_in_order]

    if label_top_accessibility > 0 and not label_rows:
        row_vars_acc = accessibility_scaled.var(axis=1)
        top_var_peaks = row_vars_acc.nlargest(label_top_accessibility).index.tolist()
        label_indices_acc = [i for i, pk in enumerate(row_order_idx) if pk in top_var_peaks]
    else:
        label_indices_acc = None
    
    if label_top_expression > 0 and not label_rows:
        row_vars_exp = expression_scaled.var(axis=1)
        top_var_genes = row_vars_exp.nlargest(label_top_expression).index.tolist()
        label_indices_exp = [i for i, gene in enumerate(row_order_idx) if gene in top_var_genes]
    else:
        label_indices_exp = None
    
    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)
    
    colors_solar = ['#3d155f', '#5a187b', '#801f95', '#a82ba9', '#cf3db9', 
                    '#f25cb4', '#fc82a8', '#ffa598', '#ffc58b', '#ffe085']
    cmap_solar = LinearSegmentedColormap.from_list('solarExtra', colors_solar)
    
    colors_horizon = ['#000436', '#011271', '#021e90', '#042aa1', '#0636b0', 
                      '#0942bd', '#0c4fc8', '#105dd2', '#156bdb', '#1a79e3']
    cmap_horizon = LinearSegmentedColormap.from_list('horizonExtra', colors_horizon)
    
    im1 = axes[0].imshow(
        accessibility_ordered.values,
        aspect='auto',
        cmap=cmap_solar,
        vmin=limits[0],
        vmax=limits[1],
        interpolation='nearest'
    )
    axes[0].set_title(name1, fontsize=12, fontweight='bold')
    axes[0].set_xlabel('Pseudotime', fontsize=10)
    axes[0].set_ylabel('Peaks', fontsize=10)
    
    n_bins = len(accessibility_ordered.columns)
    n_ticks = min(10, n_bins)
    tick_positions = np.linspace(0, n_bins-1, n_ticks, dtype=int)
    axes[0].set_xticks(tick_positions)
    axes[0].set_xticklabels([f'{int(i*100//(n_bins-1))}%' for i in tick_positions])
    
    cbar1 = plt.colorbar(im1, ax=axes[0])
    cbar1.set_label(name1, rotation=270, labelpad=15)
    
    if label_rows:
        axes[0].set_yticks(range(len(row_order_idx)))
        axes[0].set_yticklabels(row_order_idx.tolist(), fontsize=4)
    elif label_indices_acc is not None and len(label_indices_acc) > 0:
        axes[0].set_yticks(label_indices_acc)
        axes[0].set_yticklabels([row_order_idx[i] for i in label_indices_acc], fontsize=6)
    else:
        axes[0].set_yticks([])
    
    im2 = axes[1].imshow(
        expression_ordered.values,
        aspect='auto',
        cmap=cmap_horizon,
        vmin=limits[0],
        vmax=limits[1],
        interpolation='nearest'
    )
    axes[1].set_title(name2, fontsize=12, fontweight='bold')
    axes[1].set_xlabel('Pseudotime', fontsize=10)
    axes[1].set_ylabel('Genes', fontsize=10)
    
    axes[1].set_xticks(tick_positions)
    axes[1].set_xticklabels([f'{int(i*100//(n_bins-1))}%' for i in tick_positions])
    
    cbar2 = plt.colorbar(im2, ax=axes[1])
    cbar2.set_label(name2, rotation=270, labelpad=15)
    
    if label_rows:
        axes[1].set_yticks(range(len(row_order_idx)))
        axes[1].set_yticklabels(expression_ordered.index.tolist(), fontsize=4)
    elif label_indices_exp is not None and len(label_indices_exp) > 0:
        axes[1].set_yticks(label_indices_exp)
        axes[1].set_yticklabels([expression_ordered.index[i] for i in label_indices_exp], fontsize=6)
    else:
        axes[1].set_yticks([])
    
    if title:
        fig.suptitle(title, fontsize=14, fontweight='bold', y=0.98)
    else:
        fig.suptitle('Chromatin Accessibility and Gene Expression Along Trajectory',
                     fontsize=14, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig
