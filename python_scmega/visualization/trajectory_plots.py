"""
Trajectory Visualization Module

Implements comprehensive trajectory visualization functions equivalent to scMEGA's
TrajectoryPlot and PseudotimePlot functions, providing high-quality visualization
for pseudotime trajectories and developmental paths.

Key functions:
- plot_trajectory: Main trajectory visualization (R TrajectoryPlot equivalent)
- plot_pseudotime: TF dynamics plot (R PseudotimePlot equivalent)
- plot_trajectory_genes: Gene expression along trajectories
- plot_trajectory_heatmap: Heatmap visualization

References:
- Original R code: TrajectoryPlot, PseudotimePlot functions in scMEGA
- ArchR documentation for trajectory visualization
- scMEGA R package: https://github.com/CostaLab/scMEGA
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Optional, Dict, Any, List, Union, Tuple, TYPE_CHECKING
import warnings
from scipy.signal import savgol_filter

# Import from trajectory_analysis module
from ..trajectory_analysis.pseudotime_analysis import get_trajectory_data

# MuData support
try:
    import muon as mu
    MUON_AVAILABLE = True
except ImportError:
    mu = None
    MUON_AVAILABLE = False

try:
    import anndata as ad
except ImportError:
    ad = None


def plot_trajectory(
    data: Union['ad.AnnData', 'mu.MuData'],
    trajectory: str = "Trajectory",
    reduction: str = "X_umap",
    assay: str = "rna",
    size: float = 0.2,
    quantile_cut: Optional[Tuple[float, float]] = (0.01, 0.99),
    continuous_set: str = "viridis",
    discrete_set: Optional[str] = None,
    randomize: bool = True,
    keep_axis: bool = False,
    add_arrow: bool = False,
    smooth_window: int = 5,
    figsize: Tuple[int, int] = (8, 8),
    title: Optional[str] = None,
    rastr: bool = False,
    base_size: int = 6,
    save_path: Optional[str] = None,
    **kwargs
) -> plt.Figure:
    """
    Visualize cell trajectories in embedding space.
    
    Python implementation of R's TrajectoryPlot function from scMEGA package.
    
    R equivalent:
    ```R
    TrajectoryPlot(
        object = obj,
        trajectory = "Trajectory",
        reduction = "umap",
        size = 0.2,
        quantCut = c(0.01, 0.99),
        continuousSet = "horizonExtra",
        discreteSet = "stallion",
        randomize = TRUE,
        keepAxis = FALSE,
        addArrow = FALSE,
        smoothWindow = 5
    )
    ```
    
    Args:
        data: MuData or AnnData object with trajectory in obs
        trajectory: Column name in obs for trajectory values (pseudotime)
        reduction: Key in obsm for embedding coordinates (e.g., "X_umap", "X_pca")
        assay: Which modality to use for MuData ("rna", "atac")
        size: Point size (default 0.2 to match R)
        quantile_cut: Quantile cutoffs for color scaling to handle outliers
        continuous_set: Color palette name for continuous values
        discrete_set: Color palette name for discrete values
        randomize: Whether to randomize point plotting order
        keep_axis: Whether to keep axis ticks and labels
        add_arrow: Whether to add trajectory arrow overlay
        smooth_window: Smoothing window for arrow (R default: 5)
        figsize: Figure size
        title: Plot title
        rastr: Whether to rasterize points (for large datasets)
        base_size: Base font size
        save_path: Path to save the plot
        **kwargs: Additional parameters
        
    Returns:
        matplotlib.Figure: Trajectory visualization
        
    Example:
        ```python
        # Basic usage
        fig = plot_trajectory(
            mdata,
            trajectory="Trajectory",
            reduction="X_umap",
            add_arrow=True
        )
        
        # With custom styling
        fig = plot_trajectory(
            mdata,
            trajectory="Trajectory",
            reduction="X_umap",
            size=0.5,
            quantile_cut=(0.025, 0.975),
            continuous_set="viridis",
            add_arrow=True,
            smooth_window=7
        )
        ```
    """
    # Extract data from MuData/AnnData
    if MUON_AVAILABLE and isinstance(data, mu.MuData):
        if assay not in data.mod:
            raise ValueError(f"Assay '{assay}' not found in MuData. Available: {list(data.mod.keys())}")
        adata = data.mod[assay]
    else:
        adata = data
    
    # Get trajectory values
    if trajectory not in adata.obs.columns:
        raise ValueError(f"Trajectory '{trajectory}' not found in obs columns")
    
    pseudotime = adata.obs[trajectory].values.copy()
    
    # Get embedding coordinates
    if reduction not in adata.obsm.keys():
        raise ValueError(f"Reduction '{reduction}' not found in obsm. Available: {list(adata.obsm.keys())}")
    
    embedding = adata.obsm[reduction][:, :2].copy()
    
    # Filter NA values
    valid_mask = ~np.isnan(pseudotime)
    pseudotime_valid = pseudotime[valid_mask]
    embedding_valid = embedding[valid_mask]
    
    if len(pseudotime_valid) == 0:
        raise ValueError("No valid trajectory values found (all NaN)")
    
    print(f"Plotting {len(pseudotime_valid)} cells with valid trajectory values")
    
    # Apply quantile cutting to color values (like R's quantCut)
    if quantile_cut is not None:
        lower_q = np.quantile(pseudotime_valid, quantile_cut[0])
        upper_q = np.quantile(pseudotime_valid, quantile_cut[1])
        color_values = np.clip(pseudotime_valid, lower_q, upper_q)
        print(f"Applied quantile cut: [{lower_q:.2f}, {upper_q:.2f}]")
    else:
        color_values = pseudotime_valid
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Randomize plotting order if requested (prevents clustering artifacts)
    if randomize:
        random_idx = np.random.permutation(len(embedding_valid))
        embedding_plot = embedding_valid[random_idx]
        color_values_plot = color_values[random_idx]
        pseudotime_plot = pseudotime_valid[random_idx]
    else:
        embedding_plot = embedding_valid
        color_values_plot = color_values
        pseudotime_plot = pseudotime_valid
    
    # Plot cells
    scatter = ax.scatter(
        embedding_plot[:, 0], 
        embedding_plot[:, 1],
        c=color_values_plot,
        s=size,
        cmap=continuous_set,
        alpha=0.8,
        rasterized=rastr
    )
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Pseudotime', fontsize=base_size + 2)
    
    # Add arrow overlay if requested (matches R's addArrow=TRUE)
    if add_arrow:
        _add_trajectory_arrow(
            ax, 
            embedding_valid, 
            pseudotime_valid, 
            smooth_window=smooth_window
        )
    
    # Set labels
    reduction_name = reduction.replace("X_", "").upper()
    ax.set_xlabel(f'{reduction_name} 1', fontsize=base_size + 2)
    ax.set_ylabel(f'{reduction_name} 2', fontsize=base_size + 2)
    ax.set_title(title or f'Trajectory: {trajectory}', fontsize=base_size + 4)
    
    # Remove axis ticks if keep_axis=False (match R default)
    if not keep_axis:
        ax.set_xticks([])
        ax.set_yticks([])
    else:
        ax.tick_params(labelsize=base_size)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved to {save_path}")
    
    return fig


def _add_trajectory_arrow(ax, embedding, pseudotime, smooth_window=5):
    """
    Add trajectory arrow overlay to plot.
    
    This implements R's addArrow=TRUE logic:
    1. Bin cells by pseudotime (floor(pseudotime / 1.01))
    2. Calculate mean coordinates for each bin
    3. Apply rolling mean smoothing (centerRollMean)
    4. Draw path with arrow
    
    Args:
        ax: Matplotlib axis
        embedding: Cell coordinates (n_cells, 2)
        pseudotime: Pseudotime values (n_cells,)
        smooth_window: Smoothing window size
    """
    # Sort by pseudotime
    sort_idx = np.argsort(pseudotime)
    sorted_embedding = embedding[sort_idx]
    sorted_pseudotime = pseudotime[sort_idx]
    
    # Bin pseudotime and calculate means (matches R logic)
    # R: split(dfT, floor(dfT$PseudoTime / 1.01))
    bins = np.floor(sorted_pseudotime / 1.01).astype(int)
    unique_bins = np.unique(bins)
    
    arrow_coords = []
    for bin_val in unique_bins:
        bin_mask = bins == bin_val
        if bin_mask.sum() > 0:
            mean_coords = sorted_embedding[bin_mask].mean(axis=0)
            arrow_coords.append(mean_coords)
    
    if len(arrow_coords) < 2:
        warnings.warn("Not enough trajectory points for arrow, skipping")
        return
    
    arrow_coords = np.array(arrow_coords)
    
    # Apply smoothing (rolling mean) - matches R's .centerRollMean
    if smooth_window > 0 and len(arrow_coords) >= smooth_window:
        smoothed_x = _center_roll_mean(arrow_coords[:, 0], smooth_window)
        smoothed_y = _center_roll_mean(arrow_coords[:, 1], smooth_window)
        arrow_coords = np.column_stack([smoothed_x, smoothed_y])
    
    # Plot arrow path
    ax.plot(
        arrow_coords[:, 0], 
        arrow_coords[:, 1],
        color='black',
        linewidth=2,
        alpha=0.8,
        zorder=100
    )
    
    # Add arrowhead at the end (matches R's arrow(type="open", length=0.1))
    if len(arrow_coords) > 1:
        ax.annotate(
            '',
            xy=arrow_coords[-1],
            xytext=arrow_coords[-2],
            arrowprops=dict(
                arrowstyle='->',
                lw=2,
                color='black',
                alpha=0.8
            ),
            zorder=101
        )
    
    print(f"Added trajectory arrow with {len(arrow_coords)} points")


def _center_roll_mean(x, k):
    """
    Center-aligned rolling mean.
    
    Matches R's .centerRollMean function from ArchR.
    
    Args:
        x: Input array
        k: Window size
        
    Returns:
        Smoothed array
    """
    if k <= 1:
        return x
    
    n = len(x)
    result = np.zeros(n)
    
    for i in range(n):
        # Center window around position i
        left = max(0, i - k // 2)
        right = min(n, i + k // 2 + 1)
        result[i] = np.mean(x[left:right])
    
    return result




# Export main functions
__all__ = [
    'plot_trajectory',
]