"""
Visualization Module

This module provides comprehensive visualization tools for single-cell multiome
data analysis, including network visualization, trajectory plots, heatmaps.

Key Components:
- network_plots: Gene regulatory network visualization (GRNPlot, NetCentPlot equivalents)
- trajectory_plots: Trajectory and pseudotime visualization
- heatmaps: Expression and correlation heatmaps


Main Functions:
- plot_grn_network: Main GRN network visualization function
- plot_network_centrality: Network centrality analysis plots
- plot_grn_heatmap: GRN correlation heatmaps
- plot_trajectory_heatmap: Trajectory expression heatmaps

Example Usage:
```python
import python_scmega as pymega
import matplotlib.pyplot as plt

# Network visualization
fig = pymega.plot_grn_network(
    grn_network,
    layout="spring",
    node_size_col="centrality",
    edge_width_col="confidence"
)

# Trajectory heatmap
fig = pymega.plot_trajectory_heatmap(
    trajectory_data,
    pseudotime_col="pseudotime",
    cluster_genes=True
)


```
"""

from .network_plots import (
    plot_network_centrality,
    plot_network_topology,
    plot_tf_targets,
    create_network_layout
)



from .trajectory_plots import (
    plot_trajectory
)

from .heatmaps import (
    tf_gene_correlation_heatmap,
    tf_heatmap,
    gene_heatmap
)


__all__ = [
    # Network plots
    'plot_network_centrality',
    'plot_network_topology',
    'plot_tf_targets',
    'create_network_layout',
    
    # Trajectory plots
    'plot_trajectory',
    
    # Heatmaps
    'tf_gene_correlation_heatmap',
    'gene_heatmap',
    'tf_heatmap',
    
    
]