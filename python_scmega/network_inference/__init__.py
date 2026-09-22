"""
Network Inference Module

This module provides comprehensive tools for gene regulatory network (GRN)
inference from multiome data, integrating TF-gene correlations, motif matching,
and peak-gene links to construct and validate regulatory networks.

Key Components:
- grn_inference: Core GRN inference algorithms (GetGRN, GetTFGeneCorrelation equivalents)
- network_construction: Network topology construction and integration

Main Functions:
- infer_grn: Main GRN inference function
- get_tf_gene_correlation: TF-gene correlation analysis
- construct_grn_network: Network construction from multiple data sources

Example Usage:
```python
from python_scmega.network_inference import (
    infer_grn, 
    get_tf_gene_correlation,
    plot_grn_network
)

# Step 1: Calculate TF-gene correlations
tf_gene_cor = get_tf_gene_correlation(
    data=multiome,
    tf_list=tf_list,
    gene_list=gene_list,
    method="pearson",
    trajectory_name="Trajectory",
    motif_names=motif_names
)

# Step 2: Infer GRN network
grn_network = infer_grn(
    motif_matching=motif_matrix,
    tf_gene_cor=tf_gene_cor,
    peak_gene_links=peak_gene_links,
    correlation_threshold=0.5
)

# Step 3: Visualize network
fig = plot_grn_network(
    grn_network=grn_network,
    tfs_timepoint=tfs_timepoint,
    show_tf_labels=True,
    seed=42
)
```
"""

from .grn_inference import (
    infer_grn,
    get_tf_gene_correlation
)


# Network visualization functions
from .network_construction import (
    plot_grn_network,
    create_tf_timepoint_dict,
    plot_grn_network_bipartite,
    plot_grn_network_circular
)


__all__ = [
    # GRN inference (grn_inference.py)
    'infer_grn',
    'get_tf_gene_correlation',
    
    # Network visualization (network_construction.py)
    'plot_grn_network',
    'create_tf_timepoint_dict',
    'plot_grn_network_bipartite',
    'plot_grn_network_circular',
]
