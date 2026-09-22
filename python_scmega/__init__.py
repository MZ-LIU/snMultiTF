"""
PyMEGA: Single-cell Multiomic Enhancer-based Gene Regulatory Network inference

A Python implementation of scMEGA for inferring gene regulatory networks 
from single-cell multiome data, particularly optimized for 10X Genomics 
multiome datasets.

This package provides a complete workflow for:
1. Data processing (10X multiome data loading and QC)
2. Cell type annotation 
3. Multiome data integration
4. Trajectory analysis
5. Feature selection (TFs and genes)
6. Peak-to-gene linking
7. Gene regulatory network inference
8. Comprehensive visualization

Author: PyMEGA Development Team
License: MIT

References:
Li, Z., Nagai, J.S., Kuppe, C. et al. scMEGA: single-cell multi-omic 
enhancer-based gene regulatory network inference. 
Bioinformatics Advances 3, vbad003 (2023).
"""
# # Fix rpy2 Windows cffi issue BEFORE any imports
# import os
# os.environ['RPY2_CFFI_MODE'] = 'ABI'


__version__ = "0.1.0"
__author__ = "PyMEGA Development Team"

# Core data structures (Stage 1)
# from .data_processing.multiome_io import MultiomeData, load_10x_multiome

# Utility functions (Stage 1)
from .utils.helpers import get_quantiles, center_rolling_mean

# Stage 2 implementations - Data processing module
from .data_processing.quality_control import (
    quality_control_separate,
    plot_qc_metrics,
    calculate_rna_qc_metrics,
    filter_rna_cells_features,
    filter_atac_data,
    sync_common_cells,
    get_qc_summary
)
from .data_processing.normalization import (
    normalize_rna,
    normalize_atac,
    normalize_multiome_muon,
    get_normalization_summary_muon
)

# Stage 4 implementations - Multiome integration module
from .multiome_integration.coembedding import coembed_data

# Stage 5 implementations - Trajectory analysis module
from .trajectory_analysis.trajectory_inference import add_trajectory, add_trajectory_archr
from .trajectory_analysis.pseudotime_analysis import get_trajectory_data


# Stage 6 implementations - Feature selection module
from .feature_selection.tf_selection import select_tfs
from .feature_selection.gene_selection import select_genes
# Note: select_trajectory_features removed - use select_trajectory_genes from trajectory_analysis
# R-based chromVAR (via rpy2)
from .feature_selection.motif_analysis import (
    load_jaspar_motifs_r, add_motifs_r, run_chromvar_r
)

# Stage 7 implementations - Peak-gene linking module (R-compatible)
from .peak_gene_linking.peak_to_gene import (
    link_peaks_to_genes,
    load_gene_annotation,
    parse_peak_annotation
)


# Stage 8 implementations - Network inference module
# Core GRN inference (R's GetGRN and GetTFGeneCorrelation equivalents)
from .network_inference.grn_inference import (
    infer_grn,
    get_tf_gene_correlation,
    save_switchtfi_inputs_from_baseline
)
# Add R-compatible network visualization from network_inference
from .network_inference.network_construction import (
    plot_grn_network ,  # R-compatible version
    create_tf_timepoint_dict,
    plot_grn_network_bipartite,
    plot_grn_network_circular
)

# Stage 9 implementations - Visualization module
from .visualization.network_plots import (
    plot_network_centrality, plot_network_topology,
    plot_tf_targets, create_network_layout
)
from .visualization.trajectory_plots import (
    plot_trajectory
)
from .visualization.heatmaps import (
    tf_gene_correlation_heatmap, gene_heatmap, tf_heatmap
)


def get_correlation(*args, **kwargs):
    """TF-gene correlation analysis. Alias for get_tf_gene_correlation."""
    return get_tf_gene_correlation(*args, **kwargs)

def deprecated_placeholder(*args, **kwargs):
    """Deprecated placeholder function."""
    raise NotImplementedError("This function has been replaced by stage implementations")

def deprecated_stage8_placeholder(*args, **kwargs):
    """Deprecated stage 8 placeholder function."""
    raise NotImplementedError("This function has been replaced by stage 8 implementations")

def deprecated_grn_placeholder(*args, **kwargs):
    """Deprecated GRN inference placeholder function."""
    raise NotImplementedError("This function has been replaced by stage 8 implementations")

__all__ = [
    # Core classes and data loading (Stage 1)
    # MuData 
    
    # Utility functions (Stage 1)
    "get_quantiles", 
    "center_rolling_mean",
    
    # Data processing (Stage 2)
    "quality_control_separate",
    "calculate_rna_qc_metrics",
    "filter_rna_cells_features",
    "filter_atac_data",
    "sync_common_cells",
    "get_qc_summary",
    "plot_qc_metrics",
    "normalize_rna",
    "normalize_atac",
    "normalize_multiome_muon",
    "get_normalization_summary_muon",
    
    # Multiome integration (Stage 4)
    "coembed_data",
    
    # Trajectory analysis (Stage 5)
    "add_trajectory",
    "add_trajectory_archr",
    "get_trajectory_data",
    
    # Feature selection (Stage 6)
    "select_tfs",
    "select_genes",
    
    # R-based chromVAR (via rpy2) - motif analysis functions
    "load_jaspar_motifs_r",
    "add_motifs_r",
    "run_chromvar_r",
 
    # Peak-gene linking - R-compatible (Stage 7)
    "link_peaks_to_genes",  
    "load_gene_annotation",
    "parse_peak_annotation",
   
    
    # Network inference (Stage 8)
    "infer_grn",                                    
    "get_tf_gene_correlation",                     
    "save_switchtfi_inputs_from_baseline",
    
    # Network visualization (R-compatible)
    "plot_grn_network",  
    "plot_grn_network_bipartite",  
    "plot_grn_network_circular",
    "create_tf_timepoint_dict",    

    # Visualization (Stage 9)
    # Network plots
    "plot_network_centrality",
    "plot_network_topology",
    "plot_tf_targets",
    "create_network_layout",
    
    # Trajectory plots
    "plot_trajectory",

    # Heatmaps
    "tf_gene_correlation_heatmap",
    "gene_heatmap",
    "tf_heatmap",
]
