"""
Data Processing Module

Handles 10X multiome data loading, quality control, normalization, and format conversion.
This module is specifically optimized for 10X Genomics multiome datasets where RNA and ATAC
measurements come from the same cells.

Key components:
- multiome_io.py: 10X multiome data loading (Read10X_h5 equivalent)
- quality_control.py: Quality control and filtering (PercentageFeatureSet equivalent)
- normalization.py: Data normalization (NormalizeData, ScaleData, RunPCA equivalent)
- format_conversion.py: Enhanced data format conversion (Seurat <-> AnnData)

Stage 2 Implementation Status:
 multiome_io.py - 10X multiome data loading (completed in Stage 1)
 quality_control.py - Complete QC pipeline with PercentageFeatureSet equivalent
 normalization.py - RNA/ATAC normalization with NormalizeData/ScaleData/RunPCA equivalent
 format_conversion.py - Enhanced format conversion utilities

References:
- Original R code: vignettes/pbmc_10x_multiome.Rmd
- R implementation: inputdata.10x <- Read10X_h5("data.h5")
"""

# Quality control functions (Stage 2)
from .quality_control import (
    calculate_rna_qc_metrics,
    filter_rna_cells_features,  # function name in actual file
    filter_atac_data,
    sync_common_cells,
    quality_control_separate,

    plot_qc_metrics,  # function name in actual file
    get_qc_summary
)

# Normalization functions
from .normalization import (
    normalize_rna,
    normalize_atac,
    normalize_multiome_muon,
    get_normalization_summary_muon
)


__all__ = [
    
    "calculate_rna_qc_metrics",
    "filter_rna_cells_features",
    "filter_atac_data",
    "sync_common_cells",
    "quality_control_separate",
    
    "plot_qc_metrics",
    "get_qc_summary",
    
    # Normalization
    "normalize_rna",
    "normalize_atac",
    "normalize_multiome_muon",
    "get_normalization_summary_muon",
    

]
