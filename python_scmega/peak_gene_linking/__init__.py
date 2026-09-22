"""
Peak-Gene Linking Module

This module provides tools for linking ATAC-seq peaks to their target genes
using trajectory-based correlation analysis, matching R scMEGA's PeakToGene function.

Distance Calculation (R-compatible):
    The distance calculation method matches R's PeakToGene function:
    distance = abs(peak_center - gene_TSS)
    
    Where:
    - peak_center = (peak_start + peak_end) / 2
    - gene_TSS = gene_start (for + strand) or gene_end (for - strand)
    - Only same-chromosome pairs are considered
    - Default max_distance = 250,000 bp (250kb)

Key Components:
- peak_to_gene: Peak-gene linking algorithms (R PeakToGene equivalent)

Main Functions:
- link_peaks_to_genes: R-compatible trajectory-based peak-gene linking
- fast_correlation: High-performance correlation calculation (R's rowCorCpp equivalent)
- load_gene_annotation: Load real gene annotations from .rda files
- parse_peak_annotation: Parse peak coordinates from peak names

Example Usage:
```python
import python_scmega as pymega

# Get trajectory matrices
traj_rna = pymega.get_trajectory(multiome, assay="RNA", ...)
traj_atac = pymega.get_trajectory(multiome, assay="ATAC", ...)

# Link peaks to genes (R-compatible)
df_p2g = pymega.link_peaks_to_genes(
    peak_mat=traj_atac['smooth_matrix'],
    gene_mat=traj_rna['smooth_matrix'],
    genome="hg38",
    max_distance=250000  # 250kb default, same as R
)
```
"""

from .peak_to_gene import (
    link_peaks_to_genes,  # Equivalent to PeakToGene
)



__all__ = [
    # Peak-to-gene linking (R-compatible)
    'link_peaks_to_genes',

]



