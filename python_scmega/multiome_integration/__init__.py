"""
Multiome Integration Module

This module provides comprehensive tools for integrating multi-modal single-cell data,
specifically designed for scenarios involving RNA and ATAC data from different sources
or experimental conditions.

Key Components:
- coembedding: Multi-modal co-embedding and paired object creation
- cell_pairing: Cell pairing algorithms for non-paired datasets
- batch_correction: Batch effect correction across modalities
- integration_qc: Quality control and evaluation of integration results

Main Functions:
- coembed_data: Create joint embedding of RNA and ATAC data
- pair_cells: Pair cells between RNA and ATAC datasets
- correct_batch_effects: Correct batch effects in multiome data
- evaluate_integration_quality: Assess integration quality

Example Usage:
```python
import python_scmega as pymega

# Co-embedding for paired data
coembed_result = pymega.coembed_data(rna_adata, atac_adata)

# Cell pairing for unpaired data
paired_multiome = pymega.pair_cells(
    rna_adata, atac_adata,
    method="optimal_transport"
)

# Batch correction
corrected_data = pymega.correct_batch_effects(
    multiome_data,
    batch_key="batch",
    method="harmony"
)

# Integration quality assessment
quality_metrics = pymega.evaluate_integration_quality(
    integrated_data,
    batch_key="batch",
    cell_type_key="celltype"
)
```
"""

from .coembedding import (
    coembed_data
)

__all__ = [
    # Co-embedding
    'coembed_data',
    
]