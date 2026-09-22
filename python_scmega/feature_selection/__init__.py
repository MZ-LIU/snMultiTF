"""
Feature Selection Module

This module provides comprehensive tools for feature selection in single-cell
multiome data analysis, including TF selection, gene selection, trajectory-based
feature selection, and motif analysis.

Key Components:
- tf_selection: TF selection algorithms (SelectTFs equivalent)
- gene_selection: Gene selection algorithms (SelectGenes equivalent)  
- trajectory_features: Trajectory-specific feature selection
- motif_analysis: ChromVAR-equivalent motif analysis

Main Functions:
- select_tfs: Select transcription factors based on various criteria
- select_genes: Select genes based on variance, trajectory, or marker analysis
- select_trajectory_features: Comprehensive trajectory feature selection
- run_chromvar: ChromVAR-equivalent motif activity analysis

Example Usage:
```python
import python_scmega as pymega

# TF selection
selected_tfs = pymega.select_tfs(
    multiome,
    trajectory_name="Trajectory",
    method="activity",
    n_tfs=50
)

# Gene selection (R SelectGenes equivalent)
gene_results = pymega.select_genes(
    multiome,
    var_cutoff_gene=0.9,
    trajectory_name="Trajectory",
    distance_cutoff=2000,
    cor_cutoff=0.0,
    fdr_cutoff=1e-04
)
selected_genes = gene_results['selected_genes']
p2g_links = gene_results['p2g']

# Comprehensive trajectory feature selection
trajectory_features = pymega.select_trajectory_features(
    multiome,
    trajectory_name="Trajectory",
    n_tfs=30,
    n_genes=500
)

# ChromVAR motif activity analysis (calculate motif activity, use GC content matching)
# R equivalent: 
#   obj <- addMotifAnnotations(obj, motif.set = motif.set)
#   obj <- runChromVAR(obj, genome = BSgenome.Hsapiens.UCSC.hg38)
#
# Note: A real motif database must be provided, simulated data cannot be used
# from pyjaspar import jaspardb
# motif_db = jaspardb.load_jaspar_motifs('JASPAR2020_CORE')
#
# # Optional: Compute GC content for more accurate background matching (recommended)
# gc_content = pymega.compute_gc_content(multiome.atac.var, genome_fasta='hg38.fa')
# if gc_content is not None:
#     multiome.atac.var['gc_content'] = gc_content
#
# chromvar_results = pymega.run_chromvar(
#     multiome,
# motif_annotations=motif_db, # REQUIRED: real motif database
#     assay="ATAC"
# )
# # chromVAR will automatically use GC content if available in adata.var
```
"""

from .tf_selection import (
    select_tfs,
)

from .gene_selection import (
    select_genes,
    #SelectGenes  # R-compatible alias
)

# Note: trajectory_features module has been removed
# Trajectory-specific gene selection is available in trajectory_analysis module

# R-based chromVAR functions (via rpy2)
from .motif_analysis import (
    load_jaspar_motifs_r,
    add_motifs_r,
    run_chromvar_r
)

__all__ = [
    # TF selection
    'select_tfs',
    
    # Gene selection
    'select_genes',
    #'SelectGenes',  # R-compatible alias
    
    # Note: Trajectory-specific features moved to trajectory_analysis module
    # Use pymega.select_trajectory_genes from trajectory_analysis instead
    
    # R-based chromVAR (via rpy2) - motif analysis functions
    'load_jaspar_motifs_r',
    'add_motifs_r',
    'run_chromvar_r',
]
