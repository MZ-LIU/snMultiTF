Feature Selection Module
========================

The feature selection module provides comprehensive tools for selecting informative 
features in single-cell multiome data analysis, including TF selection, gene selection,
trajectory-based feature selection, and motif analysis.

This module handles the complete feature selection workflow from basic variance-based
selection through sophisticated trajectory-specific and motif-based approaches.

TF Selection
------------

.. automodule:: python_scmega.feature_selection.tf_selection
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
~~~~~~~~~~~~~

.. autofunction:: python_scmega.feature_selection.tf_selection.select_tfs
.. autofunction:: python_scmega.feature_selection.tf_selection.rank_tfs_by_activity
.. autofunction:: python_scmega.feature_selection.tf_selection.select_trajectory_tfs
.. autofunction:: python_scmega.feature_selection.tf_selection.filter_tfs_by_expression
.. autofunction:: python_scmega.feature_selection.tf_selection.validate_tf_selection

Example Usage
~~~~~~~~~~~~~

.. code-block:: python

    import python_scmega as pymega
    
    # Activity-based TF selection
    selected_tfs = pymega.select_tfs(
        multiome,
        method="activity",
        n_tfs=50,
        trajectory_name="Trajectory"
    )
    
    # Trajectory-specific TF selection
    trajectory_tfs = pymega.select_tfs(
        multiome,
        method="trajectory",
        trajectory_name="Trajectory",
        n_tfs=30
    )
    
    # Variance-based TF selection
    variable_tfs = pymega.select_tfs(
        multiome,
        method="variance",
        n_tfs=40
    )
    
  
Gene Selection
--------------

.. automodule:: python_scmega.feature_selection.gene_selection
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
~~~~~~~~~~~~~

.. autofunction:: python_scmega.feature_selection.gene_selection.select_genes

.. note::
   The ``select_genes()`` function provides R-compatible peak-to-gene based gene selection,
   equivalent to the R scMEGA ``SelectGenes()`` function.
   
   For trajectory-specific gene selection, see :func:`python_scmega.trajectory_analysis.trajectory_inference.select_trajectory_genes`.


Example Usage
~~~~~~~~~~~~~

.. code-block:: python

    import python_scmega as pymega
    
    # R SelectGenes equivalent: Peak-to-gene based gene selection
    gene_results = pymega.select_genes(
        multiome,
        atac_assay="ATAC",
        rna_assay="RNA",
        var_cutoff_gene=0.9,          # Top 10% variable genes
        trajectory_name="Trajectory",
        distance_cutoff=2000,          # Minimum distance (bp)
        cor_cutoff=0.0,                # Correlation threshold
        fdr_cutoff=1e-04,              # FDR threshold
        genome="hg38"
    )
    
    # Access results
    selected_genes = gene_results['selected_genes']
    p2g_links = gene_results['p2g']
    print(f"Selected {len(selected_genes)} genes with {len(p2g_links)} peak-gene links")
    
    # For trajectory-specific gene selection, use trajectory_analysis module
    from python_scmega.trajectory_analysis import select_trajectory_genes
    
    trajectory_genes = select_trajectory_genes(
        multiome,
        trajectory_name="Trajectory",
        n_genes=500,
        method="variance"  # Options: "variance", "correlation", "differential"
    )

Trajectory-Specific Gene Selection
-----------------------------------

.. deprecated:: 0.1.0
   The ``trajectory_features`` module has been removed. Trajectory-specific gene 
   selection is now available in the trajectory analysis module.

.. note::
   **Migration Guide:**
   
   - Old: ``pymega.select_trajectory_features()``
   - New: Use ``pymega.select_trajectory_genes()`` from trajectory_analysis module
   
   For comprehensive feature selection, combine TF and gene selection:

.. code-block:: python

    import python_scmega as pymega
    
    # Step 1: Select TFs
    selected_tfs = pymega.select_tfs(
        multiome,
        trajectory_name="Trajectory",
        method="activity",
        n_tfs=50
    )
    
    # Step 2: Select trajectory-specific genes
    trajectory_genes = pymega.select_trajectory_genes(
        multiome,
        trajectory_name="Trajectory",
        assay="RNA",
        n_genes=500,
        method="variance"  # or "correlation", "differential"
    )
    
    # Step 3: Or use peak-to-gene based gene selection
    gene_results = pymega.select_genes(
        multiome,
        var_cutoff_gene=0.9,
        trajectory_name="Trajectory",
        distance_cutoff=2000,
        cor_cutoff=0.0,
        fdr_cutoff=1e-04
    )
    p2g_selected_genes = gene_results['selected_genes']
    
    # Step 4: Combine features
    final_features = list(set(selected_tfs + trajectory_genes))
    
    print(f"Feature selection results:")
    print(f"  TFs: {len(selected_tfs)}")
    print(f"  Trajectory genes: {len(trajectory_genes)}")
    print(f"  Peak-to-gene genes: {len(p2g_selected_genes)}")
    print(f"  Combined features: {len(final_features)}")

See :doc:`trajectory_analysis` for detailed documentation on ``select_trajectory_genes()``.


Motif Analysis
--------------

.. automodule:: python_scmega.feature_selection.motif_analysis
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
~~~~~~~~~~~~~

.. autofunction:: python_scmega.feature_selection.motif_analysis.run_chromvar
.. autofunction:: python_scmega.feature_selection.motif_analysis.compute_gc_content
.. autofunction:: python_scmega.feature_selection.motif_analysis.calculate_motif_activity
.. autofunction:: python_scmega.feature_selection.motif_analysis.identify_variable_motifs
.. autofunction:: python_scmega.feature_selection.motif_analysis.motif_enrichment_analysis
.. autofunction:: python_scmega.feature_selection.motif_analysis.select_motif_features

Example Usage
~~~~~~~~~~~~~

.. code-block:: python

    # ChromVAR analysis with GC-matched background (like R version)
    
    # Step 1: Compute GC content (recommended)
    gc_content = pymega.compute_gc_content(
        multiome.atac.var,
        genome_fasta='hg38.fa'
    )
    if gc_content is not None:
        multiome.atac.var['gc_content'] = gc_content
    
    # Step 2: Run chromVAR (will use GC content for background matching)
    chromvar_results = pymega.run_chromvar(
        multiome,
        assay="ATAC",
        motif_annotations=motif_database,
        n_components=50
    )
    
    # Select variable motifs
    motif_features = pymega.select_motif_features(
        chromvar_results,
        selection_method="variable",
        n_features=100
    )
    
    # Motif enrichment analysis
    enrichment_results = pymega.motif_enrichment_analysis(
        accessibility_data,
        motif_annotations,
        cell_groups
    )
    
    # Create motif activity matrix
    motif_adata = pymega.create_motif_activity_matrix(
        multiome,
        chromvar_results,
        assay="ATAC"
    )

Complete Feature Selection Workflow
------------------------------------

.. code-block:: python

    import python_scmega as pymega
    
    # 1. Load and preprocess data
    multiome = pymega.load_10x_multiome("data.h5")
    multiome = pymega.quality_control_multiome(multiome)
    multiome = pymega.normalize_multiome(multiome)
    
    # 2. Trajectory inference (R-compatible ArchR algorithm)
    multiome = pymega.add_trajectory(
        multiome,
        trajectory=["0", "1", "2"],     # R: trajectory = c("0", "1", "2")
        group_by="leiden",              # R: group.by = "leiden"
        reduction="pca",
        dims=list(range(30)),
        pre_filter_quantile=0.9,
        post_filter_quantile=0.9,
        dof=250,
        spar=1.0,
        name="Trajectory"
    )
    
    # 3. TF selection
    selected_tfs = pymega.select_tfs(
        multiome,
        method="activity",
        trajectory_name="Trajectory",
        n_tfs=50
    )
    
    # 4. Gene selection (R SelectGenes equivalent)
    gene_results = pymega.select_genes(
        multiome,
        var_cutoff_gene=0.9,
        trajectory_name="Trajectory",
        distance_cutoff=2000,
        cor_cutoff=0.0,
        fdr_cutoff=1e-04,
        genome="hg38"
    )
    selected_genes = gene_results['selected_genes']
    p2g_links = gene_results['p2g']
    
    # 5. Trajectory-specific gene selection (alternative approach)
    trajectory_genes = pymega.select_trajectory_genes(
        multiome,
        trajectory_name="Trajectory",
        n_genes=500,
        method="variance"
    )
    
    # 6. ChromVAR motif analysis (with GC-matched background)
    # Step 6.1: Compute GC content (recommended for accurate background matching)
    gc_content = pymega.compute_gc_content(
        multiome.atac.var,
        genome_fasta='hg38.fa'  # or mm10.fa for mouse
    )
    if gc_content is not None:
        multiome.atac.var['gc_content'] = gc_content
    
    # Step 6.2: Run chromVAR (will use GC content for background matching)
    chromvar_results = pymega.run_chromvar(
        multiome,
        assay="ATAC",
        motif_annotations=motif_db
    )
    
    print("Feature selection completed!")
    print(f"  TFs: {len(selected_tfs)}")
    print(f"  Peak-to-gene genes: {len(selected_genes)}")
    print(f"  Trajectory genes: {len(trajectory_genes)}")

R Equivalence
-------------

This module provides Python equivalents for key R functions:

.. list-table:: R to Python Function Mapping
   :widths: 40 40 20
   :header-rows: 1

   * - R Function
     - PyMEGA Equivalent
     - Module
   * - ``SelectTFs()``
     - ``select_tfs()``
     - tf_selection
   * - ``SelectGenes()``
     - ``select_genes()``
     - gene_selection
   * - ``SelectGenes()`` (trajectory-specific)
     - ``select_trajectory_genes()``
     - trajectory_analysis.trajectory_inference
   * - ``runChromVAR()``
     - ``run_chromvar()``
     - motif_analysis
   * - ``addMotifAnnotations()``
     - ``calculate_motif_activity()``
     - motif_analysis

Algorithm Details
-----------------

TF Selection Methods
~~~~~~~~~~~~~~~~~~~~

1. **Activity-based**: Combines expression level, variability, and trajectory correlation
2. **Trajectory-based**: Selects TFs with highest correlation to pseudotime
3. **Variance-based**: Selects TFs with highest coefficient of variation
4. **Expression-based**: Selects TFs with highest mean expression

Gene Selection Methods
~~~~~~~~~~~~~~~~~~~~~~

1. **Peak-to-gene based** (R SelectGenes): Links peaks to genes based on trajectory and correlation
2. **Trajectory-based**: Selects genes correlated with pseudotime (use select_trajectory_genes)
3. **Variance-based**: Selects highly variable genes across trajectory
4. **Differential**: Uses differential expression analysis

ChromVAR Implementation
~~~~~~~~~~~~~~~~~~~~~~~

1. **Normalization**: TF-IDF or log normalization of accessibility data
2. **Motif Activity**: Background-corrected mean accessibility across motif sites
3. **Variability**: Median absolute deviation, standard deviation, or CV
4. **Enrichment**: Mann-Whitney U tests with multiple testing correction

Performance Considerations
--------------------------

- **Memory Efficiency**: Batch processing for large gene sets
- **Scalability**: Optimized for datasets with 100K+ cells
- **Speed**: Vectorized operations using numpy and scipy
- **Flexibility**: Multiple methods and customizable parameters

Dependencies
------------

Core dependencies:
- ``numpy``
- ``pandas``
- ``scipy``
- ``scikit-learn``
- ``scanpy``
- ``anndata``

Optional dependencies:
- ``statsmodels`` (for multiple testing correction)
- ``matplotlib`` (for plotting functions)

Method Comparison
-----------------

TF Selection Methods
~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 20 30 25 25
   :header-rows: 1

   * - Method
     - Best For
     - Advantages
     - Limitations
   * - Activity
     - General TF discovery
     - Combines multiple criteria
     - May miss context-specific TFs
   * - Trajectory
     - Developmental analysis
     - Trajectory-specific
     - Requires trajectory data
   * - Variance
     - Variable TF identification
     - Simple and robust
     - May select noisy genes
   * - Expression
     - Highly expressed TFs
     - Identifies active TFs
     - Biased toward abundant TFs

Gene Selection Methods
~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 20 30 25 25
   :header-rows: 1

   * - Method
     - Best For
     - Advantages
     - Limitations
   * - Variance (HVG)
     - General analysis
     - Well-established
     - May miss trajectory genes
   * - Trajectory
     - Temporal analysis
     - Trajectory-specific
     - Requires trajectory
   * - Marker
     - Cell type analysis
     - Biologically relevant
     - Requires clustering
   * - Dispersion
     - Overdispersed genes
     - Accounts for mean-variance
     - Sensitive to outliers

Use Cases
---------

Developmental Biology
~~~~~~~~~~~~~~~~~~~~
- Select TFs and genes involved in cell differentiation
- Identify trajectory-specific regulatory programs
- Analyze motif activity changes during development

Cancer Research
~~~~~~~~~~~~~~~
- Identify oncogenic TFs and target genes
- Select features for tumor progression analysis
- Analyze chromatin accessibility changes

Immunology
~~~~~~~~~~
- Select immune response genes and TFs
- Identify activation-specific features
- Analyze T cell differentiation trajectories

Stem Cell Biology
~~~~~~~~~~~~~~~~~
- Select pluripotency and differentiation factors
- Identify lineage-specific features
- Analyze reprogramming trajectories

Best Practices
--------------

Data Preprocessing
~~~~~~~~~~~~~~~~~~
1. Ensure proper quality control before feature selection
2. Normalize data appropriately for each assay
3. Remove low-quality cells and features

Feature Selection Strategy
~~~~~~~~~~~~~~~~~~~~~~~~~~
1. Start with variance-based methods for exploration
2. Use trajectory-based methods for temporal analysis
3. Combine multiple methods for robust selection
4. Validate selections with known biology

Feature Combination Strategy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1. Combine TF and gene selections for comprehensive analysis
2. Use select_tfs() for transcription factor selection
3. Use select_genes() for peak-to-gene based gene selection
4. Use select_trajectory_genes() for trajectory-specific genes

Validation
~~~~~~~~~~
1. Always validate selected features
2. Check correlation with known markers
3. Assess biological relevance
4. Test robustness across parameters

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**No features selected**
- Check expression filters (min_expression, min_cells)
- Verify trajectory information is present
- Adjust correlation thresholds

**Poor trajectory correlation**
- Check trajectory inference quality
- Consider different trajectory methods
- Verify pseudotime ordering

**ChromVAR fails**
- Ensure ATAC data is properly normalized
- Check motif annotation format
- Verify peak annotations

**Memory errors**
- Reduce batch size for large datasets
- Use sparse matrices where possible
- Process features in chunks

Performance Optimization
~~~~~~~~~~~~~~~~~~~~~~~~

**For large datasets**
- Use batch processing
- Subset to highly variable features first
- Consider dimensionality reduction

**For many features**
- Pre-filter by expression
- Use parallel processing where available
- Cache intermediate results

See Also
--------

- :doc:`trajectory_analysis`: Trajectory inference and analysis
- :doc:`multiome_integration`: Multi-modal data integration
- :doc:`../user_guide/installation`: Installation instructions
