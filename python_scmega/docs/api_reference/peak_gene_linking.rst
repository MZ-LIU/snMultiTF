Peak-Gene Linking Module
========================

The peak-gene linking module provides comprehensive tools for connecting ATAC-seq 
peaks to their target genes through distance-based, correlation-based, and 
regulatory element-based approaches.

This module handles the complete peak-gene linking workflow from basic distance 
calculations through sophisticated regulatory element annotation and functional 
scoring.

Peak-to-Gene Linking
--------------------

.. automodule:: python_scmega.peak_gene_linking.peak_to_gene
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
~~~~~~~~~~~~~

.. autofunction:: python_scmega.peak_gene_linking.peak_to_gene.link_peaks_to_genes
.. autofunction:: python_scmega.peak_gene_linking.peak_to_gene.calculate_peak_gene_distance
.. autofunction:: python_scmega.peak_gene_linking.peak_to_gene.filter_links_by_distance
.. autofunction:: python_scmega.peak_gene_linking.peak_to_gene.score_peak_gene_links
.. autofunction:: python_scmega.peak_gene_linking.peak_to_gene.validate_peak_gene_links

Example Usage
~~~~~~~~~~~~~

.. code-block:: python

    import python_scmega as pymega
    
    # Basic peak-gene linking
    peak_gene_links = pymega.link_peaks_to_genes(
        multiome,
        genome="hg38",
        max_distance=250000,
        correlation_threshold=0.1
    )
    
    # Distance-based linking
    distance_links = pymega.link_peaks_to_genes(
        multiome,
        method="distance",
        max_distance=100000
    )
    
    # Combined correlation and distance linking
    combined_links = pymega.link_peaks_to_genes(
        multiome,
        method="combined",
        max_distance=200000,
        correlation_threshold=0.15,
        p_value_threshold=0.05
    )
    
    # With custom annotations
    peak_gene_links = pymega.link_peaks_to_genes(
        peak_data=multiome.atac,
        gene_data=multiome.rna,
        peak_annotation=peak_annotations,
        gene_annotation=gene_annotations,
        genome="hg38"
    )

Correlation Analysis
--------------------

.. automodule:: python_scmega.peak_gene_linking.correlation_analysis
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
~~~~~~~~~~~~~

.. autofunction:: python_scmega.peak_gene_linking.correlation_analysis.calculate_peak_gene_correlation
.. autofunction:: python_scmega.peak_gene_linking.correlation_analysis.batch_correlation_analysis
.. autofunction:: python_scmega.peak_gene_linking.correlation_analysis.correlation_significance_test
.. autofunction:: python_scmega.peak_gene_linking.correlation_analysis.partial_correlation_analysis
.. autofunction:: python_scmega.peak_gene_linking.correlation_analysis.temporal_correlation_analysis

Example Usage
~~~~~~~~~~~~~

.. code-block:: python

    # Basic correlation analysis
    correlations = pymega.calculate_peak_gene_correlation(
        multiome,
        method="pearson",
        min_cells=20,
        batch_size=1000
    )
    
    # Spearman correlation
    spearman_corrs = pymega.calculate_peak_gene_correlation(
        multiome,
        method="spearman",
        min_expression=0.01
    )
    
    # Correlation with specific features
    correlations = pymega.calculate_peak_gene_correlation(
        multiome,
        peak_list=["peak_001", "peak_002"],
        gene_list=["gene_001", "gene_002"],
        method="pearson"
    )
    
    # Statistical significance testing
    significant_corrs = pymega.correlation_significance_test(
        correlations,
        multiple_testing_method="bonferroni",
        alpha=0.05
    )
    
    # Partial correlation controlling for confounders
    partial_corrs = pymega.partial_correlation_analysis(
        peak_data=multiome.atac,
        gene_data=multiome.rna,
        confounding_factors=["cell_cycle_phase", "batch"]
    )
    
    # Temporal correlation analysis
    temporal_corrs = pymega.temporal_correlation_analysis(
        peak_data=multiome.atac,
        gene_data=multiome.rna,
        time_points=multiome.rna.obs["pseudotime"]
    )

Genomic Distance
----------------

.. automodule:: python_scmega.peak_gene_linking.genomic_distance
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
~~~~~~~~~~~~~

.. autofunction:: python_scmega.peak_gene_linking.genomic_distance.calculate_genomic_distance
.. autofunction:: python_scmega.peak_gene_linking.genomic_distance.linear_distance
.. autofunction:: python_scmega.peak_gene_linking.genomic_distance.tad_aware_distance
.. autofunction:: python_scmega.peak_gene_linking.genomic_distance.estimate_3d_distance
.. autofunction:: python_scmega.peak_gene_linking.genomic_distance.distance_based_scoring

Example Usage
~~~~~~~~~~~~~

.. code-block:: python

    # Linear genomic distance
    linear_distances = pymega.calculate_genomic_distance(
        peak_annotations,
        gene_annotations,
        distance_type="linear",
        max_distance=500000
    )
    
    # TAD-aware distance
    tad_distances = pymega.calculate_genomic_distance(
        peak_annotations,
        gene_annotations,
        distance_type="tad_aware",
        tad_boundaries=tad_data,
        tad_penalty=2.0
    )
    
    # 3D distance estimation
    distance_3d = pymega.calculate_genomic_distance(
        peak_annotations,
        gene_annotations,
        distance_type="3d",
        scaling_factor=0.75
    )
    
    # Contact-based distance using Hi-C
    contact_distances = pymega.calculate_genomic_distance(
        peak_annotations,
        gene_annotations,
        distance_type="contact",
        contact_matrix=hic_matrix,
        resolution=10000
    )
    
    # Distance binning and analysis
    binned_distances = pymega.calculate_distance_bins(
        distances,
        bin_method="quantile",
        n_bins=10
    )
    
    # Distance distribution analysis
    dist_analysis = pymega.analyze_distance_distribution(
        distances,
        distance_column="distance"
    )

Regulatory Elements
-------------------

.. automodule:: python_scmega.peak_gene_linking.regulatory_elements
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
~~~~~~~~~~~~~

.. note::
   **Regulatory element annotation features removed**
   
   The following functions were removed as they don't exist in the R version of scMEGA:
   
   - ``annotate_regulatory_elements()``
   - ``identify_promoters()``
   - ``identify_enhancers()``
   - ``functional_element_scoring()``
   
   Only distance-based analysis is supported to match the R version's PeakToGene function.

Complete Peak-Gene Linking Workflow
------------------------------------

.. code-block:: python

    import python_scmega as pymega
    
    # 1. Load and preprocess data
    multiome = pymega.load_10x_multiome("data.h5")
    multiome = pymega.quality_control_multiome(multiome)
    multiome = pymega.normalize_multiome(multiome)
    
    # 2. Create genomic annotations (or load existing ones)
    peak_annotations = create_peak_annotations()  # Your annotation function
    gene_annotations = create_gene_annotations()  # Your annotation function
    
    # 3. Basic peak-gene linking
    peak_gene_links = pymega.link_peaks_to_genes(
        multiome,
        genome="hg38",
        max_distance=250000,
        method="combined",
        peak_annotation=peak_annotations,
        gene_annotation=gene_annotations
    )
    
    # 4. Detailed correlation analysis
    correlations = pymega.calculate_peak_gene_correlation(
        multiome,
        method="pearson",
        min_cells=20,
        batch_size=1000
    )
    
    # Statistical significance testing
    significant_corrs = pymega.correlation_significance_test(
        correlations,
        multiple_testing_method="bonferroni"
    )
    
    # 5. Distance analysis (matches R version)
    distances = pymega.calculate_genomic_distance(
        peak_annotations,
        gene_annotations,
        max_distance=250000  # Default: 250kb, same as R version
    )
    
    # 6. Integration and validation
    validation_results = pymega.validate_peak_gene_links(
        peak_gene_links,
        validation_method="distance_distribution"
    )

R Equivalence
-------------

This module provides Python equivalents for key R functions:

.. list-table:: R to Python Function Mapping
   :widths: 40 40 20
   :header-rows: 1

   * - R Function
     - PyMEGA Equivalent
     - Module
   * - ``PeakToGene()``
     - ``link_peaks_to_genes()``
     - peak_to_gene
   * - ``CalculatePeakGeneCorrelation()``
     - ``calculate_peak_gene_correlation()``
     - correlation_analysis
   * - ``CalculateGenomicDistance()``
     - ``calculate_genomic_distance()``
     - genomic_distance

Algorithm Details
-----------------

Peak-Gene Linking Methods
~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Distance-based**: Links peaks to genes within a maximum genomic distance
2. **Correlation-based**: Links peaks to genes based on expression correlation
3. **Combined**: Integrates distance and correlation with weighted scoring

Correlation Methods
~~~~~~~~~~~~~~~~~~~

1. **Pearson**: Linear correlation coefficient
2. **Spearman**: Rank-based correlation coefficient  
3. **Kendall**: Tau rank correlation coefficient
4. **Partial**: Correlation controlling for confounding factors

Distance Types
~~~~~~~~~~~~~~

1. **Linear**: Simple genomic distance (bp)
2. **TAD-aware**: Distance with TAD boundary penalties
3. **3D**: Estimated 3D nuclear distance using polymer scaling
4. **Contact-based**: Distance derived from Hi-C contact frequencies

Regulatory Element Classification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Promoters**: Peaks within promoter windows (≤2kb from TSS)
2. **Enhancers**: Distal peaks with regulatory potential
3. **Silencers**: Repressive regulatory elements
4. **Insulators**: Boundary elements

Performance Considerations
--------------------------

- **Memory Efficiency**: Batch processing for large datasets
- **Scalability**: Optimized for datasets with 100K+ peaks and genes
- **Speed**: Vectorized operations and optional parallel processing
- **Flexibility**: Multiple linking methods and customizable parameters

Dependencies
------------

Core dependencies:
- ``numpy``
- ``pandas``
- ``scipy``
- ``scikit-learn``
- ``anndata``

Optional dependencies:
- ``statsmodels`` (for multiple testing correction)
- ``concurrent.futures`` (for parallel processing)

Method Comparison
-----------------

Linking Methods
~~~~~~~~~~~~~~~

.. list-table::
   :widths: 20 30 25 25
   :header-rows: 1

   * - Method
     - Best For
     - Advantages
     - Limitations
   * - Distance
     - Promoter identification
     - Fast and simple
     - Ignores expression
   * - Correlation
     - Functional relationships
     - Expression-based
     - Requires shared cells
   * - Combined
     - Comprehensive analysis
     - Integrates multiple signals
     - More complex parameters

Distance Types
~~~~~~~~~~~~~~

.. list-table::
   :widths: 20 30 25 25
   :header-rows: 1

   * - Type
     - Best For
     - Advantages
     - Limitations
   * - Linear
     - General analysis
     - Simple and interpretable
     - Ignores 3D structure
   * - TAD-aware
     - Chromatin domains
     - Biologically relevant
     - Requires TAD data
   * - 3D
     - Nuclear organization
     - Realistic distances
     - Approximate estimation
   * - Contact-based
     - Hi-C integration
     - Experimental evidence
     - Requires Hi-C data

Use Cases
---------

Enhancer-Promoter Linking
~~~~~~~~~~~~~~~~~~~~~~~~~~
- Identify distal regulatory elements
- Map enhancer-target gene relationships
- Analyze tissue-specific regulation

Promoter Analysis
~~~~~~~~~~~~~~~~~
- Annotate promoter-associated peaks
- Identify core vs. extended promoters
- Analyze promoter accessibility

Disease Studies
~~~~~~~~~~~~~~~
- Link disease-associated variants to target genes
- Identify regulatory disruptions
- Map pathogenic mechanisms

Development and Differentiation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Track regulatory changes during development
- Identify lineage-specific regulatory elements
- Map trajectory-associated peak-gene links

Best Practices
--------------

Data Preprocessing
~~~~~~~~~~~~~~~~~~
1. Ensure proper quality control before linking
2. Normalize both ATAC and RNA data appropriately
3. Filter low-quality peaks and genes

Parameter Selection
~~~~~~~~~~~~~~~~~~~
1. Start with default parameters for exploration
2. Adjust distance thresholds based on biology
3. Use multiple methods for robust results
4. Validate key findings experimentally

Annotation Quality
~~~~~~~~~~~~~~~~~~
1. Use high-quality genomic annotations
2. Include multiple annotation sources
3. Validate regulatory element classifications
4. Consider cell type-specific annotations

Statistical Considerations
~~~~~~~~~~~~~~~~~~~~~~~~~~
1. Apply multiple testing correction
2. Set appropriate significance thresholds
3. Consider effect sizes, not just p-values
4. Validate with independent datasets

Troubleshooting
---------------

Common Issues
~~~~~~~~~~~~~

**No links identified**
- Check distance thresholds (may be too restrictive)
- Verify data preprocessing quality
- Ensure shared cells between ATAC and RNA data

**Low correlation values**
- Check data normalization
- Consider batch effects
- Verify cell type homogeneity

**Memory errors**
- Reduce batch size for correlation analysis
- Use feature subsets for initial exploration
- Consider parallel processing options

**Slow performance**
- Use smaller batch sizes
- Filter features before analysis
- Enable parallel processing

Performance Optimization
~~~~~~~~~~~~~~~~~~~~~~~~

**For large datasets**
- Use batch processing for correlations
- Subset to high-quality features first
- Consider distance pre-filtering

**For many peak-gene pairs**
- Use parallel processing
- Implement progressive filtering
- Cache intermediate results

See Also
--------

- :doc:`feature_selection`: Feature selection for regulatory analysis
- :doc:`trajectory_analysis`: Trajectory-based regulatory dynamics
- :doc:`../user_guide/installation`: Installation instructions
