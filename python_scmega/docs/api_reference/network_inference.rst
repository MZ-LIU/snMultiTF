Network Inference Module
=========================

The network inference module provides comprehensive tools for gene regulatory 
network (GRN) inference from single-cell multiome data, including TF-gene 
correlation analysis, network construction, validation, and target prediction.

This module implements the core network inference algorithms from the original 
scMEGA R package, providing Python equivalents for GetGRN, GetTFGeneCorrelation, 
AddTargetAssay, and other network inference functions.

GRN Inference
-------------

.. automodule:: python_scmega.network_inference.grn_inference
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
~~~~~~~~~~~~~

.. autofunction:: python_scmega.network_inference.grn_inference.infer_grn
.. autofunction:: python_scmega.network_inference.grn_inference.get_tf_gene_correlation
.. autofunction:: python_scmega.network_inference.grn_inference.calculate_tf_activity
.. autofunction:: python_scmega.network_inference.grn_inference.integrate_regulatory_evidence
.. autofunction:: python_scmega.network_inference.grn_inference.rank_regulatory_interactions
.. autofunction:: python_scmega.network_inference.grn_inference.filter_grn_by_confidence

Example Usage
~~~~~~~~~~~~~

Basic GRN inference workflow:

.. code-block:: python

    import python_scmega as pymega
    
    # Step 1: Calculate TF-gene correlations
    tf_gene_cor = pymega.get_tf_gene_correlation(
        data=multiome_data,
        tf_list=tf_list,
        gene_list=gene_list,
        method="pearson",
        min_cells=50
    )
    
    # Step 2: Infer GRN
    grn_network = pymega.infer_grn(
        motif_matching=motif_matrix,
        tf_gene_cor=tf_gene_cor,
        peak_gene_links=peak_gene_links,
        score_threshold=0.1,
        fdr_threshold=0.25
    )
    
    # Step 3: Filter high-confidence edges
    high_conf_network = pymega.filter_grn_by_confidence(
        grn_network,
        min_score=0.2,
        max_p_value=0.01,
        min_evidence=2
    )

Correlation Methods
-------------------

.. automodule:: python_scmega.network_inference.correlation_methods
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
~~~~~~~~~~~~~

.. autofunction:: python_scmega.network_inference.correlation_methods.fast_correlation
.. autofunction:: python_scmega.network_inference.correlation_methods.batch_correlation_analysis
.. autofunction:: python_scmega.network_inference.correlation_methods.sparse_correlation
.. autofunction:: python_scmega.network_inference.correlation_methods.partial_correlation_matrix
.. autofunction:: python_scmega.network_inference.correlation_methods.correlation_significance_test

Example Usage
~~~~~~~~~~~~~

High-performance correlation calculation:

.. code-block:: python

    import python_scmega as pymega
    import numpy as np
    
    # Fast correlation between two matrices
    matrix1 = np.random.normal(0, 1, (1000, 500))  # 1000 features x 500 samples
    matrix2 = np.random.normal(0, 1, (1200, 500))  # 1200 features x 500 samples
    
    # Calculate correlations using optimized algorithm
    correlations = pymega.fast_correlation(matrix1, matrix2)
    
    # Batch processing for very large matrices
    batch_correlations = pymega.batch_correlation_analysis(
        matrix1,
        matrix2,
        batch_size=200,
        n_jobs=4
    )
    
    # Sparse matrix correlation
    sparse_correlations = pymega.sparse_correlation(
        sparse_matrix1,
        sparse_matrix2,
        min_overlap=10
    )

Network Construction
--------------------

.. automodule:: python_scmega.network_inference.network_construction
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
~~~~~~~~~~~~~

.. autofunction:: python_scmega.network_inference.network_construction.construct_grn_network
.. autofunction:: python_scmega.network_inference.network_construction.build_network_topology
.. autofunction:: python_scmega.network_inference.network_construction.integrate_multimodal_evidence
.. autofunction:: python_scmega.network_inference.network_construction.score_regulatory_edges
.. autofunction:: python_scmega.network_inference.network_construction.prune_network_by_confidence
.. autofunction:: python_scmega.network_inference.network_construction.create_network_summary

Example Usage
~~~~~~~~~~~~~

Advanced network construction with multiple evidence integration:

.. code-block:: python

    import python_scmega as pymega
    
    # Construct network with weighted ensemble integration
    grn_network = pymega.construct_grn_network(
        tf_gene_correlations=tf_cor_results,
        motif_matching=motif_matrix,
        peak_gene_links=peak_links,
        integration_method="weighted_ensemble",
        confidence_threshold=0.15,
        network_type="directed",
        max_targets_per_tf=100
    )
    
    # Score regulatory edges
    scored_network = pymega.score_regulatory_edges(
        grn_network,
        evidence_sources=evidence_dict,
        scoring_method="composite"
    )
    
    # Prune network by confidence
    pruned_network = pymega.prune_network_by_confidence(
        scored_network,
        confidence_threshold=0.2,
        statistical_threshold=0.01
    )
    
    # Generate network summary
    network_summary = pymega.create_network_summary(
        pruned_network,
        include_metrics=True
    )

Network Validation
------------------

.. automodule:: python_scmega.network_inference.network_validation
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
~~~~~~~~~~~~~

.. autofunction:: python_scmega.network_inference.network_validation.validate_grn_network
.. autofunction:: python_scmega.network_inference.network_validation.calculate_network_metrics
.. autofunction:: python_scmega.network_inference.network_validation.benchmark_against_reference
.. autofunction:: python_scmega.network_inference.network_validation.assess_biological_coherence
.. autofunction:: python_scmega.network_inference.network_validation.cross_validate_network
.. autofunction:: python_scmega.network_inference.network_validation.generate_validation_report

Example Usage
~~~~~~~~~~~~~

Comprehensive network validation:

.. code-block:: python

    import python_scmega as pymega
    
    # Comprehensive network validation
    validation_results = pymega.validate_grn_network(
        grn_network=inferred_network,
        reference_networks=["string", "reactome"],
        expression_data=expression_matrix,
        validation_metrics=["topology", "biological", "statistical"],
        cross_validation=True,
        n_folds=5
    )
    
    # Calculate detailed network metrics
    network_metrics = pymega.calculate_network_metrics(
        grn_network,
        metric_types=["basic", "centrality", "clustering", "connectivity"]
    )
    
    # Benchmark against reference networks
    benchmark_results = pymega.benchmark_against_reference(
        grn_network,
        reference_networks=[reference_network_df],
        comparison_metrics=["overlap", "precision", "recall", "f1"]
    )
    
    # Assess biological coherence
    coherence_results = pymega.assess_biological_coherence(
        grn_network,
        coherence_metrics=["go_enrichment", "pathway_enrichment"]
    )
    
    # Generate validation report
    report = pymega.generate_validation_report(
        validation_results,
        report_format="text"
    )

Target Prediction
-----------------

.. automodule:: python_scmega.network_inference.target_prediction
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
~~~~~~~~~~~~~

.. autofunction:: python_scmega.network_inference.target_prediction.predict_targets
.. autofunction:: python_scmega.network_inference.target_prediction.add_target_assay
.. autofunction:: python_scmega.network_inference.target_prediction.score_tf_target_pairs
.. autofunction:: python_scmega.network_inference.target_prediction.validate_target_predictions
.. autofunction:: python_scmega.network_inference.target_prediction.prioritize_regulatory_targets
.. autofunction:: python_scmega.network_inference.target_prediction.create_target_summary

Example Usage
~~~~~~~~~~~~~

Target prediction and analysis workflow:

.. code-block:: python

    import python_scmega as pymega
    
    # Predict TF targets using multiple methods
    target_predictions = pymega.predict_targets(
        multiome_data=multiome,
        tf_list=tf_list,
        prediction_method="integrated",
        motif_data=motif_matrix,
        peak_gene_links=peak_links,
        confidence_threshold=0.15,
        max_targets_per_tf=50
    )
    
    # Add target assay to multiome data
    multiome_with_targets = pymega.add_target_assay(
        multiome_data=multiome,
        tf_targets=target_predictions,
        assay_name="predicted_targets"
    )
    
    # Score TF-target pairs
    scored_pairs = pymega.score_tf_target_pairs(
        expression_data=expression_matrix,
        tf_target_pairs=target_predictions,
        scoring_method="correlation"
    )
    
    # Validate predictions
    validation_results = pymega.validate_target_predictions(
        target_predictions,
        known_targets=known_tf_targets,
        validation_method="cross_validation"
    )
    
    # Prioritize targets
    prioritized_targets = pymega.prioritize_regulatory_targets(
        target_predictions,
        prioritization_criteria=["confidence", "expression_level"],
        expression_data=expression_matrix
    )
    
    # Create summary
    target_summary = pymega.create_target_summary(
        prioritized_targets,
        include_statistics=True
    )

Complete Workflow Example
-------------------------

Here's a complete example demonstrating the full network inference workflow:

.. code-block:: python

    import python_scmega as pymega
    import numpy as np
    import pandas as pd
    
    # Load multiome data
    multiome_data = pymega.load_10x_multiome("multiome_data.h5")
    
    # Define TFs and genes of interest
    tf_list = ["TF1", "TF2", "TF3", "TF4", "TF5"]
    gene_list = multiome_data.rna.var.index.tolist()
    
    # Step 1: Calculate TF-gene correlations
    print("Calculating TF-gene correlations...")
    tf_gene_correlations = pymega.get_tf_gene_correlation(
        data=multiome_data,
        tf_list=tf_list,
        gene_list=gene_list,
        method="pearson",
        min_cells=50,
        min_expression=0.1
    )
    
    # Step 2: Prepare motif and peak-gene data
    motif_matrix = np.random.random((len(tf_list), 1000))  # TFs x peaks
    peak_gene_links = pd.DataFrame({
        'peak': [f'Peak_{i}' for i in range(1000)],
        'gene': np.random.choice(gene_list, 1000),
        'score': np.random.uniform(0.1, 1.0, 1000)
    })
    
    # Step 3: Infer GRN
    print("Inferring gene regulatory network...")
    grn_network = pymega.infer_grn(
        motif_matching=motif_matrix,
        tf_gene_cor=tf_gene_correlations,
        peak_gene_links=peak_gene_links,
        integration_method="weighted_average",
        score_threshold=0.1,
        fdr_threshold=0.25
    )
    
    # Step 4: Validate network
    print("Validating network...")
    validation_results = pymega.validate_grn_network(
        grn_network=grn_network,
        validation_metrics=["topology", "statistical"],
        cross_validation=False
    )
    
    # Step 5: Predict targets
    print("Predicting TF targets...")
    target_predictions = pymega.predict_targets(
        multiome_data=multiome_data,
        tf_list=tf_list,
        prediction_method="correlation_based",
        confidence_threshold=0.1
    )
    
    # Step 6: Add target assay
    print("Adding target assay...")
    multiome_with_targets = pymega.add_target_assay(
        multiome_data=multiome_data,
        tf_targets=target_predictions,
        assay_name="inferred_targets"
    )
    
    # Results summary
    print(f"✓ Inferred {len(grn_network)} regulatory edges")
    print(f"✓ Predicted {len(target_predictions)} TF-target relationships")
    print(f"✓ Network validation completed")
    print("✓ Analysis complete!")

Performance Considerations
--------------------------

The network inference module is designed for high performance and scalability:

**Memory Management:**
- Use batch processing for large correlation matrices
- Implement sparse matrix operations where applicable
- Optimize memory usage with chunked processing

**Computational Efficiency:**
- Numba JIT compilation for correlation calculations
- Parallel processing for independent computations
- Vectorized operations for matrix calculations

**Scalability:**
- Support for datasets with 10K+ cells and genes
- Efficient algorithms for large network inference
- Memory-efficient data structures

**Best Practices:**
- Filter low-quality cells and genes before analysis
- Use appropriate correlation thresholds
- Validate results with multiple metrics
- Cross-validate with independent datasets

Troubleshooting
---------------

Common issues and solutions:

**Empty Results:**
- Check data quality and filtering parameters
- Verify TF and gene lists are present in data
- Adjust correlation and confidence thresholds

**Memory Errors:**
- Reduce batch size for correlation analysis
- Use sparse matrix operations
- Process data in chunks

**Performance Issues:**
- Enable parallel processing
- Use optimized correlation methods
- Consider data subsampling for initial analysis

**Validation Failures:**
- Check network structure and connectivity
- Verify reference data compatibility
- Adjust validation parameters

References
----------

1. **Original scMEGA Package**: R implementation of multiome network inference
2. **Aibar et al. (2017)**: SCENIC: single-cell regulatory network inference and clustering
3. **Pliner et al. (2018)**: Cicero predicts cis-regulatory DNA interactions from single-cell chromatin accessibility data
4. **Granja et al. (2021)**: ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis
5. **Marbach et al. (2012)**: Wisdom of crowds for robust gene network inference
