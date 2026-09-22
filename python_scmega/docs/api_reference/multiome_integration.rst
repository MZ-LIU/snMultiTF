Multiome Integration Module
==========================

The multiome integration module provides comprehensive tools for integrating multi-modal 
single-cell data, specifically designed for scenarios involving RNA and ATAC data from 
different sources or experimental conditions.

This module handles both paired data (same cells measured for both modalities) and unpaired 
data (different cells for RNA and ATAC), providing sophisticated algorithms for co-embedding, 
cell pairing, batch correction, and integration quality assessment.

Co-embedding
------------

.. automodule:: python_scmega.multiome_integration.coembedding
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
~~~~~~~~~~~~~

.. autofunction:: python_scmega.multiome_integration.coembedding.coembed_data
.. autofunction:: python_scmega.multiome_integration.coembedding.create_paired_object
.. autofunction:: python_scmega.multiome_integration.coembedding.weighted_nearest_neighbors
.. autofunction:: python_scmega.multiome_integration.coembedding.joint_pca

Example Usage
~~~~~~~~~~~~~

.. code-block:: python

    import python_scmega as pymega
    
    # Co-embedding for paired data (same cells in both modalities)
    coembed_result = pymega.coembed_data(
        rna_adata, atac_adata,
        rna_dims=list(range(30)),
        atac_dims=list(range(30)),
        lambda_param=0.6
    )
    
    # Create paired object with different methods
    paired_pca = pymega.create_paired_object(
        rna_adata, atac_adata,
        pairing_method="pca",
        n_components=30
    )
    
    paired_cca = pymega.create_paired_object(
        rna_adata, atac_adata,
        pairing_method="cca",
        n_components=25
    )

Cell Pairing
------------

.. automodule:: python_scmega.multiome_integration.cell_pairing
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
~~~~~~~~~~~~~

.. autofunction:: python_scmega.multiome_integration.cell_pairing.pair_cells
.. autofunction:: python_scmega.multiome_integration.cell_pairing.get_pair_list
.. autofunction:: python_scmega.multiome_integration.cell_pairing.refine_cell_pairing
.. autofunction:: python_scmega.multiome_integration.cell_pairing.validate_cell_pairing

Example Usage
~~~~~~~~~~~~~

.. code-block:: python

    # Pair cells between unpaired RNA and ATAC datasets
    paired_multiome = pymega.pair_cells(
        rna_adata, atac_adata,
        method="optimal_transport",
        max_pairs=1000
    )
    
    # Get pairing list
    pair_df = pymega.get_pair_list(paired_multiome)
    print(pair_df.head())
    
    # Validate pairing quality
    validation = pymega.validate_cell_pairing(
        paired_multiome,
        validation_genes=["CD3D", "CD8A", "CD19"]
    )
    
    # Compare different pairing methods
    comparison = pymega.compare_pairing_methods(
        rna_adata, atac_adata,
        methods=["optimal_transport", "correlation", "nearest_neighbor"]
    )

Batch Correction
----------------

.. automodule:: python_scmega.multiome_integration.batch_correction
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
~~~~~~~~~~~~~

.. autofunction:: python_scmega.multiome_integration.batch_correction.correct_batch_effects
.. autofunction:: python_scmega.multiome_integration.batch_correction.evaluate_batch_correction
.. autofunction:: python_scmega.multiome_integration.batch_correction.multiome_batch_correction

Example Usage
~~~~~~~~~~~~~

.. code-block:: python

    # Correct batch effects in single modality
    corrected_rna = pymega.correct_batch_effects(
        rna_adata,
        batch_key="batch",
        method="harmony"
    )
    
    # Correct batch effects in multiome data
    corrected_multiome = pymega.correct_batch_effects(
        multiome_data,
        batch_key="batch",
        method="harmony",
        assay="both"  # Correct both RNA and ATAC
    )
    
    # Comprehensive multiome batch correction
    corrected_multiome = pymega.multiome_batch_correction(
        multiome_data,
        batch_key="batch",
        method="harmony",
        rna_params={'theta': 2.0},
        atac_params={'theta': 1.5},
        joint_correction=True
    )
    
    # Evaluate batch correction quality
    evaluation = pymega.evaluate_batch_correction(
        corrected_rna,
        batch_key="batch",
        corrected_embedding_key="X_pca_harmony",
        cell_type_key="celltype"
    )


Complete Integration Workflow
-----------------------------

.. code-block:: python

    import python_scmega as pymega
    
    # 1. Load and preprocess data
    rna_adata = pymega.load_data("rna_data.h5ad")
    atac_adata = pymega.load_data("atac_data.h5ad")
    
    # 2. For paired data: Co-embedding
    if have_paired_data:
        integrated_data = pymega.coembed_data(
            rna_adata, atac_adata,
            rna_dims=list(range(30)),
            atac_dims=list(range(30))
        )
    
    # 3. For unpaired data: Cell pairing
    else:
        integrated_data = pymega.pair_cells(
            rna_adata, atac_adata,
            method="optimal_transport",
            max_pairs=min(len(rna_adata), len(atac_adata))
        )
    
    # 4. Batch correction (if needed)
    if "batch" in integrated_data.rna.obs.columns:
        integrated_data = pymega.correct_batch_effects(
            integrated_data,
            batch_key="batch",
            method="harmony",
            assay="both"
        )
    
    # 5. Quality assessment
    quality_metrics = pymega.evaluate_integration_quality(
        integrated_data,
        batch_key="batch",
        cell_type_key="celltype"
    )
    
    # 6. Visualization
    pymega.plot_integration_diagnostics(
        integrated_data,
        batch_key="batch",
        cell_type_key="celltype"
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
   * - ``CoembedData()``
     - ``coembed_data()``
     - coembedding
   * - ``CreatePairedObject()``
     - ``create_paired_object()``
     - coembedding
   * - ``PairCells()``
     - ``pair_cells()``
     - cell_pairing
   * - ``GetPairList()``
     - ``get_pair_list()``
     - cell_pairing
   * - ``RunHarmony()``
     - ``correct_batch_effects(method="harmony")``
     - batch_correction
   * - ``RunComBat()``
     - ``correct_batch_effects(method="combat")``
     - batch_correction

Algorithm Details
-----------------

Co-embedding Algorithms
~~~~~~~~~~~~~~~~~~~~~~~~

1. **Weighted Combination**: Combines RNA and ATAC embeddings using weighted average
2. **Canonical Correlation Analysis (CCA)**: Finds linear combinations that maximize correlation
3. **Joint PCA**: Performs PCA on concatenated feature matrices
4. **UMAP Integration**: Applies UMAP to joint embeddings

Cell Pairing Algorithms
~~~~~~~~~~~~~~~~~~~~~~~

1. **Optimal Transport**: Uses Hungarian algorithm to find optimal cell pairs
2. **Nearest Neighbor**: Pairs cells based on embedding similarity
3. **Correlation-based**: Pairs cells based on feature correlations
4. **Mutual Nearest Neighbors**: Finds mutually nearest cells

Batch Correction Methods
~~~~~~~~~~~~~~~~~~~~~~~~

1. **Harmony**: Fast integration using iterative clustering and correction
2. **ComBat**: Empirical Bayes batch correction
3. **Mutual Nearest Neighbors (MNN)**: Correction based on MNN pairs
4. **Scanorama**: Panoramic stitching of datasets
5. **Simple Centering**: Basic batch mean centering

Quality Metrics
~~~~~~~~~~~~~~~

1. **Batch Mixing**: Measures how well batches are mixed in integrated space
2. **Neighborhood Preservation**: Evaluates preservation of local neighborhoods
3. **Cell Type Preservation**: Assesses preservation of biological signal
4. **Silhouette Scores**: Measures separation of batches and cell types
5. **kBET Score**: k-nearest neighbor batch effect test
6. **Label Transfer Accuracy**: Cross-batch label transfer performance

Performance Considerations
--------------------------

- **Memory Efficiency**: Supports sparse matrices and large datasets
- **Scalability**: Optimized algorithms for datasets with 100K+ cells
- **Speed**: Uses efficient implementations from scikit-learn and scanpy
- **Parallelization**: Many algorithms support parallel computation

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
- ``harmonypy`` (for Harmony batch correction)
- ``scanorama`` (for Scanorama integration)
- ``umap-learn`` (for UMAP integration)
- ``matplotlib`` (for plotting)

See Also
--------

- :doc:`data_processing`: Data preprocessing and normalization
- :doc:`../user_guide/installation`: Installation instructions
