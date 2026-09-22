Trajectory Analysis Module
==========================

The trajectory analysis module provides comprehensive tools for single-cell trajectory 
inference and analysis, including pseudotime analysis, TF-gene correlation analysis, 
diffusion mapping, and trajectory validation.

This module handles the complete trajectory analysis workflow from initial trajectory 
inference through detailed analysis of gene expression dynamics along developmental 
or differentiation trajectories.

Trajectory Inference
--------------------

.. automodule:: python_scmega.trajectory_analysis.trajectory_inference
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
~~~~~~~~~~~~~

.. autofunction:: python_scmega.trajectory_analysis.trajectory_inference.add_trajectory
.. autofunction:: python_scmega.trajectory_analysis.trajectory_inference.select_trajectory_genes
.. autofunction:: python_scmega.trajectory_analysis.trajectory_inference.get_trajectory_summary

Example Usage
~~~~~~~~~~~~~

.. code-block:: python

    import python_scmega as pymega
    
    # Infer trajectory using R-compatible ArchR algorithm
    multiome = pymega.add_trajectory(
        multiome,
        trajectory=["0", "1", "2"],      # R: trajectory = c("0", "1", "2")
        group_by="leiden",               # R: group.by = "leiden"
        reduction="pca",                 # R: reduction = "pca"
        dims=list(range(30)),           # R: dims = 1:30
        pre_filter_quantile=0.9,        # R: pre.filter.quantile = 0.9
        post_filter_quantile=0.9,       # R: post.filter.quantile = 0.9
        dof=250,                        # R: dof = 250
        spar=1.0,                       # R: spar = 1
        name="Trajectory"               # R: name = "Trajectory"
    )
    
    # Alternative: Use cell type progression
    multiome = pymega.add_trajectory(
        multiome,
        trajectory=["early", "middle", "late"],
        group_by="cell_type",
        reduction="pca",
        dims=list(range(30)),
        name="CellTypeTrajectory"
    )
    
    # Alternative: Use stricter filtering parameters
    multiome = pymega.add_trajectory(
        multiome,
        trajectory=["0", "1", "2"],
        group_by="leiden",
        reduction="pca",
        dims=list(range(30)),
        pre_filter_quantile=0.95,       # Stricter filtering
        post_filter_quantile=0.95,
        dof=200,                        # Lower DOF for faster computation
        name="StrictTrajectory"
    )
    
    # Select genes that vary along trajectory
    trajectory_genes = pymega.select_trajectory_genes(
        multiome,
        trajectory_name="Trajectory",
        n_genes=1000,
        method="correlation"
    )

Pseudotime Analysis
-------------------

.. automodule:: python_scmega.trajectory_analysis.pseudotime_analysis
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
~~~~~~~~~~~~~

.. autofunction:: python_scmega.trajectory_analysis.pseudotime_analysis.get_trajectory_data
.. autofunction:: python_scmega.trajectory_analysis.pseudotime_analysis.analyze_pseudotime_dynamics
.. autofunction:: python_scmega.trajectory_analysis.pseudotime_analysis.fit_trajectory_curves
.. autofunction:: python_scmega.trajectory_analysis.pseudotime_analysis.identify_trajectory_markers
.. autofunction:: python_scmega.trajectory_analysis.pseudotime_analysis.plot_trajectory_heatmap

Example Usage
~~~~~~~~~~~~~

.. code-block:: python

    # Extract trajectory data for analysis
    trajectory_df = pymega.get_trajectory_data(
        multiome,
        trajectory_name="Trajectory",
        genes=["CD34", "CD38", "CD14"]
    )
    
    # Analyze gene expression dynamics along pseudotime
    dynamics = pymega.analyze_pseudotime_dynamics(
        trajectory_df,
        genes=["CD34", "CD38", "CD14"],
        n_bins=50,
        smooth_method="spline"
    )
    
    # Fit different curve types to gene expression
    curves = pymega.fit_trajectory_curves(
        trajectory_df,
        genes=["CD34", "CD14"],
        curve_types=["linear", "polynomial", "sigmoid"]
    )
    
    # Identify genes with significant trajectory dynamics
    markers = pymega.identify_trajectory_markers(
        trajectory_df,
        method="correlation",
        min_correlation=0.4
    )
    
    # Create heatmap of gene expression along trajectory
    heatmap_info = pymega.plot_trajectory_heatmap(
        trajectory_df,
        genes=markers.head(50)['gene'].tolist(),
        n_bins=100,
        scale_genes=True
    )

TF-Gene Correlation Analysis
----------------------------

.. automodule:: python_scmega.trajectory_analysis.getcorrelation
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
~~~~~~~~~~~~~

.. autofunction:: python_scmega.trajectory_analysis.getcorrelation.calculate_trajectory_correlation
.. autofunction:: python_scmega.trajectory_analysis.getcorrelation.get_tf_gene_correlations
.. autofunction:: python_scmega.trajectory_analysis.getcorrelation.identify_trajectory_tfs
.. autofunction:: python_scmega.trajectory_analysis.getcorrelation.analyze_tf_target_networks
.. autofunction:: python_scmega.trajectory_analysis.getcorrelation.plot_correlation_heatmap

Example Usage
~~~~~~~~~~~~~

.. code-block:: python

    # Calculate TF-gene correlations along trajectory
    correlations = pymega.calculate_trajectory_correlation(
        multiome,
        trajectory_name="Trajectory",
        tf_genes=["GATA1", "PU1", "CEBPA"],
        target_genes=["CD34", "CD14", "CD16"],
        method="pearson"
    )
    
    # Get significant TF-gene correlations
    tf_correlations = pymega.get_tf_gene_correlations(
        multiome,
        tf_list=["GATA1", "PU1"],
        correlation_threshold=0.4
    )
    
    # Identify TFs with significant regulatory activity
    trajectory_tfs = pymega.identify_trajectory_tfs(
        multiome,
        min_targets=10,
        correlation_threshold=0.3
    )
    
    # Analyze TF-target gene regulatory networks
    networks = pymega.analyze_tf_target_networks(
        multiome,
        top_tfs=["GATA1", "PU1", "CEBPA"],
        correlation_threshold=0.4
    )
    
    # Create correlation heatmap
    plot_info = pymega.plot_correlation_heatmap(
        correlations,
        top_n=50
    )

Diffusion Mapping
-----------------

.. automodule:: python_scmega.trajectory_analysis.diffusion_mapping
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
~~~~~~~~~~~~~

.. autofunction:: python_scmega.trajectory_analysis.diffusion_mapping.run_diffusion_map
.. autofunction:: python_scmega.trajectory_analysis.diffusion_mapping.calculate_diffusion_pseudotime
.. autofunction:: python_scmega.trajectory_analysis.diffusion_mapping.diffusion_components_analysis
.. autofunction:: python_scmega.trajectory_analysis.diffusion_mapping.plot_diffusion_map

Example Usage
~~~~~~~~~~~~~

.. code-block:: python

    # Run diffusion mapping
    multiome = pymega.run_diffusion_map(
        multiome,
        reduction="pca",
        dims=list(range(30)),
        n_neighbors=30,
        n_components=10
    )
    
    # Calculate diffusion pseudotime
    multiome = pymega.calculate_diffusion_pseudotime(
        multiome,
        start_cluster="0",
        diffusion_component=1
    )
    
    # Analyze diffusion components
    component_analysis = pymega.diffusion_components_analysis(
        multiome,
        n_components=5
    )
    
    # Plot diffusion map
    plot_info = pymega.plot_diffusion_map(
        multiome,
        components=(1, 2),
        color_by="leiden"
    )

Trajectory Validation
---------------------

.. automodule:: python_scmega.trajectory_analysis.trajectory_validation
   :members:
   :undoc-members:
   :show-inheritance:

Key Functions
~~~~~~~~~~~~~

.. autofunction:: python_scmega.trajectory_analysis.trajectory_validation.validate_trajectory_quality
.. autofunction:: python_scmega.trajectory_analysis.trajectory_validation.assess_trajectory_robustness
.. autofunction:: python_scmega.trajectory_analysis.trajectory_validation.validate_pseudotime_ordering

Example Usage
~~~~~~~~~~~~~

.. code-block:: python

    # Comprehensive trajectory quality validation
    validation = pymega.validate_trajectory_quality(
        multiome,
        trajectory_name="Trajectory",
        validation_genes=["CD34", "CD38", "CD14"],
        known_markers={
            "early": ["CD34", "KIT"],
            "late": ["CD14", "CD16"]
        }
    )
    
    # Assess trajectory robustness
    robustness = pymega.assess_trajectory_robustness(
        multiome,
        trajectory_function=pymega.add_trajectory,
        parameter_ranges={
            'method': ['monocle', 'slingshot'],
            'n_neighbors': [15, 30, 50]
        },
        n_bootstrap=10
    )
    
    # Validate pseudotime ordering
    ordering_validation = pymega.validate_pseudotime_ordering(
        multiome,
        reference_genes=["CD34", "CD38", "CD14"],
        time_points=["day0", "day3", "day7"]
    )

Complete Trajectory Analysis Workflow
-------------------------------------

.. code-block:: python

    import python_scmega as pymega
    
    # 1. Load and preprocess data
    multiome = pymega.load_10x_multiome("data.h5")
    multiome = pymega.quality_control_multiome(multiome)
    multiome = pymega.normalize_multiome(multiome)
    
    # 2. Supply cluster labels from an upstream workflow
    if "leiden" not in multiome.obs:
        raise KeyError("Expected a 'leiden' column in multiome.obs")
    
    # 3. Trajectory inference (R-compatible ArchR algorithm)
    multiome = pymega.add_trajectory(
        multiome,
        trajectory=["0", "1", "2"],
        group_by="leiden",
        reduction="pca",
        dims=list(range(30)),
        name="Trajectory"
    )
    
    # 4. Extract and analyze trajectory data
    trajectory_df = pymega.get_trajectory_data(multiome)
    
    dynamics = pymega.analyze_pseudotime_dynamics(
        trajectory_df,
        n_bins=50
    )
    
    markers = pymega.identify_trajectory_markers(
        trajectory_df,
        method="correlation"
    )
    
    # 5. TF-gene correlation analysis
    correlations = pymega.calculate_trajectory_correlation(
        multiome,
        tf_genes=["GATA1", "PU1", "CEBPA"]
    )
    
    # 6. Diffusion mapping
    multiome = pymega.run_diffusion_map(
        multiome,
        n_neighbors=30,
        n_components=10
    )
    
    # 7. Validate trajectory quality
    validation = pymega.validate_trajectory_quality(
        multiome,
        known_markers={
            "early": ["CD34", "KIT"],
            "late": ["CD14", "CD16"]
        }
    )
    
    # 8. Visualization
    pymega.plot_trajectory_heatmap(trajectory_df, genes=markers.head(50)['gene'])
    pymega.plot_correlation_heatmap(correlations)
    pymega.plot_diffusion_map(multiome, color_by="leiden")

R Equivalence
-------------

This module provides Python equivalents for key R functions:

.. list-table:: R to Python Function Mapping
   :widths: 40 40 20
   :header-rows: 1

   * - R Function
     - PyMEGA Equivalent
     - Module
   * - ``AddTrajectory()``
     - ``add_trajectory()``
     - trajectory_inference
   * - ``GetTrajectory()``
     - ``get_trajectory_data()``
     - pseudotime_analysis
   * - ``GetCorrelation()``
     - ``calculate_trajectory_correlation()``
     - getcorrelation
   * - ``RunDiffusionMap()``
     - ``run_diffusion_map()``
     - diffusion_mapping

Algorithm Details
-----------------

Trajectory Inference Methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Monocle3-style**: Builds minimum spanning tree between cluster centers
2. **Slingshot-style**: Fits principal curves through ordered clusters  
3. **Diffusion-based**: Uses diffusion distances for trajectory ordering

Pseudotime Analysis
~~~~~~~~~~~~~~~~~~~

1. **Dynamics Analysis**: Fits smooth curves to gene expression along pseudotime
2. **Curve Fitting**: Supports linear, polynomial, and sigmoid curve types
3. **Marker Identification**: Uses correlation and trend analysis methods

TF-Gene Correlations
~~~~~~~~~~~~~~~~~~~~

1. **Correlation Methods**: Pearson, Spearman, and mutual information
2. **Binned Analysis**: Correlation analysis within pseudotime bins
3. **Network Analysis**: TF-target network reconstruction and analysis

Diffusion Mapping
~~~~~~~~~~~~~~~~~

1. **Graph Construction**: k-nearest neighbor graph with Gaussian kernel
2. **Eigendecomposition**: Spectral analysis of diffusion operator
3. **Pseudotime Calculation**: Diffusion distance from start point

Validation Metrics
~~~~~~~~~~~~~~~~~~

1. **Topology**: Path smoothness, straightness, density variation
2. **Smoothness**: Gene expression smoothness along trajectory
3. **Expression**: Correlation with pseudotime, monotonicity
4. **Markers**: Known marker gene validation
5. **Clustering**: Consistency with clustering results

Performance Considerations
--------------------------

- **Memory Efficiency**: Supports sparse matrices and large datasets
- **Scalability**: Optimized algorithms for datasets with 100K+ cells
- **Speed**: Efficient implementations using numpy and scipy
- **Robustness**: Multiple methods and validation approaches

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
- ``matplotlib`` (for plotting)
- ``seaborn`` (for enhanced visualization)
- ``statsmodels`` (for statistical tests)

See Also
--------

- :doc:`data_processing`: Data preprocessing and normalization
- :doc:`multiome_integration`: Multi-modal data integration
- :doc:`../user_guide/installation`: Installation instructions
