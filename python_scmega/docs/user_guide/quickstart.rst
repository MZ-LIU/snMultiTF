Quick Start Guide
=================

This guide will get you up and running with PyMEGA in just a few minutes. We'll walk through a complete analysis workflow using example data.

Overview
--------

A typical PyMEGA analysis follows these steps:

1. **Data Loading** - Load 10X multiome data or create MultiomeData objects
2. **Quality Control** - Filter low-quality cells and features
3. **Normalization** - Normalize RNA and ATAC data
4. **Integration** - Integrate RNA and ATAC modalities
5. **Trajectory Analysis** - Infer developmental trajectories
6. **Network Inference** - Build gene regulatory networks
7. **Visualization** - Create publication-ready plots

Example Analysis
----------------

Let's walk through a complete analysis using simulated data:

Loading Data
~~~~~~~~~~~~

.. code-block:: python

   import pymega
   import numpy as np
   import pandas as pd
   
   # For this example, we'll create synthetic data
   # In practice, you would load real 10X data:
   # multiome = pymega.load_10x_multiome("path/to/10x/data")
   
   # Create synthetic multiome data
   n_cells, n_genes, n_peaks = 2000, 3000, 5000
   
   # RNA data (log-normalized counts)
   rna_data = np.random.lognormal(0, 1, (n_cells, n_genes))
   rna_genes = [f"Gene_{i:04d}" for i in range(n_genes)]
   
   # ATAC data (binary accessibility)
   atac_data = np.random.binomial(1, 0.1, (n_cells, n_peaks))
   atac_peaks = [f"Peak_{i:04d}" for i in range(n_peaks)]
   
   # Cell metadata
   cell_ids = [f"Cell_{i:04d}" for i in range(n_cells)]
   cell_types = np.random.choice(['TypeA', 'TypeB', 'TypeC'], n_cells)
   
   # Create MultiomeData object
   multiome = pymega.MultiomeData(
       rna_data=pd.DataFrame(rna_data, index=cell_ids, columns=rna_genes),
       atac_data=pd.DataFrame(atac_data, index=cell_ids, columns=atac_peaks),
       obs=pd.DataFrame({'cell_type': cell_types}, index=cell_ids)
   )
   
   print(f"Loaded data: {multiome.n_obs} cells, {multiome.n_rna_vars} genes, {multiome.n_atac_vars} peaks")

Quality Control
~~~~~~~~~~~~~~~

.. code-block:: python

   # Calculate QC metrics
   multiome = pymega.quality_control_multiome(
       multiome,
       min_genes=100,        # Minimum genes per cell
       min_cells=10,         # Minimum cells per gene
       max_mito_pct=20,      # Maximum mitochondrial percentage
       min_peaks=50,         # Minimum peaks per cell
       max_peaks=10000       # Maximum peaks per cell
   )
   
   print(f"After QC: {multiome.n_obs} cells, {multiome.n_rna_vars} genes, {multiome.n_atac_vars} peaks")
   
   # Visualize QC metrics
   pymega.plot_qc_metrics(multiome)

Normalization and Dimensionality Reduction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Normalize multiome data
   multiome = pymega.normalize_multiome(
       multiome,
       rna_method="log1p",           # Log-normalize RNA
       atac_method="tfidf",          # TF-IDF normalize ATAC
       rna_n_components=50,          # RNA PCA components
       atac_n_components=50          # ATAC SVD components
   )
   
   print("✓ Normalization completed")
   
   # The normalized data and dimensionality reductions are now stored in multiome

Cell-Type Metadata
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Cell-type labels must be supplied by the input data or an upstream
   # annotation workflow. The legacy annotation API is no longer bundled.
   if 'cell_type' not in multiome.obs:
       raise KeyError("Expected a 'cell_type' column in multiome.obs")

Multiome Integration
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Co-embed RNA and ATAC data
   multiome = pymega.coembed_data(
       multiome,
       method="cca",              # Canonical Correlation Analysis
       n_components=30,           # Number of co-embedding components
       resolution=0.5             # Clustering resolution
   )
   
   # Pair cells between modalities (for unpaired data)
   if not multiome.is_paired:
       multiome = pymega.pair_cells(multiome, method="knn", k=5)
   
   print("✓ Multiome integration completed")

Trajectory Analysis
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Infer developmental trajectories using R-compatible ArchR algorithm
   multiome = pymega.add_trajectory(
       multiome,
       trajectory=["TypeA", "TypeB", "TypeC"],  # Trajectory path (R: trajectory = c(...))
       group_by="cell_type",                     # Group column (R: group.by = "cell_type")
       reduction="pca",                          # Dimensionality reduction
       dims=list(range(30)),                     # Use first 30 PCs
       pre_filter_quantile=0.9,                  # Pre-filtering quantile
       post_filter_quantile=0.9,                 # Post-filtering quantile
       dof=250,                                  # Degrees of freedom for spline
       spar=1.0,                                 # Spline smoothing parameter
       name="Trajectory"                         # Trajectory name
   )
   
   # Get trajectory data (R GetTrajectory equivalent)
   trajectory_result = pymega.get_trajectory_data(
       multiome,
       trajectory_name="Trajectory",
       group_every=1,
       log2_norm=True,
       smooth_window=11
   )
   
   print(f"✓ Trajectory inferred with {trajectory_result['group_matrix'].shape[1]} bins")
   
   # Visualize trajectories
   pymega.plot_trajectory(multiome, color_by='pseudotime')

Gene Regulatory Network Inference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Define transcription factors and target genes
   # In practice, you would use known TF lists
   tf_genes = [gene for gene in multiome.rna_var_names if 'TF' in gene][:50]
   target_genes = multiome.rna_var_names[:500].tolist()
   
   if len(tf_genes) == 0:
       # For demo, use first 50 genes as TFs
       tf_genes = multiome.rna_var_names[:50].tolist()
   
   # Infer gene regulatory network
   grn_network = pymega.infer_grn(
       multiome,
       tf_list=tf_genes,
       target_genes=target_genes,
       method="correlation",       # Network inference method
       min_correlation=0.3,        # Minimum correlation threshold
       fdr_threshold=0.05          # False discovery rate
   )
   
   print(f"✓ Inferred GRN with {len(grn_network)} TF-gene interactions")
   
   # Link peaks to genes
   peak_gene_links = pymega.link_peaks_to_genes(
       multiome,
       distance_threshold=100000,   # 100kb distance threshold
       correlation_threshold=0.2    # Minimum correlation
   )
   
   print(f"✓ Identified {len(peak_gene_links)} peak-gene links")

Visualization
~~~~~~~~~~~~~

.. code-block:: python

   # Create comprehensive visualizations
   
   # 1. Network visualization
   fig_network = pymega.plot_grn_network(
       grn_network.head(100),  # Top 100 interactions
       layout="spring",
       node_size_col="degree",
       edge_width_col="correlation",
       interactive=False
   )
   
   # 2. Trajectory heatmap
   fig_trajectory = pymega.plot_trajectory_heatmap(
       multiome.rna_data,
       trajectory_data,
       gene_list=marker_results['gene'].head(20).tolist(),
       n_bins=50
   )
   
   # 3. Expression heatmap
   fig_expression = pymega.plot_expression_heatmap(
       multiome.rna_data.iloc[:100, :50],  # Subset for visualization
       cluster_genes=True,
       cluster_cells=True
   )
   
   # 4. Cell type proportions
   fig_proportions = pymega.plot_cell_proportions(
       multiome.obs,
       group_col='cell_type',
       cell_type_col='cell_type'
   )
   
   print("✓ Generated comprehensive visualizations")

Interactive Analysis
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Create interactive visualizations for exploration
   
   # Interactive network plot
   interactive_network = pymega.plot_grn_network(
       grn_network.head(50),
       interactive=True,
       layout="force"
   )
   # interactive_network.show()  # Uncomment to display
   
   # Interactive trajectory plot
   interactive_trajectory = pymega.create_interactive_trajectory(
       trajectory_data,
       multiome.obsm['X_umap']  # UMAP coordinates
   )
   # interactive_trajectory.show()  # Uncomment to display
   
   print("✓ Created interactive visualizations")

Save Results
~~~~~~~~~~~~

.. code-block:: python

   # Save analysis results
   
   # Save the multiome object
   multiome.write_h5ad("analysis_results.h5ad")
   
   # Save network data
   grn_network.to_csv("grn_network.csv", index=False)
   peak_gene_links.to_csv("peak_gene_links.csv", index=False)
   
   # Save plots
   fig_network.savefig("grn_network.png", dpi=300, bbox_inches='tight')
   fig_trajectory.savefig("trajectory_heatmap.png", dpi=300, bbox_inches='tight')
   
   print("✓ Results saved successfully")

Complete Example Script
-----------------------

Here's the complete example as a single script:

.. code-block:: python

   #!/usr/bin/env python3
   """
   PyMEGA Quick Start Example
   
   A complete single-cell multiome analysis workflow.
   """
   
   import pymega
   import numpy as np
   import pandas as pd
   import warnings
   warnings.filterwarnings('ignore')
   
   def main():
       print("PyMEGA Quick Start Analysis")
       print("=" * 40)
       
       # 1. Create/Load Data
       print("\\n1. Loading data...")
       multiome = create_example_data()
       
       # 2. Quality Control
       print("\\n2. Quality control...")
       multiome = pymega.quality_control_multiome(multiome)
       
       # 3. Normalization
       print("\\n3. Normalization...")
       multiome = pymega.normalize_multiome(multiome)
       
       # 4. Integration
       print("\\n4. Multiome integration...")
       multiome = pymega.coembed_data(multiome)
       
       # 5. Trajectory Analysis
       print("\\n5. Trajectory analysis...")
       multiome = pymega.add_trajectory(multiome)
       
       # 6. Network Inference
       print("\\n6. Network inference...")
       grn_network = pymega.infer_grn(
           multiome,
           tf_list=multiome.rna_var_names[:50].tolist(),
           target_genes=multiome.rna_var_names[:200].tolist()
       )
       
       # 7. Visualization
       print("\\n7. Creating visualizations...")
       create_visualizations(multiome, grn_network)
       
       print("\\n✓ Analysis completed successfully!")
       return multiome, grn_network
   
   def create_example_data():
       """Create synthetic multiome data for demonstration."""
       np.random.seed(42)
       n_cells, n_genes, n_peaks = 1000, 2000, 3000
       
       # Generate data with some structure
       rna_data = np.random.lognormal(0, 1, (n_cells, n_genes))
       atac_data = np.random.binomial(1, 0.1, (n_cells, n_peaks))
       
       # Create identifiers
       cell_ids = [f"Cell_{i:04d}" for i in range(n_cells)]
       rna_genes = [f"Gene_{i:04d}" for i in range(n_genes)]
       atac_peaks = [f"Peak_{i:04d}" for i in range(n_peaks)]
       cell_types = np.random.choice(['Stem', 'Progenitor', 'Differentiated'], n_cells)
       
       return pymega.MultiomeData(
           rna_data=pd.DataFrame(rna_data, index=cell_ids, columns=rna_genes),
           atac_data=pd.DataFrame(atac_data, index=cell_ids, columns=atac_peaks),
           obs=pd.DataFrame({'cell_type': cell_types}, index=cell_ids)
       )
   
   def create_visualizations(multiome, grn_network):
       """Create example visualizations."""
       try:
           # Network plot
           pymega.plot_grn_network(grn_network.head(50), interactive=False)
           
           # Expression heatmap
           pymega.plot_expression_heatmap(
               multiome.rna_data.iloc[:50, :30],
               cluster_genes=True
           )
           
           print("  ✓ Visualizations created")
       except Exception as e:
           print(f"  ⚠ Visualization error: {e}")
   
   if __name__ == "__main__":
       multiome, grn_network = main()

Next Steps
----------

Now that you've completed your first PyMEGA analysis, here are some next steps:

**Explore Your Data:**
- Use interactive visualizations to explore your results
- Try different parameter settings for network inference
- Experiment with different trajectory methods

**Learn More:**
- :doc:`../tutorials/basic_workflow` - More detailed tutorial

**Optimize Performance:**
- :doc:`performance_optimization` - Speed up your analysis
- Use parallel processing for large datasets
- Optimize memory usage for better performance

**Get Help:**
- Check the :doc:`../api_reference/index` for function details
- Visit our `GitHub discussions <https://github.com/pymega/pymega/discussions>`_ for questions

Happy analyzing! 🧬
