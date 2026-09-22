Basic Workflow Tutorial
=======================

This tutorial walks through a complete PyMEGA analysis workflow using real-world single-cell multiome data. We'll cover all the essential steps from data loading to network inference and visualization.

Overview
--------

In this tutorial, you'll learn how to:

1. Load and explore 10X multiome data
2. Perform quality control and filtering
3. Normalize RNA and ATAC data
4. Annotate cell types using reference data
5. Integrate RNA and ATAC modalities
6. Infer developmental trajectories
7. Build gene regulatory networks
8. Create publication-ready visualizations

Dataset
-------

We'll use a publicly available dataset of mouse brain development containing:
- **Cells**: ~5,000 cells
- **Genes**: ~20,000 genes
- **Peaks**: ~50,000 ATAC peaks
- **Cell Types**: Neural stem cells, neurons, glia

You can download the example dataset from: https://example.com/mouse_brain_multiome.h5

Step 1: Environment Setup
-------------------------

First, let's set up our analysis environment:

.. code-block:: python

   import pymega
   import numpy as np
   import pandas as pd
   import matplotlib.pyplot as plt
   import seaborn as sns
   from pathlib import Path
   import warnings
   
   # Configure display options
   warnings.filterwarnings('ignore')
   plt.style.use('seaborn-v0_8')
   pd.set_option('display.max_columns', None)
   
   # Set random seed for reproducibility
   np.random.seed(42)
   
   print(f"PyMEGA version: {pymega.__version__}")
   print(f"Analysis started at: {pd.Timestamp.now()}")

Step 2: Data Loading
--------------------

Load the 10X multiome data:

.. code-block:: python

   # Define data paths
   data_dir = Path("data/mouse_brain_multiome/")
   
   # Load 10X multiome data
   multiome = pymega.load_10x_multiome(
       data_dir,
       rna_h5="filtered_feature_bc_matrix.h5",
       atac_fragments="atac_fragments.tsv.gz",
       atac_peaks="peaks.bed",
       genome="mm10"
   )
   
   # Basic information about the dataset
   print(f"Dataset overview:")
   print(f"  Cells: {multiome.n_obs:,}")
   print(f"  RNA features: {multiome.n_rna_vars:,}")
   print(f"  ATAC features: {multiome.n_atac_vars:,}")
   print(f"  Total features: {multiome.n_vars:,}")
   
   # Examine the data structure
   print(f"\\nRNA data shape: {multiome.rna_data.shape}")
   print(f"ATAC data shape: {multiome.atac_data.shape}")
   print(f"Observation metadata: {list(multiome.obs.columns)}")

Let's explore the raw data:

.. code-block:: python

   # Plot basic statistics
   fig, axes = plt.subplots(2, 2, figsize=(12, 10))
   
   # Total RNA counts per cell
   total_rna_counts = multiome.rna_data.sum(axis=1)
   axes[0, 0].hist(total_rna_counts, bins=50, alpha=0.7)
   axes[0, 0].set_xlabel('Total RNA counts')
   axes[0, 0].set_ylabel('Number of cells')
   axes[0, 0].set_title('RNA Count Distribution')
   
   # Number of detected genes per cell
   n_genes_detected = (multiome.rna_data > 0).sum(axis=1)
   axes[0, 1].hist(n_genes_detected, bins=50, alpha=0.7, color='orange')
   axes[0, 1].set_xlabel('Number of genes detected')
   axes[0, 1].set_ylabel('Number of cells')
   axes[0, 1].set_title('Gene Detection per Cell')
   
   # Total ATAC counts per cell
   total_atac_counts = multiome.atac_data.sum(axis=1)
   axes[1, 0].hist(total_atac_counts, bins=50, alpha=0.7, color='green')
   axes[1, 0].set_xlabel('Total ATAC counts')
   axes[1, 0].set_ylabel('Number of cells')
   axes[1, 0].set_title('ATAC Count Distribution')
   
   # Number of accessible peaks per cell
   n_peaks_accessible = (multiome.atac_data > 0).sum(axis=1)
   axes[1, 1].hist(n_peaks_accessible, bins=50, alpha=0.7, color='red')
   axes[1, 1].set_xlabel('Number of accessible peaks')
   axes[1, 1].set_ylabel('Number of cells')
   axes[1, 1].set_title('Peak Accessibility per Cell')
   
   plt.tight_layout()
   plt.show()

Step 3: Quality Control
-----------------------

Now let's perform quality control to filter low-quality cells and features:

.. code-block:: python

   # Calculate QC metrics
   print("Calculating QC metrics...")
   
   # Add mitochondrial gene information
   multiome.var['mt'] = multiome.rna_var_names.str.startswith('mt-')
   
   # Calculate QC metrics for each cell
   multiome.obs['total_rna_counts'] = multiome.rna_data.sum(axis=1)
   multiome.obs['n_genes_detected'] = (multiome.rna_data > 0).sum(axis=1)
   multiome.obs['pct_mt'] = (
       multiome.rna_data.loc[:, multiome.var['mt']].sum(axis=1) / 
       multiome.obs['total_rna_counts'] * 100
   )
   
   multiome.obs['total_atac_counts'] = multiome.atac_data.sum(axis=1)
   multiome.obs['n_peaks_accessible'] = (multiome.atac_data > 0).sum(axis=1)
   
   # Display QC summary
   print("\\nQC Metrics Summary:")
   print(multiome.obs[['total_rna_counts', 'n_genes_detected', 'pct_mt', 
                      'total_atac_counts', 'n_peaks_accessible']].describe())

Visualize QC metrics:

.. code-block:: python

   # Create QC plots
   fig, axes = plt.subplots(2, 3, figsize=(15, 10))
   
   # Violin plots for key metrics
   metrics = ['total_rna_counts', 'n_genes_detected', 'pct_mt', 
              'total_atac_counts', 'n_peaks_accessible']
   
   for i, metric in enumerate(metrics):
       row, col = i // 3, i % 3
       sns.violinplot(y=multiome.obs[metric], ax=axes[row, col])
       axes[row, col].set_title(f'{metric}')
   
   # Scatter plot: RNA vs ATAC counts
   axes[1, 2].scatter(multiome.obs['total_rna_counts'], 
                      multiome.obs['total_atac_counts'], 
                      alpha=0.5, s=1)
   axes[1, 2].set_xlabel('Total RNA counts')
   axes[1, 2].set_ylabel('Total ATAC counts')
   axes[1, 2].set_title('RNA vs ATAC counts')
   
   plt.tight_layout()
   plt.show()

Apply quality control filters:

.. code-block:: python

   # Define QC thresholds
   qc_params = {
       'min_genes': 500,           # Minimum genes per cell
       'max_genes': 7000,          # Maximum genes per cell
       'min_rna_counts': 1000,     # Minimum RNA counts per cell
       'max_rna_counts': 50000,    # Maximum RNA counts per cell
       'max_mito_pct': 20,         # Maximum mitochondrial percentage
       'min_peaks': 200,           # Minimum peaks per cell
       'max_peaks': 20000,         # Maximum peaks per cell
       'min_cells_gene': 10,       # Minimum cells per gene
       'min_cells_peak': 10        # Minimum cells per peak
   }
   
   print(f"Before QC: {multiome.n_obs} cells, {multiome.n_rna_vars} genes, {multiome.n_atac_vars} peaks")
   
   # Apply QC filters
   multiome = pymega.quality_control_multiome(multiome, **qc_params)
   
   print(f"After QC: {multiome.n_obs} cells, {multiome.n_rna_vars} genes, {multiome.n_atac_vars} peaks")
   
   # Calculate filtering statistics
   cells_retained = multiome.n_obs / len(multiome.obs) * 100
   genes_retained = multiome.n_rna_vars / multiome.rna_data.shape[1] * 100
   peaks_retained = multiome.n_atac_vars / multiome.atac_data.shape[1] * 100
   
   print(f"\\nFiltering summary:")
   print(f"  Cells retained: {cells_retained:.1f}%")
   print(f"  Genes retained: {genes_retained:.1f}%")
   print(f"  Peaks retained: {peaks_retained:.1f}%")

Step 4: Normalization and Dimensionality Reduction
---------------------------------------------------

Normalize the data and perform dimensionality reduction:

.. code-block:: python

   # Normalize multiome data
   print("Normalizing data...")
   
   multiome = pymega.normalize_multiome(
       multiome,
       rna_method="log1p",              # Log-normalize RNA data
       atac_method="tfidf",             # TF-IDF normalize ATAC data
       rna_scale=True,                  # Scale RNA data
       rna_n_components=50,             # Number of PCA components for RNA
       atac_n_components=50,            # Number of SVD components for ATAC
       rna_highly_variable=True,        # Find highly variable genes
       atac_highly_variable=True       # Find highly variable peaks
   )
   
   print("✓ Normalization completed")
   
   # Check the results
   print(f"\\nRNA PCA components: {multiome.obsm['X_rna_pca'].shape}")
   print(f"ATAC SVD components: {multiome.obsm['X_atac_svd'].shape}")
   print(f"Highly variable genes: {multiome.var['highly_variable_rna'].sum()}")
   print(f"Highly variable peaks: {multiome.var['highly_variable_atac'].sum()}")

Visualize the normalized data:

.. code-block:: python

   # Create UMAP embeddings for visualization
   import umap
   
   # RNA UMAP
   rna_umap = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
   multiome.obsm['X_rna_umap'] = rna_umap.fit_transform(multiome.obsm['X_rna_pca'])
   
   # ATAC UMAP
   atac_umap = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
   multiome.obsm['X_atac_umap'] = atac_umap.fit_transform(multiome.obsm['X_atac_svd'])
   
   # Plot UMAPs
   fig, axes = plt.subplots(1, 2, figsize=(12, 5))
   
   # RNA UMAP
   scatter = axes[0].scatter(multiome.obsm['X_rna_umap'][:, 0], 
                            multiome.obsm['X_rna_umap'][:, 1],
                            c=multiome.obs['total_rna_counts'], 
                            cmap='viridis', s=1)
   axes[0].set_title('RNA UMAP')
   axes[0].set_xlabel('UMAP 1')
   axes[0].set_ylabel('UMAP 2')
   plt.colorbar(scatter, ax=axes[0], label='Total RNA counts')
   
   # ATAC UMAP
   scatter = axes[1].scatter(multiome.obsm['X_atac_umap'][:, 0], 
                            multiome.obsm['X_atac_umap'][:, 1],
                            c=multiome.obs['total_atac_counts'], 
                            cmap='plasma', s=1)
   axes[1].set_title('ATAC UMAP')
   axes[1].set_xlabel('UMAP 1')
   axes[1].set_ylabel('UMAP 2')
   plt.colorbar(scatter, ax=axes[1], label='Total ATAC counts')
   
   plt.tight_layout()
   plt.show()

Step 5: Cell Type Annotation
-----------------------------

Perform clustering and cell type annotation:

.. code-block:: python

   # Perform clustering on RNA data
   print("Performing clustering...")
   
   # Cluster labels must be supplied by an upstream workflow.
   if 'clusters' not in multiome.obs:
       raise KeyError("Expected a 'clusters' column in multiome.obs")
   
   print(f"Identified {len(multiome.obs['clusters'].unique())} clusters")
   
   # Marker discovery and cell-type assignment are expected to be performed by
   # an upstream workflow; the legacy annotation API is no longer bundled.

Visualize clusters and marker genes:

.. code-block:: python

   # Plot clusters on UMAP
   fig, axes = plt.subplots(1, 3, figsize=(18, 5))
   
   # Clusters on RNA UMAP
   scatter = axes[0].scatter(multiome.obsm['X_rna_umap'][:, 0], 
                            multiome.obsm['X_rna_umap'][:, 1],
                            c=pd.Categorical(multiome.obs['clusters']).codes, 
                            cmap='tab20', s=1)
   axes[0].set_title('RNA UMAP - Clusters')
   axes[0].set_xlabel('UMAP 1')
   axes[0].set_ylabel('UMAP 2')
   
   # Clusters on ATAC UMAP
   scatter = axes[1].scatter(multiome.obsm['X_atac_umap'][:, 0], 
                            multiome.obsm['X_atac_umap'][:, 1],
                            c=pd.Categorical(multiome.obs['clusters']).codes, 
                            cmap='tab20', s=1)
   axes[1].set_title('ATAC UMAP - Clusters')
   axes[1].set_xlabel('UMAP 1')
   axes[1].set_ylabel('UMAP 2')
   
   # Example marker gene expression
   if 'Neurod1' in multiome.rna_var_names:
       gene_expr = multiome.rna_data['Neurod1']
       scatter = axes[2].scatter(multiome.obsm['X_rna_umap'][:, 0], 
                                multiome.obsm['X_rna_umap'][:, 1],
                                c=gene_expr, cmap='Reds', s=1)
       axes[2].set_title('Neurod1 Expression')
       plt.colorbar(scatter, ax=axes[2])
   else:
       axes[2].text(0.5, 0.5, 'Neurod1 not found\\nin dataset', 
                    ha='center', va='center', transform=axes[2].transAxes)
       axes[2].set_title('Marker Gene Expression')
   
   plt.tight_layout()
   plt.show()

Assign cell type labels based on marker genes:

.. code-block:: python

   # Manual cell type annotation based on known markers
   # Cell-type labels can also be supplied by an upstream annotation workflow.
   
   cell_type_mapping = {
       '0': 'Neural_Stem_Cells',
       '1': 'Neurons',
       '2': 'Oligodendrocytes', 
       '3': 'Astrocytes',
       '4': 'Microglia',
       '5': 'Neural_Progenitors',
       '6': 'Endothelial_Cells'
   }
   
   # Map clusters to cell types
   multiome.obs['cell_type'] = multiome.obs['clusters'].map(cell_type_mapping)
   multiome.obs['cell_type'] = multiome.obs['cell_type'].fillna('Unknown')
   
   print("Cell type distribution:")
   print(multiome.obs['cell_type'].value_counts())
   
   # Plot cell types
   fig, ax = plt.subplots(figsize=(10, 8))
   scatter = ax.scatter(multiome.obsm['X_rna_umap'][:, 0], 
                       multiome.obsm['X_rna_umap'][:, 1],
                       c=pd.Categorical(multiome.obs['cell_type']).codes, 
                       cmap='Set3', s=1)
   ax.set_title('Cell Types')
   ax.set_xlabel('UMAP 1')
   ax.set_ylabel('UMAP 2')
   
   # Add legend
   unique_types = multiome.obs['cell_type'].unique()
   colors = plt.cm.Set3(np.linspace(0, 1, len(unique_types)))
   for i, cell_type in enumerate(unique_types):
       ax.scatter([], [], c=[colors[i]], label=cell_type, s=20)
   ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
   
   plt.tight_layout()
   plt.show()

Step 6: Multiome Integration
----------------------------

Integrate RNA and ATAC modalities:

.. code-block:: python

   # Co-embed RNA and ATAC data
   print("Performing multiome integration...")
   
   multiome = pymega.coembed_data(
       multiome,
       method="cca",                    # Canonical Correlation Analysis
       n_components=30,                 # Number of co-embedding components
       rna_use_rep="X_rna_pca",        # RNA representation to use
       atac_use_rep="X_atac_svd"       # ATAC representation to use
   )
   
   print("✓ Co-embedding completed")
   
   # Create integrated UMAP
   integrated_umap = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
   multiome.obsm['X_integrated_umap'] = integrated_umap.fit_transform(
       multiome.obsm['X_integrated']
   )
   
   # Compare individual vs integrated embeddings
   fig, axes = plt.subplots(1, 3, figsize=(18, 5))
   
   cell_type_codes = pd.Categorical(multiome.obs['cell_type']).codes
   
   # RNA UMAP
   axes[0].scatter(multiome.obsm['X_rna_umap'][:, 0], 
                   multiome.obsm['X_rna_umap'][:, 1],
                   c=cell_type_codes, cmap='Set3', s=1)
   axes[0].set_title('RNA UMAP')
   
   # ATAC UMAP
   axes[1].scatter(multiome.obsm['X_atac_umap'][:, 0], 
                   multiome.obsm['X_atac_umap'][:, 1],
                   c=cell_type_codes, cmap='Set3', s=1)
   axes[1].set_title('ATAC UMAP')
   
   # Integrated UMAP
   axes[2].scatter(multiome.obsm['X_integrated_umap'][:, 0], 
                   multiome.obsm['X_integrated_umap'][:, 1],
                   c=cell_type_codes, cmap='Set3', s=1)
   axes[2].set_title('Integrated UMAP')
   
   for ax in axes:
       ax.set_xlabel('UMAP 1')
       ax.set_ylabel('UMAP 2')
   
   plt.tight_layout()
   plt.show()

Step 7: Trajectory Analysis
---------------------------

Infer developmental trajectories:

.. code-block:: python

   # Perform trajectory inference
   print("Inferring developmental trajectories using R-compatible ArchR algorithm...")
   
   # Define trajectory path (cell type progression)
   trajectory_path = ["Neural_Stem_Cells", "Progenitors", "Neurons"]
   
   multiome = pymega.add_trajectory(
       multiome,
       trajectory=trajectory_path,      # R: trajectory = c("Neural_Stem_Cells", ...)
       group_by="cell_type",            # R: group.by = "cell_type"
       reduction="pca",                 # R: reduction = "pca"
       dims=list(range(30)),           # R: dims = 1:30
       pre_filter_quantile=0.9,        # R: pre.filter.quantile = 0.9
       post_filter_quantile=0.9,       # R: post.filter.quantile = 0.9
       dof=250,                        # R: dof = 250
       spar=1.0,                       # R: spar = 1
       name="Trajectory"               # R: name = "Trajectory"
   )
   
   # Get trajectory data (R GetTrajectory equivalent)
   trajectory_result = pymega.get_trajectory_data(
       multiome,
       trajectory_name="Trajectory",
       group_every=1,
       log2_norm=True,
       smooth_window=11
   )
   
   print(f"✓ Trajectory inferred: {' → '.join(trajectory_path)}")
   print(f"  Trajectory bins: {trajectory_result['group_matrix'].shape[1]}")
   
   # Pseudotime is stored in 'Trajectory' column (0-100 range)
   if 'Trajectory' in multiome.obs.columns:
       print(f"  Cells in trajectory: {multiome.obs['Trajectory'].notna().sum()}")
       print(f"  Trajectory range: {multiome.obs['Trajectory'].min():.1f} - {multiome.obs['Trajectory'].max():.1f}")
   
   # Visualize trajectories
   fig, axes = plt.subplots(1, 2, figsize=(12, 5))
   
   # Trajectory on integrated UMAP
   scatter = axes[0].scatter(multiome.obsm['X_integrated_umap'][:, 0], 
                            multiome.obsm['X_integrated_umap'][:, 1],
                            c=multiome.obs['pseudotime'], 
                            cmap='viridis', s=1)
   axes[0].set_title('Developmental Trajectory')
   axes[0].set_xlabel('UMAP 1')
   axes[0].set_ylabel('UMAP 2')
   plt.colorbar(scatter, ax=axes[0], label='Pseudotime')
   
   # Pseudotime distribution by cell type
   sns.boxplot(data=multiome.obs, x='cell_type', y='pseudotime', ax=axes[1])
   axes[1].set_title('Pseudotime by Cell Type')
   axes[1].tick_params(axis='x', rotation=45)
   
   plt.tight_layout()
   plt.show()

Step 8: Gene Regulatory Network Inference
------------------------------------------

Build gene regulatory networks:

.. code-block:: python

   # Prepare data for network inference
   print("Preparing for network inference...")
   
   # Get transcription factors (TFs)
   # In practice, you would load a curated TF list
   tf_list = [gene for gene in multiome.rna_var_names 
              if any(tf_keyword in gene.lower() 
                    for tf_keyword in ['sox', 'pou', 'tbr', 'neurog', 'ascl'])][:50]
   
   if len(tf_list) < 10:
       # Use highly variable genes as proxy TFs
       tf_list = multiome.var[multiome.var['highly_variable_rna']].index[:50].tolist()
   
   # Get target genes (highly variable genes)
   target_genes = multiome.var[multiome.var['highly_variable_rna']].index[:500].tolist()
   
   print(f"Using {len(tf_list)} TFs and {len(target_genes)} target genes")
   
   # Infer gene regulatory network
   print("Inferring gene regulatory network...")
   
   grn_network = pymega.infer_grn(
       multiome,
       tf_list=tf_list,
       target_genes=target_genes,
       method="correlation",            # Network inference method
       min_correlation=0.3,             # Minimum correlation threshold
       fdr_threshold=0.05,              # False discovery rate threshold
       use_pseudotime=True              # Use pseudotime for dynamic correlation
   )
   
   print(f"✓ Inferred GRN with {len(grn_network)} TF-gene interactions")
   
   # Display network summary
   print("\\nNetwork summary:")
   print(f"  Unique TFs: {grn_network['TF'].nunique()}")
   print(f"  Unique targets: {grn_network['target_gene'].nunique()}")
   print(f"  Average correlation: {grn_network['correlation'].mean():.3f}")
   print(f"  Strongest interaction: {grn_network.loc[grn_network['correlation'].idxmax(), 'TF']} → "
         f"{grn_network.loc[grn_network['correlation'].idxmax(), 'target_gene']} "
         f"(r={grn_network['correlation'].max():.3f})")

Link peaks to genes:

.. code-block:: python

   # Link ATAC peaks to genes
   print("Linking peaks to genes...")
   
   peak_gene_links = pymega.link_peaks_to_genes(
       multiome,
       distance_threshold=100000,       # 100kb distance threshold
       correlation_threshold=0.2,       # Minimum correlation
       use_pseudotime=True             # Use pseudotime for dynamic correlation
   )
   
   print(f"✓ Identified {len(peak_gene_links)} peak-gene links")
   
   # Display peak-gene linking summary
   print("\\nPeak-gene linking summary:")
   print(f"  Unique peaks: {peak_gene_links['peak'].nunique()}")
   print(f"  Unique genes: {peak_gene_links['gene'].nunique()}")
   print(f"  Average distance: {peak_gene_links['distance'].mean()/1000:.1f} kb")
   print(f"  Average correlation: {peak_gene_links['correlation'].mean():.3f}")

Step 9: Visualization and Results
---------------------------------

Create comprehensive visualizations:

.. code-block:: python

   # 1. Gene regulatory network visualization
   print("Creating network visualization...")
   
   # Select top interactions for visualization
   top_grn = grn_network.nlargest(100, 'correlation')
   
   fig_network = pymega.plot_grn_network(
       top_grn,
       layout="spring",
       node_size_col="degree",
       edge_width_col="correlation",
       node_color_col="node_type",
       title="Top 100 Gene Regulatory Interactions",
       figsize=(12, 10),
       interactive=False
   )
   plt.show()
   
   # 2. Trajectory gene expression heatmap
   print("Creating trajectory heatmap...")
   
   # Select genes that change along pseudotime
   trajectory_genes = grn_network.nlargest(50, 'correlation')['target_gene'].unique()[:20]
   
   fig_heatmap = pymega.plot_trajectory_heatmap(
       multiome.rna_data,
       trajectory_data,
       gene_list=trajectory_genes,
       n_bins=50,
       title="Gene Expression Along Trajectory"
   )
   plt.show()
   
   # 3. Cell type proportions
   print("Creating cell proportion plots...")
   
   fig_props = pymega.plot_cell_proportions(
       multiome.obs,
       group_col='cell_type',
       title="Cell Type Distribution"
   )
   plt.show()
   
   # 4. Peak-gene correlation plot
   print("Creating peak-gene correlation plot...")
   
   # Select top peak-gene links
   top_links = peak_gene_links.nlargest(20, 'correlation')
   
   fig, ax = plt.subplots(figsize=(10, 6))
   scatter = ax.scatter(top_links['distance']/1000, 
                       top_links['correlation'],
                       alpha=0.7, s=50)
   ax.set_xlabel('Distance (kb)')
   ax.set_ylabel('Peak-Gene Correlation')
   ax.set_title('Peak-Gene Links: Distance vs Correlation')
   
   # Add trend line
   z = np.polyfit(top_links['distance']/1000, top_links['correlation'], 1)
   p = np.poly1d(z)
   ax.plot(sorted(top_links['distance']/1000), 
           p(sorted(top_links['distance']/1000)), 
           "r--", alpha=0.8)
   
   plt.tight_layout()
   plt.show()

Generate summary statistics:

.. code-block:: python

   # Generate analysis summary
   print("\\n" + "="*60)
   print("ANALYSIS SUMMARY")
   print("="*60)
   
   print(f"Dataset: Mouse brain multiome data")
   print(f"Final cell count: {multiome.n_obs:,}")
   print(f"Final gene count: {multiome.n_rna_vars:,}")
   print(f"Final peak count: {multiome.n_atac_vars:,}")
   
   print(f"\\nCell type annotation:")
   cell_type_counts = multiome.obs['cell_type'].value_counts()
   for cell_type, count in cell_type_counts.items():
       pct = count / len(multiome.obs) * 100
       print(f"  {cell_type}: {count:,} cells ({pct:.1f}%)")
   
   print(f"\\nTrajectory analysis:")
   print(f"  Number of lineages: {len(trajectory_data['lineages'])}")
   print(f"  Pseudotime range: {multiome.obs['pseudotime'].min():.2f} - {multiome.obs['pseudotime'].max():.2f}")
   
   print(f"\\nGene regulatory network:")
   print(f"  TF-gene interactions: {len(grn_network):,}")
   print(f"  Unique TFs: {grn_network['TF'].nunique()}")
   print(f"  Unique target genes: {grn_network['target_gene'].nunique()}")
   print(f"  Average correlation: {grn_network['correlation'].mean():.3f}")
   
   print(f"\\nPeak-gene links:")
   print(f"  Total links: {len(peak_gene_links):,}")
   print(f"  Average distance: {peak_gene_links['distance'].mean()/1000:.1f} kb")
   print(f"  Average correlation: {peak_gene_links['correlation'].mean():.3f}")

Step 10: Save Results
---------------------

Save your analysis results:

.. code-block:: python

   # Create output directory
   output_dir = Path("results/mouse_brain_analysis/")
   output_dir.mkdir(parents=True, exist_ok=True)
   
   print(f"Saving results to {output_dir}")
   
   # Save the processed multiome object
   multiome.write_h5ad(output_dir / "processed_multiome.h5ad")
   
   # Save network data
   grn_network.to_csv(output_dir / "grn_network.csv", index=False)
   peak_gene_links.to_csv(output_dir / "peak_gene_links.csv", index=False)
   
   # Save trajectory data
   trajectory_df = pd.DataFrame({
       'cell_id': multiome.obs.index,
       'pseudotime': multiome.obs['pseudotime'],
       'cell_type': multiome.obs['cell_type']
   })
   trajectory_df.to_csv(output_dir / "trajectory_data.csv", index=False)
   
   # Save embeddings
   embeddings_df = pd.DataFrame({
       'cell_id': multiome.obs.index,
       'rna_umap1': multiome.obsm['X_rna_umap'][:, 0],
       'rna_umap2': multiome.obsm['X_rna_umap'][:, 1],
       'atac_umap1': multiome.obsm['X_atac_umap'][:, 0],
       'atac_umap2': multiome.obsm['X_atac_umap'][:, 1],
       'integrated_umap1': multiome.obsm['X_integrated_umap'][:, 0],
       'integrated_umap2': multiome.obsm['X_integrated_umap'][:, 1]
   })
   embeddings_df.to_csv(output_dir / "embeddings.csv", index=False)
   
   # Save analysis parameters
   analysis_params = {
       'qc_params': qc_params,
       'n_tf_genes': len(tf_list),
       'n_target_genes': len(target_genes),
       'grn_correlation_threshold': 0.3,
       'peak_gene_distance_threshold': 100000,
       'analysis_date': pd.Timestamp.now().isoformat()
   }
   
   import json
   with open(output_dir / "analysis_parameters.json", 'w') as f:
       json.dump(analysis_params, f, indent=2)
   
   print("✓ All results saved successfully!")
   
   # Print file list
   print("\\nSaved files:")
   for file_path in sorted(output_dir.glob("*")):
       file_size = file_path.stat().st_size / 1024 / 1024  # MB
       print(f"  {file_path.name} ({file_size:.1f} MB)")

Conclusion
----------

Congratulations! You've completed a comprehensive single-cell multiome analysis using PyMEGA. In this tutorial, you learned how to:

✅ Load and explore 10X multiome data  
✅ Perform quality control and normalization  
✅ Annotate cell types and find marker genes  
✅ Integrate RNA and ATAC modalities  
✅ Infer developmental trajectories  
✅ Build gene regulatory networks  
✅ Link peaks to genes  
✅ Create publication-ready visualizations  
✅ Save and organize results  

Next Steps
----------

Now that you've mastered the basic workflow, consider exploring:

1. **Advanced Analysis**: 
   - Custom network inference methods
   - Spatial analysis if you have spatial data

2. **Performance Optimization**:
   - :doc:`../user_guide/performance_optimization` - Speed up your analysis
   - Use parallel processing for large datasets
   - Memory optimization techniques

3. **Visualization**:
   - Interactive plots with Plotly
   - Custom visualization functions
   - Integration with other plotting libraries

4. **Integration with Other Tools**:
   - Export to Cytoscape for network analysis
   - Integration with R/Seurat workflows
   - Connection to pathway analysis tools

Happy analyzing! 🧬✨
