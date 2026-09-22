Data Processing Module
======================

The data processing module provides comprehensive tools for loading, quality control,
normalization, and format conversion of single-cell multiome data.

.. currentmodule:: python_scmega.data_processing

Core Classes
------------

.. autosummary::
   :toctree: generated/

   MultiomeData
   load_10x_multiome

Quality Control
---------------

.. autosummary::
   :toctree: generated/

   calculate_qc_metrics
   filter_cells_and_features
   quality_control_multiome
   plot_qc_metrics
   get_qc_summary

Normalization
-------------

.. autosummary::
   :toctree: generated/

   normalize_rna
   normalize_atac
   normalize_multiome
   run_umap
   find_neighbors_and_clusters
   get_normalization_summary

Format Conversion
-----------------

.. autosummary::
   :toctree: generated/

   convert_seurat_to_multiome
   convert_multiome_to_seurat
   convert_h5ad_to_multiome
   convert_multiome_to_h5ad
   validate_format_conversion

Overview
--------

The data processing module implements the complete preprocessing pipeline for 10X multiome data,
equivalent to the R scMEGA workflow:

**R Workflow:**

.. code-block:: r

   # Load data
   inputdata.10x <- Read10X_h5("data.h5")
   rna_counts <- inputdata.10x$`Gene Expression`
   atac_counts <- inputdata.10x$Peaks
   
   # Quality control
   obj.rna <- CreateSeuratObject(counts = rna_counts)
   obj.rna[["percent.mt"]] <- PercentageFeatureSet(obj.rna, pattern = "^MT-")
   obj.rna <- subset(obj.rna, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 20)
   
   # Normalization
   obj.rna <- obj.rna %>%
       NormalizeData(verbose=FALSE) %>%
       FindVariableFeatures(nfeatures=3000, verbose=F) %>%
       ScaleData() %>%
       RunPCA(verbose = FALSE)
   
   obj.atac <- obj.atac %>%
       RunTFIDF() %>%
       FindTopFeatures() %>%
       RunSVD()

**Python Equivalent:**

.. code-block:: python

   import python_scmega as pymega
   
   # Load data
   multiome = pymega.load_10x_multiome("data.h5")
   
   # Quality control
   multiome_qc = pymega.quality_control_multiome(
       multiome,
       rna_params={'mt_pattern': '^MT-'},
       filter_params={
           'min_features': 200,
           'max_features': 3000,
           'max_percent_mito': 20
       }
   )
   
   # Normalization
   multiome_norm = pymega.normalize_multiome(
       multiome_qc,
       rna_params={
           'target_sum': 1e4,
           'n_top_genes': 3000,
           'n_pcs': 50
       },
       atac_params={
           'n_top_features': 50000,
           'n_components': 50
       }
   )

Key Features
------------

**10X Multiome Optimization**
   - Specialized handling for 10X Genomics multiome data
   - Automatic cell barcode matching between RNA and ATAC
   - Chromosome filtering for ATAC peaks

**Quality Control**
   - Mitochondrial gene percentage calculation (PercentageFeatureSet equivalent)
   - Comprehensive QC metrics for both RNA and ATAC data
   - Flexible filtering based on multiple criteria
   - QC visualization and summary statistics

**Normalization**
   - RNA: Log-normalization, scaling, highly variable genes, PCA
   - ATAC: TF-IDF normalization, top feature selection, LSI/SVD
   - Integrated pipeline preserving multiome structure

**Format Conversion**
   - Bidirectional Seurat <-> MultiomeData conversion
   - H5AD format support for scanpy ecosystem
   - Validation tools for conversion accuracy

Examples
--------

**Basic Quality Control:**

.. code-block:: python

   from python_scmega.data_processing import calculate_qc_metrics, filter_cells_and_features
   
   # Calculate QC metrics
   rna_qc = calculate_qc_metrics(adata, data_type="RNA", mt_pattern="^MT-")
   
   # Filter based on QC metrics
   rna_filtered = filter_cells_and_features(
       rna_qc,
       data_type="RNA",
       min_features=200,
       max_features=3000,
       max_percent_mito=20
   )

**RNA Normalization Pipeline:**

.. code-block:: python

   from python_scmega.data_processing import normalize_rna
   
   # Complete RNA normalization
   rna_norm = normalize_rna(
       adata,
       target_sum=1e4,
       n_top_genes=3000,
       n_pcs=50
   )

**ATAC Normalization Pipeline:**

.. code-block:: python

   from python_scmega.data_processing import normalize_atac
   
   # Complete ATAC normalization
   atac_norm = normalize_atac(
       adata,
       n_top_features=50000,
       n_components=50
   )

**Format Conversion:**

.. code-block:: python

   from python_scmega.data_processing import convert_multiome_to_h5ad
   
   # Convert to H5AD format
   saved_files = convert_multiome_to_h5ad(
       multiome,
       output_dir="output/",
       save_separate=True,
       save_combined=True
   )
