Performance Tips and Troubleshooting Guide
==========================================

This guide provides tips for working with large datasets and troubleshooting common issues.

Hardware Recommendations
------------------------

**For Small Datasets (< 10k cells):**
- RAM: 8GB minimum
- CPU: 4+ cores
- Storage: SSD recommended

**For Medium Datasets (10k-100k cells):**
- RAM: 32GB minimum
- CPU: 8+ cores
- Storage: Fast SSD (NVMe)

**For Large Datasets (> 100k cells):**
- RAM: 64GB+ recommended
- CPU: 16+ cores
- Storage: High-speed NVMe SSD
- Consider cluster computing

General Performance Tips
-------------------------

1. **Data Preprocessing**

   - Filter low-quality cells early to reduce dataset size
   - Use appropriate quality control thresholds
   - Remove unnecessary features before analysis

   .. code-block:: python

      # Filter cells and features early
      mdata = pymega.quality_control_separate(
          mdata,
          min_rna_genes=200,
          min_atac_peaks=500
      )

2. **Memory Management**

   - Work with sparse matrices when possible
   - Process data in chunks for very large datasets
   - Clear intermediate results when no longer needed

   .. code-block:: python

      import gc
      
      # Clear memory after intensive operations
      del intermediate_result
      gc.collect()

3. **Efficient Data Structures**

   - Use MuData for multiome data instead of separate objects
   - Leverage AnnData's efficient sparse matrix storage
   - Use appropriate data types (e.g., float32 instead of float64 when precision allows)

4. **Analysis Optimization**

   - Select highly variable features before network inference
   - Use subset of cells for parameter tuning
   - Cache intermediate results

   .. code-block:: python

      # Work with highly variable features
      selected_genes = pymega.select_genes(
          mdata,
          var_cutoff=0.9,
          n_genes=2000
      )

Troubleshooting Common Issues
------------------------------

Memory Errors
~~~~~~~~~~~~~

**Issue:** Out of memory errors during analysis

**Solutions:**

1. Reduce dataset size through quality control
2. Select subset of features
3. Process in smaller batches
4. Use a machine with more RAM or cloud computing

.. code-block:: python

   # Example: Process in batches
   batch_size = 5000
   n_cells = mdata['rna'].n_obs
   
   for i in range(0, n_cells, batch_size):
       batch_data = mdata[i:i+batch_size]
       # Process batch
       ...

Slow Performance
~~~~~~~~~~~~~~~~

**Issue:** Analysis takes too long

**Solutions:**

1. Use highly variable features only
2. Reduce number of TFs and genes for network inference
3. Use appropriate correlation thresholds to filter edges
4. Consider cloud computing for very large datasets

.. code-block:: python

   # Reduce computational burden
   tf_gene_cor = pymega.get_tf_gene_correlation(
       mdata,
       tf_list=selected_tfs[:50],  # Limit TFs
       gene_list=hvg[:500],         # Limit genes
       min_correlation=0.3          # Filter weak correlations
   )

Network Inference Issues
~~~~~~~~~~~~~~~~~~~~~~~~

**Issue:** Network inference produces too many/too few edges

**Solutions:**

1. Adjust correlation threshold
2. Filter by FDR (false discovery rate)
3. Check trajectory quality
4. Verify peak-gene links

.. code-block:: python

   # Filter network by multiple criteria
   grn_filtered = grn_network[
       (grn_network['correlation'].abs() > 0.4) &
       (grn_network['fdr'] < 0.05)
   ]

Visualization Problems
~~~~~~~~~~~~~~~~~~~~~~

**Issue:** Plots are slow or crash

**Solutions:**

1. Reduce number of cells/features shown
2. Use static plots instead of interactive
3. Downsample cells for visualization
4. Export data and visualize in smaller chunks

.. code-block:: python

   # Downsample for visualization
   import scanpy as sc
   
   # Downsample to 10k cells
   sc.pp.subsample(mdata['rna'], n_obs=10000)
   
   # Then create plots
   pymega.plot_trajectory(mdata, ...)

Best Practices for Large Datasets
----------------------------------

1. **Progressive Analysis**

   - Start with small subset for parameter tuning
   - Validate pipeline on subset
   - Apply to full dataset once optimized

2. **Checkpoint Results**

   - Save intermediate results frequently
   - Use HDF5 format for efficient storage
   - Can resume from checkpoints if interrupted

   .. code-block:: python

      # Save checkpoints
      mdata.write('checkpoint_normalized.h5mu')
      
      # Load checkpoint
      import muon as mu
      mdata = mu.read('checkpoint_normalized.h5mu')

3. **Resource Monitoring**

   - Monitor RAM and CPU usage
   - Estimate resource needs before large runs
   - Use job schedulers for cluster computing

4. **Parallel Processing**

   - Most PyMEGA functions use optimized NumPy/SciPy
   - These libraries already use multi-threading
   - Set number of threads via environment variables if needed

   .. code-block:: bash

      # Set number of threads (bash)
      export OMP_NUM_THREADS=8
      export MKL_NUM_THREADS=8

Cloud Computing Options
-----------------------

For very large datasets, consider:

1. **Google Colab**: Free tier with GPU (limited resources)
2. **AWS/Azure/GCP**: Scalable resources
3. **HPC Clusters**: If available at your institution

Example Colab setup:

.. code-block:: python

   # Install PyMEGA in Colab
   !pip install python-scmega
   
   # Mount Google Drive for data access
   from google.colab import drive
   drive.mount('/content/drive')

Getting Help
------------

If you encounter issues:

1. Check this troubleshooting guide
2. Review example workflows in ``examples/``
3. Open an issue on GitHub with:
   - Dataset size (cells, features)
   - System resources (RAM, CPU)
   - Complete error message
   - Minimal reproducible example

Additional Resources
--------------------

- **scanpy performance tips**: https://scanpy.readthedocs.io/en/stable/tutorials.html
- **AnnData documentation**: https://anndata.readthedocs.io/
- **MuData documentation**: https://muon.readthedocs.io/
