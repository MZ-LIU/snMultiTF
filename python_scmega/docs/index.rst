PyMEGA: Python Single-Cell Multiome Gene Regulatory Network Analysis
==================================================================

.. image:: _static/pymega_logo.png
   :width: 200px
   :align: center
   :alt: PyMEGA Logo

**PyMEGA** (Python Multiome Epigenome and Gene Analysis) is a comprehensive Python package for single-cell multiome data analysis and gene regulatory network inference. It provides a complete toolkit for analyzing paired scRNA-seq and scATAC-seq data to uncover gene regulatory mechanisms.

Key Features
------------

✨ **Comprehensive Analysis Pipeline**
   - Complete single-cell multiome data processing workflow
   - From raw data loading to network inference and visualization
   - Modular design for flexible analysis customization

🧬 **Multiome Data Integration**
   - Paired scRNA-seq and scATAC-seq analysis
   - Cell type annotation and trajectory inference
   - Peak-gene linking and regulatory element annotation

⚡ **High Performance Computing**
   - Numba JIT acceleration for computational bottlenecks
   - Parallel processing for multi-core systems
   - Memory-efficient algorithms for large datasets
   - Sparse matrix optimizations

🎨 **Rich Visualizations**
   - Static and interactive network visualizations
   - Trajectory and pseudotime analysis plots
   - Expression heatmaps and correlation analysis
   - Spatial and cell proportion visualizations

🔬 **Gene Regulatory Networks**
   - Advanced network inference algorithms
   - TF-gene correlation analysis
   - Peak-gene linking with genomic distance
   - Network validation and quality assessment

📊 **Benchmarking and Validation**
   - Comprehensive performance benchmarking
   - Accuracy validation against reference implementations
   - Memory and speed profiling tools

Quick Start
-----------

Installation
~~~~~~~~~~~~

.. code-block:: bash

   # Install from PyPI (when available)
   pip install pymega
   
   # Or install from source
   git clone https://github.com/pymega/pymega.git
   cd pymega
   pip install -e .

Basic Usage
~~~~~~~~~~~

.. code-block:: python

   import pymega
   
   # Load 10X multiome data
   multiome = pymega.load_10x_multiome(
       "path/to/10x/data",
       rna_h5="filtered_feature_bc_matrix.h5",
       atac_h5="atac_fragments.h5"
   )
   
   # Quality control and normalization
   multiome = pymega.quality_control_multiome(multiome)
   multiome = pymega.normalize_multiome(multiome)
   
   # Trajectory analysis
   multiome = pymega.add_trajectory(multiome)
   
   # Gene regulatory network inference
   grn_network = pymega.infer_grn(
       multiome,
       tf_list=tf_genes,
       target_genes=target_genes
   )
   
   # Visualization
   pymega.plot_grn_network(grn_network, interactive=True)

Table of Contents
-----------------

.. toctree::
   :maxdepth: 2
   :caption: User Guide
   
   user_guide/installation
   user_guide/quickstart
   user_guide/performance_optimization

.. toctree::
   :maxdepth: 2
   :caption: Tutorials
   
   tutorials/basic_workflow

.. toctree::
   :maxdepth: 2
   :caption: API Reference
   
   api_reference/index

.. toctree::
   :maxdepth: 2
   :caption: Developer Guide
   
   developer/contributing

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Citation
========

If you use PyMEGA in your research, please cite:

.. code-block:: bibtex

   @software{pymega2024,
     title={PyMEGA: Python Single-Cell Multiome Gene Regulatory Network Analysis},
     author={PyMEGA Development Team},
     year={2024},
     url={https://github.com/pymega/pymega},
     version={1.0.0}
   }

License
=======

PyMEGA is released under the MIT License. See the `LICENSE <https://github.com/pymega/pymega/blob/main/LICENSE>`_ file for details.

Support
=======

- **Documentation**: https://pymega.readthedocs.io
- **Issues**: https://github.com/pymega/pymega/issues
- **Discussions**: https://github.com/pymega/pymega/discussions
- **Email**: pymega-support@example.com

Acknowledgments
===============

PyMEGA builds upon excellent packages in the Python scientific ecosystem:

- `scanpy <https://scanpy.readthedocs.io/>`_ for single-cell analysis
- `anndata <https://anndata.readthedocs.io/>`_ for data structures
- `numba <https://numba.pydata.org/>`_ for high-performance computing
- `plotly <https://plotly.com/python/>`_ for interactive visualizations
- `networkx <https://networkx.org/>`_ for network analysis

We thank the developers and communities of these projects for their contributions to scientific computing.
