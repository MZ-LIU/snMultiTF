Installation Guide
==================

This guide covers the installation of PyMEGA and its dependencies.

System Requirements
-------------------

**Operating Systems:**
- Linux (Ubuntu 18.04+, CentOS 7+)
- macOS (10.14+)
- Windows (10+)

**Python Version:**
- Python 3.8 or higher
- Python 3.9-3.11 recommended

**Hardware Recommendations:**
- **RAM**: 16GB minimum, 32GB+ recommended for large datasets
- **CPU**: Multi-core processor (4+ cores recommended)
- **Storage**: 10GB+ free space for data and temporary files

Installation Methods
--------------------

Method 1: Install from PyPI (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once PyMEGA is published to PyPI, you can install it using pip:

.. code-block:: bash

   # Install latest stable version
   pip install pymega
   
   # Install with all optional dependencies
   pip install pymega[all]
   
   # Install specific optional dependencies
   pip install pymega[visualization,genomics,performance]

Method 2: Install from Source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the latest development version or to contribute to PyMEGA:

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/pymega/pymega.git
   cd pymega
   
   # Create a virtual environment (recommended)
   python -m venv pymega_env
   source pymega_env/bin/activate  # On Windows: pymega_env\Scripts\activate
   
   # Install in development mode
   pip install -e .
   
   # Or install with development dependencies
   pip install -e .[dev]

Method 3: Using Conda/Mamba
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Create a new conda environment
   conda create -n pymega python=3.9
   conda activate pymega
   
   # Install dependencies from conda-forge
   conda install -c conda-forge scanpy anndata h5py numba plotly networkx
   
   # Install PyMEGA
   pip install pymega

Dependencies
------------

Core Dependencies
~~~~~~~~~~~~~~~~~

These are automatically installed with PyMEGA:

.. code-block:: text

   numpy>=1.20.0
   pandas>=1.3.0
   scipy>=1.7.0
   scikit-learn>=1.0.0
   matplotlib>=3.5.0
   seaborn>=0.11.0

Single-Cell Analysis
~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   scanpy>=1.8.0
   anndata>=0.8.0
   h5py>=3.0.0

Optional Dependencies
~~~~~~~~~~~~~~~~~~~~~

**Performance (install with ``pip install pymega[performance]``):**

.. code-block:: text

   numba>=0.56.0
   joblib>=1.1.0
   psutil>=5.8.0

**Visualization (install with ``pip install pymega[visualization]``):**

.. code-block:: text

   plotly>=5.0.0
   networkx>=2.6.0
   ipywidgets>=7.6.0

**Genomics (install with ``pip install pymega[genomics]``):**

.. code-block:: text

   pyranges>=0.0.111
   pybigwig>=0.3.18
   pysam>=0.19.0

**Machine Learning (install with ``pip install pymega[ml]``):**

.. code-block:: text

   umap-learn>=0.5.0
   scanpy[leiden]>=1.8.0

**Development (install with ``pip install pymega[dev]``):**

.. code-block:: text

   pytest>=6.0.0
   pytest-cov>=3.0.0
   black>=22.0.0
   flake8>=4.0.0
   mypy>=0.910
   sphinx>=4.0.0
   sphinx-rtd-theme>=1.0.0

Verification
------------

Test Your Installation
~~~~~~~~~~~~~~~~~~~~~~~

After installation, verify that PyMEGA is working correctly:

.. code-block:: python

   import python_scmega as pymega
   print(f"PyMEGA version: {pymega.__version__}")
   
   # Test basic functionality
   import numpy as np
   test_data = np.random.randn(100, 50)
   
   # Test correlation calculation
   corr_matrix = np.corrcoef(test_data.T)
   print(f"Correlation matrix shape: {corr_matrix.shape}")
   
   print("✓ PyMEGA installed successfully!")

Run Example Analysis
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Run examples to test functionality
   # See examples/ folder for complete workflows
   # Example: examples/01_basic_usage.py
   
   import python_scmega as pymega
   print("PyMEGA is ready to use!")

Troubleshooting
---------------

Common Installation Issues
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Issue: "No module named 'numba'"**

Solution:

.. code-block:: bash

   pip install numba

**Issue: "Failed building wheel for numba"**

This usually occurs on systems without a C compiler. Solutions:

.. code-block:: bash

   # On Ubuntu/Debian
   sudo apt-get install build-essential
   
   # On CentOS/RHEL
   sudo yum groupinstall "Development Tools"
   
   # On macOS
   xcode-select --install
   
   # Then reinstall
   pip install numba

**Issue: "Microsoft Visual C++ 14.0 is required" (Windows)**

Solution:
1. Install Microsoft C++ Build Tools from https://visualstudio.microsoft.com/visual-cpp-build-tools/
2. Or install Visual Studio Community

**Issue: Memory errors during installation**

Solution:

.. code-block:: bash

   # Increase pip timeout and use no cache
   pip install --no-cache-dir --timeout 1000 pymega

**Issue: "Failed to import scanpy"**

Solution:

.. code-block:: bash

   # Install scanpy separately
   pip install scanpy
   
   # Or use conda
   conda install -c conda-forge scanpy

Performance Issues
~~~~~~~~~~~~~~~~~~

**Slow correlation calculations:**

1. Ensure numba is installed: ``pip install numba``
2. Check if JIT compilation is working: ``python -c "import numba; print(numba.__version__)"`
3. For very large datasets, use chunked processing

**Memory errors:**

1. Increase system memory if possible
2. Use chunked processing: ``pymega.chunked_data_processor()``
3. Optimize data types: ``pymega.optimize_data_types()``
4. Use sparse matrices for ATAC-seq data

**Slow parallel processing:**

1. Check number of available cores: ``pymega.setup_parallel_backend()``
2. Adjust chunk sizes for your system
3. Consider using ``backend="threading"`` vs ``backend="multiprocessing"``

Environment-Specific Issues
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Jupyter Notebook Integration:**

.. code-block:: bash

   # Install jupyter extensions
   pip install ipywidgets
   jupyter nbextension enable --py widgetsnbextension
   
   # For JupyterLab
   jupyter labextension install @jupyter-widgets/jupyterlab-manager

**Cluster/HPC Environments:**

.. code-block:: bash

   # Load required modules (example for SLURM)
   module load python/3.9
   module load gcc/9.3.0
   
   # Create virtual environment
   python -m venv --system-site-packages pymega_env
   source pymega_env/bin/activate
   
   # Install with minimal dependencies
   pip install pymega --no-deps
   pip install numpy pandas scipy matplotlib

Getting Help
------------

If you encounter issues not covered here:

1. **Check the documentation**: https://pymega.readthedocs.io
2. **Search existing issues**: https://github.com/pymega/pymega/issues
3. **Ask questions**: https://github.com/pymega/pymega/discussions
4. **Report bugs**: https://github.com/pymega/pymega/issues/new

When reporting issues, please include:

- Operating system and version
- Python version
- PyMEGA version
- Complete error message
- Minimal code example that reproduces the issue

Next Steps
----------

After successful installation, proceed to:

- :doc:`quickstart` - Quick introduction to PyMEGA
- :doc:`../tutorials/basic_workflow` - Complete analysis tutorial
