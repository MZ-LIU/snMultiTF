Visualization Module
===================

The visualization module provides comprehensive plotting functions for single-cell multiome data analysis.

.. currentmodule:: python_scmega.visualization

Network Visualization
---------------------

.. automodule:: python_scmega.visualization.network_plots
    :members:
    :undoc-members:
    :show-inheritance:

Trajectory Visualization
------------------------

.. automodule:: python_scmega.visualization.trajectory_plots
    :members:
    :undoc-members:
    :show-inheritance:

Heatmap Visualization
--------------------

.. automodule:: python_scmega.visualization.heatmaps
    :members:
    :undoc-members:
    :show-inheritance:


Cell Proportion Visualization
-----------------------------

.. automodule:: python_scmega.visualization.cell_proportion_plots
    :members:
    :undoc-members:
    :show-inheritance:




Examples
--------

Basic Network Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import python_scmega as pymega
    import pandas as pd
    
    # Create GRN network data
    grn_network = pd.DataFrame({
        'TF': ['TF1', 'TF1', 'TF2'],
        'target_gene': ['Gene1', 'Gene2', 'Gene1'],
        'confidence': [0.8, 0.6, 0.9]
    })
    
    # Plot static network
    fig = pymega.plot_grn_network(
        grn_network,
        layout="spring",
        interactive=False
    )
    
    # Plot interactive network
    interactive_fig = pymega.plot_grn_network(
        grn_network,
        interactive=True
    )

Trajectory Visualization
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import numpy as np
    
    # Create trajectory data
    n_cells = 100
    trajectory_data = pd.DataFrame({
        'pseudotime': np.linspace(0, 1, n_cells),
        'cell_type': np.random.choice(['Type1', 'Type2'], n_cells)
    })
    
    embedding = pd.DataFrame({
        'UMAP1': np.random.normal(0, 1, n_cells),
        'UMAP2': np.random.normal(0, 1, n_cells)
    })
    
    # Plot trajectory
    fig = pymega.plot_trajectory(
        trajectory_data,
        embedding,
        color_by="pseudotime"
    )

Heatmap Visualization
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    # Create expression data
    n_cells, n_genes = 100, 50
    expression_data = pd.DataFrame(
        np.random.normal(0, 1, (n_cells, n_genes)),
        columns=[f'Gene{i}' for i in range(n_genes)]
    )
    
    # Plot trajectory heatmap
    fig = pymega.plot_trajectory_heatmap(
        expression_data,
        trajectory_data,
        gene_list=expression_data.columns[:20],
        n_bins=30,
        cluster_genes=True
    )




