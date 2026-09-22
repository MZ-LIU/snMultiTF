"""
GRN Network Visualization - R-Compatible Implementation

This module provides Python implementations of scMEGA-style GRN plots.
TF nodes are drawn with a fixed light-blue color.

Completely based on the R version of GRNPlot implementation (grn.R:253-439)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import networkx as nx
from typing import Optional, Dict, Any, List, Union, Tuple
import warnings
from scipy import stats
from sklearn.preprocessing import MinMaxScaler

TF_FIXED_COLOR = "#4C9BD6"  # deep blue
TARGET_FIXED_COLOR = "#D9D9D9"  # light gray
EDGE_FIXED_COLOR = "gray"  # ordinary edges


def plot_grn_network(
    grn_network: pd.DataFrame,
    tfs_timepoint: Dict[str, float],
    tfs_use: Optional[List[str]] = None,
    genes_use: Optional[List[str]] = None,
    genes_cluster: Optional[pd.DataFrame] = None,
    genes_highlight: Optional[List[str]] = None,
    colors_highlight: str = "#984ea3",
    show_tf_labels: bool = True,
    seed: int = 42,
    plot_importance: bool = True,
    min_importance: float = 2.0,
    remove_isolated: bool = False,
    layout_algorithm: str = "fruchterman_reingold",
    figsize: Tuple[int, int] = (14, 12),
    title: Optional[str] = None,
    save_path: Optional[str] = None
) -> plt.Figure:
    """
    Draw a gene regulatory network diagram.
    
    This is a Python implementation of the R version of the GRNPlot function, retaining some interface compatibility.
    
    **Core rules** (compare with R language implementation):
    
    1. **TF node color definition**:
       - TF nodes use fixed light blue uniformly
       - No more sorting and continuous colormap mapping based on pseudo-time or time_point
       - Show label of TF node (gene name, orange #ff7f00)
    
    2. **TF node size definition**:
       - Treat normalized PageRank and betweenness centrality as coordinates on a two-dimensional plane
       - Calculate the Euclidean distance from each TF to the "least important TF" (minimum point of the two indicators)
       - This distance is used as the comprehensive importance score of TF to control the node size
       - Formula: importance = sqrt((pagerank_scaled - min_page)² + (betweenness_scaled - min_bet)²)
       - Corresponding R code: grn.R:289-313
    
    3. **TF time sequence definition**:
       - Compare R language rules: R/select_tf_gene.R:65-66
       - TF activity matrix sorted based on TrajectoryHeatmap
       - Time points are evenly distributed from 1 to 100
    
    4. **Other features**:
       - Gene nodes are colored according to clustering or unified to gray
       - The gene node size is the TF minimum importance value
       - Node transparency: TF=1.0, Gene=0.5
    
    Corresponding R code: grn.R:253-439
    
    Args:
        grn_network: GRN network data frame, including columns: tf, gene, correlation, p_value, fdr, etc.
        tfs_timepoint: reserved for compatibility with older calls; does not currently participate in TF sorting or color mapping.
        tfs_use: list of TFs to be displayed (optional, all TFs are displayed by default)
        genes_use: list of genes to display (optional, default shows all genes)
        genes_cluster: Gene clustering information, DataFrame contains columns: gene, cluster
        genes_highlight: list of genes to highlight
        colors_highlight: Color of highlighted genes
        show_tf_labels: whether to display TF labels
        seed: random seed (for layout)
        plot_importance: whether to draw a TF importance scatter plot
        min_importance: Minimum importance threshold (used to filter displayed TF tags)
        remove_isolated: whether to remove isolated nodes
        layout_algorithm: layout algorithm
        figsize: graphic size
        title: graphic title
        save_path: save path
        
    Returns:
        matplotlib.Figure: drawn graphics object
        
    Example:
        ```python
        # Prepare TF pseudo-time points (obtained from SelectTFs results or other methods)
        tfs_timepoint = {
            "GATA1": 20.5,
            "PU1": 45.3,
            "CEBPA": 78.9
        }
        
        # Draw the network
        fig = plot_grn_network_r_style(
            grn_network=df_grn,
            tfs_timepoint=tfs_timepoint,
            show_tf_labels=True,
            plot_importance=True
        )
        ```
    """
    
    print(' Draw GRN network diagram (R-compatible version)...')
    
    df_grn = grn_network.copy()
    
    if tfs_use is not None:
        df_grn = df_grn[df_grn['tf'].isin(tfs_use)]
        print(f"  Filter TF: {len(tfs_use)}")
    
    if genes_use is not None:
        df_grn = df_grn[df_grn['gene'].isin(genes_use)]
        print(f"  Screening genes: {len(genes_use)}")
    
    if df_grn.empty:
        raise ValueError('There are no edges after filtering! Please check the tfs_use and genes_use parameters.')
    
    tf_list = df_grn['tf'].unique().tolist()
    gene_list = df_grn['gene'].unique().tolist()
    # Remove genes that are also TF from the gene list (to avoid duplication)
    gene_list = [g for g in gene_list if g not in tf_list]
    
    print(f"  The network contains {len(tf_list)} TFs, {len(gene_list)} genes")
    print(f"  Total {len(df_grn)} regulatory relationships")
    
    G = nx.DiGraph()  # directed graph
    
    for _, row in df_grn.iterrows():
        tf = row['tf']
        gene = row['gene']
        
        if 'weights' in df_grn.columns:
            w = row['weights']
        elif 'weight' in df_grn.columns:
            w = row['weight']
        else:
            w = row['correlation']
        # Write 'weight' and 'weights' together to be compatible with networkx and R semantics
        G.add_edge(tf, gene, weight=w, weights=w)
    
    print(f"  The graph contains {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")
    
    if remove_isolated:
        isolated = list(nx.isolates(G))
        G.remove_nodes_from(isolated)
        print(f"  Remove {len(isolated)} orphan nodes")
    
    # === Step 4: Calculate node importance ===
    # 
    # 1. Treat normalized PageRank and betweenness centrality as coordinates on a two-dimensional plane
    # 2. Calculate the Euclidean distance from each TF to the "least important TF" (minimum point of the two indicators)
    # 3. This distance is used as the comprehensive importance score of TF to control the node size.
    
    print('  Calculate node centrality index...')
    
    pagerank = nx.pagerank(G, weight='weight')
    pagerank_series = pd.Series(pagerank)
    
    betweenness = nx.betweenness_centrality(G, weight='weight', normalized=True)
    betweenness_series = pd.Series(betweenness)
    
    # Only calculate the importance of TF
    tf_nodes_in_graph = [node for node in G.nodes() if node in tf_list]
    
    df_measure = pd.DataFrame({
        'tf': tf_nodes_in_graph,
        'pagerank': [pagerank[tf] for tf in tf_nodes_in_graph],
        'betweenness': [betweenness[tf] for tf in tf_nodes_in_graph]
    })
    
    df_measure['pagerank'] = (df_measure['pagerank'] - df_measure['pagerank'].mean()) / df_measure['pagerank'].std()
    df_measure['betweenness'] = (df_measure['betweenness'] - df_measure['betweenness'].mean()) / df_measure['betweenness'].std()
    
    # Step 4: Calculate the Euclidean distance to the "least important TF" (minimum value point)
    # This distance is the overall importance score
    min_page = df_measure['pagerank'].min()  # Find the minimum PageRank (x-coordinate of the least important TF)
    min_bet = df_measure['betweenness'].min()  # Find the minimum Betweenness (y coordinate of the least important TF)
    
    df_measure['importance'] = np.sqrt(
        (df_measure['pagerank'] - min_page) ** 2 +  # The square of the distance in the x direction
        (df_measure['betweenness'] - min_bet) ** 2  # The square of the distance in the y direction
    )
    
    print(f"  Calculated importance scores for {len(df_measure)} TFs")
    
    # === Step 5: Draw TF importance graph (optional) ===
    if plot_importance:
        fig_importance, ax_importance = plt.subplots(figsize=(10, 6))
        
        df_measure_sorted = df_measure.sort_values('importance', ascending=False)
        
        ax_importance.scatter(range(len(df_measure_sorted)), 
                            df_measure_sorted['importance'].values,
                            alpha=0.6, s=50, color='steelblue')
        
        ax_importance.set_xlabel('TFs (sorted by importance)', fontsize=12)
        ax_importance.set_ylabel('Importance Score', fontsize=12)
        ax_importance.set_title('TF Importance (PageRank + Betweenness)', fontsize=14, fontweight='bold')
        ax_importance.grid(True, alpha=0.3)
        
        for i in range(min(10, len(df_measure_sorted))):
            tf_name = df_measure_sorted.iloc[i]['tf']
            importance = df_measure_sorted.iloc[i]['importance']
            ax_importance.annotate(tf_name, (i, importance), 
                                  xytext=(0, 5), textcoords='offset points',
                                  fontsize=8, ha='center', rotation=45)
        
        plt.tight_layout()
        
        if save_path:
            importance_path = save_path.replace('.png', '_importance.png')
            plt.savefig(importance_path, dpi=300, bbox_inches='tight')
            print(f"  Save TF importance map: {importance_path}")
        
        plt.show()
    
    # Filter important TF for label display
    df_measure_important = df_measure[df_measure['importance'] > min_importance]
    print(f"  Importance score of {len(df_measure_important)} TF > {min_importance}")
    
    
    # TF node size is determined by importance
    tf_size = {}
    for _, row in df_measure.iterrows():
        tf_size[row['tf']] = row['importance']
    
    min_tf_importance = df_measure['importance'].min()
    gene_size = {gene: min_tf_importance for gene in gene_list}
    
    node_size_dict = {**tf_size, **gene_size}
    
    # 
    # 2. No longer sorting and continuous colormap mapping based on pseudo-time or time_point
    
    print('  Assign node color...')
    
    # 7.1 TF color: fixed to light blue; tfs_timepoint is only retained for compatibility with old calls.
    cols_tf = {tf: TF_FIXED_COLOR for tf in tf_list}
    print(f"  Assign fixed light blue color to {len(tf_list)} TFs")
        
    if genes_cluster is not None and not genes_cluster.empty:
        unique_clusters = genes_cluster['cluster'].unique()
        n_clusters = len(unique_clusters)
        
        if n_clusters <= 10:
            cluster_colors = plt.cm.tab10(np.linspace(0, 1, n_clusters))
        else:
            cluster_colors = plt.cm.tab20(np.linspace(0, 1, n_clusters))
        
        cluster_color_map = {
            cluster: f"#{int(color[0]*255):02x}{int(color[1]*255):02x}{int(color[2]*255):02x}"
            for cluster, color in zip(unique_clusters, cluster_colors)
        }
        
        cols_gene = {}
        for _, row in genes_cluster.iterrows():
            if row['gene'] in gene_list:
                cols_gene[row['gene']] = cluster_color_map[row['cluster']]
        
        for gene in gene_list:
            if gene not in cols_gene:
                cols_gene[gene] = TARGET_FIXED_COLOR
        
        print(f"  Assign colors to genes based on {n_clusters} clusters")
    else:
        cols_gene = {gene: TARGET_FIXED_COLOR for gene in gene_list}
        print(f"  Color {len(gene_list)} genes using default gray")
    
    node_color_dict = {**cols_tf, **cols_gene}
    
    node_alpha_dict = {}
    for node in G.nodes():
        if node in tf_list:
            node_alpha_dict[node] = 1.0  # TF is completely opaque
        else:
            node_alpha_dict[node] = 1.0
    
    
    print(f"  Calculate network layout ({layout_algorithm})...")
    np.random.seed(seed)
    
    if layout_algorithm == "fruchterman_reingold" or layout_algorithm == "spring":
        pos = nx.spring_layout(
            G, 
            k=1.0,  # optimal distance between nodes
            iterations=1000,  # Corresponding to R's niter=1000
            weight='weight',
            seed=seed
        )
    elif layout_algorithm == "kamada_kawai":
        pos = nx.kamada_kawai_layout(G, weight='weight')
    elif layout_algorithm == "circular":
        pos = nx.circular_layout(G)
    else:
        pos = nx.spring_layout(G, seed=seed)
    
    print(f"  Layout calculation completed")
    
    print('  Draw network diagram...')
    
    fig, ax = plt.subplots(figsize=figsize)
    
    edges = list(G.edges())
    edge_colors = [EDGE_FIXED_COLOR] * len(edges)
    edge_alphas = [0.25] * len(edges)
    
    for i, (u, v) in enumerate(edges):
        x = [pos[u][0], pos[v][0]]
        y = [pos[u][1], pos[v][1]]
        ax.plot(x, y, color=edge_colors[i], alpha=edge_alphas[i], linewidth=0.5, zorder=1)
    
    
    node_x = []
    node_y = []
    node_sizes = []
    node_colors = []
    node_alphas = []
    
    for node in G.nodes():
        node_x.append(pos[node][0])
        node_y.append(pos[node][1])
        
        size = node_size_dict.get(node, 1.0)
        node_sizes.append(size)
        
        node_colors.append(node_color_dict.get(node, TARGET_FIXED_COLOR))
        
        node_alphas.append(node_alpha_dict.get(node, 0.5))
    
    if node_sizes:
        node_sizes_array = np.array(node_sizes).reshape(-1, 1)
        
        if np.isnan(node_sizes_array).all():
            node_sizes_scaled = [300] * len(node_x)
        elif np.nanmin(node_sizes_array) == np.nanmax(node_sizes_array):
            node_sizes_scaled = [500] * len(node_x)
        else:
            scaler = MinMaxScaler(feature_range=(100, 1000))
            node_sizes_scaled = scaler.fit_transform(node_sizes_array).flatten()
    else:
        node_sizes_scaled = [300] * len(node_x)
    
    # Draw nodes with different transparency respectively (because scatter does not support point-by-point transparency)
    gene_indices = [i for i, node in enumerate(G.nodes()) if node in gene_list]
    if gene_indices:
        ax.scatter(
            [node_x[i] for i in gene_indices],
            [node_y[i] for i in gene_indices],
            s=[node_sizes_scaled[i] for i in gene_indices],
            c=[node_colors[i] for i in gene_indices],
            alpha=0.5,
            edgecolors='none',
            zorder=2
        )
    
    tf_indices = [i for i, node in enumerate(G.nodes()) if node in tf_list]
    if tf_indices:
        ax.scatter(
            [node_x[i] for i in tf_indices],
            [node_y[i] for i in tf_indices],
            s=[node_sizes_scaled[i] for i in tf_indices],
            c=[node_colors[i] for i in tf_indices],
            alpha=1.0,
            edgecolors='black',
            linewidths=0.5,
            zorder=3
        )
    
    # **TF Node Color Definition Rules Article 3**: Display the label of the TF node (gene name)
    if show_tf_labels:
        print(f"  Add TF tag...")
        
        # Show tags for all TFs (or only important TFs)
        tfs_to_label = tf_list  # Show all TF
        # If you want to display only important TF, uncomment the following:
        # tfs_to_label = df_measure_important['tf'].tolist()
        
        for node in G.nodes():
            if node in tfs_to_label:
                ax.annotate(
                    node,  # Gene name
                    xy=(pos[node][0], pos[node][1]),
                    xytext=(5, 5),
                    textcoords='offset points',
                    fontsize=8,
                    color='#ff7f00',  # Orange (consistent with R version)
                    fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                             edgecolor='none', alpha=0.7),
                    zorder=4
                )
    
    if genes_highlight:
        print(f"  Highlight {len(genes_highlight)} genes...")
        
        for gene in genes_highlight:
            if gene in G.nodes():
                ax.annotate(
                    gene,
                    xy=(pos[gene][0], pos[gene][1]),
                    xytext=(5, 5),
                    textcoords='offset points',
                    fontsize=8,
                    color=colors_highlight,
                    fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                             edgecolor='none', alpha=0.7),
                    zorder=4
                )
    
    ax.set_aspect('equal')
    ax.axis('off')
    
    if title:
        ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    else:
        ax.set_title('Gene Regulatory Network', fontsize=16, fontweight='bold', pad=20)
    
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=TF_FIXED_COLOR, 
               markersize=10, label='Transcription Factor', markeredgecolor='black'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=TARGET_FIXED_COLOR,
               markersize=8, label='Target Gene', alpha=0.5)
    ]
    
    ax.legend(handles=legend_elements, loc='upper right', frameon=True, 
             facecolor='white', edgecolor='black', framealpha=0.9)
    
    # Use try-except to catch tight_layout warnings
    try:
        plt.tight_layout()
    except Exception:
        pass
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f" Save network diagram: {save_path}")
    
    print(' GRN network plot completed!')
    
    return fig


def create_tf_timepoint_dict(
    tfs_df: pd.DataFrame,
    tf_col: str = 'tfs',
    timepoint_col: str = 'time_point'
) -> Dict[str, float]:
    """
    Create a dictionary of TF pseudo-time points from SelectTFs results.
    
    Args:
        tfs_df: TF data frame returned by SelectTFs
        tf_col: TF name column
        timepoint_col: pseudo time point column
        
    Returns:
        Dict[str, float]: mapping from TF to pseudo time points
        
    Example:
        ```python
        #Create from SelectTFs results
        tfs_timepoint = create_tf_timepoint_dict(sel_tfs['tfs'])
        
        # then used for drawing
        fig = plot_grn_network_r_style(grn, tfs_timepoint)
        ```
    """
    if tf_col not in tfs_df.columns or timepoint_col not in tfs_df.columns:
        raise ValueError(f"The data frame must contain columns: {tf_col} and {timepoint_col}")
    
    return dict(zip(tfs_df[tf_col], tfs_df[timepoint_col]))



def plot_grn_network_circular(
    grn_network: pd.DataFrame,
    tfs_timepoint: Dict[str, float],
    tf_node_size: float = 500,
    gene_node_size: float = 300,
    inner_radius: float = 0.3,  # TF inner radius
    outer_radius: float = 0.6,  # Gene outer radius
    tf_radius_spread: float = 0.05,  # The dispersion range of TF nodes in the radial direction
    gene_radius_spread: float = 0.1,  # Radial dispersion range of Gene nodes
    show_tf_labels: bool = True,
    figsize: Tuple[int, int] = (15, 15),
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    edge_alpha: float = 0.15,
    edge_width: float = 0.5,
    seed: int = 42
) -> plt.Figure:
    """
    Draw a circular structure gene regulatory network diagram (inner and outer layer layout).
    
    Inner layer: TF nodes (distributed on the inner circle, fixed light blue)
    Outer layer: Gene nodes (distributed on the outer circle, unified gray circle)
    
    Args:
        grn_network: GRN network data, including 'tf', 'gene', 'correlation' columns
        tfs_timepoint: reserved for compatibility with older calls; does not currently participate in TF sorting or color mapping.
        tf_node_size: TF node size (unified)
        gene_node_size: Gene node size (unified)
        inner_radius: the radius of the inner circle of TF
        outer_radius: Radius of Gene's outer circle
        tf_radius_spread: Random dispersion range of TF nodes in the radial direction
        gene_radius_spread: Random dispersion range of Gene nodes in the radial direction
        show_tf_labels: whether to display TF labels
        figsize: graphic size
        title: graphic title
        save_path: save path
        edge_alpha: edge transparency
        edge_width: width of edge
        seed: random seed
        
    Returns:
        matplotlib.Figure: drawn graphics object
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import networkx as nx
    from matplotlib.patches import Circle
    
    print('Draw a circular inner and outer layout GRN network diagram...')
    
    np.random.seed(seed)
    
    df_grn = grn_network.copy()
    tf_list = df_grn['tf'].unique().tolist()
    gene_list = df_grn['gene'].unique().tolist()
    gene_list = [g for g in gene_list if g not in tf_list]
    
    print(f"  The network contains {len(tf_list)} TFs, {len(gene_list)} genes")
    print(f"  Total {len(df_grn)} regulatory relationships")
    
    G = nx.DiGraph()
    for _, row in df_grn.iterrows():
        tf = row['tf']
        gene = row['gene']
        if tf in tf_list and gene in gene_list:
            weight = row.get('correlation', row.get('weights', row.get('weight', 1.0)))
            G.add_edge(tf, gene, weight=weight)
    print(f"  Added {G.number_of_edges()} TF->Gene regulatory edges")

    # tfs_timepoint is ignored so TF layout and colors are not pseudotime-driven.
    tf_list_sorted = list(tf_list)
    
    pos = {}
    center_x, center_y = 0.0, 0.0  # Circle center coordinates
    n_tfs = len(tf_list_sorted)
    for i, tf in enumerate(tf_list_sorted):
        angle = 2 * np.pi * i / n_tfs if n_tfs > 0 else 0
        angle_jitter = np.random.uniform(-np.pi / (2 * max(n_tfs, 1)), 
                                         np.pi / (2 * max(n_tfs, 1)))
        angle += angle_jitter 
        radius = inner_radius + np.random.uniform(-tf_radius_spread, tf_radius_spread)
        x = center_x + radius * np.cos(angle)
        y = center_y + radius * np.sin(angle)
        pos[tf] = (x, y)
    
    n_genes = len(gene_list) 
    for i, gene in enumerate(gene_list):
        angle = 2 * np.pi * i / n_genes if n_genes > 0 else 0   
        angle_jitter = np.random.uniform(-np.pi / (2 * max(n_genes, 1)), 
                                         np.pi / (2 * max(n_genes, 1)))
        angle += angle_jitter   
        radius = outer_radius + np.random.uniform(-gene_radius_spread, gene_radius_spread)  
        x = center_x + radius * np.cos(angle)
        y = center_y + radius * np.sin(angle)   
        pos[gene] = (x, y)
    
    print(f"  The circular layout is completed: inner radius={inner_radius}, outer radius={outer_radius}")

    # pos = {}

    # n_tfs = len(tf_list_sorted)
    # for i, tf in enumerate(tf_list_sorted):
    #     angle = 2 * np.pi * i / n_tfs if n_tfs > 0 else 0
    #     angle_jitter = np.random.uniform(-np.pi / (2 * max(n_tfs, 1)), 
    #                                     np.pi / (2 * max(n_tfs, 1)))
    #     angle += angle_jitter 
    #     radius = inner_radius * np.sqrt(np.random.uniform(0, 1))
    #     x = center_x + radius * np.cos(angle)
    #     y = center_y + radius * np.sin(angle)
    #     pos[tf] = (x, y)

    # n_genes = len(gene_list)
    # for i, gene in enumerate(gene_list):
    #     angle = 2 * np.pi * i / n_genes if n_genes > 0 else 0
    #     angle_jitter = np.random.uniform(-np.pi / (2 * max(n_genes, 1)), 
    #                                     np.pi / (2 * max(n_genes, 1)))
    #     angle += angle_jitter
    #     r_squared = np.random.uniform(inner_radius**2, outer_radius**2)
    #     radius = np.sqrt(r_squared)
    #     x = center_x + radius * np.cos(angle)
    #     y = center_y + radius * np.sin(angle)
    #     pos[gene] = (x, y)


    tf_colors = {tf: TF_FIXED_COLOR for tf in tf_list_sorted}
    
    print('  Draw network diagram...')
    fig, ax = plt.subplots(figsize=figsize)
    
    # tf_circle = Circle(
    #     (center_x, center_y),
    #     inner_radius,
    #     fill=False,
    #     edgecolor='lightgreen',
    #     linestyle='--',
    #     linewidth=2,
    #     alpha=0.5,
    #     zorder=0
    # )
    # ax.add_patch(tf_circle)
    
    # gene_circle = Circle(
    #     (center_x, center_y),
    #     outer_radius,
    #     fill=False,
    #     edgecolor='lightgreen',
    #     linestyle='--',
    #     linewidth=2,
    #     alpha=0.5,
    #     zorder=0
    # )
    # ax.add_patch(gene_circle)
    tf_circle = Circle(
        (center_x, center_y),
        inner_radius,
        fill=False,
        edgecolor='lightblue',
        linestyle='--',
        linewidth=2,
        alpha=0.5,
        zorder=0
    )
    ax.add_patch(tf_circle)

    gene_circle = Circle(
        (center_x, center_y),
        outer_radius,
        fill=False,
        edgecolor='lightgreen',
        linestyle='--',
        linewidth=2,
        alpha=0.5,
        zorder=0
    )
    ax.add_patch(gene_circle)

    for u, v in G.edges():
        if u in pos and v in pos:
            x = [pos[u][0], pos[v][0]]
            y = [pos[u][1], pos[v][1]]
            ax.plot(x, y, color=EDGE_FIXED_COLOR, alpha=edge_alpha, linewidth=edge_width, zorder=1)
    
    print(f"  Plotted {G.number_of_edges()} edges")
    
    gene_x = [pos[g][0] for g in gene_list if g in pos]
    gene_y = [pos[g][1] for g in gene_list if g in pos]
    
    if gene_x:
        ax.scatter(gene_x, gene_y, s=gene_node_size, c=TARGET_FIXED_COLOR,
                  alpha=1.0, edgecolors='white', linewidths=0.5, zorder=2)
        print(f"  Plotted {len(gene_x)} gene nodes (gray circles)")
    
    tf_x = [pos[tf][0] for tf in tf_list_sorted if tf in pos]
    tf_y = [pos[tf][1] for tf in tf_list_sorted if tf in pos]
    tf_c = [tf_colors[tf] for tf in tf_list_sorted if tf in pos]
    
    if tf_x:
        ax.scatter(tf_x, tf_y, s=tf_node_size, c=tf_c,
                  alpha=1.0, edgecolors='black', linewidths=1.5, zorder=3)
        print(f"  Plotted {len(tf_x)} TF nodes")
    
    if show_tf_labels:
        for tf in tf_list_sorted:
            if tf in pos:
                node_x, node_y = pos[tf]
                angle = np.arctan2(node_y - center_y, node_x - center_x)
                label_offset = 15
                offset_x = label_offset * np.cos(angle)
                offset_y = label_offset * np.sin(angle)
                
                ax.annotate(
                    tf,
                    xy=(node_x, node_y),
                    xytext=(offset_x, offset_y),
                    textcoords='offset points',
                    fontsize=8,
                    color='#ff7f00',
                    fontweight='bold',
                    ha='center',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                             edgecolor='none', alpha=0.8),
                    zorder=4
                )
    
    plot_margin = 0.2
    ax.set_xlim(center_x - outer_radius - plot_margin, 
                center_x + outer_radius + plot_margin)
    ax.set_ylim(center_y - outer_radius - plot_margin, 
                center_y + outer_radius + plot_margin)
    ax.set_aspect('equal')  # Keep the round shape without deformation
    ax.axis('off')
    
    if title:
        ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    else:
        ax.set_title('Gene Regulatory Network', 
                    fontsize=18, fontweight='bold', pad=20)
    
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=TF_FIXED_COLOR,
               markersize=14, label='Transcription Factor',
               markeredgecolor='black', markeredgewidth=1.5),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=TARGET_FIXED_COLOR,
               markersize=12, label='Target Gene',
               markeredgecolor='black', markeredgewidth=0.5)
    ]
    ax.legend(handles=legend_elements, loc='upper left', frameon=True, 
         facecolor='white', edgecolor='black', framealpha=0.9, fontsize=12,
         markerscale=1.2, handlelength=2, handleheight=1.5)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
        print(f" Save network diagram: {save_path}")
    
    print(' Circular inner/outer GRN network plot completed!')
    
    return fig




def plot_grn_network_bipartite(
    grn_network: pd.DataFrame,
    tfs_timepoint: Dict[str, float],
    tf_node_size: float = 500,
    gene_node_size: float = 300,
    tf_layer_y_center: float = 1.5,
    gene_layer_y_center: float = 0.5,
    tf_layer_height: float = 0.2,  # TF layer vertical range
    tf_layer_width: float = 1.5,  # TF layer horizontal range
    gene_layer_height: float = 0.3,  # Gene layer vertical range (can be larger)
    gene_layer_width: float = 2.0,  # Gene layer horizontal range (can be larger)
    show_tf_labels: bool = True,
    figsize: Tuple[int, int] = (20, 10),
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    edge_alpha: float = 0.15,
    edge_width: float = 0.5,
    seed: int = 42
) -> plt.Figure:
    """
    Draw a gene regulatory network diagram of a two-layer structure (organic layout).
    
    Upper layer: TF nodes (randomly distributed within the elliptical area, fixed light blue)
    Lower layer: Gene nodes (randomly distributed within the elliptical area, uniform gray circles)
    
    Args:
        grn_network: GRN network data, including 'tf', 'gene', 'correlation' columns
        tfs_timepoint: reserved for compatibility with older calls; does not currently participate in TF sorting or color mapping.
        tf_node_size: TF node size (unified)
        gene_node_size: Gene node size (unified)
        tf_layer_y_center: center y coordinate of TF layer
        gene_layer_y_center: the center y coordinate of the Gene layer
        tf_layer_height: vertical dispersion range within the TF layer
        tf_layer_width: horizontal dispersion range within the TF layer
        gene_layer_height: vertical dispersion range within the Gene layer
        gene_layer_width: horizontal dispersion range within the Gene layer
        show_tf_labels: whether to display TF labels
        figsize: graphic size
        title: graphic title
        save_path: save path
        edge_alpha: edge transparency
        edge_width: width of edge
        seed: random seed
        
    Returns:
        matplotlib.Figure: drawn graphics object
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import networkx as nx
    from matplotlib.patches import Ellipse
    
    print('Draw a double-layer organic layout GRN network diagram...')
    
    np.random.seed(seed)
    
    df_grn = grn_network.copy()
    tf_list = df_grn['tf'].unique().tolist()
    gene_list = df_grn['gene'].unique().tolist()
    gene_list = [g for g in gene_list if g not in tf_list]
    
    print(f"  The network contains {len(tf_list)} TFs, {len(gene_list)} genes")
    print(f"  Total {len(df_grn)} regulatory relationships")
    
    G = nx.DiGraph()
    for _, row in df_grn.iterrows():
        tf = row['tf']
        gene = row['gene']
        # weight = row.get('correlation', row.get('weights', row.get('weight', 1.0)))
        # G.add_edge(tf, gene, weight=weight)
        if tf in tf_list and gene in gene_list:
            weight = row.get('correlation', row.get('weights', row.get('weight', 1.0)))
            G.add_edge(tf, gene, weight=weight)
    print(f"  Added {G.number_of_edges()} TF->Gene regulatory edges")

    # tfs_timepoint is ignored so TF layout and colors are not pseudotime-driven.
    tf_list_sorted = list(tf_list)
    
    pos = {}
    
    n_tfs = len(tf_list_sorted)
    
    for i, tf in enumerate(tf_list_sorted):
        if n_tfs == 1:
            base_x = 0.5
        else:
            base_x = 0.5 + (i / (n_tfs - 1) - 0.5) * tf_layer_width * 0.85
        
        x_jitter = np.random.uniform(-0.03, 0.03)
        y_jitter = np.random.uniform(-tf_layer_height/2, tf_layer_height/2)
        
        x = base_x + x_jitter
        y = tf_layer_y_center + y_jitter
        
        pos[tf] = (x, y)
    
        # 4.2 Lower layer: Gene nodes are randomly distributed within the elliptical area of the Gene layer (similar to the TF layer)
    n_genes = len(gene_list)
    
    if n_genes == 1:
        pos[gene_list[0]] = (0.5, gene_layer_y_center)
    else:
        a = gene_layer_width / 2  # Gene layer horizontal semi-axis
        b = gene_layer_height / 2  # Gene layer vertical semi-axis
        
        # Randomly assign positions to each gene node (evenly distributed within the ellipse)
        for i, gene in enumerate(gene_list):
            # This ensures that the nodes are evenly distributed within the ellipse and are not crowded.
            
            base_angle = 2 * np.pi * i / n_genes  # Base angle (uniform distribution)
            angle_jitter = np.random.uniform(-np.pi/n_genes*0.8, np.pi/n_genes*0.8) if n_genes > 1 else 0
            angle = base_angle + angle_jitter
            
            # Radius: random, but biased toward the periphery of the ellipse (to avoid overcrowding in the center)
            r = np.sqrt(np.random.uniform(0.0, 1.0))  # Square root distribution of 0.3-1.0
            
            x = 0.5 + r * a * np.cos(angle)
            y = gene_layer_y_center + r * b * np.sin(angle)
            
            pos[gene] = (x, y)
    
    print(f"  Organic layout completed: TF layer ({tf_layer_width}x{tf_layer_height}), Gene layer ({gene_layer_width}x{gene_layer_height})")

    tf_colors = {tf: TF_FIXED_COLOR for tf in tf_list_sorted}
    
    print('  Draw network diagram...')
    fig, ax = plt.subplots(figsize=figsize)
    
    tf_ellipse = Ellipse(
        xy=(0.5, tf_layer_y_center),
        width=tf_layer_width * 1.05,  # Use TF layer width
        height=tf_layer_height * 2.0,  # Use TF layer height
        fill=False,
        edgecolor='lightblue',
        linestyle='--',
        linewidth=1.5,
        alpha=0.4,
        zorder=0
    )
    ax.add_patch(tf_ellipse)
    
    gene_ellipse = Ellipse(
        xy=(0.5, gene_layer_y_center),
        width=gene_layer_width * 1.05,  # Use Gene layer width
        height=gene_layer_height * 2.0,  # Use Gene layer height
        fill=False,
        edgecolor='lightgreen',
        linestyle='--',
        linewidth=1.5,
        alpha=0.4,
        zorder=0
    )
    ax.add_patch(gene_ellipse)
    
    for u, v in G.edges():
        if u in pos and v in pos:
            x = [pos[u][0], pos[v][0]]
            y = [pos[u][1], pos[v][1]]
            ax.plot(x, y, color=EDGE_FIXED_COLOR, alpha=edge_alpha, linewidth=edge_width, zorder=1)
    
    print(f"  Plotted {G.number_of_edges()} edges")
    
    gene_x = [pos[g][0] for g in gene_list if g in pos]
    gene_y = [pos[g][1] for g in gene_list if g in pos]
    
    if gene_x:
        ax.scatter(gene_x, gene_y, s=gene_node_size, c=TARGET_FIXED_COLOR,
                  alpha=1.0, edgecolors='white', linewidths=0.5, zorder=2)
        print(f"  Plotted {len(gene_x)} gene nodes (gray circles)")
    
    tf_x = [pos[tf][0] for tf in tf_list_sorted if tf in pos]
    tf_y = [pos[tf][1] for tf in tf_list_sorted if tf in pos]
    tf_c = [tf_colors[tf] for tf in tf_list_sorted if tf in pos]
    
    if tf_x:
        ax.scatter(tf_x, tf_y, s=tf_node_size, c=tf_c,
                  alpha=1.0, edgecolors='black', linewidths=1.5, zorder=3)
        print(f"  Plotted {len(tf_x)} TF nodes")
    
    if show_tf_labels:
        for tf in tf_list_sorted:
            if tf in pos:
                ax.annotate(
                    tf,
                    xy=(pos[tf][0], pos[tf][1]),
                    xytext=(0, 8),
                    textcoords='offset points',
                    fontsize=8,
                    color='#ff7f00',
                    fontweight='bold',
                    ha='center',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                             edgecolor='none', alpha=0.8),
                    zorder=4
                )
    
    max_width = max(tf_layer_width, gene_layer_width)
    ax.set_xlim(0.5 - max_width/2 - 0.3, 0.5 + max_width/2 + 0.3)

    ax.set_ylim(gene_layer_y_center - 0.4, tf_layer_y_center + 0.4)
    # y_min = gene_layer_y_center - gene_layer_height - y_margin
    # y_max = tf_layer_y_center + tf_layer_height + y_margin
    # ax.set_ylim(y_min, y_max)
    # #--------------------------------
    ax.set_aspect('auto')
    ax.axis('off')
    
    if title:
        ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    else:
        ax.set_title('Gene Regulatory Network (Bipartite Organic Layout)', 
                    fontsize=16, fontweight='bold', pad=20)
    
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=TF_FIXED_COLOR,
               markersize=14, label='Transcription Factor',
               markeredgecolor='black', markeredgewidth=1.5),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=TARGET_FIXED_COLOR,
               markersize=12, label='Target Gene',
               markeredgecolor='black', markeredgewidth=0.5)
    ]
    ax.legend(handles=legend_elements, loc='upper left', frameon=True, 
         facecolor='white', edgecolor='black', framealpha=0.9, fontsize=12,
         markerscale=1.2, handlelength=2, handleheight=1.5)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
        print(f" Save network diagram: {save_path}")
    
    print(' Two-layer organic-layout GRN network plot completed!')
    
    return fig



# def plot_grn_network_bipartite(
#     grn_network: pd.DataFrame,
#     tfs_timepoint: Dict[str, float],
#     tf_node_size: float = 500,
#     gene_node_size: float = 300,
#     tf_layer_y: float = 1.0,
#     gene_layer_y: float = 0.0,
#     vertical_spacing: float = 0.5,
#     show_tf_labels: bool = True,
#     show_gene_labels: bool = False,
#     figsize: Tuple[int, int] = (16, 10),
#     title: Optional[str] = None,
#     save_path: Optional[str] = None,
#     edge_alpha: float = 0.2,
#     edge_width: float = 0.5
# ) -> plt.Figure:
#     """
    
    
#     Args:
        
#     Returns:
#     """
#     import numpy as np
#     import pandas as pd
#     import matplotlib.pyplot as plt
#     import networkx as nx
#     from matplotlib.colors import LinearSegmentedColormap, Normalize
    
    
#     df_grn = grn_network.copy()
#     tf_list = df_grn['tf'].unique().tolist()
#     gene_list = df_grn['gene'].unique().tolist()
# # Remove genes that are also TF from the gene list (to avoid duplication)
#     gene_list = [g for g in gene_list if g not in tf_list]
    
    
#     G = nx.DiGraph()
#     for _, row in df_grn.iterrows():
#         tf = row['tf']
#         gene = row['gene']
#         weight = row.get('correlation', row.get('weights', row.get('weight', 1.0)))
#         G.add_edge(tf, gene, weight=weight)
    
    
# # Only keep TF with time point information
#     tfs_with_time = [tf for tf in tf_list if tf in tfs_timepoint]
#     tfs_without_time = [tf for tf in tf_list if tf not in tfs_timepoint]
    
#     tfs_with_time_sorted = sorted(tfs_with_time, key=lambda x: tfs_timepoint[x])
    
#     tf_list_sorted = tfs_with_time_sorted + tfs_without_time
    
#     if tfs_with_time_sorted:
#         time_range = [tfs_timepoint[tfs_with_time_sorted[0]], 
#                      tfs_timepoint[tfs_with_time_sorted[-1]]]
    
#     pos = {}
    
#     n_tfs = len(tf_list_sorted)
#     if n_tfs == 1:
#         tf_x_positions = [0.5]
#     else:
#         tf_x_positions = np.linspace(0, 1, n_tfs)
    
#     for i, tf in enumerate(tf_list_sorted):
#         pos[tf] = (tf_x_positions[i], tf_layer_y)
    
#     n_genes = len(gene_list)
#     if n_genes == 0:
#         gene_x_positions = []
#     elif n_genes == 1:
#         gene_x_positions = [0.5]
#     else:
#         gene_x_positions = np.linspace(0, 1, n_genes)
    
#     for i, gene in enumerate(gene_list):
#         pos[gene] = (gene_x_positions[i], gene_layer_y)
    
    
#     colors_base = [
#     ]
#     cmap_blueyellow = LinearSegmentedColormap.from_list('blueyellow', colors_base, N=256)
    
#     tf_colors = {}
#     if tfs_with_time_sorted:
#         timepoint_values = [tfs_timepoint[tf] for tf in tfs_with_time_sorted]
#         min_time = min(timepoint_values)
#         max_time = max(timepoint_values)
#         norm = Normalize(vmin=min_time, vmax=max_time)
        
#         for tf in tfs_with_time_sorted:
#             time_val = tfs_timepoint[tf]
#             normalized_val = norm(time_val)
#             rgba = cmap_blueyellow(normalized_val)
#             hex_color = f"#{int(rgba[0]*255):02x}{int(rgba[1]*255):02x}{int(rgba[2]*255):02x}"
#             tf_colors[tf] = hex_color
    
#     for tf in tfs_without_time:
#         tf_colors[tf] = "#808080"
    
    
#     fig, ax = plt.subplots(figsize=figsize)
    
#     for u, v in G.edges():
#         x = [pos[u][0], pos[v][0]]
#         y = [pos[u][1], pos[v][1]]
#         ax.plot(x, y, color='gray', alpha=edge_alpha, linewidth=edge_width, zorder=1)
    
    
#     gene_x = [pos[g][0] for g in gene_list if g in pos]
#     gene_y = [pos[g][1] for g in gene_list if g in pos]
    
#     if gene_x:
#         ax.scatter(gene_x, gene_y, s=gene_node_size, c='#D3D3D3', 
#                   alpha=1.0, edgecolors='black', linewidths=0.5, zorder=2, label='Target Gene')
    
#     tf_x = [pos[tf][0] for tf in tf_list_sorted if tf in pos]
#     tf_y = [pos[tf][1] for tf in tf_list_sorted if tf in pos]
#     tf_c = [tf_colors[tf] for tf in tf_list_sorted if tf in pos]
    
#     if tf_x:
#         ax.scatter(tf_x, tf_y, s=tf_node_size, c=tf_c,
#                   alpha=1.0, edgecolors='black', linewidths=1.0, zorder=3, label='Transcription Factor')
    
#     if show_tf_labels:
#         for tf in tf_list_sorted:
#             if tf in pos:
#                 ax.annotate(
#                     tf,
#                     xy=(pos[tf][0], pos[tf][1]),
#                     xytext=(0, 8),
#                     textcoords='offset points',
#                     fontsize=8,
#                     color='#ff7f00',
#                     fontweight='bold',
#                     ha='center',
#                     bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
#                              edgecolor='none', alpha=0.7),
#                     zorder=4
#                 )
    
#     if show_gene_labels:
#         for gene in gene_list:
#             if gene in pos:
#                 ax.annotate(
#                     gene,
#                     xy=(pos[gene][0], pos[gene][1]),
#                     xytext=(0, -12),
#                     textcoords='offset points',
#                     fontsize=6,
#                     color='black',
#                     ha='center',
#                     alpha=0.7,
#                     zorder=4
#                 )
    
#     ax.set_xlim(-0.05, 1.05)
#     ax.set_ylim(gene_layer_y - 0.2, tf_layer_y + 0.2)
#     ax.set_aspect('auto')
#     ax.axis('off')
    
#     if title:
#         ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
#     else:
#         ax.set_title('Gene Regulatory Network (Bipartite Layout)', 
#                     fontsize=16, fontweight='bold', pad=20)
    
#     if tfs_with_time_sorted:
#         from matplotlib.colorbar import ColorbarBase
        
#         cax = fig.add_axes([0.92, 0.3, 0.02, 0.4])
        
#         norm = Normalize(vmin=min_time, vmax=max_time)
#         cb = ColorbarBase(cax, cmap=cmap_blueyellow, norm=norm, orientation='vertical')
#         cb.set_label('TF Pseudotime', rotation=270, labelpad=20, fontsize=10)
    
#     from matplotlib.lines import Line2D
#     legend_elements = [
#         Line2D([0], [0], marker='o', color='w', markerfacecolor='#ff7f00', 
#                markersize=10, label='Transcription Factor', markeredgecolor='black'),
#         Line2D([0], [0], marker='o', color='w', markerfacecolor='#D3D3D3', 
#                markersize=8, label='Target Gene', markeredgecolor='black')
#     ]
#     ax.legend(handles=legend_elements, loc='upper left', frameon=True, 
#              facecolor='white', edgecolor='black', framealpha=0.9)
    
#     plt.tight_layout()
    
#     if save_path:
#         plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    
#     return fig
