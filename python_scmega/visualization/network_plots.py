"""
Network Analysis and Visualization Tools

This module provides supplementary network analysis and visualization functions
for gene regulatory networks. For the main GRN visualization (R-compatible),
use plot_grn_network from network_inference.network_construction module.

Provided functions:
- plot_network_centrality: Centrality analysisNetwork centrality analysis
- plot_network_topology: Topology analysisNetwork topology analysis
- plot_tf_targets: TF-specific subnetworkTF-specific subnetwork
- create_network_layout: Layout utilities layout tool
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from typing import Optional, Dict, Any, List, Union, Tuple



def plot_network_centrality(
    grn_network: pd.DataFrame,
    centrality_metrics: List[str] = ["degree", "betweenness", "closeness", "eigenvector"],
    top_k: int = 20,
    figsize: Tuple[int, int] = (15, 10),
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    **kwargs
) -> plt.Figure:
    """
    Visualize network centrality analysis.
    
    Args:
        grn_network: GRN network DataFrame
        centrality_metrics: List of centrality metrics to compute
        top_k: Number of top nodes to display
        figsize: Figure size
        title: Plot title
        save_path: Path to save figure
        **kwargs: Additional parameters
        
    Returns:
        matplotlib.Figure: Centrality analysis visualization
    """
    
    # Create graph
    G = nx.DiGraph()
    for _, row in grn_network.iterrows():
        G.add_edge(row['TF'], row['target_gene'], 
                  weight=row.get('correlation', 1.0))
    
    # Calculate centrality metrics
    centrality_data = {}
    
    for metric in centrality_metrics:
        if metric == "degree":
            centrality_data[metric] = pd.Series(dict(G.degree()))
        elif metric == "betweenness":
            centrality_data[metric] = pd.Series(
                nx.betweenness_centrality(G, weight='weight', normalized=True)
            )
        elif metric == "closeness":
            centrality_data[metric] = pd.Series(
                nx.closeness_centrality(G, distance='weight')
            )
        elif metric == "eigenvector":
            try:
                centrality_data[metric] = pd.Series(
                    nx.eigenvector_centrality(G, weight='weight', max_iter=1000)
                )
            except:
                centrality_data[metric] = pd.Series(
                    nx.degree_centrality(G)
                )
        elif metric == "pagerank":
            centrality_data[metric] = pd.Series(
                nx.pagerank(G, weight='weight')
            )
    
    df_centrality = pd.DataFrame(centrality_data)
    
    # Create subplots
    n_metrics = len(centrality_metrics)
    n_cols = min(n_metrics, 2)
    n_rows = (n_metrics + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    if n_metrics == 1:
        axes = [axes]
    else:
        axes = axes.flatten() if n_rows > 1 else axes
    
    # Plot each metric
    for i, metric in enumerate(centrality_metrics):
        ax = axes[i] if n_metrics > 1 else axes[0]
        
        # Get top-k nodes
        top_nodes = df_centrality[metric].nlargest(top_k)
        
        # Bar plot
        ax.barh(range(len(top_nodes)), top_nodes.values, color='steelblue', alpha=0.7)
        ax.set_yticks(range(len(top_nodes)))
        ax.set_yticklabels(top_nodes.index, fontsize=9)
        ax.set_xlabel(f'{metric.capitalize()} Centrality', fontsize=10)
        ax.set_title(f'Top {top_k} Nodes by {metric.capitalize()}', 
                    fontsize=11, fontweight='bold')
        ax.grid(axis='x', alpha=0.3)
    
    # Hide unused subplots
    for i in range(n_metrics, len(axes)):
        axes[i].set_visible(False)
    
    if title:
        fig.suptitle(title, fontsize=14, fontweight='bold', y=0.995)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig


def plot_network_topology(
    grn_network: pd.DataFrame,
    analysis_type: str = "degree_distribution",
    figsize: Tuple[int, int] = (12, 8),
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    **kwargs
) -> plt.Figure:
    """
    Visualize network topology properties.
    
    Args:
        grn_network: GRN network DataFrame
        analysis_type: Type of analysis ('degree_distribution', 'clustering', 'path_length')
        figsize: Figure size
        title: Plot title
        save_path: Path to save figure
        **kwargs: Additional parameters
        
    Returns:
        matplotlib.Figure: Topology visualization
    """
    
    # Create graph
    G = nx.DiGraph()
    for _, row in grn_network.iterrows():
        G.add_edge(row['TF'], row['target_gene'])
    
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    axes = axes.flatten()
    
    # 1. Degree distribution
    degrees = [d for n, d in G.degree()]
    axes[0].hist(degrees, bins=30, alpha=0.7, color='steelblue', edgecolor='black')
    axes[0].set_xlabel('Degree', fontsize=10)
    axes[0].set_ylabel('Frequency', fontsize=10)
    axes[0].set_title('Degree Distribution', fontsize=11, fontweight='bold')
    axes[0].grid(alpha=0.3)
    
    # 2. In-degree vs Out-degree
    in_degrees = [d for n, d in G.in_degree()]
    out_degrees = [d for n, d in G.out_degree()]
    axes[1].scatter(in_degrees, out_degrees, alpha=0.5, s=30)
    axes[1].set_xlabel('In-degree', fontsize=10)
    axes[1].set_ylabel('Out-degree', fontsize=10)
    axes[1].set_title('In-degree vs Out-degree', fontsize=11, fontweight='bold')
    axes[1].grid(alpha=0.3)
    
    # 3. Clustering coefficient distribution
    G_undirected = G.to_undirected()
    clustering_coeffs = list(nx.clustering(G_undirected).values())
    axes[2].hist(clustering_coeffs, bins=30, alpha=0.7, color='coral', edgecolor='black')
    axes[2].set_xlabel('Clustering Coefficient', fontsize=10)
    axes[2].set_ylabel('Frequency', fontsize=10)
    axes[2].set_title('Clustering Coefficient Distribution', fontsize=11, fontweight='bold')
    axes[2].grid(alpha=0.3)
    
    # 4. Network statistics
    axes[3].axis('off')
    stats_text = f"""
    Network Statistics:
    
    Nodes: {G.number_of_nodes()}
    Edges: {G.number_of_edges()}
    Density: {nx.density(G):.4f}
    Avg Degree: {np.mean(degrees):.2f}
    Avg Clustering: {np.mean(clustering_coeffs):.4f}
    """
    axes[3].text(0.1, 0.5, stats_text, fontsize=11, 
                verticalalignment='center',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    if title:
        fig.suptitle(title, fontsize=14, fontweight='bold', y=0.995)
    else:
        fig.suptitle('Network Topology Analysis', fontsize=14, fontweight='bold', y=0.995)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig


def plot_tf_targets(
    grn_network: pd.DataFrame,
    tf_name: str,
    layout: str = "circular",
    top_k: int = 20,
    figsize: Tuple[int, int] = (10, 10),
    title: Optional[str] = None,
    save_path: Optional[str] = None,
    **kwargs
) -> plt.Figure:
    """
    Visualize a TF and its target genes (ego network).
    
    Args:
        grn_network: GRN network DataFrame
        tf_name: TF name to visualize
        layout: Layout algorithm
        top_k: Number of top targets to show
        figsize: Figure size
        title: Plot title
        save_path: Path to save figure
        **kwargs: Additional parameters
        
    Returns:
        matplotlib.Figure: TF-target visualization
    """
    
    # Filter edges for this TF
    tf_edges = grn_network[grn_network['TF'] == tf_name].copy()
    
    if tf_edges.empty:
        raise ValueError(f"No targets found for TF: {tf_name}")
    
    # Get top-k targets by correlation
    if 'correlation' in tf_edges.columns:
        tf_edges = tf_edges.nlargest(top_k, 'correlation')
    else:
        tf_edges = tf_edges.head(top_k)
    
    # Create graph
    G = nx.DiGraph()
    for _, row in tf_edges.iterrows():
        G.add_edge(row['TF'], row['target_gene'], 
                  weight=row.get('correlation', 1.0))
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Layout
    if layout == 'circular':
        pos = nx.circular_layout(G)
    elif layout == 'spring':
        pos = nx.spring_layout(G, weight='weight', iterations=100)
    else:
        pos = nx.spring_layout(G)
    
    # Draw
    # TF node
    nx.draw_networkx_nodes(G, pos, nodelist=[tf_name], 
                          node_color='red', node_size=1000, 
                          alpha=1.0, ax=ax)
    
    # Target nodes
    targets = [n for n in G.nodes() if n != tf_name]
    nx.draw_networkx_nodes(G, pos, nodelist=targets,
                          node_color='lightblue', node_size=500,
                          alpha=0.7, ax=ax)
    
    # Edges
    nx.draw_networkx_edges(G, pos, alpha=0.5, edge_color='gray',
                          width=2, arrows=True, arrowsize=20, ax=ax)
    
    # Labels
    nx.draw_networkx_labels(G, pos, font_size=9, font_color='black', ax=ax)
    
    ax.set_title(title or f'{tf_name} Regulatory Targets (Top {top_k})',
                fontsize=14, fontweight='bold')
    ax.axis('off')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig


def create_network_layout(
    grn_network: pd.DataFrame,
    layout_type: str = "spring",
    **kwargs
) -> Dict[str, Tuple[float, float]]:
    """
    Create network layout coordinates.
    
    Args:
        grn_network: GRN network DataFrame
        layout_type: Layout algorithm
        **kwargs: Additional parameters for layout algorithm
        
    Returns:
        Dict[str, Tuple[float, float]]: Node positions
    """
    
    G = nx.DiGraph()
    for _, row in grn_network.iterrows():
        G.add_edge(row['TF'], row['target_gene'], 
                  weight=row.get('correlation', 1.0))
    
    seed = kwargs.get('seed', 42)
    np.random.seed(seed)
    
    if layout_type == 'spring' or layout_type == 'fruchterman_reingold':
        pos = nx.spring_layout(G, weight='weight', 
                              iterations=kwargs.get('iterations', 1000),
                              seed=seed)
    elif layout_type == 'kamada_kawai':
        pos = nx.kamada_kawai_layout(G, weight='weight')
    elif layout_type == 'circular':
        pos = nx.circular_layout(G)
    elif layout_type == 'spectral':
        pos = nx.spectral_layout(G)
    else:
        pos = nx.spring_layout(G, seed=seed)
    
    return pos


# Export main functions
__all__ = [
    'plot_network_centrality',
    'plot_network_topology',
    'plot_tf_targets',
    'create_network_layout',
]
