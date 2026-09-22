"""
4_celltype.py
================
Cell type annotation process (using RAG system)

Function description:
  1. Read the differential gene list and RAG database
  2. Automatic cell type annotation using the RAG system
  3. Get annotation results

How to run:
  python 4_celltype.py
"""

import os
import sys
import warnings
warnings.filterwarnings('ignore')

from pathlib import Path

import pandas as pd

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import scanpy as sc
import muon as mu
import anndata as ad

import colorcet as cc
import matplotlib.colors as mcolors
from dotenv import load_dotenv

# ============================================================
# ============================================================
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, ".."))

if project_root not in sys.path:
    sys.path.insert(0, project_root)
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

from python_scmega.cell_annotation.cell_typ import CellTypeAnnotator


# ============================================================
# ============================================================
CLUSTER_BASE_DIR = os.path.join(current_dir, "results_sc", "sc_multiome", "Clustering_results")
K_VALUE = 8  # K value used

JSON_DB_PATH = os.path.join(project_root, "Maker_dataset", "Cell_marker_JSON.json")
MODEL_NAME = 'gpt-5'  # Modify based on actual available models
load_dotenv(os.path.join(project_root, "python_scmega", ".env"))

OUTPUT_BASE_DIR = os.path.join(current_dir, "results_sc", "sc_multiome", "celltype_annotation")
os.makedirs(OUTPUT_BASE_DIR , exist_ok=True)

Tissue = "Pancreas"

SHOW_PLOTS = False  # Set to True to automatically close after 5 seconds of display, and False to only save without displaying.


# ============================================================
# ============================================================
def make_glasbey_palette(categories):
    """
    Use Glasbey Color Matching - Highly Distinguished Color Matching Designed for a Large Range of Categories
    """
    k = len(categories)
    glasbey_colors = cc.glasbey[:k]
    return [mcolors.to_hex(color) for color in glasbey_colors]


def show_plot_if_enabled():
    """Determine whether to display images based on configuration"""
    if SHOW_PLOTS:
        plt.show(block=False)
        plt.pause(5)
        plt.close()


# ============================================================
# ============================================================
def main():
    print("="*70)
    print('6. Cell type annotation (RAG-based automatic annotation)')
    print("="*70)

    # timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    # output_dir = Path(OUTPUT_BASE_DIR) / timestamp
    # output_dir.mkdir(parents=True, exist_ok=True)
    output_dir = Path(OUTPUT_BASE_DIR)

    k_dir = Path(CLUSTER_BASE_DIR) / f"k_{K_VALUE}"

    # ============================================================
    # ============================================================
    print('\nStep 1: Prepare marker gene data for each cluster...')

    marker_file = k_dir / 'cluster_marker_genes_10.csv'
    marker_df = pd.read_csv(marker_file)
    required_cols = ["cluster", "gene"]
    missing_cols = [c for c in required_cols if c not in marker_df.columns]
    if missing_cols:
        raise ValueError(f"cluster_marker_genes_10.csv is missing a required column: {missing_cols}")

    annotation_input = []
    for cluster in sorted(marker_df['cluster'].unique(), key=lambda x: int(x)):
        cluster_genes = marker_df[marker_df['cluster'] == cluster]['gene'].tolist()
        genes_str = ', '.join(cluster_genes[:10])  # Take the first 15

        annotation_input.append({
            "cluster": str(cluster),
            "tissue": Tissue,
            "marker": genes_str,
        })

    annotation_df = pd.DataFrame(annotation_input)
    print(f"  Prepared marker data for {len(annotation_df)} clusters")

    input_file = output_dir / 'cluster_markers_for_annotation.xlsx'
    annotation_df.to_excel(input_file, index=False)
    print(f"  Saved to: {input_file}")

    # ============================================================
    # ============================================================
    print('\nStep 2: Initialize the RAG cell annotation system...')

    annotator = CellTypeAnnotator(
        json_db_path=JSON_DB_PATH,
        default_model=MODEL_NAME,
        verbose=True
    )

    # ============================================================
    # ============================================================
    print('\nStep 3: Perform batch cell type annotation...')

    output_file = output_dir / 'cluster_annotation_results.xlsx'
    try:
        annotation_results = annotator.batch_annotate(
            input_file=str(input_file),
            output_file=str(output_file),
            tissue_col='tissue',
            genes_col='marker',
            top_k=3,  # Return top3 candidates
            use_api=True,  # Annotate using API
            model=MODEL_NAME,
            save_interval=5  # Save progress every 5 clusters
        )

        print(f"\n Batch annotation completed!")
        print(f"  Results saved to: {output_file}")

    except Exception as e:
        print(f"\n Batch annotation failed: {e}")
        print('  If it is an API related error, please check the environment variable configuration')
        print('  Continue using manual annotation results...')

        annotation_results = annotation_df.copy()
        annotation_results['Predict_cell_type'] = 'Unknown'

    # ============================================================
    # Step 4: Save annotation result only; write-back is handled by 4.2_celltype_writeback.py
    # ============================================================
    if not os.path.exists(output_file):
        annotation_results.to_excel(output_file, index=False)
        print(f"  Fallback annotation result saved to: {output_file}")

    print("\nAnnotation results have been generated.")
    print(f"  Annotation result file: {output_file}")
    print('  To write back to MuData and generate mdata_with_celltype.h5mu, run: workflow_unlabeled/4.2_celltype_writeback.py')
    print("="*70)


if __name__ == "__main__":
    main()
