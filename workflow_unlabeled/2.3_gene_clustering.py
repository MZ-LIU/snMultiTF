import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

# sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'jullei'))
# from _01_05_01_05_ACM_WEP import cluster_demo_pro
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, ".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from ensemble_clustering.jullei._01_05_01_05_ACM_WEP import  cluster_demo_pro

BASE_DIR   = os.path.join(current_dir, "results_sc", "sc_multiome", "Clustering_results", "k_8")
CSV_PATH   = os.path.join(BASE_DIR, "gene_stability_result.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, 'Gene_clustering_results')
CLUSTER_COUNT = 5  # Number of clusters (i.e. number of edges removed from MST)

df = pd.read_csv(CSV_PATH)
print(f"{len(df)} genes were read")

scores = df['mean_score_when_detected'].values

data_2d = scores.reshape(-1, 1)

print(f"Start clustering, number of clusters = {CLUSTER_COUNT}...")
labels, beilv = cluster_demo_pro(data_2d, CLUSTER_COUNT)

df['cluster'] = labels

df = df.sort_values(by=['cluster', 'mean_score_when_detected'], ascending=[True, False])

print('\n===== Clustering result statistics =====')
for cid in sorted(df['cluster'].unique()):
    cluster_genes = df[df['cluster'] == cid]
    scores_in_cluster = cluster_genes['mean_score_when_detected']
    print(f"Cluster {cid}: {len(cluster_genes)} genes, "
          f"mean_score range [{scores_in_cluster.min():.4f}, {scores_in_cluster.max():.4f}]")

os.makedirs(OUTPUT_DIR, exist_ok=True)
output_path = os.path.join(OUTPUT_DIR, 'gene_clustering_result.csv')
df[['gene', 'mean_score_when_detected', 'cluster']].to_csv(output_path, index=False, encoding='utf-8-sig')
print(f"\n clustering results have been saved to: {output_path}")

plt.figure(figsize=(10, 6))
plt.plot(range(1, len(beilv) + 1), beilv, marker='o', markersize=4)
plt.xlabel('Iteration')
plt.ylabel('Beilv')
plt.title('Beilv Curve')
plt.grid(True, alpha=0.3)
plt.tight_layout()
beilv_path = os.path.join(OUTPUT_DIR, f"beilv_curve_{CLUSTER_COUNT}.png")
plt.savefig(beilv_path, dpi=150)
plt.close()
print(f"beilv curve graph has been saved to: {beilv_path}")
