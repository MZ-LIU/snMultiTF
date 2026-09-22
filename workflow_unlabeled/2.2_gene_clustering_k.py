import sys
import os
import pandas as pd
import matplotlib.pyplot as plt

current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, ".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from ensemble_clustering.jullei._01_05_01_05_ACM_WEP import cluster_demo_pro

BASE_DIR   = os.path.join(current_dir, "results_sc", "sc_multiome", "Clustering_results", "k_8")
CSV_PATH   = os.path.join(BASE_DIR, "gene_stability_result.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, 'Gene clustering results')
CLUSTER_COUNT = 10  # Number of clusters (i.e. number of edges removed from MST)

df = pd.read_csv(CSV_PATH)
print(f"{len(df)} genes were read")

scores = df['mean_score_when_detected'].values

data_2d = scores.reshape(-1, 1)

print(f"Start clustering, number of clusters = {CLUSTER_COUNT}...")
labels, beilv = cluster_demo_pro(data_2d, CLUSTER_COUNT)

os.makedirs(OUTPUT_DIR, exist_ok=True)
plt.figure(figsize=(10, 6))
plt.plot(range(1, len(beilv) + 1), beilv, marker='o', markersize=4)
plt.xlabel('Iteration')
plt.ylabel('Beilv')
plt.title('Beilv Curve')
plt.grid(True, alpha=0.3)
plt.tight_layout()
beilv_path = os.path.join(OUTPUT_DIR, f"gene_beilv_curve_{CLUSTER_COUNT}.png")
plt.savefig(beilv_path, dpi=150)
plt.close()
print(f"beilv curve graph has been saved to: {beilv_path}")
