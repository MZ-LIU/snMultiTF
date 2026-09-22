import os
os.environ["OMP_NUM_THREADS"] = "2"

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from ensemble_clustering.tool.metric_calculation import calculate_metrics


def csp_seurat_from_array(data_embed, resolution=0.8, n_neighbors=15, method='leiden', random_state=42):
    if method not in ['leiden', 'louvain']:
        raise ValueError("method can only be 'leiden' or 'louvain'")

    X_embed = np.asarray(data_embed)
    if X_embed.ndim != 2:
        raise ValueError('The input data must be a two-dimensional matrix (n_cells, n_features)')

    adata = AnnData(X=np.zeros((X_embed.shape[0], 1), dtype=np.float32))
    adata.obsm['X_pca'] = X_embed
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X_pca', random_state=random_state)

    cluster_key = f'{method}_r{resolution:g}'
    if method == 'leiden':
        sc.tl.leiden(adata,resolution=resolution,key_added=cluster_key,random_state=random_state)
    else:
        sc.tl.louvain(adata, resolution=resolution, key_added=cluster_key, random_state=random_state)

    sc.tl.umap(adata, random_state=random_state)

    labels_pred = adata.obs[cluster_key].astype(int).to_numpy()
    embed2 = adata.obsm['X_umap']
    return labels_pred, cluster_key, embed2


def infer_space_name(file_stem):
    name = file_stem.lower()
    if 'pca' in name:
        digits = ''.join([c for c in name if c.isdigit()])
        return f'PCA{digits}' if digits else 'PCA'
    if 'umap' in name:
        digits = ''.join([c for c in name if c.isdigit()])
        return f'UMAP{digits}' if digits else 'UMAP'
    return 'EMBED'


def append_metrics_to_excel(excel_path, row_dict):
    columns_order = [
        'input space',
        'algorithm',
        'K value',
        "DBI ↓",
        "SC ↑",
        "CH ↑",
        "Dunn Index ↑",
        "S_Dbw ↓",
        'Maximum cluster proportion (%)',
    ]

    if excel_path.exists():
        df_old = pd.read_excel(excel_path)
    else:
        df_old = pd.DataFrame(columns=columns_order)

    for c in columns_order:
        if c not in df_old.columns:
            df_old[c] = np.nan
    df_old = df_old[columns_order]

    df_new = pd.concat([df_old, pd.DataFrame([row_dict], columns=columns_order)], ignore_index=True)
    try:
        df_new.to_excel(excel_path, index=False)
        return excel_path
    except PermissionError:
        fallback = excel_path.with_name(excel_path.stem + "_seurat_append.xlsx")
        df_new.to_excel(fallback, index=False)
        return fallback


if __name__ == '__main__':
    script_dir = Path(__file__).resolve().parent
    data_dir = script_dir / 'data'
    results_dir = script_dir / 'results'
    results_dir.mkdir(parents=True, exist_ok=True)

    # data_path = data_dir / 'data_pca50.npy'
    data_path = data_dir / 'data_umap2.npy'
    cell_index_path = data_dir / 'cell_index.npy'
    metrics_xlsx = results_dir / 'metrics_comparison_leiden.xlsx'

    method = 'leiden'
    resolution = 0.8
    n_neighbors = 15
    random_state = 42

    time0 = time.time()

    if data_path.exists():
        data_eval = np.load(data_path)
        if data_eval.ndim != 2:
            raise ValueError('The input data must be a two-dimensional matrix (n_cells, n_features)')

        labels, cluster_key, embed2 = csp_seurat_from_array(
            data_eval,
            resolution=resolution,
            n_neighbors=n_neighbors,
            method=method,
            random_state=random_state,
        )

        out_csv = data_dir / f'{data_path.stem}_{cluster_key}_labels.csv'
        out_embed = data_dir / f'{data_path.stem}_{cluster_key}_umap2.npy'

        if cell_index_path.exists():
            cell_ids = np.load(cell_index_path, allow_pickle=True)
            if len(cell_ids) == len(labels):
                df = pd.DataFrame({'cell_index': cell_ids, 'cluster': labels})
            else:
                df = pd.DataFrame({'cluster': labels})
        else:
            df = pd.DataFrame({'cluster': labels})

        df.to_csv(out_csv, index=False, encoding='utf-8-sig')
        np.save(out_embed, embed2)

        k, dbi, scv, ch, dunn, s_dbw, noise_ratio, max_cluster_ratio = calculate_metrics(data_eval, labels)

        space_name = infer_space_name(data_path.stem)
        algo_name = f"Seurat-{method.capitalize()}"
        row = {
            'input space': space_name,
            'algorithm': algo_name,
            'K value': int(k),
            "DBI ↓": float(dbi) if np.isfinite(dbi) else np.nan,
            "SC ↑": float(scv) if np.isfinite(scv) else np.nan,
            "CH ↑": float(ch) if np.isfinite(ch) else np.nan,
            "Dunn Index ↑": float(dunn) if np.isfinite(dunn) else np.nan,
            "S_Dbw ↓": float(s_dbw) if np.isfinite(s_dbw) else np.nan,
            'Maximum cluster proportion (%)': float(max_cluster_ratio),
        }
        metrics_saved_path = append_metrics_to_excel(metrics_xlsx, row)

        print(type(labels))
        print(labels[:10])
        print('cluster_key:', cluster_key)
        print(
            f"params -> resolution={resolution}, n_neighbors={n_neighbors}, method={method}, random_state={random_state}"
        )
        print(
            f"metrics -> K={k}, DBI={dbi:.4f}, SC={scv:.4f}, CH={ch:.2f}, "
            f"Dunn={dunn:.4f}, S_Dbw={s_dbw:.4f}, Noise={noise_ratio:.2f}%, MaxCluster={max_cluster_ratio:.2f}%"
        )
        print('saved:', out_csv)
        print('saved:', out_embed)
        print('saved:', metrics_saved_path)
    else:
        print('Input file not found:', data_path)

    time1 = time.time()
    print(time1 - time0)
