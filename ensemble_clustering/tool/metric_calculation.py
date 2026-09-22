import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist, pdist
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, silhouette_score


def dunn_index(data: np.ndarray, labels: np.ndarray) -> float:
    clusters = [data[labels == c] for c in np.unique(labels)]
    if len(clusters) < 2:
        return float("nan")

    max_intra = 0.0
    for cluster in clusters:
        if cluster.shape[0] >= 2:
            diameter = float(np.max(pdist(cluster, metric="euclidean")))
            max_intra = max(max_intra, diameter)
    if max_intra <= 0:
        return float("nan")

    min_inter = float("inf")
    for i in range(len(clusters)):
        for j in range(i + 1, len(clusters)):
            inter_dist = float(np.min(cdist(clusters[i], clusters[j], metric="euclidean")))
            min_inter = min(min_inter, inter_dist)

    if not np.isfinite(min_inter):
        return float("nan")
    return float(min_inter / (max_intra + 1e-12))


def s_dbw_index(data: np.ndarray, labels: np.ndarray) -> float:
    unique_labels = np.unique(labels)
    k = len(unique_labels)
    if k < 2:
        return float("nan")

    clusters = [data[labels == c] for c in unique_labels]
    centers = [np.mean(c, axis=0) for c in clusters]
    sigmas = [np.linalg.norm(np.std(c, axis=0)) for c in clusters]

    sigma_data = float(np.linalg.norm(np.std(data, axis=0)))
    if sigma_data <= 0:
        return float("nan")
    scat = float(np.mean(sigmas) / (sigma_data + 1e-12))

    stdev = float(np.sqrt(np.mean(np.square(sigmas))))
    if stdev <= 0:
        return float("nan")

    def density(points: np.ndarray, center: np.ndarray) -> int:
        d = np.linalg.norm(points - center, axis=1)
        return int(np.sum(d <= stdev))

    pair_sum = 0.0
    pairs = 0
    for i in range(k):
        for j in range(i + 1, k):
            union_ij = np.vstack([clusters[i], clusters[j]])
            mid_ij = 0.5 * (centers[i] + centers[j])
            dens_mid = density(union_ij, mid_ij)
            dens_i = density(clusters[i], centers[i])
            dens_j = density(clusters[j], centers[j])
            denom = max(dens_i, dens_j, 1)
            pair_sum += float(dens_mid / denom)
            pairs += 1

    if pairs == 0:
        return float("nan")
    dens_bw = float(pair_sum / pairs)
    return float(scat + dens_bw)


def calculate_metrics(data_eval: np.ndarray, labels_pred: np.ndarray):
    """
    Compute internal clustering metrics in the same coordinate space as model input.

    Returns:
        (k, dbi, sc, ch, dunn, s_dbw, noise_ratio, max_cluster_ratio)
    """
    labels_pred = np.asarray(labels_pred)
    total_cells = len(labels_pred)
    noise_ratio = float(np.sum(labels_pred == -1) / total_cells * 100.0)

    mask = labels_pred != -1
    labels_eval = labels_pred[mask]
    data_eval_clean = data_eval[mask]
    k = int(len(np.unique(labels_eval)))

    if k < 2:
        max_ratio = 100.0 if k == 1 else 0.0
        return k, np.nan, np.nan, np.nan, np.nan, np.nan, noise_ratio, max_ratio

    dbi = float(davies_bouldin_score(data_eval_clean, labels_eval))
    sc = float(silhouette_score(data_eval_clean, labels_eval))
    ch = float(calinski_harabasz_score(data_eval_clean, labels_eval))
    dunn = dunn_index(data_eval_clean, labels_eval)
    s_dbw = s_dbw_index(data_eval_clean, labels_eval)
    max_cluster_size = int(pd.Series(labels_eval).value_counts().max())
    max_cluster_ratio = float(max_cluster_size / total_cells * 100.0)
    return k, dbi, sc, ch, dunn, s_dbw, noise_ratio, max_cluster_ratio
