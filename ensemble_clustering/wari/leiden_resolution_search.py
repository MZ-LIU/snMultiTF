"""
Utilities for choosing a Leiden/Seurat resolution that matches a target K.

Search strategy:
    coarse: [0.3, 1.5] with step 0.1

If multiple tested resolutions produce the target K exactly, sort those
resolutions from small to large and choose the middle tested resolution.
For an even number of exact matches, the lower middle tested resolution is used.
"""
from __future__ import annotations

from typing import Any

import numpy as np


def find_leiden_resolution_for_k(
    data_eval: np.ndarray,
    target_k: int,
    n_neighbors: int = 20,
    random_state: int = 42,
    coarse_range: tuple[float, float] = (0.3, 1.5),
    coarse_step: float = 0.1,
    verbose: bool = False,
) -> tuple[float, int, list[dict[str, Any]]]:
    """
    Search for the Leiden resolution whose cluster count is closest to target_k.

    Returns
    -------
    best_res, actual_k, search_log
        Best resolution, its resulting cluster count, and all trial records.
    """
    from jullei._01_04_01_04_seurat import csp_seurat_from_array

    if target_k <= 0:
        raise ValueError(f"target_k must be positive, got {target_k}")

    search_log: list[dict[str, Any]] = []

    def _eval_resolution(res: float) -> int:
        try:
            labels, _, _ = csp_seurat_from_array(
                data_eval,
                resolution=round(float(res), 4),
                n_neighbors=n_neighbors,
                random_state=random_state,
            )
            labels = np.asarray(labels)
            return int(len(np.unique(labels[labels != -1])))
        except Exception as exc:
            if verbose:
                print(f"    resolution={res:.4f} failed: {exc}")
            return -1

    def _record(phase: str, res: float, actual_k: int) -> dict[str, Any]:
        diff = abs(actual_k - target_k) if actual_k != -1 else float("inf")
        return {
            "phase": phase,
            "resolution": round(float(res), 4),
            "actual_k": actual_k,
            "target_k": target_k,
            "diff": diff,
        }

    def _select_middle_exact_match(records: list[dict[str, Any]]) -> dict[str, Any] | None:
        exact_records = [record for record in records if record["actual_k"] == target_k]
        if not exact_records:
            return None
        exact_records.sort(key=lambda record: record["resolution"])
        return exact_records[(len(exact_records) - 1) // 2]

    def _best_by_diff(records: list[dict[str, Any]]) -> dict[str, Any]:
        def _sort_key(record: dict[str, Any]) -> tuple:
            diff = record["diff"]
            actual_k = record["actual_k"]
            prefer_under = (diff == 1 and actual_k != target_k - 1)
            return (diff, prefer_under, record["resolution"])
        return min(records, key=_sort_key)

    coarse_list = np.arange(
        coarse_range[0],
        coarse_range[1] + coarse_step * 0.5,
        coarse_step,
    ).tolist()

    print(
        f"  [resolution search] coarse candidates={len(coarse_list)} "
        f"[{coarse_range[0]:.1f} -> {coarse_range[1]:.1f}, step={coarse_step}], "
        f"target K={target_k}, n_neighbors={n_neighbors}"
    )

    coarse_records: list[dict[str, Any]] = []
    for res in coarse_list:
        actual_k = _eval_resolution(res)
        record = _record("coarse", res, actual_k)
        search_log.append(record)
        coarse_records.append(record)

        if verbose:
            print(
                f"    [coarse] resolution={record['resolution']:.4f}: "
                f"K={record['actual_k']}, diff={record['diff']}"
            )

    coarse_exact = _select_middle_exact_match(coarse_records)
    if coarse_exact is not None:
        n_exact = sum(record["actual_k"] == target_k for record in coarse_records)
        best_res = float(coarse_exact["resolution"])
        print(
            f"  [resolution search] coarse exact matches={n_exact}; "
            f"selected middle resolution={best_res:.4f}, K={target_k}"
        )
        return round(best_res, 4), target_k, search_log

    coarse_best = _best_by_diff(coarse_records)
    best_res = float(coarse_best["resolution"])
    best_k = int(coarse_best["actual_k"])
    best_diff = coarse_best["diff"]

    print(
        f"  [resolution search] coarse best: "
        f"resolution={best_res:.4f}, K={best_k}, diff={best_diff}"
    )

    final_diff = abs(best_k - target_k) if best_k != -1 else float("inf")
    print(
        f"  [resolution search] final selected: "
        f"resolution={best_res:.4f}, actual K={best_k}, target K={target_k}, diff={final_diff}"
    )

    return round(float(best_res), 4), best_k, search_log


def get_leiden_labels_auto_resolution(
    data_eval: np.ndarray,
    n_clusters: int,
    cluster_params: dict[str, Any] | None = None,
    iter_idx: int | None = None,
    verbose: bool = False,
) -> np.ndarray:
    """
    Run Leiden after searching for the resolution closest to n_clusters.
    """
    from jullei._01_04_01_04_seurat import csp_seurat_from_array

    cluster_params = cluster_params or {}
    n_neighbors = cluster_params.get("n_neighbors", 20)
    random_state = cluster_params.get("random_state", 42)

    if verbose and iter_idx is not None:
        print(f"  [resolution search] iteration={iter_idx}")

    best_res, _, _ = find_leiden_resolution_for_k(
        data_eval=data_eval,
        target_k=n_clusters,
        n_neighbors=n_neighbors,
        random_state=random_state,
        verbose=verbose,
    )

    return get_leiden_labels_fixed_resolution(
        data_eval=data_eval,
        resolution=best_res,
        n_neighbors=n_neighbors,
        random_state=random_state,
    )


def get_leiden_labels_fixed_resolution(
    data_eval: np.ndarray,
    resolution: float,
    n_neighbors: int = 20,
    random_state: int = 42,
) -> np.ndarray:
    """Run Leiden/Seurat with a fixed resolution."""
    from jullei._01_04_01_04_seurat import csp_seurat_from_array

    labels, _, _ = csp_seurat_from_array(
        data_eval,
        resolution=resolution,
        n_neighbors=n_neighbors,
        random_state=random_state,
    )
    return np.asarray(labels)


def get_leiden_labels_with_params(
    data_eval: np.ndarray,
    n_clusters: int,
    n_neighbors: int = 20,
    random_state: int = 42,
    verbose: bool = False,
) -> np.ndarray:
    """
    Convenience wrapper for call sites that pass n_neighbors directly.
    """
    return get_leiden_labels_auto_resolution(
        data_eval=data_eval,
        n_clusters=n_clusters,
        cluster_params={
            "n_neighbors": n_neighbors,
            "random_state": random_state,
        },
        verbose=verbose,
    )
