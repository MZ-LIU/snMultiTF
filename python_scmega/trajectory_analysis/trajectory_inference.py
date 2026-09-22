"""
Trajectory Inference Module

Implements trajectory inference algorithms strictly following the R version of scMEGA,
specifically the AddTrajectory function based on ArchR methodology.

This is a faithful Python implementation of the R AddTrajectory function:
```R
AddTrajectory(object, trajectory, group.by, reduction, dims, 
              pre.filter.quantile, post.filter.quantile, use.all, dof, spar, seed)
```

Key functions:
- add_trajectory: Main trajectory inference function (R AddTrajectory equivalent)
- get_quantiles: Quantile calculation utility (R getQuantiles equivalent)
- center_roll_mean: Rolling mean calculation (R centerRollMean equivalent)

References:
- Original R code: AddTrajectory function in scMEGA R package
- ArchR: Granja et al. (2021) ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis
- ArchR addTrajectory: https://www.archrproject.com/reference/addTrajectory.html
"""

import json
import subprocess
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
from typing import Optional, List, Union
from scipy import interpolate
from sklearn.neighbors import NearestNeighbors

# Try to import MuData, set to None if not available
try:
    import muon as mu
    MUON_AVAILABLE = True
except ImportError:
    mu = None
    MUON_AVAILABLE = False


def _is_mudata(data) -> bool:
    """Return True for MuData-like objects across muon/mudata version aliases."""
    if MUON_AVAILABLE and mu is not None:
        mudata_cls = getattr(mu, "MuData", None)
        if mudata_cls is not None and isinstance(data, mudata_cls):
            return True

    return all(hasattr(data, attr) for attr in ("mod", "obs", "obsm"))


def _default_dims_for_reduction(reduction: str) -> List[int]:
    if reduction == "X_umap":
        return [0, 1]
    if reduction in ["X_pca", "X_joint", "X_integrated"]:
        return list(range(30))
    return list(range(10))


def _sync_trajectory_to_modalities(mdata, name: str) -> None:
    for mod in ("rna", "atac"):
        if mod in mdata.mod and name in mdata.obs.columns:
            mdata[mod].obs[name] = mdata.obs.loc[mdata[mod].obs_names, name].to_numpy()


def _write_archr_trajectory_inputs(
    mdata,
    group_by: str,
    obs_source: pd.DataFrame,
    reduction: str,
    dims: List[int],
    output_dir: Path,
) -> tuple:
    if reduction not in mdata.obsm:
        raise ValueError(
            f"Reduction '{reduction}' not found in MuData. "
            f"Available reductions: {list(mdata.obsm.keys())}"
        )

    embedding = np.asarray(mdata.obsm[reduction])
    if embedding.ndim != 2:
        raise ValueError(f"mdata.obsm['{reduction}'] must be a 2D matrix.")
    if max(dims) >= embedding.shape[1] or min(dims) < 0:
        raise ValueError(
            f"Requested dims {dims} exceed reduction '{reduction}' "
            f"with shape {embedding.shape}."
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    metadata_path = output_dir / "archr_cell_metadata.csv"
    embedding_path = output_dir / "archr_embedding.csv"

    if group_by not in obs_source.columns:
        raise ValueError(f"group_by column '{group_by}' not found in observation metadata.")

    metadata = pd.DataFrame({
        "mdata_cell_id": mdata.obs_names.astype(str),
        group_by: obs_source.loc[mdata.obs_names, group_by].astype(str).to_numpy(),
    })
    metadata.to_csv(metadata_path, index=False)

    embedding_df = pd.DataFrame(
        embedding[:, dims],
        columns=[f"Dim{i + 1}" for i in range(len(dims))],
    )
    embedding_df.insert(0, "mdata_cell_id", mdata.obs_names.astype(str))
    embedding_df.to_csv(embedding_path, index=False)

    return metadata_path, embedding_path


def add_trajectory_archr(data: 'mu.MuData',
                         trajectory: List[str],
                         group_by: str,
                         assay: str = "RNA",
                         reduction: str = "X_umap",
                         dims: Optional[List[int]] = None,
                         pre_filter_quantile: float = 0.9,
                         post_filter_quantile: float = 0.9,
                         use_all: bool = False,
                         dof: int = 250,
                         spar: float = 1.0,
                         name: str = "Trajectory",
                         seed: int = 42,
                         archr_project_dir: Optional[Union[str, Path]] = None,
                         archr_output_dir: Optional[Union[str, Path]] = None,
                         archr_rscript_path: Optional[Union[str, Path]] = None,
                         archr_r_executable: str = "Rscript",
                         archr_embedding_name: str = "pymegaTrajectoryEmbedding",
                         archr_fragments_path: Optional[Union[str, Path]] = None,
                         archr_genome: str = "hg38",
                         archr_sample_name: str = "sample1",
                         archr_min_tss: float = 0,
                         archr_min_frags: int = 0,
                         archr_add_tile_mat: bool = False,
                         archr_add_gene_score_mat: bool = False,
                         archr_threads: int = 1,
                         archr_force: bool = True,
                         archr_save_project: bool = False) -> 'mu.MuData':
    """Run ArchR::addTrajectory() via Rscript and write pseudotime back to MuData."""
    if not _is_mudata(data):
        raise ValueError(
            "ArchR trajectory requires a MuData object with mod/obs/obsm. "
            f"Got {type(data).__module__}.{type(data).__name__}."
        )
    if archr_project_dir is None and archr_fragments_path is None:
        raise ValueError(
            "backend='archr' requires either archr_project_dir or archr_fragments_path. "
            "Current workflow can pass FRAGMENTS_PATH to build a temporary ArchRProject."
        )
    if not trajectory:
        raise ValueError("Please specify the trajectory.")
    if not group_by:
        raise ValueError("Please specify group_by.")

    mdata = data
    if assay == "RNA" and 'rna' in mdata.mod:
        adata = mdata['rna']
    elif assay == "ATAC" and 'atac' in mdata.mod:
        adata = mdata['atac']
    else:
        raise ValueError(f"Assay '{assay}' not found in MuData modalities: {list(mdata.mod.keys())}")

    if group_by in mdata.obs.columns:
        obs_source = mdata.obs
    elif group_by in adata.obs.columns:
        obs_source = adata.obs
        mdata.obs[group_by] = adata.obs.loc[mdata.obs_names, group_by].to_numpy()
    else:
        raise ValueError(f"group_by column '{group_by}' not found in mdata.obs or {assay}.obs.")

    if dims is None:
        dims = _default_dims_for_reduction(reduction)

    if archr_output_dir is None:
        output_dir = Path(tempfile.mkdtemp(prefix="pymega_archr_trajectory_"))
    else:
        output_dir = Path(archr_output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    if archr_project_dir is not None:
        archr_project_dir = Path(archr_project_dir)
        if not archr_project_dir.exists():
            raise FileNotFoundError(f"ArchRProject directory not found: {archr_project_dir}")
    if archr_fragments_path is not None:
        archr_fragments_path = Path(archr_fragments_path)
        if not archr_fragments_path.exists():
            raise FileNotFoundError(f"fragments file not found: {archr_fragments_path}")

    if archr_rscript_path is None:
        archr_rscript_path = Path(__file__).resolve().parent / "r_scripts" / "run_archr_trajectory.R"
    else:
        archr_rscript_path = Path(archr_rscript_path)
    if not archr_rscript_path.exists():
        raise FileNotFoundError(f"ArchR trajectory R script not found: {archr_rscript_path}")

    metadata_path, embedding_path = _write_archr_trajectory_inputs(
        mdata=mdata,
        group_by=group_by,
        obs_source=obs_source,
        reduction=reduction,
        dims=dims,
        output_dir=output_dir,
    )

    output_csv = output_dir / "archr_trajectory_values.csv"
    config = {
        "archr_project_dir": str(archr_project_dir) if archr_project_dir is not None else None,
        "fragments_path": str(archr_fragments_path) if archr_fragments_path is not None else None,
        "genome": archr_genome,
        "sample_name": archr_sample_name,
        "min_tss": archr_min_tss,
        "min_frags": int(archr_min_frags),
        "add_tile_mat": bool(archr_add_tile_mat),
        "add_gene_score_mat": bool(archr_add_gene_score_mat),
        "metadata_csv": str(metadata_path),
        "embedding_csv": str(embedding_path),
        "output_csv": str(output_csv),
        "trajectory_name": name,
        "trajectory": [str(x) for x in trajectory],
        "group_by": group_by,
        "embedding_name": archr_embedding_name,
        "pre_filter_quantile": pre_filter_quantile,
        "post_filter_quantile": post_filter_quantile,
        "use_all": bool(use_all),
        "dof": int(dof),
        "spar": float(spar),
        "seed": int(seed),
        "threads": int(archr_threads),
        "force": bool(archr_force),
        "save_project": bool(archr_save_project),
    }
    config_path = output_dir / "archr_trajectory_config.json"
    config_path.write_text(json.dumps(config, indent=2), encoding="utf-8")

    cmd = [archr_r_executable, str(archr_rscript_path), str(config_path)]
    print("Running ArchR trajectory via Rscript...")
    print(f"  ArchRProject: {archr_project_dir if archr_project_dir is not None else '[temporary from fragments]'}")
    print(f"  Output dir: {output_dir}")
    completed = subprocess.run(
        cmd,
        cwd=str(output_dir),
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "ArchR trajectory Rscript failed.\n"
            f"Command: {' '.join(cmd)}\n"
            f"STDOUT:\n{completed.stdout}\n"
            f"STDERR:\n{completed.stderr}"
        )
    if not output_csv.exists():
        raise FileNotFoundError(f"ArchR trajectory output was not created: {output_csv}")

    trajectory_df = pd.read_csv(output_csv)
    required_cols = {"mdata_cell_id", name}
    missing_cols = required_cols.difference(trajectory_df.columns)
    if missing_cols:
        raise ValueError(f"ArchR output missing required columns: {sorted(missing_cols)}")

    trajectory_series = pd.Series(np.nan, index=mdata.obs_names, dtype="float64")
    values = pd.to_numeric(trajectory_df[name], errors="coerce")
    mapped = pd.Series(values.to_numpy(), index=trajectory_df["mdata_cell_id"].astype(str))
    overlapping = trajectory_series.index.intersection(mapped.index)
    trajectory_series.loc[overlapping] = mapped.loc[overlapping].to_numpy()

    mdata.obs[name] = trajectory_series
    _sync_trajectory_to_modalities(mdata, name)

    valid_values = trajectory_series.dropna()
    if len(valid_values) == 0:
        raise ValueError("ArchR returned no non-NA trajectory values for MuData cells.")

    mdata.uns.setdefault("trajectory_analysis", {})
    mdata.uns["trajectory_analysis"][name] = {
        "method": "ArchR::addTrajectory",
        "trajectory_groups": list(trajectory),
        "group_by": group_by,
        "reduction": reduction,
        "dims": dims,
        "pre_filter_quantile": pre_filter_quantile,
        "post_filter_quantile": post_filter_quantile,
        "dof": dof,
        "spar": spar,
        "use_all": use_all,
        "seed": seed,
        "archr_project_dir": str(archr_project_dir) if archr_project_dir is not None else None,
        "archr_fragments_path": str(archr_fragments_path) if archr_fragments_path is not None else None,
        "archr_genome": archr_genome,
        "archr_output_dir": str(output_dir),
        "archr_embedding_name": archr_embedding_name,
        "n_trajectory_cells": int(len(valid_values)),
        "trajectory_range": [float(valid_values.min()), float(valid_values.max())],
    }

    print(f" ArchR trajectory '{name}' added successfully")
    print(f"Trajectory cells: {len(valid_values)} / {mdata.n_obs}")
    print(f"Trajectory range: {valid_values.min():.1f} - {valid_values.max():.1f}")
    return mdata


def add_trajectory(data: 'mu.MuData',
                   trajectory: List[str],
                   group_by: str,
                   assay: str = "RNA",
                   reduction: str = "X_umap",  # Change to use UMAP by default (integrated)
                   dims: Optional[List[int]] = None,
                   pre_filter_quantile: float = 0.9,
                   post_filter_quantile: float = 0.9,
                   use_all: bool = False,
                   dof: int = 250,
                   spar: float = 1.0,
                   name: str = "Trajectory",
                   seed: int = 42,
                   backend: str = "python",
                   archr_project_dir: Optional[Union[str, Path]] = None,
                   archr_output_dir: Optional[Union[str, Path]] = None,
                   archr_rscript_path: Optional[Union[str, Path]] = None,
                   archr_r_executable: str = "Rscript",
                   archr_embedding_name: str = "pymegaTrajectoryEmbedding",
                   archr_fragments_path: Optional[Union[str, Path]] = None,
                   archr_genome: str = "hg38",
                   archr_sample_name: str = "sample1",
                   archr_min_tss: float = 0,
                   archr_min_frags: int = 0,
                   archr_add_tile_mat: bool = False,
                   archr_add_gene_score_mat: bool = False,
                   archr_threads: int = 1,
                   archr_force: bool = True,
                   archr_save_project: bool = False) -> 'mu.MuData':
    """
    Add trajectory information to single-cell data using ArchR methodology.
    
    This function is a faithful Python implementation of the R AddTrajectory function:
    ```R
    obj <- AddTrajectory(
        object = obj,
        trajectory = c("group1", "group2", "group3"),
        group.by = "leiden",
        reduction = "MOJITOO_UMAP", # For multi-omics data, use integrated UMAP
        dims = 1:2,                  # R uses 1-based indexing
        pre.filter.quantile = 0.9,
        post.filter.quantile = 0.9,
        use.all = FALSE,
        dof = 250,
        spar = 1,
        seed = 42
    )
    ```
    
    **IMPORTANT - Dimension Indexing Difference**:
    - **R version**: Uses 1-based indexing, e.g., `dims = 1:3` means dimensions 1, 2, 3
    - **Python version**: Uses 0-based indexing, e.g., `dims = [0, 1, 2]` means dimensions 0, 1, 2
    - To replicate R's `dims = 1:30`, use Python's `dims = list(range(30))` (0-29)
    - To replicate R's `dims = 1:3`, use Python's `dims = [0, 1, 2]`
    
    Args:
        data: MultiomeData or AnnData object
        trajectory: The order of cell groups to be used for constraining the initial 
                   supervised fitting procedure (corresponds to R 'trajectory' parameter)
        group_by: The column name in metadata that contains the cell group definitions 
                 used in trajectory (corresponds to R 'group.by' parameter)
        assay: Which assay to use ("RNA" or "ATAC")
        reduction: Name of dimension reduction used to infer the trajectory.
          **For multiome data, ONLY integrated reductions are supported**: 
          - "X_umap": integrated UMAP (REQUIRED, matches R's MOJITOO_UMAP)
          - "X_integrated": integrated embedding from CCA
          - "X_joint": joint embedding (alternative name)
          Single-modality reductions (X_pca, X_lsi) are NOT supported for MuData.
        dims: **0-BASED** list of dimension indices to use (corresponds to R 'dims' parameter).
              **R→Python conversion**: R's `dims = 1:n` → Python's `dims = list(range(n))`
              Default: [0, 1] for UMAP (2D), [0-29] for PCA/joint embeddings
        pre_filter_quantile: Prior to the initial supervised trajectory fitting,
                           cells whose euclidean distance from the cell-grouping center is
                           above the provided quantile will be excluded
        post_filter_quantile: After initial supervised trajectory fitting, cells
                            whose euclidean distance from the cell-grouping center is above 
                            the provided quantile will be excluded
        use_all: Whether or not to use cells outside of trajectory groups for post-fitting procedure
        dof: The number of degrees of freedom to be used in the spline fit
        spar: The sparsity to be used in the spline fit
        name: A string indicating the name of the fitted trajectory to be added in metadata
        seed: Random seed
        
    Returns:
        Data object with trajectory information added (pseudotime in 0-100 range)
        
    Example:
        ```python
        # Recommended method: use multi-omics integrated UMAP (consistent with the R version MOJITOO_UMAP)
        # R version: AddTrajectory(obj, reduction = "MOJITOO_UMAP", dims = 1:2)
        # Python equivalent:
        multiome = add_trajectory(
            multiome,
            trajectory=["CD4 Naive", "CD4 TCM"],
            group_by="predicted.id",
            reduction="umap", # Use integrated UMAP
            dims=[0, 1]           # 2D UMAP (R's 1:2 → Python's [0,1])
        )
        
        # Or use integrated embedding (high-dimensional)
        # R version: AddTrajectory(obj, reduction = "MOJITOO", dims = 1:30)
        # Python equivalent:
        multiome = add_trajectory(
            multiome,
            trajectory=["CD4 Naive", "CD4 TCM"],
            group_by="predicted.id",
            reduction="joint", # or "integrated"
            dims=list(range(30))      # 0-29 in Python = 1-30 in R
        )
        
        # Single modal analysis (not recommended for multi-omics data)
        multiome = add_trajectory(
            multiome,
            trajectory=["0", "1", "2"],
            group_by="leiden",
            reduction="pca",
            dims=list(range(30))
        )
        ```
    """
    if backend.lower() == "archr":
        return add_trajectory_archr(
            data=data,
            trajectory=trajectory,
            group_by=group_by,
            assay=assay,
            reduction=reduction,
            dims=dims,
            pre_filter_quantile=pre_filter_quantile,
            post_filter_quantile=post_filter_quantile,
            use_all=use_all,
            dof=dof,
            spar=spar,
            name=name,
            seed=seed,
            archr_project_dir=archr_project_dir,
            archr_output_dir=archr_output_dir,
            archr_rscript_path=archr_rscript_path,
            archr_r_executable=archr_r_executable,
            archr_embedding_name=archr_embedding_name,
            archr_fragments_path=archr_fragments_path,
            archr_genome=archr_genome,
            archr_sample_name=archr_sample_name,
            archr_min_tss=archr_min_tss,
            archr_min_frags=archr_min_frags,
            archr_add_tile_mat=archr_add_tile_mat,
            archr_add_gene_score_mat=archr_add_gene_score_mat,
            archr_threads=archr_threads,
            archr_force=archr_force,
            archr_save_project=archr_save_project,
        )
    if backend.lower() != "python":
        raise ValueError("backend must be 'python' or 'archr'.")

    print(f"Inferring trajectory using ArchR methodology...")
    print(f"Trajectory path: {' → '.join(trajectory)}")
    
    # Set random seed
    np.random.seed(seed)
    
    if _is_mudata(data):
        is_multiome = True
        mdata = data
        if assay == "RNA" and 'rna' in data.mod:
            adata = data['rna']  # Only used to get obs
        elif assay == "ATAC" and 'atac' in data.mod:
            adata = data['atac']  # Only used to get obs
        else:
            raise ValueError(f"Assay '{assay}' not found in MuData modalities: {list(data.mod.keys())}")
        
        # Search logic for dimensionality reduction data: Prioritize the MuData layer (integrated dimensionality reduction), followed by the sub-object layer
        embedding_key = f'{reduction}'
        
        if embedding_key in mdata.obsm:
            use_mdata_layer = True
            print(f"  Integrated dimensionality reduction using MuData layer: {embedding_key}")
        else:
            raise ValueError(
                f"Reduction '{embedding_key}' not found in MuData. Please compute it first.\n"
                f"Available reductions: {list(mdata.obsm.keys())}"
            )
            
    else:
        raise ValueError("Please provide a MuData object!")
        

    # Validate parameters
    if not trajectory:
        raise ValueError("Please specify the trajectory!")

    if not group_by:
        raise ValueError("Please specify the group based on which the trajectory is inferred!")
    if is_multiome and use_mdata_layer:
        # When using integrated data, the grouping information is read from the global layer first, and if it does not exist, it falls back from the modal layer.
        if group_by in mdata.obs.columns:
            obs_source = mdata.obs  # Use global layer
            print(f" Use the grouping information of the global layer: {group_by}")
        elif group_by in adata.obs.columns:
            obs_source = adata.obs  # Fall back to modal layer
            print(f" The group information '{group_by}' is in the modal layer. It is recommended to synchronize it to the global layer.")
            print(f" Tip: mdata.obs['{group_by}'] = mdata['{assay}'].obs['{group_by}'].copy()")
        else:
            raise ValueError(
                f"Grouping column '{group_by}' not found! \n"
                f"MuData global layer: {list(mdata.obs.columns)}\n"
                f"{assay} Modal layer: {list(adata.obs.columns)}"
            )

    data_source = mdata  # Data sources always use the global layer
    
    if not reduction:
        raise ValueError("Please provide the dimensional reduction!")

    # Set default dimensions based on reduction type
    if dims is None:
        if reduction == "X_umap":
            dims = [0, 1]  # UMAP is typically 2D (corresponding to R dims = 1:2)
        elif reduction in ["X_pca", "X_joint", "X_integrated"]:
            dims = list(range(30))  # High-dimensional embeddings (corresponds to R dims = 1:30)
        else:
            dims = list(range(10))  # Conservative default for unknown reductions

    
    if len(dims) == 1:
        data_all = mdata.obsm[embedding_key][:, dims].reshape(-1, 1)
    else:
        data_all = mdata.obsm[embedding_key][:, dims]

    print(f"Using {data_all.shape[1]} dimensions from {reduction}")
    
    # Get group information
    df_group = obs_source[[group_by]].copy()
    df_group.columns = ['group']
    
    # Filter to trajectory groups only
    trajectory_mask = df_group['group'].astype(str).isin([str(t) for t in trajectory])
    df_group_filtered = df_group[trajectory_mask].copy()
    
    if len(df_group_filtered) == 0:
        raise ValueError("Cannot find the specified trajectory groups!")
    
    unique_groups = df_group_filtered['group'].astype(str).unique()
    missing_groups = set([str(t) for t in trajectory]) - set(unique_groups)
    if missing_groups:
        print(f"Warning: Trajectory groups not found in data: {missing_groups}")
    
    print(f"Trajectory groups found: {list(unique_groups)}")
    print(f"Total trajectory cells: {len(df_group_filtered)}")
    
    # Extract data for trajectory groups only
    # R equivalent: data.use <- data.use[rownames(df.group), , drop = FALSE]
    data_use = data_all[trajectory_mask]
    trajectory_cell_indices = np.where(trajectory_mask)[0]  # Store original indices
    
    print(f"Trajectory cells: {data_use.shape[0]} / {data_all.shape[0]}")
    
    ######################################################
    # Filter Outliers (Pre-filtering)
    # R equivalent: filterObj <- lapply(...) in lines 95-110
    ######################################################
    print("Filtering outliers (pre-filtering)...")
    
    filter_objects = []
    
    for group in trajectory:
        group_str = str(group)
        if group_str not in unique_groups:
            continue
            
        # Get cells in this group (within trajectory cells)
        # R: groupsx <- rownames(df.group)[df.group[, 1] == trajectory[x]]
        group_mask_in_trajectory = df_group_filtered['group'].astype(str) == group_str
        
        # Get data for this group from data_use
        # R: matx <- data.use[groupsx, , drop = FALSE]
        group_data = data_use[group_mask_in_trajectory.values]
        
        if len(group_data) == 0:
            continue
        
        # Calculate distances from group center
        # R: matMeanx <- colMeans(matx)
        # R: diffx <- sqrt(colSums((t(matx) - matMeanx) ^ 2))
        group_center = np.mean(group_data, axis=0)
        distances = np.sqrt(np.sum((group_data - group_center) ** 2, axis=1))
        
        # Filter by quantile
        # R: idxKeep <- which(diffx <= quantile(diffx, pre.filter.quantile))
        distance_threshold = np.quantile(distances, pre_filter_quantile)
        keep_mask = distances <= distance_threshold
        
        filtered_data = group_data[keep_mask]
        # Store indices within data_use (trajectory cells)
        filtered_indices_in_trajectory = np.where(group_mask_in_trajectory.values)[0][keep_mask]
        
        filter_objects.append({
            'data': filtered_data,
            'indices': filtered_indices_in_trajectory,  # Indices within data_use
            'group': group_str
        })
        
        print(f"  Group {group}: {len(group_data)} → {len(filtered_data)} cells")
    
    if len(filter_objects) == 0:
        raise ValueError("No cells remaining after pre-filtering")
    
    # Combine filtered data
    # R: matFilter <- lapply(...) %>% Reduce("rbind", .)
    mat_filter = np.vstack([obj['data'] for obj in filter_objects])
    filter_indices_in_trajectory = np.concatenate([obj['indices'] for obj in filter_objects])
    groups_filter = []
    
    for obj in filter_objects:
        groups_filter.extend([obj['group']] * len(obj['data']))
    
    groups_filter = np.array(groups_filter)
    
    print(f"After pre-filtering: {len(mat_filter)} cells")
    
    ######################################################
    ######################################################
    print("Computing initial time alignment...")
    
    initial_times = []
    cell_times = {}
    
    for i, group in enumerate(trajectory):
        group_str = str(group)
        group_mask = groups_filter == group_str
        
        if not np.any(group_mask):
            continue
            
        group_data = mat_filter[group_mask]
        group_indices_in_trajectory = filter_indices_in_trajectory[group_mask]
        
        if i < len(trajectory) - 1:
            # Not the last group - compute distance to next group
            next_group = str(trajectory[i + 1])
            next_group_mask = groups_filter == next_group
            
            if np.any(next_group_mask):
                next_group_data = mat_filter[next_group_mask]
                next_group_center = np.mean(next_group_data, axis=0)
                
                # Calculate distances to next group center
                distances = np.sqrt(np.sum((group_data - next_group_center) ** 2, axis=1))
                # Convert to quantiles and flip (1 - quantiles) then add group index
                time_values = (1 - get_quantiles(distances)) + i
            else:
                time_values = np.full(len(group_data), i)
        else:
            # Last group - compute distance to previous group
            prev_group = str(trajectory[i - 1])
            prev_group_mask = groups_filter == prev_group
            
            if np.any(prev_group_mask):
                prev_group_data = mat_filter[prev_group_mask]
                prev_group_center = np.mean(prev_group_data, axis=0)
                
                # Calculate distances to previous group center
                distances = np.sqrt(np.sum((group_data - prev_group_center) ** 2, axis=1))
                # Convert to quantiles and add group index
                time_values = get_quantiles(distances) + i
            else:
                time_values = np.full(len(group_data), i)
        
        # Store times for each cell (using indices within trajectory cells)
        for j, cell_idx_in_traj in enumerate(group_indices_in_trajectory):
            cell_times[cell_idx_in_traj] = time_values[j]
            initial_times.append(time_values[j])
    
    initial_times = np.array(initial_times)
    
    print(f"Initial time range: {initial_times.min():.3f} - {initial_times.max():.3f}")
    
    ######################################################
    ######################################################
    print("Fitting cubic splines...")
    
    # Sort by initial time
    time_order = np.argsort(initial_times)
    sorted_times = initial_times[time_order]
    sorted_data = mat_filter[time_order]
    
    # Fit spline for each dimension
    spline_data = np.zeros_like(sorted_data)
    
    for dim in range(sorted_data.shape[1]):
        y_values = sorted_data[:, dim]
        
        try:
            # Use UnivariateSpline which is closest to R's smooth.spline
            # Convert dof to smoothing parameter s
            if len(np.unique(sorted_times)) > 3:  # Need at least 4 unique points
                spline = interpolate.UnivariateSpline(
                    sorted_times, y_values, 
                    s=len(sorted_times) / dof * spar,  # Approximate conversion
                    k=3  # Cubic spline
                )
                spline_data[:, dim] = spline(sorted_times)
            else:
                # Fallback to linear interpolation
                spline_data[:, dim] = np.interp(sorted_times, sorted_times, y_values)
                
        except Exception as e:
            print(f"Spline fitting failed for dimension {dim}: {e}")
            spline_data[:, dim] = y_values
    
    print("Cubic spline fitting completed")
    
    ######################################################
    # R equivalent: knnObj <- nabor::knn(data = matSpline, query = data.use, k = 3)
    ######################################################
    print("Performing KNN fitting (k=3)...")
    
    # Build KNN with k=3 (matching R version exactly)
    nbrs = NearestNeighbors(n_neighbors=3, algorithm='auto')
    nbrs.fit(spline_data)
    
    # Query against trajectory cells only (data_use)
    # R: query = data.use (which is already filtered to trajectory groups)
    distances, indices = nbrs.kneighbors(data_use)
    
    # Extract KNN results (for trajectory cells only)
    knn_idx = indices  # Shape: (n_trajectory_cells, 3)
    knn_dist = distances  # Shape: (n_trajectory_cells, 3)
    
    # Calculate direction correction (matching R logic exactly)
    # R: knnDiff <- ifelse(knnIdx[, 2] > knnIdx[, 3], 1, -1)
    knn_diff = np.where(knn_idx[:, 1] > knn_idx[:, 2], 1, -1)
    
    # Convert distances to quantiles
    # R: knnDistQ <- getQuantiles(knnDist[, 1])
    knn_dist_q = get_quantiles(knn_dist[:, 0])  # Use distance to 1st neighbor
    
    # Calculate trajectory position index
    # R: DistanceIdx = knnIdx[, 1] + knnDiff * knnDistQ
    distance_idx = knn_idx[:, 0] + knn_diff * knn_dist_q
    
    ######################################################
    # Post-filtering and Final Trajectory Assignment
    # R equivalent: lines 191-210
    ######################################################
    print("Applying post-filtering...")
    
    # Filter outlier cells based on distance to trajectory (within trajectory cells)
    # R: idxKeep <- which(knnDist[, 1] <= quantile(knnDist[, 1], post.filter.quantile))
    distance_threshold = np.quantile(knn_dist[:, 0], post_filter_quantile)
    keep_mask_in_trajectory = knn_dist[:, 0] <= distance_threshold
    
    print(f"Post-filtering: {np.sum(keep_mask_in_trajectory)} / {len(keep_mask_in_trajectory)} trajectory cells kept")
    
    # Convert distance index to 0-100 trajectory values (matching R exactly)
    # R: nas <- rep(NA, dim(object)[2])
    trajectory_values = np.full(len(adata), np.nan)
    
    if np.sum(keep_mask_in_trajectory) > 0:
        # R: dfTrajectory3$Trajectory <- 100 * getQuantiles(dfTrajectory3[, 2])
        kept_distance_idx = distance_idx[keep_mask_in_trajectory]
        trajectory_quantiles = get_quantiles(kept_distance_idx)
        trajectory_scaled = trajectory_quantiles * 100
        
        # Map back to original cell indices
        # R: nas[rownames(dfTrajectory3)] <- dfTrajectory3$Trajectory
        kept_cell_indices = trajectory_cell_indices[keep_mask_in_trajectory]
        trajectory_values[kept_cell_indices] = trajectory_scaled
    
    # Add trajectory to observations
    data_source.obs[name] = trajectory_values
    
    # Store trajectory metadata
    if 'trajectory_analysis' not in data_source.uns:
        data_source.uns['trajectory_analysis'] = {}
    
    data_source.uns['trajectory_analysis'][name] = {
        'method': 'yangtiaopinghua',
        'trajectory_groups': trajectory,
        'group_by': group_by,
        'reduction': reduction,
        'dims': dims,
        'pre_filter_quantile': pre_filter_quantile,
        'post_filter_quantile': post_filter_quantile,
        'dof': dof,
        'spar': spar,
        'use_all': use_all,
        'seed': seed,
        'n_trajectory_cells': np.sum(~np.isnan(trajectory_values)),
        'trajectory_range': [np.nanmin(trajectory_values), np.nanmax(trajectory_values)]
    }
    
    print(f" Trajectory '{name}' added successfully")
    print(f"Trajectory cells: {np.sum(~np.isnan(trajectory_values))} / {len(data_source)}")
    
    if np.sum(~np.isnan(trajectory_values)) > 0:
        valid_values = trajectory_values[~np.isnan(trajectory_values)]
        print(f"Trajectory range: {valid_values.min():.1f} - {valid_values.max():.1f}")
    
    # Handle MultiomeData output
    if is_multiome:
        # If integrated data is used, the trajectory information is already in the global layer and needs to be synchronized to the modal layer.
        if use_mdata_layer:
            if 'rna' in mdata.mod and 'atac' in mdata.mod:
                shared_cells = list(set(mdata['rna'].obs.index) & set(mdata['atac'].obs.index))
                if shared_cells and name in mdata.obs.columns:
                    _sync_trajectory_to_modalities(mdata, name)
                    print(f" Synchronized '{name}' to both RNA and ATAC modalities ({len(shared_cells)} cells)")
            
    return mdata  # Return MuData object
   

def get_quantiles(v: np.ndarray, length: Optional[int] = None) -> np.ndarray:
    """
    Calculate quantiles for a vector, matching R getQuantiles function exactly.
    
    This is a faithful Python implementation of the R getQuantiles function:
    ```R
    getQuantiles <- function(v = NULL, len = length(v)){
        p <- trunc(rank(v2))/length(v2)
        return(p)
    }
    ```
    
    Args:
        v: Input vector
        length: Target length (defaults to length of v)
        
    Returns:
        Quantile values
    """
    if length is None:
        length = len(v)
    
    if len(v) < length:
        v2 = np.zeros(length)
        v2[:len(v)] = v
    else:
        v2 = v.copy()
    
    # Calculate quantiles using truncated ranks (matching R exactly)
    # R's rank() uses ties.method="average" by default, then trunc()
    from scipy.stats import rankdata
    ranks = rankdata(v2, method='average')  # average rank for ties, matching R default
    p = np.trunc(ranks) / len(v2)  # Truncate ranks and normalize
    
    if len(v) < length:
        p = p[:len(v)]
    
    return p


def center_roll_mean(v: np.ndarray, k: int) -> np.ndarray:
    """
    Calculate center-aligned rolling mean, matching R centerRollMean function.
    
    This is a faithful Python implementation of the R centerRollMean function:
    ```R
    centerRollMean <- function(v = NULL, k = NULL){
        o1 <- data.table::frollmean(v, k, align = "right", na.rm = FALSE)
        # ... (padding logic)
    }
    ```
    
    Args:
        v: Input vector
        k: Window size
        
    Returns:
        Center-aligned rolling mean
    """
    if k <= 1:
        return v.copy()
    
    # Calculate right-aligned rolling mean
    rolled = pd.Series(v).rolling(window=k, min_periods=k).mean().values
    
    # Handle padding to center-align
    if k % 2 == 0:
        # Even window
        pad_left = k // 2 - 1
        pad_right = k // 2
    else:
        # Odd window
        pad_left = k // 2
        pad_right = k // 2
    
    # Remove initial NaN values and pad
    valid_start = k - 1
    valid_rolled = rolled[valid_start:]
    
    if len(valid_rolled) == 0:
        return v.copy()
    
    # Create output with padding
    result = np.zeros_like(v)
    
    # Left padding
    if pad_left > 0:
        result[:pad_left] = valid_rolled[0]
    
    # Main values
    end_main = min(len(result) - pad_right, len(valid_rolled))
    result[pad_left:pad_left + end_main] = valid_rolled[:end_main]
    
    # Right padding
    if pad_right > 0:
        result[-pad_right:] = valid_rolled[-1]
    
    return result
