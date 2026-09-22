"""
Gene Regulatory Network (GRN) Inference Module

Implements core GRN inference algorithms equivalent to scMEGA's GetGRN and 
GetTFGeneCorrelation functions, providing comprehensive regulatory network 
reconstruction from multiome data.

Key functions:
- infer_grn: Main GRN inference function (GetGRN equivalent)
- get_tf_gene_correlation: TF-gene correlation analysis (GetTFGeneCorrelation equivalent)
- calculate_tf_activity: Calculate TF activity from motif and expression data
- integrate_regulatory_evidence: Integrate multiple sources of regulatory evidence
- rank_regulatory_interactions: Rank and filter regulatory interactions

References:
- Original R code: GetGRN and GetTFGeneCorrelation functions in scMEGA
- Aibar et al. (2017) SCENIC: single-cell regulatory network inference and clustering
- Pliner et al. (2018) Cicero predicts cis-regulatory DNA interactions
- Granja et al. (2021) ArchR is a scalable software package for integrative analysis
"""

import numpy as np
import pandas as pd
import anndata as ad
from typing import Optional, Dict, Any, List, Union, Tuple
from pathlib import Path
import json
import warnings


# Try to import MuData, set to None if not available
try:
    import muon as mu
    MUON_AVAILABLE = True
except ImportError:
    mu = None
    MUON_AVAILABLE = False


def _ensure_dense_array(x: Any) -> np.ndarray:
    """Convert sparse/dense matrix-like objects to a 2D numpy array."""
    try:
        from scipy import sparse
        if sparse.issparse(x):
            return x.toarray()
    except Exception:
        pass
    return np.asarray(x)


def _aggregate_activity_by_tf_names(
    activity_matrix: Any,
    tf_names: List[str],
) -> Tuple[np.ndarray, List[str], pd.Series]:
    """Average duplicated motif activity columns that map to the same TF name."""
    tf_index = pd.Index([str(name).strip() for name in tf_names])
    if len(tf_index) == 0:
        raise ValueError("TF activity aggregation requires at least one TF name.")
    if tf_index.isna().any() or (tf_index == "").any():
        raise ValueError("TF activity aggregation found empty TF names.")

    activity = _ensure_dense_array(activity_matrix)
    if activity.ndim != 2:
        raise ValueError(f"TF activity matrix must be 2D, got shape {activity.shape}.")
    if activity.shape[1] != len(tf_index):
        raise ValueError(
            "TF activity matrix column count does not match TF names: "
            f"{activity.shape[1]} columns vs {len(tf_index)} names."
        )

    activity_df = pd.DataFrame(activity, columns=tf_index)
    aggregated = activity_df.T.groupby(level=0, sort=False).mean().T
    source_counts = pd.Series(tf_index).value_counts(sort=False).reindex(aggregated.columns)
    return aggregated.to_numpy(), aggregated.columns.astype(str).tolist(), source_counts


def _aggregate_chromvar_activity_by_tf(chromvar_adata: ad.AnnData) -> ad.AnnData:
    """Return chromVAR activity with duplicated motifs collapsed by TF gene name."""
    if "motif_name" in chromvar_adata.var.columns:
        tf_names = chromvar_adata.var["motif_name"].astype(str).tolist()
    else:
        tf_names = chromvar_adata.var_names.astype(str).tolist()

    aggregated_x, unique_tfs, source_counts = _aggregate_activity_by_tf_names(
        chromvar_adata.X,
        tf_names,
    )

    var = pd.DataFrame(index=pd.Index(unique_tfs, name=chromvar_adata.var_names.name))
    var["tf_name"] = unique_tfs
    var["n_motifs_aggregated"] = source_counts.to_numpy(dtype=int)

    chromvar_out = ad.AnnData(
        X=aggregated_x,
        obs=chromvar_adata.obs.copy(),
        var=var,
    )
    chromvar_out.uns["tf_activity_aggregation"] = {
        "method": "mean_by_tf_name",
        "input_motifs": int(chromvar_adata.n_vars),
        "output_tfs": int(chromvar_out.n_vars),
        "duplicated_tf_count": int((source_counts > 1).sum()),
    }
    return chromvar_out


def _apply_motif_names_from_atac(chromvar_adata: ad.AnnData, data: Any) -> ad.AnnData:
    """Attach ATAC motif names to chromVAR features when available."""
    if not (hasattr(data, "mod") and "atac" in data.mod):
        return chromvar_adata

    motif_info = data["atac"].uns.get("motifs", {})
    motif_names = motif_info.get("motif_names", None)
    if motif_names is None or len(motif_names) != chromvar_adata.n_vars:
        return chromvar_adata

    chromvar_adata = chromvar_adata.copy()
    chromvar_adata.var["motif_name"] = [str(name).strip() for name in motif_names]
    return chromvar_adata


def infer_grn(motif_matching: Union[np.ndarray, pd.DataFrame],
             tf_gene_cor: pd.DataFrame,
             peak_gene_links: pd.DataFrame,
             correlation_threshold: Optional[float] = None,  # Make optional
             return_dict: bool = False,  # New: whether to return dictionary format
             **kwargs) -> Dict[str, Any]:
    """
    Infer gene regulatory network - R-compatible (GetGRN) behavior:
      - Left join TF-gene correlation only if TF-peak-gene path exists
      - Filter only by TF-gene correlation threshold
    
    This function replicates scMEGA's GetGRN function (R/grn.R:130-172):
    R logic:
    1. Build TF-peak relationships (df.p2m) from motif_matching
    2. Join with peak-gene links to get TF-peak-gene paths (df.m2g)
    3. LEFT JOIN df.m2g with tf_gene_cor (correlation is optional)
    4. Return unfiltered results (filtering done in workflow)
    
    **KEY PRINCIPLE**: Only TF-gene pairs with peak-gene linkage support are retained,
    even if they have high correlation but no peak support.
    
    Args:
        motif_matching: Peak x TF motif matching matrix
        tf_gene_cor: TF-gene correlation DataFrame
        peak_gene_links: Peak-gene linkage DataFrame
        correlation_threshold: Optional correlation threshold for filtering.
                              If None (default), returns all pairs (R behavior).
                              If specified, filters after construction.
        return_dict: If True, returns dict with network graph. 
                    If False (default), returns DataFrame only (R behavior).
        **kwargs: Additional parameters
        
    Returns:
        If return_dict=False (default): pd.DataFrame (R-compatible)
        If return_dict=True: Dict with 'grn_network' and 'grn_edges'
    """
    print("=" * 70)
    print("Building GRN using R-compatible GetGRN logic")
    print("=" * 70)
    
    # ========================================================================
    # Step 0: Input validation and format conversion
    # ========================================================================
    # --- Input validation (match R behavior: all three inputs required) ---
    if motif_matching is None:
        raise ValueError("Please provide a motif matching matrix (peak × TF).")
    if tf_gene_cor is None or tf_gene_cor.empty:
        raise ValueError("Please provide a non-empty TF-gene correlation table (df.cor).")
    if peak_gene_links is None or peak_gene_links.empty:
        raise ValueError("Please provide peak-to-gene links (df.p2g).")
    
    if isinstance(motif_matching, np.ndarray):
        motif_df = pd.DataFrame(motif_matching)
        motif_df.index = motif_df.index.astype(str)
        motif_df.columns = motif_df.columns.astype(str)
    elif isinstance(motif_matching, pd.DataFrame):
        motif_df = motif_matching.copy()
        motif_df.index = motif_df.index.astype(str)
        motif_df.columns = motif_df.columns.astype(str)
    else:
        raise TypeError('motif_matching must be an ndarray or DataFrame.')
    if motif_df.shape[0] == 0 or motif_df.shape[1] == 0:
        raise ValueError('The motif_matching matrix cannot be empty.')

    df_cor = tf_gene_cor.rename(columns={'TF': 'tf', 'Gene': 'gene'}).copy()
    # Only the two columns tf/gene are required to exist; other statistical columns (correlation/p_value/fdr/t_stat) are retained as they are.
    required_cor_cols = {'tf', 'gene'}
    if not required_cor_cols.issubset(df_cor.columns):
        raise ValueError(f"tf_gene_cor needs to contain at least column {required_cor_cols}, current column: {list(df_cor.columns)}")

    # Ensure correct format: rows = peaks, columns = TFs
    # motif_df.index = peaks, motif_df.columns = TFs
    print(f"Input data:")
    print(f"  - Motif matching: {motif_df.shape[0]} peaks x {motif_df.shape[1]} TFs")
    print(f"  - Peak-gene links: {len(peak_gene_links)} links")
    print(f"  - TF-gene correlations: {len(tf_gene_cor) if not tf_gene_cor.empty else 0} pairs")
    



    vals = motif_df.values
    nz_r, nz_c = np.where(vals > 0)
    if len(nz_r) == 0:
        extra_cols = [c for c in df_cor.columns if c not in ('tf', 'gene')]
        empty_df = pd.DataFrame(columns=['tf', 'gene', 'n_peaks'] + extra_cols)
        return {'grn_edges': empty_df} if return_dict else empty_df

    df_p2m = pd.DataFrame({
        'peak': motif_df.index.to_numpy()[nz_r],
        'tf': motif_df.columns.to_numpy()[nz_c],
        'is_bound': vals[nz_r, nz_c]
    })
    df_p2m['peak'] = df_p2m['peak'].astype(str)
    df_p2m['tf'] = df_p2m['tf'].astype(str)

    print(f"  Found {len(df_p2m)} peak-TF binding events")
    print(f"  Unique TFs with binding sites: {df_p2m['tf'].nunique()}")
    print(f"  Unique peaks with TF binding: {df_p2m['peak'].nunique()}")
    

    df_p2g = peak_gene_links[['peak', 'gene']].copy()
    df_p2g['peak'] = df_p2g['peak'].astype(str)
    df_p2g['gene'] = df_p2g['gene'].astype(str)

    df_m2g_raw = pd.merge(df_p2m, df_p2g, on='peak', how='left').dropna(subset=['gene'])
    if df_m2g_raw.empty:
        extra_cols = [c for c in df_cor.columns if c not in ('tf', 'gene')]
        empty_df = pd.DataFrame(columns=['tf', 'gene', 'n_peaks'] + extra_cols)
        return {'grn_edges': empty_df} if return_dict else empty_df

    df_m2g = (
        df_m2g_raw
        .groupby(['tf', 'gene'], as_index=False)
        .size()
        .rename(columns={'size': 'n_peaks'})
    )

    print(f"  Found {len(df_m2g)} TF-gene pairs with peak support")
    print(f"  TFs with target genes: {df_m2g['tf'].nunique()}")
    print(f"  Genes regulated by TFs: {df_m2g['gene'].nunique()}")
    print(f"  Average peaks per TF-gene pair: {df_m2g['n_peaks'].mean():.2f}")
    print(f"  Peak support range: {df_m2g['n_peaks'].min()}-{df_m2g['n_peaks'].max()}")

    # --- Step 3: Left join with correlation (retain only (tf,gene) of df.m2g) ---
    print("\nStep 3: Adding TF-gene correlation (left join)...")
    df_grn = pd.merge(df_m2g, df_cor, on=['tf', 'gene'], how='left')

    # --- Optional filtering: compatible with R workflow (not filtered here by default) ---
    if correlation_threshold is not None and 'correlation' in df_grn.columns:
        df_grn = df_grn[df_grn['correlation'].abs() >= float(correlation_threshold)].reset_index(drop=True)

    # ========================================================================
    # Summary
    # ========================================================================
    print("\n" + "=" * 70)
    print("GRN construction completed!")
    print("=" * 70)
    print(f"Results:")
    print(f"  - Total TF-gene pairs: {len(df_grn)}")
    print(f"  - Unique TFs: {df_grn['tf'].nunique()}")
    print(f"  - Unique genes: {df_grn['gene'].nunique()}")
    print(f"  - Average targets per TF: {len(df_grn) / df_grn['tf'].nunique():.1f}")
    print(f"  - Average regulators per gene: {len(df_grn) / df_grn['gene'].nunique():.1f}")
    print("=" * 70)

    # R style: directly return the complete connection table (do not do any threshold filtering/mapping here)
    # Optional: save SwitchTFI-required artifacts during baseline GRN construction.
    if bool(kwargs.get("save_switchtfi_inputs", False)):
        try:
            save_info = save_switchtfi_inputs_from_baseline(
                data=kwargs.get("switchtfi_data", kwargs.get("data", None)),
                grn_df=df_grn,
                output_dir=kwargs.get("switchtfi_output_dir", None),
                trajectory_key=kwargs.get("switchtfi_trajectory_key", "Trajectory"),
                start_state=kwargs.get("switchtfi_start_state", None),
                end_state=kwargs.get("switchtfi_end_state", None),
                label_col_name=kwargs.get("switchtfi_label_col_name", "switch_label"),
                save_multimodal_h5ad=kwargs.get("save_switchtfi_multiome_h5ad", True),
                multiome_h5ad_filename=kwargs.get(
                    "switchtfi_multiome_h5ad_filename",
                    "switchtfi_start_end_multiome.h5ad",
                ),
                save_multimodal_h5mu=kwargs.get("save_switchtfi_multiome_h5mu", False),
                multiome_h5mu_filename=kwargs.get(
                    "switchtfi_multiome_h5mu_filename",
                    "trajectory_multiome.h5mu",
                ),
                switchtfi_network_filename=kwargs.get(
                    "switchtfi_network_filename",
                    "grn_switchTFI.csv",
                ),
                save_supporting_tables=kwargs.get("save_switchtfi_supporting_tables", True),
                tf_gene_cor_df=tf_gene_cor,
                peak_gene_links_df=peak_gene_links,
                motif_matching=motif_matching,
                tf_gene_cor_filename=kwargs.get(
                    "switchtfi_tf_gene_cor_filename",
                    "tf_gene_correlation.csv",
                ),
                peak_gene_links_filename=kwargs.get(
                    "switchtfi_peak_gene_links_filename",
                    "peak_gene_links.csv",
                ),
                motif_matching_filename=kwargs.get(
                    "switchtfi_motif_matching_filename",
                    "motif_matching.csv",
                ),
                full_grn_filename=kwargs.get(
                    "switchtfi_full_grn_filename",
                    "grn_full_baseline.csv",
                ),
                require_multimodal=kwargs.get("switchtfi_require_multimodal", True),
            )
            print(f"  - SwitchTFI artifacts saved to: {save_info['output_dir']}")
        except Exception as exc:
            import traceback
            traceback.print_exc()
            if bool(kwargs.get("switchtfi_raise_on_error", False)):
                raise RuntimeError(
                    "GRN built successfully, but saving SwitchTFI artifacts failed."
                ) from exc
            warnings.warn(
                f"GRN built successfully, but saving SwitchTFI artifacts failed: {exc}",
                RuntimeWarning,
            )

    return df_grn

def _resolve_start_end_states(
    state_series: pd.Series,
    start_state: Optional[str] = None,
    end_state: Optional[str] = None
) -> Tuple[str, str]:
    """Resolve start/end states from a trajectory state series."""
    states = state_series.dropna().astype(str).unique().tolist()
    if len(states) < 2:
        raise ValueError(
            f"Need at least 2 trajectory states, found {len(states)}: {states}"
        )

    resolved_start = str(states[0]) if start_state is None else str(start_state)
    resolved_end = str(states[-1]) if end_state is None else str(end_state)

    if resolved_start not in states:
        raise ValueError(f"start_state '{resolved_start}' not found in trajectory states: {states}")
    if resolved_end not in states:
        raise ValueError(f"end_state '{resolved_end}' not found in trajectory states: {states}")
    if resolved_start == resolved_end:
        raise ValueError("start_state and end_state must be different.")

    return resolved_start, resolved_end


def _resolve_start_end_mask_from_context(
    obs_df: pd.DataFrame,
    uns_dict: Dict[str, Any],
    trajectory_key: str,
    start_state: Optional[str],
    end_state: Optional[str]
) -> Tuple[np.ndarray, str, str, str]:
    """
    Resolve start/end cells using trajectory metadata first.

    Strict behavior:
    - Must resolve from ``uns['trajectory_analysis'][trajectory_key]`` with
      ``group_by`` + ``trajectory_groups``.
    - No fallback to direct ``obs[trajectory_key]`` matching, because trajectory
      values are often continuous pseudotime and not cell-type labels.
    """

    try:
        traj_meta = uns_dict.get('trajectory_analysis', {}).get(trajectory_key, {})
    except Exception:
        traj_meta = {}

    group_by = traj_meta.get('group_by', None)
    trajectory_groups = traj_meta.get('trajectory_groups', None)

    if (
        isinstance(group_by, str)
        and group_by in obs_df.columns
        and isinstance(trajectory_groups, (list, tuple))
        and len(trajectory_groups) >= 2
    ):
        group_series = obs_df[group_by].astype(str)
        states = group_series.dropna().unique().tolist()

        resolved_start = str(trajectory_groups[0]) if start_state is None else str(start_state)
        resolved_end = str(trajectory_groups[-1]) if end_state is None else str(end_state)

        if resolved_start not in states:
            raise ValueError(f"start_state '{resolved_start}' not found in '{group_by}' states: {states}")
        if resolved_end not in states:
            raise ValueError(f"end_state '{resolved_end}' not found in '{group_by}' states: {states}")
        if resolved_start == resolved_end:
            raise ValueError("start_state and end_state must be different.")

        mask = group_series.isin([resolved_start, resolved_end]).to_numpy()
        return mask, resolved_start, resolved_end, group_by

    raise ValueError(
        f"Cannot resolve start/end cell types for trajectory_key='{trajectory_key}'. "
        "Expected trajectory metadata at uns['trajectory_analysis'][trajectory_key] with "
        "'group_by' and 'trajectory_groups'."
    )


def _build_start_end_multiome_h5ad(
    data: Union[ad.AnnData, 'mu.MuData'],
    trajectory_key: str,
    start_state: Optional[str],
    end_state: Optional[str],
    label_col_name: str = "switch_label"
) -> Tuple[ad.AnnData, Dict[str, str]]:
    """
    Build a start/end-cell AnnData snapshot for SwitchTFI.

    For MuData input:
    - Uses RNA as primary AnnData matrix
    - Stores ATAC matrix in obsm['X_atac'] when available
    - Stores ATAC feature names in uns['atac_var_names']
    """
    if hasattr(data, "mod"):
        if 'rna' not in data.mod:
            raise ValueError("MuData missing 'rna' modality.")
        rna = data['rna']
        mask, resolved_start, resolved_end, state_key_used = _resolve_start_end_mask_from_context(
            obs_df=rna.obs,
            uns_dict=data.uns if hasattr(data, 'uns') else {},
            trajectory_key=trajectory_key,
            start_state=start_state,
            end_state=end_state
        )
        if mask.sum() == 0:
            raise ValueError("No cells found for start/end states in MuData.")

        rna_subset = rna[mask].copy()
        selected_obs = rna_subset.obs_names
        if not selected_obs.is_unique:
            raise ValueError("RNA obs_names must be unique to align ATAC by cell name.")
        adata_out = rna_subset.copy()
        if "counts" not in rna_subset.layers:
            raise ValueError(
                "RNA raw counts layer is required to export unnormalized SwitchTFI input, "
                "but rna.layers['counts'] was not found."
            )
        adata_out.X = rna_subset.layers["counts"].copy()
        adata_out.obs[label_col_name] = (
            rna_subset.obs[state_key_used].astype(str)
            .map({resolved_start: "start", resolved_end: "end"})
            .astype("category")
        )

        if 'atac' in data.mod:
            missing_atac_obs = selected_obs.difference(data['atac'].obs_names)
            if len(missing_atac_obs) > 0:
                raise ValueError(
                    "ATAC modality is missing cells required for SwitchTFI export: "
                    f"{list(missing_atac_obs[:5])}"
                )
            atac_subset = data['atac'][selected_obs].copy()
            if not np.array_equal(atac_subset.obs_names, selected_obs):
                raise ValueError("RNA and ATAC cells are not aligned after name-based subsetting.")
            if "counts" not in atac_subset.layers:
                raise ValueError(
                    "ATAC raw counts layer is required to export unnormalized SwitchTFI input, "
                    "but atac.layers['counts'] was not found."
                )
            adata_out.obsm['X_atac'] = atac_subset.layers["counts"].copy()
            adata_out.uns['atac_var_names'] = atac_subset.var_names.astype(str).tolist()

        # Expose TF activity for downstream improved SwitchTFI ranking.
        if 'chromvar' in data.mod:
            chromvar = data['chromvar']
            if chromvar.n_obs == rna.n_obs and np.array_equal(chromvar.obs_names, rna.obs_names):
                chromvar_subset = chromvar[mask].copy()
            else:
                missing_chromvar_obs = selected_obs.difference(chromvar.obs_names)
                if len(missing_chromvar_obs) > 0:
                    raise ValueError(
                        "chromVAR activity is missing cells required for SwitchTFI export: "
                        f"{list(missing_chromvar_obs[:5])}"
                    )
                chromvar_subset = chromvar[selected_obs].copy()
            chromvar_subset = _apply_motif_names_from_atac(chromvar_subset, data)
            chromvar_agg = _aggregate_chromvar_activity_by_tf(chromvar_subset)
            adata_out.obsm['tf_activity'] = _ensure_dense_array(chromvar_agg.X)
            adata_out.uns['tf_activity_genes'] = chromvar_agg.var_names.astype(str).tolist()
            adata_out.uns['tf_activity_aggregation'] = chromvar_agg.uns.get(
                'tf_activity_aggregation',
                {},
            )

        # Keep a unified pseudotime column when possible.
        if 'pseudotime' in rna_subset.obs.columns:
            adata_out.obs['pseudotime'] = rna_subset.obs['pseudotime'].to_numpy()
        else:
            traj_vals = pd.to_numeric(rna_subset.obs[trajectory_key], errors='coerce')
            valid = ~traj_vals.isna()
            if valid.any():
                min_v = float(traj_vals[valid].min())
                max_v = float(traj_vals[valid].max())
                if max_v > min_v:
                    adata_out.obs['pseudotime'] = (traj_vals - min_v) / (max_v - min_v)
                else:
                    adata_out.obs['pseudotime'] = 0.0

        return adata_out, {
            "start_state": resolved_start,
            "end_state": resolved_end,
            "state_key_used": state_key_used
        }

    if not isinstance(data, ad.AnnData):
        raise TypeError("data must be AnnData or MuData.")
    mask, resolved_start, resolved_end, state_key_used = _resolve_start_end_mask_from_context(
        obs_df=data.obs,
        uns_dict=data.uns if hasattr(data, 'uns') else {},
        trajectory_key=trajectory_key,
        start_state=start_state,
        end_state=end_state
    )
    if mask.sum() == 0:
        raise ValueError("No cells found for start/end states in AnnData.")

    adata_out = data[mask].copy()
    adata_out.obs[label_col_name] = (
        adata_out.obs[state_key_used].astype(str)
        .map({resolved_start: "start", resolved_end: "end"})
        .astype("category")
    )

    # Normalize ATAC obsm key to X_atac when possible.
    if 'X_atac' not in adata_out.obsm:
        for candidate in ('atac', 'ATAC', 'X_ATAC'):
            if candidate in adata_out.obsm:
                adata_out.obsm['X_atac'] = adata_out.obsm[candidate]
                break

    if 'tf_activity' not in adata_out.obsm and 'X_chromvar' in adata_out.obsm:
        adata_out.obsm['tf_activity'] = _ensure_dense_array(adata_out.obsm['X_chromvar'])
        if 'tf_activity_genes' not in adata_out.uns:
            adata_out.uns['tf_activity_genes'] = adata_out.var_names.astype(str).tolist()

    if 'pseudotime' not in adata_out.obs.columns:
        traj_vals = pd.to_numeric(adata_out.obs[trajectory_key], errors='coerce')
        valid = ~traj_vals.isna()
        if valid.any():
            min_v = float(traj_vals[valid].min())
            max_v = float(traj_vals[valid].max())
            if max_v > min_v:
                adata_out.obs['pseudotime'] = (traj_vals - min_v) / (max_v - min_v)
            else:
                adata_out.obs['pseudotime'] = 0.0

    return adata_out, {
        "start_state": resolved_start,
        "end_state": resolved_end,
        "state_key_used": state_key_used
    }


def _make_mudata(modalities: Dict[str, ad.AnnData]):
    """Create a MuData object using whichever MuData implementation is available."""
    try:
        import mudata as md  # type: ignore
        return md.MuData(modalities)
    except Exception:
        pass

    if MUON_AVAILABLE and mu is not None:
        return mu.MuData(modalities)

    raise ImportError("Saving trajectory_multiome.h5mu requires the 'mudata' or 'muon' package.")


def _build_start_end_multiome_h5mu(
    data: Union[ad.AnnData, 'mu.MuData'],
    trajectory_key: str,
    start_state: Optional[str],
    end_state: Optional[str],
    label_col_name: str = "switch_label"
) -> Tuple[Any, Dict[str, str]]:
    """
    Build a start/end-cell MuData snapshot matching the R export layout:
    modalities are RNA, ATAC, and required chromvar_activity.
    """
    if not hasattr(data, "mod"):
        raise TypeError("trajectory_multiome.h5mu export requires MuData input.")
    if 'rna' not in data.mod:
        raise ValueError("MuData missing 'rna' modality.")

    rna = data['rna']
    mask, resolved_start, resolved_end, state_key_used = _resolve_start_end_mask_from_context(
        obs_df=rna.obs,
        uns_dict=data.uns if hasattr(data, 'uns') else {},
        trajectory_key=trajectory_key,
        start_state=start_state,
        end_state=end_state
    )
    if mask.sum() == 0:
        raise ValueError("No cells found for start/end states in MuData.")

    rna_subset = rna[mask].copy()
    selected_obs = rna_subset.obs_names
    if not selected_obs.is_unique:
        raise ValueError("RNA obs_names must be unique to align ATAC by cell name.")
    if "counts" not in rna_subset.layers:
        raise ValueError(
            "RNA raw counts layer is required to export trajectory_multiome.h5mu, "
            "but rna.layers['counts'] was not found."
        )

    rna_out = rna_subset.copy()
    rna_out.X = rna_subset.layers["counts"].copy()
    rna_out.obs[label_col_name] = (
        rna_subset.obs[state_key_used].astype(str)
        .map({resolved_start: "start", resolved_end: "end"})
        .astype("category")
    )

    if 'pseudotime' in rna_subset.obs.columns:
        rna_out.obs['pseudotime'] = rna_subset.obs['pseudotime'].to_numpy()
    else:
        traj_vals = pd.to_numeric(rna_subset.obs[trajectory_key], errors='coerce')
        valid = ~traj_vals.isna()
        if valid.any():
            min_v = float(traj_vals[valid].min())
            max_v = float(traj_vals[valid].max())
            if max_v > min_v:
                rna_out.obs['pseudotime'] = (traj_vals - min_v) / (max_v - min_v)
            else:
                rna_out.obs['pseudotime'] = 0.0
    if 'pseudotime' not in rna_out.obs.columns:
        raise ValueError(
            "trajectory_multiome.h5mu requires obs['pseudotime']. "
            f"Provide rna.obs['pseudotime'] or numeric obs['{trajectory_key}']."
        )
    if trajectory_key not in rna_out.obs.columns:
        raise ValueError(
            f"trajectory_multiome.h5mu requires obs['{trajectory_key}']."
        )

    modalities: Dict[str, ad.AnnData] = {"RNA": rna_out}

    if 'atac' in data.mod:
        missing_atac_obs = selected_obs.difference(data['atac'].obs_names)
        if len(missing_atac_obs) > 0:
            raise ValueError(
                "ATAC modality is missing cells required for trajectory_multiome.h5mu: "
                f"{list(missing_atac_obs[:5])}"
            )
        atac_subset = data['atac'][selected_obs].copy()
        if not np.array_equal(atac_subset.obs_names, selected_obs):
            raise ValueError("RNA and ATAC cells are not aligned after name-based subsetting.")
        if "counts" not in atac_subset.layers:
            raise ValueError(
                "ATAC raw counts layer is required to export trajectory_multiome.h5mu, "
                "but atac.layers['counts'] was not found."
            )
        atac_out = atac_subset.copy()
        atac_out.X = atac_subset.layers["counts"].copy()
        atac_out.obs = rna_out.obs.copy()
        modalities["ATAC"] = atac_out

    if 'chromvar' not in data.mod:
        raise ValueError(
            "trajectory_multiome.h5mu requires TF activity modality "
            "data.mod['chromvar'] to export mod['chromvar_activity'].X."
        )
    chromvar = data['chromvar']
    if chromvar.n_obs == rna.n_obs and np.array_equal(chromvar.obs_names, rna.obs_names):
        chromvar_subset = chromvar[mask].copy()
    else:
        missing_obs = selected_obs.difference(chromvar.obs_names)
        if len(missing_obs) > 0:
            raise ValueError(
                "chromVAR activity is missing cells required for trajectory_multiome.h5mu: "
                f"{list(missing_obs[:5])}"
            )
        chromvar_subset = chromvar[selected_obs].copy()
    chromvar_subset = _apply_motif_names_from_atac(chromvar_subset, data)
    chromvar_out = _aggregate_chromvar_activity_by_tf(chromvar_subset)
    chromvar_out.obs = pd.DataFrame(index=chromvar_out.obs_names.copy())
    modalities["chromvar_activity"] = chromvar_out

    mdata_out = _make_mudata(modalities)
    for obs_key in (trajectory_key, label_col_name, 'pseudotime'):
        mdata_out.obs[obs_key] = rna_out.obs[obs_key].to_numpy()

    return mdata_out, {
        "start_state": resolved_start,
        "end_state": resolved_end,
        "state_key_used": state_key_used
    }


def _to_switchtfi_network(grn_df: pd.DataFrame) -> pd.DataFrame:
    """Convert baseline GRN columns to SwitchTFI naming convention."""
    network = grn_df.copy()
    rename_map = {}
    if "tf" in network.columns:
        rename_map["tf"] = "TF"
    elif "TF" not in network.columns:
        raise ValueError("GRN table must contain 'tf' or 'TF' column.")

    if "gene" in network.columns:
        rename_map["gene"] = "target"
    elif "target" not in network.columns:
        raise ValueError("GRN table must contain 'gene' or 'target' column.")

    if rename_map:
        network = network.rename(columns=rename_map)
    return network


def save_switchtfi_inputs_from_baseline(
    data: Union[ad.AnnData, 'mu.MuData'],
    grn_df: pd.DataFrame,
    output_dir: Union[str, Path],
    trajectory_key: str = "Trajectory",
    start_state: Optional[str] = None,
    end_state: Optional[str] = None,
    label_col_name: str = "switch_label",
    save_multimodal_h5ad: bool = True,
    multiome_h5ad_filename: str = "switchtfi_start_end_multiome.h5ad",
    save_multimodal_h5mu: bool = False,
    multiome_h5mu_filename: str = "trajectory_multiome.h5mu",
    switchtfi_network_filename: str = "grn_switchTFI.csv",
    save_supporting_tables: bool = True,
    tf_gene_cor_df: Optional[pd.DataFrame] = None,
    peak_gene_links_df: Optional[pd.DataFrame] = None,
    motif_matching: Optional[Union[np.ndarray, pd.DataFrame]] = None,
    tf_gene_cor_filename: str = "tf_gene_correlation.csv",
    peak_gene_links_filename: str = "peak_gene_links.csv",
    motif_matching_filename: str = "motif_matching.csv",
    full_grn_filename: str = "grn_full_baseline.csv",
    require_multimodal: bool = True
) -> Dict[str, Any]:
    """
    Save SwitchTFI-required artifacts at baseline network construction stage.

    This function performs:
    1) Save start/end cell multiome snapshot as h5ad and/or h5mu
    2) Save baseline GRN in SwitchTFI network format (TF/target)
    3) Save full baseline GRN and optional supporting evidence tables
    4) Save manifest with key names and file paths
    """
    if output_dir is None:
        raise ValueError("output_dir is required to save SwitchTFI artifacts.")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths: Dict[str, str] = {}
    state_meta: Dict[str, Optional[str]] = {
        "start_state": start_state,
        "end_state": end_state,
    }

    if save_multimodal_h5ad:
        adata_switchtfi, resolved = _build_start_end_multiome_h5ad(
            data=data,
            trajectory_key=trajectory_key,
            start_state=start_state,
            end_state=end_state,
            label_col_name=label_col_name
        )

        if require_multimodal and 'X_atac' not in adata_switchtfi.obsm:
            raise ValueError(
                "Saved start/end h5ad is missing ATAC in obsm['X_atac']. "
                "This workflow requires RNA+ATAC in one .h5ad."
            )

        h5ad_path = output_dir / multiome_h5ad_filename
        adata_switchtfi.write_h5ad(h5ad_path)
        paths["multiome_h5ad"] = str(h5ad_path)
        state_meta.update(resolved)

    if save_multimodal_h5mu:
        mdata_switchtfi, resolved = _build_start_end_multiome_h5mu(
            data=data,
            trajectory_key=trajectory_key,
            start_state=start_state,
            end_state=end_state,
            label_col_name=label_col_name
        )

        if require_multimodal and 'ATAC' not in mdata_switchtfi.mod:
            raise ValueError(
                "Saved trajectory h5mu is missing ATAC modality. "
                "This workflow requires RNA+ATAC in one .h5mu."
            )
        if require_multimodal and 'chromvar_activity' not in mdata_switchtfi.mod:
            raise ValueError(
                "Saved trajectory h5mu is missing chromvar_activity modality. "
                "This workflow requires TF activity in one .h5mu."
            )
        if require_multimodal and 'pseudotime' not in mdata_switchtfi.obs.columns:
            raise ValueError(
                "Saved trajectory h5mu is missing obs['pseudotime']. "
                "This workflow requires pseudotime in one .h5mu."
            )

        h5mu_path = output_dir / multiome_h5mu_filename
        mdata_switchtfi.write_h5mu(h5mu_path)
        paths["trajectory_multiome_h5mu"] = str(h5mu_path)
        state_meta.update(resolved)

    # Save full baseline GRN (all columns) for reproducibility.
    full_grn_path = output_dir / full_grn_filename
    grn_df.to_csv(full_grn_path, index=False)
    paths["baseline_grn_full_csv"] = str(full_grn_path)

    grn_switchtfi = _to_switchtfi_network(grn_df)
    network_path = output_dir / switchtfi_network_filename
    grn_switchtfi.to_csv(network_path, index=False)
    paths["switchtfi_network_csv"] = str(network_path)

    if save_supporting_tables:
        if tf_gene_cor_df is not None:
            tf_cor_path = output_dir / tf_gene_cor_filename
            tf_gene_cor_df.to_csv(tf_cor_path, index=False)
            paths["tf_gene_correlation_csv"] = str(tf_cor_path)

        if peak_gene_links_df is not None:
            p2g_path = output_dir / peak_gene_links_filename
            peak_gene_links_df.to_csv(p2g_path, index=False)
            paths["peak_gene_links_csv"] = str(p2g_path)

        if motif_matching is not None:
            motif_path = output_dir / motif_matching_filename
            if isinstance(motif_matching, pd.DataFrame):
                motif_matching.to_csv(motif_path)
            else:
                pd.DataFrame(motif_matching).to_csv(motif_path, index=False)
            paths["motif_matching_csv"] = str(motif_path)

    h5ad_keys: Dict[str, str] = {}
    if save_multimodal_h5ad:
        h5ad_keys = {
            "rna_matrix": "X",
            "atac_matrix": "obsm['X_atac']",
            "trajectory_key": trajectory_key,
            "label_key": label_col_name,
            "tf_activity_key": "obsm['tf_activity'] (if present)",
            "pseudotime_key": "obs['pseudotime'] (if present)",
            "tf_activity_genes_key": "uns['tf_activity_genes'] (if present)",
            "atac_var_names_key": "uns['atac_var_names'] (if present)"
        }

    h5mu_keys: Dict[str, str] = {}
    if save_multimodal_h5mu:
        h5mu_keys = {
            "rna_modality": "RNA",
            "atac_modality": "ATAC",
            "tf_activity_modality": "chromvar_activity",
            "rna_matrix": "mod['RNA'].X",
            "atac_matrix": "mod['ATAC'].X",
            "tf_activity_matrix": "mod['chromvar_activity'].X",
            "trajectory_key": trajectory_key,
            "label_key": label_col_name,
            "pseudotime_key": "obs['pseudotime']"
        }

    manifest = {
        "trajectory_key": trajectory_key,
        "start_state": state_meta.get("start_state"),
        "end_state": state_meta.get("end_state"),
        "state_key_used": state_meta.get("state_key_used"),
        "label_col_name": label_col_name,
        "require_multimodal": require_multimodal,
        "saved_supporting_tables": bool(save_supporting_tables),
        "h5ad_keys": h5ad_keys,
        "h5mu_keys": h5mu_keys,
        "paths": paths
    }
    manifest_path = output_dir / "switchtfi_inputs_manifest.json"
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, ensure_ascii=False, indent=2)
    paths["manifest_json"] = str(manifest_path)

    return {
        "output_dir": str(output_dir),
        "start_state": manifest["start_state"],
        "end_state": manifest["end_state"],
        "paths": paths
    }


def get_tf_gene_correlation(data: Union[ad.AnnData, 'mu.MuData'],
                           tf_list: Optional[List[str]] = None,
                           gene_list: Optional[List[str]] = None,
                           method: str = "pearson",
                           trajectory_name: str = "Trajectory",
                           group_every: int = 1,
                           smooth_window: int = 7,
                           tf_assay: str = "chromvar",
                           gene_assay: str = "RNA",
                           atac_assay: str = "ATAC",
                           motif_names: Optional[Dict] = None,
                           return_heatmap: bool = False,
                           **kwargs) -> pd.DataFrame:
    """
    Calculate TF-gene correlations along trajectory - EXACT R-compatible implementation.
    
    This function replicates scMEGA's GetTFGeneCorrelation (R/grn.R:19-112).
    
    R logic (step by step):
    1. GetTrajectory for TF activity (chromvar, log2Norm=FALSE)
    2. GetTrajectory for gene expression (RNA, log2Norm=TRUE)
    3. CRITICAL: Replace trajMM rownames with motif names
    4. TrajectoryHeatmap: variance filter + Z-score + clip to [-2, 2]
    5. Filter TFs and genes if specified
    6. Calculate correlation ALONG TIME BINS (not cells)
    7. Compute t-statistics and p-values
    
    Args:
        data: MultiomeData or AnnData object with trajectory information
        tf_list: List of TFs to analyze
        gene_list: List of target genes to analyze
        method: Correlation method ("pearson" matches R default)
        trajectory_name: Name of the trajectory
        group_every: Trajectory binning (default: 1)
        smooth_window: Smoothing window (default: 7)
        tf_assay: Assay for TF activity (default: "chromvar")
         gene_assay: Assay for gene expression (default: "RNA")
        atac_assay: ATAC assay name
        motif_names: **REQUIRED** Dictionary mapping feature indices/names to TF/motif names
        **kwargs: Additional parameters
        
    Returns:
        pd.DataFrame with columns [tf, gene, correlation, t_stat, p_value, fdr]
        
    Example:
        ```python
        # Must provide motif_names mapping!
        motif_dict = {
            'MA0098.3': 'ETS1',
            'MA0080.5': 'SPI1',
            # ... more mappings
        }
        
        tf_gene_cor = get_tf_gene_correlation(
            multiome,
            tf_list=['ETS1', 'SPI1', 'GATA1'],
            gene_list=selected_genes,
            trajectory_name="Trajectory",
            motif_names=motif_dict
        )
        ```
    """
    print("=" * 70)
    print("Computing TF-gene correlation (R GetTFGeneCorrelation equivalent)")
    print("=" * 70)
    
    from ..trajectory_analysis.pseudotime_analysis import get_trajectory_data
    
    
    # Deprecated parameter warning
    if motif_names is not None:
        warnings.warn(
            "The 'motif_names' parameter is deprecated. "
            "Motif names are now automatically extracted from data['atac'].uns['chromvar']['motif_names']. "
            "This parameter will be ignored.",
            DeprecationWarning
        )
    
    # Step 1: Get TF activity trajectory (chromvar)
    # Now automatically extracts from obsm['X_chromvar'] and applies motif names
    print("\nStep 1: Getting TF activity trajectory (chromvar)...")
    
    traj_mm = get_trajectory_data(
        data=data,
        trajectory_name=trajectory_name,
        assay=tf_assay,  # "chromvar" - automatically extract and apply motif names from obsm
        slot="X",
        group_every=group_every,
        log2_norm=False,
        scale_to=None,
        smooth_window=smooth_window,
        return_matrix=True
    )
    
    traj_mm_mat = traj_mm  # traj_mm is already a matrix, not a dictionary
    #traj_mm_mat = traj_mm['smooth_matrix'] if 'smooth_matrix' in traj_mm else traj_mm['group_matrix']
    print(f"  TF trajectory: {traj_mm_mat.shape[0]} TFs x {traj_mm_mat.shape[1]} time bins")
    print(f"   Motif names automatically applied")
    
    # Step 2: Get gene expression trajectory (RNA)
    print("\nStep 2: Getting gene expression trajectory (RNA)...")
    
    traj_rna = get_trajectory_data(
        data=data,
        trajectory_name=trajectory_name,
        assay=gene_assay,  # "RNA"
        slot="counts",     # read raw counts, then apply scale_to + log2_norm
        group_every=group_every,
        log2_norm=True,
        scale_to=10000,
        smooth_window=smooth_window,
        return_matrix=True
    )

    traj_rna_mat = traj_rna  # traj_rna is already a matrix, not a dictionary
    #traj_rna_mat = traj_rna['smooth_matrix'] if 'smooth_matrix' in traj_rna else traj_rna['group_matrix']
    print(f"  Gene trajectory: {traj_rna_mat.shape[0]} genes x {traj_rna_mat.shape[1]} time bins")
    
    # ========================================================================
    # Step 4: Apply TrajectoryHeatmap transformation
    # ========================================================================
    # R code lines 51-71
    
    print("\nStep 4: Applying TrajectoryHeatmap transformation...")
    
    # TF activity: varCutOff=0 (keep all), limits=(-2, 2)
    print("  Processing TF activity...")
    tf_activity = _trajectory_heatmap_r_exact(
        mat=traj_mm_mat,
        var_cutoff=0.0,
        max_features=25000,
        scale_rows=True,
        limits=(-2, 2)
    )
    print(f"    Result: {tf_activity.shape[0]} TFs x {tf_activity.shape[1]} time bins")
    
    # Gene expression: varCutOff=0.9 (top 10% by variance), limits=(-2, 2)
    print("  Processing gene expression...")
    gene_expression = _trajectory_heatmap_r_exact(
        mat=traj_rna_mat,
        var_cutoff=0.9,  # Retain the top ten percent of hypervariable genes
        max_features=25000,
        scale_rows=True,
        limits=(-2, 2)
    )
    print(f"    Result: {gene_expression.shape[0]} genes x {gene_expression.shape[1]} time bins")
    
    # ========================================================================
    # Step 5: Filter TFs and genes
    # ========================================================================
    # R code lines 74-84
    
    print("\nStep 5: Filtering features...")
    
    if tf_list is not None:
        # R: tf_activity <- tf_activity[tf.use, ]
        available_tfs = [tf for tf in tf_list if tf in tf_activity.index]
        if len(available_tfs) == 0:
            raise ValueError(f"None of the specified TFs found in TF activity matrix. "
                           f"Available TFs: {list(tf_activity.index[:10])}")
        tf_activity = tf_activity.loc[available_tfs]
        print(f"  Filtered to {len(available_tfs)} specified TFs")
    else:
        print(f"  Using all {len(tf_activity)} TFs")
    
    if gene_list is not None:
        # R: sel_genes <- intersect(rownames(gene_expression), gene.use)
        available_genes = list(set(gene_list) & set(gene_expression.index))
        if len(available_genes) == 0:
            raise ValueError(f"None of the specified genes found in gene expression matrix")
        gene_expression = gene_expression.loc[available_genes]
        print(f"  Filtered to {len(available_genes)} specified genes")
    else:
        print(f"  Using all {len(gene_expression)} genes")
    
    # ========================================================================
    # Step 6: Calculate correlations along trajectory time bins
    # ========================================================================
    # R code lines 86-99
    
    print("\nStep 6: Calculating correlations...")
    print(f"  Method: {method}")
    print(f"  Dimension: along {tf_activity.shape[1]} time bins")
    
    # R: df.cor <- t(cor(t(tf_activity), t(gene_expression)))
    # Breakdown:
    #   t(tf_activity): time_bins x TFs
    #   t(gene_expression): time_bins x genes  
    #   cor(...): computes correlation between columns -> TFs x genes
    #   t(...): transpose to genes x TFs
    #   Then convert to long format
    
    n_time_bins = tf_activity.shape[1]
    
    if n_time_bins < 3:
        raise ValueError(f"Insufficient time bins: {n_time_bins}. Need at least 3.")
    
    # Compute correlation matrix efficiently
    if method == "pearson":
        # ============================================================
        # ============================================================
        print("  Using np.corrcoef() for vectorized correlation calculation...")
        
        # combined shape: (n_genes + n_tfs, n_time_bins)
        combined = np.vstack([gene_expression.values, tf_activity.values])
        
        # cor_matrix shape: (n_genes + n_tfs, n_genes + n_tfs)
        cor_matrix = np.corrcoef(combined)
        
        n_genes = gene_expression.shape[0]
        n_tfs = tf_activity.shape[0]
        
        # cor_values shape: (n_genes, n_tfs)
        cor_values = cor_matrix[:n_genes, n_genes:]
        
        correlations_list = []
        for i, gene_name in enumerate(gene_expression.index):
            for j, tf_name in enumerate(tf_activity.index):
                correlation = cor_values[i, j]
                
                if not np.isnan(correlation):
                    correlations_list.append({
                        'tf': tf_name,
                        'gene': gene_name,
                        'correlation': correlation
                    })
        
        print(f"   Calculated {len(correlations_list)} correlations using vectorized operation")

    else:
        raise ValueError(f"Unknown method: {method}")
            

    
    if len(correlations_list) == 0:
        raise ValueError("No valid correlations computed")
    
    df_cor = pd.DataFrame(correlations_list)
    
    print(f"  Computed {len(df_cor)} TF-gene pairs")
    print(f"  Mean |correlation|: {df_cor['correlation'].abs().mean():.3f}")
    
    # ========================================================================
    # Step 7: Calculate statistics
    # ========================================================================
    # R code lines 101-108
    
    print("\nStep 7: Computing statistics...")
    
    # Degrees of freedom = number of time bins - 2
    df_stat = n_time_bins - 2
    
    # T-statistic: t = r * sqrt(df / (1 - r^2))
    # R code: (df.cor$correlation / sqrt((pmax(1 - df.cor$correlation^2, 1e-16)) / (ncol(tf_activity) - 2)))
    df_cor['t_stat'] = df_cor['correlation'] / np.sqrt(
        np.maximum(1 - df_cor['correlation'] ** 2, 1e-16) / df_stat
    )
    
    # P-value: two-tailed t-test
    # R code: 2 * pt(-abs(df.cor$t_stat), ncol(tf_activity) - 2)
    from scipy.stats import t as t_dist
    df_cor['p_value'] = 2 * t_dist.cdf(-np.abs(df_cor['t_stat']), df_stat)
    
    # FDR correction
    # R code: p.adjust(df.cor$p_value, method = "fdr")
    try:
        from scipy.stats import false_discovery_control
        df_cor['fdr'] = false_discovery_control(df_cor['p_value'].values, method='bh')
    except:
        # Fallback for older scipy
        from statsmodels.stats.multitest import multipletests
        _, fdr_vals, _, _ = multipletests(df_cor['p_value'].values, method='fdr_bh')
        df_cor['fdr'] = fdr_vals
    
    # Sort by absolute correlation (descending)
    df_cor = df_cor.sort_values('correlation', key=abs, ascending=False).reset_index(drop=True)
    
    # ========================================================================
    # Summary
    # ========================================================================
    print("\n" + "=" * 70)
    print("TF-gene correlation analysis completed!")
    print("=" * 70)
    print(f"Results:")
    print(f"  - Total TF-gene pairs: {len(df_cor)}")
    print(f"  - Unique TFs: {df_cor['tf'].nunique()}")
    print(f"  - Unique genes: {df_cor['gene'].nunique()}")
    print(f"  - Correlation range: [{df_cor['correlation'].min():.3f}, {df_cor['correlation'].max():.3f}]")
    print(f"  - Mean |correlation|: {df_cor['correlation'].abs().mean():.3f}")
    print(f"  - Significant (FDR<0.05): {(df_cor['fdr'] < 0.05).sum()}")
    print(f"  - Strong correlation (|r|>0.5): {(df_cor['correlation'].abs() > 0.5).sum()}")
    print("=" * 70)
    
    # If heatmap is not required, return DataFrame directly (backwards compatible)
    if not return_heatmap:
        return df_cor
    # ========================================================================
    # ========================================================================
    print("\nGenerating GRN correlation heatmap...")
    heatmap_fig = None
    try:
        from ..visualization.heatmaps import tf_gene_correlation_heatmap
        
        heatmap_fig = tf_gene_correlation_heatmap(
            tf_gene_cor=df_cor,  # TF-gene correlationDataFrame
            n_clusters=1,  # Corresponding to the km parameter of R
            figsize=(12, 10),
            show_row_names=False,  # Corresponds to R's show_row_names = FALSE
            show_column_names=True,  # Corresponds to R's show_column_names = TRUE
            cluster_method='ward',  # ward.D2 corresponding to R
            cmap='RdBu_r'  # Correlation color palette corresponding to R
        )
        
        print(f"   GRN heatmap generated")
        print(f"    TFs: {df_cor['tf'].nunique()}")
        print(f"    Genes: {df_cor['gene'].nunique()}")
        print(f"    Total pairs: {len(df_cor)}")
        
    # except Exception as e:
    #     print(f"   Heatmap generation failed: {e}")
    #     heatmap_fig = None
    except Exception as e:
        print(f"   Heatmap generation failed in grn_inference: {e}")
        import traceback
        traceback.print_exc()  # Print complete stack information
        heatmap_fig = None
    
    # ========================================================================
    # ========================================================================
    result = {
        'correlation': df_cor,  # TF-gene correlation data
        'heatmap': heatmap_fig  # GRN heatmap or None
    }
    
    print("\nReturning R-compatible dictionary:")
    print(f"  - 'correlation': DataFrame with {len(df_cor)} TF-gene pairs")
    print(f"  - 'heatmap': {'matplotlib.Figure' if heatmap_fig is not None else 'None'}")
    print("=" * 70)
    
    return result


def _trajectory_heatmap_r_exact(mat: pd.DataFrame,
                               var_cutoff: float = 0.9,
                               max_features: int = 25000,
                               scale_rows: bool = True,
                               limits: Tuple[float, float] = (-1.5, 1.5)) -> pd.DataFrame:
    """
    Exact replication of R's TrajectoryHeatmap transformation (R/visualization.R:559-664).
    
    Args:
        mat: Input matrix (features x time_bins)
        var_cutoff: Variance quantile cutoff (0-1)
        max_features: Maximum features to keep
        scale_rows: Whether to apply Z-score normalization
        limits: Clipping limits for scaled values
        
    Returns:
        Transformed matrix
    """
    # Step 1: Remove rows with NA values (R lines 574-578)
    na_rows = mat.isna().sum(axis=1)
    if (na_rows > 0).sum() > 0:
        mat = mat[na_rows == 0]
    
    # Step 2: Remove rows with zero standard deviation (R lines 580-584)
    row_stds = mat.std(axis=1)
    zero_std_mask = row_stds != 0
    if zero_std_mask.sum() < len(mat):
        mat = mat[zero_std_mask]
    
    if len(mat) == 0:
        raise ValueError("No features remaining after filtering NA and zero-variance rows")
    
    # Step 3: Calculate variance and filter (R lines 586-604)
    row_vars = mat.var(axis=1)
    
    # Convert to quantiles (matching ArchR:::.getQuantiles)
    from scipy.stats import rankdata
    var_quantiles = rankdata(row_vars, method='ordinal') / len(row_vars)
    
    # Order by variance (descending)
    var_order = np.argsort(row_vars.values)[::-1]
    mat_ordered = mat.iloc[var_order]
    
    # Calculate number of features to keep
    if var_cutoff is None and max_features is None:
        n_keep = len(mat_ordered)
    elif var_cutoff is None:
        n_keep = max_features
    elif max_features is None:
        n_keep = int((1 - var_cutoff) * len(mat_ordered))
    else:
        n_keep = min(int((1 - var_cutoff) * len(mat_ordered)), max_features)
    
    n_keep = min(n_keep, len(mat_ordered))
    mat_filtered = mat_ordered.iloc[:n_keep]
    
    # Step 4: Scale rows (Z-score) and clip (R lines 618-622)
    if scale_rows:
        # Row-wise Z-score: (x - mean) / std
        row_means = mat_filtered.mean(axis=1)
        row_stds = mat_filtered.std(axis=1)
        
        # Avoid division by zero
        row_stds[row_stds == 0] = 1
        
        mat_scaled = mat_filtered.sub(row_means, axis=0).div(row_stds, axis=0)
        
        # Clip to limits
        mat_scaled = mat_scaled.clip(lower=limits[0], upper=limits[1])
        
        return mat_scaled
    else:
        return mat_filtered
