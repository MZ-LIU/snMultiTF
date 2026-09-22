"""
4.2_celltype_writeback.py
====================
Add cell type information based on original tag metadata
Read cell-by-cell annotations in cell_metadata.csv,
Write cell_type back to MuData and save the annotated data.
"""

import os
import sys
from pathlib import Path

import pandas as pd
import muon as mu
import anndata as ad


current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, ".."))

if project_root not in sys.path:
    sys.path.insert(0, project_root)
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)


# ============================================================
# ============================================================
RESULT_BASE_DIR = r'E:\日志\snMultiTF\workflow_labeled\results\pbmc_multiome'

ANNOTATION_RESULT_PATH = r"E:/data/10x_pbmc_celltype/cell_metadata.csv"

CELLTYPE_ANNOTATION_DIR = str(Path(RESULT_BASE_DIR) / "celltype_annotation")

MUDATA_PATH = str(Path(RESULT_BASE_DIR) / "mdata_qc.h5mu")

MDATA_WITH_CELLTYPE_PATH = str(Path(RESULT_BASE_DIR) / "mdata_with_celltype.h5mu")

def find_annotation_dir():
    result_file = Path(ANNOTATION_RESULT_PATH)
    if not result_file.exists():
        raise FileNotFoundError(f"Annotation result file not found: {result_file}")
    output_dir = Path(CELLTYPE_ANNOTATION_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def load_cell_metadata(annotation_file):
    metadata = pd.read_csv(annotation_file)

    id_candidates = ["cell_id", "cell_name", "barcode", "Unnamed: 0"]
    id_col = next((col for col in id_candidates if col in metadata.columns), None)
    if id_col is None:
        id_col = metadata.columns[0]

    if "celltype" not in metadata.columns:
        raise ValueError('Required celltype column not found in cell_metadata.csv')

    metadata[id_col] = metadata[id_col].astype(str).str.strip()
    metadata = metadata.drop_duplicates(subset=id_col, keep="first").set_index(id_col)
    metadata["celltype"] = metadata["celltype"].fillna("Unknown").astype(str).str.strip()
    metadata.loc[metadata["celltype"].eq(""), "celltype"] = "Unknown"

    if "broad_celltype" in metadata.columns:
        metadata["broad_cell_type"] = (
            metadata["broad_celltype"].fillna("Unknown").astype(str).str.strip()
        )

    return metadata


def map_metadata_to_obs(obs_names, metadata, column, default="Unknown"):
    obs_index = pd.Index(obs_names.astype(str).str.strip())
    mapped = obs_index.to_series(index=obs_names).map(metadata[column])

    if mapped.isna().any():
        normalized_metadata = metadata.copy()
        normalized_metadata.index = normalized_metadata.index.astype(str).str.replace(
            r"-1$", "", regex=True
        )
        normalized_metadata = normalized_metadata[~normalized_metadata.index.duplicated(keep="first")]
        normalized_obs = obs_index.str.replace(r"-1$", "", regex=True)
        normalized_mapped = pd.Series(normalized_obs, index=obs_names).map(
            normalized_metadata[column]
        )
        mapped = mapped.fillna(normalized_mapped)

    return mapped.fillna(default).astype(str)


def main():
    print("=" * 70)
    print('4.2 Cell type annotations written back to MuData')
    print("=" * 70)

    output_dir = find_annotation_dir()
    annotation_file = Path(ANNOTATION_RESULT_PATH)

    print(f"Annotation result directory: {output_dir}")
    print(f"Annotation result file: {annotation_file}")

    cell_metadata = load_cell_metadata(annotation_file)

    print('\nSummary of annotation results:')
    print(f"  - cell_metadata Cell number: {len(cell_metadata)}")
    print(f"  - Number of cell types: {cell_metadata['celltype'].nunique()} types")
    print('  - Distribution of each cell type in cell_metadata:')
    for celltype, count in cell_metadata["celltype"].value_counts().head(10).items():
        print(f"    · {celltype}: {count} cells")

    file_ext = os.path.splitext(MUDATA_PATH)[1].lower()
    if file_ext == ".h5mu":
        mdata = mu.read_h5mu(MUDATA_PATH)
    elif file_ext == ".h5ad":
        adata = ad.read_h5ad(MUDATA_PATH)
        mdata = mu.MuData({"rna": adata})
    else:
        raise ValueError(f"Unsupported file format: {file_ext}")

    mdata["rna"].obs["celltype"] = map_metadata_to_obs(
        mdata["rna"].obs_names, cell_metadata, "celltype"
    )
    mdata["rna"].obs["cell_type"] = mdata["rna"].obs["celltype"]
    if "broad_cell_type" in cell_metadata.columns:
        mdata["rna"].obs["broad_cell_type"] = map_metadata_to_obs(
            mdata["rna"].obs_names, cell_metadata, "broad_cell_type"
        )

    unmatched_rna = mdata["rna"].obs_names[mdata["rna"].obs["celltype"].eq("Unknown")]
    print(f"  - RNA cell number: {mdata['rna'].n_obs}")
    print(f"  - Number of RNA cells not matched to cell_metadata: {len(unmatched_rna)}")
    if len(unmatched_rna) > 0:
        pd.Series(unmatched_rna, name="rna_cell_id").to_csv(
            output_dir / "unmatched_rna_cells.csv",
            index=False,
        )
        print(f"  - Unmatched RNA cell list saved: {output_dir / 'unmatched_rna_cells.csv'}")

    if "atac" in mdata.mod:
        mdata["atac"].obs["celltype"] = map_metadata_to_obs(
            mdata["atac"].obs_names, cell_metadata, "celltype"
        )
        mdata["atac"].obs["cell_type"] = mdata["atac"].obs["celltype"]
        if "broad_cell_type" in cell_metadata.columns:
            mdata["atac"].obs["broad_cell_type"] = map_metadata_to_obs(
                mdata["atac"].obs_names, cell_metadata, "broad_cell_type"
            )

        unmatched_atac = mdata["atac"].obs_names[mdata["atac"].obs["celltype"].eq("Unknown")]
        print(f"  - ATAC cell number: {mdata['atac'].n_obs}")
        print(f"  - Number of ATAC cells not matched to cell_metadata: {len(unmatched_atac)}")
        if len(unmatched_atac) > 0:
            preview = unmatched_atac[:10].tolist()
            print(f"  - Examples of unmatched ATAC cells (up to 10): {preview}")
            pd.Series(unmatched_atac, name="atac_cell_id").to_csv(
                output_dir / "unmatched_atac_cells.csv",
                index=False,
            )
            print(f"  - Unmatched ATAC cell list saved: {output_dir / 'unmatched_atac_cells.csv'}")

    mdata.obs["celltype"] = map_metadata_to_obs(
        mdata.obs_names, cell_metadata, "celltype"
    )
    mdata.obs["cell_type"] = mdata.obs["celltype"]
    if "broad_cell_type" in cell_metadata.columns:
        mdata.obs["broad_cell_type"] = map_metadata_to_obs(
            mdata.obs_names, cell_metadata, "broad_cell_type"
        )

    print('\n Cell type annotation write back completed')
    print(f"  - Number of identified cell types: {mdata['rna'].obs['cell_type'].nunique()} species")
    print('  - Distribution of each cell type:')
    for celltype, count in mdata["rna"].obs["cell_type"].value_counts().head(10).items():
        print(f"    · {celltype}: {count} cells ({count / len(mdata['rna'].obs) * 100:.1f}%)")

    export_cols = [
        col for col in ["celltype", "cell_type", "broad_cell_type"]
        if col in mdata["rna"].obs.columns
    ]
    celltype_results = mdata["rna"].obs[export_cols].copy()
    celltype_results.index.name = "cell_id"
    celltype_results.to_csv(output_dir / "cell_type_annotation_results.csv")
    print(f"   Cell type annotation results have been saved: {output_dir / 'cell_type_annotation_results.csv'}")

    for mod in ["rna", "atac"]:
        if mod in mdata.mod:
            if "counts" not in mdata[mod].layers:
                mdata[mod].layers["counts"] = mdata[mod].X.copy()
                print(f"  - Added layers['counts'] for {mod.upper()}")
            else:
                print(f"  - {mod.upper()} layers['counts'] already exists and will not be overwritten repeatedly.")

    mdata_output_path = Path(MDATA_WITH_CELLTYPE_PATH)
    mdata_output_path.parent.mkdir(parents=True, exist_ok=True)
    mdata.write_h5mu(mdata_output_path)
    print(f"   Annotated MuData saved: {mdata_output_path}")

    # 4.2 The script is only responsible for writing the cell-by-cell annotations in cell_metadata.csv back to MuData after quality control.
    # Integrated, UMAP plot colored by cell_type generated by 5_network.py after multi-omics integration.
    print('  - RNA/ATAC UMAP plots in 4.2 have been skipped.')
    print("=" * 70)
    print(f"The write-back process is completed and the results are saved to: {mdata_output_path.parent}")
    print("=" * 70)


if __name__ == "__main__":
    main()
