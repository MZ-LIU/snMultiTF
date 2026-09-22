"""
1_Quality control index calculation.py
================
Load label-free multi-omics data (RNA + ATAC), calculate QC metrics and generate QC visualization charts before filtering.
Provide 1_data quality control.py with reference to the chart to determine the filtering threshold.
"""

import os
import sys
import argparse
import json
from pathlib import Path

import muon as mu

# Compatible with package import paths when running scripts directly
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_DIR, ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import python_scmega as pymega


def parse_args():
    """Parse command line arguments. Supports two modes: online tool calling and direct operation."""
    parser = argparse.ArgumentParser(description='Calculation of quality control indicators')
    parser.add_argument("--config", type=str, help='parameters.json path')
    parser.add_argument("--run_dir", type=str, help='Run output directory')
    parser.add_argument("--result_dir", type=str, help='Results directory')
    return parser.parse_args()


def resolve_paths(args):
    """Determine the data path and output directory based on parameters."""
    if args.config and args.run_dir:
        with open(args.config, "r", encoding="utf-8") as f:
            params = json.load(f)
        data_path = params.get("h5mu_path", "")
        run_dir = Path(args.run_dir)
        output_dir = run_dir / "plots"
        output_dir.mkdir(parents=True, exist_ok=True)
        return data_path, None, output_dir
    else:
        data_path = r"E:\data\Pancreas\GSM5979673_week2.filtered_feature_bc_matrix.h5"
        fragments_path = r"E:\data\Pancreas\GSM5979672_week2.atac_fragments.tsv.gz"
        output_dir = Path(CURRENT_DIR) / "results_sc" / "sc_multiome"
        output_dir.mkdir(parents=True, exist_ok=True)
        return data_path, fragments_path, output_dir


# =========================
# =========================
def load_mudata(data_path: str, fragments_path: str) -> mu.MuData:
    """Load 10X multiome h5 data and record fragments path."""
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"The data file {data_path} does not exist, please check whether the path is correct.")
    if not os.path.exists(fragments_path):
        raise FileNotFoundError(f"The fragment file {fragments_path} does not exist, please check whether the path is correct.")

    mdata = mu.read_10x_h5(data_path)
    print(f"Included modalities: {mdata.mod.keys()}")

    for mod in mdata.mod:
        if mdata[mod].var_names.duplicated().any():
            n_dup = mdata[mod].var_names.duplicated().sum()
            print(f"  {mod}: {n_dup} duplicate variable names were found and are being removed...")
            mdata[mod].var_names_make_unique()

    mdata["atac"].uns["fragments"] = {
        "path": str(Path(fragments_path).absolute()),
        "genome": "hg38",
    }
    print(f"  Fragments path has been recorded: {mdata['atac'].uns['fragments']['path']}")

    print(f"Raw data overview:")
    print(f"  RNA data shape: {mdata['rna'].shape}")
    print(f"  ATAC data shape: {mdata['atac'].shape}")
    return mdata


# =========================
# =========================
def main() -> None:
    args = parse_args()
    data_path, fragments_path, output_dir = resolve_paths(args)

    print("=" * 70)
    print('Quality control index calculation (RNA + ATAC)')
    print("=" * 70)
    print(f"Data path: {data_path}")
    print(f"Output directory: {output_dir}")

    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    print('\n[1/3] Loading 10X multiome data...')
    mdata = load_mudata(data_path, fragments_path)
    rna = mdata["rna"]

    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    print('\n[2/3] Calculate RNA quality control indicators...')
    pymega.calculate_rna_qc_metrics(rna, mt_pattern="^MT-", inplace=True)

    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    print('\n[3/3] Generate QC visualization charts...')
    pymega.plot_qc_metrics(
        rna,
        save_path=str(output_dir / "QC_metrics.png"),
    )
    print(f"  QC chart saved: {output_dir / 'QC_metrics.png'}")
    print('\nPlease determine the filtering threshold based on the QC chart, and then run 1_dataQC.py.')


if __name__ == "__main__":
    main()
