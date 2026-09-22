# snMultiTF

snMultiTF is an integrated framework for single-nucleus multi-omics analysis, combining ensemble consensus clustering, LLM-assisted cell-type annotation, trajectory-associated gene regulatory network (GRN) inference, and regulatory transcription factor (TF) identification.

The `workflow_unlabeled` pipeline is designed for datasets without predefined cell-type labels, where cell clusters are identified by ensemble consensus clustering and annotated with snTypeGPT. The `workflow_labeled` pipeline is designed for datasets with predefined cell-type labels, which are incorporated directly into the downstream workflow for multi-omics integration, trajectory-associated GRN inference, and regulatory TF identification.

## Overview

The framework contains four main components:

1. **Ensemble consensus clustering** for stable cell clustering.
2. **snTypeGPT** for marker- and LLM-assisted cell-type annotation.
3. **Trajectory-associated GRN inference** by integrating RNA expression and ATAC chromatin accessibility.
4. **NACR** for regulatory TF identification.

## Repository Structure

```text
snMultiTF/
├── workflow_unlabeled/        # Main workflow for datasets without predefined cell-type labels
├── workflow_labeled/          # Main workflow for datasets with predefined cell-type labels
├── ensemble_clustering/       # Ensemble consensus clustering
├── python_scmega/             # Multi-omics integration and GRN inference
├── tf_driver/                 # NACR-related scripts
├── Maker_dataset/             # Marker-gene database
├── requirements.txt           # Python dependencies
└── requirements_R.txt         # Environment record
```

## Requirements

The main analysis was implemented using Python 3.10 and R 4.5.x.

```bash
pip install -r requirements.txt
Rscript r_requirements.R
```

An available LLM API is required for snTypeGPT annotation.

## Data

The associated study uses publicly available PBMC, human pancreas, and stem cell-derived islet datasets. Dataset information and preprocessing details are provided in the manuscript.
