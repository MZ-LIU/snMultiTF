# snMultiTF

snMultiTF is an integrated framework for single-nucleus multi-omics analysis, combining ensemble consensus clustering, LLM-assisted cell-type annotation, trajectory-associated gene regulatory network (GRN) inference, and regulatory transcription factor (TF) identification.

The `workflow_unlabeled` pipeline is designed for datasets without predefined cell-type labels, where cell clusters are identified by ensemble consensus clustering and annotated with snTypeGPT. The `workflow_labeled` pipeline is designed for datasets with predefined cell-type labels, which are incorporated directly into the downstream workflow for multi-omics integration, trajectory-associated GRN inference, and regulatory TF identification.
