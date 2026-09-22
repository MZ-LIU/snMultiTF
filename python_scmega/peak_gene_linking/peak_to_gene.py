"""
Peak-to-Gene Linking Module

Implements peak-to-gene linking algorithms equivalent to R scMEGA's PeakToGene function.
This module provides trajectory-based peak-to-gene linking for GRN inference.

Key functions:
- link_peaks_to_genes: Main peak-to-gene linking function (R PeakToGene equivalent)
- load_gene_annotation: Load real gene annotations from .rda files
- parse_peak_annotation: Parse peak coordinates from peak names (chr-start-end format)

Distance Calculation (R-compatible):
    distance = abs(peak_center - gene_TSS)
    
Where:
    - peak_center = (peak_start + peak_end) / 2
    - gene_TSS = gene_start (for + strand genes) or gene_end (for - strand genes)
    - Only peak-gene pairs on the same chromosome are considered
    - Default maximum distance: 250kb (250,000 bp)

References:
- Original R code: PeakToGene function in scMEGA R/peak_to_gene.R
- ArchR: Granja et al. (2021) ArchR is a scalable software package for integrative analysis
"""

import numpy as np
import pandas as pd
from typing import Optional
import warnings
#from scipy.stats import pearsonr
from scipy.stats import t as t_dist  
from pathlib import Path

from functools import lru_cache

# Import rpy2 for loading .rda files (REQUIRED)
try:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    R_AVAILABLE = True
except ImportError:
    raise ImportError(
        "rpy2 is required for loading gene annotations from .rda files. "
        "Please install: pip install rpy2>=3.6.2"
    )


def load_gene_annotation(genome: str = "hg38") -> pd.DataFrame:
    """
    Load real gene annotations from .rda files (R-compatible).
    
    This function loads the same gene annotation data used in R scMEGA,
    matching the R implementation in peak_to_gene.R lines 38-46.
    
    Args:
        genome: Reference genome ("hg38", "hg19", "mm10", "mm9")
        
    Returns:
        DataFrame with columns: chr, start, end, strand, symbol
        Index: gene symbols
        
    Raises:
        ValueError: If genome is not supported
        FileNotFoundError: If annotation file is not found
        RuntimeError: If failed to load annotation data
        
    Example:
        ```python
        gene_anno = load_gene_annotation(genome="hg38")
        # Returns real gene annotations from geneAnnoHg38.rda
        ```
    """
    if genome not in ["hg19", "hg38", "mm9", "mm10"]:
        raise ValueError(f"Available genomes are: hg19, hg38, mm9, mm10. Got: {genome}")
    
    # Get path to annotation file
    data_dir = Path(__file__).parent.parent / "data"
    anno_file = data_dir / f"geneAnno{genome.capitalize()}.rda"
    
    if not anno_file.exists():
        raise FileNotFoundError(
            f"Gene annotation file not found: {anno_file}\n"
            f"Expected location: python_scmega/data/geneAnno{genome.capitalize()}.rda"
        )
    
    try:
        # Load .rda file (use as_posix() for Windows path compatibility)
        r_path = anno_file.as_posix()  # Convert to forward slash format
        print(f"  Loading gene annotation from: {r_path}")  # Add debugging information
        ro.r(f'load("{r_path}")')
        
        # Load the GenomicRanges library (required! for seqnames, start, end, strand functions)
        ro.r('library(GenomicRanges)')
        
        # Get the gene annotation object (R: geneAnnoHg38$genes, etc.)
        gene_anno_name = f"geneAnno{genome.capitalize()}"
        ro.r(f'genes <- {gene_anno_name}$genes')
        
        # Extract gene information (matching R lines 49-58)
        ro.r('''
        gene_df <- data.frame(
            symbol = genes$symbol,
            chr = as.character(seqnames(genes)),
            start = start(genes),
            end = end(genes),
            strand = as.character(strand(genes)),
            stringsAsFactors = FALSE
        )
        ''')
        
        # Convert to pandas DataFrame
        with localconverter(ro.default_converter + pandas2ri.converter):
            gene_df = ro.r('gene_df')

        
        # Remove duplicate gene symbols (keep the first occurrence)
        # This is necessary because some genes may have multiple transcripts or locations
        if gene_df['symbol'].duplicated().any():
            n_duplicates = gene_df['symbol'].duplicated().sum()
            print(f"  Warning: Found {n_duplicates} duplicate gene symbols, keeping first occurrence")
            gene_df = gene_df.drop_duplicates(subset='symbol', keep='first')
        
        
        # Set index to gene symbol
        gene_df = gene_df.set_index('symbol')
        
        print(f"  Loaded {len(gene_df)} genes from {genome} annotation")
        
        return gene_df
        
    except Exception as e:
        raise RuntimeError(f"Failed to load gene annotation from {anno_file}: {str(e)}")


def parse_peak_annotation(peak_names: pd.Index) -> pd.DataFrame:
    """
    Parse peak coordinates from peak names.
    
    Two formats are supported:
    - chr:start-end (e.g. "chr1:1000-2000") - 10X Genomics standard format
    - chr-start-end (e.g. "chr1-1000-2000") - alternative format
    
    Args:
        peak_names: index of peak identifier
        
    Returns:
        DataFrame containing columns: chr, start, end
        Index: peak name
        
    Raises:
        ValueError: if the peak name is not in the correct format
    """
    annotations = []
    
    for peak_name in peak_names:
        try:
            peak_str = str(peak_name)
            
            if ':' in peak_str and '-' in peak_str:
                chr_name, coords = peak_str.split(':', 1)
                start_str, end_str = coords.split('-', 1)
            elif peak_str.count('-') == 2:
                parts = peak_str.split('-')
                if len(parts) == 3:
                    chr_name, start_str, end_str = parts
                else:
                    raise ValueError(f"Peak name '{peak_name}' format is abnormal")
            else:
                raise ValueError(
                    f"Peak name '{peak_name}' does not conform to a known format."
                    f"Expected format 'chr:start-end' or 'chr-start-end'"
                )
            
            if not chr_name.startswith('chr'):
                raise ValueError(
                    f"The chromosome name '{chr_name}' in peak '{peak_name}' is invalid."
                    f"Expected format: chr1, chr2, chrX, etc."
                )
            
            try:
                start = int(start_str)
                end = int(end_str)
            except ValueError:
                raise ValueError(
                    f"The coordinates in Peak '{peak_name}' are invalid."
                    f"The start and end positions must be integers."
                )
            
            if start >= end:
                raise ValueError(
                    f"The coordinates in Peak '{peak_name}' are invalid."
                    f"The start position ({start}) must be smaller than the end position ({end})."
                )
            if start < 0 or end < 0:
                raise ValueError(
                    f"The coordinates in Peak '{peak_name}' are invalid."
                    f"Coordinates must be non-negative."
                )
            
            annotations.append({
                'chr': chr_name,
                'start': start,
                'end': end
            })
            
        except Exception as e:
            raise ValueError(
                f"Parsing peak '{peak_name}' failed: {str(e)}\n"
                f"Peak names must conform to the following format: 'chr:start-end' (e.g. 'chr1:1000-2000') "
                f"or 'chr-start-end' (e.g. 'chr1-1000-2000')"
            )
    
    df = pd.DataFrame(annotations, index=peak_names)
    print(f"  {len(df)} peaks have been resolved, spanning {df['chr'].nunique()} chromosomes")
    
    return df

            
@lru_cache(maxsize=4)
def _cached_gene_annotation(genome: str) -> pd.DataFrame:
    """LRU caches gene annotations to avoid repeated rpy2 loading."""
    return load_gene_annotation(genome)


def _enumerate_pairs_fast(
    peak_annotation: pd.DataFrame,
    gene_annotation: pd.DataFrame,
    peak_mat: pd.DataFrame,
    gene_mat_filtered: pd.DataFrame,
    max_distance: int
):
    """
    Quickly enumerate candidate pairs according to the chromosome+TSS distance window to avoid O(P×G) full scan.
    """
    pairs_info, peak_indices, gene_indices = [], [], []
    peak_name_to_idx = {name: idx for idx, name in enumerate(peak_mat.index)}
    gene_name_to_idx = {name: idx for idx, name in enumerate(gene_mat_filtered.index)}

    peaks_by_chr = {}
    for chr_name, dfp in peak_annotation.groupby('chr'):
        centers = ((dfp['start'].values + dfp['end'].values) * 0.5).astype(np.float64)
        order = np.argsort(centers)
        peaks_by_chr[chr_name] = (centers[order], dfp.index.values[order])

    for chr_name, dfg in gene_annotation.groupby('chr'):
        if chr_name not in peaks_by_chr:
            continue
        centers_sorted, peak_names_sorted = peaks_by_chr[chr_name]
        if centers_sorted.size == 0:
            continue

        for gene_name, row in dfg.iterrows():
            if gene_name not in gene_name_to_idx:
                continue  # Filtered genes only
            tss = row['start'] if row['strand'] == '+' else row['end']
            lo, hi = tss - max_distance, tss + max_distance

            L = np.searchsorted(centers_sorted, lo, side='left')
            R = np.searchsorted(centers_sorted, hi, side='right')
            if L >= R:
                continue

            centers_sel = centers_sorted[L:R]
            peak_names_sel = peak_names_sorted[L:R]
            distances = np.abs(centers_sel - tss)

            for pk, dist, pc in zip(peak_names_sel, distances, centers_sel):
                pairs_info.append({
                    'peak': pk,
                    'gene': gene_name,
                    'chr': chr_name,
                    'distance': float(dist),
                    'peak_center': float(pc),
                    'gene_tss': float(tss)
                })
                peak_indices.append(peak_name_to_idx[pk])
                gene_indices.append(gene_name_to_idx[gene_name])

    return pairs_info, np.array(peak_indices, dtype=np.int64), np.array(gene_indices, dtype=np.int64)        

def link_peaks_to_genes(
    peak_mat: pd.DataFrame,
    gene_mat: pd.DataFrame,
    genome: str = "hg38",
    peak_annotation: Optional[pd.DataFrame] = None,
    gene_annotation: Optional[pd.DataFrame] = None,
    max_distance: int = 250000
) -> pd.DataFrame:
    """
    Link peaks to genes using trajectory matrices (R PeakToGene equivalent).
    
    This is the R-compatible version that works on trajectory matrices
    rather than single-cell data. It implements the exact logic from R's
    PeakToGene function:
    
    ```R
    df.p2g <- PeakToGene(
        peak.mat = groupMatATAC,
        gene.mat = groupMatRNA,
        genome = genome
    )
    ```
    
    Algorithm (matching R implementation in R/peak_to_gene.R):
    1. Load real gene annotations from .rda files (R lines 38-46)
    2. Calculate TSS for each gene (considering strand direction)
    3. Find peak-gene pairs within max_distance
    4. Calculate Pearson correlation across trajectory time points
    5. Compute t-statistic and p-value
    6. Apply FDR correction (Benjamini-Hochberg)
    """
    print(f"Linking peaks to genes from trajectory matrices (R PeakToGene equivalent)...")
    print(f"  Peak matrix: {peak_mat.shape} (peaks × time_bins)")
    print(f"  Gene matrix: {gene_mat.shape} (genes × time_bins)")
    
    # Load genomic annotations
    if peak_annotation is None:
        peak_annotation = parse_peak_annotation(peak_mat.index)
    
    # Unified handling of gene annotation and matrix filtering (whether provided externally or not)
    if gene_annotation is None:
        print("DEBUG: About to load gene annotation...")
        gene_annotation = _cached_gene_annotation(genome)
        print("DEBUG: Gene annotation loaded successfully")
    
    # Only keep the genes that overlap with each other
    gene_use = list(set(gene_mat.index).intersection(set(gene_annotation.index)))
    if len(gene_use) == 0:
        raise ValueError(
            f"No genes in gene_mat found in {genome} annotations. "
            f"Gene matrix has {len(gene_mat)} genes, but none match the annotation. "
            f"Please check that gene names are in the correct format (gene symbols)."
        )
    gene_annotation = gene_annotation.loc[gene_use]
    gene_mat_filtered = gene_mat.loc[gene_use].copy()

    print(f"  Using {len(gene_annotation)} genes with annotations")
    print(f"  Filtered gene_mat to {len(gene_mat_filtered)} genes")
    
    print(f"  Identifying peak-gene pairs within {max_distance}bp...")
    pairs_info, peak_indices_arr, gene_indices_arr = _enumerate_pairs_fast(
        peak_annotation=peak_annotation,
        gene_annotation=gene_annotation,
        peak_mat=peak_mat,
        gene_mat_filtered=gene_mat_filtered,
        max_distance=max_distance
    )
    
    if len(pairs_info) == 0:
        warnings.warn("No peak-gene pairs found within distance threshold")
        return pd.DataFrame(columns=['peak', 'gene', 'chr', 'distance', 'peak_center', 
                                     'gene_tss', 'Correlation', 'TStat', 'pvalue', 'FDR'])
    print(f"  Found {len(pairs_info)} peak-gene pairs to test")
    
    # print(f"  Calculating correlations for {len(pairs_info)} pairs (vectorized GEMM)...")
    # used_peak_idx = np.unique(peak_indices_arr)
    # used_gene_idx = np.unique(gene_indices_arr)

    # peak_map = -np.ones(len(peak_mat), dtype=np.int64)
    # gene_map = -np.ones(len(gene_mat_filtered), dtype=np.int64)
    # peak_map[used_peak_idx] = np.arange(used_peak_idx.size, dtype=np.int64)
    # gene_map[used_gene_idx] = np.arange(used_gene_idx.size, dtype=np.int64)
    # peak_idx_comp = peak_map[peak_indices_arr]
    # gene_idx_comp = gene_map[gene_indices_arr]

    # X = peak_mat.values[used_peak_idx, :]
    # Y = gene_mat_filtered.values[used_gene_idx, :]

    # Xc = X - X.mean(axis=1, keepdims=True)
    # Yc = Y - Y.mean(axis=1, keepdims=True)

    # eps = 1e-12
    # sX = np.sqrt((Xc ** 2).sum(axis=1, keepdims=True))
    # sY = np.sqrt((Yc ** 2).sum(axis=1, keepdims=True))
    # sX = np.maximum(sX, eps)
    # sY = np.maximum(sY, eps)


    # correlations = Corr_full[peak_idx_comp, gene_idx_comp]
    # print(f"   Correlations calculated")
    # Step 2: Correlation (calculate in batches to avoid memory overflow)
    print(f"  Calculating correlations for {len(pairs_info)} pairs (batched)...")
    used_peak_idx = np.unique(peak_indices_arr)
    used_gene_idx = np.unique(gene_indices_arr)

    peak_map = -np.ones(len(peak_mat), dtype=np.int64)
    gene_map = -np.ones(len(gene_mat_filtered), dtype=np.int64)
    peak_map[used_peak_idx] = np.arange(used_peak_idx.size, dtype=np.int64)
    gene_map[used_gene_idx] = np.arange(used_gene_idx.size, dtype=np.int64)
    peak_idx_comp = peak_map[peak_indices_arr]
    gene_idx_comp = gene_map[gene_indices_arr]

    X = peak_mat.values[used_peak_idx, :]
    Y = gene_mat_filtered.values[used_gene_idx, :]

    Xc = X - X.mean(axis=1, keepdims=True)
    Yc = Y - Y.mean(axis=1, keepdims=True)
    eps = 1e-12
    sX = np.sqrt((Xc ** 2).sum(axis=1, keepdims=True))
    sY = np.sqrt((Yc ** 2).sum(axis=1, keepdims=True))
    sX = np.maximum(sX, eps)
    sY = np.maximum(sY, eps)

    Xc_norm = Xc / sX
    Yc_norm = Yc / sY

    # Calculate correlation coefficients in batches (avoid creating full matrix)
    batch_size = 10000  # Process 10,000 peaks per batch
    n_peaks = len(used_peak_idx)
    correlations = np.zeros(len(pairs_info), dtype=np.float64)

    print(f"  Processing in batches of {batch_size} peaks...")
    for batch_start in range(0, n_peaks, batch_size):
        batch_end = min(batch_start + batch_size, n_peaks)
        
        Corr_batch = Xc_norm[batch_start:batch_end] @ Yc_norm.T
        
        batch_mask = (peak_idx_comp >= batch_start) & (peak_idx_comp < batch_end)
        batch_pairs = np.where(batch_mask)[0]
        
        for pair_idx in batch_pairs:
            local_peak_idx = peak_idx_comp[pair_idx] - batch_start
            gene_idx = gene_idx_comp[pair_idx]
            correlations[pair_idx] = Corr_batch[local_peak_idx, gene_idx]
        
        if (batch_start // batch_size + 1) % 5 == 0:
            print(f"    Processed {batch_end}/{n_peaks} peaks...")

    print(f"   Correlations calculated")
    print(f"  Computing statistics...")
    n_timepoints = peak_mat.shape[1]
    results = []
    
    for i, pair_info in enumerate(pairs_info):
        corr = correlations[i]
        if not np.isnan(corr):
            # t-stat
            denom = np.sqrt(max(1 - corr**2, 1e-16) / (n_timepoints - 2))
            t_stat = corr / denom
            pval = 2 * t_dist.cdf(-abs(t_stat), df=n_timepoints-2)
            results.append({
                **pair_info,
                'Correlation': corr,
                'TStat': t_stat,
                'pvalue': pval
            })
    
    if len(results) == 0:
        warnings.warn("No valid peak-gene correlations found")
        return pd.DataFrame(columns=['peak', 'gene', 'chr', 'distance', 'peak_center', 
                                     'gene_tss', 'Correlation', 'TStat', 'pvalue', 'FDR'])
    
    df_p2g = pd.DataFrame(results)
    
    # FDR
    try:
        from statsmodels.stats.multitest import multipletests
        _, fdr_values, _, _ = multipletests(df_p2g['pvalue'], method='fdr_bh')
        df_p2g['FDR'] = fdr_values
    except ImportError:
        warnings.warn("statsmodels not available, FDR will be same as p-value")
        df_p2g['FDR'] = df_p2g['pvalue']
    
    df_p2g = df_p2g[~df_p2g['FDR'].isna()]
    df_p2g = df_p2g.sort_values('Correlation', ascending=False)
    
    print(f"   Found {len(df_p2g)} peak-gene links within {max_distance}bp")
    return df_p2g

