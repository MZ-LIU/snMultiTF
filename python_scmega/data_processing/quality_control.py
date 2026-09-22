"""
Quality Control Module - MuData/AnnData format support

Support workflow:
1. Use muon to read 10X data: mdata = mu.io.read_10x_h5("data.h5")
2. Extract RNA and ATAC respectively: rna = mdata['rna'].copy(), atac = mdata['atac'].copy()
3. Separate quality control: calculate_rna_qc_metrics(), filter_rna_cells(), etc.
4. Synchronize cells: sync_common_cells()

RNA filtering standards:
- nFeature_RNA: 200-5000 (number of genes)
- percent.mt: ≤20% (mitochondrial gene proportion)

ATAC filtering standards:
- Only peaks with standard chromosome naming (chr1-chr22, chrX, chrY, chrM) are retained
- min.cells: filter low quality peaks
- Add gene annotation information to establish RNA-ATAC correspondence

Cell synchronization:
- Find public cells through cell barcodes to ensure that multi-omics data comes from the same cell
"""

import numpy as np
import pandas as pd
import anndata as ad
from typing import Optional, Dict, Any, Union, Tuple
from pathlib import Path
import warnings
import re
import scanpy as sc
from typing import Union, Optional
from muon import MuData

try:
    import muon as mu
    MUON_AVAILABLE = True
except ImportError:
    MUON_AVAILABLE = False
    warnings.warn('The muon library is not installed and cannot read 10X multiome data. Please install: pip install muon')



def _filter_standard_chromosomes(peak_names: list) -> np.ndarray:
    """
    Filter peaks with non-standard chromosome naming
    
    Standard format: chr1:10000-10500 or chr1-10000-10500
    Reserved: chr1-chr22, chrX, chrY, chrM
    """
    standard_chrs = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']
    
    def is_standard_peak(peak_name):
        match = re.match(r'(chr[^:\-]+)', peak_name)
        if match:
            chr_name = match.group(1)
            return chr_name in standard_chrs
        return False
    
    mask = np.array([is_standard_peak(peak) for peak in peak_names])
    return mask


def ensure_ensembl_installed(release=86, cache_dir=None):
    """Make sure Ensembl data is installed"""
    import pyensembl
    import os
    from pathlib import Path
    
    if cache_dir is None:
        cache_dir =  "E:/scMEGA9.25/python_scmega/data/.pyensembl"
    else:
        cache_dir = Path(cache_dir)
    
    os.environ['PYENSEMBL_CACHE_DIR'] = str(cache_dir)

    ensembl = pyensembl.EnsemblRelease(release)
    if not ensembl.required_local_files_exist():
        print(f"Ensembl {release} data is not installed, automatic download starts...")

        # from .install_ensembl import download_ensembl_annotation
        # download_ensembl_annotation(cache_dir=cache_dir, release=release)

    # Make sure the database index is built
    try:
        ensembl.genes()  # Test if available
    except (ValueError, Exception):
        print(f"Ensembl {release} index not built, start building...")
        ensembl.index()
    
    return ensembl

def _add_gene_annotation(atac_adata, genome="hg38"):
    """Add gene annotation information"""
    try:
        if genome == "hg38":
            # Make sure the data is installed
            ensembl = ensure_ensembl_installed(86, cache_dir =  "E:/scMEGA9.25/python_scmega/data/.pyensembl")
            
            genes = ensembl.genes()  # Returns a list of gene objects
            
            gene_df = pd.DataFrame([
                {
                    "gene_id": gene.id,
                    "chr": gene.contig,  # chromosome
                    "start": gene.start,
                    "end": gene.end,
                    "strand": gene.strand,
                    "gene_name": gene.name,
                    "gene_biotype": gene.biotype  # Gene type (e.g. protein coding, lncRNA)
                }
                for gene in genes
            ])
            
            atac_adata.uns["gene_annotation"] = gene_df
            print(f"   Annotation information for {len(gene_df)} genes has been added (Ensembl 86)")
            return atac_adata
        else:
            warnings.warn(f"The annotation of genome {genome} is not supported yet, skip it")
            return atac_adata
            
    except Exception as e:
        warnings.warn(f"Failed to add gene annotation: {str(e)}")
        return atac_adata



def calculate_rna_qc_metrics(adata: ad.AnnData, 
                             mt_pattern: str = "^MT-",
                             inplace: bool = True) -> Optional[ad.AnnData]:
    """
    Calculate quality control indicators for single RNA AnnData (used after extraction from MuData)
    
    Equivalent to R code:
    ```R
    obj.rna <- CreateSeuratObject(counts = rna_counts)
    obj.rna[["percent.mt"]] <- PercentageFeatureSet(obj.rna, pattern = "^MT-")
    ```
    
    Args:
        adata: RNA AnnData object
        mt_pattern: mitochondrial gene pattern
        inplace: whether to modify in place
        
    Returns:
        AnnData with QC indicators added (if inplace=False)
    """
    if inplace:
        target = adata
    else:
        target = adata.copy()

    
    target.var['mt'] = target.var_names.str.match(mt_pattern, case=False)
    sc.pp.calculate_qc_metrics(
        target, 
        qc_vars=['mt'],
        percent_top=None, 
        log1p=False, 
        inplace=True
    )
    # Add Seurat compatible column names
    target.obs['nFeature_RNA'] = target.obs['n_genes_by_counts']
    target.obs['nCount_RNA'] = target.obs['total_counts']
    target.obs['percent.mt'] = target.obs['pct_counts_mt']
    
    
    print(f"RNA QC indicator calculation completed:")
    print(f"  Cell number: {target.n_obs}")
    print(f"  Number of genes: {target.n_vars}")
    if 'nFeature_RNA' in target.obs.columns:
        print(f"  Average genes/cell: {target.obs['nFeature_RNA'].mean():.0f}")
    if 'nCount_RNA' in target.obs.columns:
        print(f"  Average UMI/cell: {target.obs['nCount_RNA'].mean():.0f}")
    if 'percent.mt' in target.obs.columns:
        print(f"  Average mitochondria %: {target.obs['percent.mt'].mean():.2f}%")
    
    if inplace:
        return None
    else:
        return target


def filter_rna_cells_features(
    adata,
    min_nFeature_RNA: int = 200,
    max_nFeature_RNA: int = 5000,
    min_nCount_RNA: int = 0,
    max_nCount_RNA: int = np.inf,
    max_percent_mt: float = 20.0,
    min_cells_per_gene: int = 0,
    exclude_gene_patterns: list = None,
    add_qc_flag: bool = False  # Optional: whether to retain QC marks
):
    """
    Filtering RNA cells and genes (revised version)
    """

    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    # ============================================================
    # ============================================================
    cell_filter = (
        (adata.obs['nFeature_RNA'] >= min_nFeature_RNA) &
        (adata.obs['nFeature_RNA'] <= max_nFeature_RNA) &
        (adata.obs['nCount_RNA'] >= min_nCount_RNA) &
        (adata.obs['nCount_RNA'] <= max_nCount_RNA) &
        (adata.obs['percent.mt'] <= max_percent_mt)
    )

    # Optional: Keep QC marks (on original data)
    if add_qc_flag:
        adata.obs['pass_rnaQC'] = cell_filter

    # ============================================================
    # ============================================================
    filtered = adata[cell_filter, :].copy()
    n_cells_after = filtered.n_obs

    print(f"[Cell filter]")
    print(f"  Condition: nFeature [{min_nFeature_RNA}, {max_nFeature_RNA}], "
          f"nCount [{min_nCount_RNA}, {max_nCount_RNA if max_nCount_RNA != np.inf else '∞'}], "
          f"percent.mt ≤ {max_percent_mt}%")
    print(f"  Cell: {n_cells_before} → {n_cells_after} (-{n_cells_before - n_cells_after})")

    # ============================================================
    # ============================================================
    n_genes_current = filtered.n_vars

    if min_cells_per_gene > 0:
        sc.pp.filter_genes(filtered, min_cells=min_cells_per_gene)
        print(f"[Gene filter] min_cells={min_cells_per_gene}: "
              f"{n_genes_current} → {filtered.n_vars}")
        n_genes_current = filtered.n_vars

    if exclude_gene_patterns:
        for pattern in exclude_gene_patterns:
            mask = ~filtered.var_names.str.contains(pattern, case=False, regex=True)
            filtered = filtered[:, mask]

        print(f"[Gene Pattern Exclusion] {exclude_gene_patterns}: "
              f"{n_genes_current} → {filtered.n_vars}")
        n_genes_current = filtered.n_vars

    print(f"[Final number of genes] {n_genes_before} → {n_genes_current}")

    # # ============================================================
    # # ============================================================
    # if 'pass_rnaQC' in filtered.obs.columns:
    #     n_before_marked = filtered.n_obs
    #     filtered = filtered[filtered.obs['pass_rnaQC'] == True].copy()
    #     n_after_marked = filtered.n_obs
    #     if n_before_marked > n_after_marked:
    
    return filtered


def filter_atac_data(adata: ad.AnnData,
                            min_cells: int = 1,
                            add_gene_annotation: bool = False,
                            genome: str = "hg38") -> ad.AnnData:
    """
    Filter ATAC data (used after extracting from MuData)
    According to the R version logic:
    1. Filter non-standard chromosome naming
    2. Filter low-quality peaks (min.cells)
    3. Add gene annotation
    
    Equivalent to R code:
    ```R
    atac_counts <- atac_counts[grep("chr", rownames(atac_counts)), ]
    chrom_assay <- CreateChromatinAssay(
        counts = atac_counts,
        min.cells = 1, # Do not filter peaks
        genome = 'hg38'
    )
    ```
    
     Args:
        adata: ATAC AnnData object
        min_cells: At least how many cells each peak is counted in (default 1, matches R's min.cells=1)
        add_gene_annotation: whether to add gene annotation
        genome: genome version
        
    Returns:
        Filtered AnnData
    """
    filtered = adata.copy()
    
    n_peaks_before = filtered.n_vars
    n_cells_before = filtered.n_obs

    
    print(f"Filter ATAC chromosomes...")
    peak_filter = _filter_standard_chromosomes(filtered.var_names.tolist())
    filtered = filtered[:, peak_filter]
    n_peaks_after_chr = filtered.n_vars
    print(f"  [1/3] Chromosome filter: {n_peaks_before} → {n_peaks_after_chr} peaks")
    
    if min_cells > 0:
        n_peaks_before_filter = filtered.n_vars
        sc.pp.filter_genes(filtered, min_cells=min_cells)
        n_peaks_after_filter = filtered.n_vars
        print(f"  [2/3] Peak filtering (min_cells={min_cells}): {n_peaks_before_filter} → {n_peaks_after_filter} peaks")
    # else:
    
    if add_gene_annotation:
        print(f"  [3/3] Add gene annotation...")
        filtered = _add_gene_annotation(filtered, genome=genome)
    else:
        print('  [3/3] Skip Python pyensembl gene annotation; R/EnsDb.Hsapiens.v86 annotation will be performed in 5.2_network.py')
    
    print(f"  Final result: {n_cells_before} cells × {filtered.n_vars} peaks")

    
    return filtered


def sync_common_cells(rna_adata: ad.AnnData, 
                     atac_adata: ad.AnnData) -> tuple:
    """
    Synchronize common cells for RNA and ATAC (used after extraction from MuData)
    
    Equivalent to R code:
    ```R
    cell.sel <- intersect(colnames(obj.rna), colnames(obj.atac))
    obj.rna <- subset(obj.rna, cells = cell.sel)
    obj.atac <- subset(obj.atac, cells = cell.sel)
    ```
    
    Args:
        rna_adata: RNA AnnData object
        atac_adata: ATAC AnnData object
        
    Returns:
        (rna_filtered, atac_filtered): tuples filtered to common cells
    """
    rna_cells = set(rna_adata.obs.index)
    atac_cells = set(atac_adata.obs.index)
    shared_cells = sorted(list(rna_cells.intersection(atac_cells)))
    
    if len(shared_cells) == 0:
        raise ValueError('Error: No public cells! Please check if the cell barcodes match.')
    
    print(f"\n cell synchronization:")
    print(f"  RNA cell number: {len(rna_cells)}")
    print(f"  ATAC cell number: {len(atac_cells)}")
    print(f"  Public cell number: {len(shared_cells)}")
    
    rna_filtered = rna_adata[shared_cells, :].copy()
    atac_filtered = atac_adata[shared_cells, :].copy()
    
    return rna_filtered, atac_filtered


def quality_control_separate(rna_adata: ad.AnnData,
                             atac_adata: ad.AnnData,
                             min_nFeature_RNA: int = 200,
                             max_nFeature_RNA: int = 5000,
                             min_nCount_RNA: int = 0,
                             max_nCount_RNA: int = np.inf,
                             max_percent_mt: float = 20.0,
                             min_cells_per_peak: int = 1,
                             min_cells_per_gene: int = 1,
                             exclude_gene_patterns: list = None,
                             mt_pattern: str = "^MT-",
                             add_gene_annotation: bool = True,
                             genome: str = "hg38",
                             return_mudata: bool = True) -> Union[MuData, tuple]:  
    """
    Complete separate quality control process (one-stop function extracted from MuData)
    
    Equivalent to the full R workflow:
    ```R
    #RNA QC - Filter cells and genes
    obj.rna <- CreateSeuratObject(counts = rna_counts)
    obj.rna[["percent.mt"]] <- PercentageFeatureSet(obj.rna, pattern = "^MT-")
    obj.rna <- subset(obj.rna, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
    obj.rna <- obj.rna[rowSums(obj.rna) >= min_counts, ] # Filter genes
    
    # ATAC QC - filter chromosomes, peaks and add annotations
    atac_counts <- atac_counts[grep("chr", rownames(atac_counts)), ]
    chrom_assay <- CreateChromatinAssay(counts = atac_counts, min.cells = 1, genome = 'hg38')
    obj.atac <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC")
    Annotation(obj.atac) <- annotations
    
    # synchronize cells
    cell.sel <- intersect(colnames(obj.rna), colnames(obj.atac))
    obj.rna <- subset(obj.rna, cells = cell.sel)
    obj.atac <- subset(obj.atac, cells = cell.sel)
    ```
    
    Args:
        rna_adata: RNA AnnData object
        atac_adata: ATAC AnnData object
        min_nFeature_RNA: Minimum number of RNA genes
        max_nFeature_RNA: Maximum number of RNA genes
        min_nCount_RNA: Minimum UMI count (default 0, no filtering)
        max_nCount_RNA: Maximum UMI count (default is infinity, no filtering)
        max_percent_mt: Maximum mitochondrial ratio
        min_cells_per_peak: Minimum number of cells per peak (default 1, matches R's min.cells=1)
        min_cells_per_gene: At least how many cells the gene is expressed in (default 0, no filtering)
        exclude_gene_patterns: list of excluded gene patterns, such as ["MT-", "RP"] (default None)
        mt_pattern: mitochondrial gene matching pattern
        add_gene_annotation: whether to add gene annotation
        genome: genome version
        return_mudata: whether to return MuData object (True) or tuple (False)
        
    Returns:
        MuData object or (rna_qc, atac_qc) tuple: data after quality control
        
    Example:
        ```python
        import muon as mu
        from python_scmega.data_processing import quality_control_separate
        
        #Read data
        mdata = mu.io.read_10x_h5("data.h5")
        rna_adata = mdata['rna'].copy()
        atac_adata = mdata['atac'].copy()
        
        #Basic quality control (cell filtration only)
        mdata_qc = quality_control_separate(
            rna_adata, atac_adata,
            min_nFeature_RNA=200,
            max_nFeature_RNA=5000,
            max_percent_mt=20,
            min_cells_per_peak=1
        )
        
        # Advanced quality control (including gene filtering and exclusion)
        mdata_qc = quality_control_separate(
            rna_adata, atac_adata,
            min_nFeature_RNA=200,
            max_nFeature_RNA=5000,
            min_nCount_RNA=800, # Minimum UMI count
            max_nCount_RNA=70000, # Maximum UMI count
            max_percent_mt=20,
            min_cells_per_peak=1,
            min_cells_per_gene=10, # The gene is expressed in at least 10 cells
            exclude_gene_patterns=["MT-", "RP"] # Exclude mitochondrial/ribosomal genes
        )
        ```
    """
    print("=" * 70)
    print('Separate quality control process (R style)')
    print("=" * 70)
    
    print('\n[Step 1/4] Calculate RNA quality control indicators...')
    calculate_rna_qc_metrics(rna_adata, mt_pattern=mt_pattern, inplace=True)
    
    print('\n[Step 2/4] Filter RNA cells and genes...')
    rna_filtered = filter_rna_cells_features(
        rna_adata,
        min_nFeature_RNA=min_nFeature_RNA,
        max_nFeature_RNA=max_nFeature_RNA,
        min_nCount_RNA=min_nCount_RNA,
        max_nCount_RNA=max_nCount_RNA,
        max_percent_mt=max_percent_mt,
        min_cells_per_gene=min_cells_per_gene,
        exclude_gene_patterns=exclude_gene_patterns
    )
    # rna_qc, atac_qc = sync_common_cells(rna_filtered, atac_filtered)

    
    print('\n[Step 3/4] Filter ATAC (chromosome + peaks + annotation)...')
    atac_filtered = filter_atac_data(
        atac_adata,
        min_cells=min_cells_per_peak,
        add_gene_annotation=add_gene_annotation,
        genome=genome
    )
    print('\n[Step 4/4] Synchronize public cells...')
    rna_qc, atac_qc = sync_common_cells(rna_filtered, atac_filtered)
    
    print("\n" + "=" * 70)
    print('Separate quality control completed!')
    print(f"  Final number of public cells: {rna_qc.n_obs}")
    print(f"  - RNA:  {rna_qc.n_vars} genes")
    print(f"  - ATAC: {atac_qc.n_vars} peaks")
    print("=" * 70)
    
    if return_mudata: 
        print('  Construct MuData object...')
        from muon import MuData
        mdata_qc = MuData({'rna': rna_qc, 'atac': atac_qc})
        
        # Keep ATAC fragments information (if it exists)
        if 'fragments' in atac_adata.uns:
            mdata_qc['atac'].uns['fragments'] = atac_adata.uns['fragments']
            print(f"   Fragments path information has been retained")
        print(' MuData object construction completed')
        return mdata_qc
    else:
        return rna_qc, atac_qc


def plot_qc_metrics(adata: ad.AnnData, save_path: Optional[str] = None):
    """
    Draw RNA quality control indicators
    
    Args:
        adata: RNA AnnData object
        save_path: save path (optional)
        
    Example:
        ```python
        import muon as mu
        mdata = mu.io.read_10x_h5("data.h5")
        rna_adata = mdata['rna'].copy()
        calculate_rna_qc_metrics(rna_adata, inplace=True)
        plot_qc_metrics_rna(rna_adata, save_path="qc_rna.png")
        ```
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
        matplotlib.rcParams['axes.unicode_minus'] = False
    except ImportError:
        warnings.warn('Requires matplotlib for plotting. Please install: pip install matplotlib')
        return

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    metrics = ['nFeature_RNA', 'nCount_RNA', 'percent.mt']
    labels = ['Number of genes (nFeature_RNA)', 'UMI Count (nCount_RNA)', 'Mitochondria % (percent.mt)']
    
    for ax, metric, label in zip(axes, metrics, labels):
        if metric in adata.obs.columns:
            values = adata.obs[metric].dropna()
            parts = ax.violinplot([values], positions=[0], widths=0.7, 
                                 showmeans=True, showmedians=True)
            for pc in parts['bodies']:
                pc.set_facecolor('#1f77b4')
                pc.set_alpha(0.7)
            ax.set_ylabel(label, fontsize=12)
            ax.set_title(label, fontsize=13, fontweight='bold')
            ax.set_xticks([])
            ax.grid(True, alpha=0.3, axis='y')
    
    plt.suptitle('RNA quality control indicators', fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Picture saved: {save_path}")
    
    plt.show()
    
def get_qc_summary(rna_adata: Optional[ad.AnnData] = None,
                  atac_adata: Optional[ad.AnnData] = None) -> Dict[str, Any]:
    """
    Get QC statistical summary of RNA and/or ATAC (unified function)
    
    Args:
        rna_adata: RNA AnnData object (optional)
        atac_adata: ATAC AnnData object (optional)
        
    Returns:
        Dictionary containing RNA and/or ATAC statistics
        
    Example:
        ```python
        # RNA summary alone
        summary = get_qc_summary(rna_adata=rna_qc)
        print(f"RNA cell number: {summary['rna']['n_cells']}")
        
        # Separate ATAC summary
        summary = get_qc_summary(atac_adata=atac_qc)
        
        #RNA and ATAC full summary
        rna_qc, atac_qc = quality_control_separate(rna, atac)
        summary = get_qc_summary(rna_qc, atac_qc)
        print(f"RNA cell number: {summary['rna']['n_cells']}")
        print(f"Number of ATAC peaks: {summary['atac']['n_peaks']}")
        ```
    """
    summary = {}
    
    if rna_adata is not None:
        rna_summary = {
            'n_cells': rna_adata.n_obs,
            'n_genes': rna_adata.n_vars,
        }
        
        if 'nFeature_RNA' in rna_adata.obs.columns:
            rna_summary['features_per_cell'] = {
                'mean': float(rna_adata.obs['nFeature_RNA'].mean()),
                'median': float(rna_adata.obs['nFeature_RNA'].median()),
                'std': float(rna_adata.obs['nFeature_RNA'].std())
            }
        
        if 'nCount_RNA' in rna_adata.obs.columns:
            rna_summary['counts_per_cell'] = {
                'mean': float(rna_adata.obs['nCount_RNA'].mean()),
                'median': float(rna_adata.obs['nCount_RNA'].median()),
                'std': float(rna_adata.obs['nCount_RNA'].std())
            }
        
        if 'percent.mt' in rna_adata.obs.columns:
            rna_summary['percent_mito'] = {
                'mean': float(rna_adata.obs['percent.mt'].mean()),
                'median': float(rna_adata.obs['percent.mt'].median()),
                'std': float(rna_adata.obs['percent.mt'].std())
            }
        
        summary['rna'] = rna_summary
    
    if atac_adata is not None:
        atac_summary = {
            'n_cells': atac_adata.n_obs,
            'n_peaks': atac_adata.n_vars
        }
        
        if 'nFeature_ATAC' in atac_adata.obs.columns:
            atac_summary['features_per_cell'] = {
                'mean': float(atac_adata.obs['nFeature_ATAC'].mean()),
                'median': float(atac_adata.obs['nFeature_ATAC'].median()),
                'std': float(atac_adata.obs['nFeature_ATAC'].std())
            }
        
        if 'nCount_ATAC' in atac_adata.obs.columns:
            atac_summary['counts_per_cell'] = {
                'mean': float(atac_adata.obs['nCount_ATAC'].mean()),
                'median': float(atac_adata.obs['nCount_ATAC'].median()),
                'std': float(atac_adata.obs['nCount_ATAC'].std())
            }
        
        summary['atac'] = atac_summary
    
    return summary
