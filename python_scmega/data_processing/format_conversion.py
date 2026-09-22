"""
Format Conversion Module

Provides conversion functions between Python AnnData and R Seurat formats,
with special support for ATAC data with fragment association.

Key functions:
- _atac_anndata_to_seurat_r: Convert ATAC AnnData to R ChromatinAssay with fragments
- _diagnose_fragment_compatibility: Diagnose fragment file compatibility
- _diagnose_failure_reason: Analyze R error reasons
- _suggest_solutions: Provide troubleshooting suggestions

References:
- Signac CreateChromatinAssay for fragment association
- Seurat multiome object structure
"""

import anndata as ad
from typing import Optional, Dict, Any
import warnings
from pathlib import Path
from scipy import sparse


# R/Python interoperability
try:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri, numpy2ri
    from rpy2.robjects.conversion import localconverter
    R_AVAILABLE = True
except ImportError:
    R_AVAILABLE = False
    warnings.warn(
        "rpy2 is not available. R conversion functions will not work. "
        "Install rpy2 with: pip install rpy2"
    )


def _atac_anndata_to_seurat_r(atac_adata: ad.AnnData,
                              fragments_path: Optional[str] = None,
                              genome: str = "hg38",
                              verbose: bool = True) -> bool:
    """
    Convert ATAC AnnData to R Seurat ChromatinAssay WITH fragment association.
    
    This function performs the complete workflow:
    1. Transfer ATAC counts matrix to R
    2. Create ChromatinAssay with fragment file association
    3. Create Seurat object with the ChromatinAssay
    
    This is NOT just format conversion - it actually associates the fragments file.
    
    Args:
        atac_adata: AnnData object containing ATAC data
        fragments_path: Path to fragments.tsv.gz file (required for association)
        genome: Genome assembly (e.g., "hg38", "mm10")
        verbose: Whether to print progress messages
        
    Returns:
        True if successful, False otherwise
        
    Raises:
        ImportError: If rpy2 not available
        ValueError: If fragments_path not provided
        FileNotFoundError: If fragments file doesn't exist
    """
    if not R_AVAILABLE:
        raise ImportError("rpy2 is required for R conversion and fragment association")
    
    if fragments_path is None:
        raise ValueError("fragments_path is required for fragment association")
    
    fragments_path = Path(fragments_path)
    if not fragments_path.exists():
        raise FileNotFoundError(f"Fragment file not found: {fragments_path}")
    
    # Configure rpy2 to handle non-UTF-8 console output (e.g., Chinese error messages)
    from rpy2.rinterface_lib import callbacks
    import sys
    
    def console_write_ex(s):
        """Handle R console output with encoding error tolerance"""
        try:
            # Try UTF-8 first
            if isinstance(s, bytes):
                s = s.decode('utf-8')
            sys.stdout.write(s)
        except UnicodeDecodeError:
            # UTF-8 failed, try common Windows encodings
            if isinstance(s, bytes):
                for encoding in ['gbk', 'gb2312', 'latin1', 'cp1252']:
                    try:
                        s = s.decode(encoding)
                        break
                    except (UnicodeDecodeError, LookupError):
                        continue
                else:
                    # All encodings failed, use replace strategy
                    s = s.decode('utf-8', errors='replace')  # Replace bad bytes with 
            sys.stdout.write(s)
        sys.stdout.flush()
    
    # Set the custom console callbacks for both normal output and warnings/errors
    original_console_write = callbacks.consolewrite_print
    original_console_warnerror = callbacks.consolewrite_warnerror
    callbacks.consolewrite_print = console_write_ex
    callbacks.consolewrite_warnerror = console_write_ex  # Critical: handle R warnings/errors
    
    if verbose:
        print(f"Converting ATAC data to R ChromatinAssay with fragment association...")
        print(f"  ATAC: {atac_adata.n_obs} cells × {atac_adata.n_vars} peaks")
        print(f"  Fragment file: {fragments_path.name}")
        print(f"  Genome: {genome}")
    
    # Pre-diagnose fragment compatibility
    diagnostic_result = _diagnose_fragment_compatibility(
        atac_adata, 
        str(fragments_path), 
        verbose=verbose
    )
    
    try:
        # Load R libraries
        ro.r('suppressMessages(library(Seurat))')
        ro.r('suppressMessages(library(Signac))')
        ro.r('suppressMessages(library(Matrix))')

        # Prefer raw counts (emulates R’s ChromatinAssay@counts slot)
        if 'counts' in atac_adata.layers:
            if verbose:
                print("  Using raw counts from layers['counts']")
            counts_data = atac_adata.layers['counts']
        elif atac_adata.raw is not None:
            if verbose:
                print("  Using raw counts from .raw.X")
            counts_data = atac_adata.raw.X
        else:
            if verbose:
                print("  WARNING: No raw counts found in layers['counts'] or .raw")
                print("  Using .X directly (may be normalized TF-IDF data!)")
                print("  This may cause incorrect results in AddMotifs/RunChromVAR")
            counts_data = atac_adata.X
        
        if sparse.issparse(counts_data):
            counts_matrix = counts_data.T
            if not isinstance(counts_matrix, sparse.csc_matrix):
                counts_matrix = counts_matrix.tocsc()
        else:
            counts_matrix = sparse.csc_matrix(counts_data.T)

        if counts_matrix.shape[0] == 0 or counts_matrix.shape[1] == 0:
            raise ValueError(f"Invalid matrix dimensions: {counts_matrix.shape}")

        if counts_matrix.nnz == 0:
            raise ValueError('The matrix has no non-zero elements')

        if verbose:
            print(f"  Matrix verification passed: {counts_matrix.shape}, non-zero elements: {counts_matrix.nnz}")

        # Transfer sparse matrix to R using Matrix package
        # Extract sparse matrix components (CSC format matches R's dgCMatrix)
        data = counts_matrix.data
        indices = counts_matrix.indices + 1  # R uses 1-based indexing
        indptr = counts_matrix.indptr + 1
        shape = counts_matrix.shape

        # Prepare path string with forward slashes (R supports this on Windows)
        fragments_path_r = str(fragments_path.absolute()).replace('\\', '/')

        with localconverter(ro.default_converter + numpy2ri.converter + pandas2ri.converter):
            ro.globalenv['sparse_data'] = data
            ro.globalenv['sparse_i'] = indices
            ro.globalenv['sparse_p'] = indptr
            ro.globalenv['nrow'] = shape[0]
            ro.globalenv['ncol'] = shape[1]
            ro.globalenv['peak_names'] = atac_adata.var.index.tolist()
            ro.globalenv['cell_names'] = atac_adata.obs.index.tolist()
            
            # Transfer cell metadata to preserve data consistency
            if not atac_adata.obs.empty and len(atac_adata.obs.columns) > 0:
                ro.globalenv['cell_metadata'] = atac_adata.obs
                ro.globalenv['has_cell_metadata'] = True
            else:
                ro.globalenv['has_cell_metadata'] = False

        if verbose:
            print(f"  Creating sparse dgCMatrix...")
            sparsity = 100 * (1 - len(data) / (shape[0] * shape[1]))
            print(f"  Sparsity: {sparsity:.2f}% ({len(data):,} non-zero elements)")
            print(f"  Fragment path for R: {fragments_path_r}")
            if ro.globalenv['has_cell_metadata'][0]:
                print(f"  Cell metadata: {len(atac_adata.obs.columns)} columns")

        # Create ChromatinAssay WITH fragment association.
        # Match the scMEGA PBMC R vignette:
        # CreateChromatinAssay(..., sep = c(":", "-"), min.cells = 1,
        #                       genome = "hg38", fragments = fragments_path)
        r_code = f'''
        library(Signac)
        library(Seurat)
        library(GenomicRanges)
        library(GenomeInfoDb)

        # Create dgCMatrix from sparse components
        atac_counts <- new("dgCMatrix",
                        i = as.integer(sparse_i - 1),
                        p = as.integer(sparse_p - 1), 
                        x = as.numeric(sparse_data),
                        Dim = as.integer(c(nrow, ncol)))

        rownames(atac_counts) <- peak_names
        colnames(atac_counts) <- cell_names

        # Initialize variables
        association_success <- FALSE
        creation_error <- NULL

        # Create ChromatinAssay with error handling
        tryCatch({{
            chromatin_assay <- CreateChromatinAssay(
                counts = atac_counts,
                sep = c(":", "-"),
                min.cells = 1,
                genome = "{genome}",
                fragments = "{fragments_path_r}"
            )
            association_success <- TRUE
        }}, error = function(e) {{
            creation_error <<- conditionMessage(e)
            association_success <<- FALSE
        }})

         # Create Seurat object only if ChromatinAssay succeeded
        if (association_success) {{
            seurat_obj <- CreateSeuratObject(
                counts = chromatin_assay,
                assay = "ATAC"
            )
            
            # Transfer cell metadata to preserve data consistency
            # FIX: Add columns individually to avoid dimension mismatch
            if (exists("has_cell_metadata") && has_cell_metadata && exists("cell_metadata")) {{
                shared_cells <- intersect(colnames(seurat_obj), rownames(cell_metadata))
                if(length(shared_cells) > 0) {{
                    # Add each column individually instead of wholesale replacement
                    for (col_name in colnames(cell_metadata)) {{
                        # Convert to character to avoid factor level issues
                        seurat_obj@meta.data[shared_cells, col_name] <- as.character(cell_metadata[shared_cells, col_name])
                    }}
                }}
            }}

            # Match the scMEGA PBMC vignette gene annotation source.
            # annotations <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)
            # seqlevelsStyle(annotations) <- "UCSC"
            # Annotation(obj.atac) <- annotations
            annotation_success <- FALSE
            annotation_error <- NULL
            if ("{genome}" == "hg38") {{
                tryCatch({{
                    suppressMessages(library(EnsDb.Hsapiens.v86))
                    annotations <- GetGRangesFromEnsDb(
                        ensdb = EnsDb.Hsapiens.v86,
                        verbose = FALSE
                    )
                    seqlevelsStyle(annotations) <- "UCSC"
                    Annotation(seurat_obj[["ATAC"]]) <- annotations
                    annotation_success <- TRUE
                }}, error = function(e) {{
                    annotation_error <<- conditionMessage(e)
                    annotation_success <<- FALSE
                }})
            }}
            
            # Verify fragment association
            has_fragments <- !is.null(Fragments(seurat_obj[["ATAC"]]))
        }} else {{
            has_fragments <- FALSE
            annotation_success <- FALSE
            annotation_error <- NULL
        }}
        # # Restore output
        # sink()
        # options(warn = 0)
        '''
        
        try:
            ro.r(r_code)
            
            try:
                association_success = ro.globalenv['association_success'][0]
                has_fragments = ro.globalenv['has_fragments'][0]
                annotation_success = bool(ro.globalenv['annotation_success'][0])
                
                creation_error_r = ro.globalenv['creation_error']
                if creation_error_r == ro.NULL or creation_error_r is ro.NULL:
                    creation_error = None
                else:
                    creation_error = str(creation_error_r[0])

                annotation_error_r = ro.globalenv['annotation_error']
                if annotation_error_r == ro.NULL or annotation_error_r is ro.NULL:
                    annotation_error = None
                else:
                    annotation_error = str(annotation_error_r[0])
                
            except Exception as e:
                if verbose:
                    print(f"  [FAIL] Failed to retrieve R results: {str(e)}")
                association_success = False
                has_fragments = False
                annotation_success = False
                creation_error = str(e)
                annotation_error = str(e)
            
            if verbose:
                if association_success and has_fragments:
                    print(f"  [OK] Fragment association successful")
                    if annotation_success:
                        print(f"  [OK] Gene annotation added from EnsDb.Hsapiens.v86")
                    elif annotation_error:
                        print(f"  [WARN] Gene annotation failed: {annotation_error}")
                elif creation_error:
                    print(f"  [FAIL] ChromatinAssay creation failed: {creation_error}")
                    _diagnose_failure_reason(creation_error, diagnostic_result)
                else:
                    print(f"  [WARN] ChromatinAssay created but fragments not associated")
                    _suggest_solutions(diagnostic_result)
            
            return association_success and has_fragments

        except Exception as e:
            if verbose:
                print(f"  [FAIL] R execution failed: {str(e)}")
                # Provide direct diagnostic info based on exception type
                if 'UnicodeDecodeError' in str(type(e).__name__):
                    print(f"    - Encoding error in R console output")
                    print(f"    - R may have output non-UTF-8 characters")
                    print(f"    - This is a known issue with rpy2 on non-English systems")
                elif diagnostic_result and not diagnostic_result.get('is_compatible', True):
                    print(f"    - {diagnostic_result.get('warning', 'Unknown compatibility issue')}")
            return False
    
    except Exception as e:
        if verbose:
            print(f"   Error: {str(e)}")
        warnings.warn(f"Failed to create ChromatinAssay with fragments: {str(e)}")
        return False
    
    finally:
        # Restore original console callbacks
        try:
            callbacks.consolewrite_print = original_console_write
            callbacks.consolewrite_warnerror = original_console_warnerror
        except:
            pass


def _diagnose_fragment_compatibility(atac_adata: ad.AnnData, 
                                   fragments_path: str, 
                                   verbose: bool = False) -> Dict[str, Any]:
    """
    Pre-diagnose the compatibility of fragment files and ATAC data on the Python side
    """
    fragments_path = Path(fragments_path)
    result = {
        'is_compatible': True,
        'warning': None,
        'cell_overlap_rate': 0.0,
        'fragment_cells': 0,
        'atac_cells': 0,
        'matching_cells': 0
    }
    
    try:
        # According to the R version logic, ASCII encoding format is preferred.
        if fragments_path.suffix == '.gz':
            import gzip
            with gzip.open(fragments_path, 'rt', encoding='ascii', errors='ignore') as f:
                sample_lines = []
                for _ in range(1000):
                    line = f.readline()
                    if not line:  # end of file
                        break
                    sample_lines.append(line.strip())
        else:
            with open(fragments_path, 'r', encoding='ascii', errors='ignore') as f:
                sample_lines = []
                for _ in range(1000):
                    line = f.readline()
                    if not line:  # end of file
                        break
                    sample_lines.append(line.strip())
        fragment_barcodes = set()
        for line in sample_lines:
            if line and len(line.split('\t')) >= 4:
                fragment_barcodes.add(line.split('\t')[3])
        atac_barcodes = set(atac_adata.obs.index)
        matching_barcodes = fragment_barcodes.intersection(atac_barcodes)
        overlap_rate = len(matching_barcodes) / len(atac_barcodes) if atac_barcodes else 0
        
        result.update({
            'cell_overlap_rate': overlap_rate,
            'fragment_cells': len(fragment_barcodes),
            'atac_cells': len(atac_barcodes),
            'matching_cells': len(matching_barcodes)
        })
        
        # Determine compatibility
        if overlap_rate < 0.1:  # Overlap rate is less than 10%
            result['is_compatible'] = False
            result['warning'] = f"Low cell overlap rate: {overlap_rate:.1%} ({len(matching_barcodes)}/{len(atac_barcodes)})"
        
        if verbose:
            print(f"    Fragment cells (sample): {len(fragment_barcodes)}")
            print(f"    ATAC cells: {len(atac_barcodes)}")
            print(f"    Matching cells: {len(matching_barcodes)}")
            print(f"    Overlap rate: {overlap_rate:.1%}")
        
    except Exception as e:
        result['is_compatible'] = False
        result['warning'] = f"Cannot read fragment file: {str(e)}"
    
    return result

def _diagnose_fragment_compatibility(atac_adata: ad.AnnData,
                                   fragments_path: str,
                                   verbose: bool = False) -> Dict[str, Any]:
    """
    Diagnose fragment compatibility by scanning the full fragments file.

    A head-only sample can strongly underestimate overlap after trajectory
    filtering, because the selected cells may appear later in the fragments
    file. The scan stops early if all ATAC cells are already matched.
    """
    fragments_path = Path(fragments_path)
    result = {
        'is_compatible': True,
        'warning': None,
        'cell_overlap_rate': 0.0,
        'fragment_cells': 0,
        'atac_cells': 0,
        'matching_cells': 0,
    }

    try:
        atac_barcodes = set(map(str, atac_adata.obs.index))
        fragment_barcodes = set()
        matching_barcodes = set()

        if fragments_path.suffix == '.gz':
            import gzip
            opener = gzip.open
        else:
            opener = open

        with opener(fragments_path, 'rt', encoding='ascii', errors='ignore') as f:
            for line in f:
                if not line or line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 4:
                    continue
                barcode = parts[3]
                fragment_barcodes.add(barcode)
                if barcode in atac_barcodes:
                    matching_barcodes.add(barcode)

                if len(matching_barcodes) == len(atac_barcodes):
                    break

        overlap_rate = len(matching_barcodes) / len(atac_barcodes) if atac_barcodes else 0
        result.update({
            'cell_overlap_rate': overlap_rate,
            'fragment_cells': len(fragment_barcodes),
            'atac_cells': len(atac_barcodes),
            'matching_cells': len(matching_barcodes),
        })

        if overlap_rate < 0.1:
            result['is_compatible'] = False
            result['warning'] = (
                f"Low cell overlap rate: {overlap_rate:.1%} "
                f"({len(matching_barcodes)}/{len(atac_barcodes)})"
            )

        if verbose:
            print(f"    Fragment cells: {len(fragment_barcodes)}")
            print(f"    ATAC cells: {len(atac_barcodes)}")
            print(f"    Matching cells: {len(matching_barcodes)}")
            print(f"    Overlap rate: {overlap_rate:.1%}")

    except Exception as e:
        result['is_compatible'] = False
        result['warning'] = f"Cannot read fragment file: {str(e)}"

    return result


def _diagnose_failure_reason(r_error: str, diagnostic_result: Dict[str, Any]):
    """Analyze failure causes based on R errors and pre-diagnostic information"""
    print(f"  Diagnostic analysis:")
    
    if "No such file" in r_error or "cannot open" in r_error:
        print(f"    - Fragment file access issue")
        print(f"    - Check file path and permissions")
    
    elif diagnostic_result['cell_overlap_rate'] < 0.1:
        print(f"    - Cell barcode mismatch (overlap: {diagnostic_result['cell_overlap_rate']:.1%})")
        print(f"    - Fragment cells: {diagnostic_result['fragment_cells']}")
        print(f"    - ATAC cells: {diagnostic_result['atac_cells']}")
        print(f"    - Matching: {diagnostic_result['matching_cells']}")
    
    elif "Invalid" in r_error or "malformed" in r_error:
        print(f"    - Fragment file format issue")
        print(f"    - Expected: chr\tstart\tend\tbarcode\tdup_count")
    
    else:
        print(f"    - Unknown R error: {r_error}")

def _suggest_solutions(diagnostic_result: Dict[str, Any]):
    """Provide solution suggestions based on diagnosis results"""
    print(f"  Suggested solutions:")
    
    if diagnostic_result['cell_overlap_rate'] < 0.1:
        print(f"    1. Check if fragment file matches ATAC data")
        print(f"    2. Verify cell barcode format consistency")
        print(f"    3. Consider using subset of cells with fragments")
    
    print(f"    4. Check fragment file format (5-column TSV)")
    print(f"    5. Verify file permissions and accessibility")
