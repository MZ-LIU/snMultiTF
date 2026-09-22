"""
ChromVAR Analysis using rpy2

This module provides R-based chromVAR analysis through rpy2, ensuring complete
compatibility with the R version of scMEGA.

Key functions:
- load_jaspar_motifs_r: Load JASPAR motif database using R
- add_motifs_r: Add motif information using R's AddMotifs
- run_chromvar_r: Run chromVAR analysis using R

References:
- ChromVAR: Schep et al. (2017) Nature Methods
- Signac: Stuart et al. (2021) Nature Methods
- JASPAR: Fornes et al. (2020) Nucleic Acids Research
"""

import numpy as np
import pandas as pd
import anndata as ad
from typing import Optional, Dict, Any, List, Union
import warnings

from ..data_processing.format_conversion import _atac_anndata_to_seurat_r

# Try to import MuData, set to None if not available
try:
    import muon as mu
    MUON_AVAILABLE = True
except ImportError:
    mu = None
    MUON_AVAILABLE = False




def load_jaspar_motifs_r(collection: str = "CORE",
                        species: Optional[str] = "Homo sapiens",
                        version: str = "JASPAR2020",
                        tax_group: Optional[str] = None,
                        **kwargs) -> Any:
    """
    Load JASPAR motif database using R.
    
    R equivalent:
    ```R
    library(JASPAR2020)
    library(TFBSTools)
    
    pfm <- getMatrixSet(
        x = JASPAR2020,
        opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
    )
    ```
    
    Args:
        collection: JASPAR collection ("CORE", "UNVALIDATED", "PBM", etc.)
        species: Species name ("Homo sapiens", "Mus musculus", etc.).
            Set to None to omit species from getMatrixSet opts.
        version: JASPAR version ("JASPAR2020", "JASPAR2022", etc.)
        tax_group: JASPAR taxonomic group (e.g. "vertebrates"). When provided,
            it is added to getMatrixSet opts.
        **kwargs: Additional parameters for getMatrixSet
        
    Returns:
        R PFMatrixList object containing motif data
        
    Example:
        ```python
        # Load human core motifs
        motifs = load_jaspar_motifs_r(
            collection="CORE",
            species="Homo sapiens"
        )
        ```
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.conversion import localconverter
        
        print(f"Loading JASPAR motifs using R ({version})...")
        
        # Import required R packages
        try:
            jaspar_pkg = importr(version)
            tfbstools = importr('TFBSTools')
        except Exception as e:
            raise ImportError(
                f"Failed to import R packages. Please install them in R:\n"
                f"  install.packages('BiocManager')\n"
                f"  BiocManager::install(c('{version}', 'TFBSTools'))\n"
                f"\nError: {e}"
            )
        
        # Use context manager for conversion
        with localconverter(ro.default_converter + pandas2ri.converter):
            # Get JASPAR database object.
            #
            # Older JASPAR packages exposed the DB object directly as JASPAR20xx.
            # Newer packages such as JASPAR2024 expose JASPAR2024 as a constructor
            # function, so passing ro.r[version] directly to getMatrixSet() sends a
            # function and triggers:
            #   unable to find an inherited method for getMatrixSet, x = "function"
            jaspar_db = ro.r(
                f"if (is.function({version})) {version}() else {version}"
            )
            
            # Create opts list for getMatrixSet
            opts_dict = {
                'collection': collection,
                'all_versions': False
            }
            if species is not None:
                opts_dict['species'] = species
            if tax_group is not None:
                opts_dict['tax_group'] = tax_group
            opts = ro.ListVector(opts_dict)
            
            # Add any additional options
            if kwargs:
                opts_dict = dict(opts.items())
                opts_dict.update(kwargs)
                opts = ro.ListVector(opts_dict)
            
            # Call getMatrixSet
            opts_preview = ", ".join(f"{k}={v!r}" for k, v in dict(opts.items()).items())
            print(f"  Calling getMatrixSet({opts_preview})...")
            pfm = tfbstools.getMatrixSet(jaspar_db, opts=opts)
            
            # Get motif count
            n_motifs = int(ro.r.length(pfm)[0])
            print(f"   Loaded {n_motifs} motifs from JASPAR")
        
        return pfm
        
    except ImportError as e:
        raise ImportError(
            f"rpy2 is not installed or configured correctly.\n"
            f"Please install rpy2: pip install rpy2\n"
            f"And ensure R is installed with BiocManager, {version}, and TFBSTools packages.\n"
            f"\nError: {e}"
        )
    except Exception as e:
        raise RuntimeError(f"Failed to load JASPAR motifs: {e}")


def add_motifs_r(data: Union[ad.AnnData, 'mu.MuData'],
                motifs: Any,
                genome: str = "hg38",
                assay: str = "ATAC",
                fragments_path: Optional[str] = None,
                **kwargs) -> Dict[str, Any]:
    """
    Add motif information using R's AddMotifs function from Signac.
    
    R equivalent:
    ```R
    library(Signac)
    library(BSgenome.Hsapiens.UCSC.hg38)
    
    obj <- AddMotifs(
        object = obj,
        genome = BSgenome.Hsapiens.UCSC.hg38,
        pfm = pfm,
        assay = "ATAC"
    )
    ```
    
    Args:
        data: MultiomeData or AnnData object with ATAC-seq data
        motifs: R PFMatrixList object from load_jaspar_motifs_r
        genome: Genome build ("hg38", "hg19", "mm10", "mm9")
        assay: Which assay to use (default: "ATAC")
        fragments_path: Path to fragments file. If provided, this will be used
                       instead of extracting from data object metadata.
        **kwargs: Additional parameters for AddMotifs

    Returns:
        Modified MultiomeData or AnnData object with motif information added
            
    Example:
        ```python
        # Load motifs
        motifs = load_jaspar_motifs_r()
        
        # Add motifs to data - returns modified data object
        multiome_data = add_motifs_r(
            multiome_data,
            motifs=motifs,
            genome="hg38",
            fragments_path="/path/to/fragments.tsv.gz"
        )
        
        # Now multiome_data contains motif information
        # Can proceed to run chromVAR
        multiome_data = run_chromvar_r(multiome_data, genome="hg38")
        ```
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr
        from rpy2.robjects import pandas2ri, numpy2ri
        from rpy2.robjects.conversion import localconverter
        
        print(f"\nAdding motifs using R's AddMotifs (genome={genome})...")
        
        # Import required R packages
        try:
            signac = importr('Signac')
            seurat = importr('Seurat')
            
            # Load appropriate BSgenome
            genome_map = {
                'hg38': 'BSgenome.Hsapiens.UCSC.hg38',
                'hg19': 'BSgenome.Hsapiens.UCSC.hg19',
                'mm10': 'BSgenome.Mmusculus.UCSC.mm10',
                'mm9': 'BSgenome.Mmusculus.UCSC.mm9'
            }
            
            if genome not in genome_map:
                raise ValueError(f"Unsupported genome: {genome}. Use one of {list(genome_map.keys())}")
            
            bsgenome_name = genome_map[genome]
            print(f"  Loading {bsgenome_name}...")
            bsgenome = importr(bsgenome_name)
            genome_obj = ro.r[bsgenome_name]
            
        except Exception as e:
            raise ImportError(
                f"Failed to import R packages. Please install them in R:\n"
                f"  install.packages('BiocManager')\n"
                f"  BiocManager::install(c('Signac', 'Seurat', '{genome_map.get(genome, 'BSgenome')}'))\n"
                f"\nError: {e}"
            )
        
        # Use context manager for conversion
        with localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
            if MUON_AVAILABLE and isinstance(data, mu.MuData):
                print("  Detected MuData object")
                atac_adata = data['atac']
                if fragments_path is None:
                    print(f"  No fragments_path provided")
                else:
                    print(f"  Using fragments_path from function parameter: {fragments_path}")
            else:
                print("  Detected AnnData object")
                atac_adata = data
                # Prefer fragments_path in function parameters, or extract from uns if not provided
                if fragments_path is None:
                    fragments_path = atac_adata.uns.get('fragments', {}).get('path', None)
                    print(f"  Using fragments_path from uns: {fragments_path}")
                else:
                    print(f"  Using fragments_path from function parameter: {fragments_path}")
            
            print("  Converting ATAC data to Seurat and associating fragments...")
            success = _atac_anndata_to_seurat_r(
                atac_adata=atac_adata,
                fragments_path=fragments_path,
                genome=genome,
                verbose=True
            )
            
            if not success:
                warnings.warn("Fragment association failed, but continuing with AddMotifs...")
            
            seurat_obj = ro.r['seurat_obj']
            
            # Call AddMotifs
            print(f"  Running AddMotifs...")
            ro.r.assign('seurat_obj', seurat_obj)
            ro.r.assign('pfm', motifs)
            ro.r.assign('genome_obj', genome_obj)
            
            # Build R command
            r_cmd = f"""
            library(Signac)
            library(Seurat)
            set.seed(42)
            
            seurat_obj <- AddMotifs(
                object = seurat_obj,
                genome = genome_obj,
                pfm = pfm,
                assay = "{assay}"
            )
            """

            # Execute R command
            ro.r(r_cmd)
            seurat_obj_with_motifs = ro.r['seurat_obj']
            
            # Extract motif matching matrix
            print("  Extracting motif matching matrix...")
            motif_data = _extract_motif_matrix_r(seurat_obj_with_motifs, assay=assay)
            
            print(f"   AddMotifs completed")
            print(f"   Motif matching matrix: {motif_data['motif_matching'].shape}")

            motif_info = {
                'motif_matching': motif_data['motif_matching'],
                'motif_names': motif_data['motif_names'],
                'seurat_obj': seurat_obj_with_motifs,  # Storing R object references
                'genome': genome,
                'assay': assay,
                'summary': motif_data.get('summary', {})
            }
            
            if MUON_AVAILABLE and isinstance(data, mu.MuData):
                data['atac'].uns['motifs'] = motif_info
                print(f"   Motif data stored in MuData['atac'].uns['motifs']")
                print(f"     - motif_matching: pandas DataFrame (for GRN inference)")
                print(f"     - seurat_obj: R object reference (temporary, for chromVAR)")
            else:
                data.uns['motifs'] = motif_info
                print(f"   Motif data stored in AnnData.uns['motifs']")
        
        return data
        
    except ImportError as e:
        raise ImportError(f"rpy2 or R packages not available: {e}")
    except Exception as e:
        raise RuntimeError(f"AddMotifs failed: {e}")
        


def run_chromvar_r(data: Union[ad.AnnData, 'mu.MuData'],
                  genome: str = "hg38",
                  assay: str = "ATAC",
                  fragments_path: Optional[str] = None,
                  **kwargs) -> Union[ad.AnnData, 'mu.MuData']:
    """
    Run chromVAR analysis using R's RunChromVAR function from Signac.
    
    For MuData objects, creates an independent 'chromvar' modality following
    the standard MuData structure specification.
    
    R equivalent:
    ```R
    library(Signac)
    library(chromVAR)
    library(BSgenome.Hsapiens.UCSC.hg38)
    
    obj <- RunChromVAR(
        object = obj,
        genome = BSgenome.Hsapiens.UCSC.hg38,
        assay = "ATAC"
    )
    ```
    
    Args:
        data: MuData or AnnData object with motif annotations (from add_motifs_r)
        genome: Genome build ("hg38", "hg19", "mm10", "mm9")
        assay: Which assay to use (default: "ATAC")
        fragments_path: Path to fragments file (optional, extracted from metadata if not provided)
        **kwargs: Additional parameters for RunChromVAR
        
    Returns:
    Modified MuData object with chromVAR results
    - Creates independent chromvar modality at mdata['chromvar']
    - chromVAR results stored as: mdata['chromvar'].X, .layers, .var_names
            
    Example:
        ```python
        # For MuData (recommended)
        motifs = load_jaspar_motifs_r()
        mdata = add_motifs_r(mdata, motifs=motifs, genome="hg38")
        mdata = run_chromvar_r(mdata, genome="hg38")
        
        # Access chromVAR results
        chromvar_modality = mdata['chromvar']  # Independent AnnData modality
        tf_activities = chromvar_modality.X     # TF activity matrix
        z_scores = chromvar_modality.layers['z_scores']
        
        # For GRN inference
        motif_matching = mdata['atac'].uns['motifs']['motif_matching']
        ```
        
    Note:
        The chromVAR results are stored following MuData standard structure:
        - mdata['chromvar']: Independent AnnData modality
        - mdata['chromvar'].X: TF activity matrix (n_cells × n_TFs)
        - mdata['chromvar'].layers: z_scores, deviations
        - mdata['atac'].uns['motifs']['motif_matching']: Preserved for GRN inference
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr
        from rpy2.robjects import pandas2ri, numpy2ri
        from rpy2.robjects.conversion import localconverter
        
        print(f"\nRunning chromVAR using R (genome={genome})...")
        
        if MUON_AVAILABLE and isinstance(data, mu.MuData):
            print("  Detected MuData object")
            atac_adata = data['atac']
            motif_info = data['atac'].uns.get('motifs', None)
        else:
            print("  Detected AnnData object")
            atac_adata = data
            motif_info = data.uns.get('motifs', None)
        
        print(f"  Checking for motif data and Seurat object...")
        if motif_info is None:
            raise ValueError(
                "\n" + "="*70 + "\n"
                " Motif data not found!\n"
                "="*70 + "\n"
                "Please run add_motifs_r() before run_chromvar_r().\n"
                + "="*70
            )
        
        if 'seurat_obj' not in motif_info:
            raise ValueError(
                "\n" + "="*70 + "\n"
                " R Seurat object not found in motif data!\n"
                "="*70 + "\n"
                "Please re-run add_motifs_r().\n"
                + "="*70
            )
        
        seurat_obj_with_motifs = motif_info['seurat_obj']
        print(f"   Found motif data: {len(motif_info['motif_names'])} motifs")
        print(f"   Retrieved R Seurat object from data object")
        
        # Import R packages
        try:
            signac = importr('Signac')
            chromvar = importr('chromVAR')
            seurat = importr('Seurat')
            
            genome_map = {
                'hg38': 'BSgenome.Hsapiens.UCSC.hg38',
                'hg19': 'BSgenome.Hsapiens.UCSC.hg19',
                'mm10': 'BSgenome.Mmusculus.UCSC.mm10',
                'mm9': 'BSgenome.Mmusculus.UCSC.mm9'
            }
            
            if genome not in genome_map:
                raise ValueError(f"Unsupported genome: {genome}")
            
            bsgenome_name = genome_map[genome]
            print(f"  Loading {bsgenome_name}...")
            bsgenome = importr(bsgenome_name)
            genome_obj = ro.r[bsgenome_name]
            
        except Exception as e:
            raise ImportError(
                f"Failed to import R packages:\n"
                f"  BiocManager::install(c('Signac', 'chromVAR','BiocParallel', '{genome_map.get(genome, 'BSgenome')}'))\n"
                f"\nError: {e}"
            )
        
        # SerialParam() critical fix: force single core to avoid Windows parallel serialization errors
        # Use context manager for conversion
        with localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
            print(f"  Running RunChromVAR on existing Seurat object...")
            ro.r.assign('seurat_obj_with_motifs', seurat_obj_with_motifs)
            ro.r.assign('genome_obj', genome_obj)
            
            r_cmd = f"""
            library(Signac)
            library(chromVAR)
            library(Seurat)
            library(BiocParallel)

            set.seed(42)
            BiocParallel::register(BiocParallel::SerialParam())
            
            seurat_obj_with_chromvar <- RunChromVAR(
                object = seurat_obj_with_motifs,
                genome = genome_obj,
                assay = "{assay}"
            )
            """
            
            ro.r(r_cmd)
            seurat_obj_with_chromvar = ro.r['seurat_obj_with_chromvar']
            
            print("  Extracting chromVAR results to Python format...")
            chromvar_data = _extract_chromvar_results_r(seurat_obj_with_chromvar)

            # Make sure to include motif_names
            if 'motif_names' not in chromvar_data:
                chromvar_data['motif_names'] = motif_info['motif_names']
            print(f"   chromVAR completed")
            print(f"   Motif activities: {chromvar_data['motif_activities'].shape}")

    
            if MUON_AVAILABLE and isinstance(data, mu.MuData):
                chromvar_adata = ad.AnnData(
                    X=chromvar_data['motif_activities'].values,  # TF activity matrix (cells × TFs)
                    obs=data['rna'].obs.copy(),  # Share cell metadata
                    var=pd.DataFrame(
                        index=chromvar_data['motif_names'],
                        data={'motif_name': chromvar_data['motif_names']}
                    )
                )
                
                # if 'z_scores' in chromvar_data and chromvar_data['z_scores'] is not None:
                #     chromvar_adata.layers['z_scores'] = chromvar_data['z_scores']
                # if 'deviations' in chromvar_data and chromvar_data['deviations'] is not None:
                #     chromvar_adata.layers['deviations'] = chromvar_data['deviations']
                if 'motif_annotations' in chromvar_data:
                    chromvar_adata.uns['motif_annotations'] = chromvar_data['motif_annotations']
                
                chromvar_adata.uns['summary'] = chromvar_data.get('summary', {
                    'n_motifs': len(chromvar_data['motif_names']),
                    'n_cells': chromvar_data['motif_activities'].shape[0]
                })
                
                data.mod['chromvar'] = chromvar_adata
                # NOTE: There is no need to call data.update() because chromvar is a standalone modal

                print(f"   chromvar modality created: {chromvar_adata.shape}")
                print(f"    - X: TF activity matrix ({chromvar_adata.n_obs} cells × {chromvar_adata.n_vars} TFs)")
                if chromvar_adata.layers:
                    print(f"    - layers: {', '.join(chromvar_adata.layers.keys())}")
                print(f"    - uns: {', '.join(chromvar_adata.uns.keys())}")
                
                # -----------------------------------------------------------
                # -----------------------------------------------------------
                if 'seurat_obj' in data['atac'].uns.get('motifs', {}):
                    del data['atac'].uns['motifs']['seurat_obj']
                    print(f"   Cleaned up temporary R Seurat object reference")
                
                # -----------------------------------------------------------
                # Keep motif_matching (required for GRN inference)
                # -----------------------------------------------------------
                print(f"   Preserved motif_matching in mdata['atac'].uns['motifs']")
                
                print("\n" + "="*70)
                print("chromVAR Results Stored (MuData Standard Format)")
                print("="*70)
                print(f" Independent modality: mdata['chromvar']")
                print(f"  - Type: AnnData")
                print(f"  - Shape: {chromvar_adata.shape}")
                print(f"  - Access TF activities: mdata['chromvar'].X")
                # print(f"  - Access Z-scores: mdata['chromvar'].layers['z_scores']")
                print(f"\n For GRN inference:")
                print(f"  - Motif matching: mdata['atac'].uns['motifs']['motif_matching']")
                print("="*70)
                
            else:
                error_msg = (
                    "\n" + "="*70 + "\n"
                    "chromVAR requires MuData object!\n"
                    "="*70 + "\n"
                    "This function only supports MuData objects following the standard\n"
                    "multi-omics data structure. AnnData objects are no longer supported.\n"
                    "\n"
                    "Please convert your data to MuData format:\n"
                    "  import muon as mu\n"
                    "  mdata = mu.MuData({'rna': rna_adata, 'atac': atac_adata})\n"
                    "\n"
                    "Then run chromVAR:\n"
                    "  mdata = run_chromvar_r(mdata, genome='hg38')\n"
                    "="*70
                )
                raise TypeError(error_msg)              
        
        return data
        
    except ImportError as e:
        raise ImportError(f"rpy2 or R packages not available: {e}")
    except Exception as e:
        raise RuntimeError(f"RunChromVAR failed: {e}")
                


def _extract_motif_matrix_r(seurat_obj: Any, assay: str = "ATAC") -> Dict[str, Any]:
    """Extract motif matching matrix from Seurat object."""
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    
    # Use context manager for conversion
    with localconverter(ro.default_converter + pandas2ri.converter):
        ro.r.assign('seurat_obj', seurat_obj)
        
        # Extract motif data
        r_cmd = f"""
        library(Signac)
        library(Seurat)
        
        # Get motif object
        motif_obj <- Motifs(seurat_obj[["{assay}"]])
        
        # Get motif matching matrix
        motif_matrix <- motif_obj@data
        motif_names <- motif_obj@motif.names
        peak_names <- rownames(motif_matrix)
        
        # Convert to data frame
        motif_df <- as.data.frame(as.matrix(motif_matrix))
        colnames(motif_df) <- motif_names
        rownames(motif_df) <- peak_names
        """
        
        ro.r(r_cmd)
        
        # Get results
        motif_df = ro.r['motif_df']

        def convert_peak_format(peak_name):
            """Convert peak name from 'chr-start-end' to 'chr:start-end'"""
            parts = str(peak_name).split('-')
            if len(parts) == 3:
                chr_name, start, end = parts
                return f"{chr_name}:{start}-{end}"
            return peak_name

        motif_df.index = [convert_peak_format(p) for p in motif_df.index]
        print(f"  Peak format converted: {list(motif_df.index[:3])}")

        import re
        motif_names_raw = list(ro.r['motif_names'])
        motif_names = []
        for name in motif_names_raw:
            name_str = str(name)
            match = re.search(r'"([^"]+)"', name_str)
            if match:
                motif_names.append(match.group(1))
            else:
                motif_names.append(name_str.strip())
    
    return {
        'motif_matching': motif_df,
        'motif_names': motif_names,
        'summary': {
            'n_peaks': motif_df.shape[0],
            'n_motifs': motif_df.shape[1],
            'total_matches': int(motif_df.sum().sum())
        }
    }


def _extract_chromvar_results_r(seurat_obj: Any) -> Dict[str, Any]:
    """Extract chromVAR results from Seurat object."""
    '\n    Extract chromVAR results from Seurat object.\n    \n    Extracts all layers from chromVAR assay:\n    - data: TF activity scores (motif_activities)\n    - scale.data: Z-scores (z_scores) - Normalized activity scores\n    - counts: Deviation scores (deviations) - original deviations\n    \n    Args:\n        seurat_obj: R Seurat object with chromVAR results\n        \n    Returns:\n        Dictionary with:\n        - motif_activities: DataFrame (cells × TFs)\n        - z_scores: DataFrame (cells × TFs) if available\n        - deviations: DataFrame (cells × TFs) if available\n        - summary: dict with metadata\n    '
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    
    # Use context manager for conversion
    with localconverter(ro.default_converter + pandas2ri.converter):
        ro.r.assign('seurat_obj', seurat_obj)
        
        # Extract chromVAR data
        r_cmd = """
        library(Signac)
        library(Seurat)
        
        # Get chromVAR assay
        chromvar_data <- GetAssayData(seurat_obj, assay = "chromvar", layer = "data")
        
        # Convert to data frame (cells x motifs)
        chromvar_df <- as.data.frame(t(as.matrix(chromvar_data)))
        cell_names <- colnames(chromvar_data)
        motif_names <- rownames(chromvar_data)
        
        rownames(chromvar_df) <- cell_names
        colnames(chromvar_df) <- motif_names
        """
        
        ro.r(r_cmd)
        
        # Get results
        chromvar_df = ro.r['chromvar_df']
    
    return {
        'motif_activities': chromvar_df,
        'summary': {
            'n_cells': chromvar_df.shape[0],
            'n_motifs': chromvar_df.shape[1]
        }
    }


