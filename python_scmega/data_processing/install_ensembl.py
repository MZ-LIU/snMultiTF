# """
# """
# import os
# import urllib.request
# from pathlib import Path
# import gzip
# import shutil


# def download_ensembl_annotation(cache_dir=None, genome="hg38", release=86):
#     """
    
#     Args:
#     """
    
#     if cache_dir is None:
#         cache_dir =  "E:/scMEGA9.25/python_scmega/data/.pyensembl"
#     else:
#         cache_dir = Path(cache_dir)
    
#     if genome == "hg38":
#         genome_name = "GRCh38"
#     else:
    
#     target_dir = cache_dir / genome_name / f"ensembl{release}"
#     target_dir.mkdir(parents=True, exist_ok=True)
    
    
#     base_url = f"http://ftp.ensembl.org/pub/release-{release}"
    
#     files_to_download = [
#         {
#             "url": f"{base_url}/gtf/homo_sapiens/Homo_sapiens.{genome_name}.{release}.gtf.gz",
#             "filename": f"Homo_sapiens.{genome_name}.{release}.gtf.gz",
#             "required": True
#         },
#         {
#             "url": f"{base_url}/fasta/homo_sapiens/cdna/Homo_sapiens.{genome_name}.cdna.all.fa.gz",
#             "filename": f"Homo_sapiens.{genome_name}.cdna.all.fa.gz",
#             "required": False
#         },
#         {
#             "url": f"{base_url}/fasta/homo_sapiens/pep/Homo_sapiens.{genome_name}.pep.all.fa.gz",
#             "filename": f"Homo_sapiens.{genome_name}.pep.all.fa.gz",
#             "required": False
#         }
#     ]
    
#     for file_info in files_to_download:
#         target_file = target_dir / file_info["filename"]
        
#         if target_file.exists():
#             continue
        
#         print(f"  URL: {file_info['url']}")
        
#         try:
#             def progress_hook(block_num, block_size, total_size):
#                 if total_size > 0:
#                     percent = block_num * block_size / total_size * 100
#                     mb_downloaded = block_num * block_size / (1024 * 1024)
#                     mb_total = total_size / (1024 * 1024)
            
#             urllib.request.urlretrieve(
#                 file_info['url'], 
#                 target_file,
#                 reporthook=progress_hook
#             )
            
#         except Exception as e:
#             if file_info['required']:
#             else:
    
#     print("\n" + "="*70)
#     print("="*70)
    
#     try:
#         import pyensembl
        
#         os.environ['PYENSEMBL_CACHE_DIR'] = str(cache_dir)
        
#         ensembl = pyensembl.EnsemblRelease(release)
        
#         if not ensembl.is_installed():
#             ensembl.index()
#         else:
        
#         genes = ensembl.genes()
        
#         return True
        
#     except ImportError:
#         print("    import pyensembl")
#         print(f"    ensembl = pyensembl.EnsemblRelease({release})")
#         print("    ensembl.index()")
#         return False
    
#     except Exception as e:
#         return False


# if __name__ == "__main__":
#     import sys
#     project_dir = "E:/scMEGA9.25/python_scmega/data/.pyensembl"
#     download_ensembl_annotation(cache_dir=project_dir / "data")
    
    
