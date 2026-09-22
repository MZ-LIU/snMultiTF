"""
RAG cell annotation system - based on two-layer JSON database
Data source: Cell_marker_JSON.json
Data structure:
  {
    "tissue_name": {
      "gene_index": {"GENE1": ["CellType1", "CellType2"], ...},
      "cell_metadata": {"CellType1": {"total_markers": 10, "cancer_type": "No"}, ...}
    }
  }

Search process:
1. Physical positioning: lock data shards based on Tissue
2. Inverted index: traverse the user gene list and count the number of hits for each cell type
3. Dual scoring:
   - Recall = number of hits / total number of markers for this cell in the database
   - Jaccard = number of hits / (user input length + total number of database markers - number of hits)
4. Candidate screening: Sort by Jaccard score and return to Top-3
5. LLM generation: Pass the list of candidate cells and hit genes to LLM for final judgment.

Modification instructions:
- generate_annotation: Single call to LLM, prompt is consistent with GPT.py, only cell type name is returned
- batch_annotate: save Predict_cell_type and Matched_Genes
-Default model: gpt-5
"""

import re
import time 
import os
import json
import pandas as pd
from collections import defaultdict
from typing import List, Dict, Tuple, Optional, Union
import warnings

try:
    from openai import OpenAI
    OPENAI_AVAILABLE = True
except ImportError:
    OPENAI_AVAILABLE = False
    warnings.warn(' The openai library is not installed and the API function is not available')

try:
    from dotenv import load_dotenv
    load_dotenv('E:\\\u65e5\u5fd7\\snMultiTF\\python_scmega\\.env')
except ImportError:
    warnings.warn(' python-dotenv is not installed, you need to manually set environment variables')


class CellTypeAnnotator:
    """
    RAG-based cell type annotator
    
    parameter:
        json_db_path: JSON database path
        api_key: OpenAI API key (optional, use environment variables first)
        base_url: API base URL (optional, use environment variables first)
        default_model: the name of the model used by default
        verbose: whether to print detailed debugging information
    """
    
    def __init__(self, 
                 json_db_path: str = None,
                 api_key: str = None,
                 base_url: str = None,
                 default_model: str = "gpt-5",
                 verbose: bool = True):
        """Initialize the cell type annotator"""
        
        self.verbose = verbose
        
        if json_db_path is None:
            current_dir = os.path.dirname(os.path.abspath(__file__))
            potential_paths = [
                'E:\\\u65e5\u5fd7\\snMultiTF\\Maker_dataset\\Cell_marker_JSON.json',
                os.path.join(current_dir, '..', 'data', 'Cell_marker_JSON.json')
            ]
            for path in potential_paths:
                if os.path.exists(path):
                    json_db_path = path
                    break
        
        self.json_db_path = json_db_path
        self.default_model = default_model
        self.knowledge_base = {}
        self.client = None
        
        if json_db_path and os.path.exists(json_db_path):
            self._load_database()
        else:
            if self.verbose:
                warnings.warn(f" Database file not found: {json_db_path}")
        
        if OPENAI_AVAILABLE:
            self._init_client(api_key, base_url)
    
    def _load_database(self):
        """Load JSON database"""
        try:
            if self.verbose:
                print("="*70)
                print('Loading JSON database...')
                print("="*70)
            
            with open(self.json_db_path, 'r', encoding='utf-8') as f:
                self.knowledge_base = json.load(f)
            
            if self.verbose:
                print(f" Database loaded successfully")
                print(f"  Number of tissues: {len(self.knowledge_base)}")
                print(f"  Tissue list: {list(self.knowledge_base.keys())[:10]}...")
                print("="*70 + "\n")
        except Exception as e:
            warnings.warn(f" Database loading failed: {e}")
            self.knowledge_base = {}
    
    def _init_client(self, api_key: str = None, base_url: str = None):
        """Initialize OpenAI client"""
        try:
            self.client = OpenAI(
                api_key=api_key or os.getenv("OPENAI_API_KEY"),
                base_url=base_url or os.getenv("OPENAI_BASE_URL"),
            )
            if self.verbose:
                print(' OpenAI Client initialized successfully')
        except Exception as e:
            warnings.warn(f" OpenAI Client initialization failed: {e}")
            self.client = None
    
    def retrieve_candidates(self, 
                          user_genes: Union[str, List[str]], 
                          user_tissue: str = None,
                          top_k: int = 3) -> List[Dict]:
        """
        Candidate cell type retrieval based on two-layer JSON database - supports cross-tissue retrieval
        """
        if isinstance(user_genes, str):
            query_genes = set([g.strip().upper() for g in user_genes.split(',') if g.strip()])
        else:
            query_genes = set([g.strip().upper() for g in user_genes if g.strip()])
        
        tissue_key = user_tissue.strip().capitalize() if user_tissue and user_tissue.strip() else None
        
        if self.verbose:
            print(f"Query tissue: {tissue_key if tissue_key else 'Not specified (cross-tissue search)'}")
        
        cell_hits = defaultdict(lambda: {"genes": []})
        
        if self.verbose:
            total_tissues = len(self.knowledge_base)
        
        for current_tissue, tissue_data in self.knowledge_base.items():
            gene_index = tissue_data.get("gene_index", {})
            
            for db_gene, matched_cells in gene_index.items():
                match = re.match(r'^([^\(]+)(?:\s*\(([^\)]+)\))?$', db_gene)
                if match:
                    main_gene = match.group(1).strip().upper()
                    alias_gene = match.group(2).strip().upper() if match.group(2) else None
                    
                    matched_gene_name = None
                    if main_gene in query_genes:
                        matched_gene_name = main_gene
                    elif alias_gene and alias_gene in query_genes:
                        matched_gene_name = alias_gene
                    
                    if matched_gene_name:
                        for cell_type in matched_cells:
                            composite_key = (cell_type, current_tissue)
                            cell_hits[composite_key]["genes"].append(matched_gene_name)
        
        if self.verbose:
            print(f" Cross-tissue search completed, {len(cell_hits)} (cell type, tissue) hit entries found")
        
        candidates = []
        user_gene_count = len(query_genes)
        
        for (cell_type, current_tissue), hit_data in cell_hits.items():
            matched_genes = list(set(hit_data["genes"]))
            
            if current_tissue not in self.knowledge_base:
                continue
            
            tissue_data = self.knowledge_base[current_tissue]
            cell_metadata = tissue_data.get("cell_metadata", {})
            metadata = cell_metadata.get(cell_type, {})
            total_markers = metadata.get("total_markers", 0)
            cancer_type = metadata.get("cancer_type", "Unknown")
            
            hits = len(matched_genes)
            recall = hits / total_markers if total_markers > 0 else 0
            jaccard_denom = user_gene_count + total_markers - hits
            jaccard = hits / jaccard_denom if jaccard_denom > 0 else 0
            
            if tissue_key is None:
                tissue_weight = 1.0
            else:
                tissue_weight = 1.5 if current_tissue == tissue_key else 1.0
            
            weighted_jaccard = jaccard * tissue_weight
            
            candidates.append({
                "cell_type": cell_type,
                "tissue": current_tissue,
                "hits": hits,
                "recall": round(recall, 4),
                "jaccard": round(jaccard, 4),
                "weighted_jaccard": round(weighted_jaccard, 4),
                "matched_genes": matched_genes,
                "total_markers": total_markers,
                "cancer_type": cancer_type,
                "tissue_weight": tissue_weight
            })
        
        candidates.sort(key=lambda x: x["weighted_jaccard"], reverse=True)
        top_candidates = candidates[:top_k]
        
        if self.verbose:
            if tissue_key is None:
                print(f"  Cross-tissue search: {len(top_candidates)} candidates in total (tissue not specified)")
            else:
                same_tissue_count = sum(1 for c in top_candidates if c["tissue"] == tissue_key)
                diff_tissue_count = len(top_candidates) - same_tissue_count
                print(f"  Same-tissue candidates: {same_tissue_count} | Cross-tissue candidates: {diff_tissue_count}")
        
        return top_candidates
    
    @staticmethod
    def _count_tokens(text: str) -> int:
        """Estimate the number of tokens in a text"""
        if not text:
            return 0
        chinese_chars = len([c for c in text if 0x4e00 <= ord(c) <= 0x9fff])
        english_chars = len([c for c in text if c.isalpha()])
        numbers = len([c for c in text if c.isdigit()])
        symbols = len([c for c in text if not c.isalpha() and not c.isdigit() and not c.isspace()])
        tokens = chinese_chars + english_chars * 0.5 + numbers * 0.5 + symbols * 0.3
        return int(tokens)
    
    @staticmethod
    def _count_message_tokens(messages: List[Dict]) -> int:
        """Calculate the total number of tokens in the message list"""
        total_tokens = 0
        for message in messages:
            total_tokens += 2
            total_tokens += CellTypeAnnotator._count_tokens(message.get("role", ""))
            total_tokens += CellTypeAnnotator._count_tokens(message.get("content", ""))
        total_tokens += 3
        return int(total_tokens)
    
    def _call_api(self,
                  system_prompt: str,
                  user_content: str,
                  model: str = None,
                  max_retries: int = 3) -> Tuple[Optional[str], Dict]:
        """Call GPT API"""
        if self.client is None:
            return None, {"error": ' API Client is not initialized'}
        
        model = model or self.default_model
        messages = [
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_content}
        ]
        
        for attempt in range(1, max_retries + 1):
            try:
                if self.verbose:
                    print(f"   Calling API... ({attempt}/{max_retries})")
                
                start_time = time.time()
                
                response = self.client.chat.completions.create(
                    model=model,
                    messages=messages,
                    temperature=0.3,
                    max_tokens=10000
                )
                
                end_time = time.time()
                text = response.choices[0].message.content.strip()
                
                if hasattr(response, "usage") and response.usage:
                    stats = {
                        "model": model,
                        "prompt_tokens": response.usage.prompt_tokens,
                        "completion_tokens": response.usage.completion_tokens,
                        "total_tokens": response.usage.total_tokens,
                        "duration": round(end_time - start_time, 2),
                        "attempts": attempt
                    }
                else:
                    stats = {
                        "model": model,
                        "prompt_tokens": self._count_message_tokens(messages),
                        "completion_tokens": self._count_tokens(text),
                        "total_tokens": None,
                        "duration": round(end_time - start_time, 2),
                        "attempts": attempt
                    }
                
                if self.verbose:
                    print(f"   API successful | {stats['duration']}s | tokens: {stats.get('total_tokens','?')}")
                
                return text, stats
            
            except Exception as e:
                if self.verbose:
                    print(f"   API failed ({attempt}/{max_retries}): {str(e)}")
                
                if attempt == max_retries:
                    return None, {"error": str(e), "attempts": attempt, "model": model}
                
                wait_time = 2 ** (attempt - 1)
                if self.verbose:
                    print(f"   Wait for {wait_time}s and try again...")
                time.sleep(wait_time)
        
        return None, {"error": 'unknown error'}
    
    def generate_annotation(self,
                          user_genes: Union[str, List[str]],
                          user_tissue: str,
                          candidates: List[Dict],
                          use_api: bool = True,
                          model: str = None) -> Dict:
        """
        Generate cell annotations based on search results (single call, logic consistent with GPT.py).
        """
        if not candidates:
            return {
                "cell_type": 'Unable to judge',
                "jaccard_score": 0,
                "recall_score": 0,
                "matched_genes": [],
                "stats": None
            }

        if isinstance(user_genes, str):
            gene_tokens = [g for g in re.split(r"[,;\s]+", user_genes) if g]
        else:
            gene_tokens = [str(g) for g in user_genes]
        gene_tokens = [g.strip().upper() for g in gene_tokens if str(g).strip()]
        genes_str = ", ".join(gene_tokens)

        user_tissue_normalized = user_tissue.strip().capitalize() if user_tissue and user_tissue.strip() else None

        context_str = ""
        for i, cand in enumerate(candidates):
            matched_genes_str = ', '.join(cand['matched_genes'])
            if user_tissue_normalized is None:
                tissue_match_indicator = 'Across organizations'
            else:
                tissue_match_indicator = 'Same organization' if cand['tissue'] == user_tissue_normalized else 'Across organizations'
            context_str += f"""
    [Candidate {i+1}]: {cand['cell_type']} (Tissue: {cand['tissue']} | {tissue_match_indicator})
    - Hit marker gene: {matched_genes_str}
    -------------------------------------------
    """

        system_prompt = """You are an expert in single-cell cell type annotation. Your task is to identify the cell type based on DEGs (Differentially Expressed Genes) and tissue context.

    ### Guidelines:
    1. **Filter Noise**: Prioritize specific marker genes. Ignore housekeeping, ribosomal, or cell-cycle genes.
    2. **Evaluate Candidates**: Use the provided candidate list as a baseline, focusing on "same-tissue" matches.
    3. **Override Rule**: If the DEGs strongly indicate a cell type NOT in the candidate list, output the correct biological name instead.
    4. **Standard**: Only output a cell type in precise Cell Ontology (CL) Label, no explanations.

    Examples of correct outputs:
    - "Macrophage"
    - "CD8+ T cell"
    - "Endothelial cell"
    - "Unknown" (if evidence is insufficient)
    """

        tissue_display = user_tissue if (user_tissue and user_tissue.strip()) else "Unspecified"
        user_content = f"""
    ### Data for Annotation:
    - **Tissue**: Human {tissue_display}
    - **Top DEGs**: {genes_str}
    - **Reference Candidates**: {context_str}
    ### Instruction:
    Identify and output the cell type name (English only).
    """

        if not use_api:
            if self.verbose:
                print("\n" + "="*60)
                print('User input:')
                print("="*60)
                print(user_content)
                print("="*60 + "\n")
            return {
                "cell_type": '(API not called)',
                "jaccard_score": candidates[0]['weighted_jaccard'],
                "recall_score": candidates[0]['recall'],
                "matched_genes": candidates[0]['matched_genes'],
                "stats": None
            }

        if self.verbose:
            print(f"\n is calling AI for cell annotation...")

        response, stats = self._call_api(system_prompt, user_content, model=model)

        if self.verbose:
            print(f"Debug information:")
            print(f"   API return value type: {type(response)}")
            print(f"   API return value content: '{response}'")
            print(f"   Statistics: {stats}")

        final_cell_type = response.strip() if response and response.strip() else 'Annotation failed'

        return {
            "cell_type": final_cell_type,
            "jaccard_score": candidates[0]['weighted_jaccard'],
            "recall_score": candidates[0]['recall'],
            "matched_genes": candidates[0]['matched_genes'],
            "stats": stats
        }
    
    def annotate_cell_type(self,
                          genes: Union[str, List[str]],
                          tissue: str,
                          top_k: int = 3,
                          use_api: bool = True,
                          model: str = None) -> Tuple[Dict, List[Dict], Optional[Dict]]:
        """Complete cell annotation process"""
        if self.verbose:
            print("\n" + "="*80)
            print('Start the cell type annotation process')
            print("="*80)
            print('\n[Step 1] Search the database...')

        candidates = self.retrieve_candidates(genes, tissue, top_k=top_k)

        if self.verbose and candidates:
            print(f"\n retrieved {len(candidates)} candidate cell types:")
            tissue_normalized = tissue.strip().capitalize() if tissue and tissue.strip() else None
            for i, cand in enumerate(candidates):
                if tissue_normalized is None:
                    tissue_indicator = "-"
                else:
                    tissue_indicator = "" if cand['tissue'] == tissue_normalized else "×"
                print(f"\n{i+1}. {cand['cell_type']} [{tissue_indicator} {cand['tissue']}]")
                print(f"   Weighted Jaccard: {cand['weighted_jaccard']} (original: {cand['jaccard']}, weight: {cand['tissue_weight']}x)")
                print(f"   Recall: {cand['recall']} | Number of hits: {cand['hits']}")
                print(f"   Hit gene: {', '.join(cand['matched_genes'][:10])}...")

            print("\n" + "="*80)
            print('[Step 2] Generate cell annotations...')
            print("="*80)

        annotation = self.generate_annotation(
            genes, tissue, candidates,
            use_api=use_api, model=model
        )

        return annotation, candidates, annotation.get('stats')
    
    def batch_annotate(self,
                      input_file: str,
                      output_file: str = None,
                      tissue_col: str = 'Tissue',
                      genes_col: str = 'Markers',
                      top_k: int = 3,
                      use_api: bool = True,
                      model: str = None,
                      save_interval: int = 10) -> pd.DataFrame:
        """
        Batch annotate Excel files (single call, logic consistent with GPT.py)

        Excel column structure:
            Predict_cell_type (prediction result)
            Matched_Genes (database hit marker)
        """
        if output_file is None:
            base_name = os.path.splitext(input_file)[0]
            output_file = f"{base_name}_annotated.xlsx"

        if self.verbose:
            print(f"\n{'='*60}")
            print(f"Read file: {input_file}")

        try:
            df = pd.read_excel(input_file)
        except Exception as e:
            raise ValueError(f"Failed to read file: {e}")

        if tissue_col not in df.columns or genes_col not in df.columns:
            raise ValueError(f"File must contain '{tissue_col}' and '{genes_col}' columns")

        if self.verbose:
            print(f"A total of {len(df)} rows of data need to be processed")

        df['Predict_cell_type'] = ''
        df['Matched_Genes'] = ''

        total = len(df)
        success_count = 0
        fail_count = 0
        skip_count = 0

        if self.verbose:
            print(f"\n{'='*60}")
            print('Start batch annotation...')

        for idx, row in df.iterrows():
            tissue = str(row[tissue_col]).strip() if pd.notna(row[tissue_col]) else ''
            genes = str(row[genes_col]).strip() if pd.notna(row[genes_col]) else ''

            if not genes:
                df.at[idx, 'Predict_cell_type'] = 'null value'
                skip_count += 1
                if self.verbose:
                    print(f"[{idx+1}/{total}] Skip empty genes")
                continue

            try:
                candidates = self.retrieve_candidates(genes, tissue, top_k=top_k)

                if candidates:
                    result = self.generate_annotation(
                        genes, tissue, candidates,
                        use_api=use_api, model=model
                    )

                    df.at[idx, 'Predict_cell_type'] = result['cell_type']
                    df.at[idx, 'Matched_Genes'] = ', '.join(result['matched_genes'])
                    success_count += 1

                    if self.verbose:
                        tissue_display = tissue[:15] if tissue else 'not specified'
                        print(f"[{idx+1}/{total}] {tissue_display:15s} → {result['cell_type'][:40]}")
                else:
                    df.at[idx, 'Predict_cell_type'] = 'No candidate found'
                    fail_count += 1

                    if self.verbose:
                        tissue_display = tissue[:15] if tissue else 'not specified'
                        print(f"[{idx+1}/{total}] {tissue_display:15s} → No candidate found")

                if (idx + 1) % save_interval == 0:
                    df.to_excel(output_file, index=False)
                    if self.verbose:
                        print(f"   Progress saved ({idx+1}/{total})")

            except Exception as e:
                fail_count += 1
                df.at[idx, 'Predict_cell_type'] = f"ERROR: {str(e)[:50]}"
                if self.verbose:
                    print(f"[{idx+1}/{total}] Processing failed: {str(e)[:50]}")

            if use_api and idx < total - 1 and OPENAI_AVAILABLE:
                time.sleep(0.5)

        df.to_excel(output_file, index=False)

        if self.verbose:
            print(f"\n{'='*60}")
            print('Batch annotation completed!')
            print(f"Success: {success_count} | Skip: {skip_count} | Failure: {fail_count}")
            print(f"Result saved: {output_file}")
            print(f"{'='*60}\n")
        
        return df


def annotate_cell_type(genes: Union[str, List[str]],
                       tissue: str,
                       json_db_path: str = None,
                       top_k: int = 3,
                       use_api: bool = True,
                       api_key: str = None,
                       base_url: str = None,
                       model: str = "gpt-5") -> Tuple[Dict, List[Dict], Optional[Dict]]:
    """Convenience function: one-shot cell type annotation"""
    annotator = CellTypeAnnotator(
        json_db_path=json_db_path,
        api_key=api_key,
        base_url=base_url,
        default_model=model
    )

    return annotator.annotate_cell_type(
        genes=genes,
        tissue=tissue,
        top_k=top_k,
        use_api=use_api,
        model=model
    )


def batch_annotate_from_excel(input_file: str,
                              output_file: str = None,
                              json_db_path: str = None,
                              tissue_col: str = 'Tissue',
                              genes_col: str = 'Markers',
                              top_k: int = 3,
                              use_api: bool = True,
                              api_key: str = None,
                              base_url: str = None,
                              model: str = "gpt-5",
                              save_interval: int = 10) -> pd.DataFrame:
    """Convenient function: batch annotate Excel files"""
    annotator = CellTypeAnnotator(
        json_db_path=json_db_path,
        api_key=api_key,
        base_url=base_url,
        default_model=model
    )

    return annotator.batch_annotate(
        input_file=input_file,
        output_file=output_file,
        tissue_col=tissue_col,
        genes_col=genes_col,
        top_k=top_k,
        use_api=use_api,
        model=model,
        save_interval=save_interval
    )
