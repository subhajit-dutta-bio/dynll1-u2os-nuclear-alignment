"""Verify original files, target extraction and offline deterministic reproduction."""
from pathlib import Path
import json,hashlib,runpy,requests,socket,contextlib
import pandas as pd
import numpy as np
P=Path(__file__).resolve().parent;D=P/'data';R=P/'results';Q=P/'qa';Q.mkdir(exist_ok=True)
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def blocked(*a,**k):raise RuntimeError('Network disabled in verification')
def main():
    # Every retained target count row must exactly match the original full table.
    specs=[('GSE232592_CRISPR_p53_counts.txt.gz','\t','Gene'),('GSE290538_U2OS_HRE_Count_table.txt.gz','\t','Gene'),('GSE283147_combined-guide-level-counts.csv.gz',',','GENE'),('GSE233195_JSMO.count.txt.gz','\t','Gene')];counts={}
    for file,sep,gene in specs:
        a=pd.read_csv(D/file,sep=sep,low_memory=False);expected=a[a[gene].eq('DYNLL1')].reset_index(drop=True);saved=pd.read_csv(R/(file.split('_')[0]+'_DYNLL1_guide_counts.csv'));pd.testing.assert_frame_equal(expected,saved,check_dtype=False);counts[file.split('_')[0]]=len(saved)
    published=pd.read_csv(R/'GSE232592_DYNLL1_published_gene_results.csv');fdr=[float(v) for c in published if c.endswith('_FDR') for v in published[c].dropna()];assert len(fdr)==8 and min(fdr)>.05
    for record in json.loads((P/'depmap_matrix_provenance.json').read_text()):assert sha(P/record['retained_file'].replace('\\','/'))==record['retained_sha256']
    guide=pd.read_csv(R/'public_CRISPR_guide_overlap_with_JUMP.csv');assert len(guide[(guide.source=='GSE290538')&guide.identical_JUMP_sequence.notna()])==1;assert len(guide[(guide.source=='DepMap24Q4')&guide.used_in_primary_score.eq(True)&guide.identical_JUMP_sequence.notna()])==1
    # No external requests are used by the retained-data extraction.
    names=['GSE232592_DYNLL1_published_gene_results.csv','GSE232592_DYNLL1_guide_counts.csv','GSE290538_DYNLL1_guide_counts.csv','GSE283147_DYNLL1_guide_counts.csv','GSE233195_DYNLL1_guide_counts.csv','GSE226994_DYNLL1_guide_counts.csv','verified_U2OS_CRISPR_sources.csv','DepMap24Q4_DYNLL1_Avana_guides.csv','DepMap24Q4_U2OS_DYNLL1.json','public_CRISPR_guide_overlap_with_JUMP.csv','U2OS_APMS_DYNLL1_DYNLL2_associations.csv','U2OS_LC8_localization_attributes.json','U2OS_APMS_summary.json','qChIP_U2OS_DYNLL1_DYNLL2_protein_group.csv','qChIP_U2OS_LC8_KO_control_rows.csv','qChIP_U2OS_control_summary.json','LINCS_DYNLL1_instances_all_cell_metadata.csv','LINCS_DYNLL1_signatures_all_cell_metadata.csv','perturbation_catalog_scope.json','GEO_U2OS_CRISPR_series.csv','EuropePMC_DYNLL1_U2OS_candidates.csv','summary.json','U2OS_DYNLL1_validation_sources_ranked.csv','search_scope.json'];outputs=[R/name for name in names];before={f.name:sha(f) for f in outputs};requests.sessions.Session.request=blocked;socket.create_connection=blocked
    with (Q/'offline_reproduction.log').open('w',encoding='utf-8') as out,contextlib.redirect_stdout(out):
        runpy.run_path(str(P/'analyze_sources.py'),run_name='source_audit_verification')['main']();runpy.run_path(str(P/'curate_evidence.py'),run_name='curated_source_verification')['main']()
    after={f.name:sha(f) for f in outputs};assert before==after
    result={'all_pass':True,'DYNLL1_guide_count_rows_exact':counts,'published_GSE232592_FDR_values_checked':len(fdr),'minimum_published_FDR':min(fdr),'DepMap_retained_row_hashes_match':True,'guide_overlap_counts_match':True,'result_files_recomputed_and_reproduced_byte_identically':len(outputs),'verified_output_files':names,'HTTP_disabled':True,'connection_creation_disabled':True,'maximum_numeric_difference':0.0,'scope':'The 24 named quantitative and eligibility exports were recomputed. Additional exploratory/source-index exports are retained but are not counted as recomputed outputs. No new statistical inference or independent spatial validation.'};(Q/'verification.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':main()
