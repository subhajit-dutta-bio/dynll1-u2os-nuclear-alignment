from pathlib import Path
import hashlib,json,runpy,socket,contextlib
import numpy as np
import pandas as pd
import requests
P=Path(__file__).resolve().parent;R=P/'results';Q=P/'qa';Q.mkdir(exist_ok=True)
def block(*a,**k):raise RuntimeError('Networking disabled for reproduction')
def main():
    # Independent direct calculation of all hTERT target-guide ratios from XLSX.
    sheets={n:pd.read_excel(P.parent/'phase9/data/GSE226994_MAGeCK_count_normalized.xlsx',sheet_name=n) for n in ['U2OS-rep1.count_normalized','U2OS-rep2.count_normalized','KO-rep1.count_normalized','KO-rep2.count_normalized']};audit=pd.read_csv(R/'hTERT_DYNLL1_DYNLL2_reagent_and_interaction_audit.csv');largest=0.
    for row in audit.itertuples():
        for context,prefix in [('WT','U2OS'),('KO','KO')]:
            for i in [1,2]:
                tab=sheets[f'{prefix}-rep{i}.count_normalized'].set_index('sgRNA');x=tab.loc[row.guide];value=float(np.log2((x[f'{prefix}-d21-rep{i}']+1)/(x[f'{prefix}-d0-rep{i}']+1)));delta=abs(value-getattr(row,f'LFC_r{i}_{context}'));assert delta<1e-12;largest=max(largest,delta)
    outputs={f.name:pd.read_csv(f) for f in R.glob('*.csv')};outputs.update({f.name:pd.read_csv(f) for f in R.glob('*gene_effects.csv.gz')});outputs['hTERT_all_gene_interactions.csv.gz']=pd.read_csv(R/'hTERT_all_gene_interactions.csv.gz');before={f.name:pd.read_parquet(f) for f in (R/'guides').glob('*.parquet')}
    summaries={name:json.loads((R/name).read_text(encoding='utf-8')) for name in ['summary.json','hTERT_candidate_audit.json','APMS_deep_summary.json','published_score_scale_clarification.json']}
    requests.sessions.Session.request=block;socket.create_connection=block
    with (Q/'offline_run.log').open('w',encoding='utf-8') as log,contextlib.redirect_stdout(log):
        runpy.run_path(str(P/'analyze.py'),run_name='offline')['main']();runpy.run_path(str(P/'followup.py'),run_name='offline_followup')['main']()
    for name,old in outputs.items():pd.testing.assert_frame_equal(old,pd.read_csv(R/name),check_exact=False,atol=1e-12,rtol=0)
    for name,old in before.items():pd.testing.assert_frame_equal(old,pd.read_parquet(R/'guides'/name),check_exact=True)
    for name,old in summaries.items():assert old==json.loads((R/name).read_text(encoding='utf-8')),name
    s=summaries['hTERT_candidate_audit.json'];assert s['nonoverlapping_guides']==1 and s['overlapping_JUMP_intervals']==3 and s['exact_JUMP_shared_sequences']==1;assert s['nonoverlapping_interactions'][0]['delta_r1']>0 and s['nonoverlapping_interactions'][0]['delta_r2']<0
    assert len(before)==13
    result={'all_pass':True,'HTTP_disabled':True,'connection_creation_disabled':True,'screen_contrasts_recomputed':13,'guide_tables_reproduced_exactly':len(before),'CSV_tables_reproduced':len(outputs),'JSON_summaries_reproduced_exactly':len(summaries),'independent_raw_XLSX_target_ratio_maximum_difference':largest,'official_guide_identity_and_overlap_checks':True,'scope':'Computational reproduction and extraction checks, not independent biological validation.'};(Q/'verification.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':main()
