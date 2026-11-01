"""Recompute every Phase11 analysis offline and compare all numeric outputs."""
from pathlib import Path
import hashlib,json,runpy,sys,socket
import numpy as np
import pandas as pd
import requests
P=Path(__file__).resolve().parent;ROOT=P.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def blocked(*a,**k):raise RuntimeError('Network disabled during reproducibility verification')
def main():
    oldfiles=[ROOT/'phase6/results/adaptive_conditional_contrasts.csv',ROOT/'phase6/results/module_tests.json',ROOT/'phase10/results/hTERT_candidate_audit.json']
    frozen={str(p.relative_to(ROOT)):sha(p) for p in oldfiles}
    before={p.name:pd.read_csv(p) for p in (P/'results').glob('*.csv*')}
    before_json={p.name:json.loads(p.read_text(encoding='utf-8')) for p in (P/'results').glob('*.json') if p.name!='field_descriptor_provenance.json'}
    requests.sessions.Session.request=blocked;socket.create_connection=blocked
    for name,args in [('spatial.py',['--workers','1','--refresh']),('context.py',[]),('prediction.py',[]),('angle_QA.py',[])]:
        sys.argv=[name]+args;runpy.run_path(str(P/name),run_name='__main__')
    checked=[]
    for name,df in before.items():
        after=pd.read_csv(P/'results'/name);pd.testing.assert_frame_equal(df,after,check_exact=False,atol=1e-12,rtol=1e-12);checked.append(name)
    for name,obj in before_json.items():assert obj==json.loads((P/'results'/name).read_text(encoding='utf-8')),name
    for name,digest in frozen.items():assert sha(ROOT/name)==digest,name
    reg=json.loads((P/'panel_registry.json').read_text(encoding='utf-8'));assert len(reg)==66
    assert len({(x['figure'],x['panel']) for x in reg})==66
    for x in reg:assert len(pd.read_csv(ROOT/x['source_data']))>0 and not x['independent_spatial_replication']
    c=pd.read_csv(P/'results/new_descriptor_contrasts.csv');assert len(c)==452 and c.evaluable.all()
    s=json.loads((P/'results/spatial_summary.json').read_text());assert s['eligible_directors']==615037 and s['descriptor_eligible_fields']==1807
    angle=json.loads((P/'qa/angle_convention.json').read_text());assert angle['all_pass']
    result={'all_pass':True,'network_disabled':True,'geometry_cache_bypassed':True,'csv_outputs_reproduced':len(checked),'json_outputs_reproduced':len(before_json),'numeric_tolerance':1e-12,'panel_sources_checked':len(reg),'metric_gene_date_contrasts':len(c),'unchanged_prior_inputs':frozen,'output_files':checked,'scope':'Recomputes Phase11 from retained geometry and Phase10 processed guide tables. Original segmentation and raw screen reconstruction have separate earlier verification records.'}
    (P/'qa/verification.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
