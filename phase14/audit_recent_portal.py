"""Eligibility-only lookup in a newly checked public perturbation registry.

No expression direction or pathway is selected. Candidate signatures require
an exact target and separately verified U2OS context before any analysis.
"""
from pathlib import Path
import json,hashlib,re,sys,requests,time
import pandas as pd
P=Path(__file__).resolve().parent;D=P/'data/recent_portal';D.mkdir(parents=True,exist_ok=True)
NAMES=['lincs-l1000-xpr-cell_line_metadata.parquet','lincs-l1000-xpr-gene_target_metadata.parquet','perturb-atlas-human-cell_line_metadata.parquet','perturb-atlas-human-gene_target_metadata.parquet','rummageo-gene-gene_target_metadata.parquet']
def main():
    sys.stdout.reconfigure(encoding='utf-8');records=[]
    for name in NAMES:
        url='https://perturbseqr.s3.amazonaws.com/'+name;p=D/name
        if not p.exists():r=requests.get(url,timeout=90);r.raise_for_status();p.write_bytes(r.content)
        d=pd.read_parquet(p);text=d.astype(str).agg('|'.join,axis=1);target=d[text.str.contains(r'(?i)\bDYNLL1\b',regex=True)];u2os=d[text.str.contains(r'(?i)\bU[- ]?2[- ]?OS\b',regex=True)];q={'file':name,'url':url,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size,'rows':len(d),'columns':list(d),'DYNLL1_rows':target.to_dict('records'),'U2OS_rows':u2os.to_dict('records')};records.append(q);print(name,len(d),'DYNLL1',len(target),'U2OS',len(u2os),d.columns.tolist(),flush=True)
        if len(target):print(target.to_string(index=False),flush=True)
        if len(u2os):print(u2os.to_string(index=False),flush=True)
    (D/'metadata_audit.json').write_text(json.dumps({'checked_UTC':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'registry':'https://perturbseqr.maayanlab.cloud/download','purpose':'Eligibility discovery only. Marginal target and cell-line coverage does not establish a joint U2OS DYNLL1 signature. No expression outcome or pathway selected.','records':records},indent=2,default=str),encoding='utf-8')
if __name__=='__main__':main()
