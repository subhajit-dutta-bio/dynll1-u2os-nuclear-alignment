"""Audit covariate distances from the unchanged target-focal matching rule."""
from pathlib import Path
import os
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[name]='2'
import json,hashlib
import numpy as np,pandas as pd
from scipy.spatial import cKDTree
P=Path(__file__).resolve().parent;R=P/'results'

def main():
    a=pd.read_parquet(R/'analysis_focal_inputs.parquet');a=a[a.adaptive_conditional.notna()].copy();counts=a.groupby('field',observed=True).size();a=a[a.field.isin(counts[counts.ge(20)].index)]
    manifest=pd.read_csv(P/'frozen_contrasts.csv');ret=pd.read_csv(R/'matching_retention.csv').set_index(['gene','batch']);parts=[]
    columns=['local_density','ar','area'];calipers=np.array([.2,.2,.3])
    for (batch,noguide),g in a.groupby(['batch',a.jcp2022.eq('JCP2022_800001')],observed=True):
        ctrl=g[g.role.ne('target')];target=g[g.role.eq('target')]
        if not len(target):continue
        cx=np.log(ctrl[columns].to_numpy(float))/calipers;gx=np.log(target[columns].to_numpy(float))/calipers;tree=cKDTree(cx);cf=ctrl.field.to_numpy();tf=target.field.to_numpy();diff=np.full((len(target),3),np.nan);matched=np.full((len(target),3),np.nan)
        for start in range(0,len(target),4096):
            stop=min(start+4096,len(target));ids=tree.query(gx[start:stop],k=64)[1]
            for j,js in enumerate(ids,start):
                js=js[(cf[js]!=tf[j])&(np.abs(cx[js]-gx[j]).max(axis=1)<=1)][:3]
                if len(js)==3:matched[j]=(cx[js]*calipers).mean(axis=0);diff[j]=gx[j]*calipers-matched[j]
        base=target[['gene','batch']].reset_index(drop=True).copy();base['retained']=np.isfinite(diff).all(axis=1)
        for k,column in enumerate(columns):base['target_'+column]=gx[:,k]*calipers[k];base['matched_'+column]=matched[:,k];base['difference_'+column]=diff[:,k]
        for (gene,date),q in base.groupby(['gene','batch'],observed=True):
            kept=q[q.retained];record={'gene':gene,'batch':date,'primary_focals':len(q),'matched_focals':len(kept)}
            assert len(q)==ret.loc[(gene,date),'primary_focals'] and len(kept)==ret.loc[(gene,date),'matched_focals']
            for k,column in enumerate(columns):
                values=kept['difference_'+column];absolute=values.abs();record.update({column+'_target_mean_log':kept['target_'+column].mean(),column+'_matched_control_mean_log':kept['matched_'+column].mean(),column+'_mean_paired_log_difference':values.mean(),column+'_median_absolute_log_difference':absolute.median(),column+'_p95_absolute_log_difference':absolute.quantile(.95),column+'_maximum_absolute_log_difference':absolute.max()})
                assert not len(values) or absolute.max()<=calipers[k]+1e-12
            parts.append(record)
        print('Audited unchanged target matching',batch,len(target),flush=True)
    result=manifest[['gene','panel','batch']].merge(pd.DataFrame(parts),on=['gene','batch'],how='left',validate='one_to_one');result[['primary_focals','matched_focals']]=result[['primary_focals','matched_focals']].fillna(0).astype(int);assert len(result)==335
    result.to_csv(R/'matching_covariate_balance.csv',index=False)
    record={'all_pass':True,'rows':len(result),'exact_original_retention_reproduces':True,'difference_definition':'Target log covariate minus the mean log covariate of its three matched controls; absolute-distance summaries use these target-specific differences.','unchanged_calipers':dict(zip(columns,calipers)),'scope':'Post-outcome diagnostic of the existing matching estimator, prompted by unmatched gene-level area/AR means. No rematching rule, outcome, eligibility threshold or primary estimate was changed.','focal_input_sha256':hashlib.sha256((R/'analysis_focal_inputs.parquet').read_bytes()).hexdigest()};(P/'qa/matching_covariate_balance.json').write_text(json.dumps(record,indent=2));print('Matching covariate audit complete; all335 target/date rows retained')
if __name__=='__main__':main()
