"""Repeat exact native inference on outcome-independent fixed example fields."""
from pathlib import Path
import json,hashlib,os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import numpy as np,pandas as pd,tifffile
P=Path(__file__).resolve().parent;ROOT=P.parent;R=P/'results'
def selection():
    d=pd.read_csv(P/'frozen_fields.csv').fillna({'gene':'','panel':''});rows=[]
    for key,g in d[d.role.eq('target')].groupby('panel'):rows.append(dict(g.sort_values('field').iloc[0],check_group='target/'+key))
    for key,g in d[d.role.ne('target')].groupby('role'):rows.append(dict(g.sort_values('field').iloc[0],check_group=key))
    return pd.DataFrame(rows)
def main():
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--plan-only',action='store_true');ap.add_argument('--weights',type=Path,help='Exact archived checkpoint; SHA-256 must match');a=ap.parse_args();plan=P/'pixel_repeat_plan.csv'
    if a.plan_only:selection().to_csv(plan,index=False);print('Saved fixed first-lexicographic group/role repeat plan');return
    assert json.loads((P/'qa/analysis_provenance.json').read_text())['all_fields_completed'];d=pd.read_csv(plan).fillna({'gene':'','panel':''});assert d.field.tolist()==selection().field.tolist()
    import torch
    from cellpose import models
    torch.set_num_threads(4);torch.manual_seed(2026103012);np.random.seed(2026103012);policy=json.loads((ROOT/'phase12/frozen_policy.json').read_text());weights=a.weights if a.weights is not None else Path(policy['model']['local_weights']);assert hashlib.sha256(weights.read_bytes()).hexdigest()==policy['model']['sha256'];m=models.CellposeModel(gpu=True,pretrained_model=str(weights),use_bfloat16=True);records=[]
    for row in d.itertuples():
        raw=tifffile.imread(P/'data/images'/f'{row.field}.tif');reference=np.load(R/'masks'/f'{row.field}.npz')['labels'];mask,_,_=m.eval(raw,batch_size=16,normalize=True,diameter=None,flow_threshold=.4,cellprob_threshold=0.,min_size=15,max_size_fraction=.4,resample=True,augment=False);rec={'field':row.field,'check_group':row.check_group,'raw_sha256':hashlib.sha256((P/'data/images'/f'{row.field}.tif').read_bytes()).hexdigest(),'equal_pixels':bool(np.array_equal(reference,mask)),'different_pixels':int((reference!=mask).sum()),'objects':int(mask.max())};records.append(rec);print(rec,flush=True)
    out={'all_pass':all(x['equal_pixels'] for x in records),'plan_sha256':hashlib.sha256(plan.read_bytes()).hexdigest(),'fields':records,'scope':'Exact inference repeatability only; no additional biological replicate or independent reagent.'};(P/'qa/pixel_repeats.json').write_text(json.dumps(out,indent=2));assert out['all_pass']
if __name__=='__main__':main()
