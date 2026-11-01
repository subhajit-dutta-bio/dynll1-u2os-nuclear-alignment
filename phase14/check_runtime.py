"""Accept optimization only on pixel-identical existing reference fields."""
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[k]='4'
import sys,json,time,hashlib
import torch,numpy as np,pandas as pd,tifffile
P=Path(__file__).resolve().parent;ROOT=P.parent;sys.path.insert(0,str(ROOT/'phase12'));from second_segmentation import infer
from cellpose import models
from graph_runtime import install
def main():
    (P/'qa').mkdir(exist_ok=True);policy=json.loads((ROOT/'phase12/frozen_policy.json').read_text());weights=Path(policy['model']['local_weights']);assert hashlib.sha256(weights.read_bytes()).hexdigest()==policy['model']['sha256'];torch.set_num_threads(4);torch.manual_seed(2026103012);np.random.seed(2026103012);m=models.CellposeModel(gpu=True,pretrained_model=str(weights),use_bfloat16=True)
    d=pd.read_csv(ROOT/'phase5/jump_raw_frozen_fields.csv');d=d[d.site.eq(0)&d.role.isin(['target','source_NT'])&d.well.isin(['A02','L20'])].sort_values(['batch','role']);assert len(d)==10
    inputs=[];records=[]
    for r in d.itertuples():
        f=ROOT/'phase5/data/jump_raw/images'/f'{r.field}.tif';a=tifffile.imread(f);mask,seconds=infer(m,a);old=np.load(ROOT/'phase12/results/cellpose_crispr/masks'/f'{r.field}.npz')['labels'];assert np.array_equal(mask,old);inputs.append((r.field,a,mask));records.append({'field':r.field,'native_seconds':seconds,'native_equals_Phase12':True});print('Native',r.field,round(seconds,3),flush=True)
    error=None
    try:
        install(m.net)
        for (field,a,native),rec in zip(inputs,records):
            mask,seconds=infer(m,a);rec.update(graph_seconds=seconds,identical_pixels=bool(np.array_equal(mask,native)),different_pixels=int((mask!=native).sum()));print('Graph',field,round(seconds,3),rec['identical_pixels'],flush=True)
    except Exception as e:error=repr(e)
    ok=error is None and all(r.get('identical_pixels',False) for r in records);native=float(np.median([r['native_seconds'] for r in records]));graph=float(np.median([r['graph_seconds'] for r in records])) if error is None else None
    report={'all_pixels_identical':ok,'accept_graph':bool(ok and graph<native),'native_median_seconds':native,'graph_median_seconds':graph,'error':error,'fields':records,'new_comparator_outcomes_used':False};(P/'qa/runtime_check.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='fields'},indent=2),flush=True)
if __name__=='__main__':main()
