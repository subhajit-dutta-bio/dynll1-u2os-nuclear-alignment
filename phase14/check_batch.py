from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[k]='4'
import json,time,hashlib
import torch,numpy as np,pandas as pd,tifffile
from cellpose import models
P=Path(__file__).resolve().parent;ROOT=P.parent
def main():
    policy=json.loads((ROOT/'phase12/frozen_policy.json').read_text());weights=Path(policy['model']['local_weights']);assert hashlib.sha256(weights.read_bytes()).hexdigest()==policy['model']['sha256'];torch.set_num_threads(4);torch.manual_seed(2026103012);np.random.seed(2026103012);m=models.CellposeModel(gpu=True,pretrained_model=str(weights),use_bfloat16=True)
    d=pd.read_csv(ROOT/'phase5/jump_raw_frozen_fields.csv');d=d[d.site.eq(0)&d.role.isin(['target','source_NT'])&d.well.isin(['A02','L20'])].sort_values(['batch','role']);assert len(d)==10;records=[];error=None
    try:
        for r in d.itertuples():
            raw=tifffile.imread(ROOT/'phase5/data/jump_raw/images'/f'{r.field}.tif');t=time.monotonic();mask,_,_=m.eval(raw,batch_size=64,normalize=True,diameter=None,flow_threshold=.4,cellprob_threshold=0.,min_size=15,max_size_fraction=.4,resample=True,augment=False);elapsed=time.monotonic()-t;old=np.load(ROOT/'phase12/results/cellpose_crispr/masks'/f'{r.field}.npz')['labels'];rec={'field':r.field,'seconds':elapsed,'identical_pixels':bool(np.array_equal(mask,old)),'different_pixels':int((mask!=old).sum())};records.append(rec);print(rec,flush=True)
    except Exception as e:error=repr(e)
    baseline=json.loads((P/'qa/runtime_check.json').read_text())['native_median_seconds'];median=float(np.median([r['seconds'] for r in records])) if records else None;ok=error is None and len(records)==10 and all(r['identical_pixels'] for r in records);report={'accept_batch64':bool(ok and median<baseline),'all_pixels_identical':ok,'batch64_median_seconds':median,'native_batch16_median_seconds':baseline,'peak_GPU_allocated_bytes':torch.cuda.max_memory_allocated(),'error':error,'fields':records};(P/'qa/batch_runtime_check.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='fields'},indent=2),flush=True)
if __name__=='__main__':main()
