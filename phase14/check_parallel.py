from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[k]='2'
import time,json,multiprocessing as mp
import numpy as np,pandas as pd,tifffile
P=Path(__file__).resolve().parent;ROOT=P.parent
def init(barrier):
    global MODEL,BARRIER
    import torch
    from cellpose import models
    torch.set_num_threads(2);torch.manual_seed(2026103012);np.random.seed(2026103012);pol=json.loads((ROOT/'phase12/frozen_policy.json').read_text());MODEL=models.CellposeModel(gpu=True,pretrained_model=pol['model']['local_weights'],use_bfloat16=True);BARRIER=barrier
def worker(index):
    d=pd.read_csv(ROOT/'phase5/jump_raw_frozen_fields.csv');d=d[d.site.eq(0)&d.role.isin(['target','source_NT'])&d.well.isin(['A02','L20'])].sort_values(['batch','role']);BARRIER.wait();start=time.monotonic();rows=[]
    for r in d.itertuples():
        raw=tifffile.imread(ROOT/'phase5/data/jump_raw/images'/f'{r.field}.tif');t=time.monotonic();mask,_,_=MODEL.eval(raw,batch_size=16,normalize=True,diameter=None,flow_threshold=.4,cellprob_threshold=0.,min_size=15,max_size_fraction=.4,resample=True,augment=False);old=np.load(ROOT/'phase12/results/cellpose_crispr/masks'/f'{r.field}.npz')['labels'];rows.append({'field':r.field,'identical_pixels':bool(np.array_equal(mask,old)),'seconds':time.monotonic()-t})
    return {'worker':index,'start':start,'end':time.monotonic(),'fields':rows}
def main():
    ctx=mp.get_context('spawn');barrier=ctx.Barrier(2)
    with ctx.Pool(2,initializer=init,initargs=(barrier,)) as pool:r=pool.map(worker,[0,1])
    seconds=max(x['end'] for x in r)-min(x['start'] for x in r);native=json.loads((P/'qa/runtime_check.json').read_text())['native_median_seconds'];ratio=20/seconds*native;identical=all(y['identical_pixels'] for x in r for y in x['fields']);record={'accept_two_workers':bool(identical and ratio>=1.1),'all_pixels_identical':identical,'elapsed_seconds_for_20_fields':seconds,'speed_ratio':ratio,'workers':r};(P/'qa/parallel_runtime_check.json').write_text(json.dumps(record,indent=2));print(json.dumps({k:v for k,v in record.items() if k!='workers'},indent=2),flush=True)
if __name__=='__main__':main()
