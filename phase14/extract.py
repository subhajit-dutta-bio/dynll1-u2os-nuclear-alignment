"""Frozen native Cellpose inference with parallel geometry measurement."""
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import sys,json,hashlib,time,argparse,importlib.util
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
import numpy as np,pandas as pd,tifffile
from skimage import measure as sm
P=Path(__file__).resolve().parent;ROOT=P.parent;R=P/'results';SCALES=json.loads((ROOT/'phase5/jump_raw_policy.json').read_text())['source_control_scales']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def raw_path(field):
    return next((f for f in [P/'data/images'/f'{field}.tif',ROOT/'phase6/data/raw/images'/f'{field}.tif',ROOT/'phase5/data/jump_raw/images'/f'{field}.tif'] if f.exists()),None)
def init_worker():
    global measure,sp
    (ROOT/'phase5/results').mkdir(parents=True,exist_ok=True);sys.path.insert(0,str(ROOT/'phase5'));from harmonize import measure
    spec=importlib.util.spec_from_file_location('spatial_frozen',ROOT/'phase11/spatial.py');sp=importlib.util.module_from_spec(spec);spec.loader.exec_module(sp);sp.path_for=lambda field:R/'geometry'/f'{field}.parquet'
def postprocess(row,mask,seconds,reused):
    raw=tifffile.imread(raw_path(row['field']));assert raw.shape==(1080,1280)
    g=pd.DataFrame(sm.regionprops_table(mask,intensity_image=raw,properties=['label','area','centroid','axis_major_length','axis_minor_length','orientation','bbox','solidity','intensity_mean']))
    g=g.rename(columns={'label':'nucleus_id','centroid-0':'y','centroid-1':'x','axis_major_length':'major','axis_minor_length':'minor','orientation':'theta','bbox-0':'y0','bbox-1':'x0','bbox-2':'y1','bbox-3':'x1'});g['ar']=g.major/g.minor.replace(0,np.nan);bg=float(np.median(raw[mask==0]));g['DNA_background_subtracted_sum']=(g.intensity_mean-bg)*g.area;g['DNA_background']=bg;g.to_parquet(R/'geometry'/f"{row['field']}.parquet",index=False);np.savez_compressed(R/'masks'/f"{row['field']}.npz",labels=mask.astype(np.uint16))
    sc=SCALES[row['batch']];f=measure(g,1280,1080,sc['area'],sc['major'])
    for k in ['field','batch','plate','well','site','role','jcp2022','gene','panel']:f[k]=row[k]
    f.to_parquet(R/'focals'/f"{row['field']}.parquet",index=False);rec=sp.one(dict(row,multiscale=False));rec.update(objects=len(g),focals=len(f),primary_focals=int(f.adaptive_conditional.notna().sum()) if 'adaptive_conditional' in f else 0,inference_seconds=seconds,reused_Phase12_mask=reused,raw_sha256=sha(raw_path(row['field'])));(R/'fields'/f"{row['field']}.json").write_text(json.dumps(rec,indent=2,allow_nan=True));return row['field']
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=6);ap.add_argument('--weights',type=Path,help='Path to the exact archived checkpoint; SHA-256 must match');a=ap.parse_args();plan=json.loads((P/'frozen_plan.json').read_text())
    for name,digest in plan['hashes'].items():assert sha(P/name)==digest,name
    for name,digest in plan['dependencies'].items():assert sha(ROOT/name)==digest,name
    for child in ['fields','geometry','masks','focals']:(R/child).mkdir(parents=True,exist_ok=True)
    runtime=json.loads((P/'qa/batch_runtime_check.json').read_text());batch_size=64 if runtime['accept_batch64'] else 16
    import torch
    from cellpose import models
    torch.set_num_threads(4);torch.manual_seed(2026103012);np.random.seed(2026103012);policy=json.loads((ROOT/'phase12/frozen_policy.json').read_text());weights=a.weights if a.weights is not None else Path(policy['model']['local_weights']);assert sha(weights)==policy['model']['sha256'];assert torch.cuda.is_available();model=models.CellposeModel(gpu=True,pretrained_model=str(weights),use_bfloat16=True)
    (P/'qa/extraction_runtime.json').write_text(json.dumps({'batch_size':batch_size,'checkpoint_sha256':sha(weights),'device':torch.cuda.get_device_name(),'torch':torch.__version__,'graph_used':False,'process_id':os.getpid()},indent=2));d=pd.read_csv(P/'frozen_fields.csv').fillna('');complete=sum((R/'fields'/f'{f}.json').exists() for f in d.field);pending=d[~d.field.map(lambda f:(R/'fields'/f'{f}.json').exists())];errors=[];start=time.monotonic();jobs={};inferred=0;cached=0;print('Starting',len(d),'fields; completed',complete,'batch',batch_size,flush=True)
    def collect(block):
        nonlocal complete,jobs
        if not jobs:return
        done,_=wait(jobs,return_when=FIRST_COMPLETED,timeout=None if block else 0)
        for fut in done:
            field=jobs.pop(fut)
            try:fut.result();complete+=1
            except Exception as e:errors.append({'field':field,'error':repr(e)})
    with ProcessPoolExecutor(max_workers=a.workers,initializer=init_worker) as pool:
        for row in pending.to_dict('records'):
            try:
                old=ROOT/'phase12/results/cellpose_crispr/masks'/f"{row['field']}.npz";reused=old.exists()
                if reused:mask=np.load(old)['labels'];seconds=0;cached+=1
                else:
                    while not (P/'data/provenance'/f"{row['field']}.json").exists():
                        collect(False);time.sleep(.5)
                    raw=tifffile.imread(raw_path(row['field']));assert raw.shape==(1080,1280);t=time.monotonic();mask,_,_=model.eval(raw,batch_size=batch_size,normalize=True,diameter=None,flow_threshold=.4,cellprob_threshold=0.,min_size=15,max_size_fraction=.4,resample=True,augment=False);seconds=time.monotonic()-t;inferred+=1
                jobs[pool.submit(postprocess,row,mask,seconds,reused)]=row['field'];collect(False)
                while len(jobs)>=4*a.workers:collect(True)
                if (inferred+cached)%25==0:
                    status={'completed':complete,'planned':len(d),'GPU_inferred':inferred,'cached_masks':cached,'queued_geometry':len(jobs),'technical_errors':len(errors),'elapsed_seconds':time.monotonic()-start};(P/'qa/extraction_progress.json').write_text(json.dumps(status,indent=2));print(json.dumps(status),flush=True)
            except Exception as e:errors.append({'field':row['field'],'error':repr(e)});print('ERROR',row['field'],repr(e),flush=True)
        while jobs:collect(True)
    (P/'qa/extraction_errors.json').write_text(json.dumps(errors,indent=2));assert not errors,errors[:5];records=pd.DataFrame([json.loads((R/'fields'/f'{f}.json').read_text()) for f in d.field]);assert len(records)==len(d);records.to_csv(R/'field_descriptors.csv',index=False);print('Extraction complete',len(records),'native fields',flush=True)
if __name__=='__main__':main()
