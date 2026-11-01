"""Frozen default Cellpose-SAM sensitivity; no model fitting or target tuning."""
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='4'
import sys,json,hashlib,time,zipfile,io,argparse
import numpy as np
import pandas as pd
import tifffile,torch
from PIL import Image
from scipy.optimize import linear_sum_assignment
from skimage import measure as sm
from cellpose import models
P=Path(__file__).resolve().parent;ROOT=P.parent
sys.path.insert(0,str(ROOT/'phase5'));from harmonize import measure,NAMES
sys.path.insert(0,str(ROOT/'phase7'));from benchmark import summary
POL=json.loads((P/'frozen_policy.json').read_text());SCALES=json.loads((ROOT/'phase5/jump_raw_policy.json').read_text())['source_control_scales']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def geom(mask,raw):
    g=pd.DataFrame(sm.regionprops_table(mask,intensity_image=raw,properties=['label','area','centroid','axis_major_length','axis_minor_length','orientation','bbox','solidity','intensity_mean']))
    g=g.rename(columns={'label':'nucleus_id','centroid-0':'y','centroid-1':'x','axis_major_length':'major','axis_minor_length':'minor','orientation':'theta','bbox-0':'y0','bbox-1':'x0','bbox-2':'y1','bbox-3':'x1'});g['ar']=g.major/g.minor.replace(0,np.nan);bg=float(np.median(raw[mask==0]));g['DNA_background_subtracted_sum']=(g.intensity_mean-bg)*g.area;g['DNA_background']=bg;return g
def infer(model,raw):
    t=time.monotonic();mask,_,_=model.eval(raw,batch_size=16,normalize=True,diameter=None,flow_threshold=.4,cellprob_threshold=0.,min_size=15,max_size_fraction=.4,resample=True,augment=False)
    return mask.astype(np.uint16),time.monotonic()-t
def benchmark_record(gt,pred,gg,gp,sc,row):
    ng,npred=int(gt.max()),int(pred.max());ov=np.bincount((gt.astype(np.int64)*(npred+1)+pred).ravel(),minlength=(ng+1)*(npred+1)).reshape(ng+1,npred+1);ag=np.bincount(gt.ravel(),minlength=ng+1)[1:];ap=np.bincount(pred.ravel(),minlength=npred+1)[1:];inter=ov[1:,1:];iou=inter/(ag[:,None]+ap[None,:]-inter);gi,pi=linear_sum_assignment(iou,maximize=True);ok=iou[gi,pi]>=.5;gi,pi=gi[ok],pi[ok]
    a=gg.copy();mapping={int(g)+1:int(p)+1 for g,p in zip(gi,pi)};a['pred_id']=pd.array(a.nucleus_id.map(mapping),dtype='Int64');pg=gp[['nucleus_id','ar','theta']].rename(columns={'nucleus_id':'pred_id','ar':'pred_ar','theta':'pred_theta'});pg.pred_id=pg.pred_id.astype('Int64');a=a.merge(pg,on='pred_id',how='left',validate='many_to_one');h,w=gt.shape
    a['manual_eligible']=a.area.between(.1*sc['area'],10*sc['area'])&a.ar.ge(1.2)&a.x0.gt(0)&a.y0.gt(0)&a.x1.lt(w)&a.y1.lt(h);a['usable_axis']=a.manual_eligible&a.pred_ar.ge(1.2)&a.pred_id.notna();delta=(a.theta-a.pred_theta+np.pi/2)%np.pi-np.pi/2;a['axial_agreement']=np.where(a.usable_axis,np.cos(2*delta),np.nan);a['absolute_axis_error_deg']=np.where(a.usable_axis,np.abs(np.rad2deg(delta)),np.nan);tp=len(gi);n=int(a.manual_eligible.sum());u=int(a.usable_axis.sum())
    rec={'field':row['field'],'split':row['split'],'manual_objects':ng,'predicted_objects':npred,'TP':tp,'FP':npred-tp,'FN':ng-tp,'F1':2*tp/(ng+npred) if ng+npred else None,'manual_eligible_axes':n,'usable_axes':u,'axis_coverage':u/n if n else None,'axial_mean':float(a.axial_agreement.mean()) if u else None,'axial_sum':float(a.axial_agreement.sum()),'axis_error_sum':float(a.absolute_axis_error_deg.sum()),'above_15deg':int(a.absolute_axis_error_deg.gt(15).sum()),'height':h,'width':w}
    for name,g in [('manual',gg),('automatic',gp)]:
        f=measure(g,w,h,sc['area'],sc['major'])
        for endpoint in NAMES:
            v=f[endpoint].dropna() if endpoint in f else [];rec[name+'_'+endpoint+'_focals']=len(v);rec[name+'_'+endpoint]=float(np.mean(v)) if len(v)>=20 else None
    return rec,a
def run_benchmark(model,limit):
    folder=P/'results/cellpose_benchmark'
    for child in ['fields','objects','geometry','masks']:(folder/child).mkdir(parents=True,exist_ok=True)
    plan=json.loads((ROOT/'phase7/benchmark_plan.json').read_text());sc={k:v*plan['pixel_ratio']**(2 if k=='area' else 1) for k,v in SCALES[plan['primary_setting']].items()};design=pd.read_csv(ROOT/'phase7/benchmark_fields.csv').sort_values('mask_file');design=design.head(limit) if limit else design
    iz=zipfile.ZipFile(ROOT/'phase7/data/bbbc039_images.zip');mz=zipfile.ZipFile(ROOT/'phase7/data/bbbc039_masks.zip')
    for i,row in enumerate(design.to_dict('records'),1):
        field=Path(row['mask_file']).stem;row['field']=field;out=folder/'fields'/f'{field}.json'
        if out.exists():continue
        raw=tifffile.imread(io.BytesIO(iz.read('images/'+field+'.tif')));gt=sm.label(np.asarray(Image.open(io.BytesIO(mz.read('masks/'+row['mask_file']))))[:,:,0],background=0,connectivity=2);pred,elapsed=infer(model,raw);gg=geom(gt,raw);gp=geom(pred,raw);rec,a=benchmark_record(gt,pred,gg,gp,sc,row);rec['inference_seconds']=elapsed;np.savez_compressed(folder/'masks'/f'{field}.npz',labels=pred);gp.to_parquet(folder/'geometry'/f'{field}.parquet',index=False);a.to_parquet(folder/'objects'/f'{field}.parquet',index=False);out.write_text(json.dumps(rec,indent=2,allow_nan=False));print('Benchmark',i,len(design),field,'seconds',round(elapsed,2),'F1',rec['F1'],flush=True)
    frame=pd.DataFrame([json.loads(x.read_text()) for x in sorted((folder/'fields').glob('*.json'))]);frame.to_csv(folder/'all_fields.csv',index=False)
    if not limit:
        rng=np.random.default_rng(2026103012);result={split:summary(g,rng) for split,g in frame.groupby('split')};result['scope']='Reference diagnostic for pretrained default model. Prior training overlap with BBBC039 is unresolved; these are not guaranteed unseen images.';(folder/'summary.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
def run_crispr(model,limit,all_fields):
    folder=P/'results/cellpose_crispr'
    for child in ['fields','geometry','masks','focals']:(folder/child).mkdir(parents=True,exist_ok=True)
    design=pd.read_csv(ROOT/('phase6/frozen_fields.csv' if all_fields else 'phase5/jump_raw_frozen_fields.csv')).sort_values('field');design=design.head(limit) if limit else design
    for i,row in enumerate(design.itertuples(),1):
        out=folder/'fields'/f'{row.field}.json'
        if out.exists():continue
        rawfile=next(x for x in [ROOT/'phase6/data/raw/images'/f'{row.field}.tif',ROOT/'phase5/data/jump_raw/images'/f'{row.field}.tif'] if x.exists());raw=tifffile.imread(rawfile);pred,elapsed=infer(model,raw);g=geom(pred,raw);g.to_parquet(folder/'geometry'/f'{row.field}.parquet',index=False);np.savez_compressed(folder/'masks'/f'{row.field}.npz',labels=pred);sc=SCALES[row.batch];f=measure(g,raw.shape[1],raw.shape[0],sc['area'],sc['major'])
        for k in ['field','batch','plate','well','site','role','jcp2022']:f[k]=getattr(row,k)
        if all_fields:f['gene']=getattr(row,'gene','')
        f.to_parquet(folder/'focals'/f'{row.field}.parquet',index=False);rec={k:getattr(row,k) for k in ['field','batch','plate','well','site','role','jcp2022']};rec.update({'objects':len(g),'focals':len(f),'inference_seconds':elapsed,'raw_sha256':sha(rawfile)});out.write_text(json.dumps(rec,indent=2));
        if i%5==0:print('CRISPR',i,len(design),row.field,'seconds',round(elapsed,2),'objects',len(g),flush=True)
    pd.DataFrame([json.loads(x.read_text()) for x in sorted((folder/'fields').glob('*.json'))]).to_csv(folder/'all_fields.csv',index=False)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['benchmark','crispr','both']);ap.add_argument('--limit',type=int);ap.add_argument('--all',action='store_true');a=ap.parse_args();torch.set_num_threads(4);torch.manual_seed(2026103012);np.random.seed(2026103012);assert torch.cuda.is_available();assert sha(Path(POL['model']['local_weights']))==POL['model']['sha256'];model=models.CellposeModel(gpu=True,pretrained_model=POL['model']['local_weights'],use_bfloat16=True);print('Loaded frozen model',torch.cuda.get_device_name(),flush=True)
    if a.mode in ['benchmark','both']:run_benchmark(model,a.limit)
    if a.mode in ['crispr','both']:run_crispr(model,a.limit,a.all)
if __name__=='__main__':main()
