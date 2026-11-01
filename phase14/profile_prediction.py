"""Nested linear/nonlinear profile prediction on the fixed comparison panel."""
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ[k]='2'
import json,hashlib,itertools,argparse
import numpy as np,pandas as pd,pyarrow as pa,pyarrow.parquet as pq
from sklearn.metrics import r2_score,mean_absolute_error
P=Path(__file__).resolve().parent;ROOT=P.parent;R=P/'results/prediction';R.mkdir(parents=True,exist_ok=True);POL=json.loads((P/'prediction_policy.json').read_text());ALPHAS=np.asarray(POL['ridge_alpha'],float);GAMMAS=POL['RBF_gamma_multiplier']
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def contrast_inputs(join_effects=True):
    manifest=pd.read_csv(P/'frozen_contrasts.csv');needed=set()
    for row in manifest.itertuples():
        needed.update([(row.target_plate,row.target_well),(row.negative_plate,row.negative_well)]);needed.update((row.target_plate,w) for w in row.source_NT_wells.split(';'));needed.update((row.negative_plate,w) for w in row.anchor_wells.split(';'))
    cache=P/'data/prediction_wells.parquet'
    if not cache.exists():
        source=ROOT/'phase2/data/crispr_interpretable.parquet';parts=[]
        for batch in pq.ParquetFile(source).iter_batches(batch_size=1024):
            plates=batch.column(batch.schema.get_field_index('Metadata_Plate')).to_pylist();wells=batch.column(batch.schema.get_field_index('Metadata_Well')).to_pylist();ids=[i for i,key in enumerate(zip(plates,wells)) if key in needed]
            if ids:parts.append(pa.Table.from_batches([batch.take(pa.array(ids))]))
        pa.concat_tables(parts).to_pandas().to_parquet(cache,index=False);(P/'data/prediction_profile_provenance.json').write_text(json.dumps({'source':str(source.relative_to(ROOT)),'source_sha256':sha(source),'source_rows':51185,'subset_sha256':sha(cache),'required_wells':len(needed),'scope':'Original interpretable globally preprocessed profile stage, prior to proximity-bias correction'},indent=2))
    wells=pd.read_parquet(cache);assert not wells.duplicated(['Metadata_Plate','Metadata_Well']).any();features=[c for c in wells if not c.startswith('Metadata_')];wells=wells.set_index(['Metadata_Plate','Metadata_Well']);records=[];coverage=[]
    for row in manifest.itertuples():
        keys=[(row.target_plate,row.target_well)]+[(row.target_plate,w) for w in row.source_NT_wells.split(';')]+[(row.negative_plate,row.negative_well)]+[(row.negative_plate,w) for w in row.anchor_wells.split(';')];missing=[p+'/'+w for p,w in keys if (p,w) not in wells.index];coverage.append({'gene':row.gene,'batch':row.batch,'complete':not missing,'missing':';'.join(missing)})
        if missing:continue
        v=wells.loc[keys,features].to_numpy(float);delta=v[0]-v[1:3].mean(axis=0)-v[3]+v[4:6].mean(axis=0);records.append({'gene':row.gene,'batch':row.batch,**dict(zip(features,delta))})
    pd.DataFrame(coverage).to_csv(R/'profile_coverage.csv',index=False);d=pd.DataFrame(records);d.to_parquet(R/'feature_contrasts.parquet',index=False)
    if not join_effects:return d,features
    effects=pd.read_csv(P/'results/alignment_contrasts.csv');effects=effects[effects.model.eq('adaptive_conditional')&effects.effect.notna()];d=d.merge(effects[['gene','batch','effect']],on=['gene','batch'],validate='one_to_one');return d,features
def prepare(train,test):
    keep=np.isfinite(train).any(axis=0);a=train[:,keep].astype(float);b=test[:,keep].astype(float);a[~np.isfinite(a)]=np.nan;b[~np.isfinite(b)]=np.nan;median=np.nanmedian(a,axis=0);a=np.where(np.isnan(a),median,a);b=np.where(np.isnan(b),median,b);mean=a.mean(axis=0);scale=a.std(axis=0);keep=scale>1e-12;assert keep.any();return (a[:,keep]-mean[keep])/scale[keep],(b[:,keep]-mean[keep])/scale[keep]
def predict_grid(xtr,ytr,xte,kind):
    a,b=prepare(xtr,xte);center=ytr.mean();gram=a@a.T;cross=b@a.T;grids=[(gram,cross,None)]
    if kind=='RBF':
        dt=np.maximum((a*a).sum(axis=1)[:,None]+(a*a).sum(axis=1)[None,:]-2*gram,0);dv=np.maximum((b*b).sum(axis=1)[:,None]+(a*a).sum(axis=1)[None,:]-2*cross,0);grids=[(np.exp(-g*dt/a.shape[1]),np.exp(-g*dv/a.shape[1]),g) for g in GAMMAS]
    predictions=[];settings=[]
    for k,kcross,gamma in grids:
        eigen,vectors=np.linalg.eigh(k);eigen=np.maximum(eigen,0);z=vectors.T@(ytr-center);coef=vectors@(z[:,None]/(eigen[:,None]+ALPHAS[None,:]));p=kcross@coef+center
        for j,alpha in enumerate(ALPHAS):predictions.append(p[:,j]);settings.append({'alpha':float(alpha),'gamma_multiplier':gamma,'retained_features':a.shape[1]})
    return np.asarray(predictions),settings
def inner_splits(groups,grouping):
    unique=np.unique(groups)
    if grouping=='gene':
        rng=np.random.default_rng(2026103014);folds=np.array_split(rng.permutation(unique),min(5,len(unique)))
    else:folds=[np.asarray([x]) for x in unique]
    for fold in folds:
        test=np.isin(groups,fold);yield np.flatnonzero(~test),np.flatnonzero(test)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--prepare-only',action='store_true');args=ap.parse_args()
    if args.prepare_only:
        data,features=contrast_inputs(False);print('Prepared',len(data),'feature contrasts and',len(features),'features without reading new spatial outcomes');return
    data,features=contrast_inputs();assert len(data)>50;configs=[('full_ridge',features,'ridge'),('shape_count_ridge',[f for f in features if '_AreaShape_' in f or f=='Cells_Count_Count'],'ridge'),('without_neighbor_orientation_ridge',[f for f in features if '_Neighbors_' not in f and not f.endswith('_AreaShape_Orientation')],'ridge'),('full_RBF',features,'RBF')];outputs=[];folds=[];summaries=[]
    for name,cols,kind in configs:
        X=data[cols].to_numpy(float);Y=data.effect.to_numpy(float)
        for grouping in ['gene','batch']:
            groups=data[grouping].to_numpy();pred=np.full(len(Y),np.nan);baseline=np.full(len(Y),np.nan)
            for held in np.unique(groups):
                test=np.flatnonzero(groups==held);train=np.flatnonzero(groups!=held);total=None;count=0
                for ti,vi in inner_splits(groups[train],grouping):
                    estimates,settings=predict_grid(X[train][ti],Y[train][ti],X[train][vi],kind);sse=((estimates-Y[train][vi][None,:])**2).sum(axis=1);total=sse if total is None else total+sse;count+=len(vi)
                chosen=int(np.argmin(total/count));estimates,settings=predict_grid(X[train],Y[train],X[test],kind);pred[test]=estimates[chosen];baseline[test]=Y[train].mean();folds.append({'model':name,'grouping':grouping,'held_out':held,'training':len(train),'test':len(test),'inner_MSE':float(total[chosen]/count),**settings[chosen]})
            out=data[['gene','batch','effect']].copy();out['model']=name;out['grouping']=grouping;out['predicted']=pred;out['training_mean']=baseline;outputs.append(out);sse=((Y-pred)**2).sum();base=((Y-baseline)**2).sum();record={'model':name,'grouping':grouping,'observations':len(Y),'genes':data.gene.nunique(),'features':len(cols),'R2':float(r2_score(Y,pred)),'MAE':float(mean_absolute_error(Y,pred)),'training_mean_R2':float(r2_score(Y,baseline)),'SSE_improvement_vs_training_mean':float(1-sse/base)};summaries.append(record);print(record,flush=True)
    pd.concat(outputs,ignore_index=True).to_csv(R/'OOF_predictions.csv',index=False);pd.DataFrame(folds).to_csv(R/'folds.csv',index=False);pd.DataFrame(summaries).to_csv(R/'summary.csv',index=False);(P/'qa/prediction_provenance.json').write_text(json.dumps({'policy_sha256':sha(P/'prediction_policy.json'),'script_sha256':sha(Path(__file__)),'features':len(features),'observations':len(data),'all_models_reported':True,'scope':'Targeted-panel prediction, not a genome-wide irreducibility test'},indent=2))
if __name__=='__main__':main()
