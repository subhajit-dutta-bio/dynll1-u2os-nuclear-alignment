"""New exploratory field descriptors and common-focal multiscale analysis."""
from pathlib import Path
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import json,itertools,argparse,hashlib
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import HuberRegressor

P=Path(__file__).resolve().parent;ROOT=P.parent;R=P/'results';R.mkdir(exist_ok=True)
SCALES=json.loads((ROOT/'phase5/jump_raw_policy.json').read_text())['source_control_scales']
WIDTH=1280;HEIGHT=1080
METRICS=['log_area','log_AR','solidity','log_density','scaled_NND','global_pair_order','position_axis_coupling','log_DNA_proxy']

def path_for(field):
    paths=[ROOT/'phase6/data/raw/geometry'/f'{field}.parquet',ROOT/'phase5/data/jump_raw/geometry'/f'{field}.parquet']
    return next(f for f in paths if f.exists())

def one(record):
    row=dict(record);field=row['field'];g=pd.read_parquet(path_for(field));s=SCALES[row['batch']];v=g.area.between(.1*s['area'],10*s['area'])&np.isfinite(g[['x','y','theta','ar']]).all(axis=1)&g.minor.gt(0)&g.x0.gt(0)&g.y0.gt(0)&g.x1.lt(WIDTH)&g.y1.lt(HEIGHT);g=g[v].copy();g['border']=np.minimum.reduce([g.x,g.y,WIDTH-1-g.x,HEIGHT-1-g.y]);direct=g[g.ar.ge(1.2)].copy();n=len(direct);rec={k:row[k] for k in ['field','batch','plate','well','site','role','gene','jcp2022']};rec.update({'valid_nuclei':len(g),'directors':n,'field_area':WIDTH*HEIGHT,'reference_area':s['area'],'reference_major':s['major']})
    if len(g)<7 or n<100:return rec
    xy=g[['x','y']].to_numpy();tree=cKDTree(xy);dist,ids=tree.query(xy,k=7);density=5/(np.pi*dist[:,6]**2);g['rho']=density;direct=g[g.ar.ge(1.2)].copy();didx=g.index.get_indexer(direct.index);z=np.exp(2j*direct.theta.to_numpy());n=len(z);nullz=(z.sum()-z)/(n-1);near=ids[didx,1];delta=xy[near]-direct[['x','y']].to_numpy();phi=np.arctan2(delta[:,1],delta[:,0]);w=np.exp(2j*phi)
    # regionprops theta is measured from the row axis. Displacement phi is
    # measured from x, so this new positional statistic needs an explicit
    # conversion. Pairwise nuclear angle differences are frame invariant.
    z_position=np.exp(2j*(np.pi/2-direct.theta.to_numpy()));null_position=(z_position.sum()-z_position)/(n-1);coupling=np.real(w*np.conj(z_position-null_position));valid=direct.border.to_numpy()>=np.maximum(3*s['major'],dist[didx,1]+s['major']);dna=g.loc[g.DNA_background_subtracted_sum.gt(0),'DNA_background_subtracted_sum']
    rec.update({'rho':float(np.median(density)),'log_area':float(np.log(g.area.median()/s['area'])),'log_AR':float(np.log(g.ar.median())),'solidity':float(g.solidity.mean()),'log_density':float(np.log(len(g)/(WIDTH*HEIGHT))),'scaled_NND':float(np.median(dist[:,1])*np.sqrt(len(g)/(WIDTH*HEIGHT))),'global_pair_order':float((abs(z.sum())**2-n)/(n*(n-1))),'position_axis_focals':int(valid.sum()),'position_axis_coupling':float(coupling[valid].mean()) if valid.sum()>=20 else np.nan,'log_DNA_proxy':float(np.log(dna.median())) if len(dna) else np.nan})
    if not row['multiscale']:return rec
    dtree=cKDTree(direct[['x','y']]);dxy=direct[['x','y']].to_numpy();base=1/np.sqrt(direct.rho.to_numpy());obs=[];counts=[]
    for c in [1.5,2.,2.5,3.]:
        lists=dtree.query_ball_point(dxy,c*base);values=[];nn=[]
        for i,js in enumerate(lists):
            js=np.asarray(js,int);js=js[js!=i];nn.append(len(js));values.append(float(np.real(z[i]*np.conj(z[js].mean()))) if len(js) else np.nan)
        obs.append(np.asarray(values));counts.append(np.asarray(nn))
    common=(np.vstack(counts)>=3).all(axis=0)&(direct.border.to_numpy()>=np.maximum(3*s['major'],3*base+s['major']));rec['multiscale_common_focals']=int(common.sum());null=np.real(z*np.conj(nullz))
    for c,arr in zip([1.5,2.,2.5,3.],obs):rec['radius_'+str(c).replace('.','p')]=float((arr[common]-null[common]).mean()) if common.sum()>=20 else np.nan
    return rec

def contrast(values,manifest):
    out=[]
    for row in manifest.itertuples():
        b=values[values.batch.eq(row.batch)];tg=b[b.plate.eq(row.target_plate)&b.well.eq(row.target_well)];source=b[b.plate.eq(row.target_plate)&b.well.isin(row.source_NT_wells.split(';'))];pos=b[b.plate.eq(row.negative_plate)&b.well.eq(row.negative_well)];anchor=b[b.plate.eq(row.negative_plate)&b.well.isin(row.anchor_wells.split(';'))];joined=pd.concat([tg,source,pos,anchor]);good=len(tg)==1 and len(source)==2 and len(pos)==1 and len(anchor)==2 and joined.sites.ge(6).all();effect=float(tg.value.mean()-source.value.mean()-pos.value.mean()+anchor.value.mean()) if good else np.nan
        out.append({'gene':row.gene,'panel':row.panel,'batch':row.batch,'evaluable':bool(good),'effect':effect,'target_sites':int(tg.sites.sum()),'source_sites_min':int(source.sites.min()) if len(source) else 0})
    return pd.DataFrame(out)

def fit_metric(d,metric,manifest):
    a=d[d[metric].notna()&d.directors.ge(100)].copy();ld=a.log_density;features={'noguide':a.jcp2022.eq('JCP2022_800001').astype(float)}
    if metric!='log_density':features.update({'ld':ld,'ld2':ld**2})
    if metric in ['global_pair_order','position_axis_coupling','log_DNA_proxy'] or metric.startswith('radius_'):features['area']=a.log_area
    if metric in ['global_pair_order','position_axis_coupling'] or metric.startswith('radius_'):features['AR']=a.log_AR
    x=pd.DataFrame(features,index=a.index);x=pd.concat([x,pd.get_dummies(a.batch,dtype=float)],axis=1);ctrl=a.role.ne('target');sc=StandardScaler().fit(x[ctrl]);xx=sc.transform(x);model=HuberRegressor(max_iter=3000).fit(xx[ctrl],a.loc[ctrl,metric]);assert model.n_iter_<3000;a['value']=a[metric]-model.predict(xx);w=a.groupby(['batch','plate','well','role'],as_index=False).agg(value=('value','mean'),sites=('field','nunique'));c=contrast(w,manifest);sd=float(a.loc[ctrl,'value'].std(ddof=1));c['control_field_residual_SD']=sd;c['effect_control_SD']=c.effect/sd;c['metric']=metric
    modelinfo={'metric':metric,'features':x.columns.tolist(),'control_fields':int(ctrl.sum()),'total_fields':len(a),'coef':model.coef_.tolist(),'intercept':float(model.intercept_),'feature_mean':sc.mean_.tolist(),'feature_scale':sc.scale_.tolist(),'control_residual_SD':sd,'iterations':int(model.n_iter_)}
    return a,c,modelinfo

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=6);ap.add_argument('--refresh',action='store_true');args=ap.parse_args();design=pd.read_csv(ROOT/'phase6/frozen_fields.csv').fillna({'gene':''});m=pd.read_csv(ROOT/'phase6/frozen_contrasts.csv').fillna({'RNF8_nearby_NT_wells':''});dfields=set(pd.read_csv(ROOT/'phase5/jump_raw_frozen_fields.csv').field) if (ROOT/'phase5/jump_raw_frozen_fields.csv').exists() else set()
    # The same270 original DYNLL1 contrast fields, identified from its manifest.
    dyn=m[m.gene.eq('DYNLL1')];wanted=set()
    for x in dyn.itertuples():
        for plate,wells in [(x.target_plate,[x.target_well]+x.source_NT_wells.split(';')),(x.negative_plate,[x.negative_well]+x.anchor_wells.split(';'))]:
            for well in wells:
                for site in range(9):wanted.add(f'{plate}__{well}__s{site}')
    assert len(wanted)==270;design['multiscale']=design.field.isin(wanted)
    cache=R/'field_descriptors.csv';provenance=R/'field_descriptor_provenance.json';algorithm_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if not args.refresh and cache.exists() and provenance.exists() and json.loads(provenance.read_text())['algorithm_sha256']==algorithm_hash:fields=pd.read_csv(cache).fillna({'gene':''});assert set(fields.field)==set(design.field)
    else:
        records=[]
        if args.workers==1:
            for i,record in enumerate(map(one,design.to_dict('records')),1):
                records.append(record)
                if i%300==0:print('Geometry summaries',i,len(design),flush=True)
        else:
            with ProcessPoolExecutor(max_workers=args.workers) as pool:
                for i,record in enumerate(pool.map(one,design.to_dict('records'),chunksize=15),1):
                    records.append(record)
                    if i%300==0:print('Geometry summaries',i,len(design),flush=True)
        fields=pd.DataFrame(records).sort_values('field');fields.to_csv(cache,index=False);provenance.write_text(json.dumps({'algorithm_sha256':algorithm_hash,'angle_frame':'regionprops row-axis angles converted to image x-axis only for displacement coupling and display. Pairwise nematic differences invariant.','all_geometry_fields':len(fields)},indent=2))
    results=[];models=[];fieldparts=[]
    for metric in METRICS+['radius_1p5','radius_2p0','radius_2p5','radius_3p0']:
        inputfields=fields[fields.field.isin(wanted)] if metric.startswith('radius_') else fields;manifest=dyn if metric.startswith('radius_') else m;a,c,model=fit_metric(inputfields,metric,manifest);results.append(c);models.append(model);fieldparts.append(a[['field','batch','plate','well','site','role','gene','value']].assign(metric=metric));print('Calibrated',metric,'evaluable',int(c.evaluable.sum()),'of',len(c),flush=True)
    contrasts=pd.concat(results,ignore_index=True);contrasts.to_csv(R/'new_descriptor_contrasts.csv',index=False);pd.concat(fieldparts,ignore_index=True).to_csv(R/'new_descriptor_adjusted_fields.csv',index=False);(R/'new_descriptor_models.json').write_text(json.dumps(models,indent=2));contrasts.groupby(['metric','gene'],as_index=False).agg(mean_effect=('effect','mean'),mean_effect_in_control_SD=('effect_control_SD','mean'),dates=('effect','count'),positive_dates=('effect',lambda x:int((x>0).sum()))).to_csv(R/'new_descriptor_gene_summary.csv',index=False)
    primary=pd.read_csv(ROOT/'phase6/results/adaptive_conditional_contrasts.csv');dv=primary[primary.gene.eq('DYNLL1')];pd.DataFrame([{'omitted_date':x.batch,'remaining_dates':len(dv)-1,'mean_effect':float(dv[dv.batch.ne(x.batch)].effect.mean())} for x in dv.itertuples()]).to_csv(R/'leave_one_date_out.csv',index=False)
    original=pd.read_csv(ROOT/'phase6/results/adaptive_conditional_fields.csv');original['site']=original.field.str.extract(r'__s(\d+)$').astype(int);site_results=[]
    for site in range(9):
        subset=original[original.site.ne(site)];w=subset.groupby(['batch','plate','well','role'],as_index=False).agg(value=('residual','mean'),sites=('field','nunique'));c=contrast(w,dyn);c['omitted_site']=site;site_results.append(c)
    pd.concat(site_results,ignore_index=True).to_csv(R/'leave_one_site_out.csv',index=False);summary={'fields':len(fields),'wells':int(design[['plate','well']].drop_duplicates().shape[0]),'valid_nuclei_all_geometry_fields':int(fields.valid_nuclei.sum()),'directors_all_geometry_fields':int(fields.directors.sum()),'descriptor_eligible_fields':int(fields.directors.ge(100).sum()),'eligible_directors':int(fields.loc[fields.directors.ge(100),'directors'].sum()),'low_director_yield_fields':int(fields.directors.lt(100).sum()),'metrics':METRICS,'multiscale_fields':len(wanted),'all_metric_gene_date_contrasts':len(contrasts),'unevaluable_contrasts':int((~contrasts.evaluable).sum()),'scope':'Exploratory field-level descriptors with existing segmentation and controls; not a replacement for original spatial inference or independent validation.'};(R/'spatial_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
