"""Common dimensionless focal definitions in JUMP and A549; diagnostic only."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

P=Path(__file__).resolve().parent;ROOT=P.parent;R=P/'results/harmonization';R.mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT/'phase2'));from analyze_discovery import calibrate
sys.path.insert(0,str(ROOT/'phase4'));from analyze import adjust,aggregate
NAMES=['adaptive_conditional','adaptive_unconditional','fixed_conditional','fixed_unconditional']

def measure(g,width,height,refarea,refmajor,selected=None):
    g=g.copy();v=g.area.between(.1*refarea,10*refarea)&np.isfinite(g[['x','y','theta','ar']]).all(axis=1)&(g.minor>0)&(g.x0>0)&(g.y0>0)&(g.x1<width)&(g.y1<height)
    alln=g[v].copy()
    if len(alln)<7:return pd.DataFrame()
    dd=cKDTree(alln[['x','y']]).query(alln[['x','y']],k=7)[0][:,6]
    alln['local_density']=5/(np.pi*dd**2)
    directors=alln[alln.ar>=1.2].copy()
    if len(directors)<100:return pd.DataFrame()
    directors['border']=np.minimum.reduce([directors.x,directors.y,width-1-directors.x,height-1-directors.y])
    if selected is None:f=directors.copy()
    else:f=selected.merge(directors,on='nucleus_id',how='inner',validate='many_to_one')
    if not len(f):return f
    xy=directors[['x','y']].to_numpy();z=np.exp(2j*directors.theta.to_numpy());ids=directors.nucleus_id.to_numpy();tree=cKDTree(xy);zz=np.exp(2j*f.theta.to_numpy());n=len(z)
    globalnull=(abs(z.sum())**2-n)/(n*(n-1));conditional=np.real(zz*np.conj((z.sum()-zz)/(n-1)))
    for mode in ['adaptive','fixed']:
        radii=2/np.sqrt(f.local_density.to_numpy()) if mode=='adaptive' else np.repeat(2/np.sqrt(len(alln)/(width*height)),len(f))
        neighborhoods=tree.query_ball_point(f[['x','y']].to_numpy(),radii)
        neighbors=[];obs=[]
        for i,js in enumerate(neighborhoods):
            q=np.asarray(js,dtype=int);q=q[ids[q]!=f.nucleus_id.iloc[i]];neighbors.append(len(q));obs.append(float(np.real(zz[i]*np.conj(z[q].mean()))) if len(q) else np.nan)
        valid=(np.asarray(neighbors)>=3)&(f.border.to_numpy()>=np.maximum(3*refmajor,radii+refmajor))
        observed=np.asarray(obs);f[mode+'_conditional']=np.where(valid,observed-conditional,np.nan);f[mode+'_unconditional']=np.where(valid,observed-globalnull,np.nan)
        f[mode+'_neighbors']=neighbors
    return f

def jump_geometry(row):
    name=f'{row.Metadata_Plate}__{row.Metadata_Well}__s{int(row.Metadata_Site)}_Nuclei.parquet'
    paths=[ROOT/'phase3/data/fields'/name,ROOT/'phase2/data/fields'/name]
    path=next((f for f in paths if f.exists()),None)
    if path is None:raise FileNotFoundError(name)
    d=pd.read_parquet(path);prefix='Nuclei_AreaShape_'
    ren={'ObjectNumber':'nucleus_id',prefix+'Area':'area',prefix+'MajorAxisLength':'major',prefix+'MinorAxisLength':'minor',prefix+'Center_X':'x',prefix+'Center_Y':'y',prefix+'BoundingBoxMinimum_X':'x0',prefix+'BoundingBoxMinimum_Y':'y0',prefix+'BoundingBoxMaximum_X':'x1',prefix+'BoundingBoxMaximum_Y':'y1'}
    g=d.rename(columns=ren);g['theta']=np.deg2rad(g[prefix+'Orientation']);g['ar']=g.major/g.minor.replace(0,np.nan)
    return g

def jump():
    out=R/'jump_common_fields.csv'
    old=pd.concat([pd.read_csv(f) for f in (ROOT/'phase2/results').glob('CP-*_fields.csv')],ignore_index=True);new=pd.read_csv(ROOT/'phase3/results/all_measured_fields.csv')
    d=pd.concat([old,new],ignore_index=True).drop_duplicates(['Metadata_Plate','Metadata_Well','Metadata_Site'],keep='last')
    d=d[d.Metadata_Batch.str.contains('Run4|Run5') & (d.role.isin(['control','position_control'])|(d.Metadata_Symbol=='DYNLL1'))].copy();d['analysis_split']='holdout'
    d['row']=d.Metadata_Well.str[0].map({chr(65+i):i for i in range(16)})/15;d['col']=(d.Metadata_Well.str[1:].astype(int)-1)/23;d['noguide']=d.Metadata_JCP2022.eq('JCP2022_800001').astype(float)
    if not out.exists():
        control=d[(d.role=='control')&(d.Metadata_JCP2022=='JCP2022_800002')];sc=[]
        for row in control.itertuples():
            g=jump_geometry(row);sc.append(g[['area','major']].replace([np.inf,-np.inf],np.nan).dropna())
        scales=pd.concat(sc).median();(R/'jump_control_scales.json').write_text(json.dumps(scales.to_dict(),indent=2));records=[]
        for k,row in enumerate(d.itertuples(),1):
            m=measure(jump_geometry(row),row.Width_OrigDNA,row.Height_OrigDNA,scales.area,scales.major)
            rec={'Metadata_Plate':row.Metadata_Plate,'Metadata_Well':row.Metadata_Well,'Metadata_Site':row.Metadata_Site}
            for name in NAMES:
                z=m[name].dropna() if name in m else pd.Series(dtype=float);rec[name]=float(z.mean()) if len(z)>=20 else np.nan;rec[name+'_n_focals']=len(z)
            records.append(rec)
            if k%100==0:print('JUMP common definitions',k,len(d),flush=True)
        pd.DataFrame(records).to_csv(out,index=False)
    merged=d.merge(pd.read_csv(out),on=['Metadata_Plate','Metadata_Well','Metadata_Site'],validate='one_to_one');tests=[]
    for name in NAMES:
        q=merged.copy();q['nuclear_excess_1p2']=q[name];q['nuclear_spatial_qc_1p2']=q[name].notna()
        for shape in [False,True]:
            a,_,model=calibrate(q,'nuclear',shape=shape);a.to_csv(R/f'jump_{name}_shape{int(shape)}_wells.csv',index=False)
            z=a[(a.role=='target')&(a.Metadata_Symbol=='DYNLL1')&a.reference_available]
            for batch,g in z.groupby('Metadata_Batch'):tests.append({'platform':'JUMP','definition':name,'shape_adjusted':shape,'batch':batch,'effect':float(g.adjusted.mean()),'raw_plate_delta':float(g.raw_plate_delta.mean()),'density_supported':bool(g.density_supported.all()),'valid_sites':float(g.valid_sites.mean())})
    pd.DataFrame(tests).to_csv(R/'jump_definition_effects.csv',index=False);print(pd.DataFrame(tests).to_string(index=False),flush=True)

def a549():
    cache=R/'a549_common_focal_cells.parquet'
    if not cache.exists():
        base=pd.concat([pd.read_parquet(ROOT/'phase4/results/all_focal_cells.parquet'),pd.read_parquet(ROOT/'phase4/holdout/results/all_focal_cells.parquet')],ignore_index=True)
        scales=base[base.gene=='nontargeting'][['area','major']].median();(R/'a549_control_scales.json').write_text(json.dumps(scales.to_dict(),indent=2));parts=[]
        for k,(field,subset) in enumerate(base.groupby('field'),1):
            g=pd.read_parquet(ROOT/'phase4/data/geometry'/f'{field}.parquet')
            keep=['field','plate','well','gene','guide','cell_id','nucleus_id','quality','registered']
            m=measure(g,5500,5500,scales.area,scales.major,subset[subset.registered][keep])
            if len(m):parts.append(m)
            if k%50==0:print('A549 common definitions',k,270,flush=True)
        pd.concat(parts,ignore_index=True).to_parquet(cache,index=False)
    d=pd.read_parquet(cache);records=[]
    for name in NAMES:
        a,_=adjust(d,name);gp,g=aggregate(a);g.to_csv(R/f'a549_{name}_guide_effects.csv',index=False);gp.to_csv(R/f'a549_{name}_guide_plate_effects.csv',index=False)
        t=g[g.gene=='DYNLL1'];c=g[(g.gene=='nontargeting')&g.eligible]
        records.append({'platform':'A549','definition':name,'effect':float(t.effect.mean()-c.effect.mean()),'guides':len(t),'negative_guides':int((t.effect<0).sum()),'target_cells':int(t.cells.sum())})
    pd.DataFrame(records).to_csv(R/'a549_definition_effects.csv',index=False);print(pd.DataFrame(records).to_string(index=False),flush=True)

if __name__=='__main__':
    policy=json.loads((P/'harmonization_policy.json').read_text());jump();a549()
