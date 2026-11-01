"""Frozen guide-level inference; do not call until every field is measured."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.linear_model import HuberRegressor
from sklearn.preprocessing import StandardScaler
from scipy.stats import spearmanr

P=Path(__file__).resolve().parent;RES=P/'results'
SEED=20261029

def adjust(d,endpoint,solidity=False,position=False):
    d=d[d[endpoint].notna()].copy()
    ld=np.log(d.local_density.to_numpy());la=np.log(d.ar.to_numpy());lz=np.log(d.area.to_numpy())
    x=np.column_stack([ld,ld**2,la,lz]+([d.solidity.to_numpy()] if solidity else []))
    if position:
        xx=d.x.to_numpy()/5500;yy=d.y.to_numpy()/5500
        x=np.column_stack([x,xx,yy,xx*xx,yy*yy,xx*yy])
    control=d.gene.eq('nontargeting').to_numpy()
    sc=StandardScaler().fit(x[control]);xs=sc.transform(x)
    model=HuberRegressor(max_iter=2000).fit(xs[control],d.loc[control,endpoint])
    d['raw_residual']=d[endpoint]-model.predict(xs)
    refs=d[control].groupby('field').raw_residual.mean()
    d['adjusted']=d.raw_residual-d.field.map(refs)
    # A field with no eligible NT after QC cannot supply a within-field baseline.
    d=d[d.adjusted.notna()].copy()
    return d,{'scaler_mean':sc.mean_.tolist(),'scaler_scale':sc.scale_.tolist(),'coef':model.coef_.tolist(),'intercept':model.intercept_,'n_controls':int(control.sum())}

def aggregate(d):
    gp=d.groupby(['gene','guide','plate'],as_index=False).agg(effect=('adjusted','mean'),cells=('adjusted','size'),fields=('field','nunique'))
    g=gp.groupby(['gene','guide'],as_index=False).agg(effect=('effect','mean'),cells=('cells','sum'),plates=('plate','nunique'))
    g['eligible']=(g.cells>=10)&(g.plates>=3)
    return gp,g

def welch_t(a,b):
    return (a.mean(axis=-1)-b.mean(axis=-1))/np.sqrt(a.var(axis=-1,ddof=1)/a.shape[-1]+b.var(axis=-1,ddof=1)/b.shape[-1])

def infer(g):
    rng=np.random.default_rng(SEED)
    t=g[(g.gene=='DYNLL1')&g.eligible].effect.to_numpy()
    c=g[(g.gene=='nontargeting')&g.eligible].effect.to_numpy()
    if len(t)<4 or len(c)<12:return {'testable':False,'target_guides':len(t),'nt_guides':len(c)}
    obs=float(welch_t(t,c));pool=np.r_[t,c];n=len(pool)
    total=pool.sum();squares=(pool**2).sum();null=[]
    # Draw uniform four-guide subsets without replacement using random keys.
    # NT guides remain clustered; no cell-level permutation is used.
    for k in range(20):
        ix=np.argpartition(rng.random((5000,n)),3,axis=1)[:,:4]
        a=pool[ix];am=a.mean(axis=1);av=a.var(axis=1,ddof=1)
        cm=(total-a.sum(axis=1))/(n-4)
        cv=(squares-(a*a).sum(axis=1)-(n-4)*cm**2)/(n-5)
        null.extend(((am-cm)/np.sqrt(av/4+cv/(n-4))).tolist())
    null=np.asarray(null)
    draws=rng.choice(t,(50000,4),replace=True).mean(axis=1)
    # Include uncertainty of the mean NT-guide reference.
    cb=[]
    for k in range(10):cb.extend(rng.choice(c,(5000,len(c)),replace=True).mean(axis=1).tolist())
    ci=np.quantile(draws-np.asarray(cb),[.025,.975])
    return {'testable':True,'target_guides':len(t),'nt_guides':len(c),'effect':float(t.mean()-c.mean()),'target_guide_mean':float(t.mean()),'nt_guide_mean':float(c.mean()),'welch_t':obs,'p_one_sided':float((1+(null>=obs).sum())/(len(null)+1)),'p_two_sided':float((1+(abs(null)>=abs(obs)).sum())/(len(null)+1)),'guide_bootstrap_95ci':ci.tolist(),'positive_target_guides':int((t>0).sum()),'leave_one_guide_out_effects':[(float(np.delete(t,i).mean()-c.mean())) for i in range(len(t))],'permutations':len(null)}

def main():
    d=pd.read_parquet(RES/'all_focal_cells.parquet')
    design=pd.read_csv(P/'frozen_fields_v2.csv')
    assert set(d.field)==set(design.field) and len(d)==len(pd.read_csv(P/'frozen_cells_v2.csv'))
    assert d[d.gene=='nontargeting'].groupby('plate').size().index.nunique()==7
    br=json.loads((P/'biological_replicate_amendment.json').read_text())['map']
    d['biorep']=d.plate.map(br);assert d.biorep.notna().all()
    qc=d.groupby(['gene','guide','plate'],as_index=False).agg(selected=('cell_id','size'),registered=('registered','sum'),eligible=('primary',lambda x:x.notna().sum()))
    qc.to_csv(RES/'guide_plate_qc.csv',index=False)
    results={};models={};sensitivity=[]
    for endpoint,sol in [('primary',False),('primary',True),('ar11',False),('ar13',False),('r15',False),('r25',False),('fixed100',False),('raw_alignment',False)]:
        name=endpoint+('_solidity_adjusted' if sol else '')
        a,model=adjust(d,endpoint,sol);models[name]=model
        gp,g=aggregate(a);test=infer(g);results[name]=test
        gp.to_csv(RES/f'{name}_guide_plate.csv',index=False);g.to_csv(RES/f'{name}_guides.csv',index=False)
        if endpoint=='primary' and not sol:
            a.to_parquet(RES/'primary_adjusted_cells.parquet',index=False)
            plates=gp.groupby('plate',as_index=False).agg(effect=('effect','mean'))
            # Gene effects only: control guide effects are not part of target plate means.
            plates=gp[gp.gene=='DYNLL1'].groupby('plate',as_index=False).agg(effect=('effect','mean'),cells=('cells','sum'),guides=('guide','nunique'))
            plates['eligible']=plates.cells>=3;plates.to_csv(RES/'primary_target_plates.csv',index=False)
            b=gp[gp.gene=='DYNLL1'].copy();b['biorep']=b.plate.map(br)
            bg=b.groupby(['biorep','guide'],as_index=False).agg(effect=('effect','mean'),cells=('cells','sum'))
            bt=bg.groupby('biorep',as_index=False).agg(effect=('effect','mean'),cells=('cells','sum'),guides=('guide','nunique'))
            bg.to_csv(RES/'primary_guide_bioreplicates.csv',index=False);bt.to_csv(RES/'primary_bioreplicates.csv',index=False)
            test['biological_replicates']=bt.to_dict('records')
            test['positive_eligible_plates']=int(((plates.effect>0)&plates.eligible).sum());test['eligible_plates']=int(plates.eligible.sum())
            test['bioreplicate_consistency_pass']=bool(len(bt)==3 and (bt.effect>0).all() and (bt.cells>=10).all())
            loo=[]
            for rep in sorted(a.biorep.unique()):
                _,lg=aggregate(a[a.biorep!=rep]);tg=lg[(lg.gene=='DYNLL1')&lg.eligible];cg=lg[(lg.gene=='nontargeting')&lg.eligible]
                loo.append({'excluded_biorep':rep,'target_guides':len(tg),'effect':float(tg.effect.mean()-cg.effect.mean()) if len(tg) else None})
            test['leave_one_bioreplicate_out']=loo
            allgp=gp.copy();allgp['biorep']=allgp.plate.map(br)
            brg=allgp.groupby(['gene','guide','biorep'],as_index=False).agg(effect=('effect','mean'),cells=('cells','sum'))
            eq=brg.groupby(['gene','guide'],as_index=False).agg(effect=('effect','mean'),cells=('cells','sum'),bioreplicates=('biorep','nunique'))
            eq['plates']=eq.guide.map(g.set_index('guide').plates);eq['eligible']=(eq.cells>=10)&(eq.plates>=3)&(eq.bioreplicates==3)
            eq.to_csv(RES/'equal_bioreplicate_guides.csv',index=False)
            test['equal_bioreplicate_weight_sensitivity']=infer(eq)
            mat=bg.pivot(index='biorep',columns='guide',values='effect')
            if mat.shape==(3,4) and mat.notna().all().all():
                rng2=np.random.default_rng(SEED);m=mat.to_numpy();v=[]
                for k in range(10):
                    ri=rng2.integers(0,3,(5000,3));ci=rng2.integers(0,4,(5000,4));v.extend(m[ri[:,:,None],ci[:,None,:]].mean(axis=(1,2)).tolist())
                test['two_way_bioreplicate_guide_bootstrap_95ci']=np.quantile(v,[.025,.975]).tolist()
                test['two_way_bootstrap_scope']='Descriptive equal-biological-replicate mean; only three transductions and four guides, including an overlapping pair.'
            # Descriptive field-paired bootstrap, distinct from the guide primary test.
            f=a[a.gene=='DYNLL1'].groupby('field').adjusted.mean().to_numpy()
            rng=np.random.default_rng(SEED);bs=[]
            for k in range(10):bs.extend(rng.choice(f,(5000,len(f)),replace=True).mean(axis=1).tolist())
            test['field_bootstrap_95ci']=np.quantile(bs,[.025,.975]).tolist();test['field_mean']=float(f.mean());test['eligible_target_fields']=len(f)
        sensitivity.append({'endpoint':name,'role':'primary' if name=='primary' else 'sensitivity; cannot rescue primary failure',**{k:v for k,v in test.items() if not isinstance(v,(list,dict))}})
    primary=results['primary'];sol=results['primary_solidity_adjusted']
    primary['full_support_rule_pass']=bool(primary.get('testable',False) and primary['p_one_sided']<=.05 and primary['positive_target_guides']>=3 and min(primary['leave_one_guide_out_effects'])>0 and primary['positive_eligible_plates']/max(1,primary['eligible_plates'])>=.7 and sol.get('effect',-np.inf)>0 and primary.get('bioreplicate_consistency_pass',False))
    # Rank-concordance is a registration check, not another biological validation.
    check=[]
    for f in sorted((P.parent/'phase2/data/periscope').glob('*_DYNLL1.csv.gz')):
        x=pd.read_csv(f,usecols=['Metadata_Foci_site','Metadata_Nuclei_ObjectNumber_x','Nuclei_AreaShape_Solidity','Nuclei_AreaShape_Eccentricity'])
        x=x.rename(columns={'Metadata_Foci_site':'field','Metadata_Nuclei_ObjectNumber_x':'nucleus_id'})
        x=x.merge(d[d.gene=='DYNLL1'][['field','nucleus_id','plate','solidity','ar']],on=['field','nucleus_id'],validate='one_to_one')
        for plate,z in x.groupby('plate'):
            z=z.dropna()
            if len(z)>=4:check.append({'guide':f.name.split('__')[-1].split('_')[0],'plate':plate,'cells':len(z),'solidity_spearman':float(spearmanr(z.solidity,z.Nuclei_AreaShape_Solidity).statistic),'ar_eccentricity_spearman':float(spearmanr(z.ar,z.Nuclei_AreaShape_Eccentricity).statistic)})
    pd.DataFrame(check).to_csv(RES/'physical_normalized_registration_check.csv',index=False)
    pd.DataFrame(sensitivity).to_csv(RES/'spatial_sensitivity.csv',index=False)
    (RES/'spatial_tests.json').write_text(json.dumps(results,indent=2,allow_nan=False))
    (RES/'control_models.json').write_text(json.dumps(models,indent=2))
    print(json.dumps(primary,indent=2))

if __name__=='__main__':main()
