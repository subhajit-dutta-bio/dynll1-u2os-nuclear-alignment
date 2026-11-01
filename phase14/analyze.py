"""Prespecified calibration, matching and gene-comparator reporting."""
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import sys,json,itertools,importlib.util,hashlib,argparse
import numpy as np,pandas as pd
from scipy.spatial import cKDTree
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import HuberRegressor
P=Path(__file__).resolve().parent;ROOT=P.parent;R=P/'results'
spec=importlib.util.spec_from_file_location('spatial_p14',ROOT/'phase11/spatial.py');sp=importlib.util.module_from_spec(spec);spec.loader.exec_module(sp)
ENDPOINTS=['adaptive_conditional','adaptive_unconditional','fixed_conditional','fixed_unconditional']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ready():
    plan=json.loads((P/'frozen_plan.json').read_text())
    for f,h in plan['hashes'].items():assert sha(P/f)==h,f
    design=pd.read_csv(P/'frozen_fields.csv').fillna({'gene':'','panel':''})
    assert all((R/'fields'/f'{f}.json').exists() for f in design.field),'Extraction incomplete; no partial gene ranking'
    return design,pd.read_csv(P/'frozen_contrasts.csv')
def aggregate(a,value,manifest):
    f=a[a[value].notna()].groupby(['batch','plate','well','role','field'],as_index=False,observed=True).agg(value=(value,'mean'),focals=(value,'size'));f=f[f.focals.ge(20)]
    w=f.groupby(['batch','plate','well','role'],as_index=False,observed=True).agg(value=('value','mean'),sites=('field','nunique'));return sp.contrast(w,manifest),f,w
def fit_alignment(d,endpoint,manifest,quadratic=False):
    a=d[d[endpoint].notna()].copy();counts=a.groupby('field',observed=True).size();a=a[a.field.isin(counts[counts.ge(20)].index)];ld=np.log(a.local_density);la=np.log(a.area)
    x=pd.DataFrame({'ld':ld,'ld2':ld**2,'lar':np.log(a.ar),'la':la,'noguide':a.jcp2022.eq('JCP2022_800001').astype(float)},index=a.index)
    if quadratic:x['la2']=la**2
    x=pd.concat([x,pd.get_dummies(a.batch,dtype=float)],axis=1);ctrl=a.role.ne('target');sc=StandardScaler().fit(x[ctrl]);xx=sc.transform(x);m=HuberRegressor(max_iter=3000).fit(xx[ctrl],a.loc[ctrl,endpoint]);assert m.n_iter_<3000;a['adjusted']=a[endpoint]-m.predict(xx);c,f,w=aggregate(a,'adjusted',manifest)
    info={'endpoint':endpoint,'quadratic_area':quadratic,'control_focals':int(ctrl.sum()),'focals':len(a),'features':x.columns.tolist(),'coef':m.coef_.tolist(),'intercept':float(m.intercept_),'means':sc.mean_.tolist(),'scales':sc.scale_.tolist(),'iterations':int(m.n_iter_)};return c,f,w,info,a
def support(a,manifest):
    out=[]
    for row in manifest.itertuples():
        b=a[a.batch.eq(row.batch)];t=b[b.plate.eq(row.target_plate)&b.well.eq(row.target_well)];ctrl=b[b.plate.eq(row.target_plate)&b.well.isin(row.source_NT_wells.split(';'))];rec={'gene':row.gene,'panel':row.panel,'batch':row.batch,'target_focals':len(t),'source_NT_focals':len(ctrl)}
        for col in ['local_density','area','ar']:
            low,high=ctrl[col].quantile([.01,.99]) if len(ctrl) else [np.nan,np.nan];rec[col+'_lower']=low;rec[col+'_upper']=high;rec[col+'_fraction_inside']=float(t[col].between(low,high).mean()) if len(t) and len(ctrl) else np.nan
        out.append(rec)
    return pd.DataFrame(out)
def direct_match(a,manifest):
    a=a.copy();a['matched']=np.nan;ctrlmask=a.role.ne('target')
    for (batch,noguide),g in a.groupby(['batch',a.jcp2022.eq('JCP2022_800001')],observed=True):
        ctrl=g[ctrlmask.loc[g.index]];cx=np.column_stack([np.log(ctrl.local_density)/.2,np.log(ctrl.ar)/.2,np.log(ctrl.area)/.3]);tree=cKDTree(cx);cy=ctrl.adaptive_conditional.to_numpy();cf=ctrl.field.to_numpy();values=np.full(len(g),np.nan);gx=np.column_stack([np.log(g.local_density)/.2,np.log(g.ar)/.2,np.log(g.area)/.3]);fields=g.field.to_numpy();ys=g.adaptive_conditional.to_numpy()
        for start in range(0,len(g),4096):
            stop=min(start+4096,len(g));ids=tree.query(gx[start:stop],k=64)[1]
            for j,js in enumerate(ids,start):
                js=js[(cf[js]!=fields[j])&(np.abs(cx[js]-gx[j]).max(axis=1)<=1)][:3]
                if len(js)==3:values[j]=ys[j]-cy[js].mean()
        a.loc[g.index,'matched']=values;print('Matched',batch,bool(noguide),len(g),flush=True)
    retention=a[a.role.eq('target')].groupby(['gene','batch'],observed=True).matched.agg(primary_focals='size',matched_focals='count').reset_index();retention=manifest[['gene','panel','batch']].merge(retention,on=['gene','batch'],how='left',validate='one_to_one');retention[['primary_focals','matched_focals']]=retention[['primary_focals','matched_focals']].fillna(0).astype(int);retention['retained_fraction']=retention.matched_focals/retention.primary_focals.replace(0,np.nan);c,f,w=aggregate(a,'matched',manifest);return c,f,w,retention
def fit_spacing(fields,manifest):
    results=[];models=[];plans=[('density_only',[]),('plus_area',['log_area']),('plus_area_AR',['log_area','log_AR']),('quadratic_area_AR',['log_area','log_area_squared','log_AR'])]
    a=fields[fields.directors.ge(100)&np.isfinite(fields[['scaled_NND','log_density','log_area','log_AR']]).all(axis=1)].copy()
    for name,extra in plans:
        x=pd.DataFrame({'noguide':a.jcp2022.eq('JCP2022_800001').astype(float),'ld':a.log_density,'ld2':a.log_density**2},index=a.index)
        for key in extra:x[key]=a.log_area**2 if key=='log_area_squared' else a[key]
        x=pd.concat([x,pd.get_dummies(a.batch,dtype=float)],axis=1);ctrl=a.role.ne('target');sc=StandardScaler().fit(x[ctrl]);xx=sc.transform(x);m=HuberRegressor(max_iter=3000).fit(xx[ctrl],a.loc[ctrl,'scaled_NND']);assert m.n_iter_<3000;a['value']=a.scaled_NND-m.predict(xx);w=a.groupby(['batch','plate','well','role'],as_index=False).agg(value=('value','mean'),sites=('field','nunique'));c=sp.contrast(w,manifest);c['metric']='scaled_NND';c['model']=name;c['control_residual_SD']=a.loc[ctrl,'value'].std(ddof=1);c['effect_control_SD']=c.effect/c.control_residual_SD;results.append(c);models.append({'model':name,'features':x.columns.tolist(),'coef':m.coef_.tolist(),'intercept':float(m.intercept_),'means':sc.mean_.tolist(),'scales':sc.scale_.tolist(),'iterations':int(m.n_iter_),'fields':len(a),'control_fields':int(ctrl.sum())})
    return pd.concat(results,ignore_index=True),models
def signflip(v):
    v=np.asarray(v,float)
    if len(v)!=5 or not np.isfinite(v).all():return {'evaluable':False}
    null=np.asarray(list(itertools.product([-1,1],repeat=5)))@v/5;obs=v.mean();return {'evaluable':True,'mean':float(obs),'positive_dates':int((v>0).sum()),'positive_tail':float((null>=obs-1e-12).mean()),'two_sided_tail':float((np.abs(null)>=abs(obs)-1e-12).mean()),'sign_assignments':32}
def summarize(primary,density):
    summary=primary.groupby(['gene','panel'],as_index=False).agg(mean=('effect','mean'),dates=('effect','count'),positive_dates=('effect',lambda v:int(v.gt(0).sum())),SD=('effect','std'));summary['complete_five_dates']=summary.dates.eq(5);summary['descriptive_date_SE']=summary.SD/np.sqrt(summary.dates);summary.to_csv(R/'gene_summary.csv',index=False)
    d=summary[summary.gene.eq('DYNLL1')].iloc[0];groups={};datecomparisons=[];wide=primary.pivot(index='gene',columns='batch',values='effect');dyn=wide.loc['DYNLL1'];dates=wide.columns.tolist()
    for group in ['matched_other_chromosome','proximal_12q','distal_12q','motor_context']:
        q=summary[summary.panel.eq(group)];complete=q[q.complete_five_dates];means=complete['mean'];g={'planned':len(q),'complete_genes':len(complete),'complete_gene_names':complete.gene.tolist(),'median':float(means.median()) if len(means) else None,'q95':float(means.quantile(.95)) if len(means) else None,'DYNLL1_rank_within_set_plus_anchor':int(means.ge(d['mean']).sum()+1) if d.complete_five_dates else None,'DYNLL1_mean':float(d['mean']) if d.complete_five_dates else None}
        if len(complete):
            median=wide.loc[complete.gene].median(axis=0);delta=dyn-median
            for batch in dates:datecomparisons.append({'group':group,'batch':batch,'complete_comparator_genes':len(complete),'DYNLL1_effect':dyn[batch],'group_median':median[batch],'DYNLL1_minus_median':delta[batch]})
            g['date_contrast']=signflip(delta.to_numpy())
        groups[group]=g
    pd.DataFrame(datecomparisons).to_csv(R/'DYNLL1_comparator_date_contrasts.csv',index=False)
    motor=primary[primary.panel.eq('motor_context')];motor_dates=[]
    for batch,g in motor.groupby('batch'):motor_dates.append({'batch':batch,'genes_evaluable':int(g.effect.notna().sum()),'value':float(g.effect.mean()) if g.effect.notna().sum()==3 else np.nan})
    md=pd.DataFrame(motor_dates);md.to_csv(R/'motor_module_dates.csv',index=False);mtest=signflip(md.value);mtest['all_three_gene_means_positive']=bool(summary[summary.panel.eq('motor_context')]['mean'].gt(0).all());mtest['scope']='Exploratory motor-context association; not a complete NuMA module or causal mechanism.'
    q=groups['matched_other_chromosome'];dd=density[density.gene.eq('DYNLL1')];densityok=len(dd)==5 and dd.local_density_fraction_inside.ge(.8).all();positive=q.get('date_contrast',{}).get('positive_dates',0);rule=bool(d.complete_five_dates and q['complete_genes']>=24 and d['mean']>q['q95'] and positive>=4 and densityok)
    result={'scope':'Targeted comparison, not a genome-wide null; ranks are descriptive','DYNLL1_same_cohort':signflip(dyn.to_numpy()),'groups':groups,'motor_module':mtest,'DYNLL1_specificity_descriptive_rule_met':rule,'DYNLL1_all_dates_density_support_80percent':bool(densityok),'minimum_DYNLL1_density_support':float(dd.local_density_fraction_inside.min()),'complete_genes':int(summary.complete_five_dates.sum()),'planned_genes':67};(R/'comparison_summary.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
    # Repeatability is measured across genes and dates, not inferred from cell N.
    complete=wide.dropna();long=[]
    for a,b in itertools.combinations(dates,2):long.append({'date_a':a,'date_b':b,'genes':len(complete),'Spearman':float(stats.spearmanr(complete[a],complete[b]).statistic)})
    pd.DataFrame(long).to_csv(R/'date_repeatability.csv',index=False)
def main(offline=False):
    if offline:
        design=pd.read_csv(P/'frozen_fields.csv').fillna({'gene':'','panel':''});manifest=pd.read_csv(P/'frozen_contrasts.csv');input_record=json.loads((R/'offline_input_provenance.json').read_text());assert len(design)==7875 and len(manifest)==335
        for name,digest in input_record['inputs'].items():assert sha(R/name)==digest,name
        for name,digest in json.loads((P/'frozen_plan.json').read_text())['hashes'].items():assert sha(P/name)==digest,name
    else:design,manifest=ready()
    files=[R/'focals'/f'{f}.parquet' for f in design.field];parts=[]
    keep=['field','batch','plate','well','role','gene','jcp2022','area','ar','local_density']+ENDPOINTS
    if offline:d=pd.read_parquet(R/'analysis_focal_inputs.parquet')
    else:
        for f in files:
            q=pd.read_parquet(f)
            if len(q):parts.append(q[keep])
        d=pd.concat(parts,ignore_index=True)
    rawfields=d.groupby('field',as_index=False).agg(raw_alignment=('adaptive_conditional','mean'),primary_focals=('adaptive_conditional','count'),median_local_density=('local_density','median'),median_focal_area=('area','median'))
    rawfields.to_csv(R/'raw_alignment_field_summary.csv',index=False)
    for key in ['field','batch','plate','well','role','gene','jcp2022']:d[key]=d[key].astype('category')
    d.to_parquet(R/'analysis_focal_inputs.parquet',index=False);models=[];contrasts=[];fieldparts=[];primary=None
    for endpoint in ENDPOINTS:
        c,f,w,m,a=fit_alignment(d,endpoint,manifest);c['model']=endpoint;contrasts.append(c);models.append(m);f['model']=endpoint;fieldparts.append(f);w.to_csv(R/(endpoint+'_wells.csv'),index=False)
        if endpoint=='adaptive_conditional':primary=c.copy();a_primary=a.copy();density=support(a,manifest);density.to_csv(R/'density_area_support.csv',index=False)
        print('Calibrated',endpoint,'evaluable',int(c.evaluable.sum()),len(c),flush=True)
    c,f,w,m,_=fit_alignment(d,'adaptive_conditional',manifest,quadratic=True);c['model']='adaptive_conditional_quadratic_area';contrasts.append(c);models.append(m);f['model']=c.model.iloc[0];fieldparts.append(f)
    c,f,w,retention=direct_match(a_primary,manifest);c['model']='direct_matching';contrasts.append(c);f['model']='direct_matching';fieldparts.append(f);w.to_csv(R/'direct_matching_wells.csv',index=False);retention.to_csv(R/'matching_retention.csv',index=False)
    pd.concat(contrasts,ignore_index=True).to_csv(R/'alignment_contrasts.csv',index=False);pd.concat(fieldparts,ignore_index=True).to_csv(R/'adjusted_alignment_fields.csv',index=False);(R/'alignment_models.json').write_text(json.dumps(models,indent=2))
    fields=pd.read_parquet(R/'analysis_field_inputs.parquet') if offline else pd.DataFrame([json.loads((R/'fields'/f'{f}.json').read_text()) for f in design.field]);fields['panel']=fields.field.map(design.set_index('field').panel);fields.to_csv(R/'field_descriptors.csv',index=False);fields.to_parquet(R/'analysis_field_inputs.parquet',index=False);metrics=[];modelrecords=[]
    for metric in sp.METRICS:
        if metric=='scaled_NND':continue
        _,c,m=sp.fit_metric(fields,metric,manifest);c['model']='standard';metrics.append(c);modelrecords.append(m)
    spacing,spacingmodels=fit_spacing(fields,manifest);metrics.append(spacing);pd.concat(metrics,ignore_index=True).to_csv(R/'descriptor_contrasts.csv',index=False);(R/'descriptor_models.json').write_text(json.dumps(modelrecords+spacingmodels,indent=2))
    summarize(primary,density);qc=fields.groupby(['batch','plate','well','role','gene','panel'],as_index=False,dropna=False).agg(fields=('field','size'),valid_nuclei=('valid_nuclei','sum'),directors=('directors','sum'),primary_eligible_fields=('primary_focals',lambda v:int(v.ge(20).sum())));qc.to_csv(R/'well_QC.csv',index=False)
    (P/'qa/analysis_provenance.json').write_text(json.dumps({'all_fields_completed':True,'fields':len(design),'focal_input_rows':len(d),'policy_sha256':sha(P/'analysis_policy.json'),'script_sha256':sha(Path(__file__)),'focal_input_sha256':sha(R/'analysis_focal_inputs.parquet'),'partial_ranking_disallowed':True},indent=2))
    (R/'offline_input_provenance.json').write_text(json.dumps({'scope':'Complete measured focal and field inputs for offline statistical reproduction; does not rerun pixel inference.','fields':len(design),'inputs':{name:sha(R/name) for name in ['analysis_focal_inputs.parquet','analysis_field_inputs.parquet']}},indent=2))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--offline-inputs',action='store_true');args=ap.parse_args();main(args.offline_inputs)
