"""All-site nuclear count and area context without collective coverage selection."""
from pathlib import Path
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='2'
import json,argparse,hashlib
import pandas as pd,numpy as np
P=Path(__file__).resolve().parent;ROOT=P.parent;R=P/'results';SCALES=json.loads((ROOT/'phase5/jump_raw_policy.json').read_text())['source_control_scales']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def raw_fields(design):
    rows=[]
    for i,row in enumerate(design.itertuples(),1):
        g=pd.read_parquet(R/'geometry'/f'{row.field}.parquet');s=SCALES[row.batch];valid=g.area.between(.1*s['area'],10*s['area'])&np.isfinite(g[['x','y','theta','ar']]).all(axis=1)&g.minor.gt(0)&g.x0.gt(0)&g.y0.gt(0)&g.x1.lt(1280)&g.y1.lt(1080);g=g[valid];rows.append({'field':row.field,'batch':row.batch,'plate':row.plate,'well':row.well,'site':row.site,'role':row.role,'gene':row.gene,'panel':row.panel,'nuclei':len(g),'directors':int(g.ar.ge(1.2).sum()),'nuclei_per_mm2':len(g)/(1280*1080*.65**2)*1e6,'median_area_um2':g.area.median()*.65**2,'median_AR':g.ar.median()})
        if i%1000==0:print('All-site count/area QC',i,len(design),flush=True)
    return pd.DataFrame(rows)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--offline-inputs',action='store_true');a=ap.parse_args();design=pd.read_csv(P/'frozen_fields.csv').fillna({'gene':'','panel':''});manifest=pd.read_csv(P/'frozen_contrasts.csv');assert len(design)==7875 and json.loads((P/'qa/analysis_provenance.json').read_text())['all_fields_completed'];source=R/'analysis_raw_covariate_inputs.parquet'
    if a.offline_inputs:
        provenance=json.loads((P/'qa/coverage_covariates.json').read_text());assert sha(source)==provenance['measured_input_sha256'];d=pd.read_parquet(source)
    else:d=raw_fields(design);d.to_parquet(source,index=False)
    expected=pd.read_csv(R/'field_descriptors.csv').set_index('field');observed=d.set_index('field');assert len(d)==7875 and d.field.is_unique and np.array_equal(observed.loc[expected.index,'nuclei'],expected.valid_nuclei) and np.array_equal(observed.loc[expected.index,'directors'],expected.directors)
    d.to_csv(R/'all_site_covariate_fields.csv',index=False);w=d.groupby(['batch','plate','well','role','gene','panel'],as_index=False,dropna=False).agg(sites=('field','size'),mean_nuclei=('nuclei','mean'),total_nuclei=('nuclei','sum'),mean_nuclei_per_mm2=('nuclei_per_mm2','mean'),mean_field_median_area_um2=('median_area_um2','mean'),area_defined_sites=('median_area_um2','count'),mean_field_median_AR=('median_AR','mean'),collective_director_sites=('directors',lambda v:int(v.ge(100).sum())));assert len(w)==875 and w.sites.eq(9).all();w.to_csv(R/'all_site_covariate_wells.csv',index=False);records=[]
    for row in manifest.itertuples():
        b=w[w.batch.eq(row.batch)];t=b[b.plate.eq(row.target_plate)&b.well.eq(row.target_well)];sourceNT=b[b.plate.eq(row.target_plate)&b.well.isin(row.source_NT_wells.split(';'))];p=b[b.plate.eq(row.negative_plate)&b.well.eq(row.negative_well)];anchor=b[b.plate.eq(row.negative_plate)&b.well.isin(row.anchor_wells.split(';'))];assert len(t)==1 and len(sourceNT)==2 and len(p)==1 and len(anchor)==2;rec={'gene':row.gene,'panel':row.panel,'batch':row.batch}
        for label,column in [('count','mean_nuclei'),('area_um2','mean_field_median_area_um2'),('AR','mean_field_median_AR')]:
            vals=[float(x[column].mean()) for x in [t,sourceNT,p,anchor]]
            for name,v in zip(['target','source_NT','position','anchor'],vals):rec[label+'_'+name]=v
            rec[label+'_target_source_ratio']=vals[0]/vals[1] if vals[1]>0 else np.nan;rec[label+'_raw_position_contrast']=vals[0]-vals[1]-vals[2]+vals[3]
        for name,frame in [('target',t),('source_NT',sourceNT),('position',p),('anchor',anchor)]:rec['area_defined_sites_'+name]=int(frame.area_defined_sites.min())
        rec['area_all_required_sites_defined']=all(frame.area_defined_sites.eq(9).all() for frame in [t,sourceNT,p,anchor])
        if not rec['area_all_required_sites_defined']:rec['area_um2_target_source_ratio']=np.nan;rec['area_um2_raw_position_contrast']=np.nan;rec['AR_target_source_ratio']=np.nan;rec['AR_raw_position_contrast']=np.nan
        rec['target_director_eligible_sites']=int(t.collective_director_sites.iloc[0]);records.append(rec)
    context=pd.DataFrame(records);assert len(context)==335;context.to_csv(R/'all_site_covariate_contrasts.csv',index=False);(P/'qa/coverage_covariates.json').write_text(json.dumps({'all_pass':True,'fields':len(d),'wells':len(w),'target_date_rows':len(context),'all_nine_sites_retained':True,'zero_count_fields':int(d.nuclei.eq(0).sum()),'measured_counts_equal_original_geometry':True,'measured_input_sha256':sha(R/'analysis_raw_covariate_inputs.parquet'),'addendum_sha256':sha(P/'coverage_QC_addendum.json'),'scope':'Raw count/area coverage context, not a relaxation of spatial eligibility or a proliferation assay'},indent=2));print('All-site covariate context complete',len(context),'planned target/date rows')
if __name__=='__main__':main()
