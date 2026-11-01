"""Discovery ranking with control calibration, empirical nulls and a frozen shortlist.

The resulting q values quantify evidence relative to the sampled screen controls.
They do not resolve gene/position confounding or replace independent guide tests.
"""
from pathlib import Path
import argparse
import json
import hashlib
import warnings
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.linear_model import HuberRegressor
from sklearn.model_selection import GroupKFold,cross_val_predict
from sklearn.metrics import r2_score
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=Path(__file__).resolve().parent;R=P/'results'
ENDPOINTS={'cell':'cell_excess_1p2','nuclear':'nuclear_excess_1p2'}

def load_fields(split):
    frames=[pd.read_csv(path) for path in R.glob('CP-*_fields.csv')]
    data=pd.concat(frames,ignore_index=True)
    data=data[data.metrics_version.eq(3)&data.analysis_split.eq(split)].copy()
    if data.duplicated(['Metadata_Plate','Metadata_Well','Metadata_Site']).any():
        raise RuntimeError('Duplicate plate/well/site observations')
    data['row']=data.Metadata_Well.str[0].map({chr(65+i):i for i in range(16)})/15
    data['col']=(data.Metadata_Well.str[1:].astype(int)-1)/23
    data['noguide']=data.Metadata_JCP2022.eq('JCP2022_800001').astype(float)
    return data

def features(data,shape=False,position=True):
    # Calibrate batch with contemporaneous controls. Only controls fit this model.
    x=pd.DataFrame(index=data.index)
    logdensity=np.log(np.maximum(data.density_including_edges,1))
    loglocal=np.log(np.maximum(data.cell_local_density_median,1))
    x['log_density']=logdensity;x['log_density_squared']=logdensity**2
    x['log_local_density']=loglocal;x['noguide']=data.noguide
    if position:
        x['row']=data.row;x['column']=data.col;x['row_column']=data.row*data.col
    x=pd.concat([x,pd.get_dummies(data.Metadata_Batch,prefix='batch',dtype=float)],axis=1)
    if shape:
        for col in ['cell_ar_median','nuclear_ar_median','cell_q_median','nuclear_q_median','cell_area_fraction_retained']:
            x[col]=data[col]
    return x

def calibrate(data,endpoint,shape=False,position=True):
    name=ENDPOINTS[endpoint]
    eligible=data[endpoint+'_spatial_qc_1p2'].astype(bool)&np.isfinite(data[name])
    d=data.loc[eligible].copy()
    d['valid_sites']=d.groupby(['Metadata_Plate','Metadata_Well']).Metadata_Site.transform('nunique')
    d=d[d.analysis_split.ne('holdout')|d.valid_sites.ge(2)].copy()
    # Equal-area fields are averaged within a well before statistical inference.
    # This is a no-op for the initial single-field discovery pass.
    if d.duplicated(['Metadata_Plate','Metadata_Well']).any():
        metadata=['Metadata_Source','Metadata_Batch','Metadata_Plate','Metadata_Well','Metadata_JCP2022','Metadata_Symbol','panel','role','analysis_split']
        numbers=[c for c in d.select_dtypes(include=[np.number,bool]).columns if c not in metadata]
        agg={c:'mean' for c in numbers};agg.update({c:'first' for c in d.columns if c not in numbers and c not in ['Metadata_Plate','Metadata_Well']})
        d=d.groupby(['Metadata_Plate','Metadata_Well'],as_index=False).agg(agg)
    controls=d.role.isin(['control','position_control'])
    x=features(d,shape,position)
    model=make_pipeline(SimpleImputer(strategy='median'),StandardScaler(),HuberRegressor(epsilon=1.5,alpha=1,max_iter=2000,tol=1e-6))
    model.fit(x[controls],d.loc[controls,name])
    prediction=model.predict(x)
    residual=d[name].to_numpy()-prediction
    d['model_prediction']=prediction;d['residual']=residual
    # Use assay-matched non-targeting controls as the primary plate reference.
    nt=d[d.Metadata_JCP2022.eq('JCP2022_800002')]
    source_nt=nt[nt.role.eq('control')]
    plate_reference=source_nt.groupby('Metadata_Plate').residual.median()
    d['before_position_reference']=d.residual-d.Metadata_Plate.map(plate_reference)
    reference=d[d.role.eq('position_control')].copy()
    anchors=reference[reference.Metadata_Well.isin(['A02','I23'])].groupby('Metadata_Batch').residual.mean()
    reference['position_offset']=reference.residual-reference.Metadata_Batch.map(anchors)
    offsets=reference.groupby(['Metadata_Batch','Metadata_Well']).position_offset.mean()
    d['position_reference_offset']=[offsets.get((b,w),np.nan) for b,w in zip(d.Metadata_Batch,d.Metadata_Well)]
    d['adjusted']=d.before_position_reference-d.position_reference_offset
    raw_reference=source_nt.groupby('Metadata_Plate')[name].median()
    d['raw_plate_delta']=d[name]-d.Metadata_Plate.map(raw_reference)
    d['density_control_min']=d.Metadata_Batch.map(nt.groupby('Metadata_Batch').density_including_edges.min())
    d['density_control_max']=d.Metadata_Batch.map(nt.groupby('Metadata_Batch').density_including_edges.max())
    d['density_supported']=(d.density_including_edges>=d.density_control_min)&(d.density_including_edges<=d.density_control_max)
    # Separate negative-control classes are not assumed biologically identical.
    reference_count=source_nt.groupby('Metadata_Plate').Metadata_Well.nunique()
    d['reference_control_wells']=d.Metadata_Plate.map(reference_count)
    d['reference_available']=np.isfinite(d.adjusted)&d.reference_control_wells.ge(2)
    ctrl=d[d.role.eq('control')].copy()
    pairdiff=[]
    for (batch,plate,jcp),g in ctrl.groupby(['Metadata_Batch','Metadata_Plate','Metadata_JCP2022']):
        if len(g)>=2:
            # Repeated fields are averaged within well first; present run uses s0.
            values=g.groupby('Metadata_Well').adjusted.mean().to_numpy()
            if len(values)>=2:pairdiff.append((batch,plate,jcp,float(values[0]-values[1])))
    pairs=pd.DataFrame(pairdiff,columns=['batch','plate','jcp','difference'])
    all_controls=d[controls].copy()
    cv=r2_score(all_controls[name],cross_val_predict(model,x.loc[all_controls.index],all_controls[name],groups=all_controls.Metadata_Plate,cv=GroupKFold(5)))
    stats=dict(endpoint=endpoint,shape_adjusted=shape,position_adjusted=position,control_fields=len(ctrl),
               control_held_plate_R2=float(cv),control_pair_difference_sd=float(pairs.difference.std(ddof=1)),
               matched_position_control_fields=len(reference),
               density_supported_target_fields=int((d.role.eq('target')&d.density_supported).sum()))
    return d,pairs,stats

def null_distribution(pairs,batches,n=40000):
    rng=np.random.default_rng(20261029)
    result=np.zeros(n)
    for batch in batches:
        pool=pairs[pairs.batch.eq(batch)].difference.dropna().to_numpy()
        if not len(pool):return np.full(n,np.nan)
        # In the matched-position difference-in-differences contrast, a target
        # and its reference, each minus two anchors, have variance 3*sigma^2.
        # A matched-control pair difference has variance 4*sigma^2.
        result += rng.choice(pool,n,replace=True)*rng.choice([-1.,1.],n)*np.sqrt(.75)
    return result/len(batches)

def rank_genes(calibrated,pairs,endpoint):
    treated=calibrated[calibrated.role.eq('target')&calibrated.reference_available].copy()
    records=[]
    for gene,g in treated.groupby('Metadata_Symbol'):
        batches=g.groupby('Metadata_Batch').adjusted.mean()
        supported=g[g.density_supported].groupby('Metadata_Batch').adjusted.mean()
        effect=float(batches.mean())
        null=null_distribution(pairs,batches.index.tolist())
        p=(1+np.sum(np.abs(null)>=abs(effect)))/(len(null)+1) if np.isfinite(null).any() else np.nan
        records.append(dict(gene=gene,jcp2022=g.Metadata_JCP2022.iloc[0],panel=g.panel.iloc[0],endpoint=endpoint,
                            n_wells=g.Metadata_Well.groupby(g.Metadata_Plate).size().sum(),n_batches=len(batches),
                            effect=effect,median_effect=float(batches.median()),
                            raw_effect=float(g.groupby('Metadata_Batch').raw_plate_delta.mean().mean()),
                            density_supported_effect=float(supported.mean()) if len(supported) else np.nan,
                            density_supported_batches=len(supported),direction_consistency=float(np.mean(np.sign(batches)==np.sign(effect))),
                            between_batch_sd=float(batches.std(ddof=1)),empirical_p=float(p),
                            control_null_sd=float(np.nanstd(null,ddof=1)),
                            cell_density_mean=float(g.density_including_edges.mean()),
                            q_cell=float(g.cell_q_median.mean()),ar_cell=float(g.cell_ar_median.mean())))
    return pd.DataFrame(records)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--allow-partial',action='store_true');args=ap.parse_args()
    warnings.filterwarnings('ignore',category=RuntimeWarning)
    d=load_fields('discovery')
    expected=pd.concat([pd.read_csv(P/file) for file in ['design.csv','extension_design.csv','position_control_design.csv']],ignore_index=True)
    expected=expected.query("analysis_split == 'discovery'").drop_duplicates(['Metadata_Plate','Metadata_Well'])
    have=set(zip(d.Metadata_Plate,d.Metadata_Well));need=set(zip(expected.Metadata_Plate,expected.Metadata_Well))
    if not need.issubset(have) and not args.allow_partial:raise RuntimeError(f'Incomplete discovery extraction: {len(need&have)}/{len(need)} wells')
    if args.allow_partial:raise RuntimeError('Partial candidate ranking is deliberately disabled; use QC summaries while extracting')
    ranks=[];calibration=[];calibrated=[]
    for endpoint in ENDPOINTS:
        basic,pairs,stats=calibrate(d,endpoint,shape=False,position=True)
        ranks.append(rank_genes(basic,pairs,endpoint));calibration.append(stats)
        shape,_,s=calibrate(d,endpoint,shape=True,position=True);calibration.append(s)
        no_position,_,s=calibrate(d,endpoint,shape=False,position=False);calibration.append(s)
        effects=shape[shape.role.eq('target')].groupby('Metadata_Symbol').adjusted.mean()
        no_pos_effect=no_position[no_position.role.eq('target')].groupby('Metadata_Symbol').adjusted.mean()
        ranks[-1]['shape_adjusted_effect']=ranks[-1].gene.map(effects)
        ranks[-1]['no_position_effect']=ranks[-1].gene.map(no_pos_effect)
        ranks[-1]['shape_adjusted_same_direction']=np.sign(ranks[-1].effect)==np.sign(ranks[-1].shape_adjusted_effect)
        basic['endpoint']=endpoint;calibrated.append(basic)
    ranking=pd.concat(ranks,ignore_index=True)
    ranking['cohort']=np.where(ranking.panel.eq('standard_profile_extension'),'exploratory_extension','frozen_primary_panel')
    for cohort,indices in ranking.groupby('cohort').groups.items():
        valid=ranking.loc[indices].empirical_p.dropna().index
        ranking.loc[valid,'screen_fdr']=multipletests(ranking.loc[valid,'empirical_p'],method='fdr_bh')[1]
    ranking['consistent_candidate']=(ranking.n_batches>=3)&ranking.direction_consistency.eq(1)&ranking.density_supported_batches.ge(2)&ranking.shape_adjusted_same_direction
    ranking=ranking.sort_values(['screen_fdr','empirical_p','effect'],ascending=[True,True,False])
    ranking.to_csv(R/'discovery_gene_ranking.csv',index=False)
    pd.concat(calibrated).to_csv(R/'discovery_calibrated_fields.csv',index=False)
    (R/'control_calibration.json').write_text(json.dumps(calibration,indent=2))
    shortlist=ranking[ranking.consistent_candidate].copy()
    shortlist['strength']=abs(shortlist.effect)/shortlist.control_null_sd
    # A fixed maximum of ten genes limits downstream cherry-picking. Even if no
    # FDR threshold is passed, top consistent associations are tested as tentative.
    shortlist=shortlist.sort_values(['screen_fdr','strength'],ascending=[True,False]).drop_duplicates('gene').head(10)
    shortlist['status']=np.where(shortlist.screen_fdr.le(.1),'discovery_screen_support','tentative_below_screen_threshold')
    path=R/'frozen_candidates.csv';shortlist.to_csv(path,index=False)
    freeze=dict(date='2026-10-29',discovery_fields=len(d),candidate_genes=shortlist.gene.tolist(),
                candidate_file_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                primary_endpoints=list(ENDPOINTS),screen_tests=len(ranking),
                screen_fdr_supported_gene_endpoint_tests=int(ranking.screen_fdr.le(.1).sum()),
                rule='Up to ten genes with three-batch sign consistency, density support in at least two batches, and shape-adjusted sign agreement. Screen FDR ranks first; standardized magnitude breaks ties.',
                caveat='Gene and well position are partially confounded. These are screen associations, not independent experimental validation.')
    freeze['analysis_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    freeze['input_design_sha256']={file:hashlib.sha256((P/file).read_bytes()).hexdigest() for file in ['design.csv','extension_design.csv','position_control_design.csv']}
    (R/'candidate_freeze.json').write_text(json.dumps(freeze,indent=2))
    # Extract only frozen candidates and contemporaneous controls in holdout batches.
    design=pd.concat([pd.read_csv(P/'design.csv'),pd.read_csv(P/'extension_design.csv')],ignore_index=True).drop_duplicates(['Metadata_Plate','Metadata_Well'])
    targets=design[design.analysis_split.eq('holdout')&design.Metadata_Symbol.isin(shortlist.gene)]
    controls=design[design.analysis_split.eq('holdout')&design.role.eq('control')&design.Metadata_Plate.isin(targets.Metadata_Plate)]
    position=pd.read_csv(P/'position_control_design.csv')
    positions=set(targets.Metadata_Well)|{'A02','I23','B02','J23'}|set(controls.Metadata_Well)
    position=position[position.analysis_split.eq('holdout')&position.Metadata_Well.isin(positions)&position.Metadata_Batch.isin(targets.Metadata_Batch)]
    pd.concat([targets,controls,position]).to_csv(P/'holdout_design.csv',index=False)
    print(json.dumps(freeze,indent=2));print(ranking.head(24).to_string(index=False))
    # Replication across discovery batches is assessed across genes, not cells.
    checks=[]
    for endpoint,g in pd.concat(calibrated).query("role == 'target'").groupby('endpoint'):
        pivot=g.pivot_table(index='Metadata_Symbol',columns='Metadata_Batch',values='adjusted',aggfunc='mean')
        columns=pivot.columns.tolist()
        for i in range(len(columns)):
            for j in range(i+1,len(columns)):
                f=pivot[[columns[i],columns[j]]].dropna();rho,p=spearmanr(f.iloc[:,0],f.iloc[:,1])
                checks.append(dict(endpoint=endpoint,batch1=columns[i],batch2=columns[j],genes=len(f),spearman=float(rho),p=float(p)))
    pd.DataFrame(checks).to_csv(R/'discovery_replicability.csv',index=False)
    fig,axs=plt.subplots(1,2,figsize=(11,5),layout='constrained')
    for ax,endpoint in zip(axs,ENDPOINTS):
        f=ranking[ranking.endpoint.eq(endpoint)].sort_values('effect')
        top=f.loc[f.gene.isin(shortlist.gene)].copy()
        ax.scatter(f.effect,-np.log10(f.empirical_p),c=np.where(f.panel.eq('random_background'),'#B3B3B3','#24798A'),alpha=.8)
        for _,row in top.iterrows():ax.annotate(row.gene,(row.effect,-np.log10(row.empirical_p)),fontsize=8,xytext=(3,3),textcoords='offset points')
        ax.axvline(0,c='gray',lw=.7);ax.set(xlabel='Density and position calibrated alignment excess',ylabel='−log10 empirical screen p',title=endpoint.capitalize()+' spatial order')
        ax.spines[['top','right']].set_visible(False)
    fig.savefig(R/'discovery_spatial_hits.png',dpi=180);plt.close(fig)

if __name__=='__main__':main()
