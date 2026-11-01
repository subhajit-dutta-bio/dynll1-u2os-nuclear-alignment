"""Check inferential units, complete coverage and every reported new summary."""
from pathlib import Path
import json,itertools,hashlib
import numpy as np,pandas as pd
from pandas.testing import assert_frame_equal
P=Path(__file__).resolve().parent;ROOT=P.parent;R=P/'results'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    checks=[]
    def check(name,truth):assert truth,name;checks.append(name)
    plan=json.loads((P/'frozen_plan.json').read_text());design=pd.read_csv(P/'frozen_fields.csv').fillna({'gene':'','panel':''});manifest=pd.read_csv(P/'frozen_contrasts.csv');wells=pd.read_csv(P/'frozen_wells.csv');fields=pd.read_csv(R/'field_descriptors.csv');effects=pd.read_csv(R/'alignment_contrasts.csv');summary=json.loads((R/'comparison_summary.json').read_text());genes=pd.read_csv(R/'gene_summary.csv')
    check('All frozen design hashes preserved',all(sha(P/n)==v for n,v in plan['hashes'].items()))
    check('All frozen dependency hashes preserved',all(sha(ROOT/n)==v for n,v in plan['dependencies'].items()))
    check('7875 exact fields,875 wells,335 targets,67 genes and5 dates',len(fields)==7875 and fields.field.is_unique and set(fields.field)==set(design.field) and len(wells)==875 and len(manifest)==335 and manifest.gene.nunique()==67 and manifest.batch.nunique()==5)
    key=['gene','batch'];expected=set(map(tuple,manifest[key].to_numpy()));models=effects.model.unique()
    check('All six endpoint specifications present',set(models)=={'adaptive_conditional','adaptive_unconditional','fixed_conditional','fixed_unconditional','adaptive_conditional_quadratic_area','direct_matching'})
    for model,q in effects.groupby('model'):
        check(model+' preserves all335 planned contrasts',len(q)==335 and not q.duplicated(key).any() and set(map(tuple,q[key].to_numpy()))==expected)
        check(model+' missingness equals evaluability',q.effect.notna().equals(q.evaluable))
        check(model+' target and source coverage respected',q.loc[q.evaluable,'target_sites'].ge(6).all() and q.loc[q.evaluable,'source_sites_min'].ge(6).all())
    primary=effects[effects.model.eq('adaptive_conditional')];recalc=primary.groupby(['gene','panel'],as_index=False).agg(mean=('effect','mean'),dates=('effect','count'),positive_dates=('effect',lambda v:int(v.gt(0).sum())),SD=('effect','std'));assert_frame_equal(recalc,genes[recalc.columns],check_exact=False,rtol=1e-12,atol=1e-12);checks.append('Every gene mean/date-count/sign/SD reproduces')
    dyn=genes[genes.gene.eq('DYNLL1')].iloc[0];wide=primary.pivot(index='gene',columns='batch',values='effect');datecontrasts=pd.read_csv(R/'DYNLL1_comparator_date_contrasts.csv')
    for group,rec in summary['groups'].items():
        q=genes[genes.panel.eq(group)&genes.complete_five_dates];check(group+' complete membership preserved',set(q.gene)==set(rec['complete_gene_names']))
        if len(q):
            check(group+' ranks and quantiles reproduce',np.isclose(q['mean'].median(),rec['median'],atol=1e-12) and np.isclose(q['mean'].quantile(.95),rec['q95'],atol=1e-12) and int(q['mean'].ge(dyn['mean']).sum()+1)==rec['DYNLL1_rank_within_set_plus_anchor'])
            for row in datecontrasts[datecontrasts.group.eq(group)].itertuples():check(group+' '+row.batch+' paired median',np.isclose(row.DYNLL1_minus_median,wide.loc['DYNLL1',row.batch]-wide.loc[q.gene,row.batch].median(),atol=1e-12))
    sd=pd.read_csv(R/'density_area_support.csv');matched=summary['groups']['matched_other_chromosome'];rule=bool(dyn.complete_five_dates and matched['complete_genes']>=24 and dyn['mean']>matched['q95'] and matched.get('date_contrast',{}).get('positive_dates',0)>=4 and sd[sd.gene.eq('DYNLL1')].local_density_fraction_inside.ge(.8).all());check('Specificity criterion recomputes without threshold changes',rule==summary['DYNLL1_specificity_descriptive_rule_met'])
    for record in [summary['DYNLL1_same_cohort'],summary['motor_module']]+[x.get('date_contrast',{}) for x in summary['groups'].values()]:
        if record.get('evaluable'):check('Five-date exact fractions respect32 assignments',record['sign_assignments']==32 and record['two_sided_tail']>=.0625 and np.isclose(record['positive_tail']*32,round(record['positive_tail']*32)))
    d=pd.read_csv(R/'descriptor_contrasts.csv');check('All335 rows in every descriptor/model retained',all(len(q)==335 and set(map(tuple,q[key].to_numpy()))==expected for _,q in d.groupby(['metric','model'])))
    context=pd.read_csv(R/'gene_context_and_covariates.csv').set_index('gene')
    for (metric,model),q in d.groupby(['metric','model']):
        column=metric+'_'+model if metric=='scaled_NND' else metric;means=q.groupby('gene').effect.mean();check(metric+'/'+model+' context means reproduce',np.allclose(means,context.loc[means.index,column],rtol=1e-12,atol=1e-12,equal_nan=True))
    ranks={metric:{'rank':int(context[metric].gt(context.loc['DYNLL1',metric]).sum()+1) if descending else int(context[metric].lt(context.loc['DYNLL1',metric]).sum()+1),'genes':int(context[metric].count()),'order':'descending' if descending else 'ascending'} for metric,descending in [('mean',True),('log_AR',False),('position_axis_coupling',False),('solidity',True)]}
    check('Reported DYNLL1 multi-descriptor extrema hold in all67 targets',all(v['rank']==1 and v['genes']==67 for v in ranks.values()))
    for metric,positive in [('log_AR',False),('position_axis_coupling',False),('solidity',True)]:
        q=d[d.gene.eq('DYNLL1')&d.metric.eq(metric)].effect;check(metric+' reported five-date direction verified',len(q)==5 and (q.gt(0).all() if positive else q.lt(0).all()))
    raw=pd.read_csv(R/'all_site_covariate_contrasts.csv');rawfields=pd.read_csv(R/'all_site_covariate_fields.csv');rawwells=pd.read_csv(R/'all_site_covariate_wells.csv');check('All-site covariate QC retains every field, well and planned contrast',len(rawfields)==7875 and len(rawwells)==875 and rawwells.sites.eq(9).all() and len(raw)==335 and set(map(tuple,raw[key].to_numpy()))==expected)
    retention=pd.read_csv(R/'matching_retention.csv');check('Matching retention retains every planned target/date, including absent focals',len(retention)==335 and set(map(tuple,retention[key].to_numpy()))==expected and retention.loc[retention.primary_focals.eq(0),'retained_fraction'].isna().all())
    balance=pd.read_csv(R/'matching_covariate_balance.csv');check('Post-outcome matching-balance audit retains all335 rows and original retention',len(balance)==335 and set(map(tuple,balance[key].to_numpy()))==expected and balance.set_index(key).primary_focals.equals(retention.set_index(key).primary_focals) and balance.set_index(key).matched_focals.equals(retention.set_index(key).matched_focals))
    for column,caliper in [('local_density',.2),('ar',.2),('area',.3)]:check(column+' matching-balance distances respect original caliper',balance[column+'_maximum_absolute_log_difference'].dropna().le(caliper+1e-12).all())
    check('All-site raw count ratios reproduce',np.allclose(raw.count_target_source_ratio,raw['count_target']/raw['count_source_NT'].replace(0,np.nan),equal_nan=True,rtol=1e-12,atol=1e-12))
    pred=pd.read_csv(R/'prediction/OOF_predictions.csv');stats=pd.read_csv(R/'prediction/summary.csv');check('All eight prediction configurations retained',len(stats)==8 and stats.groupby('model').size().eq(2).all())
    for row in stats.itertuples():
        q=pred[pred.model.eq(row.model)&pred.grouping.eq(row.grouping)];sse=((q.effect-q.predicted)**2).sum();sst=((q.effect-q.effect.mean())**2).sum();check(row.model+'/'+row.grouping+' R2 and complete predictions',len(q)==row.observations and q.predicted.notna().all() and np.isclose(1-sse/sst,row.R2,atol=1e-12) and np.isclose((q.effect-q.predicted).abs().mean(),row.MAE,atol=1e-12))
    old={'DYNLL1_U2OS_main_manuscript_v2.docx':'3079b3b45af46a5c0dd15e958e355c19e08af0d1f694f7f8c90afaf42e6c69ed','DYNLL1_U2OS_supplementary_information_v2.docx':'8ad2f76a8f2ede737861f0596a81c571e03fa5229d09d240e07089925e3aab04'};check('Both earlier version2 documents remain byte-identical',all(sha(ROOT/n)==h for n,h in old.items()))
    report={'all_pass':True,'checks':checks,'new_comparator_genes':67,'planned_gene_date_estimates_per_model':335,'complete_primary_genes':int(genes.complete_five_dates.sum()),'DYNLL1_descriptive_ranks':ranks,'nucleus_level_significance_used':False,'all_failed_or_missing_estimates_retained':True};(P/'qa/numeric_verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
