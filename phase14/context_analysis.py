"""Observed covariate balance and genomic context of the fixed panel."""
from pathlib import Path
import json,re
import pandas as pd,numpy as np
from scipy.stats import spearmanr
P=Path(__file__).resolve().parent;ROOT=P.parent;R=P/'results'
def main():
    g=pd.read_csv(R/'gene_summary.csv');assert len(g)==67;context=pd.read_csv(ROOT/'phase13/results/all_gene_coordinate_coverage.csv.gz');g=g.merge(context[['gene','chrom','position','location','distance_bp','profile_similarity_to_DYNLL1']],on='gene',validate='one_to_one')
    p=pd.read_parquet(ROOT/'phase13/data/matching_profiles.parquet').rename(columns={'Metadata_Symbol':'gene'});g=g.merge(p,on='gene',validate='one_to_one');d=pd.read_csv(ROOT/'phase9/data/DepMap24Q4_U2OS_CRISPRGeneEffect.csv');row=d.iloc[0].to_dict();fitness={re.sub(r' \(\d+\)$','',c):float(row[c]) for c in d.columns[1:]};g['U2OS_DepMap_gene_effect']=g.gene.map(fitness)
    metrics=pd.read_csv(R/'descriptor_contrasts.csv');means=metrics.groupby(['gene','metric','model'],as_index=False).agg(value=('effect','mean'),dates=('effect','count'))
    for (metric,model),q in means.groupby(['metric','model']):
        key=metric if metric!='scaled_NND' else metric+'_'+model;g=g.merge(q[['gene','value','dates']].rename(columns={'value':key,'dates':key+'_dates'}),on='gene',how='left',validate='one_to_one')
    g.to_csv(R/'gene_context_and_covariates.csv',index=False);rows=[]
    for variable in ['Cells_Count_Count','Nuclei_AreaShape_Area','U2OS_DepMap_gene_effect','log_density','log_area','log_AR','position_axis_coupling','scaled_NND_density_only','scaled_NND_plus_area','scaled_NND_plus_area_AR','scaled_NND_quadratic_area_AR']:
        use=g[g.complete_five_dates&np.isfinite(g[['mean',variable]]).all(axis=1)]
        rows.append({'variable':variable,'genes':len(use),'Spearman':float(spearmanr(use['mean'],use[variable]).statistic),'interpretation':'Selected-panel descriptive correlation; no causal or genome-wide inference'})
    pd.DataFrame(rows).to_csv(R/'alignment_covariate_correlations.csv',index=False)
    target=g[g.gene.eq('DYNLL1')].iloc[0];matched=g[g.panel.eq('matched_other_chromosome')];balance=[]
    for variable in ['Cells_Count_Count','Nuclei_AreaShape_Area','U2OS_DepMap_gene_effect','log_density','log_area','log_AR']:
        q=matched[variable].dropna();balance.append({'variable':variable,'DYNLL1':target[variable],'matched_genes':len(q),'matched_mean':q.mean(),'matched_median':q.median(),'matched_q25':q.quantile(.25),'matched_q75':q.quantile(.75),'matched_min':q.min(),'matched_max':q.max(),'DYNLL1_inside_matched_IQR':bool(q.quantile(.25)<=target[variable]<=q.quantile(.75))})
    pd.DataFrame(balance).to_csv(R/'matched_comparator_balance.csv',index=False)
    effects=pd.read_csv(R/'alignment_contrasts.csv');summary=effects.groupby(['model','gene','panel'],as_index=False).agg(mean=('effect','mean'),dates=('effect','count'),positive_dates=('effect',lambda x:int(x.gt(0).sum())));summary.to_csv(R/'all_definition_gene_summary.csv',index=False)
    leave=[];primary=effects[effects.model.eq('adaptive_conditional')];complete=g[g.complete_five_dates&g.panel.eq('matched_other_chromosome')].gene;wide=primary[primary.gene.isin(list(complete)+['DYNLL1'])].pivot(index='gene',columns='batch',values='effect')
    for omitted in sorted(primary.batch.unique()):
        q=wide.drop(columns=omitted).mean(axis=1);other=q.drop('DYNLL1');leave.append({'omitted_date':omitted,'DYNLL1_four_date_mean':q['DYNLL1'],'matched_median_four_date_mean':other.median(),'difference':q['DYNLL1']-other.median(),'rank_plus_anchor':int(other.ge(q['DYNLL1']).sum()+1),'matched_complete_genes':len(other)})
    pd.DataFrame(leave).to_csv(R/'leave_one_date_out_comparator_sensitivity.csv',index=False)
    print('Context and measured comparator balance complete; all67 genes retained')
if __name__=='__main__':main()
