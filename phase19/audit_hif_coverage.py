"""Trace every displayed HIF guide/replicate to the full released count table.

The original guide eligibility and reported effects are kept unchanged.
Replicate-level denominator flags and pseudocount/coverage diagnostics are
additional descriptive checks, not a replacement analysis or validation.
"""
from pathlib import Path
import hashlib,json,re
import numpy as np
import pandas as pd
from numpy.testing import assert_allclose
P=Path(__file__).resolve().parent;ROOT=P.parent;R=P/'results'
def main():
    R.mkdir(exist_ok=True);(P/'qa').mkdir(exist_ok=True)
    source=ROOT/'phase9/data/GSE290538_U2OS_HRE_Count_table.txt.gz'
    a=pd.read_csv(source,sep='\t',low_memory=False)
    old=pd.read_csv(ROOT/'phase17/panel_data/supplementary_figure5_F.csv')
    num=a[['High1','High2','High3']].to_numpy(float)
    den=a[['Low1','Low2','Low3']].to_numpy(float)
    depths=np.r_[num.sum(axis=0),den.sum(axis=0)];scale=float(np.median(depths))
    numerator=num/num.sum(axis=0)*scale;denominator=den/den.sum(axis=0)*scale
    eligible=(denominator.mean(axis=1)>=30)&((denominator>=10).sum(axis=1)>=2)
    ntc=a.Gene.astype(str).str.match(r'^(NonTargetingControlGuideForHuman_|Non-TargetingControl_|NONT|Control$)',case=False).to_numpy()
    uncentered=np.log2((numerator+1)/(denominator+1))
    centers=np.median(uncentered[ntc&eligible],axis=0);lfc=uncentered-centers
    idx=a.index[a.Gene.eq('DYNLL1')].to_numpy();assert len(idx)==4
    target=a.loc[idx].set_index('sgRNA').reindex(old.guide)
    target_idx=a.index[a.sgRNA.isin(old.guide)].to_numpy()
    # Align explicitly to the retained figure order rather than relying on rows.
    target_idx=a.set_index('sgRNA').index.get_indexer(old.guide)
    assert (target_idx>=0).all()
    assert_allclose(lfc[target_idx],old[['LFC_r1','LFC_r2','LFC_r3']],rtol=1e-12,atol=1e-12)
    assert np.array_equal(eligible[target_idx],old.eligible.to_numpy())
    assert_allclose(denominator[target_idx].mean(axis=1),old.denominator_mean,rtol=1e-12)
    rows=[];wide=old.copy();sensitivity=[]
    for i,k in enumerate(target_idx):
        for j in range(3):
            low=bool(denominator[k,j]<10)
            wide.loc[i,f'raw_denominator_r{j+1}']=den[k,j]
            wide.loc[i,f'raw_numerator_r{j+1}']=num[k,j]
            wide.loc[i,f'replicate_denominator_below_10_r{j+1}']=low
            rows.append({'guide':a.sgRNA.iloc[k],'guide_prefix':a.sgRNA.iloc[k].split('_',1)[1][:7],
                'replicate':j+1,'raw_low_count':int(den[k,j]),'raw_high_count':int(num[k,j]),
                'normalized_low_count':denominator[k,j],'normalized_high_count':numerator[k,j],
                'guide_mean_normalized_low':float(denominator[k].mean()),
                'guide_replicates_low_ge_10':int((denominator[k]>=10).sum()),
                'original_guide_eligible':bool(eligible[k]),
                'replicate_denominator_below_10':low,'zero_denominator':bool(den[k,j]==0),
                'both_fractions_zero':bool(den[k,j]==0 and num[k,j]==0),
                'pseudocount_normalized_units':1,'uncentered_log2_ratio':uncentered[k,j],
                'NTC_median_center':centers[j],'displayed_centered_log2_ratio':lfc[k,j]})
    long=pd.DataFrame(rows);long.to_csv(R/'HIF_DYNLL1_all_12_guide_replicate_coverage.csv',index=False)
    wide.to_csv(R/'HIF_DYNLL1_figure_source_with_coverage.csv',index=False)
    for pseudo in [.1,1.,10.]:
        u=np.log2((numerator+pseudo)/(denominator+pseudo));c=np.median(u[ntc&eligible],axis=0);v=u-c
        for i,k in enumerate(target_idx):
            for j in range(3):sensitivity.append({'guide':a.sgRNA.iloc[k],'replicate':j+1,'pseudocount':pseudo,'NTC_median_center':c[j],'centered_log2_ratio':v[k,j],'replicate_denominator_below_10':bool(denominator[k,j]<10)})
    pd.DataFrame(sensitivity).to_csv(R/'HIF_DYNLL1_pseudocount_sensitivity.csv',index=False)
    retained=eligible[target_idx];diagnostics=[]
    for j in range(3):
        covered=retained&(denominator[target_idx,j]>=10)
        diagnostics.append({'replicate':j+1,'original_guide_count':int(retained.sum()),
            'original_gene_median':float(np.median(lfc[target_idx[retained],j])),
            'covered_guide_replicate_count':int(covered.sum()),
            'covered_only_gene_median':float(np.median(lfc[target_idx[covered],j])) if covered.any() else None})
    pd.DataFrame(diagnostics).to_csv(R/'HIF_DYNLL1_replicate_coverage_diagnostic.csv',index=False)
    x=long[long.guide_prefix.eq('ACAAGAA')&long.replicate.eq(2)].iloc[0]
    assert x.raw_low_count==0 and x.raw_high_count==1480 and x.original_guide_eligible
    assert long.replicate_denominator_below_10.sum()==4
    report={'all_pass':True,'raw_source':source.relative_to(ROOT).as_posix(),
        'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'all_original_12_display_values_reproduced':True,'original_guide_filter_unchanged':True,
        'original_guide_rule':'Mean normalized low-fraction count >=30 and low-fraction count >=10 in at least 2 of 3 replicates; no numerator filter.',
        'replicate_warning_rule':'Normalized low-fraction count <10 in this particular replicate; display flag only.',
        'flagged_cells':long.loc[long.replicate_denominator_below_10,['guide_prefix','replicate','raw_low_count','raw_high_count','original_guide_eligible']].to_dict('records'),
        'ACAAGAA_Rep2':x.to_dict(),
        'original_gene_mean':float(np.mean([r['original_gene_median'] for r in diagnostics])),
        'post_hoc_covered_only_gene_mean':float(np.mean([r['covered_only_gene_median'] for r in diagnostics])),
        'post_hoc_diagnostic_scope':'Replicate-wise coverage mask within originally eligible guides, with unchanged original NTC centering; no new inference or primary estimate replacement.',
        'interpretation':'The +10.9827 value has a zero denominator, depends on the pseudocount and NTC centering, and is not a well-measured approximately 2000-fold biological enrichment.'}
    (R/'HIF_coverage_audit.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
