"""Transparent exploratory U2OS screen effects and paralog-aware AP-MS audit."""
from pathlib import Path
import hashlib,itertools,json,zipfile
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

P=Path(__file__).resolve().parent;ROOT=P.parent;D=ROOT/'phase9/data';R=P/'results'
R.mkdir(exist_ok=True);(R/'guides').mkdir(exist_ok=True)
POLICY=json.loads((P/'analysis_policy.json').read_text())
CONTEXT=['DYNLL1','DYNLL2']+POLICY['fixed_context_genes']

def control_mask(genes):
    return genes.astype(str).str.match(r'^(NonTargetingControlGuideForHuman_|Non-TargetingControl_|NONT|Control$)',case=False).to_numpy()

def load_count(acc,file,sep,gene,idcol):
    a=pd.read_csv(D/file,sep=sep,low_memory=False).rename(columns={gene:'gene',idcol:'guide'})
    assert a.guide.is_unique and a.gene.notna().all()
    a['gene_original']=a['gene'];a['gene']=a['gene'].replace({'H2AFX':'H2AX','MRE11A':'MRE11'})
    return a

def eligibility(den,threshold=30):
    return (den.mean(axis=1)>=threshold)&((den>=10).sum(axis=1)>=int(np.ceil(den.shape[1]/2)))

def signflip(v):
    v=np.asarray(v);null=np.array([np.mean(v*s) for s in itertools.product([-1,1],repeat=len(v))])
    return float(np.mean(np.abs(null)>=abs(v.mean())-1e-12))

def ntc_null(lfc,den,eligible,ntc,target,seed):
    pool=np.flatnonzero(ntc&eligible);targets=np.flatnonzero(target&eligible)
    if len(pool)<30 or not len(targets):return {'ntc_pseudogene_p_two_sided':None,'null_draws':0}
    rng=np.random.default_rng(seed);base=den.mean(axis=1);edges=np.unique(np.quantile(base[pool],[0,.2,.4,.6,.8,1]));bins=np.digitize(base,edges[1:-1]);draws=np.zeros((20000,len(targets)),dtype=int);used=np.zeros_like(draws,dtype=bool);fallback=0
    for b in np.unique(bins[targets]):
        cols=np.flatnonzero(bins[targets]==b);candidate=pool[bins[pool]==b]
        if len(candidate)<20:candidate=pool;fallback+=len(cols)
        x=rng.choice(candidate,size=(len(draws),len(cols)),replace=True)
        if len(cols)>1:
            duplicate=(np.diff(np.sort(x,axis=1),axis=1)==0).any(axis=1)
            while duplicate.any():
                x[duplicate]=rng.choice(candidate,size=(int(duplicate.sum()),len(cols)),replace=True)
                duplicate=(np.diff(np.sort(x,axis=1),axis=1)==0).any(axis=1)
        draws[:,cols]=x
    # Different bins are disjoint except a sparse-bin fallback. Ensure all
    # pseudo-gene guides are distinct even when such a fallback occurs.
    dup=(np.diff(np.sort(draws,axis=1),axis=1)==0).any(axis=1)
    while dup.any():
        rows=np.flatnonzero(dup)
        for row in rows:draws[row]=rng.choice(pool,size=len(targets),replace=False)
        dup=(np.diff(np.sort(draws,axis=1),axis=1)==0).any(axis=1)
    null=np.median(lfc[draws],axis=1).mean(axis=1);obs=np.median(lfc[targets],axis=0).mean();lo=(np.sum(null<=obs)+1)/(len(null)+1);hi=(np.sum(null>=obs)+1)/(len(null)+1)
    return {'ntc_pseudogene_p_two_sided':float(min(1,2*min(lo,hi))),'null_draws':len(null),'null_2_5':float(np.quantile(null,.025)),'null_97_5':float(np.quantile(null,.975)),'NTC_pool':len(pool),'sparse_bin_fallback_guides':fallback}

def run_contrast(acc,name,a,numerator,denominator,normalized=False):
    assert len(numerator)==len(denominator)
    num=a[numerator].to_numpy(float);den=a[denominator].to_numpy(float);assert (num>=0).all() and (den>=0).all()
    ntc=control_mask(a.gene);guideids=a.guide.astype(str).to_numpy();genes=a.gene.astype(str).to_numpy();n=len(numerator)
    depths=np.r_[num.sum(axis=0),den.sum(axis=0)];scale=1 if normalized else np.median(depths)
    if not normalized:num=num/num.sum(axis=0)*scale;den=den/den.sum(axis=0)*scale
    raw=np.log2((num+1)/(den+1));el=eligibility(den)
    center=np.median(raw[ntc&el],axis=0) if ntc.any() else np.zeros(n)
    lfc=raw-center
    rows=a[['guide','gene']].copy();rows['denominator_mean']=den.mean(axis=1);rows['eligible']=el;rows['NTC']=ntc
    for j in range(n):rows[f'denominator_r{j+1}']=den[:,j];rows[f'numerator_r{j+1}']=num[:,j];rows[f'LFC_r{j+1}']=lfc[:,j];rows[f'uncentered_LFC_r{j+1}']=raw[:,j]
    rows.to_parquet(R/'guides'/f'{name}.parquet',index=False)
    cols=[f'LFC_r{j+1}' for j in range(n)]
    gs=rows[el&~ntc].groupby('gene')[cols].median();gs['guide_count']=rows[el&~ntc].groupby('gene').size();gs['mean_LFC']=gs[cols].mean(axis=1);gs['depletion_percentile']=gs.mean_LFC.rank(pct=True)*100
    gs.to_csv(R/f'{name}_gene_effects.csv.gz',compression='gzip')
    target_results=[];sens=[]
    for gene in CONTEXT:
        target=genes==gene;valid=target&el;idx=np.flatnonzero(valid)
        if not len(idx):
            target_results.append({'contrast':name,'accession':acc,'gene':gene,'total_guides':int(target.sum()),'eligible_guides':0,'n_replicates':n,'mean_LFC':None});continue
        v=np.median(lfc[idx],axis=0);effect=float(v.mean());loo=[float(np.median(lfc[np.setdiff1d(idx,[i])],axis=0).mean()) for i in idx] if len(idx)>1 else []
        record={'contrast':name,'accession':acc,'gene':gene,'primary_contrast':name in POLICY['primary_contrasts'],'total_guides':int(target.sum()),'eligible_guides':len(idx),'n_replicates':n,'mean_LFC':effect,'replicate_min':float(v.min()),'replicate_max':float(v.max()),'positive_replicates':int((v>0).sum()),'negative_replicates':int((v<0).sum()),'all_replicates_agree':bool((v>0).all() or (v<0).all()),'exact_replicate_signflip_p_two_sided':signflip(v),'depletion_percentile':float(gs.loc[gene,'depletion_percentile']),'positive_guides':int((lfc[idx].mean(axis=1)>0).sum()),'negative_guides':int((lfc[idx].mean(axis=1)<0).sum()),'LOO_min':min(loo) if loo else None,'LOO_max':max(loo) if loo else None,'LOO_same_direction':bool(all(np.sign(x)==np.sign(effect) for x in loo)) if loo else None,'normalization':'published_normalized_no_NTC' if normalized else 'library_depth_then_NTC_median_centering'}
        for j,value in enumerate(v):record[f'LFC_r{j+1}']=float(value)
        if gene in ['DYNLL1','DYNLL2']:
            seed=POLICY['seed']+int(hashlib.sha256((name+gene).encode()).hexdigest()[:8],16)
            record.update(ntc_null(lfc,den,el,ntc,target,seed))
        target_results.append(record)
        for threshold in [0,20,30,50,100]:
            e=eligibility(den,threshold)&target
            sens.append({'contrast':name,'gene':gene,'sensitivity':f'denominator_threshold_{threshold}','guides':int(e.sum()),'mean_LFC':float(np.median(lfc[e],axis=0).mean()) if e.any() else None})
        for label,arr in [('uncentered',raw),('all_guide_median_centered',raw-np.median(raw[el],axis=0))]:
            sens.append({'contrast':name,'gene':gene,'sensitivity':label,'guides':len(idx),'mean_LFC':float(np.median(arr[idx],axis=0).mean())})
    qc={'contrast':name,'accession':acc,'numerator_columns':'|'.join(numerator),'denominator_columns':'|'.join(denominator),'replicates':n,'guides':len(a),'eligible_guides':int(el.sum()),'NTCs':int(ntc.sum()),'eligible_NTCs':int((ntc&el).sum()),'normalization_centers':center.tolist(),'input_already_normalized':normalized}
    for i in range(n):
        for j in range(i+1,n):
            qc[f'guide_LFC_spearman_r{i+1}_r{j+1}']=float(spearmanr(lfc[el,i],lfc[el,j]).statistic)
            qc[f'gene_LFC_spearman_r{i+1}_r{j+1}']=float(spearmanr(gs[f'LFC_r{i+1}'],gs[f'LFC_r{j+1}']).statistic)
    print(name,'DYNLL1',next(x for x in target_results if x['gene']=='DYNLL1').get('mean_LFC'),flush=True)
    return target_results,sens,qc

def network():
    cx=json.loads((D/'U2OS_APMS.cx').read_text());aspects={k:v for obj in cx for k,v in obj.items()};nodes={i['@id']:i['n'] for i in aspects['nodes']};edges={(nodes[e['s']],nodes[e['t']]) for e in aspects['edges']};neighbor={g:({b if a==g else a for a,b in edges if a==g or b==g}-{g}) for g in ['DYNLL1','DYNLL2']}
    with zipfile.ZipFile(P/'data/CM4AI_PPIDownloader.zip') as z:
        attrs=pd.read_csv(z.open(next(n for n in z.namelist() if n.endswith('ppi_gene_node_attributes.tsv'))),sep='\t')
        rawedges=pd.read_csv(z.open(next(n for n in z.namelist() if n.endswith('edgeList.tsv'))),sep='\t')
        clean=pd.read_csv(z.open(next(n for n in z.namelist() if n.endswith('ppi_edgelist.tsv'))),sep='\t')
        baits=pd.read_csv(z.open(next(n for n in z.namelist() if n.endswith('baitList.tsv'))),sep='\t')
    cleanpairs=set(map(tuple,clean.to_numpy()));assert {frozenset(x) for x in cleanpairs}=={frozenset(x) for x in edges}
    attrs.to_csv(R/'APMS_node_attributes.csv',index=False);baits[baits.GeneSymbol.isin(['DYNLL1','DYNLL2'])].to_csv(R/'LC8_bait_metadata.csv',index=False)
    rows=[]
    for g in sorted(neighbor['DYNLL1']|neighbor['DYNLL2']):
        rows.append({'neighbor':g,'DYNLL1':g in neighbor['DYNLL1'],'DYNLL2':g in neighbor['DYNLL2'],'class':'shared' if g in neighbor['DYNLL1']&neighbor['DYNLL2'] else ('DYNLL1_only' if g in neighbor['DYNLL1'] else 'DYNLL2_only'),'DYNLL1_as_bait':('DYNLL1',g) in cleanpairs,'DYNLL2_as_bait':('DYNLL2',g) in cleanpairs})
    pd.DataFrame(rows).to_csv(R/'APMS_paralog_neighbors.csv',index=False)
    assemblies=pd.read_excel(D/'41586_2025_8878_MOESM1_ESM.xlsx',sheet_name='Table1',header=1);out=[]
    for x in assemblies.to_dict('records'):
        proteins=set(str(x['Proteins']).split());hit=proteins&{'DYNLL1','DYNLL2'}
        if not hit:continue
        out.append({'assembly':x['Assembly ID'],'annotation':x['Annotation'],'protein_count':len(proteins),'DYNLL1':'DYNLL1' in proteins,'DYNLL2':'DYNLL2' in proteins,'DYNLL1_neighbor_overlap':len(proteins&neighbor['DYNLL1']),'DYNLL2_neighbor_overlap':len(proteins&neighbor['DYNLL2']),'context_genes':'|'.join(sorted(proteins&set(CONTEXT))),'all_proteins':'|'.join(sorted(proteins))})
    pd.DataFrame(out).sort_values('protein_count').to_csv(R/'LC8_all_published_assemblies.csv',index=False)
    result={'DYNLL1_neighbors':len(neighbor['DYNLL1']),'DYNLL2_neighbors':len(neighbor['DYNLL2']),'shared_neighbors':len(neighbor['DYNLL1']&neighbor['DYNLL2']),'Jaccard':len(neighbor['DYNLL1']&neighbor['DYNLL2'])/len(neighbor['DYNLL1']|neighbor['DYNLL2']),'DYNLL1_only':sorted(neighbor['DYNLL1']-neighbor['DYNLL2']),'DYNLL2_only':sorted(neighbor['DYNLL2']-neighbor['DYNLL1']),'both_are_affinity_baits':bool(attrs.set_index('name').loc[['DYNLL1','DYNLL2'],'bait'].all()),'clean_network_exact_undirected_agreement':True,'raw_edge_rows':len(rawedges),'clean_edge_rows':len(clean),'interpretation':'Separate tagged baits provide paralog-resolved association context. Shared localization images remain antibody-ambiguous. A missing edge does not establish specificity or absence; no quantitative differential-binding contrast or causal validation.'}
    (R/'APMS_deep_summary.json').write_text(json.dumps(result,indent=2));return result

def main():
    allresults=[];sensitivity=[];qcs=[]
    a=load_count('GSE232592','GSE232592_CRISPR_p53_counts.txt.gz','\t','Gene','ID')
    specs=[]
    for group,label in [('WT','WT'),('p53KO','p53KO')]:
        specs += [('GSE232592',f'GSE232592_{label}_growth',a,[f'{group}_Untreated_{r}' for r in 'ABC'],[f'{group}_T0_{r}' for r in 'ABC'],False),('GSE232592',f'GSE232592_{label}_etoposide',a,[f'{group}_Etoposide_{r}' for r in 'ABC'],[f'{group}_Untreated_{r}' for r in 'ABC'],False)]
    b=load_count('GSE290538','GSE290538_U2OS_HRE_Count_table.txt.gz','\t','Gene','sgRNA')
    for label,num,den in [('HIF_high_vs_low','High','Low'),('high_vs_library','High','Brunello'),('low_vs_library','Low','Brunello')]:specs.append(('GSE290538','GSE290538_'+label,b,[num+str(i) for i in [1,2,3]],[den+str(i) for i in [1,2,3]],False))
    c=load_count('GSE283147','GSE283147_combined-guide-level-counts.csv.gz',',','GENE','ID');cols={key:[next(x for x in c.columns if x.startswith(key+'_'+r+'_')) for r in 'AB'] for key in ['DEAD','LIVE','T0','UNT']}
    for label,num,den in [('dead_vs_live','DEAD','LIVE'),('growth','UNT','T0'),('live_vs_untreated','LIVE','UNT')]:specs.append(('GSE283147','GSE283147_'+label,c,cols[num],cols[den],False))
    e=load_count('GSE233195','GSE233195_JSMO.count.txt.gz','\t','Gene','sgRNA');specs.append(('GSE233195','GSE233195_screen_vs_input',e,['S1_B01','S2_B02'],['In1_A01','In2_A02'],False))
    for context,prefix in [('parental','U2OS'),('hTERT_mutant','KO')]:
        parts=[pd.read_excel(D/'GSE226994_MAGeCK_count_normalized.xlsx',sheet_name=f'{prefix}-rep{i}.count_normalized') for i in [1,2]]
        merged=parts[0].merge(parts[1],on=['sgRNA','Gene'],validate='one_to_one').rename(columns={'sgRNA':'guide','Gene':'gene'});assert len(merged)==len(parts[0])==len(parts[1]);merged['gene_original']=merged.gene;merged['gene']=merged.gene.replace({'H2AFX':'H2AX','MRE11A':'MRE11'})
        specs.append(('GSE226994',f'GSE226994_{context}_growth',merged,[f'{prefix}-d21-rep{i}' for i in [1,2]],[f'{prefix}-d0-rep{i}' for i in [1,2]],True))
    for args in specs:
        result,sens,qc=run_contrast(*args);allresults+=result;sensitivity+=sens;qcs.append(qc)
    result=pd.DataFrame(allresults)
    for gene in ['DYNLL1','DYNLL2']:
        mask=result.gene.eq(gene)&result.primary_contrast.eq(True)&result.ntc_pseudogene_p_two_sided.notna()
        result.loc[mask,'primary_family_NTC_BH_q']=multipletests(result.loc[mask,'ntc_pseudogene_p_two_sided'],method='fdr_bh')[1]
    result.to_csv(R/'all_context_gene_effects.csv',index=False);result[result.gene.isin(['DYNLL1','DYNLL2'])].to_csv(R/'DYNLL1_DYNLL2_all_contrasts.csv',index=False);pd.DataFrame(sensitivity).to_csv(R/'normalization_and_count_sensitivity.csv',index=False);pd.DataFrame(qcs).to_json(R/'screen_QC.json',orient='records',indent=2)
    # The hTERT contrast is a difference of published normalized log ratios.
    wt=pd.read_parquet(R/'guides/GSE226994_parental_growth.parquet');ko=pd.read_parquet(R/'guides/GSE226994_hTERT_mutant_growth.parquet');joint=wt.merge(ko,on=['guide','gene'],suffixes=('_WT','_KO'),validate='one_to_one');out=[]
    for gene in CONTEXT:
        q=joint[joint.gene.eq(gene)&joint.eligible_WT&joint.eligible_KO]
        vals=[float(np.median(q[f'LFC_r{i}_KO']-q[f'LFC_r{i}_WT'])) for i in [1,2]] if len(q) else []
        out.append({'gene':gene,'guides':len(q),'replicate1_delta':vals[0] if vals else None,'replicate2_delta':vals[1] if vals else None,'mean_delta':float(np.mean(vals)) if vals else None,'scope':'Engineered hTERT context minus parental growth, descriptive; not a DYNLL1 spatial endpoint or verified independent clone test.'})
    pd.DataFrame(out).to_csv(R/'hTERT_context_interaction.csv',index=False)
    overlaps=pd.read_csv(ROOT/'phase9/results/public_CRISPR_guide_overlap_with_JUMP.csv');ind=[]
    for name in ['GSE290538_HIF_high_vs_low','GSE290538_high_vs_library','GSE290538_low_vs_library','GSE283147_dead_vs_live','GSE283147_growth','GSE283147_live_vs_untreated']:
        guides=pd.read_parquet(R/'guides'/f'{name}.parquet');acc=name.split('_')[0];m=overlaps[overlaps.source.eq(acc)].set_index('id');q=guides[guides.gene.eq('DYNLL1')].copy();q['overlap_bp']=q.guide.map(m.maximum_JUMP_target_overlap_bp);q['sequence']=q.guide.map(m.sequence);q.to_csv(R/f'{name}_DYNLL1_guide_audit.csv',index=False);non=q[q.eligible&q.overlap_bp.eq(0)];cols=[x for x in q if x.startswith('LFC_r')];ind.append({'contrast':name,'DYNLL1_guides':len(q),'eligible_nonoverlapping_guides':len(non),'nonoverlapping_median_LFC':float(non[cols].median().mean()) if len(non) else None,'sufficient_for_two_guide_support':len(non)>=2})
    pd.DataFrame(ind).to_csv(R/'nonoverlapping_guide_sensitivity.csv',index=False)
    net=network();summary={'contrasts':len(qcs),'primary_or_context_contrasts':len(POLICY['primary_contrasts']),'all_guides_retained':True,'new_alignment_validation':False,'new_biological_replicate_significance':False,'network':net,'DYNLL1':result[result.gene.eq('DYNLL1')].replace({np.nan:None}).to_dict('records'),'scope':'Exploratory functional readouts; no causal mechanism for neighbor alignment established.'};(R/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False));print('COMPLETE',len(qcs),'contrasts',flush=True)

if __name__=='__main__':main()
