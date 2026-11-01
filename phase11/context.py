"""Exploratory hTERT, network, sequence and raw-count diagnostic additions."""
from pathlib import Path
import hashlib,json,itertools,zipfile
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
P=Path(__file__).resolve().parent;ROOT=P.parent;R=P/'results';R.mkdir(exist_ok=True);OLD=ROOT/'phase10/results'
SEED=20261031
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def main():
    wt=pd.read_parquet(OLD/'guides/GSE226994_parental_growth.parquet');ko=pd.read_parquet(OLD/'guides/GSE226994_hTERT_mutant_growth.parquet');joint=wt.merge(ko,on=['guide','gene'],suffixes=('_WT','_KO'),validate='one_to_one');inter=pd.read_csv(OLD/'hTERT_all_gene_interactions.csv.gz').set_index('gene');parent=pd.read_csv(OLD/'GSE226994_parental_growth_gene_effects.csv.gz').set_index('gene');mutant=pd.read_csv(OLD/'GSE226994_hTERT_mutant_growth_gene_effects.csv.gz').set_index('gene');paired=parent[['mean_LFC']].rename(columns={'mean_LFC':'parental_LFC'}).join(mutant[['mean_LFC']].rename(columns={'mean_LFC':'engineered_LFC'}),how='inner').join(inter[['delta_r1','delta_r2','mean_delta','guides']],how='inner');paired.to_csv(R/'hTERT_genome_fitness_and_interactions.csv')
    # Fixed500-neighbor parental-fitness context; the target interaction is known.
    candidates=paired.drop(index=['DYNLL1','DYNLL2'],errors='ignore').copy();candidates['parental_distance']=abs(candidates.parental_LFC-paired.loc['DYNLL1','parental_LFC']);matched=candidates.sort_values(['parental_distance'],kind='stable').head(500);matched.to_csv(R/'hTERT_parental_fitness_matched500.csv');target=float(paired.loc['DYNLL1','mean_delta']);rank=int((matched.mean_delta<target).sum())+1
    # Threshold checks repeat the same within-guide interaction estimator.
    threshold_results=[]
    for cutoff in [0,20,30,50,100]:
        mask=(joint.denominator_mean_WT>=cutoff)&(joint.denominator_mean_KO>=cutoff)
        for context in ['WT','KO']:
            cols=[f'denominator_r{i}_{context}' for i in [1,2]];mask&=(joint[cols].to_numpy()>=10).sum(axis=1)>=1
        q=joint[mask].copy()
        for i in [1,2]:q[f'delta_r{i}']=q[f'LFC_r{i}_KO']-q[f'LFC_r{i}_WT']
        gene=q.groupby('gene')[['delta_r1','delta_r2']].median();gene['mean_delta']=gene.mean(axis=1);gene['percentile']=gene.mean_delta.rank(pct=True)*100
        for g in ['DYNLL1','DYNLL2']:
            v=gene.loc[g];threshold_results.append({'gene':g,'threshold':cutoff,'eligible_guides':int(q.gene.eq(g).sum()),'replicate1':float(v.delta_r1),'replicate2':float(v.delta_r2),'mean_delta':float(v.mean_delta),'discovery_percentile':float(v.percentile),'genes_evaluated':len(gene)})
    pd.DataFrame(threshold_results).to_csv(R/'hTERT_interaction_threshold_sensitivity.csv',index=False)
    # All network edges, not a selected pathway subset.
    aspects={k:v for obj in json.loads((ROOT/'phase9/data/U2OS_APMS.cx').read_text(encoding='utf-8')) for k,v in obj.items()};names={x['@id']:x['n'] for x in aspects['nodes']};adj={g:set() for g in names.values()}
    for e in aspects['edges']:
        a,b=names[e['s']],names[e['t']]
        if a!=b:adj[a].add(b);adj[b].add(a)
    degree=pd.DataFrame({'gene':list(adj),'degree':[len(v) for v in adj.values()]});degree.to_csv(R/'APMS_all_node_degree.csv',index=False);groups=pd.read_csv(OLD/'APMS_paralog_neighbors.csv')
    # Systematic identifiers resolve legacy CRISPR symbols without choosing
    # individual genes by their measured effect. Ambiguous mappings remain out.
    official=pd.read_csv(ROOT/'phase10/data/Brunello_official_library.txt',sep='\t',low_memory=False);libids={}
    for gid,g in official.dropna(subset=['Target Gene ID']).groupby('Target Gene ID'):
        symbols=set(g['Target Gene Symbol'].replace({'H2AFX':'H2AX','MRE11A':'MRE11'}));symbols&=set(paired.index)
        if len(symbols)==1:libids[int(gid)]=next(iter(symbols))
    hgnc=pd.read_csv(P/'data/hgnc_complete_set.txt',sep='\t',low_memory=False);approved={};aliases={}
    for row in hgnc.dropna(subset=['entrez_id']).to_dict('records'):
        try:gid=int(row['entrez_id'])
        except (ValueError,TypeError):continue
        approved[row['symbol']]=gid
        for col in ['alias_symbol','prev_symbol']:
            if isinstance(row.get(col),str):
                for alias in row[col].split('|'):aliases.setdefault(alias,set()).add(gid)
    with zipfile.ZipFile(ROOT/'phase10/data/CM4AI_PPIDownloader.zip') as z:
        raw=pd.read_csv(z.open(next(n for n in z.namelist() if n.endswith('edgeList.tsv'))),sep='\t')
    rawids=pd.concat([raw[['Symbol1','GeneID1']].rename(columns={'Symbol1':'gene','GeneID1':'entrez'}),raw[['Symbol2','GeneID2']].rename(columns={'Symbol2':'gene','GeneID2':'entrez'})]).drop_duplicates();rawidmap={g:int(q.entrez.iloc[0]) for g,q in rawids.groupby('gene') if q.entrez.nunique()==1};identity=[]
    for gene in degree.gene:
        gid=rawidmap.get(gene,approved.get(gene));method='exact_symbol' if gene in paired.index else None;countgene=gene if method else None
        if gid is None and len(aliases.get(gene,set()))==1:gid=next(iter(aliases[gene]))
        if countgene is None and gid in libids:countgene=libids[gid];method='source_or_HGNC_Entrez_to_official_library'
        identity.append({'gene':gene,'entrez_id':gid,'count_gene':countgene,'mapping_method':method or 'unresolved'})
    identity=pd.DataFrame(identity);duplicate=identity.count_gene.notna()&identity.count_gene.duplicated(keep=False);identity['ambiguous_duplicate_count_gene']=duplicate;identity.to_csv(R/'APMS_gene_identity_map.csv',index=False);usable=identity[identity.count_gene.notna()&~duplicate].merge(degree,on='gene',validate='one_to_one');whole=usable.merge(paired,left_on='count_gene',right_index=True,validate='one_to_one').set_index('gene');eligibility=groups[['neighbor','class']].merge(identity,left_on='neighbor',right_on='gene',how='left',validate='one_to_one');eligibility['included_in_gene_set_analysis']=eligibility.neighbor.isin(whole.index);eligibility.to_csv(R/'APMS_group_eligibility.csv',index=False)
    neighbors=set(groups.neighbor)|{'DYNLL1','DYNLL2'};background=whole[~whole.index.isin(neighbors)].copy();features=np.column_stack([whole.parental_LFC,np.log1p(whole.degree)]);scale=np.nanstd(features,axis=0,ddof=1);group_rows=[];balance=[];membership=[];nullrows=[]
    for gi,group in enumerate(['DYNLL1_only','shared','DYNLL2_only']):
        nominated=groups.loc[groups['class'].eq(group),'neighbor'];observed=whole.loc[whole.index.intersection(nominated)].drop(index=['DYNLL1','DYNLL2'],errors='ignore');pools=[];keep=[]
        for gene,row in observed.iterrows():
            candidates=background[background.guides.eq(row.guides)]
            if len(candidates)<50:continue
            diff=np.column_stack([candidates.parental_LFC-row.parental_LFC,np.log1p(candidates.degree)-np.log1p(row.degree)])/scale;order=np.argsort((diff**2).sum(axis=1),kind='stable')[:50];pools.append(candidates.iloc[order].index.to_numpy());keep.append(gene)
        observed=observed.loc[keep];pool=np.array(pools);assert len(observed)>0;rng=np.random.default_rng(SEED+gi);draw=pool[np.arange(len(pool))[None,:],rng.integers(0,50,size=(20000,len(pool)))];sorted_draw=np.sort(draw,axis=1);duplicate=(sorted_draw[:,1:]==sorted_draw[:,:-1]).any(axis=1);iterations=0
        while duplicate.any():
            n=int(duplicate.sum());draw[duplicate]=pool[np.arange(len(pool))[None,:],rng.integers(0,50,size=(n,len(pool)))];sorted_draw=np.sort(draw,axis=1);duplicate=(sorted_draw[:,1:]==sorted_draw[:,:-1]).any(axis=1);iterations+=1;assert iterations<300
        effect_map=whole.mean_delta.to_dict();flat=np.array([effect_map[x] for x in draw.ravel()]).reshape(draw.shape);null=np.median(flat,axis=1);obs=float(observed.mean_delta.median());p2=min(1.,2*min((np.sum(null<=obs)+1)/(len(null)+1),(np.sum(null>=obs)+1)/(len(null)+1)));group_rows.append({'group':group,'released_neighbors':len(nominated),'eligible_matched_neighbors':len(observed),'observed_median_interaction':obs,'matched_null_median':float(np.median(null)),'matched_null_2p5':float(np.quantile(null,.025)),'matched_null_97p5':float(np.quantile(null,.975)),'exploratory_gene_set_p_two_sided':float(p2),'draws':len(null)});nullrows.append(pd.DataFrame({'group':group,'draw':np.arange(len(null)),'median_interaction':null}))
        membership.append(observed.reset_index().assign(group=group))
        for variable,transform in [('parental_LFC',lambda x:x),('degree',np.log1p)]:
            vv=transform(observed[variable].to_numpy());lookup=transform(whole[variable]).to_dict();mv=np.array([lookup[x] for x in draw[:2000].ravel()]);bg=transform(background[variable].to_numpy());sd=float(np.std(bg,ddof=1));balance.append({'group':group,'feature':variable,'observed_mean':float(vv.mean()),'matched_mean':float(mv.mean()),'background_mean':float(bg.mean()),'SMD_before':float((vv.mean()-bg.mean())/sd),'SMD_after':float((vv.mean()-mv.mean())/sd)})
        print('Network group',group,'eligible',len(observed),'p',p2,flush=True)
    result=pd.DataFrame(group_rows);result['exploratory_BH_q']=multipletests(result.exploratory_gene_set_p_two_sided,method='fdr_bh')[1];result.to_csv(R/'APMS_hTERT_matched_gene_sets.csv',index=False);pd.concat(nullrows).to_csv(R/'APMS_hTERT_gene_set_nulls.csv.gz',index=False,compression='gzip');pd.concat(membership).to_csv(R/'APMS_hTERT_group_members.csv',index=False);pd.DataFrame(balance).to_csv(R/'APMS_hTERT_matching_balance.csv',index=False)
    # Local-locus mismatch diagnostics, keeping exact reference bounds/provenance.
    official=pd.read_csv(ROOT/'phase10/data/Brunello_official_library.txt',sep='\t',low_memory=False);targetguides=official[official['Target Gene Symbol'].eq('DYNLL1')];local=[]
    for row in targetguides.to_dict('records'):
        query=row['sgRNA Target Sequence'];rec={'guide_sequence':query,'cut_coordinate':int(row['Position of Base After Cut (1-based)'])}
        for gene in ['DYNLL1','DYNLL2']:
            seq=(ROOT/f'phase9/data/reagent_identity/{gene}_sequence.txt').read_text(encoding='utf-8').strip().upper();arr=np.frombuffer(seq.encode(),dtype='S1');windows=np.lib.stride_tricks.sliding_window_view(arr,len(query));best=100
            for strand,q in [('+',query),('-',rc(query))]:
                distance=(windows!=np.frombuffer(q.encode(),dtype='S1')).sum(axis=1);m=int(distance.min());best=min(best,m)
            rec[gene+'_minimum_local_mismatches']=best
        local.append(rec)
    pd.DataFrame(local).to_csv(R/'Brunello_DYNLL1_local_paralog_mismatches.csv',index=False)
    # Sample-level count diagnostics are re-extracted without any gene filtering.
    countq=[]
    for acc,file,sep,noncols in [('GSE232592','GSE232592_CRISPR_p53_counts.txt.gz','\t',['ID','Gene']),('GSE290538','GSE290538_U2OS_HRE_Count_table.txt.gz','\t',['sgRNA','Gene']),('GSE283147','GSE283147_combined-guide-level-counts.csv.gz',',',['ID','GENE','SEQUENCE','GUIDE_ID','TARGETEXON']),('GSE233195','GSE233195_JSMO.count.txt.gz','\t',['sgRNA','Gene'])]:
        a=pd.read_csv(ROOT/'phase9/data'/file,sep=sep,low_memory=False)
        for col in a.columns.difference(noncols):
            x=a[col].to_numpy(float);countq.append({'accession':acc,'sample':col,'guides':len(x),'total_reads':float(x.sum()),'zero_fraction':float((x==0).mean()),'median_reads':float(np.median(x)),'p10_reads':float(np.quantile(x,.1)),'p90_reads':float(np.quantile(x,.9))})
    pd.DataFrame(countq).to_csv(R/'screen_library_QC.csv',index=False)
    # Keep published standardized scores separate from unscaled count ratios.
    published=pd.read_excel(ROOT/'phase9/data/GSE232592_CRISPR_p53_GeneLevel.xlsx',sheet_name='fold-change').rename(columns={'Gene':'gene'});new=pd.read_csv(OLD/'GSE232592_WT_growth_gene_effects.csv.gz');published=published.merge(new[['gene','mean_LFC']],on='gene',how='inner',validate='one_to_one');published[['gene','WT_DMSOvT0_L2FC','WT_DMSOvT0_FDR','mean_LFC']].to_csv(R/'published_standardized_vs_new_growth.csv.gz',index=False,compression='gzip')
    summary={'hTERT_fitness_matched_background':len(matched),'hTERT_DYNLL1_rank_within_matched_plus_target':rank,'hTERT_DYNLL1_interaction':target,'hTERT_replicate_interaction_spearman':float(spearmanr(paired.delta_r1,paired.delta_r2).statistic),'network_group_tests':result.to_dict('records'),'local_mismatch_scope':'Retained DYNLL1/DYNLL2 loci only; no genome-wide off-target claim.','screen_libraries':len(countq),'published_vs_unscaled_growth_spearman':float(spearmanr(published.WT_DMSOvT0_L2FC,published.mean_LFC).statistic),'all_exploratory':True};(R/'context_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
