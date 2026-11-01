"""Post hoc hTERT-context audit prompted by the new exploratory count result."""
from pathlib import Path
import json,re
import numpy as np
import pandas as pd
from bs4 import BeautifulSoup
P=Path(__file__).resolve().parent;R=P/'results'

def main():
    library=pd.read_csv(P/'data/Brunello_official_library.txt',sep='\t',low_memory=False)
    library=library[library['Target Gene ID'].notna()&library['Position of Base After Cut (1-based)'].notna()].copy()
    library['key']=library['Target Gene ID'].astype(int).astype(str)+'_'+library['Target Gene Symbol']+'_'+library['Position of Base After Cut (1-based)'].astype(int).astype(str)
    # Collapse repeated transcript annotations only when every relevant value agrees.
    rows=[]
    for key,g in library.groupby('key',sort=False):
        row={'guide':key,'library_records':len(g),'sequence_count':g['sgRNA Target Sequence'].nunique()}
        for old,new in [('sgRNA Target Sequence','sequence'),('Genomic Sequence','chrom_accession'),('Position of Base After Cut (1-based)','cut_position')]:
            vals=g[old].dropna().unique();row[new]=vals[0] if len(vals)==1 else None
        rows.append(row)
    lib=pd.DataFrame(rows)
    wt=pd.read_parquet(R/'guides/GSE226994_parental_growth.parquet');ko=pd.read_parquet(R/'guides/GSE226994_hTERT_mutant_growth.parquet');q=wt.merge(ko,on=['guide','gene'],suffixes=('_WT','_KO'),validate='one_to_one').merge(lib,on='guide',how='left',validate='one_to_one')
    q['eligible_joint']=q.eligible_WT&q.eligible_KO
    for i in [1,2]:q[f'delta_r{i}']=q[f'LFC_r{i}_KO']-q[f'LFC_r{i}_WT']
    genes=q[q.eligible_joint].groupby('gene')[['delta_r1','delta_r2']].median();genes['mean_delta']=genes[['delta_r1','delta_r2']].mean(axis=1);genes['guides']=q[q.eligible_joint].groupby('gene').size();genes['depletion_percentile']=genes.mean_delta.rank(pct=True)*100;genes['chrom_accession']=q.groupby('gene').chrom_accession.agg(lambda s: s.dropna().iloc[0] if len(s.dropna().unique())==1 else None);genes['cut_position']=q.groupby('gene').cut_position.median();genes.to_csv(R/'hTERT_all_gene_interactions.csv.gz',compression='gzip')
    overlap=pd.read_csv(P.parent/'phase9/results/public_CRISPR_guide_overlap_with_JUMP.csv');lookup=overlap[overlap.source.eq('GSE290538')].set_index('sequence');target=q[q.gene.isin(['DYNLL1','DYNLL2'])].copy();target['maximum_JUMP_overlap_bp']=target.sequence.map(lookup.maximum_JUMP_target_overlap_bp);target['identical_JUMP_sequence']=target.sequence.map(lookup.identical_JUMP_sequence);target.to_csv(R/'hTERT_DYNLL1_DYNLL2_reagent_and_interaction_audit.csv',index=False)
    # Reproduce the earlier all-guide gene summary with canonicalized symbols.
    tg=target[target.gene.eq('DYNLL1')&target.eligible_joint];v=tg[['delta_r1','delta_r2']].to_numpy();non=tg[tg.maximum_JUMP_overlap_bp.eq(0)];loo=[np.median(np.delete(v,i,axis=0),axis=0).mean() for i in range(len(v))]
    targetpos=float(genes.loc['DYNLL1','cut_position']);regional=genes[(genes.chrom_accession=='NC_000012.12')&(abs(genes.cut_position-targetpos)<=2e6)&(genes.index!='DYNLL1')];regional.to_csv(R/'hTERT_chr12_within_2Mb_controls.csv')
    # Retain and check known assay context genes without selecting them by outcome.
    controls=['FANCA','FANCB','FANCC','FANCD2','FANCE','FANCF','FANCG','FANCI','FANCL','FANCM','BRCA1','BRCA2','PALB2','PARP1','TERT']
    genes.reindex(controls).to_csv(R/'hTERT_published_pathway_context.csv')
    record={'status':'Post hoc candidate from exploratory reanalysis; not independent confirmation.','official_library_rows':len(library),'unique_official_keys':len(lib),'ambiguous_sequence_keys':int(lib.sequence_count.gt(1).sum()),'input_guides':len(q),'exact_official_key_matches':int(q.sequence.notna().sum()),'DYNLL1_guides':len(tg),'DYNLL1_mutant_growth_LFC':float(pd.read_csv(R/'DYNLL1_DYNLL2_all_contrasts.csv').query("gene=='DYNLL1' and contrast=='GSE226994_hTERT_mutant_growth'").mean_LFC.iloc[0]),'DYNLL1_paired_guide_interaction_mean':float(v.mean(axis=0).mean()),'DYNLL1_paired_median_interaction_mean':float(np.median(v,axis=0).mean()),'replicate_median_interactions':np.median(v,axis=0).tolist(),'interaction_depletion_percentile':float(genes.loc['DYNLL1','depletion_percentile']),'interaction_LOO_range':[float(min(loo)),float(max(loo))],'guides_with_negative_mean_interaction':int((v.mean(axis=1)<0).sum()),'guides_with_negative_interaction_in_both_replicates':int((v<0).all(axis=1).sum()),'exact_JUMP_shared_sequences':int(tg.identical_JUMP_sequence.notna().sum()),'overlapping_JUMP_intervals':int(tg.maximum_JUMP_overlap_bp.gt(0).sum()),'nonoverlapping_guides':len(non),'nonoverlapping_interactions':non[['guide','sequence','delta_r1','delta_r2']].to_dict('records'),'regional_control_genes':len(regional),'regional_control_median_delta':float(regional.mean_delta.median()),'regional_controls_below_DYNLL1':int((regional.mean_delta<genes.loc['DYNLL1','mean_delta']).sum()),'region_scope':'Descriptive within 2Mb; not copy-number correction or a formal independence test.','clone_scope':'GEO labels two replicates of hTERT-CRISPR cells but does not assign both screened replicates to separate independently edited clones.','inference':'A promising hTERT-context fitness candidate. Most guides overlap JUMP intervals, one guide does not deplete, and the sole nonoverlapping guide has mixed differential effects across replicates. No nuclear alignment assay or rescue is provided.'}
    (R/'hTERT_candidate_audit.json').write_text(json.dumps(record,indent=2))
    # Published Honeywell L2RV values are standardized, so their numerical scale
    # cannot be compared directly with newly computed log2 abundance ratios.
    soup=BeautifulSoup((P/'data/Honeywell_primary.html').read_text(encoding='utf-8'),'html.parser');text=soup.get_text(' ',strip=True);needle='Each fold change value was z-scored';i=text.find(needle);assert i>=0
    (R/'published_score_scale_clarification.json').write_text(json.dumps({'source':'https://doi.org/10.1038/s41589-024-01584-7','published_column_labels':'WT_DMSOvT0_L2FC and corresponding etoposide/p53KO columns','clarification':'The methods state that gene-level L2RV values were z-scored relative to non-targeting pseudo-genes. Do not interpret the released values as unscaled log2 fold changes. Original published FDRs remain unchanged. New transparent median-guide log2 ratios are different estimands.','supersedes':'Any earlier prose equating the magnitude of those published columns directly with unscaled log2 fold change.'},indent=2))
    print(json.dumps(record,indent=2))
if __name__=='__main__':main()
