"""Reproducible metadata, published-result and guide-overlap extraction.

No new differential-expression, CRISPR significance or pathway-enrichment test.
"""
from pathlib import Path
import json,re,gzip,hashlib
import numpy as np
import pandas as pd
P=Path(__file__).resolve().parent;D=P/'data';R=P/'results';R.mkdir(exist_ok=True)
def load(f):return json.loads((D/f).read_text(encoding='utf-8'))
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def metadata(acc):
    text=(D/(acc+'_metadata.txt')).read_text(encoding='utf-8');return {'accession':acc,'title':next(x.split(' = ',1)[1] for x in text.splitlines() if x.startswith('!Series_title = ')),'PMIDs':'|'.join(x.split(' = ',1)[1] for x in text.splitlines() if x.startswith('!Series_pubmed_id = ')),'U2OS_model_verified':True}
def main():
    # Published gene-level results, with original multiple-testing correction.
    gene_results=[]
    for sheet in ['rate-based','fold-change']:
        a=pd.read_excel(D/'GSE232592_CRISPR_p53_GeneLevel.xlsx',sheet_name=sheet);q=a[a.Gene.eq('DYNLL1')].copy();assert len(q)==1;q.insert(0,'table',sheet);gene_results+=q.to_dict('records')
    pd.DataFrame(gene_results).to_csv(R/'GSE232592_DYNLL1_published_gene_results.csv',index=False)
    # Complete processed count files are retained; export only DYNLL1 records.
    coverage=[]
    specs=[('GSE232592','GSE232592_CRISPR_p53_counts.txt.gz','\t','Gene',3),('GSE290538','GSE290538_U2OS_HRE_Count_table.txt.gz','\t','Gene',3),('GSE283147','GSE283147_combined-guide-level-counts.csv.gz',',','GENE',2),('GSE233195','GSE233195_JSMO.count.txt.gz','\t','Gene',2)]
    for acc,file,sep,gene,nrep in specs:
        a=pd.read_csv(D/file,sep=sep,low_memory=False);q=a[a[gene].eq('DYNLL1')].copy();q.to_csv(R/(acc+'_DYNLL1_guide_counts.csv'),index=False);coverage.append({**metadata(acc),'retained_table_rows':len(a),'DYNLL1_guides':len(q),'replicates_per_arm':nrep,'count_scale':'published counts','target_sequences_in_table':'SEQUENCE' in a or acc=='GSE290538','status':'Published DYNLL1 FDRs all above .05' if acc=='GSE232592' else 'Guide coverage verified; no new screen inference in this search'})
    x=pd.ExcelFile(D/'GSE226994_MAGeCK_count_normalized.xlsx');parts=[]
    for sheet in x.sheet_names:
        a=pd.read_excel(x,sheet_name=sheet);q=a[a.Gene.eq('DYNLL1')].copy();q.insert(0,'sheet',sheet);parts.append(q)
    q=pd.concat(parts,ignore_index=True);q.to_csv(R/'GSE226994_DYNLL1_guide_counts.csv',index=False);coverage.append({**metadata('GSE226994'),'retained_table_rows':len(a),'DYNLL1_guides':int(q.groupby('sheet').size().min()),'replicates_per_arm':2,'count_scale':'published normalized counts, not raw integers','target_sequences_in_table':False,'status':'U2OS and engineered hTERT arms kept separate; no new screen inference'})
    # An explicitly versioned archived reference, not a claim about current DepMap.
    effect=pd.read_csv(D/'DepMap24Q4_U2OS_CRISPRGeneEffect.csv');dep=pd.read_csv(D/'DepMap24Q4_U2OS_CRISPRGeneDependency.csv');assert len(effect)==len(dep)==1;av=pd.read_csv(D/'DepMap24Q4_AvanaGuideMap.csv',low_memory=False);av=av[av.Gene.eq('DYNLL1 (8655)')].copy();av.to_csv(R/'DepMap24Q4_DYNLL1_Avana_guides.csv',index=False);dm={'release':'DepMap 24Q4 Public','doi':'10.25452/figshare.plus.27993248','model':'ACH-000364','screen':'SC-000364.AV01','library':'Avana','days':21,'passing_replicates':2,'GeneEffect_DYNLL1':float(effect['DYNLL1 (8655)'].iloc[0]),'Dependency_probability_DYNLL1':float(dep['DYNLL1 (8655)'].iloc[0]),'guide_map_DYNLL1_records':len(av),'guides_used_by_Chronos':int(av.UsedByChronos.sum()),'interpretation':'Modest negative fitness score; dependency probability is not a p-value. No strong dependency conclusion or spatial validation. One retained Avana guide is identical to JUMP.'};(R/'DepMap24Q4_U2OS_DYNLL1.json').write_text(json.dumps(dm,indent=2));coverage.append({'accession':'DepMap24Q4','title':'Archived U2OS basal fitness screen','PMIDs':'','U2OS_model_verified':True,'DYNLL1_guides':int(av.UsedByChronos.sum()),'replicates_per_arm':2,'count_scale':'Chronos gene score','target_sequences_in_table':True,'status':'Modest negative fitness; not spatial validation'})
    pd.DataFrame(coverage).to_csv(R/'verified_U2OS_CRISPR_sources.csv',index=False)
    # Guide identity and locus overlap: report base-pair overlap, without inventing
    # a binary independence threshold or performing a genome-wide off-target test.
    jump=pd.read_csv(P.parent/'phase7/results/DYNLL1_shared_reagent_identity.csv');seq=(D/'reagent_identity/DYNLL1_sequence.txt').read_text();annotation=json.loads((D/'reagent_identity/provenance.json').read_text())[0];start=annotation['start_0based'];allguides=[{'source':'JUMP','id':x['Component Catalog Number'],'sequence':x['Sequence'],'used_in_primary_score':True} for x in jump.to_dict('records')]
    h=pd.read_csv(R/'GSE290538_DYNLL1_guide_counts.csv');allguides +=[{'source':'GSE290538','id':x.sgRNA,'sequence':x.sgRNA.split('_',1)[1],'used_in_primary_score':None} for x in h.itertuples()]
    t=pd.read_csv(R/'GSE283147_DYNLL1_guide_counts.csv');allguides +=[{'source':'GSE283147','id':x.ID,'sequence':x.SEQUENCE,'used_in_primary_score':None} for x in t.itertuples()]
    allguides +=[{'source':'DepMap24Q4','id':x.sgRNA,'sequence':x.sgRNA,'used_in_primary_score':bool(x.UsedByChronos)} for x in av.itertuples()]
    for g in allguides:
        matches=[]
        for strand,query in [('+',g['sequence']),('-',rc(g['sequence']))]:
            for hit in re.finditer('(?='+re.escape(query)+')',seq):matches.append((start+hit.start(),strand))
        g['DYNLL1_locus_exact_matches']=len(matches);g['locus_start_0based']=matches[0][0] if len(matches)==1 else None;g['strand']=matches[0][1] if len(matches)==1 else None;g['all_locus_matches']='|'.join(str(a)+b for a,b in matches)
    for g in allguides:
        g['identical_JUMP_sequence']='|'.join(x['id'] for x in allguides if x['source']=='JUMP' and g['sequence']==x['sequence']);overlaps=[]
        for j in allguides[:4]:
            if g['locus_start_0based'] is None or j['locus_start_0based'] is None:continue
            overlap=max(0,min(g['locus_start_0based']+20,j['locus_start_0based']+20)-max(g['locus_start_0based'],j['locus_start_0based']));overlaps.append((overlap,j['id']))
        g['maximum_JUMP_target_overlap_bp']=max([x[0] for x in overlaps] or [0]);g['nearest_overlapping_JUMP_constructs']='|'.join(x[1] for x in overlaps if x[0]==g['maximum_JUMP_target_overlap_bp'] and x[0]>0)
    pd.DataFrame(allguides).to_csv(R/'public_CRISPR_guide_overlap_with_JUMP.csv',index=False)
    # Direct AP-MS associations and the explicitly ambiguous shared IF image.
    aspects={k:v for obj in load('U2OS_APMS.cx') for k,v in obj.items()};nodes={x['@id']:x['n'] for x in aspects['nodes']};ids={v:k for k,v in nodes.items() if v in ['DYNLL1','DYNLL2']};edges=[];sets={}
    for gene,i in ids.items():
        ee=[e for e in aspects['edges'] if e['s']==i or e['t']==i];sets[gene]={nodes[e['t'] if e['s']==i else e['s']] for e in ee}
        edges +=[{'query_gene':gene,'source':nodes[e['s']],'target':nodes[e['t']],'edge_id':e['@id'],'association':'AP-MS network association; not proof of direct binding or perturbation causality'} for e in ee]
    pd.DataFrame(edges).to_csv(R/'U2OS_APMS_DYNLL1_DYNLL2_associations.csv',index=False);attrs=[x for x in aspects['nodeAttributes'] if x['po'] in ids.values()];(R/'U2OS_LC8_localization_attributes.json').write_text(json.dumps(attrs,indent=2));network={'network_nodes':len(nodes),'network_edges':len(aspects['edges']),'DYNLL1_unique_neighbors':len(sets['DYNLL1']),'DYNLL2_unique_neighbors':len(sets['DYNLL2']),'shared_neighbors':len(sets['DYNLL1']&sets['DYNLL2']),'DYNLL1_selected_associations':sorted(sets['DYNLL1']&{'ATMIN','TP53BP1','FAM83D','HMMR','TLK1','TLK2','ZMYM2','ZMYM3','ZMYM4','BACH1'}),'image_gene_specificity':'Both genes share HPA039954 and an ambiguous ENSG00000088986,ENSG00000264364 annotation; do not use this image as DYNLL1-specific localization validation.','scope':'U2OS protein association context using affinity baits, not a DYNLL1 loss-of-function contrast. No new interaction significance or pathway test.'};(R/'U2OS_APMS_summary.json').write_text(json.dumps(network,indent=2))
    # qChIP-MS reports a protein group, not gene-specific DYNLL1 abundance.
    q=pd.read_excel(D/'41467_2026_73609_MOESM5_ESM.xlsx');q=q[q['Gene names'].eq('DYNLL1')];assert len(q)==1;q.to_csv(R/'qChIP_U2OS_DYNLL1_DYNLL2_protein_group.csv',index=False)
    controls=[];control_summary=[]
    for number in [10,11,12]:
        x=pd.ExcelFile(D/f'41467_2026_73609_MOESM{number}_ESM.xlsx')
        for sheet in x.sheet_names:
            a=pd.read_excel(x,sheet_name=sheet);q=a[a.astype(str).apply(lambda s:s.str.contains('DYNLL|P63167|Q96FJ2',regex=True,na=False)).any(axis=1)].copy();q.insert(0,'sheet',sheet);q.insert(0,'supplementary_data',number-2);controls.append(q);control_summary.append({'supplementary_data':number-2,'sheet':sheet,'LC8_group_rows':len(q),'source':'ZBTB48 WT versus KO, not DYNLL1 perturbation','mean_log2_fold_changes':q.mean_log2_fold_change.tolist()})
    pd.concat(controls,ignore_index=True).to_csv(R/'qChIP_U2OS_LC8_KO_control_rows.csv',index=False);(R/'qChIP_U2OS_control_summary.json').write_text(json.dumps(control_summary,indent=2))
    # Catalog scope and verified exclusions.
    catalog={}
    for file in ['CREEDS_Gene_Perturbations_Gene_Matrix_Transposed.gmt.gz','PerturbAtlas_Human_Gene_Perturbations_Gene_Matrix_Transposed.gmt.gz']:
        with gzip.open(D/file,'rt',encoding='utf-8') as stream:terms=[line.split('\t',1)[0] for line in stream]
        catalog[file]={'signature_terms':len(terms),'DYNLL1_target_terms':[t for t in terms if re.search(r'(^|::)DYNLL1(::|\b)',t,re.I)],'U2OS_labeled_terms':sum(bool(re.search('U[ -]?2[ -]?OS',t,re.I)) for t in terms)}
    cp=pd.read_parquet(D/'L1000_XPR__Gene_Knockouts__Cell_Line.json');catalog['L1000_XPR_catalog']={'cell_models':len(cp),'U2OS_in_released_cell_metadata':bool(cp.astype(str).apply(lambda s:s.str.contains('U[ -]?2[ -]?OS',regex=True,case=False,na=False)).any(axis=1).any())}
    sc=pd.read_csv(D/'scperturb_dataset_info_140922.csv');catalog['scPerturb_historical_CSV']={'source':'sanderlab/scPerturb master metadata/scperturb_dataset_info_140922.csv','dataset_rows':len(sc),'U2OS_rows':int(sc.astype(str).apply(lambda s:s.str.contains('U[ -]?2[ -]?OS',regex=True,case=False,na=False)).any(axis=1).sum()),'scope':'Historical retained CSV, not a claim about all current single-cell data.'}
    inst=pd.read_csv(D/'GSE92742_Broad_LINCS_inst_info.txt.gz',sep='\t',low_memory=False);sig=pd.read_csv(D/'GSE92742_Broad_LINCS_sig_info.txt.gz',sep='\t',low_memory=False);ii=inst[inst.pert_iname.eq('DYNLL1')];ss=sig[sig.pert_iname.eq('DYNLL1')];catalog['GSE92742_DYNLL1']={'instances':len(ii),'signatures':len(ss),'cell_models':sorted(ii.cell_id.unique().tolist()),'U2OS_instances':int(ii.cell_id.eq('U2OS').sum()),'U2OS_signatures':int(ss.cell_id.eq('U2OS').sum())};ii.to_csv(R/'LINCS_DYNLL1_instances_all_cell_metadata.csv',index=False);ss.to_csv(R/'LINCS_DYNLL1_signatures_all_cell_metadata.csv',index=False)
    (R/'perturbation_catalog_scope.json').write_text(json.dumps(catalog,indent=2))
    # Full source indexes, not just selected positive studies.
    geo=load('geo_U2OS_CRISPR_summaries.json');pd.DataFrame([geo['result'][k] for k in geo['result']['uids']])[['accession','title','summary','n_samples']].to_csv(R/'GEO_U2OS_CRISPR_series.csv',index=False)
    pmc=load('europepmc_DYNLL1_U2OS.json');pd.DataFrame(pmc['resultList']['result'])[['id','pmcid','title','firstPublicationDate','doi']].to_csv(R/'EuropePMC_DYNLL1_U2OS_candidates.csv',index=False)
    summary={'primary_system':'U2OS','new_CRISPR_source_count':6,'CRISPR_sources':coverage,'DepMap':dm,'APMS':network,'GEO_U2OS_CRISPR_series_returned':len(geo['result']['uids']),'EuropePMC_candidate_hits':pmc['hitCount'],'verified_DYNLL1_U2OS_Perturb_seq_dataset':None,'no_new_statistical_tests':True,'scope':'Search and source-result extraction. New guide count resources are not positive validation. No nuclear-neighbor endpoint, rescue of that endpoint, or mechanistic link has been obtained.'};(R/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps({'CRISPR_datasets':6,'DepMap':dm,'APMS':network},indent=2))
if __name__=='__main__':main()
