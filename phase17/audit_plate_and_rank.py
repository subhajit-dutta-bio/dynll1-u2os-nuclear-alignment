"""Post hoc layout and rank diagnostics without changing primary measurements."""
from pathlib import Path
import hashlib,itertools,json,re
import numpy as np,pandas as pd
P=Path(__file__).resolve().parent;ROOT=P.parent;R=P/'results'
NAMED=['DYNLL1','ANAPC5','SPPL3','CAMKK2','P2RX7','P2RX4','ORAI1']
CENT=['HSPB8','PRKAB1','CIT']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    R.mkdir(exist_ok=True);m=pd.read_csv(ROOT/'phase14/frozen_contrasts.csv');raw=pd.read_csv(ROOT/'phase14/results/alignment_contrasts.csv');e=raw[raw.model.eq('adaptive_conditional')];dates=sorted(e.batch.unique());assert len(m)==335
    # Verify frozen assignments against the released well, plate and reagent metadata.
    plates=pd.read_csv(ROOT/'data/metadata/plate.csv.gz');wanted=set(m.target_plate);pieces=[]
    for c in pd.read_csv(ROOT/'data/metadata/well.csv.gz',chunksize=250000):
        z=c[c.Metadata_Source.eq('source_13')&c.Metadata_Plate.isin(wanted)]
        if len(z):pieces.append(z)
    public=pd.concat(pieces).merge(plates,on=['Metadata_Source','Metadata_Plate'],validate='many_to_one')
    public=public.rename(columns={'Metadata_Plate':'target_plate','Metadata_Well':'target_well','Metadata_Batch':'batch','Metadata_JCP2022':'public_jcp2022'})
    verified=m.merge(public[['batch','target_plate','target_well','public_jcp2022','Metadata_PlateType']],on=['batch','target_plate','target_well'],validate='one_to_one');assert len(verified)==335 and verified.target_jcp2022.eq(verified.public_jcp2022).all() and verified.Metadata_PlateType.eq('CRISPR').all()
    crispr=pd.read_csv(ROOT/'data/metadata/crispr.csv.gz');symbol=crispr.set_index('Metadata_JCP2022').Metadata_Symbol;assert m.target_jcp2022.map(symbol).eq(m.gene).all()
    x=m.merge(e[['gene','batch','effect']],on=['gene','batch'],validate='one_to_one');x['layout']=x.target_plate.str.extract(r'-(\d+)$')[0];x['row']=x.target_well.str[0].map(lambda s:ord(s)-64);x['column']=x.target_well.str[1:].astype(int)
    w=pd.read_csv(ROOT/'phase14/results/adaptive_conditional_wells.csv').set_index(['batch','plate','well']).value;components=[];weights=[]
    for a in x.itertuples():
        terms=[(a.target_plate,a.target_well,1.)]+[(a.target_plate,b,-.5) for b in a.source_NT_wells.split(';')]+[(a.negative_plate,a.negative_well,-1.)]+[(a.negative_plate,b,.5) for b in a.anchor_wells.split(';')]
        assert len(terms)==6
        source=w.loc[(a.batch,a.target_plate,a.target_well)]-np.mean([w.loc[(a.batch,a.target_plate,b)] for b in a.source_NT_wells.split(';')]);pos=w.loc[(a.batch,a.negative_plate,a.negative_well)]-np.mean([w.loc[(a.batch,a.negative_plate,b)] for b in a.anchor_wells.split(';')]);assert np.isclose(source-pos,a.effect,rtol=1e-12,atol=1e-12)
        components.append((source,pos,source-pos))
        for plate in sorted({t[0] for t in terms}):
            total=sum(t[2] for t in terms if t[0]==plate);assert total==0
            weights.append({'gene':a.gene,'batch':a.batch,'plate':plate,'sum_contrast_weights':total})
    x[['source_only_contrast','position_component','reconstructed_primary']]=components
    x.to_csv(R/'plate_well_all_assignments.csv',index=False);pd.DataFrame(weights).to_csv(R/'contrast_plate_weights.csv',index=False)
    rows=[]
    for gene,z in x.groupby('gene',sort=True):
        assert z.layout.nunique()==z.target_well.nunique()==1
        rows.append({'gene':gene,'panel':z.panel.iloc[0],'layout':z.layout.iloc[0],'well':z.target_well.iloc[0],'row':z.row.iloc[0],'column':z.column.iloc[0],'dates':len(z),'physical_plates':z.target_plate.nunique(),'same_layout_and_well_all_dates':True,'primary_mean':z.effect.mean(),'source_only_mean':z.source_only_contrast.mean(),'source_only_positive_dates':int(z.source_only_contrast.gt(0).sum()),'primary_positive_dates':int(z.effect.gt(0).sum())})
    layout=pd.DataFrame(rows);layout.to_csv(R/'plate_layout_by_gene.csv',index=False)
    pairs=[]
    for date,z in x[x.gene.isin(NAMED)].groupby('batch'):
        for a,b in itertools.combinations(z.itertuples(),2):
            same=a.target_plate==b.target_plate;dr=abs(a.row-b.row);dc=abs(a.column-b.column)
            pairs.append({'batch':date,'gene_a':a.gene,'gene_b':b.gene,'plate_a':a.target_plate,'plate_b':b.target_plate,'well_a':a.target_well,'well_b':b.target_well,'same_physical_plate':same,'row_distance':dr,'column_distance':dc,'Chebyshev_distance_if_same_plate':max(dr,dc) if same else np.nan,'adjacent_including_diagonal':same and max(dr,dc)<=1})
    pairs=pd.DataFrame(pairs);pairs.to_csv(R/'named_gene_pair_layout_audit.csv',index=False);assert not pairs.adjacent_including_diagonal.any()
    v=e.pivot(index='gene',columns='batch',values='effect').sort_index();follow=v[dates[-2:]].mean(axis=1);rank=follow.rank(ascending=False,method='min');g=pd.read_csv(ROOT/'phase14/results/gene_context_and_covariates.csv');match=set(g.loc[g.panel.eq('matched_other_chromosome'),'gene']);rows=[];paired=[]
    for gene in v.index:
        delta=v.loc['DYNLL1']-v.loc[gene]
        rows.append({'gene':gene,'R4_effect':v.loc[gene,dates[3]],'R5_effect':v.loc[gene,dates[4]],'R4_R5_mean':follow.loc[gene],'descriptive_point_rank':int(rank.loc[gene]),'DYNLL1_minus_gene_R4':delta.iloc[3],'DYNLL1_minus_gene_R5':delta.iloc[4],'DYNLL1_minus_gene_followup_mean':delta.iloc[3:].mean(),'between_date_SD_all_five':v.loc[gene].std(ddof=1),'is_matched_other_chromosome':gene in match})
        for date in dates:paired.append({'gene':gene,'batch':date,'DYNLL1_effect':v.loc['DYNLL1',date],'comparator_effect':v.loc[gene,date],'paired_difference':delta.loc[date]})
    gaps=pd.DataFrame(rows).sort_values(['descriptive_point_rank','gene']);gaps.to_csv(R/'rank_gap_all_genes.csv',index=False);pd.DataFrame(paired).to_csv(R/'rank_pairwise_dates.csv',index=False)
    named=x[x.gene.isin(NAMED)];loc=pd.read_csv(ROOT/'phase15/results/all_chromosome12_positions_and_date_effects.csv');tel=set(loc.loc[loc.distance_from_DYNLL1_bp.between(0,1200000,inclusive='neither'),'gene']);local=x[x.gene.isin(tel)]
    out={'public_metadata_assignments_verified':len(verified),'named_genes':NAMED,'named_layouts':sorted(named.layout.unique()),'named_physical_plates':int(named.target_plate.nunique()),'named_same_plate_pair_comparisons_across_dates':int(pairs.same_physical_plate.sum()),'named_adjacent_pair_comparisons':int(pairs.adjacent_including_diagonal.sum()),'all_cohort_layouts':int(x.layout.nunique()),'telomeric_1_2Mb_layouts':sorted(local.layout.unique()),'all_335_contrasts_reconstructed':True,'additive_plate_weight_sums_all_zero':True,'repeated_layouts_across_dates':True,'DYNLL1_followup_mean':float(follow.loc['DYNLL1']),'matched_followup_median':float(follow.loc[sorted(match)].median()),'matched_all_five_median':float(v.loc[sorted(match)].mean(axis=1).median()),'DYNLL1_minus_TAF15_R4':float((v.loc['DYNLL1']-v.loc['TAF15']).iloc[3]),'DYNLL1_minus_TAF15_R5':float((v.loc['DYNLL1']-v.loc['TAF15']).iloc[4]),'DYNLL1_minus_TAF15_followup':float(follow.loc['DYNLL1']-follow.loc['TAF15']),'median_gene_between_date_SD':float(v.std(axis=1,ddof=1).median()),'scope':'Post hoc descriptive diagnostics. No rank significance, equivalence or independent plate randomization is claimed.'}
    (R/'plate_rank_audit_summary.json').write_text(json.dumps(out,indent=2))
    inputs=[ROOT/n for n in ['phase14/frozen_contrasts.csv','phase14/results/alignment_contrasts.csv','phase14/results/adaptive_conditional_wells.csv','data/metadata/well.csv.gz','data/metadata/plate.csv.gz','data/metadata/crispr.csv.gz']];inputs+=[P/'audit_plate_and_rank.py'];(P/'input_provenance.json').write_text(json.dumps({'inputs':[{'path':str(f.relative_to(ROOT)),'sha256':sha(f)} for f in inputs],'preserved_v4_documents':[{'path':f.name,'sha256':sha(f)} for f in ROOT.glob('*_v4.docx')]},indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
