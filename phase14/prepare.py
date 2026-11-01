"""Validate every target/control identity and freeze the new raw-image cohort."""
from pathlib import Path
import json,hashlib
from datetime import datetime,timezone
import pandas as pd
P=Path(__file__).resolve().parent;ROOT=P.parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    if (P/'frozen_plan.json').exists():print((P/'frozen_plan.json').read_text());return
    targets=pd.read_csv(ROOT/'phase13/results/proposed_target_wells.csv');groups=pd.read_csv(ROOT/'phase13/results/proposed_comparator_gene_sets.csv');assert len(groups)==67 and len(targets)==335
    plates=pd.read_csv(ROOT/'data/metadata/plate.csv.gz');well=pd.read_csv(ROOT/'data/metadata/well.csv.gz').merge(plates,on=['Metadata_Source','Metadata_Plate']);well=well[well.Metadata_Source.eq('source_13')&well.Metadata_PlateType.eq('CRISPR')];old=pd.read_csv(ROOT/'phase6/frozen_contrasts.csv');negative=old[['batch','negative_plate']].drop_duplicates().set_index('batch').negative_plate.to_dict();records={};contrasts=[]
    def add(r,role,gene='',panel=''):
        key=(r.Metadata_Plate,r.Metadata_Well);d={'batch':r.Metadata_Batch,'plate':r.Metadata_Plate,'well':r.Metadata_Well,'jcp2022':r.Metadata_JCP2022,'role':role,'gene':gene,'panel':panel}
        if key in records:assert records[key]==d,(records[key],d)
        records[key]=d
    for t in targets.itertuples():
        assert t.Metadata_Batch in negative
        check=well[well.Metadata_Plate.eq(t.Metadata_Plate)&well.Metadata_Well.eq(t.Metadata_Well)];assert len(check)==1 and check.Metadata_JCP2022.iloc[0]==t.Metadata_JCP2022
        add(t,'target',t.Metadata_Symbol,t.set);source=well[well.Metadata_Plate.eq(t.Metadata_Plate)&well.Metadata_Well.isin(['A02','I23'])];assert len(source)==2 and source.Metadata_JCP2022.eq('JCP2022_800002').all()
        for r in source.itertuples():add(r,'source_NT')
        neg=well[well.Metadata_Plate.eq(negative[t.Metadata_Batch])];assert len(neg)==384 and neg.Metadata_JCP2022.isin(['JCP2022_800001','JCP2022_800002']).all();pos=neg[neg.Metadata_Well.eq(t.Metadata_Well)];assert len(pos)==1;pos=next(pos.itertuples());add(pos,'position_control');anchors=neg[neg.Metadata_Well.isin(['A02','I23'])];assert len(anchors)==2 and anchors.Metadata_JCP2022.isin(['JCP2022_800001','JCP2022_800002']).all()
        for r in anchors.itertuples():add(r,'negative_anchor')
        contrasts.append({'gene':t.Metadata_Symbol,'panel':t.set,'batch':t.Metadata_Batch,'target_plate':t.Metadata_Plate,'target_well':t.Metadata_Well,'target_jcp2022':t.Metadata_JCP2022,'source_NT_wells':'A02;I23','negative_plate':pos.Metadata_Plate,'negative_well':pos.Metadata_Well,'position_control_jcp2022':pos.Metadata_JCP2022,'anchor_wells':'A02;I23'})
    fields=[]
    for row in records.values():
        for site in range(9):
            f=dict(row,site=site);f['field']=f"{f['plate']}__{f['well']}__s{site}";f['dna_url']=f"https://cellpainting-gallery.s3.amazonaws.com/cpg0016-jump/source_13/images/{f['batch']}/images/{f['plate']}/{f['plate']}_{f['well']}_T0001F{site+1:03d}L01A01Z01C01.tif";fields.append(f)
    d=pd.DataFrame(fields).sort_values('field');d.to_csv(P/'frozen_fields.csv',index=False);pd.DataFrame(contrasts).sort_values(['gene','batch']).to_csv(P/'frozen_contrasts.csv',index=False);pd.DataFrame(records.values()).sort_values(['batch','plate','well']).to_csv(P/'frozen_wells.csv',index=False)
    frozen=['analysis_policy.json','frozen_fields.csv','frozen_contrasts.csv','frozen_wells.csv'];dependencies=['phase13/results/proposed_comparator_gene_sets.csv','phase13/results/proposed_target_wells.csv','phase12/frozen_policy.json','phase5/jump_raw_policy.json','phase5/harmonize.py','phase12/second_segmentation.py','phase11/spatial.py']
    plan={'frozen_at_UTC':datetime.now(timezone.utc).isoformat(),'fields':len(d),'wells':len(records),'genes':67,'dates':5,'target_fields':int(d.role.eq('target').sum()),'roles':d.groupby('role').field.size().to_dict(),'source_control_identities_verified':True,'hashes':{str(x):sha(P/x) for x in frozen},'dependencies':{x:sha(ROOT/x) for x in dependencies},'new_spatial_outcomes_known_at_freeze':False};(P/'frozen_plan.json').write_text(json.dumps(plan,indent=2));print(json.dumps(plan,indent=2))
if __name__=='__main__':main()
