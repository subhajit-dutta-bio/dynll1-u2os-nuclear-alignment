"""Recompute all ten CAMKK2 comparisons on the same physical plates."""
from pathlib import Path
import pandas as pd
P=Path(__file__).resolve().parent;ROOT=P.parent
def main():
    q=pd.read_csv(ROOT/'phase14/results/alignment_contrasts.csv');q=q[q.model.eq('adaptive_conditional')];m=pd.read_csv(ROOT/'phase14/frozen_contrasts.csv');q=q.merge(m[['gene','batch','target_plate','target_well']],on=['gene','batch'],validate='one_to_one');rows=[]
    for gene in ['HSPB8','CIT']:
        z=q[q.gene.eq('CAMKK2')].merge(q[q.gene.eq(gene)],on='batch',suffixes=('_CAMKK2','_comparison'),validate='one_to_one');assert len(z)==5 and z.target_plate_CAMKK2.eq(z.target_plate_comparison).all();z['comparison_gene']=gene;z['CAMKK2_minus_comparison']=z.effect_CAMKK2-z.effect_comparison;rows.append(z)
    out=pd.concat(rows);out.to_csv(P/'results/same_plate_CAMKK2_comparisons.csv',index=False);print('Verified and exported all ten same-physical-plate comparisons.')
if __name__=='__main__':main()
