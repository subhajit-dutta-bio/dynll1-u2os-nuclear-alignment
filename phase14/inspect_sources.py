from pathlib import Path
import pandas as pd,json
P=Path(__file__).resolve().parent;R=P.parent
for name in ['phase2/results/screen_results.csv','phase2/results/frozen_candidates.csv','phase2/results/heldout_validation.csv','phase12/results/orf/all_director_yield.csv','phase12/results/orf/well_coverage.csv']:
 p=R/name
 if p.exists():
  d=pd.read_csv(p);print(name,len(d),d.columns.tolist())
print('phase2 results', [p.name for p in (R/'phase2/results').glob('*.csv')])
print('benchmark dirs',[p.name for p in (R/'phase12/results/cellpose_benchmark').iterdir() if p.is_dir()])
for p in (R/'phase12/results/cellpose_benchmark').rglob('*.parquet'):
 if 'mask' not in str(p):
  d=pd.read_parquet(p);print(p.relative_to(R),len(d),d.columns.tolist());break
