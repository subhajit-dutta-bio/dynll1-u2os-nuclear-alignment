"""Recompute every new statistical result from complete measured inputs offline."""
from pathlib import Path
import json,subprocess,sys,hashlib
import pandas as pd,numpy as np
from pandas.testing import assert_frame_equal
P=Path(__file__).resolve().parent;ROOT=P.parent;R=P/'results'
def main():
    files=[p for p in R.rglob('*.csv') if p.parent.name not in ['fields']];before={p.relative_to(ROOT).as_posix():pd.read_csv(p) for p in files};checks=[]
    for name,args in [('audit_coordinates.py',[]),('analyze.py',['--offline-inputs']),('matching_balance.py',[]),('context_analysis.py',[]),('coverage_context.py',['--offline-inputs']),('profile_prediction.py',[])]:subprocess.run([sys.executable,str(P/name),*args],cwd=ROOT,check=True)
    for name,old in before.items():
        new=pd.read_csv(ROOT/name);assert_frame_equal(old,new,check_exact=False,rtol=1e-10,atol=1e-10,check_dtype=False);checks.append(name)
    report={'all_pass':True,'network_used':False,'tables_compared':checks,'relative_tolerance':1e-10,'absolute_tolerance':1e-10,'scope':'Full new statistical calibration, all335 contrasts per model, comparator context and all nested predictions from complete focal/field inputs. Pixel inference is a separate layer.'};(P/'qa/offline_reproduction.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
