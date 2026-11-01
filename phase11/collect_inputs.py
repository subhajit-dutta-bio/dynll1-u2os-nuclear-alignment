"""Trace local reads during offline analysis and plotting for portable packaging."""
from pathlib import Path
import os,sys,json,runpy,hashlib
P=Path(__file__).resolve().parent;ROOT=P.parent;seen=set()
def audit(event,args):
    if event!='open' or not isinstance(args[0],(str,bytes,os.PathLike)):return
    mode=args[1]
    if isinstance(mode,str) and not mode.startswith('r'):return
    try:
        path=Path(os.fsdecode(args[0])).resolve();rel=path.relative_to(ROOT)
        if path.is_file() and '__pycache__' not in path.parts:seen.add(rel.as_posix())
    except (ValueError,OSError):pass
def main():
    sys.path.insert(0,str(P));sys.addaudithook(audit)
    for name in ['verify.py','plot_main.py','plot_supp.py']:
        sys.argv=[name];runpy.run_path(str(P/name),run_name='__main__')
    import pandas as pd
    from spatial import path_for
    for field in pd.read_csv(ROOT/'phase6/frozen_fields.csv').field:seen.add(path_for(field).relative_to(ROOT).as_posix())
    paths=sorted(x for x in seen if (ROOT/x).exists())
    result=[{'path':x,'bytes':(ROOT/x).stat().st_size,'sha256':hashlib.sha256((ROOT/x).read_bytes()).hexdigest()} for x in paths]
    (P/'qa/observed_local_reads.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print('Retained local read dependencies',len(result))
if __name__=='__main__':main()
