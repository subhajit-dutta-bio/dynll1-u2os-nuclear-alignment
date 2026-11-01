"""Trace imports and actual plot inputs without repeating the model fits."""
from pathlib import Path
import os,sys,runpy,json,hashlib,importlib.util,subprocess
P=Path(__file__).resolve().parent;ROOT=P.parent;seen=set()
def audit(event,args):
    if event!='open' or not isinstance(args[0],(str,bytes,os.PathLike)):return
    if isinstance(args[1],str) and not args[1].startswith('r'):return
    try:
        f=Path(os.fsdecode(args[0])).resolve()
        if f.suffix=='.pyc':f=Path(importlib.util.source_from_cache(str(f)))
        rel=f.relative_to(ROOT)
        if f.is_file():seen.add(rel.as_posix())
    except (ValueError,OSError):pass
def main():
    sys.path.insert(0,str(P));sys.addaudithook(audit)
    if '--imports-only' in sys.argv:
        spec=importlib.util.spec_from_file_location('pixel_import_dependencies',P/'extract.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);module.init_worker();record=[{'path':n,'bytes':(ROOT/n).stat().st_size,'sha256':hashlib.sha256((ROOT/n).read_bytes()).hexdigest()} for n in sorted(seen)];(P/'qa/observed_image_import_reads.json').write_text(json.dumps(record,indent=2));return
    subprocess.run([sys.executable,str(Path(__file__)),'--imports-only'],cwd=ROOT,check=True)
    seen.update(row['path'] for row in json.loads((P/'qa/observed_image_import_reads.json').read_text()))
    for name,args in [('audit_coordinates.py',[]),('plot_all.py',[])]:
        sys.argv=[name,*args];runpy.run_path(str(P/name),run_name='__main__')
    record=[{'path':n,'bytes':(ROOT/n).stat().st_size,'sha256':hashlib.sha256((ROOT/n).read_bytes()).hexdigest()} for n in sorted(seen)];(P/'qa/observed_local_reads.json').write_text(json.dumps(record,indent=2));(P/'qa/observed_read_scope.json').write_text(json.dumps({'scope':'Pixel-pipeline imports, coordinate provenance and all plot inputs. Complete numerical measured inputs and document dependencies are also explicitly included by package.py. Fresh archive validation refits the models once with network access blocked.'},indent=2));print('Recorded actual project dependencies:',len(record))
if __name__=='__main__':main()
