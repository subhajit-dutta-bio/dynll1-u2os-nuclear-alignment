"""Build distinct reading and reproducibility archives without altering history."""
from pathlib import Path
import hashlib,json,zipfile,sys,importlib.metadata as md
P=Path(__file__).resolve().parent;ROOT=P.parent;DEST=ROOT.parent
DOCS=['DYNLL1_U2OS_main_manuscript_expanded.docx','DYNLL1_U2OS_supplementary_information.docx']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def allowed(p):return '__pycache__' not in p.parts and not any(x.startswith('render_') or x.startswith('renderer_') or x.startswith('pdf_') for x in p.parts) and p.suffix.lower() not in ['.log']
def main():
    ap=sys.argv[1:];verify=json.loads((P/'qa/verification.json').read_text());assert verify['all_pass']
    light={ROOT/x for x in DOCS}|{P/'README.txt',P/'panel_registry.csv',P/'panel_registry.json'}
    for folder in ['figures','panel_data','tables']:light.update(p for p in (P/folder).rglob('*') if p.is_file())
    full=set(light)
    prior=DEST/'DYNLL1_U2OS_manuscript_and_analysis.zip'
    with zipfile.ZipFile(prior) as z:
        for name in z.namelist():
            rel=Path(*name.split('/')[1:])
            if rel.suffix.lower()=='.docx':continue
            f=ROOT/rel
            if f.is_file():full.add(f)
    for x in json.loads((P/'qa/observed_local_reads.json').read_text()):
        f=ROOT/x['path']
        if f.is_file() and allowed(f):full.add(f)
    full.update(p for p in P.rglob('*') if p.is_file() and allowed(p))
    for folder in ['phase5','phase6','phase7','phase8','phase9']:
        full.update(p for p in (ROOT/folder).glob('*.json') if p.is_file())
    full.update(p for p in (ROOT/'phase6/qa').glob('*.json') if p.is_file())
    packages={name:md.version(name) for name in ['numpy','pandas','scipy','scikit-learn','statsmodels','matplotlib','networkx','tifffile','scikit-image','pyarrow','openpyxl','requests']}
    (P/'requirements_analysis.txt').write_text('\n'.join(f'{k}=={v}' for k,v in packages.items())+'\n',encoding='utf-8');full.add(P/'requirements_analysis.txt')
    metadata=[]
    for label,files in [('reading',light),('analysis',full)]:
        archive=DEST/('DYNLL1_U2OS_expanded_manuscript_figures.zip' if label=='reading' else 'DYNLL1_U2OS_expanded_analysis.zip')
        files={p for p in files if p.name not in ['package_inventory.json','package_input_inventory.json','standalone_validation.json']}
        inventory=[{'path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(files)]
        if label=='analysis':(P/'qa/package_input_inventory.json').write_text(json.dumps(inventory,indent=2),encoding='utf-8')
        with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
            for p in sorted(files):z.write(p,ROOT.name+'/'+p.relative_to(ROOT).as_posix())
            z.writestr(ROOT.name+'/DELIVERY_FILE_MANIFEST.json',json.dumps(inventory,indent=2))
            z.writestr(ROOT.name+'/README_expanded_manuscript.txt',(P/'README.txt').read_bytes())
            proof=P/'qa/standalone_validation.json'
            if proof.exists():z.write(proof,ROOT.name+'/phase11/qa/standalone_validation.json')
        with zipfile.ZipFile(archive) as z:assert z.testzip() is None
        metadata.append({'kind':label,'archive':str(archive),'files':len(files),'bytes':archive.stat().st_size,'sha256':sha(archive)})
    (P/'qa/package_inventory.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8');print(json.dumps(metadata,indent=2))
if __name__=='__main__':main()
