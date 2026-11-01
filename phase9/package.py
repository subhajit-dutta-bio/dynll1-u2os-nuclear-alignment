"""Package the source search separately from all earlier evidence deliveries."""
from pathlib import Path
import json,hashlib,zipfile,platform,importlib.metadata
P=Path(__file__).resolve().parent;ROOT=P.parent;OUT=ROOT.parent
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def main():
    versions={x:importlib.metadata.version(x) for x in ['numpy','pandas','requests','openpyxl','pyarrow','beautifulsoup4','lxml']};(P/'runtime_versions.json').write_text(json.dumps({'Python':platform.python_version(),'packages':versions},indent=2))
    files=[]
    for f in P.rglob('*'):
        if not f.is_file() or '__pycache__' in f.parts or f.suffix=='.log':continue
        rel=f.relative_to(P)
        if 'qa' in rel.parts and f.suffix in ['.png','.pdf']:continue
        if f.name in ['checksums.json','package_manifest.json','archive_validation.json']:continue
        files.append(f)
    files +=[ROOT/'DYNLL1_U2OS_public_validation_sources.docx',ROOT/'phase7/results/DYNLL1_shared_reagent_identity.csv'];files=sorted(set(files));checks={f.relative_to(ROOT).as_posix():sha(f) for f in files};(P/'checksums.json').write_text(json.dumps(checks,indent=2));files.append(P/'checksums.json')
    groups=[('DYNLL1_U2OS_public_validation_summary.zip',[f for f in files if f.parent==P/'results' or f.name in ['DYNLL1_U2OS_public_validation_sources.docx','README.txt','verification.json','visual_review.json']]),('DYNLL1_U2OS_public_validation_data.zip',files)];result=[]
    for name,items in groups:
        dest=OUT/name
        with zipfile.ZipFile(dest,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
            for f in items:z.write(f,(Path(ROOT.name)/f.relative_to(ROOT)).as_posix(),compress_type=zipfile.ZIP_STORED if f.suffix in ['.xlsx','.zip','.gz','.docx','.parquet'] else zipfile.ZIP_DEFLATED)
        with zipfile.ZipFile(dest) as z:assert z.testzip() is None;count=len(z.namelist())
        result.append({'archive':str(dest),'bytes':dest.stat().st_size,'entries':count,'sha256':sha(dest),'CRC_pass':True})
    manifest={'archives':result,'checksummed_files':len(checks),'scope':'Source finding and published-score extraction; not new biological/spatial validation.','earlier_evidence_files_changed':False};(P/'package_manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest,indent=2))
if __name__=='__main__':main()
