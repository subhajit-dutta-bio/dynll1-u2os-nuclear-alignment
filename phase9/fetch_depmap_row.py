"""Retain the U2OS row of two official archived DepMap matrices."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import requests,json,csv,hashlib,sys
sys.stdout.reconfigure(encoding='utf-8');P=Path(__file__).resolve().parent;D=P/'data'
def one(file):
    name=file['name'];destination=D/('DepMap24Q4_U2OS_'+name);url=file['download_url'];matched=False;lines=0
    with requests.get(url,stream=True,timeout=(30,120)) as r:
        r.raise_for_status();iterator=r.iter_lines(chunk_size=65536);header=next(iterator).decode('utf-8-sig');columns=next(csv.reader([header]));assert 'DYNLL1 (8655)' in columns
        for raw in iterator:
            lines+=1
            if raw.startswith(b'ACH-000364,') or raw.startswith(b'"ACH-000364",'):
                text=raw.decode('utf-8');destination.write_text(header+'\n'+text+'\n',encoding='utf-8');values=next(csv.reader([text]));matched=True;break
        assert matched,name
        values=dict(zip(columns,values));record={'source_file':name,'download_url':url,'source_declared_size':file['size'],'official_file_md5':file.get('computed_md5'),'download_scope':'Header and U2OS row retained from streamed matrix; request stopped after locating U2OS. Not a full-file download or checksum.','rows_scanned':lines,'model':'ACH-000364','DYNLL1':values['DYNLL1 (8655)'],'DYNLL2':values.get('DYNLL2 (140735)'),'retained_file':str(destination.relative_to(P)),'retained_sha256':hashlib.sha256(destination.read_bytes()).hexdigest()}
    return record
def main():
    j=json.loads((D/'DepMap24Q4_files.json').read_text(encoding='utf-8'));files=[x for x in j['files'] if x['name'] in ['CRISPRGeneEffect.csv','CRISPRGeneDependency.csv']]
    with ThreadPoolExecutor(max_workers=2) as pool:records=list(pool.map(one,files))
    (P/'depmap_matrix_provenance.json').write_text(json.dumps(records,indent=2));print(json.dumps(records,indent=2))
if __name__=='__main__':main()
