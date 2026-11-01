from pathlib import Path
import requests,json,hashlib
P=Path(__file__).resolve().parent;D=P/'data';D.mkdir(exist_ok=True)
def main():
    url='https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt';dest=D/'hgnc_complete_set.txt'
    if not dest.exists():
        r=requests.get(url,timeout=120);r.raise_for_status();dest.write_bytes(r.content)
    out={'source':'HGNC complete approved gene set','url':url,'retrieval_date':'2026-10-30','file':dest.name,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'bytes':dest.stat().st_size};(P/'annotation_provenance.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
