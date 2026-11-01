from pathlib import Path
from bs4 import BeautifulSoup
import requests,json,hashlib
P=Path(__file__).resolve().parent;D=P/'data'
def main():
    soup=BeautifulSoup((D/'hTERT_primary.html').read_text(encoding='utf-8'),'html.parser');out=[]
    for x in soup.select('a[href]'):
        url=x['href']
        if 'media.springernature.com' not in url or not url.endswith('.xlsx'):continue
        name=url.rsplit('/',1)[-1];dest=D/name
        if not dest.exists():
            r=requests.get(url,timeout=120);r.raise_for_status();dest.write_bytes(r.content)
        out.append({'file':name,'url':url,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'bytes':dest.stat().st_size});print(name,dest.stat().st_size,flush=True)
    (P/'hTERT_table_provenance.json').write_text(json.dumps(out,indent=2))
if __name__=='__main__':main()
