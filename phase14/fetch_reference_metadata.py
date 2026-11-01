"""Archive publisher-deposited bibliographic metadata for the revised references."""
from pathlib import Path
import concurrent.futures,json,re,requests,time,sys
from references import references
P=Path(__file__).resolve().parent
def one(line):
    number=int(line.split('.')[0]);doi=re.search(r'https://doi.org/(\S+)',line).group(1).rstrip('.');url='https://api.crossref.org/works/'+doi
    r=requests.get(url,timeout=50,headers={'User-Agent':'U2OS-public-data-bibliography/1.0'});r.raise_for_status();m=r.json()['message'];keep=['DOI','title','container-title','author','issued','published','volume','issue','page','article-number','URL','type'];return str(number),{'url':url,'retrieved_UTC':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'metadata':{k:m[k] for k in keep if k in m}}
def main():
    sys.stdout.reconfigure(encoding='utf-8');destination=P/'data/reference_metadata.json';out=json.loads(destination.read_text())['records'] if destination.exists() else {};refs=[s for s in references() if re.search(r'https://doi.org/',s) and int(s.split('.')[0]) not in [12,14] and s.split('.')[0] not in out];errors=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
        for future in concurrent.futures.as_completed([ex.submit(one,s) for s in refs]):
            try:k,v=future.result();out[k]=v;print(k,v['metadata'].get('title'),v['metadata'].get('page'),v['metadata'].get('article-number'),flush=True)
            except Exception as e:errors.append(str(e));print(type(e).__name__,str(e),flush=True)
    (P/'data/reference_metadata.json').write_text(json.dumps({'records':out,'errors':errors},indent=2),encoding='utf-8')
if __name__=='__main__':main()
