"""Check curated perturbation catalogs without downloading expression matrices."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
from bs4 import BeautifulSoup
import sys,json,re,requests
sys.stdout.reconfigure(encoding='utf-8');P=Path(__file__).resolve().parent;D=P/'data';sys.path.insert(0,str(P));from search_portals import get
def main():
    log=[];r=get('perturbseqr_download.html','https://perturbseqr.maayanlab.cloud/download');log.append(r);soup=BeautifulSoup((D/'perturbseqr_download.html').read_text(encoding='utf-8'),'html.parser');tasks=[]
    sections=['L1000-XPR (Gene Knockouts)','PerturbAtlas Human Gene Perturbations','CREEDS Gene Perturbations'];seen=set()
    for heading in soup.find_all(['h2','h3']):
        label=heading.get_text(' ',strip=True)
        if label not in sections:continue
        for sibling in heading.next_siblings:
            if getattr(sibling,'name',None) in ['h2','h3']:break
            for a in sibling.find_all('a',href=True) if hasattr(sibling,'find_all') else []:
                text=a.get_text(' ',strip=True);url=a['href']
                if text in ['Cell Line','Gene Targets'] or (label!='L1000-XPR (Gene Knockouts)' and text=='Gene Matrix Transposed'):
                    key=re.sub('[^a-zA-Z0-9]','_',label+'_'+text);suffix='.gmt.gz' if text=='Gene Matrix Transposed' else '.json';tasks.append((key+suffix,url,None));seen.add(url)
    tasks.extend([('scPerturb_explorer.html','https://projects.sanderlab.org/scperturb/explorer.html',None),('DepMap_figshare_search.json','https://api.figshare.com/v2/articles?resource_doi=10.25452%2Ffigshare.plus.25880521',None)])
    with ThreadPoolExecutor(max_workers=3) as pool:
        for f in as_completed([pool.submit(get,*t) for t in tasks]):r=f.result();log.append(r);print(r['file'],r.get('status',r.get('error')),r.get('bytes'),flush=True)
    # Public metadata search, not an account mutation.
    r=requests.post('https://api.figshare.com/v2/articles/search',json={'search_for':'"DepMap Public"','order':'published_date','order_direction':'desc','page_size':5},timeout=60);r.raise_for_status();(D/'DepMap_release_search.json').write_bytes(r.content);log.append({'file':'DepMap_release_search.json','url':r.url,'method':'POST read-only search','status':r.status_code});print('DepMap release metadata',len(r.json()),flush=True)
    (P/'catalog_retrieval_provenance.json').write_text(json.dumps(log,indent=2))
if __name__=='__main__':main()
