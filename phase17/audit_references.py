"""Explicit, portable reference-by-reference audit using primary web records."""
from pathlib import Path
import argparse,concurrent.futures,hashlib,json,re,time
import pandas as pd
from references import references
P=Path(__file__).resolve().parent;ROOT=P.parent
RESOURCE={6:'https://bbbc.broadinstitute.org/BBBC039',12:'https://api.figshare.com/v2/articles/27993248',14:'https://api.datacite.org/dois/10.18130/v3/164si1',18:'https://www.genenames.org/download/archive/',19:'https://www.janelia.org/node/69901',25:'https://www.addgene.org/pooled-library/broadgpp-human-knockout-brunello/'}
def retrieve(item):
    import requests
    key,url=item
    try:
        r=requests.get(url,timeout=45,headers={'User-Agent':'U2OS-manuscript-bibliography/1.0'});r.raise_for_status();typ=r.headers.get('content-type','');payload=r.json() if 'json' in typ else r.text
        return key,{'url':url,'resolved_url':r.url,'retrieved_UTC':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'status':r.status_code,'payload':payload,'sha256':hashlib.sha256(r.content).hexdigest()}
    except Exception as e:return key,{'url':url,'error':str(e)}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--fetch',action='store_true');args=ap.parse_args();dest=P/'data/reference_web_records.json';refs=references();journal=[i for i in range(1,26) if i not in RESOURCE];pub=json.loads((ROOT/'phase14/data/pubmed_reference_audit.json').read_text());pub={r['reference']:r for r in pub};old=json.loads((ROOT/'phase14/data/reference_metadata.json').read_text())['records']
    if args.fetch:
        jobs=[('crossref_'+str(i),'https://api.crossref.org/works/'+re.search(r'https://doi.org/([^\s]+)',refs[i-1]).group(1).rstrip('.')) for i in journal];jobs += [('resource_'+str(i),u) for i,u in RESOURCE.items()]
        records=json.loads(dest.read_text(encoding='utf-8')) if dest.exists() else {}
        jobs=[j for j in jobs if 'payload' not in records.get(j[0],{})]
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:records.update(dict(pool.map(retrieve,jobs)))
        dest.write_text(json.dumps(records,indent=2),encoding='utf-8')
    records=json.loads(dest.read_text(encoding='utf-8'));rows=[]
    for i,line in enumerate(refs,1):
        if i in journal:
            doi=re.search(r'https://doi.org/([^\s]+)',line).group(1).rstrip('.');a=pub[i];assert a['exact_DOI_verified'];record=records['crossref_'+str(i)];fresh='payload' in record;meta=record['payload']['message'] if fresh else old[str(i)]['metadata'];assert meta['DOI'].lower()==doi.lower()
            rows.append({'reference':i,'kind':'Journal article','source_1':'PubMed exact DOI record','source_1_URL':'https://pubmed.ncbi.nlm.nih.gov/'+a['PMIDs'][0]+'/','source_2':'Publisher-deposited Crossref web metadata','source_2_URL':'https://api.crossref.org/works/'+doi,'checked_title':meta['title'][0],'status':'Two bibliographic records agree','current_web_fetch':fresh,'limitation':'Bibliographic identity only, not validation of every biological claim.'})
        else:
            record=records['resource_'+str(i)];ok='payload' in record
            assert ok,(i,record)
            title=' '.join(re.findall(r'<title[^>]*>(.*?)</title>',record['payload'],re.I|re.S)) if isinstance(record['payload'],str) else str(record['payload'].get('title',record['payload'].get('data',{}).get('attributes',{}).get('titles',[])))
            if i==12:assert record['payload']['doi'].lower().startswith('10.25452/figshare.plus.27993248')
            if i==14:assert record['payload']['data']['id'].lower()=='10.18130/v3/164si1'
            rows.append({'reference':i,'kind':'Preprint' if i==19 else 'Resource','source_1':'Official author page' if i==19 else 'Official resource or DOI registry','source_1_URL':RESOURCE[i],'source_2':'bioRxiv DOI and author-lab record' if i==19 else 'Archived resource files and retrieval provenance','source_2_URL':'https://www.biorxiv.org/content/10.1101/2025.04.28.651001v1' if i==19 else '', 'checked_title':re.sub(r'\s+',' ',title).strip(),'status':'Official record checked','current_web_fetch':True,'limitation':'A dataset or preprint is not represented as a PubMed-indexed journal article; resource provenance is separate from two-database journal verification.'})
    frame=pd.DataFrame(rows);frame.to_csv(P/'results/reference_checks_all_25.csv',index=False);report={'references':25,'journal_references_PubMed_and_Crossref':len(journal),'fresh_Crossref_fetches':int(frame.query('kind == "Journal article"').current_web_fetch.sum()),'official_resource_or_author_web_checks':len(RESOURCE),'errors':{k:v['error'] for k,v in records.items() if 'error' in v},'scope':'Explicit identity audit. Archived PubMed records are reused; new web retrievals supplement them. No journal fit, impact factor, DEAL eligibility or scientific validation is inferred.'};(P/'qa/reference_verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
