"""Second-source bibliography check by exact DOI in NCBI PubMed."""
from pathlib import Path
import requests,json,time,re,sys,xml.etree.ElementTree as ET
from references import references
P=Path(__file__).resolve().parent
def main():
    sys.stdout.reconfigure(encoding='utf-8');session=requests.Session();session.headers['User-Agent']='U2OS-public-data-bibliography/1.0';out=[]
    for line in references():
        n=int(line.split('.')[0]);match=re.search(r'https://doi.org/(\S+)',line)
        if not match or n in [12,14]:continue
        doi=match.group(1).rstrip('.');record={'reference':n,'doi':doi}
        try:
            response=session.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi',params={'db':'pubmed','term':doi+'[doi]','retmode':'json'},timeout=40);response.raise_for_status();record['query_url']=response.url;record['PMIDs']=response.json()['esearchresult']['idlist'];print(n,record['PMIDs'],flush=True)
        except Exception as e:record['error']=repr(e);print(n,repr(e),flush=True)
        out.append(record);time.sleep(.6)
    ids=[pmid for row in out for pmid in row.get('PMIDs',[])]
    if ids:
        response=session.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',params={'db':'pubmed','id':','.join(ids),'retmode':'xml'},timeout=60);response.raise_for_status();(P/'data/pubmed_reference_records.xml').write_bytes(response.content);records={}
        for article in ET.fromstring(response.content).findall('PubmedArticle'):
            pmid=article.findtext('MedlineCitation/PMID');title=''.join(article.find('MedlineCitation/Article/ArticleTitle').itertext());dois=[x.text for x in article.findall('PubmedData/ArticleIdList/ArticleId') if x.attrib.get('IdType')=='doi'];records[pmid]={'title':title,'doi':dois,'journal':article.findtext('MedlineCitation/Article/Journal/Title'),'volume':article.findtext('MedlineCitation/Article/Journal/JournalIssue/Volume'),'pages':article.findtext('MedlineCitation/Article/Pagination/StartPage') or article.findtext('MedlineCitation/Article/Pagination/MedlinePgn')}
        for row in out:
            row['records']=[records[x] for x in row.get('PMIDs',[]) if x in records];row['exact_DOI_verified']=any(row['doi'].lower() in [d.lower() for d in r['doi']] for r in row['records'])
    (P/'data/pubmed_reference_audit.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
if __name__=='__main__':main()
