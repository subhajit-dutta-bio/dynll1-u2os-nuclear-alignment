"""Download small published tables and verify direct-perturbation metadata."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
from bs4 import BeautifulSoup
import sys,json,re,requests,time
sys.stdout.reconfigure(encoding='utf-8');P=Path(__file__).resolve().parent;D=P/'data';sys.path.insert(0,str(P));from search_portals import get
def main():
    tasks=[];log=[]
    for filename in ['GSE232592_CRISPR_p53_GeneLevel.xlsx','GSE232592_CRISPR_p53_counts.txt.gz']:
        tasks.append((filename,'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232592/suppl/'+filename,None))
    for article,keep in [('s41467-026-73609-9',lambda text:'Supplementary Data' in text and any('Data '+str(i)+' ' in text for i in [1,2,3,4]) or 'Description of Additional' in text or 'Source Data' in text),('s41467-025-64066-x',lambda text:'Supplementary Data' in text),('s41586-025-08878-3',lambda text:'Supplementary Tables' in text)]:
        soup=BeautifulSoup((D/(article+'.html')).read_text(encoding='utf-8'),'html.parser')
        for a in soup.select('a[href]'):
            text=a.get_text(' ',strip=True);url=a['href']
            if keep(text) and re.search(r'\.(xlsx|pdf)$',url):tasks.append((url.rsplit('/',1)[-1],url,None))
    for pmc in ['3270970','6726907','6025793','11879137','6614850']:
        tasks.append(('PMC'+pmc+'_efetch.xml','https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',{'db':'pmc','id':pmc}))
    tasks.append(('ena_U2OS_DYNLL1_expanded.tsv','https://www.ebi.ac.uk/ena/portal/api/search',{'result':'read_run','query':'tax_eq(9606) AND (cell_line="U2OS" OR cell_line="U-2 OS" OR cell_line="U-2OS" OR sample_title="*U2OS*") AND (description="*DYNLL1*" OR sample_description="*DYNLL1*" OR sample_title="*DYNLL1*" OR experiment_title="*DYNLL1*")','fields':'study_accession,run_accession,sample_accession,study_title,sample_title,sample_description,cell_line,library_strategy','format':'tsv','limit':10000}))
    with ThreadPoolExecutor(max_workers=3) as pool:
        for f in as_completed([pool.submit(get,*t) for t in tasks]):r=f.result();log.append(r);print(r['file'],r.get('status',r.get('error')),r.get('bytes'),flush=True)
    # Screen alternate spellings, since a literal U2OS search is not comprehensive.
    for word in ['U-2OS','U-2 OS','U2-OS']:
        for page in range(5):
            r=get(f'pride_alias_{word.replace(" ","_")}_{page}.json','https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects',{'keyword':word,'pageSize':100,'page':page});log.append(r)
            if 'error' in r:break
            a=json.loads((D/r['file']).read_text(encoding='utf-8'));print('PRIDE',word,page,len(a),flush=True)
            if len(a)<100:break
    (P/'quantitative_retrieval_provenance.json').write_text(json.dumps(log,indent=2))
if __name__=='__main__':main()
