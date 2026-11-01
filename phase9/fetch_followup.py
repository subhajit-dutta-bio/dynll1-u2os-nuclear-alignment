"""Paginated portal retrieval and primary-source data access."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import sys,json,time,requests
sys.stdout.reconfigure(encoding='utf-8')
P=Path(__file__).resolve().parent;D=P/'data';sys.path.insert(0,str(P));from search_portals import get
def main():
    log=[];tasks=[
      ('ena_study_fields.tsv','https://www.ebi.ac.uk/ena/portal/api/searchFields',{'result':'read_study'}),
      ('ena_study_return_fields.tsv','https://www.ebi.ac.uk/ena/portal/api/returnFields',{'result':'read_study'}),
      ('europepmc_DYNLL1_U2OS.json','https://www.ebi.ac.uk/europepmc/webservices/rest/search',{'query':'(DYNLL1 OR DNCL1 OR LC8) AND (U2OS OR "U-2OS" OR "U-2 OS")','format':'json','pageSize':100,'resultType':'core'}),
      ('biostudies_exact_DYNLL1_U2OS.json','https://www.ebi.ac.uk/biostudies/api/v1/search',{'query':'"DYNLL1" AND ("U2OS" OR "U-2OS")','pageSize':100,'expandSynonyms':'false','expandEfo':'false'}),
      ('biostudies_arrayexpress_DYNLL1.json','https://www.ebi.ac.uk/biostudies/api/v1/search',{'query':'"DYNLL1"','collection':'arrayexpress','pageSize':100,'expandSynonyms':'false','expandEfo':'false'}),
      ('U2OS_APMS.cx','https://www.ndexbio.org/v2/network/95bc75d5-d1d1-11ee-8a40-005056ae23aa',None),
      ('U2OS_cellmap.cx','https://www.ndexbio.org/v2/network/f693137a-d2d7-11ef-8e41-005056ae3c32',None),
    ]
    for acc in ['PXD009335','PXD052362']:
        tasks.append((acc+'_project.json',f'https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{acc}',None));tasks.append((acc+'_files.json',f'https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{acc}/files',{'pageSize':100,'page':0}))
    for acc in ['GSE284207','GSE314191','GSE232592','GSE218249','GSE183059','GSE30417']:
        tasks.append((acc+'_metadata.txt','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi',{'acc':acc,'targ':'self','form':'text','view':'full'}))
    for pmc in ['PMC3270970','PMC6726907','PMC7155769','PMC12803952']:
        tasks.append((pmc+'.xml',f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmc}/fullTextXML',None))
    with ThreadPoolExecutor(max_workers=4) as pool:
        for f in as_completed([pool.submit(get,*task) for task in tasks]):r=f.result();log.append(r);print(r['file'],r.get('status',r.get('error')),r.get('bytes'),flush=True)
    # Enumerate all pages for the literal U2OS keyword, retaining the first pass.
    for page in range(1,15):
        r=get(f'pride_U2OS_page{page}.json','https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects',{'keyword':'U2OS','pageSize':100,'page':page});log.append(r)
        if 'error' in r:break
        a=json.loads((D/r['file']).read_text(encoding='utf-8'));print('PRIDE page',page,'records',len(a),flush=True)
        if len(a)<100:break
    (P/'followup_provenance.json').write_text(json.dumps(log,indent=2))
if __name__=='__main__':main()
