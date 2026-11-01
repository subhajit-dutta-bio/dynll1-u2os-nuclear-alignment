"""Resolve ENA query fields and retrieve targeted primary supplements."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import sys,json,requests,re,time
sys.stdout.reconfigure(encoding='utf-8');P=Path(__file__).resolve().parent;D=P/'data';sys.path.insert(0,str(P));from search_portals import get
def main():
    tasks=[];log=[]
    for word in ['DYNLL1','DNCL1','LC8','U2OS']:
        tasks.append((f'ena_{word}_corrected.tsv','https://www.ebi.ac.uk/ena/portal/api/search',{'result':'read_study','query':f'(study_title="*{word}*" OR description="*{word}*" OR sample_title="*{word}*")','fields':'study_accession,secondary_study_accession,study_title,description','format':'tsv','limit':10000}))
    for acc in ['S-SCDT-EMBOR-2018-47495-T','S-SCDT-10_15252-EMBJ_2023114473']:
        tasks.append((acc+'.json',f'https://www.ebi.ac.uk/biostudies/api/v1/studies/{acc}',None))
    for acc in ['PXD009335','PXD052362']:
        tasks.append((acc+'_files_page1.json',f'https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{acc}/files',{'pageSize':100,'page':1}))
    for acc in ['PXD039627','PXD039701','PXD076520']:
        tasks.append((acc+'_project.json',f'https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{acc}',None))
        tasks.append((acc+'_files.json',f'https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{acc}/files',{'pageSize':100,'page':0}))
    for pmc in ['PMC3270970','PMC6726907','PMC7155769','PMC6025793','PMC11879137','PMC9023290','PMC12540834','PMC6614850','PMC12648486','PMC10602929']:
        tasks.append((pmc+'_article.html',f'https://pmc.ncbi.nlm.nih.gov/articles/{pmc}/',None))
    for article in ['s41467-026-73609-9','s41467-025-64066-x','s41586-025-08878-3']:
        tasks.append((article+'.html',f'https://www.nature.com/articles/{article}',None))
    with ThreadPoolExecutor(max_workers=4) as pool:
        for f in as_completed([pool.submit(get,*task) for task in tasks]):r=f.result();log.append(r);print(r['file'],r.get('status',r.get('error')),r.get('bytes'),flush=True)
    (P/'deep_retrieval_provenance.json').write_text(json.dumps(log,indent=2))
if __name__=='__main__':main()
