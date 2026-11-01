"""Finish targeted genome-wide screen and primary article retrieval."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import sys,json,re
sys.stdout.reconfigure(encoding='utf-8');P=Path(__file__).resolve().parent;D=P/'data';sys.path.insert(0,str(P));from search_portals import get
def main():
    tasks=[];log=[]
    for acc in ['GSE290538','GSE283147','GSE226994','GSE233195']:
        tasks.append((acc+'_metadata.txt','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi',{'acc':acc,'targ':'self','form':'text','view':'full'}))
    tasks.extend([
      ('West2019_primary.html','https://oup.silverchair-cdn.com/article-minimal/5458119',None),
      ('PXD013808_project.json','https://www.ebi.ac.uk/pride/ws/archive/v2/projects/PXD013808',None),
      ('PXD013808_files.json','https://www.ebi.ac.uk/pride/ws/archive/v2/projects/PXD013808/files',{'pageSize':100,'page':0}),
      ('Jurado2012_supplement.pdf','https://pmc.ncbi.nlm.nih.gov/articles/instance/3270970/bin/supp_M111.306019_jbc.M111.306019-1.pdf',None),
      ('GSM7356374_metadata.txt','https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi',{'acc':'GSM7356374','targ':'self','form':'text','view':'full'}),
      ('CM4AI_U2OS_dataverse.json','https://dataverse.lib.virginia.edu/api/datasets/:persistentId/',{'persistentId':'doi:10.18130/V3/164SI1'}),
    ])
    with ThreadPoolExecutor(max_workers=3) as pool:
        for f in as_completed([pool.submit(get,*t) for t in tasks]):r=f.result();log.append(r);print(r['file'],r.get('status',r.get('error')),r.get('bytes'),flush=True)
    for acc in ['PXD009335','PXD013808']:
        for page in range(2 if acc=='PXD009335' else 1,8):
            r=get(f'{acc}_files_page{page}.json',f'https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{acc}/files',{'pageSize':100,'page':page});log.append(r)
            if 'error' in r:break
            a=json.loads((D/r['file']).read_text(encoding='utf-8'));print(acc,'files page',page,len(a),flush=True)
            if len(a)<100:break
    (P/'remaining_retrieval_provenance.json').write_text(json.dumps(log,indent=2))
if __name__=='__main__':main()
