"""Retain complete processed count tables for all four selected U2OS screens."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import sys,json,re
sys.stdout.reconfigure(encoding='utf-8');P=Path(__file__).resolve().parent;D=P/'data';sys.path.insert(0,str(P));from search_portals import get
def main():
    tasks=[];log=[]
    for acc in ['GSE290538','GSE283147','GSE226994','GSE233195']:
        text=(D/(acc+'_metadata.txt')).read_text(encoding='utf-8')
        for line in text.splitlines():
            if not line.startswith('!Series_supplementary_file = '):continue
            url=line.split(' = ',1)[1].replace('ftp://','https://');name=url.rsplit('/',1)[-1]
            if 'HEK293T' not in name:tasks.append((name,url,None))
    with ThreadPoolExecutor(max_workers=3) as pool:
        for f in as_completed([pool.submit(get,*t) for t in tasks]):r=f.result();log.append(r);print(r['file'],r.get('status',r.get('error')),r.get('bytes'),flush=True)
    (P/'screen_count_provenance.json').write_text(json.dumps(log,indent=2))
if __name__=='__main__':main()
