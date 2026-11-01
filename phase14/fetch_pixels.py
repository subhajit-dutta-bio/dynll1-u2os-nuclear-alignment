"""Fetch only the frozen public U2OS DNA fields, with native TIFF validation."""
from pathlib import Path
import json,hashlib,time,argparse,shutil,threading
from concurrent.futures import ThreadPoolExecutor,as_completed
import requests,pandas as pd,numpy as np,tifffile
P=Path(__file__).resolve().parent;ROOT=P.parent;D=P/'data/images';Q=P/'data/provenance'
LOCAL=threading.local()
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
def one(row):
    if not hasattr(LOCAL,'session'):LOCAL.session=requests.Session()
    f=D/(row['field']+'.tif');record=Q/(row['field']+'.json')
    if f.exists() and record.exists():return 'existing'
    for source in [ROOT/'phase6/data/raw/images'/f.name,ROOT/'phase5/data/jump_raw/images'/f.name]:
        if source.exists():
            a=tifffile.imread(source);assert a.shape==(1080,1280);shutil.copy2(source,f);record.write_text(json.dumps({'field':row['field'],'url':row['dna_url'],'sha256':sha(f),'bytes':f.stat().st_size,'source':'identical existing public pixel cache','cache_path':str(source)},indent=2));return 'cache'
    for attempt in range(5):
        try:
            r=LOCAL.session.get(row['dna_url'],timeout=(15,90));r.raise_for_status();tmp=f.with_suffix('.partial.tif');tmp.write_bytes(r.content);a=tifffile.imread(tmp);assert a.shape==(1080,1280) and np.issubdtype(a.dtype,np.integer);tmp.replace(f);record.write_text(json.dumps({'field':row['field'],'url':row['dna_url'],'sha256':sha(f),'bytes':f.stat().st_size,'ETag':r.headers.get('ETag'),'source':'public Cell Painting Gallery S3'},indent=2));return 'download'
        except Exception:
            if attempt==4:raise
            time.sleep(min(2**attempt,8))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=16);args=ap.parse_args();D.mkdir(parents=True,exist_ok=True);Q.mkdir(exist_ok=True);(P/'qa').mkdir(exist_ok=True);plan=json.loads((P/'frozen_plan.json').read_text());assert sha(P/'frozen_fields.csv')==plan['hashes']['frozen_fields.csv'];d=pd.read_csv(P/'frozen_fields.csv').fillna('');errors=[];counts={};start=time.monotonic()
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        fs={ex.submit(one,row):row['field'] for row in d.to_dict('records')}
        for i,task in enumerate(as_completed(fs),1):
            try:k=task.result();counts[k]=counts.get(k,0)+1
            except Exception as e:errors.append({'field':fs[task],'error':repr(e)})
            if i%100==0:print('Fetched',i,len(d),counts,'errors',len(errors),'seconds',round(time.monotonic()-start),flush=True)
    (P/'qa/fetch_errors.json').write_text(json.dumps(errors,indent=2));assert not errors,errors[:5];pd.DataFrame([json.loads((Q/(f+'.json')).read_text()) for f in d.field]).to_csv(P/'data/pixel_provenance.csv',index=False);print('Complete',len(d),'native fields',flush=True)
if __name__=='__main__':main()
