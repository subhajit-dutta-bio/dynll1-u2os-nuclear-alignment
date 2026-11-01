"""Record public portal queries and original responses; inspect metadata first."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import requests,json,hashlib,time,gzip
P=Path(__file__).resolve().parent;D=P/'data';D.mkdir(parents=True,exist_ok=True)
def get(name,url,params=None):
    out=D/name
    if out.exists():return {'file':name,'url':url,'params':params,'status':'cached','bytes':out.stat().st_size,'sha256':hashlib.sha256(out.read_bytes()).hexdigest()}
    try:
        r=requests.get(url,params=params,timeout=90);r.raise_for_status();out.write_bytes(r.content);return {'file':name,'url':r.url,'params':params,'status':r.status_code,'bytes':len(r.content),'sha256':hashlib.sha256(r.content).hexdigest()}
    except Exception as e:return {'file':name,'url':url,'params':params,'error':repr(e)}
def main():
    log=[]
    terms={
      'geo_DYNLL1':'(DYNLL1 OR DNCL1 OR DNCLC1) AND gse[Entry Type]',
      'geo_U2OS_DYNLL1':'(U2OS OR U-2OS OR "U-2 OS") AND (DYNLL1 OR DNCL1 OR LC8 OR ASCIZ OR ATMIN) AND gse[Entry Type]',
      'geo_U2OS_singlecell':'(U2OS OR U-2OS OR "U-2 OS") AND ("perturb-seq" OR "CROP-seq" OR "single-cell CRISPR" OR "single cell perturbation") AND gse[Entry Type]',
      'geo_U2OS_CRISPR':'(U2OS OR U-2OS OR "U-2 OS") AND (CRISPR OR Cas9) AND gse[Entry Type]',
      'geo_U2OS_RNAi':'(U2OS OR U-2OS OR "U-2 OS") AND ("DYNLL1" OR "dynein light chain" OR "LC8" OR "FAM83D") AND gse[Entry Type]',
      'sra_U2OS_DYNLL1':'(U2OS OR U-2OS) AND (DYNLL1 OR DNCL1 OR LC8)',
    }
    for name,term in terms.items():
        result=get(name+'.json','https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi',{'db':'sra' if name.startswith('sra') else 'gds','term':term,'retmode':'json','retmax':300});log.append(result)
        if 'error' not in result:
            data=json.loads((D/(name+'.json')).read_text());q=data.get('esearchresult',{});print(name,'count',q.get('count'),'errors',q.get('errorlist'),flush=True);ids=q.get('idlist',[])
            if ids:
                time.sleep(.4);log.append(get(name+'_summaries.json','https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi',{'db':'sra' if name.startswith('sra') else 'gds','id':','.join(ids),'retmode':'json'}))
        time.sleep(.45)
    tasks=[]
    for word in ['DYNLL1','LC8','U2OS']:
        tasks.append((f'pride_{word}.json','https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects',{'keyword':word,'pageSize':100,'page':0}))
        tasks.append((f'biostudies_{word}.json','https://www.ebi.ac.uk/biostudies/api/v1/search',{'query':word,'pageSize':100}))
    for word in ['DYNLL1','U2OS']:
        tasks.append((f'ena_{word}.tsv','https://www.ebi.ac.uk/ena/portal/api/search',{'result':'read_study','query':f'(study_title="*{word}*" OR study_description="*{word}*")','fields':'study_accession,secondary_study_accession,study_title,study_description','format':'tsv','limit':1000}))
    tasks.extend([
        ('GSE92742_Broad_LINCS_pert_info.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_pert_info.txt.gz',None),
        ('GSE92742_Broad_LINCS_sig_info.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_sig_info.txt.gz',None),
        ('GSE92742_Broad_LINCS_inst_info.txt.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_inst_info.txt.gz',None),
    ])
    with ThreadPoolExecutor(max_workers=4) as pool:
        fs=[pool.submit(get,*x) for x in tasks]
        for f in as_completed(fs):r=f.result();log.append(r);print(r['file'],r.get('status',r.get('error')),r.get('bytes'),flush=True)
    (P/'retrieval_provenance.json').write_text(json.dumps({'retrieved_UTC':datetime.now(timezone.utc).isoformat(),'requests':log},indent=2))
if __name__=='__main__':main()
