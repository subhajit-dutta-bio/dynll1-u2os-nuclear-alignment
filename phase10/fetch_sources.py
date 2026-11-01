from pathlib import Path
import requests,json,time,hashlib
P=Path(__file__).resolve().parent;D=P/'data';D.mkdir(exist_ok=True)
def main():
    jobs=[(f'{acc}_samples.txt',f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=gsm&form=text&view=full') for acc in ['GSE226994','GSE290538','GSE283147','GSE233195']]
    jobs += [('CM4AI_PPIDownloader.zip','https://dataverse.lib.virginia.edu/api/access/datafile/120795'),('CM4AI_release_metadata.zip','https://dataverse.lib.virginia.edu/api/access/datafile/120796')]
    jobs += [('Brunello_official_library.txt','https://www.addgene.org/static/cms/filer_public/8b/4c/8b4c89d9-eac1-44b2-bb2f-8fea95672705/broadgpp-brunello-library-contents.txt'),('hTERT_primary.html','https://www.nature.com/articles/s41556-024-01427-6'),('hTERT_correction.html','https://www.nature.com/articles/s41556-026-01884-1'),('Honeywell_primary.html','https://www.nature.com/articles/s41589-024-01584-7'),('HIF_primary.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11977238/fullTextXML')]
    for pmid in ['38806647','40195317','39713309','41577826']:
        jobs.append((f'pubmed_{pmid}.json',f'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:{pmid}%20AND%20SRC:MED&format=json&resultType=core'))
    records=[]
    for name,url in jobs:
        dest=D/name
        if dest.exists(): body=dest.read_bytes();status='retained'
        else:
            r=requests.get(url,timeout=90);body=r.content;status=r.status_code
            if r.ok:dest.write_bytes(body)
        records.append({'file':name,'url':url,'status':status,'bytes':len(body),'sha256':hashlib.sha256(body).hexdigest()});print(name,status,len(body),flush=True);time.sleep(.4)
    (P/'download_provenance.json').write_text(json.dumps(records,indent=2))
if __name__=='__main__':main()
