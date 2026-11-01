from pathlib import Path
import requests,json,time,sys,re,html
sys.stdout.reconfigure(encoding='utf-8')
P=Path(__file__).resolve().parent;D=P/'data'
DOIS={2:'10.1074/jbc.M111.306019',3:'10.1093/nar/gkz263',4:'10.1038/s41586-018-0670-5',5:'10.1038/s41594-023-01074-9',7:'10.1038/s41556-024-01427-6',8:'10.1038/s41589-024-01584-7',9:'10.1038/s41467-025-58488-w',13:'10.1038/s41586-025-08878-3',15:'10.15252/embr.201847495',16:'10.1038/s41467-026-73609-9'}
def main():
    refs={}
    for number,doi in DOIS.items():
        url='https://api.crossref.org/works/'+doi;r=requests.get(url,timeout=45)
        if not r.ok:print(number,r.status_code,flush=True);continue
        a=r.json()['message'];authors=a.get('author',[]);names=[x.get('family','')+' '+''.join(w[0] for w in x.get('given','').split() if w) for x in authors[:6]];author=', '.join(names)+(' et al.' if len(authors)>6 else '.');year=a.get('published',a.get('issued',{})).get('date-parts',[[None]])[0][0];title=a.get('title',[''])[0];journal=a.get('container-title',[''])[0];loc=a.get('volume','')+(('('+a['issue']+')') if a.get('issue') else '')+(':'+a['page'] if a.get('page') else '');refs[str(number)]=f'{number}. {author} {title} {journal}. {year};{loc}. https://doi.org/{doi}';print(number,title[:100],flush=True);time.sleep(.15)
    refs={k:html.unescape(re.sub('<[^>]+>','',v)) for k,v in refs.items()}
    (D/'crossref_references.json').write_text(json.dumps(refs,indent=2),encoding='utf-8')
if __name__=='__main__':main()
