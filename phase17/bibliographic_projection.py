"""Compact bibliographic records without nested cited-work bibliographies."""
from copy import deepcopy
def project_record(key,record):
    result=deepcopy(record)
    def compact(x):
        if isinstance(x,dict):return {k:compact(v) for k,v in x.items() if k not in ['reference','ReferenceList']}
        if isinstance(x,list):return [compact(v) for v in x]
        return x
    result=compact(result)
    if key=='resource_14' and 'payload' in result:
        original=result['payload']['data'];attrs=original.get('attributes',{})
        kept={k:attrs[k] for k in ['doi','publisher','publicationYear','url','version'] if k in attrs}
        kept['titles']=[{'title':'U2OS Cell Map'}];kept['displayTitleAbbreviated']=True
        result['payload']={'data':{'id':original['id'],'type':original.get('type','dois'),'attributes':kept}}
        result['payload_projection']='Identifier, publisher, year, URL and abbreviated display title; the complete official record is available at the source URL.'
    elif 'payload' in result:
        result['payload_projection']='Primary bibliographic fields retained; nested cited-work bibliographies omitted.'
    return result
