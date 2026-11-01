"""Audit published LC8 reagent identity against both human paralog loci."""
from pathlib import Path
import requests,json,hashlib
P=Path(__file__).resolve().parent;D=P/'data/reagent_identity';D.mkdir(parents=True,exist_ok=True)
def main():
    records=[]
    for gene in ['DYNLL1','DYNLL2']:
        url=f'https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}?content-type=application/json';r=requests.get(url,timeout=60);r.raise_for_status();a=r.json();(D/f'{gene}_annotation.json').write_text(json.dumps(a,indent=2));start=max(0,a['start']-2001);end=a['end']+2000
        sequrl=f'https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr{a["seq_region_name"]};start={start};end={end}';r=requests.get(sequrl,timeout=60);r.raise_for_status();sequence=r.json()['dna'].upper();(D/f'{gene}_sequence.txt').write_text(sequence);records.append({'gene':gene,'ensembl_id':a['id'],'chromosome':a['seq_region_name'],'start_0based':start,'end_0based_exclusive':end,'annotation_url':url,'sequence_url':sequrl,'sequence_sha256':hashlib.sha256(sequence.encode()).hexdigest()})
    (D/'provenance.json').write_text(json.dumps(records,indent=2))
if __name__=='__main__':main()
