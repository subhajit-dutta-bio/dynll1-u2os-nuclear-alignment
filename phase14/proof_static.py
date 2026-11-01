"""Internal layout proof of completed fixed methods and historical context."""
from pathlib import Path
import json
from build_documents import Draft,P,ROOT,TITLE
from manuscript_content import SUPP
from references import references
def main():
    registry=json.loads((P/'qa/static_panel_registry.json').read_text());d=Draft('Internal layout proof');d.h('Nuclear organization in U2OS Supplementary information',0);d.p(TITLE);d.p('Internal typography and historical context proof. New comparator results have not been inserted.',10);d.add_content(SUPP);d.h('Reference key',1)
    for ref in references():d.p(ref,10)
    for i in [3,4,5]:d.figure('supplementary_figure'+str(i),i,registry)
    print(d.save('phase14/qa/static_methods_proof.docx','qa/static_methods_proof.txt'))
if __name__=='__main__':main()
