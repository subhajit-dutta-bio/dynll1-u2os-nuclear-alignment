"""Write the U2OS public-validation source report with bundled document Python."""
from pathlib import Path
import json,csv
from docx import Document
from docx.shared import Inches,Pt,RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.opc.constants import RELATIONSHIP_TYPE as RT
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT
P=Path(__file__).resolve().parent;ROOT=P.parent;R=P/'results'
def main():
    scope=json.loads((R/'search_scope.json').read_text());s=json.loads((R/'summary.json').read_text());v=json.loads((P/'qa/verification.json').read_text());assert v['all_pass'];dm=s['DepMap'];net=s['APMS']
    D=Document();sec=D.sections[0];sec.page_width=Inches(8.5);sec.page_height=Inches(11);sec.top_margin=sec.bottom_margin=Inches(.65);sec.left_margin=sec.right_margin=Inches(.75)
    for name in ['Normal','Title','Heading 1','Heading 2']:
        st=D.styles[name];st.font.name='Calibri';st.font.color.rgb=RGBColor(0,0,0)
    D.styles['Normal'].font.size=Pt(11);D.styles['Normal'].paragraph_format.space_after=Pt(7);D.styles['Normal'].paragraph_format.line_spacing=1.04;D.styles['Title'].font.size=Pt(23);D.styles['Title'].paragraph_format.space_after=Pt(9)
    for name,size in [('Heading 1',15),('Heading 2',12)]:D.styles[name].font.size=Pt(size);D.styles[name].paragraph_format.space_after=Pt(6);D.styles[name].paragraph_format.space_before=Pt(14 if name=='Heading 1' else 8)
    footer=sec.footer.paragraphs[0];footer.alignment=2;footer.add_run('DYNLL1 public validation in U2OS   |   ').font.size=Pt(9);field=OxmlElement('w:fldSimple');field.set(qn('w:instr'),'PAGE');footer._p.append(field)
    def p(t):return D.add_paragraph(t)
    def h(t):return D.add_heading(t,1)
    def page(t):h(t).paragraph_format.page_break_before=True
    def link(para,label,url):
        el=OxmlElement('w:hyperlink');el.set(qn('r:id'),para.part.relate_to(url,RT.HYPERLINK,is_external=True));run=OxmlElement('w:r');pr=OxmlElement('w:rPr');color=OxmlElement('w:color');color.set(qn('w:val'),'0563C1');pr.append(color);run.append(pr);t=OxmlElement('w:t');t.text=label;run.append(t);el.append(run);para._p.append(el)
    def refs(items):
        q=p('Sources   ')
        for i,(label,url) in enumerate(items):
            if i:q.add_run('   |   ')
            link(q,label,url)
        for r in q.runs:r.font.size=Pt(9)
    def table(headers,rows,widths):
        tab=D.add_table(rows=1,cols=len(headers));tab.autofit=False
        for c,t in zip(tab.rows[0].cells,headers):c.text=t
        for row in rows:
            for c,t in zip(tab.add_row().cells,row):c.text=str(t)
        for k,row in enumerate(tab.rows):
            row._tr.get_or_add_trPr().append(OxmlElement('w:cantSplit'))
            for i,c in enumerate(row.cells):
                c.width=Inches(widths[i]);c.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER;pr=c._tc.get_or_add_tcPr();border=OxmlElement('w:tcBorders')
                for side in ['top','left','bottom','right']:
                    x=OxmlElement('w:'+side);x.set(qn('w:val'),'single');x.set(qn('w:sz'),'4');x.set(qn('w:color'),'D9D9D9');border.append(x)
                pr.append(border);mar=OxmlElement('w:tcMar')
                for side in ['top','left','bottom','right']:
                    x=OxmlElement('w:'+side);x.set(qn('w:w'),'85');x.set(qn('w:type'),'dxa');mar.append(x)
                pr.append(mar)
                if k==0:sh=OxmlElement('w:shd');sh.set(qn('w:fill'),'EAEAEA');pr.append(sh)
                for q in c.paragraphs:
                    q.paragraph_format.space_after=Pt(4);q.paragraph_format.line_spacing=1
                    for r in q.runs:r.font.size=Pt(10);r.bold=k==0
        p('').paragraph_format.space_after=Pt(1)
    D.add_heading('DYNLL1 validation sources in U2OS',0);p('Public repository and publication search   |   29 to 30 October 2026')
    p('Independent U2OS experiments already support molecular functions of DYNLL1 in nuclear protein organization. This search also located six public CRISPR source families containing DYNLL1 measurements and usable U2OS protein-interaction resources. No verified U2OS DYNLL1 Perturb-seq cohort or independent measurement of nuclear-neighbor alignment was obtained. The molecular evidence and the original spatial candidate must remain distinct.')
    h('Direct functional evidence in U2OS')
    table(['Study','Verified U2OS experiment','What it supports'],[
        ['Jurado 2012','DYNLL1 RNAi; GFP-ASCIZ reporter; MMS treatment; endogenous co-IP','Loss of damage-induced ASCIZ nuclear foci and an ASCIZ–DYNLL1 association'],
        ['West 2019','LC8 KO and 53BP1 foci in Figure 3A–B; DYNLL1 siRNA and checkpoint signaling in Figure 4A–B','Nuclear foci assembly and damage signaling; U2OS panels are identifiable'],
        ['He 2018','DYNLL1 and BRCA1 codepletion in AsiSI U2OS; resection qPCR; U2OS19 tethering','A role in nuclear DNA processing under defined damage conditions'],
        ['Swift 2023','U2OS AsiSI Figure 3B; EGFP versus FHA-DYNLL1 in 53BP1-depleted cells','Engineered recruitment affects the resection readout; n = 3 biological experiments']],[1.15,3.1,2.75])
    p('Jurado describes a SMARTpool plus four individual DYNLL1 siRNAs. We verified the published nuclear-foci result, but have not recovered the original full-field screen archive or newly quantified each siRNA. Its associated GEO array, GSE30417, is mouse fibroblasts and is excluded from the U2OS molecular cohort.')
    p('West provides both positive and negative U2OS findings. Reduced 53BP1 foci and altered checkpoint signaling do not imply a general repair defect: the paper also reports negative RPA/repair or replication findings. The published LC8 checkpoint rescue in HeLa is excluded. The printed LC8 siRNA exactly matches the retained DYNLL1 reference; the printed KO guide strings have one and two mismatches to that reference, so edited-allele identity remains a separate audit point.')
    h('Implication for the current project')
    p('There is a defensible U2OS molecular evidence base. It does not establish that these functions cause the existing neighbor-alignment association. Our frozen four-gene recruitment-module test failed its proposed shared positive spatial effect. A literature connection cannot replace that failed test or supply a rescue of the spatial endpoint.')
    refs([('Jurado 2012','https://doi.org/10.1074/jbc.M111.306019'),('West 2019','https://doi.org/10.1093/nar/gkz263'),('He 2018','https://doi.org/10.1038/s41586-018-0670-5'),('Swift 2023','https://doi.org/10.1038/s41594-023-01074-9')])
    page('Public CRISPR sources containing DYNLL1')
    p('The complete processed count tables were downloaded for the five GEO screens below. DYNLL1 rows were extracted directly from the original files. The DepMap analysis retains the official U2OS rows and complete guide and screen metadata from a fixed archived release. These resources measure different endpoints and should not be pooled into one validation statistic.')
    table(['Source','Guides and replicates','Readout','Current finding'],[
        ['DepMap 24Q4','3 guides used; 2 passing replicates; 21 days','Basal fitness','Gene effect −0.332; dependency probability 0.308'],
        ['GSE232592','5 retained guides; triplicate arms','Growth and etoposide death','No reported DYNLL1 FDR below 0.05'],
        ['GSE290538','4 guides; triplicate low/high/library arms','HIF activity reporter','Coverage verified; new inference pending'],
        ['GSE226994','4 guides per sheet; duplicate reference and engineered arms','hTERT-context fitness','Normalized counts retained; new inference pending'],
        ['GSE283147','1 released guide; duplicate arms','Pol II degradation live/dead response','Single-guide coverage; limited specificity'],
        ['GSE233195','4 guides; 2 biological replicates','Interferon-response selection','Coverage verified; lower mechanistic priority']],[1.35,1.85,1.75,2.05])
    h('Two directly available quantitative results')
    p('DepMap identifies U2OS as ACH-000364 and screen SC-000364.AV01. Both Avana replicate sequences pass the released QC. DYNLL1 has a Chronos gene-effect score of −0.331559 and a dependency probability of 0.307845. This supplies modest fitness context, not evidence for a strong dependency. The probability is not a p value, and the 21-day fitness assay cannot rule out proliferation confounding in the 96-hour imaging experiment. This is the archived 24Q4 release; the latest release was not verified.')
    p('GSE232592 reports U2OS WT growth-rate FDR 0.09290 and death-rate FDR 0.86917 for DYNLL1. All eight released FDR values across its rate and fold-change tables exceed 0.05. The p53KO arms remain separate. Its RNA-seq profiles conditioned-media exposure, not DYNLL1 knockout. The result is retained as a nonconfirming test of those endpoints.')
    h('Reagent overlap must be carried forward')
    p('One HIF guide and one Avana guide used by Chronos exactly equal JUMP construct CM-005281-02. Two further HIF targets overlap JUMP intervals by 12 and 19 bases; the remaining HIF target has no interval overlap. The sole released Pol II response guide overlaps a JUMP interval by 19 bases. Study independence therefore does not mean complete reagent independence. Target-by-target overlap is supplied in the source tables; no arbitrary independence threshold was imposed.')
    refs([('DepMap release','https://doi.org/10.25452/figshare.plus.27993248'),('GEO accession records','https://www.ncbi.nlm.nih.gov/geo/')])
    page('U2OS proteomics and protein assemblies')
    table(['Resource','Retained evidence','Interpretation'],[
        ['MSV000097168 and NDEx U2OS AP-MS','Full network and DYNLL1/DYNLL2 target edges','Protein associations; no DYNLL1 loss-of-function contrast or spatial readout'],
        ['PXD013808','Endogenous FAM83D-GFP U2OS affinity proteomics; 96-file inventory','DYNLL1 occurs in asynchronous and mitotic FAM83D IPs; FAM83D rescue is not DYNLL1 rescue'],
        ['PXD009335','FAM83-family affinity proteomics; 236-file inventory','Mixed U2OS and HEK293 data require sample-level separation'],
        ['PXD052362','U2OS native SEC-MS; 3 biological replicates; assembly tables','A biochemical assembly resource; gene-specific coelution not extracted here'],
        ['qChIP-MS 2026','U2OS ZBTB48 versus IgG and WT/KO control tables','A checked lead with mixed paralog groups and weak matched-control enrichment']],[1.55,2.7,2.75])
    h('A usable interaction map')
    p(f"The U2OS AP-MS network contains {net['network_nodes']:,} proteins and {net['network_edges']:,} associations. DYNLL1 has {net['DYNLL1_unique_neighbors']} unique neighbors, including ATMIN, TP53BP1, FAM83D, HMMR, TLK1/2 and ZMYM proteins. DYNLL2 shares {net['shared_neighbors']} of those neighbors. These are associations, not proof of direct binding or a functional pathway connecting the genes to our spatial readout.")
    p('Both LC8 nodes share antibody HPA039954, the same HPA image and an explicit DYNLL1/DYNLL2 ambiguity annotation. That image cannot validate DYNLL1-specific localization. The full raw AP-MS accession is MSV000097168. The CM4AI U2OS release and original Nature study refer to this same resource and must be deduplicated.')
    h('A proteomics lead that weakens under the right controls')
    p('The qChIP study labels a ZBTB48-versus-IgG row DYNLL1, with mean log2 fold change +4.533, but its protein group includes both DYNLL1 P63167 and DYNLL2 Q96FJ2. In the matched WT/KO tables, LC8-group or resolved DYNLL1 changes are small and mixed in sign; the no-crosslinking table has no matching row. Different input amounts and crosslinking conditions are retained separately. This does not provide convincing DYNLL1-specific ZBTB48 enrichment, and averaged log-p values were not converted into new p values.')
    p('The strongest biological context remains nuclear protein assembly and the FAM83D–HMMR spindle branch. These are hypotheses to distinguish, not validated explanations for the JUMP alignment signal. We found no U2OS DYNLL1-knockout whole-proteome contrast suitable for a new pathway claim in this search.')
    refs([('U2OS cell map','https://doi.org/10.1038/s41586-025-08878-3'),('FAM83D study','https://doi.org/10.15252/embr.201847495'),('qChIP study','https://doi.org/10.1038/s41467-026-73609-9')])
    page('Search scope and data that did not qualify')
    p('The retained search covers GEO/SRA, ENA, Europe PMC, BioStudies/ArrayExpress, PRIDE/ProteomeXchange, MassIVE, iProX/NGDC leads, CM4AI, DepMap and curated perturbation catalogs. It returned 108 GEO U2OS/CRISPR series, 58 Europe PMC keyword candidates and 330 unique PRIDE records across the combined queries; 289 PRIDE records contain a U2OS spelling in retained metadata. These are search records, not counts of independent validated experiments.')
    table(['Lead','What the metadata showed'],[
        ['GSE218249','The U2OS RNA-seq targets SUGP1, XCR1 and controls. The parent paper’s DYNLL1 hit does not make this a DYNLL1 transcriptome.'],
        ['Superb-seq GSE284207','The main single-cell experiment targets four chromatin remodelers in K562; it is not a U2OS DYNLL1 cohort.'],
        ['LINCS GSE92742','69 DYNLL1 instances and 44 signatures across nine models; zero U2OS target instances or signatures.'],
        ['L1000-XPR catalog','DYNLL1 is listed as a target, but none of the 19 released model entries are U2OS.'],
        ['PerturbAtlas and CREEDS','No DYNLL1 target terms in the retained GMT files. The historical scPerturb CSV also has no U2OS record among 44 entries.'],
        ['ENA DYNLL1 records','Five returned read records belong to other assayed models. The expanded U2OS/DYNLL1 metadata query returned no records.'],
        ['TLK2 chromatin-looping study','CRA002434, HRA001605 and iProX PXD062009 are real leads, but a U2OS DYNLL1 omics contrast was not sample-verified.']],[1.65,5.35])
    p('Metadata indexing is incomplete. The broad ENA U2OS response reached its 10,000-row cap, and the scPerturb CSV is a historical snapshot. These searches cannot prove that no unindexed dataset exists. GEO, SRA and ENA mirrors are one experiment rather than separate validations. A cell line in a paper’s materials list is insufficient to assign every assay to U2OS.')
    h('Deliverables and verification')
    p(f"The ranked inventory contains {scope['curated_ledger_records']} source families, their exact accessions, eligible experiments, local files, interpretation and limitations. Full processed counts, gene-level scores, guide-overlap tables, original portal responses and provenance are retained. The 24 main quantitative and eligibility exports were recomputed and reproduced byte-identically with HTTP and connection creation disabled. No new CRISPR significance test, expression comparison or pathway-enrichment test was run in this source search.")
    p('The immediate computational priority is to analyze eligible nonoverlapping CRISPR guides with each screen’s own controls and readout, using a frozen analysis policy. The molecular publications provide context for nuclear assembly. Replication of the actual nuclear-neighbor endpoint remains a separate requirement for the manuscript. All earlier failed or nonconfirming spatial tests remain intact.')
    refs([('GEO','https://www.ncbi.nlm.nih.gov/geo/'),('Perturbation catalog','https://perturbseqr.maayanlab.cloud/download'),('ENA','https://www.ebi.ac.uk/ena/browser/home')])
    for el in [D._element,D.styles.element]:
        for border in list(el.iter(qn('w:pBdr'))):border.getparent().remove(border)
    out=ROOT/'DYNLL1_U2OS_public_validation_sources.docx';D.save(out);print(out)
if __name__=='__main__':main()
