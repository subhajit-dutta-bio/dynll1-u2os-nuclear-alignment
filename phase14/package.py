"""Build the version3 reading and portable numeric archives after QA."""
from pathlib import Path
import json,zipfile,hashlib,importlib.metadata as im,sys,platform
P=Path(__file__).resolve().parent;ROOT=P.parent;OUT=ROOT.parent
DOCS=['DYNLL1_U2OS_main_manuscript_v3.docx','DYNLL1_U2OS_supplementary_information_v3.docx']
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()
def inventory(kind):
    files={ROOT/n for n in DOCS};files.add(P/'README.txt')
    for folder in [P/'figures',P/'panel_data',P/'tables']:files.update(f for f in folder.rglob('*') if f.is_file())
    for name in ['panel_registry.json','panel_registry.csv','revision_response.txt']:files.add(P/name)
    if kind=='analysis':
        files.update(f for f in P.iterdir() if f.is_file() and f.suffix not in ['.log'])
        # Complete measured inputs and results; bulk masks/focals/geometry are
        # reconstructed from public pixels rather than copied into this layer.
        files.update(f for f in (P/'results').iterdir() if f.is_file())
        files.update(f for f in (P/'results/prediction').iterdir() if f.is_file())
        files.update(f for f in (P/'data').iterdir() if f.is_file())
        files.update(f for f in (P/'data/recent_portal').iterdir() if f.is_file())
        files.update(f for f in (P/'data/provenance').glob('*.json'))
        files.update(f for f in (P/'qa').iterdir() if f.is_file() and f.suffix=='.json' and not f.name.startswith('static_') and f.name not in ['archive_validation.json','delivery_archives.json'])
        for row in json.loads((P/'qa/observed_local_reads.json').read_text()):files.add(ROOT/row['path'])
        # Some binary readers open files below Python's audit-hook layer.
        # Include every fixed display image and geometry explicitly as well.
        import pandas as pd
        design=pd.read_csv(P/'frozen_fields.csv').fillna({'gene':'','panel':''})
        display=design[design.gene.isin(['DYNLL1','DYNC1H1','DYNC1I2','PAFAH1B1'])&design.batch.eq('20220914_Run1')&design.site.eq(0)]
        for row in display.itertuples():
            files.add(P/'data/images'/f'{row.field}.tif');files.add(P/'results/geometry'/f'{row.field}.parquet')
            source=f'{row.plate}__A02__s0';files.add(P/'data/images'/f'{source}.tif')
            if row.gene=='DYNLL1':files.add(P/'results/geometry'/f'{source}.parquet')
        for name in json.loads((P/'frozen_plan.json').read_text())['dependencies']:files.add(ROOT/name)
        # The style helper imports the retained content module; version3
        # replaces its prose and figure titles before authoring.
        for name in ['phase11/build_documents.py','phase11/manuscript_content.py','phase10/data/crossref_references.json','phase13/results/proposed_matched_comparators.csv','phase13/data/matching_profiles.parquet','phase13/results/all_gene_coordinate_coverage.csv.gz','phase12/frozen_policy.json','phase7/results/U2OS_independent_validation_ledger.csv']:
            files.add(ROOT/name)
        for name in ['DYNLL1_U2OS_main_manuscript_v2.docx','DYNLL1_U2OS_supplementary_information_v2.docx','phase2/results/frozen_candidates.csv','phase2/results/heldout_validation.csv','phase6/results/adaptive_conditional_contrasts.csv','phase7/results/benchmark/summary.json','phase12/results/cellpose_benchmark/summary.json','phase12/results/orf/all_contrasts.csv','phase10/results/hTERT_DYNLL1_DYNLL2_reagent_and_interaction_audit.csv','phase10/results/DYNLL1_DYNLL2_all_contrasts.csv','phase9/results/GSE232592_DYNLL1_published_gene_results.csv','phase11/results/APMS_hTERT_matched_gene_sets.csv','phase9/data/DepMap24Q4_U2OS_CRISPRGeneEffect.csv','phase10/data/Brunello_official_library.txt']:
            files.add(ROOT/name)
        # Retain the original public count inputs and historical analysis code
        # as well as the numeric context panels. The new offline validator
        # explicitly refits the new spatial study, not every historical assay.
        for folder in ['phase9','phase10','phase11']:
            files.update(f for f in (ROOT/folder).iterdir() if f.is_file() and f.suffix in ['.py','.json','.txt'])
        files.update(f for f in (ROOT/'phase10/data').iterdir() if f.is_file() and f.suffix not in ['.html','.xml'])
        files.update(f for f in (ROOT/'phase10/results').iterdir() if f.is_file())
        for name in ['GSE232592_CRISPR_p53_counts.txt.gz','GSE290538_U2OS_HRE_Count_table.txt.gz','GSE283147_combined-guide-level-counts.csv.gz','GSE233195_JSMO.count.txt.gz','GSE226994_MAGeCK_count_normalized.xlsx','U2OS_APMS.cx','41586_2025_8878_MOESM1_ESM.xlsx']:
            files.add(ROOT/'phase9/data'/name)
        for name in ['public_CRISPR_guide_overlap_with_JUMP.csv','qChIP_U2OS_LC8_KO_control_rows.csv']:
            files.add(ROOT/'phase9/results'/name)
    return sorted(f for f in files if f.is_file() and '__pycache__' not in f.parts and f.suffix!='.pyc' and not f.name.startswith('static_methods_proof'))
def runtime():
    names=['numpy','pandas','scipy','scikit-learn','scikit-image','matplotlib','pyarrow','tifffile','Pillow','requests','python-docx','pypdf','PyMuPDF','cellpose','torch','statsmodels','openpyxl'];versions={n:im.version(n) for n in names};(P/'runtime_versions.json').write_text(json.dumps({'python':sys.version,'platform':platform.platform(),'versions':versions,'document_authoring':'Bundled Codex artifact Python; see document_QA.json'},indent=2));(P/'requirements_numeric.txt').write_text('\n'.join(f'{n}=={versions[n]}' for n in names if n not in ['cellpose','torch','python-docx','pypdf','PyMuPDF'])+'\n');(P/'requirements_images.txt').write_text('\n'.join(f'{n}=={versions[n]}' for n in ['cellpose','torch','scikit-image','tifffile'])+'\n');(P/'requirements_documents.txt').write_text('python-docx\nPillow\npypdf\npandas\n')
def make(kind):
    target=OUT/('DYNLL1_U2OS_manuscript_v3_figures.zip' if kind=='reading' else 'DYNLL1_U2OS_v3_analysis.zip');files=inventory(kind);manifest=[{'path':f.relative_to(ROOT).as_posix(),'bytes':f.stat().st_size,'sha256':sha(f)} for f in files];prefix=ROOT.name+'/'
    with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
        for f in files:z.write(f,prefix+f.relative_to(ROOT).as_posix())
        z.writestr(prefix+'SHA256_MANIFEST.json',json.dumps(manifest,indent=2))
    with zipfile.ZipFile(target) as z:assert z.testzip() is None
    return {'archive':str(target),'files':len(files),'bytes':target.stat().st_size,'sha256':sha(target)}
def main():
    assert json.loads((P/'qa/document_QA.json').read_text())['all_pass'];runtime();out=[make(kind) for kind in ['reading','analysis']];(P/'qa/delivery_archives.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
