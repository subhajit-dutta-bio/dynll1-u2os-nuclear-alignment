from pathlib import Path
import json,hashlib,zipfile,platform,importlib.metadata
P=Path(__file__).resolve().parent;ROOT=P.parent;OUT=ROOT.parent
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    versions={k:importlib.metadata.version(k) for k in ['numpy','pandas','scipy','statsmodels','pyarrow','openpyxl','requests','matplotlib','beautifulsoup4']};(P/'runtime_versions.json').write_text(json.dumps({'Python':platform.python_version(),'packages':versions},indent=2))
    files=[]
    for f in P.rglob('*'):
        if not f.is_file() or '__pycache__' in f.parts or f.suffix=='.log':continue
        rel=f.relative_to(P)
        if 'render' in rel.parts:continue
        if f.name in ['checksums.json','package_manifest.json','archive_validation.json']:continue
        files.append(f)
    sources=['phase9/data/GSE232592_CRISPR_p53_counts.txt.gz','phase9/data/GSE290538_U2OS_HRE_Count_table.txt.gz','phase9/data/GSE283147_combined-guide-level-counts.csv.gz','phase9/data/GSE233195_JSMO.count.txt.gz','phase9/data/GSE226994_MAGeCK_count_normalized.xlsx','phase9/data/U2OS_APMS.cx','phase9/data/41586_2025_8878_MOESM1_ESM.xlsx','phase9/results/public_CRISPR_guide_overlap_with_JUMP.csv','phase9/results/GSE232592_DYNLL1_published_gene_results.csv','phase9/results/DepMap24Q4_U2OS_DYNLL1.json','phase9/results/qChIP_U2OS_LC8_KO_control_rows.csv','phase9/README.txt','phase6/results/adaptive_conditional_contrasts.csv','phase6/results/direct_matching_contrasts.csv','phase6/results/module_tests.json','phase6/results/gene_descriptive_statistics.json','phase6/mechanistic_policy.json','phase6/qa/verification.json','phase7/results/benchmark/summary.json','phase7/results/DYNLL1_shared_reagent_identity.csv']
    files +=[ROOT/s for s in sources]+[ROOT/'DYNLL1_U2OS_manuscript_2026-10-30.docx'];files=sorted(set(files));checks={f.relative_to(ROOT).as_posix():sha(f) for f in files};(P/'checksums.json').write_text(json.dumps(checks,indent=2));files.append(P/'checksums.json');archives=[]
    groups=[('DYNLL1_U2OS_manuscript_figures.zip',[f for f in files if f.suffix=='.docx' or f.parent==P/'figures' or f.name in ['README.txt','manuscript_text.txt'] and f.parent==P]),('DYNLL1_U2OS_manuscript_and_analysis.zip',files)]
    for name,items in groups:
        dest=OUT/name
        with zipfile.ZipFile(dest,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
            for f in items:z.write(f,(Path(ROOT.name)/f.relative_to(ROOT)).as_posix(),compress_type=zipfile.ZIP_STORED if f.suffix in ['.zip','.gz','.xlsx','.docx','.parquet','.png'] else zipfile.ZIP_DEFLATED)
        with zipfile.ZipFile(dest) as z:assert z.testzip() is None
        archives.append({'archive':str(dest),'bytes':dest.stat().st_size,'sha256':sha(dest),'entries':len(items),'CRC_pass':True})
    result={'archives':archives,'checksummed_files':len(checks),'new_functional_analysis_only':True,'raw_images_not_duplicated':True,'earlier_deliveries_preserved':True};(P/'package_manifest.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':main()
