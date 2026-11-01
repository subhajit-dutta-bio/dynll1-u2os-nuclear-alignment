"""Plain-language reading guide and response to the substantive review."""
from pathlib import Path
import json,sys,pandas as pd
P=Path(__file__).resolve().parent;ROOT=P.parent;R=P/'results'
def main():
    if hasattr(sys.stdout,'reconfigure'):sys.stdout.reconfigure(encoding='utf-8')
    assert json.loads((P/'qa/numeric_verification.json').read_text())['all_pass'];s=json.loads((R/'comparison_summary.json').read_text());g=pd.read_csv(R/'gene_summary.csv');d=g[g.gene.eq('DYNLL1')].iloc[0];m=s['groups']['matched_other_chromosome'];p=s['groups']['proximal_12q'];motor=s['motor_module'];e=pd.read_csv(R/'alignment_contrasts.csv');prediction=pd.read_csv(R/'prediction/summary.csv')
    response=f'''VERSION 3 REVISION AND COMPLETED TARGETED ANALYSIS
30 October 2026

NEW MEASUREMENTS
All 7,875 native DNA fields from the fixed U2OS cohort were processed.
The design contains 67 genes, 335 target wells and 875 total wells across
the original five acquisition dates. No planned gene or date was replaced.
All analyses retain their missing values and fixed eligibility thresholds.

DYNLL1 mean primary contrast: {d['mean']:+.6f}
Evaluable dates: {int(d.dates)}; positive dates: {int(d.positive_dates)}
Complete primary genes: {s['complete_genes']} / 67
Matched comparators complete on all dates: {m['complete_genes']} / 32
DYNLL1 rank including matched comparators: {m['DYNLL1_rank_within_set_plus_anchor']}
Matched median: {m['median']:+.6f}; matched 95th percentile: {m['q95']:+.6f}
Descriptive specificity rule met: {s['DYNLL1_specificity_descriptive_rule_met']}
Minimum DYNLL1 source-control density support: {s['minimum_DYNLL1_density_support']:.4f}
Proximal12q complete genes: {p['complete_genes']} / 23
DYNLL1 rank including proximal12q comparators: {p['DYNLL1_rank_within_set_plus_anchor']}
Proximal12q median: {p['median']:+.6f}
Complete exploratory motor module evaluable: {motor['evaluable']}
Motor-module mean: {motor['mean']:+.6f}; positive dates: {motor['positive_dates']}
Motor-module exact two-sided sign-flip fraction: {motor['two_sided_tail']:.4f}

The strongest descriptive pattern is joint: DYNLL1 has the highest mean
local alignment and solidity, and the lowest mean nuclear aspect ratio
and axis-position coupling among these 67 selected genes. These extrema
were observed after the panel was complete, not four independent tests.
Alignment remains positive on all five dates after quadratic-area
adjustment and the unchanged direct nucleus-level matching procedure.

WHAT THE REVISION ADDRESSES
1. Selection history. The original discovery was a selected 64-gene pilot,
not all 7,975 genes. Its 126 tests and ten-gene frozen shortlist are now
explicit. The earlier source-feature split used three discovery and two
later dates; current raw-image follow-up reused all known dates. With five
dates, the minimum exact two-sided sign-flip fraction remains 0.0625.

2. Size and fitness. All four spacing models are retained. Earlier area
adjustment substantially attenuated spacing. The new 32 matched comparator
genes use fixed source-count/area and U2OS-fitness proxies, with actual
image covariate balance reported rather than assumed. DYNLL1 raw-image
area and aspect-ratio means fall outside the matched-set interquartile
ranges. A transparently post-outcome diagnostic of the existing focal
matches finds small residual log-covariate differences; no rematching,
threshold or effect estimate changed. All-site count/area QC retains
sparse fields without relaxing collective eligibility.

3. Motor alternative. DYNC1H1, DYNC1I2 and PAFAH1B1 are measured on every
planned date. Their joint alignment pattern does not reproduce DYNLL1.
This does not establish motor independence without measured perturbation
severity. NUMA1, NUDC and DYNLL2 were absent from the CRISPR metadata.
Missing coverage is not a null result or a reason to change the module.
The earlier TP53BP1/H2AX/RNF8/RNF168 hypothesis remains nonconfirming.

4. Chromosome context. All 23 proximal and eight distal 12q targets now have
actual raw-image spatial measurements. Lazar 2024 is cited. This diagnostic
uses median intended Brunello cut positions as gene-location proxies,
not the actual JUMP cut sites. It does not measure U2OS chromosome structure or establish a validated
proximity-bias correction for the new statistic.

5. Functional context. GEO labels are corrected to parental TERT+/+ versus
heterozygous TERT-mutant U2OS. The series describes two clones, without
explicit sample-title clone mapping. The TERT screen, AP-MS and unrelated
abundance screens are supplementary context. Shared guides, filter
sensitivity, failed gene sets and original published FDRs remain visible.

6. Measurement and presentation. Native Cellpose-SAM is the new primary
method. Every fixed JUMP crop uses a 65-micrometer bar. The text defines
unoriented nuclear directors and the cos(2Δθ) scale without converting it to
a unique angle or percentage. Five main and five supplementary figures
each contain eight panels, with 80 panel sources and 17 supplementary tables.
The final Cell 2025 citation replaces the ambiguous preprint citation.

REMAINING SCIENTIFIC LIMITS
This is a targeted post-discovery comparison, not a genome-wide spatial
atlas or a validated gene-specific mechanism. Dates share one reagent pool
per gene. Independent spatial target sites, measured protein depletion,
rescue, cell dynamics and actual chromosome structure are not supplied by
these public datasets. Author identity, affiliation, correspondence,
no funding and no competing interests were supplied by Subhajit Dutta.
No HPC acknowledgment is included. Panel count is not independent
evidence count.
'''
    (P/'revision_response.txt').write_text(response,encoding='utf-8')
    readme='''NUCLEAR ORGANIZATION IN U2OS VERSION 3

READ FIRST
  DYNLL1_U2OS_main_manuscript_v3.docx
  DYNLL1_U2OS_supplementary_information_v3.docx
  phase14/revision_response.txt

The completed fixed cohort contains 67 genes, 875 wells and 7,875 native DNA
fields across five dates. Every planned contrast and missing value remains.
Five main figures and five supplementary figures have eight panels each.
The supplement contains 17 tables. Exact results and limits are in the
manuscript and revision response; no independent experimental validation
or genome-wide spatial atlas is claimed.

FIGURES AND TABLES
  phase14/figures: individual PNG, SVG and PDF figures plus two collections.
  phase14/panel_data: 80 exact numeric panel sources.
  phase14/panel_registry.json: input pointers, units and interpretation.
  phase14/tables: all 17 exported supplementary tables.

OFFLINE NUMERIC REPRODUCTION
Keep the extracted folder structure. Install requirements_numeric.txt in
your scientific Python environment, then run from the project directory:
  python phase14/reproduce_numeric.py
  python phase14/verify_numeric.py
  python phase14/plot_all.py

This refits all new control models, all 335 target/date contrasts under
each specification, matched/chromosome comparisons, all-site covariates
and all nested prediction models from the complete measured focal/field
inputs. Historical screen-context panels reproduce from retained earlier
numeric sources. The relevant original public count inputs, historical
analysis code and fixed policies are also included for inspection and
separate reanalysis; the new offline validation does not rerun every
historical assay. Figures include only the fixed display pixels needed
for regeneration, not the entire image collection.

DOCUMENTS
  python phase14/build_documents.py
DOCX authoring used the Python document runtime. Native Word rendered
the final documents for page-by-page review after the packaged renderer
reported that LibreOffice was unavailable. QA previews are internal and
are excluded from the reading deliverables.

RAW IMAGE LAYER
The full downloaded TIFF, mask and per-object geometry collection remains
locally available. Bulk pixels and the large checkpoint are excluded from
the compact numeric archive. Public URLs and hashes for every field, the
exact checkpoint identity and native inference settings are retained.
Repeating inference requires those public pixels, the exact checkpoint
and an appropriate GPU/software environment. checkpoint_source.json gives
the public checkpoint URL and its independently published matching SHA-256.
Install the retained numeric and image requirements; the CUDA PyTorch
wheel is available from the official PyTorch CUDA 12.8 wheel index.
The portable image-layer entry points are:
  python phase14/fetch_pixels.py --workers 6
  python phase14/extract.py --workers 6 --weights PATH_TO_CPSAM_V2
The original source used for this run is preserved as extract_source_used.py.
The later checkpoint-path option and output-directory initialization do not
change image settings or measurements. No target-specific training,
resizing or threshold tuning was used.

PROVENANCE
analysis_policy.json and frozen_plan.json define the unchanged primary
comparison. coverage_QC_addendum.json transparently records the later
all-site covariate diagnostic prompted by sparse fixed example fields.
matching_balance_addendum.json records the post-outcome diagnostic of
the unchanged focal matches, including its exact covariate definition.
prediction_policy.json defines every evaluated model and grouped fold.
The fresh archive validation blocks network access in every child process.
All earlier manuscripts and archives remain unchanged. Earlier version2
documents may be present as preservation-check inputs in the analysis
archive; version3 is the current manuscript.
'''
    (P/'README.txt').write_text(readme,encoding='utf-8');print(response)
if __name__=='__main__':main()
