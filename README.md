# DYNLL1 targeting and nuclear alignment in U2OS

DYNLL1 targeting and nuclear alignment in U2OS
Release 1.0.2, 1 November 2026
Author: Subhajit Dutta

Data and reproducible analysis: https://doi.org/10.5281/zenodo.23863910
Code: https://github.com/subhajit-dutta-bio/dynll1-u2os-nuclear-alignment

This self-contained distribution consolidates the complete core statistical
inputs and corrected scientific figures. Extract it and work inside
dynll1-u2os-nuclear-alignment/. Install requirements.txt, or use the retained
version records and requirements files in phase14 for the numeric workflow.
Current figures, 88 panel-source tables and 22 supplementary table exports
are in phase19. The main spatial estimates and all twelve HIF guide/replicate
ratios and original eligibility calls remain unchanged.

Recompute primary calibrations, all six alignment definitions, matching,
shape/spacing diagnostics and grouped conventional-profile predictions:
  python -X utf8 phase14/reproduce_numeric.py

Recompute the explicit HIF raw-count and replicate-coverage audit:
  python -X utf8 phase19/audit_hif_coverage.py

Recompute the 25-reference identity audit from cached bibliographic records:
  python -X utf8 phase19/audit_references.py

Rebuild all eleven figures and 88 panel-source tables:
  python -X utf8 phase19/plot_all.py

The current 25-row reference matrix is phase19/results/reference_checks_all_25.csv.
Nineteen journal references have PubMed exact-DOI records and Crossref
metadata; six resource/preprint entries have official-source checks.
Bibliographic records are distributed as projections. Nested cited-work
bibliographies are omitted. The U2OS Cell Map dataset title is abbreviated
in display, with its exact DOI and official source URL retained. These edits
do not change the identity of a cited paper or resource or validate a
biological interpretation.

ACAAGAA Rep2 in the HIF context has raw High=1480 and Low=0. Its guide passes
the across-replicate rule, while this particular denominator is zero. The
displayed +10.9827 depends on a pseudocount and nontargeting centering. All
four zero-denominator cells are explicitly flagged. Complete counts, flags
and clearly labeled post hoc sensitivities are retained in phase19/results.

This is a descriptive public-data study of one reagent pool per gene.
Acquisition dates reuse reagents and controls. Independent targeting, rescue
and the molecular mechanism remain unvalidated. All selected genes, dates
and nonconfirming outcomes are retained. This distribution is scientific
data and code, not a journal submission or new wet-lab validation.

Raw-pixel inference is a separate layer. Bulk public JUMP images and the
retained segmentation checkpoint are linked by URLs and hashes. Repeating
the full image inference needs those public inputs and a suitable runtime.
The CPU statistical and figure-source workflows use the bundled measured
inputs and require no network access or new pixel inference.

Original code is MIT licensed and original derived outputs are CC BY 4.0.
Upstream data retain their original terms; see THIRD_PARTY_NOTICES.txt.
See SHA256_MANIFEST.json for the distributed-file checksums.
