DYNLL1 U2OS correction release 1.0.1
1 November 2026

Correction DOI: https://doi.org/10.5281/zenodo.23856046
Preserved base release: https://doi.org/10.5281/zenodo.23850589
Code: https://github.com/subhajit-dutta-bio/dynll1-u2os-nuclear-alignment

Extract DYNLL1_U2OS_public_data_v1.0.0.zip, then extract
DYNLL1_U2OS_correction_v1.0.1.zip into the same parent directory. Both archives
use dynll1-u2os-nuclear-alignment/ as their top-level directory. This adds a
phase18/ directory without replacing the original phase17/ release. The base
ZIP and its prior full statistical reproduction report are unchanged.

The corrected scientific figures are in phase18/figures/. The associated
88 panel-source files and registry are in phase18/panel_data/ and phase18/.
All 22 supplementary tables retain their original values.

Corrections:
1. Supplementary Figure 2A moves the DYNLL1 nuclear label into open space,
   with a leader to its original point.
2. Supplementary Figure 3D separates the two-line Spearman annotation from
   the vertical zero line and moves target labels away from the point cloud.
3. Supplementary Figure 5F distinguishes original across-replicate guide
   eligibility (*) from an individual replicate's low denominator (dagger).

ACAAGAA Rep2 has raw High=1480 and Low=0. Its mean normalized Low across three
replicates is 92.7065, with adequate coverage in two of three replicates; it
therefore passes the original guide rule. Its Rep2 denominator remains zero.
Library normalization gives High=1284.2616 and Low=0, so pseudocount one gives
log2((1284.2616+1)/(0+1))=10.32785. Subtracting the NTC median -0.654895 gives
10.98274. This is a pseudocount-dependent zero-denominator value rather than
a well-measured approximately 2000-fold enrichment.

The dagger flags all four zero-denominator cells, including two zero/zero
cells whose displayed +0.654895 arises entirely from NTC centering. The
asterisk remains on CATTGCG, which fails the original across-replicate rule.
All 12 original ratios, original guide eligibility and all primary spatial
results are unchanged. Original and nonconfirming outcomes remain available.

Reproduce the full raw-count audit and the three corrected figures offline:
  python -X utf8 phase18/audit_hif_coverage.py
  python -X utf8 phase18/plot_corrected.py

Rebuild all 11 figures with the corrected annotations:
  python -X utf8 phase18/plot_all.py

The all-12-cell audit, pseudocount sensitivities and replicate-coverage
diagnostic are in phase18/results/HIF_*. Pseudocount 0.1, 1 and 10 are reported
with NTC centers recomputed for each choice; the original guide eligibility
is held fixed. The separate covered-only diagnostic masks guide/replicate
cells with normalized Low<10 within originally eligible guides and retains
the original NTC centering. Its three replicate medians have mean -4.69245,
compared with the preserved primary -3.94008. This is explicitly post hoc,
does not replace the original estimate and supplies no new significance test
or validation of nuclear organization.

The same MIT code and CC BY 4.0 derived-output terms apply as in the base
release; upstream source licenses are unchanged. This is a scientific data
and figure correction, not a journal submission or new biological validation.
