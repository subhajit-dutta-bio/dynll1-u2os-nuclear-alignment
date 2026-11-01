DYNLL1 targeting and nuclear alignment in U2OS
Release 1.0.0, 1 November 2026
Author: Subhajit Dutta

Data and reproducible analysis: https://doi.org/10.5281/zenodo.23850589
Code mirror: https://github.com/subhajit-dutta-bio/dynll1-u2os-nuclear-alignment

This is a descriptive reanalysis of public U2OS imaging. It measures local
nuclear alignment after one DYNLL1 reagent pool alongside a telomeric 12q
pattern. Independent targeting, rescue and the molecular mechanism remain
unvalidated. Dates share reagents and controls; nuclei are not independent
gene-level biological replicates. All selected genes, dates and nonconfirming
outcomes are retained. The narrow two-date rank lead is not statistically
resolved. Same-layout positions repeat across acquisition dates.

The Zenodo archive contains this code, complete focal/field inputs for the
primary statistical analyses, source metadata, all retained estimates,
profile-prediction inputs and folds, public-context results, all 88 panel
source tables, 22 supplementary tables, figure exports and provenance.
Editorial drafts and review correspondence are excluded. This release is
data and analysis code, not a journal submission or a claim of acceptance.

Unzip the Zenodo archive and work inside dynll1-u2os-nuclear-alignment/.
Install the dependencies from requirements.txt. CPU statistical reproduction
does not require GPU inference or downloading the bulk JUMP image collection.

Recompute all primary calibrations, six alignment definitions, shape and
spacing diagnostics, matching diagnostics and grouped profile predictions:
  python -X utf8 phase14/reproduce_numeric.py

Recompute the plate/rank audit, cached reference audit and same-plate checks:
  python -X utf8 phase17/audit_plate_and_rank.py
  python -X utf8 phase17/audit_references.py
  python -X utf8 phase17/same_plate_comparison.py

Rebuild all 11 figures and their 88 source tables:
  python -X utf8 phase17/plot_all.py

The 25-entry reference matrix is phase17/results/reference_checks_all_25.csv.
Nineteen journal papers have exact-DOI PubMed records plus publisher-deposited
Crossref web metadata. Six resource/preprint entries have official-source
checks. Cached records make the audit reproducible offline; --fetch is optional.

For a cloned GitHub checkout, use the Zenodo archive to supply the data/results
directories at the same relative paths. The archive also includes the exact
code so it can be used without GitHub. The GitHub source mirror omits large
data assets; the Zenodo archive is the canonical data distribution.

Raw-image inference is a separate layer. The public JUMP URLs, raw hashes,
checkpoint identity and settings are preserved in phase14. A small set of
fixed display images and the CC0 BBBC039 reference archive support figure
reproduction. Repeating all raw-pixel inference requires the public inputs,
the documented Cellpose-SAM checkpoint and an appropriate GPU environment.

See SHA256_MANIFEST.json for checksums and THIRD_PARTY_NOTICES.txt for source
attribution. MIT applies to original analysis code; CC BY 4.0 applies to
original derived release outputs. Upstream material retains its own terms.
No credentials, private review correspondence or API access tokens are included.
