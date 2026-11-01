# DYNLL1 targeting and nuclear alignment in U2OS

Current release: **1.0.1**, 1 November 2026. Data and correction: [Zenodo](https://doi.org/10.5281/zenodo.23856046).

This release corrects label placement in Supplementary Figures 2A and 3D and adds explicit replicate-level denominator warnings to Supplementary Figure 5F. ACAAGAA Rep2 has raw High=1480 and Low=0; its displayed +10.9827 depends on a pseudocount and NTC centering. Its guide passes the original across-replicate rule, while this individual replicate has a zero denominator. All four such cells are flagged. All 12 original ratios, guide eligibility and primary spatial results are preserved.

Extract the unchanged v1.0.0 base ZIP and overlay the v1.0.1 correction ZIP into the same parent directory. Both use `dynll1-u2os-nuclear-alignment/` as their root. Then run:

```sh
python -X utf8 phase18/audit_hif_coverage.py
python -X utf8 phase18/plot_all.py
```

See [correction details](phase18/README_correction.txt) and the [fresh offline reproduction report](phase18/qa/correction_offline_reproduction.json). The preserved [v1.0.0 dataset](https://doi.org/10.5281/zenodo.23850589) and commit cbca0c4a14e209c53ed81ad01633deefb68246a0 remain available. Original base-release instructions are in README.txt.

This is a descriptive reanalysis of public U2OS images. Independent targeting, rescue and the molecular mechanism remain unvalidated. Dates share reagents and controls. All selected genes, dates and nonconfirming outcomes are retained. The data release excludes manuscript drafts and private review correspondence.

Original code is MIT licensed; original derived outputs are CC BY 4.0. Upstream source terms remain unchanged.
