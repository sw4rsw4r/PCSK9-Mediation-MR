# PCSK9-Mediation-MR

Scripts and variant-level tables for the real-data example of the two-step method and multivariable Mendelian randomisation (MVMR): plasma PCSK9 (exposure), LDL-C (mediator) and coronary artery disease (CAD, outcome), in two settings for the LDL-C instruments.

[Click here](https://htmlpreview.github.io/?https://github.com/sw4rsw4r/PCSK9-Mediation-MR/blob/main/PCSK9-Mediation-MR.html) for the report, which regenerates Table 3, Supplementary Table S2, Figure 1 and Supplementary Figure S1.

The full GWAS summary statistics are available from the studies cited in the paper and are not included here.

## Contents

- `summary_results/`: the variant-level tables used in the paper (PCSK9 instruments, LDL-C instruments of Settings 1 and 2, 1000 Genomes correlation matrices), and `estimates.tsv`, the estimates reported in the paper.
- `scripts/estimate_effects.py`: all estimates. IVW and multivariable IVW in generalised least squares form with the LD between instruments within 10 Mb of PCSK9, multiplicative random effects, LDL-C instruments in the APOE and LPA regions excluded in the main analysis, and covariance terms for the two-step direct effect and the MVMR indirect effect.
- `scripts/figure1.R`, `scripts/figureS1.R`: the figures.
- `PCSK9-Mediation-MR.Rmd`: the report. It runs the scripts above and checks the estimates against the MendelianRandomization R package.
- `scripts/extract_*.py`, `scripts/clump_ldl_by_z.py`, `scripts/select_cis_instruments.R`: rebuild `summary_results` from the original GWAS summary statistics (see `data/README.md`).

## Reproducing the results

```bash
conda env create -f environment.yml
conda activate pcsk9-mr
python scripts/estimate_effects.py
Rscript scripts/figure1.R
Rscript scripts/figureS1.R
Rscript -e "rmarkdown::render('PCSK9-Mediation-MR.Rmd')"
```

Outputs are written to `results/`. On Windows, run `scripts/setup_env_windows.ps1` once after creating the environment (see the script for why). The environment pins Python 3.12.9, NumPy 2.5.1, SciPy 1.18.0, pandas 3.0.3, R 4.5.3, MendelianRandomization 0.10.0 and ggplot2 4.0.4. PLINK v1.90b7.7 is needed only to rebuild the tables.
