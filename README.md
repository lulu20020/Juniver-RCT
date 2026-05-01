# Juniver RCT analysis code

This repository contains the main R analysis scripts accompanying a manuscript on the Juniver randomised controlled trial.

This is a code-only repository for transparency and code availability. The analyses were run in the authors' controlled local research environment using confidential trial data, which are not included here.

## Scripts

- `scripts/missingness_diagnostics.R`: missingness diagnostics and logistic missingness models.
- `scripts/mmrm_main_analysis.R`: primary and exploratory mixed model repeated measures analyses.
- `scripts/mmrm_subgroup_analysis.R`: baseline-defined subgroup MMRM analyses.
- `scripts/mnar_sensitivity_analysis.R`: MNAR sensitivity analyses for the primary outcome.

Earlier local scripts for data cleaning, BMI quality control, and derivation of the confidential analysis dataset are not included.

## Data

Individual participant trial data are not included because they contain confidential clinical research data.

The confidential dataset is not provided in this code-only repository. In the scripts, `data/analysis_dataset.csv` is used as a generic placeholder path for the de-identified long-format analysis dataset in the authors' controlled local analysis environment.

## Dependencies

The analyses use the following R packages: `broom`, `dplyr`, `emmeans`, `forcats`, `mmrm`, `rbmi`, `readxl`, `tibble`, and `tidyr`.
