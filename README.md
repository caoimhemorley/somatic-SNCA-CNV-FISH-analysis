
This repository contains the analysis scripts used to generate the figures, tables, and summary statistics for a manuscript investigating the role of somatic SNCA copy number variations (CNVs) in multiple system atrophy (MSA) human brain tissue.

All scripts are written in R and are organised to reflect the structure of the manuscript (main figures, supplementary figures, and tables).

The repository is intended for reproducibility and transparency of the analyses underlying the manuscript.


R/
├── 00_setup/
│   ├── 00_setup_paths.R
│   └── 01_setup_environment.R
│
├── 01_assumption_checks/
│   └── Normality_and_variance_testing.R
│
├── 02_main_analyses/
│   ├── Fig1c_Average_SNCA_CNV_MSA_vs_Control.R
│   ├── Fig2b_Reference_probe_CNV_Analysis.R
│   ├── Fig3a_Regional_SOX10pos_Gains.R
│   ├── Fig3b_Regional_SOX10neg_Gains.R
│   ├── Fig4a_Regional_SOX10pos_Losses.R
│   ├── Fig4b_Regional_SOX10neg_Losses.R
│   ├── Fig5a-d_Spearman_correlation.R
│   ├── Fig6b-f_gH2AX_analysis.R
│   └── Table1_Demographics_analysis_FISH_cohort.R
│
├── 03_supplementary_analyses/
│   ├── Supp_Fig1b_aSyn_NS_vs_sections.R
│   ├── Supp_Fig2_Summarising_Cell_Counts.R
│   ├── Supp_Fig3_GCI_distribution_SND_v_OPCA.R
│   ├── Supp_Table2_overall_cnvs.R
│   ├── Supp_Table3_Spearman_Correlations.R
│   └── Supp_Table4_Demographics_analysis_gH2AX_cohort.R



## 00_setup/

  Core setup scripts required by all analyses.

  These scripts:

- Define directory paths
- Load required R packages
- Import and clean raw data
- Generate derived datasets used throughout the manuscript


## 01_assumption_checks/

  Statistical assumption testing used to justify the choice of parametric versus non-parametric tests reported in the manuscript.

  Includes:
- Normality testing (Shapiro-Wilk)
- Variance testing (Levene's test)


## 02_main_analyses/

  Scripts corresponding directly to main manuscript figures and tables, including:

- SNCA CNV burden (MSA vs Control)
- Reference probe control analyses
- SOX10⁺ and SOX10⁻ gain and loss distributions across different brain regions
- SNCA CNV Correlation analyses
- gH2AX analyses 
- Primary FISH cohort demographics (Table 1)

## 03_supplementary_analyses/

  Scripts generating supplementary figures and tables, including:

- aSyn NS versus tissue section comparisons
- Cell count summaries
- GCI distributions by MSA subtype
- Overall CNV summaries
- Expanded correlation analyses
- gH2AX cohort demographics

## Data and outputs

- Raw input data are assumed to be located in a `data/` directory (not included in this repository).
- The input datasets used for these analyses will be provided as part of the manuscript’s supplementary material.
- All outputs (figures and summary tables) are written to structured subdirectories under `outputs/`, as defined in `00_setup/00_setup_paths.R`
