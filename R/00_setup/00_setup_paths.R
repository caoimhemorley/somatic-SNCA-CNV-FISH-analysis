library(here)

# Data
data_raw_dir <- here("data", "raw")

# Outputs
outputs_dir <- here("outputs")
dir.create(outputs_dir, recursive = TRUE, showWarnings = FALSE)

# Plots
plots_dir <- here(outputs_dir, "plots")

plots_asyn_ns_dir       <- here(plots_dir, "aSyn_NS_vs_sections")
plots_cellcounts_dir    <- here(plots_dir, "FISH_cell_counts")
plots_correlations_dir  <- here(plots_dir, "FISH_correlations")
plots_gci_dir           <- here(plots_dir, "FISH_GCI_distribution")
plots_avg_cnv_dir       <- here(plots_dir, "FISH_MSA_v_Control_CNVs")
plots_ref_cnv_dir       <- here(plots_dir, "FISH_reference_probe_CNVs")

plots_reg_sox10pos_gain <- here(plots_dir, "FISH_regional_SOX10pos_Gains")
plots_reg_sox10neg_gain <- here(plots_dir, "FISH_regional_SOX10neg_Gains")
plots_reg_sox10pos_loss <- here(plots_dir, "FISH_regional_SOX10pos_Losses")
plots_reg_sox10neg_loss <- here(plots_dir, "FISH_regional_SOX10neg_Losses")

plots_spearman_dir <- here(plots_dir, "FISH_spearmans_plots")
plots_gh2ax_dir    <- here(plots_dir, "gH2AX")

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_asyn_ns_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_cellcounts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_correlations_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_gci_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_avg_cnv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_ref_cnv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_reg_sox10pos_gain, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_reg_sox10neg_gain, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_reg_sox10pos_loss, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_reg_sox10neg_loss, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_spearman_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_gh2ax_dir, recursive = TRUE, showWarnings = FALSE)

# Statistical output
stats_dir <- here(outputs_dir, "stat_output")
stats_summary_dir <- here(stats_dir, "Stats_Summaries")

stats_asyn_ns_dir     <- here(stats_summary_dir, "aSyn_NS_vs_sections")
stats_avg_cnv_dir     <- here(stats_summary_dir, "FISH_average_SNCA_CNVs")
stats_cellcounts_dir  <- here(stats_summary_dir, "FISH_cell_counts")
stats_demo_dir        <- here(stats_summary_dir, "FISH_cohort_demographics")
stats_gci_dir         <- here(stats_summary_dir, "FISH_GCI_distribution_MSA_subtypes")
stats_overall_cnv_dir <- here(stats_summary_dir, "FISH_overall_CNVs")
stats_ref_cnv_dir     <- here(stats_summary_dir, "FISH_reference_probe_CNV_analysis")

stats_reg_sox10pos_gain <- here(stats_summary_dir, "FISH_regional_SOX10pos_gains")
stats_reg_sox10neg_gain <- here(stats_summary_dir, "FISH_regional_SOX10neg_gains")
stats_reg_sox10pos_loss <- here(stats_summary_dir, "FISH_regional_SOX10pos_losses")
stats_reg_sox10neg_loss <- here(stats_summary_dir, "FISH_regional_SOX10neg_losses")

stats_spearman_dir    <- here(stats_summary_dir, "FISH_spearman_correlations")
stats_gh2ax_dir       <- here(stats_summary_dir, "gH2AX")
stats_assumptions_dir <- here(stats_summary_dir, "Normality_and_Variance_Test")

dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_asyn_ns_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_avg_cnv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_cellcounts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_demo_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_gci_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_overall_cnv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_ref_cnv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_reg_sox10pos_gain, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_reg_sox10neg_gain, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_reg_sox10pos_loss, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_reg_sox10neg_loss, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_spearman_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_gh2ax_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stats_assumptions_dir, recursive = TRUE, showWarnings = FALSE)
