source(here("R", "00_setup", "01_setup_environment.R"))

keep_existing_vars <- function(data, vars) {
  vars[vars %in% colnames(data)]
}

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.0001 ~ "p<0.0001",
    p < 0.001 ~ "p<0.001",
    p < 0.01 ~ "p<0.01",
    TRUE ~ paste0("p=", format(round(p, 3), nsmall = 3))
  )
}

shapiro_test_function <- function(data, value_col, group_cols) {
  data %>%
    dplyr::select(dplyr::all_of(c(group_cols, value_col))) %>%
    tidyr::drop_na(dplyr::all_of(value_col)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::group_modify(~ {
      values <- .x[[value_col]]
      values <- values[!is.na(values)]
      
      if (length(unique(values)) > 1 && length(values) >= 3) {
        test <- stats::shapiro.test(values)
        tibble::tibble(
          Variable = value_col,
          Test = "Shapiro",
          W = as.numeric(test$statistic),
          p_value = test$p.value,
          p_text = format_p(test$p.value),
          note = NA_character_
        )
      } else {
        tibble::tibble(
          Variable = value_col,
          Test = "Shapiro",
          W = NA_real_,
          p_value = NA_real_,
          p_text = NA_character_,
          note = "Not enough variability / n<3"
        )
      }
    }) %>%
    dplyr::ungroup()
}

levene_test_function <- function(data, value_col, rhs_terms) {
  df <- data %>%
    dplyr::select(dplyr::all_of(c(rhs_terms, value_col))) %>%
    dplyr::rename(Value = dplyr::all_of(value_col)) %>%
    tidyr::drop_na(Value)
  
  if (nrow(df) > 0 && dplyr::n_distinct(df[[rhs_terms[1]]]) > 1) {
    test <- car::leveneTest(
      stats::as.formula(paste("Value ~", paste(rhs_terms, collapse = "*"))),
      data = df
    )
    p <- as.numeric(test[1, "Pr(>F)"])
    tibble::tibble(
      Variable = value_col,
      Test = "Levene",
      Comparison = paste(rhs_terms, collapse = "*"),
      W = NA_real_,
      p_value = p,
      p_text = format_p(p),
      note = NA_character_
    )
  } else {
    tibble::tibble(
      Variable = value_col,
      Test = "Levene",
      Comparison = paste(rhs_terms, collapse = "*"),
      W = NA_real_,
      p_value = NA_real_,
      p_text = NA_character_,
      note = "Not enough groups / no data"
    )
  }
}

write_assumptions_per_variable <- function(domain_prefix, var_name, shapiro_df, levene_df, out_dir) {
  out <- dplyr::bind_rows(shapiro_df, levene_df)
  
  out_file <- file.path(
    out_dir,
    paste0(domain_prefix, "_Shapiro_and_Levene_", var_name, ".tsv")
  )
  
  readr::write_tsv(out, out_file, na = "")
}

out_dir <- stats_assumptions_dir

cnv_regional_vars <- c(
  "Percentage_SOX10pos_gain",
  "Percentage_SOX10pos_loss",
  "Percentage_SOX10pos_inclusions",
  "Percentage_SOX10neg_gain",
  "Percentage_SOX10neg_loss",
  "Percentage_SOX10neg_inclusions",
  "Percentage_Overall_gain",
  "Percentage_Overall_loss"
)
cnv_regional_vars <- keep_existing_vars(cnv_data, cnv_regional_vars)

for (v in cnv_regional_vars) {
  sh <- shapiro_test_function(cnv_data, v, c("Region", "Classification"))
  lv <- levene_test_function(cnv_data, v, c("Region", "Classification"))
  write_assumptions_per_variable(
    domain_prefix = "CNV_regional",
    var_name = v,
    shapiro_df = sh,
    levene_df = lv,
    out_dir = out_dir
  )
}

avg_cnv_vars <- c(
  "Percentage_Average_SOX10pos_gain",
  "Percentage_Average_SOX10pos_loss",
  "Percentage_Average_SOX10neg_inclusions",
  "Percentage_Average_SOX10neg_gain",
  "Percentage_Average_SOX10neg_loss",
  "Percentage_Average_Overall_gain",
  "Percentage_Average_Overall_loss"
)
avg_cnv_vars <- keep_existing_vars(avg_cnv_data, avg_cnv_vars)

clinical_vars <- c(
  "Age_of_onset",
  "Age_of_death",
  "Disease_duration",
  "PMI"
)
clinical_vars <- keep_existing_vars(avg_cnv_data, clinical_vars)

avg_cnv_msa_v_control <- avg_cnv_data %>%
  dplyr::mutate(
    Group_MSA_vs_Control = dplyr::case_when(
      Classification %in% c("SND", "OPCA") ~ "MSA",
      TRUE ~ "Control"
    ),
    Group_MSA_vs_Control = factor(Group_MSA_vs_Control, levels = c("Control", "MSA"))
  )

for (v in avg_cnv_vars) {
  sh <- shapiro_test_function(avg_cnv_msa_v_control, v, c("Group_MSA_vs_Control"))
  lv <- levene_test_function(avg_cnv_msa_v_control, v, c("Group_MSA_vs_Control"))
  write_assumptions_per_variable(
    domain_prefix = "AvgCNV_MSA_vs_Control",
    var_name = v,
    shapiro_df = sh,
    levene_df = lv,
    out_dir = out_dir
  )
}

for (v in clinical_vars) {
  sh <- shapiro_test_function(avg_cnv_data, v, c("Classification"))
  lv <- levene_test_function(avg_cnv_data, v, c("Classification"))
  write_assumptions_per_variable(
    domain_prefix = "All_groups",
    var_name = v,
    shapiro_df = sh,
    levene_df = lv,
    out_dir = out_dir
  )
}

gh2ax_vars <- c(
  "Percentage_all_gh2ax_foci",
  "Percentage_SOX10_pos_gh2ax_pos",
  "Percentage_SOX10_neg_gh2ax_pos",
  "Percentage_aSyn_pos_gh2ax_pos",
  "Percentage_aSyn_neg_gh2ax_pos"
)
gh2ax_vars <- keep_existing_vars(gh2ax_data, gh2ax_vars)

for (v in gh2ax_vars) {
  sh <- shapiro_test_function(gh2ax_data, v, c("group_msa"))
  lv <- levene_test_function(gh2ax_data, v, c("group_msa"))
  write_assumptions_per_variable(
    domain_prefix = "gH2AX",
    var_name = v,
    shapiro_df = sh,
    levene_df = lv,
    out_dir = out_dir
  )
}



asyn_vars <- c("Percent_inclusions_sections", "Percent_inclusions_NS")
asyn_vars <- keep_existing_vars(aSyn_NS_data, asyn_vars)

if (length(asyn_vars) == 2) {
  asyn_long <- aSyn_NS_data %>%
    dplyr::select(Donor_ID, Region, Subtype, dplyr::all_of(asyn_vars)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(asyn_vars),
      names_to = "Tissue",
      values_to = "Value"
    ) %>%
    dplyr::mutate(
      Tissue = dplyr::case_when(
        Tissue == "Percent_inclusions_sections" ~ "Sections",
        Tissue == "Percent_inclusions_NS" ~ "NS",
        TRUE ~ NA_character_
      ),
      Tissue = factor(Tissue, levels = c("Sections", "NS"))
    ) %>%
    tidyr::drop_na(Tissue)
  
  sh_asyn <- shapiro_test_function(
    asyn_long %>% dplyr::rename(Percent_inclusions = Value),
    "Percent_inclusions",
    c("Tissue")
  ) %>%
    dplyr::mutate(Comparison = NA_character_) %>%
    dplyr::relocate(Tissue, .before = Variable)
  
  lv_asyn <- {
    df <- asyn_long %>% tidyr::drop_na(Value)
    if (nrow(df) > 0 && dplyr::n_distinct(df$Tissue) > 1) {
      test <- car::leveneTest(Value ~ Tissue, data = df)
      p <- as.numeric(test[1, "Pr(>F)"])
      tibble::tibble(
        Tissue = "All",
        Variable = "Percent_inclusions",
        Test = "Levene",
        Comparison = "Tissue",
        W = NA_real_,
        p_value = p,
        p_text = format_p(p),
        note = NA_character_
      )
    } else {
      tibble::tibble(
        Tissue = "All",
        Variable = "Percent_inclusions",
        Test = "Levene",
        Comparison = "Tissue",
        W = NA_real_,
        p_value = NA_real_,
        p_text = NA_character_,
        note = "Not enough groups / no data"
      )
    }
  }
  
  out_asyn <- dplyr::bind_rows(sh_asyn, lv_asyn)
  
  readr::write_tsv(
    out_asyn,
    file.path(out_dir, "aSyn_NS_vs_Sections_ShapWilk_and_Levene_Percent_inclusions.tsv"),
    na = ""
  )
}

gci_var <- "Percentage_SOX10pos_inclusions"

if (gci_var %in% colnames(cnv_data)) {
  gci_data <- cnv_data %>%
    dplyr::filter(
      Classification %in% c("SND", "OPCA"),
      Region %in% c("Putamen", "Cerebellum")
    ) %>%
    dplyr::select(Donor_ID, Classification, Region, dplyr::all_of(gci_var)) %>%
    tidyr::drop_na(dplyr::all_of(gci_var))
  
  gci_shapiro <- shapiro_test_function(gci_data, gci_var, c("Classification", "Region")) %>%
    dplyr::mutate(Comparison = NA_character_)
  
  levene_within_subtype <- purrr::map_dfr(
    c("SND", "OPCA"),
    function(cls) {
      df_sub <- gci_data %>% dplyr::filter(Classification == cls)
      levene_test_function(df_sub, gci_var, c("Region")) %>%
        dplyr::mutate(
          Classification = cls,
          Region = NA_character_,
          Comparison = paste0("Region (Putamen vs Cerebellum) within ", cls)
        )
    }
  )
  
  levene_within_region <- purrr::map_dfr(
    c("Putamen", "Cerebellum"),
    function(rg) {
      df_sub <- gci_data %>% dplyr::filter(Region == rg)
      levene_test_function(df_sub, gci_var, c("Classification")) %>%
        dplyr::mutate(
          Classification = NA_character_,
          Region = rg,
          Comparison = paste0("Subtype (SND vs OPCA) within ", rg)
        )
    }
  )
  
  gci_out <- dplyr::bind_rows(
    gci_shapiro,
    dplyr::bind_rows(levene_within_subtype, levene_within_region)
  )
  
  readr::write_tsv(
    gci_out,
    file.path(out_dir, paste0("GCI_distribution_ShapWilk_and_Levene_", gci_var, ".tsv")),
    na = ""
  )
}

ref_cnv_vars <- c(
  "Percentage_SOX10pos_Ref_gain",
  "Percentage_SOX10pos_SNCA_gain",
  "Percentage_SOX10neg_Ref_gain",
  "Percentage_SOX10neg_SNCA_gain",
  "Percentage_SOX10pos_Ref_loss",
  "Percentage_SOX10pos_SNCA_loss",
  "Percentage_SOX10neg_Ref_loss",
  "Percentage_SOX10neg_SNCA_loss"
)
ref_cnv_vars <- keep_existing_vars(ref_cnv_data, ref_cnv_vars)

ref_cnv_msa_v_control <- ref_cnv_data %>%
  dplyr::mutate(
    Group_MSA_vs_Control = dplyr::case_when(
      Group %in% c("SND", "OPCA") ~ "MSA",
      TRUE ~ "Control"
    ),
    Group_MSA_vs_Control = factor(Group_MSA_vs_Control, levels = c("Control", "MSA"))
  )

for (v in ref_cnv_vars) {
  sh <- shapiro_test_function(ref_cnv_msa_v_control, v, c("Group_MSA_vs_Control")) %>%
    dplyr::mutate(Comparison = NA_character_)
  lv <- levene_test_function(ref_cnv_msa_v_control, v, c("Group_MSA_vs_Control"))
  
  write_assumptions_per_variable(
    domain_prefix = "FISH_reference_probe_CNVs_MSA_vs_Control",
    var_name = v,
    shapiro_df = sh,
    levene_df = lv,
    out_dir = out_dir
  )
}

message("Normality + Levene outputs written to: ", out_dir)
