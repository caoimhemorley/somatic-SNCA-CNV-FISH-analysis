# Helper function for IQR
iqr_bounds <- function(x) {
  q <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  list(low = q[1], high = q[2])
}


safe_sum_logical <- function(x) {
  if (all(is.na(x))) return(NA_integer_)
  sum(x, na.rm = TRUE)
}

safe_min <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}


# Main demographics summary table
demographics_summary_table <- case_info %>%
  dplyr::mutate(Sex = toupper(Sex)) %>%
  dplyr::group_by(Classification) %>%
  dplyr::summarise(
    n = dplyr::n(),
    
    mean_Age_of_onset = mean(Age_of_onset, na.rm = TRUE),
    sd_Age_of_onset = stats::sd(Age_of_onset, na.rm = TRUE),
    min_Age_of_onset = safe_min(Age_of_onset),
    max_Age_of_onset = safe_max(Age_of_onset),
    median_Age_of_onset = stats::median(Age_of_onset, na.rm = TRUE),
    iqr_low_Age_of_onset = iqr_bounds(Age_of_onset)$low,
    iqr_high_Age_of_onset = iqr_bounds(Age_of_onset)$high,
    
    mean_Disease_duration = mean(Disease_duration, na.rm = TRUE),
    sd_Disease_duration = stats::sd(Disease_duration, na.rm = TRUE),
    min_Disease_duration = safe_min(Disease_duration),
    max_Disease_duration = safe_max(Disease_duration),
    median_Disease_duration = stats::median(Disease_duration, na.rm = TRUE),
    iqr_low_Disease_duration = iqr_bounds(Disease_duration)$low,
    iqr_high_Disease_duration = iqr_bounds(Disease_duration)$high,
    
    mean_Age_of_death = mean(Age_of_death, na.rm = TRUE),
    sd_Age_of_death = stats::sd(Age_of_death, na.rm = TRUE),
    min_Age_of_death = safe_min(Age_of_death),
    max_Age_of_death = safe_max(Age_of_death),
    median_Age_of_death = stats::median(Age_of_death, na.rm = TRUE),
    iqr_low_Age_of_death = iqr_bounds(Age_of_death)$low,
    iqr_high_Age_of_death = iqr_bounds(Age_of_death)$high,
    
    mean_PMI = mean(PMI, na.rm = TRUE),
    sd_PMI = stats::sd(PMI, na.rm = TRUE),
    min_PMI = safe_min(PMI),
    max_PMI = safe_max(PMI),
    median_PMI = stats::median(PMI, na.rm = TRUE),
    iqr_low_PMI = iqr_bounds(PMI)$low,
    iqr_high_PMI = iqr_bounds(PMI)$high,
    
    males = safe_sum_logical(Sex == "M"),
    females = safe_sum_logical(Sex == "F"),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    `Age of onset (years)` = paste0(
      round(mean_Age_of_onset, 1), " +/- ", round(sd_Age_of_onset, 1),
      " (", round(min_Age_of_onset, 1), "-", round(max_Age_of_onset, 1), ")"
    ),
    `Age of onset median [IQR]` = paste0(
      round(median_Age_of_onset, 1), " [",
      round(iqr_low_Age_of_onset, 1), "-", round(iqr_high_Age_of_onset, 1), "]"
    ),
    
    `Disease duration (years)` = paste0(
      round(mean_Disease_duration, 1), " +/- ", round(sd_Disease_duration, 1),
      " (", round(min_Disease_duration, 1), "-", round(max_Disease_duration, 1), ")"
    ),
    `Disease duration median [IQR]` = paste0(
      round(median_Disease_duration, 1), " [",
      round(iqr_low_Disease_duration, 1), "-", round(iqr_high_Disease_duration, 1), "]"
    ),
    
    `Age at death (years)` = paste0(
      round(mean_Age_of_death, 1), " +/- ", round(sd_Age_of_death, 1),
      " (", round(min_Age_of_death, 1), "-", round(max_Age_of_death, 1), ")"
    ),
    `Age at death median [IQR]` = paste0(
      round(median_Age_of_death, 1), " [",
      round(iqr_low_Age_of_death, 1), "-", round(iqr_high_Age_of_death, 1), "]"
    ),
    
    `PMI (hrs)` = paste0(
      round(mean_PMI, 1), " +/- ", round(sd_PMI, 1),
      " (", round(min_PMI, 1), "-", round(max_PMI, 1), ")"
    ),
    `PMI median [IQR]` = paste0(
      round(median_PMI, 1), " [",
      round(iqr_low_PMI, 1), "-", round(iqr_high_PMI, 1), "]"
    ),
    
    `Sex (M/F)` = dplyr::case_when(
      is.na(males) | is.na(females) ~ NA_character_,
      TRUE ~ paste0(
        males, " (", round(100 * males / n, 1), "%) / ",
        females, " (", round(100 * females / n, 1), "%)"
      )
    )
  ) %>%
  dplyr::select(
    Classification,
    n,
    `Age of onset (years)`,
    `Age of onset median [IQR]`,
    `Disease duration (years)`,
    `Disease duration median [IQR]`,
    `Age at death (years)`,
    `Age at death median [IQR]`,
    `PMI (hrs)`,
    `PMI median [IQR]`,
    `Sex (M/F)`
  ) %>%
  dplyr::arrange(Classification)

# Export
readr::write_tsv(
  demographics_summary_table,
  file.path(stats_demo_dir, "stats_demographics_full_summary_table.tsv"),
  na = ""
)

# 1. Subset OPCA and SND from averaged CNV dataset
avg_msa_subset <- avg_cnv_data_msa_only %>%
  dplyr::filter(Classification %in% c("OPCA", "SND"))

# 2. Student t-test for Age_of_onset
t_age_onset <- stats::t.test(
  Age_of_onset ~ Classification,
  data = avg_msa_subset,
  var.equal = TRUE
)

# 3. Student t-test for Disease_duration
t_disease_duration <- stats::t.test(
  Disease_duration ~ Classification,
  data = avg_msa_subset,
  var.equal = FALSE
)

# 4. Summary table for both variables
summary_duration_onset <- avg_msa_subset %>%
  dplyr::group_by(Classification) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_Age_of_onset = mean(Age_of_onset, na.rm = TRUE),
    sd_Age_of_onset = stats::sd(Age_of_onset, na.rm = TRUE),
    median_Age_of_onset = stats::median(Age_of_onset, na.rm = TRUE),
    iqr_low_Age_of_onset = stats::quantile(Age_of_onset, 0.25, na.rm = TRUE),
    iqr_high_Age_of_onset = stats::quantile(Age_of_onset, 0.75, na.rm = TRUE),
    
    mean_Disease_duration = mean(Disease_duration, na.rm = TRUE),
    sd_Disease_duration = stats::sd(Disease_duration, na.rm = TRUE),
    median_Disease_duration = stats::median(Disease_duration, na.rm = TRUE),
    iqr_low_Disease_duration = stats::quantile(Disease_duration, 0.25, na.rm = TRUE),
    iqr_high_Disease_duration = stats::quantile(Disease_duration, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_tsv(
  summary_duration_onset,
  file.path(stats_demo_dir, "stats_OPCA_vs_SND_onset_and_duration_summary.tsv"),
  na = ""
)

# 5. Save t-test results as text files
sink(file.path(stats_demo_dir, "t_test_age_of_onset_OPCA_vs_SND.txt"))
print(t_age_onset)
sink()

sink(file.path(stats_demo_dir, "t_test_disease_duration_OPCA_vs_SND.txt"))
print(t_disease_duration)
sink()

# 1. Build 3x2 contingency table (Classification x Sex)
sex_table <- table(case_info$Classification, case_info$Sex)

# Save the contingency table
readr::write_tsv(
  as.data.frame.matrix(sex_table),
  file.path(stats_demo_dir, "sex_distribution_contingency_table.tsv"),
  na = ""
)

# 2. Check expected counts for chi square assumption
chisq_test_raw <- stats::chisq.test(sex_table, correct = FALSE)
expected_counts <- chisq_test_raw$expected

# 3. Decide whether to use chi square or Fisher
if (any(expected_counts < 5)) {
  test_used <- "Fisher_exact_test"
  test_result <- stats::fisher.test(sex_table)
} else {
  test_used <- "Chi_square_test"
  test_result <- chisq_test_raw
}

# 4. Convert the result into a readable table
sex_test_summary <- tibble::tibble(
  Test_used = test_used,
  Statistic = ifelse(test_used == "Chi_square_test", test_result$statistic, NA),
  df = ifelse(test_used == "Chi_square_test", test_result$parameter, NA),
  p_value = test_result$p.value
)

# 5. Save statistical results
readr::write_tsv(
  sex_test_summary,
  file.path(stats_demo_dir, "sex_distribution_test_results.tsv"),
  na = ""
)
