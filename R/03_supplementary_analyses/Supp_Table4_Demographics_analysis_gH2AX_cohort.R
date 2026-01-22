# GH2AX_demographics_from_case_info.R
# Assumes you will source your project environment + paths
source(here("R", "00_setup", "01_setup_environment.R"))

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

safe_round <- function(x, digits = 1) {
  ifelse(is.na(x), NA_real_, round(x, digits))
}

# Output directory (from setup paths)
dir.create(stats_gh2ax_dir, recursive = TRUE, showWarnings = FALSE)

# 1) Identify unique GH2AX donors (avoid repeated regional measurements)
gh2ax_donors <- gh2ax_data %>%
  dplyr::mutate(Donor = trimws(Donor)) %>%
  dplyr::distinct(Donor) %>%
  dplyr::rename(Donor_ID = Donor)

# 2) Subset case_info to GH2AX cohort only (and recode OPCA/SND -> MSA)
gh2ax_case_info <- case_info %>%
  dplyr::mutate(
    Donor_ID = trimws(Donor_ID),
    Sex = toupper(Sex),
    Classification = dplyr::if_else(Classification %in% c("OPCA", "SND"), "MSA", Classification),
    Classification = factor(Classification, levels = c("Control", "MSA"))
  ) %>%
  dplyr::semi_join(gh2ax_donors, by = "Donor_ID")

# 3) MAIN GH2AX DEMOGRAPHICS SUMMARY TABLE
gh2ax_demographics_summary <- gh2ax_case_info %>%
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
    `Age of onset (years, mean ± SD [range])` = paste0(
      safe_round(mean_Age_of_onset, 1), " +/- ", safe_round(sd_Age_of_onset, 1),
      " (", safe_round(min_Age_of_onset, 1), "-", safe_round(max_Age_of_onset, 1), ")"
    ),
    `Age of onset median [IQR]` = paste0(
      safe_round(median_Age_of_onset, 1), " [",
      safe_round(iqr_low_Age_of_onset, 1), "-", safe_round(iqr_high_Age_of_onset, 1), "]"
    ),
    
    `Disease duration (years, mean ± SD [range])` = paste0(
      safe_round(mean_Disease_duration, 1), " +/- ", safe_round(sd_Disease_duration, 1),
      " (", safe_round(min_Disease_duration, 1), "-", safe_round(max_Disease_duration, 1), ")"
    ),
    `Disease duration median [IQR]` = paste0(
      safe_round(median_Disease_duration, 1), " [",
      safe_round(iqr_low_Disease_duration, 1), "-", safe_round(iqr_high_Disease_duration, 1), "]"
    ),
    
    `Age at death (years, mean ± SD [range])` = paste0(
      safe_round(mean_Age_of_death, 1), " +/- ", safe_round(sd_Age_of_death, 1),
      " (", safe_round(min_Age_of_death, 1), "-", safe_round(max_Age_of_death, 1), ")"
    ),
    `Age at death median [IQR]` = paste0(
      safe_round(median_Age_of_death, 1), " [",
      safe_round(iqr_low_Age_of_death, 1), "-", safe_round(iqr_high_Age_of_death, 1), "]"
    ),
    
    `PMI (hrs, mean ± SD [range])` = paste0(
      safe_round(mean_PMI, 1), " +/- ", safe_round(sd_PMI, 1),
      " (", safe_round(min_PMI, 1), "-", safe_round(max_PMI, 1), ")"
    ),
    `PMI median [IQR]` = paste0(
      safe_round(median_PMI, 1), " [",
      safe_round(iqr_low_PMI, 1), "-", safe_round(iqr_high_PMI, 1), "]"
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
    `Age of onset (years, mean ± SD [range])`,
    `Age of onset median [IQR]`,
    `Disease duration (years, mean ± SD [range])`,
    `Disease duration median [IQR]`,
    `Age at death (years, mean ± SD [range])`,
    `Age at death median [IQR]`,
    `PMI (hrs, mean ± SD [range])`,
    `PMI median [IQR]`,
    `Sex (M/F)`
  ) %>%
  dplyr::arrange(Classification)

readr::write_tsv(
  gh2ax_demographics_summary,
  file.path(stats_gh2ax_dir, "stats_GH2AX_demographics_summary_table.tsv"),
  na = ""
)

# 4) SEX DISTRIBUTION TEST (MSA vs CONTROL)
sex_table_gh2ax <- table(
  gh2ax_case_info$Classification,
  gh2ax_case_info$Sex
)

readr::write_tsv(
  as.data.frame.matrix(sex_table_gh2ax),
  file.path(stats_gh2ax_dir, "sex_distribution_contingency_table_GH2AX.tsv"),
  na = ""
)

chisq_raw <- suppressWarnings(stats::chisq.test(sex_table_gh2ax, correct = FALSE))
expected_counts <- chisq_raw$expected

if (any(expected_counts < 5)) {
  test_used <- "Fisher_exact_test"
  test_result <- stats::fisher.test(sex_table_gh2ax)
} else {
  test_used <- "Chi_square_test"
  test_result <- chisq_raw
}

sex_test_summary_gh2ax <- tibble::tibble(
  Test_used = test_used,
  Statistic = ifelse(test_used == "Chi_square_test", as.numeric(test_result$statistic), NA_real_),
  df = ifelse(test_used == "Chi_square_test", as.numeric(test_result$parameter), NA_real_),
  p_value = as.numeric(test_result$p.value)
)

readr::write_tsv(
  sex_test_summary_gh2ax,
  file.path(stats_gh2ax_dir, "sex_distribution_test_results_GH2AX.tsv"),
  na = ""
)
