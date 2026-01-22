

source(here::here("R", "00_setup", "01_setup_environment.R"))

iqr_bounds <- function(x) {
  q <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  list(low = q[1], high = q[2])
}

dir.create(stats_overall_cnv_dir, recursive = TRUE, showWarnings = FALSE)

summary_overall_gain_loss_regional <- cnv_data %>%
  dplyr::filter(!is.na(Region), !is.na(Classification)) %>%
  dplyr::group_by(Classification, Region) %>%
  dplyr::summarise(
    n_gain = sum(!is.na(Percentage_Overall_gain)),
    n_loss = sum(!is.na(Percentage_Overall_loss)),
    n = max(n_gain, n_loss),
    
    mean_Overall_gain = mean(Percentage_Overall_gain, na.rm = TRUE),
    sd_Overall_gain = stats::sd(Percentage_Overall_gain, na.rm = TRUE),
    median_Overall_gain = stats::median(Percentage_Overall_gain, na.rm = TRUE),
    q1_Overall_gain = iqr_bounds(Percentage_Overall_gain)$low,
    q3_Overall_gain = iqr_bounds(Percentage_Overall_gain)$high,
    min_Overall_gain = min(Percentage_Overall_gain, na.rm = TRUE),
    max_Overall_gain = max(Percentage_Overall_gain, na.rm = TRUE),
    
    mean_Overall_loss = mean(Percentage_Overall_loss, na.rm = TRUE),
    sd_Overall_loss = stats::sd(Percentage_Overall_loss, na.rm = TRUE),
    median_Overall_loss = stats::median(Percentage_Overall_loss, na.rm = TRUE),
    q1_Overall_loss = iqr_bounds(Percentage_Overall_loss)$low,
    q3_Overall_loss = iqr_bounds(Percentage_Overall_loss)$high,
    min_Overall_loss = min(Percentage_Overall_loss, na.rm = TRUE),
    max_Overall_loss = max(Percentage_Overall_loss, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    `Overall gain mean +/- SD (min-max)` = paste0(
      round(mean_Overall_gain, 2), " +/- ", round(sd_Overall_gain, 2),
      " (", round(min_Overall_gain, 2), "-", round(max_Overall_gain, 2), ")"
    ),
    `Overall gain median [Q1-Q3]` = paste0(
      round(median_Overall_gain, 2), " [",
      round(q1_Overall_gain, 2), "-", round(q3_Overall_gain, 2), "]"
    ),
    `Overall loss mean +/- SD (min-max)` = paste0(
      round(mean_Overall_loss, 2), " +/- ", round(sd_Overall_loss, 2),
      " (", round(min_Overall_loss, 2), "-", round(max_Overall_loss, 2), ")"
    ),
    `Overall loss median [Q1-Q3]` = paste0(
      round(median_Overall_loss, 2), " [",
      round(q1_Overall_loss, 2), "-", round(q3_Overall_loss, 2), "]"
    )
  ) %>%
  dplyr::select(
    Classification,
    Region,
    n,
    `Overall gain mean +/- SD (min-max)`,
    `Overall gain median [Q1-Q3]`,
    `Overall loss mean +/- SD (min-max)`,
    `Overall loss median [Q1-Q3]`
  ) %>%
  dplyr::arrange(Classification, Region)

readr::write_tsv(
  summary_overall_gain_loss_regional,
  file.path(stats_overall_cnv_dir, "summary_overall_gain_loss_regional.tsv"),
  na = ""
)

print(summary_overall_gain_loss_regional, n = Inf)
