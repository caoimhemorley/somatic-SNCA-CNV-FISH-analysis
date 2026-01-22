source(here("R", "00_setup", "01_setup_environment.R"))

ggplot2::theme_set(ggplot2::theme_minimal(base_family = "Helvetica"))

add_consistent_pval <- function(plot_obj, data, y_var, p_val_text, group_x = c(1, 2)) {
  max_y <- max(data[[y_var]], na.rm = TRUE)
  min_y <- min(data[[y_var]], na.rm = TRUE)
  y_range <- max_y - min_y
  
  if (!is.finite(y_range) || y_range == 0) {
    y_range <- max_y * 0.1
    if (!is.finite(y_range) || y_range == 0) y_range <- 1
  }
  
  y_top <- max_y + y_range * 0.15
  bracket_y <- max_y + y_range * 0.10
  tick_length <- y_range * 0.03
  
  plot_obj +
    ggplot2::expand_limits(y = y_top + tick_length) +
    ggplot2::annotate("text", x = mean(group_x), y = y_top, label = p_val_text, size = 6) +
    ggplot2::annotate(
      "segment",
      x = group_x[1], xend = group_x[2],
      y = bracket_y, yend = bracket_y,
      linewidth = 0.7
    ) +
    ggplot2::annotate(
      "segment",
      x = group_x[1], xend = group_x[1],
      y = bracket_y, yend = bracket_y - tick_length,
      linewidth = 0.7
    ) +
    ggplot2::annotate(
      "segment",
      x = group_x[2], xend = group_x[2],
      y = bracket_y, yend = bracket_y - tick_length,
      linewidth = 0.7
    )
}

get_summary <- function(data, value_col, group_col, comparison) {
  data %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(
      n = sum(!is.na(.data[[value_col]])),
      mean = mean(.data[[value_col]], na.rm = TRUE),
      sd = stats::sd(.data[[value_col]], na.rm = TRUE),
      median = stats::median(.data[[value_col]], na.rm = TRUE),
      q1 = stats::quantile(.data[[value_col]], 0.25, na.rm = TRUE),
      q3 = stats::quantile(.data[[value_col]], 0.75, na.rm = TRUE),
      iqr = stats::IQR(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      variable = value_col,
      comparison = comparison
    ) %>%
    dplyr::relocate(comparison, variable)
}

asyn_long <- aSyn_NS_data %>%
  dplyr::select(
    Donor_ID, Region, Subtype,
    Percent_inclusions_sections,
    Percent_inclusions_NS
  ) %>%
  tidyr::pivot_longer(
    cols = c(Percent_inclusions_sections, Percent_inclusions_NS),
    names_to = "Tissue",
    values_to = "Percent_inclusions"
  ) %>%
  dplyr::mutate(
    Tissue = dplyr::case_when(
      Tissue == "Percent_inclusions_sections" ~ "Sections",
      Tissue == "Percent_inclusions_NS" ~ "NS",
      TRUE ~ NA_character_
    ),
    Tissue = factor(Tissue, levels = c("Sections", "NS")),
    Pair_ID = interaction(Donor_ID, Region, drop = TRUE)
  ) %>%
  dplyr::filter(!is.na(Tissue))

paired_complete <- asyn_long %>%
  dplyr::group_by(Pair_ID) %>%
  dplyr::summarise(
    n_tissues = dplyr::n_distinct(Tissue[!is.na(Percent_inclusions)]),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_tissues == 2) %>%
  dplyr::pull(Pair_ID)

asyn_long <- asyn_long %>%
  dplyr::filter(Pair_ID %in% paired_complete) %>%
  tidyr::drop_na(Percent_inclusions)

asyn_wide <- asyn_long %>%
  dplyr::select(Donor_ID, Region, Subtype, Pair_ID, Tissue, Percent_inclusions) %>%
  tidyr::pivot_wider(
    names_from = Tissue,
    values_from = Percent_inclusions
  ) %>%
  tidyr::drop_na(Sections, NS)

t_res <- stats::t.test(asyn_wide$Sections, asyn_wide$NS, paired = TRUE)

p_text <- dplyr::if_else(
  is.na(t_res$p.value),
  "p = NA",
  dplyr::if_else(
    t_res$p.value < 0.05,
    sprintf("p = %.3f*", t_res$p.value),
    sprintf("p = %.3f", t_res$p.value)
  )
)

colors_box_asyn <- c(
  "Sections" = scales::alpha("skyblue", 0.5),
  "NS" = scales::alpha("aquamarine4", 0.5)
)

colors_points_asyn <- c(
  "Sections" = "skyblue",
  "NS" = "aquamarine4"
)

p_asyn <- ggplot2::ggplot(
  asyn_long,
  ggplot2::aes(x = Tissue, y = Percent_inclusions, fill = Tissue)
) +
  ggplot2::geom_boxplot(
    alpha = 0.3,
    outlier.shape = NA,
    width = 0.4,
    color = "black",
    linewidth = 0.7
  ) +
  ggplot2::geom_line(
    ggplot2::aes(group = Pair_ID),
    color = "gray40",
    linewidth = 0.6,
    alpha = 0.7
  ) +
  ggplot2::geom_jitter(
    ggplot2::aes(color = Tissue),
    width = 0.05,
    size = 3
  ) +
  ggplot2::scale_fill_manual(values = colors_box_asyn) +
  ggplot2::scale_color_manual(values = colors_points_asyn) +
  ggplot2::labs(
    title = "αSyn inclusions: Sections vs NS (paired)",
    x = NULL,
    y = "% cells with inclusions"
  ) +
  ggplot2::theme(
    legend.position = "none",
    plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title.y = ggplot2::element_text(size = 18),
    axis.text = ggplot2::element_text(size = 14),
    panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
  )

p_asyn <- add_consistent_pval(
  p_asyn,
  asyn_long,
  "Percent_inclusions",
  p_text,
  group_x = c(1, 2)
)

summary_asyn <- get_summary(
  asyn_long,
  "Percent_inclusions",
  "Tissue",
  "aSyn inclusions: Sections vs NS"
) %>%
  dplyr::mutate(dplyr::across(c(mean, sd, median, q1, q3, iqr), ~ round(.x, 2)))

ttest_table <- tibble::tibble(
  comparison = "aSyn inclusions: Sections vs NS",
  test = "paired t-test",
  n_pairs = nrow(asyn_wide),
  t_statistic = as.numeric(t_res$statistic),
  df = as.numeric(t_res$parameter),
  p_value = t_res$p.value,
  mean_sections = mean(asyn_wide$Sections, na.rm = TRUE),
  mean_ns = mean(asyn_wide$NS, na.rm = TRUE),
  mean_diff_sections_minus_ns = mean(asyn_wide$Sections - asyn_wide$NS, na.rm = TRUE),
  conf_low = t_res$conf.int[1],
  conf_high = t_res$conf.int[2]
) %>%
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ signif(.x, 4)))

readr::write_tsv(
  summary_asyn,
  file.path(stats_asyn_ns_dir, "aSyn_NS_vs_Sections_SummaryStats_Percent_inclusions.tsv"),
  na = ""
)

readr::write_tsv(
  ttest_table,
  file.path(stats_asyn_ns_dir, "aSyn_NS_vs_Sections_PairedTTest_Percent_inclusions.tsv"),
  na = ""
)

readr::write_tsv(
  asyn_wide,
  file.path(stats_asyn_ns_dir, "aSyn_NS_vs_Sections_PairedDataWide_Percent_inclusions.tsv"),
  na = ""
)

ggplot2::ggsave(
  filename = file.path(plots_asyn_ns_dir, "aSyn_NS_vs_Sections_PairedPlot_Percent_inclusions.png"),
  plot = p_asyn,
  width = 5,
  height = 6,
  dpi = 600
)

print(p_asyn)
