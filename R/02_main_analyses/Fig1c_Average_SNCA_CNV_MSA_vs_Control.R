source(here("R", "00_setup", "01_setup_environment.R"))

desired_order <- c(
  "bold(\"SOX10\")^\"+\"~bolditalic(SNCA)~bold(gains)",
  "bold(\"SOX10\")^\"−\"~bolditalic(SNCA)~bold(gains)",
  "bold(\"SOX10\")^\"+\"~bolditalic(SNCA)~bold(losses)",
  "bold(\"SOX10\")^\"−\"~bolditalic(SNCA)~bold(losses)"
)



plot_data_long <- avg_cnv_data %>%
  dplyr::mutate(
    Classification = dplyr::if_else(Classification %in% c("SND", "OPCA"), "MSA", Classification),
    Classification = factor(Classification, levels = c("Control", "MSA"))
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::starts_with("Percentage_Average_"),
    names_to = "raw_name",
    values_to = "Percentage"
  ) %>%
  dplyr::mutate(
    raw_name = trimws(raw_name),
    CNV_Label = dplyr::case_when(
      raw_name == "Percentage_Average_SOX10pos_gain" ~
        "bold(\"SOX10\")^\"+\"~bolditalic(SNCA)~bold(gains)",
      raw_name == "Percentage_Average_SOX10neg_gain" ~
        "bold(\"SOX10\")^\"−\"~bolditalic(SNCA)~bold(gains)",
      raw_name == "Percentage_Average_SOX10pos_loss" ~
        "bold(\"SOX10\")^\"+\"~bolditalic(SNCA)~bold(losses)",
      raw_name == "Percentage_Average_SOX10neg_loss" ~
        "bold(\"SOX10\")^\"−\"~bolditalic(SNCA)~bold(losses)",
      TRUE ~ NA_character_
    )
    
  ) %>%
  dplyr::filter(!is.na(CNV_Label)) %>%
  dplyr::mutate(
    CNV_Label = factor(CNV_Label, levels = desired_order)
  )

mw_results <- purrr::map_dfr(
  desired_order,
  function(lbl) {
    sub <- dplyr::filter(plot_data_long, CNV_Label == lbl)
    test <- stats::wilcox.test(Percentage ~ Classification, data = sub, exact = TRUE)
    pval <- test$p.value
    
    formatted_p <- dplyr::case_when(
      is.na(pval) ~ NA_character_,
      pval < 0.0001 ~ "p < 0.0001",
      pval < 0.001 ~ "p < 0.001",
      pval < 0.01 ~ "p < 0.01",
      TRUE ~ paste0("p = ", format(round(pval, 2), nsmall = 2))
    )
    
    tibble::tibble(
      CNV_Label = lbl,
      Test = "Mann-Whitney",
      Comparison = "Control vs MSA",
      Variable_Set = "Average_CNVs",
      p_value = pval,
      p_text = formatted_p,
      label = dplyr::if_else(!is.na(pval) & pval < 0.05, paste0(formatted_p, "*"), formatted_p)
    )
  }
)

readr::write_tsv(
  mw_results,
  file.path(
    stats_avg_cnv_dir,
    "FISH_average_SNCA_CNVs_MSA_vs_Control_MannWhitney_by_CNV_Label.tsv"
  ),
  na = ""
)

ymax <- max(plot_data_long$Percentage, na.rm = TRUE)
bracket_y <- ymax * 1.18

bracket_positions <- mw_results %>%
  dplyr::mutate(
    xmin = 1,
    xmax = 2,
    y.position = bracket_y,
    CNV_Label = factor(CNV_Label, levels = desired_order)
  )

msa_control_cnv_plot <- ggplot2::ggplot(
  plot_data_long,
  ggplot2::aes(x = Classification, y = Percentage, fill = Classification)
) +
  ggplot2::geom_boxplot(
    outlier.shape = NA,
    color = "black",
    alpha = 0.3,
    width = 0.6,
    linewidth = 0.4
  ) +
  ggplot2::geom_jitter(
    ggplot2::aes(color = Classification),
    size = 1.5,
    alpha = 0.8,
    position = ggplot2::position_jitter(width = 0.15)
  ) +
  ggplot2::scale_fill_manual(values = c("Control" = "deepskyblue4", "MSA" = "slateblue4")) +
  ggplot2::scale_color_manual(values = c("Control" = "deepskyblue4", "MSA" = "slateblue4")) +
  ggplot2::facet_wrap(
    ~ CNV_Label,
    nrow = 1,
    labeller = ggplot2::label_parsed
  ) +
  ggpubr::geom_bracket(
    data = bracket_positions,
    ggplot2::aes(xmin = xmin, xmax = xmax, y.position = y.position, label = label),
    inherit.aes = FALSE,
    size = 0.7,
    tip.length = 0.02,
    label.size = 5
  ) +
  ggplot2::labs(x = "", y = "% of cells") +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(size = 14, face = "bold"),
    axis.text.x = ggplot2::element_text(size = 12),
    panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "none"
  ) +
  ggplot2::expand_limits(y = bracket_y * 1.05)

ggplot2::ggsave(
  filename = file.path(
    plots_avg_cnv_dir,
    "FISH_average_SNCA_CNVs_MSA_vs_Control_Average_CNVs_faceted_boxplot_by_CNV_Label.png"
  ),
  plot = msa_control_cnv_plot,
  width = 12,
  height = 5,
  dpi = 1200
)

iqr_bounds <- function(x) {
  q <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  list(low = q[1], high = q[2])
}

cnv_summary_table <- plot_data_long %>%
  dplyr::group_by(CNV_Label, Classification) %>%
  dplyr::summarise(
    n = sum(!is.na(Percentage)),
    mean = mean(Percentage, na.rm = TRUE),
    sd = stats::sd(Percentage, na.rm = TRUE),
    min = min(Percentage, na.rm = TRUE),
    max = max(Percentage, na.rm = TRUE),
    median = stats::median(Percentage, na.rm = TRUE),
    iqr_low = iqr_bounds(Percentage)$low,
    iqr_high = iqr_bounds(Percentage)$high,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Mean_SD_range = paste0(
      round(mean, 2), " +/- ", round(sd, 2),
      " (", round(min, 2), "-", round(max, 2), ")"
    ),
    Median_IQR = paste0(
      round(median, 2), " [",
      round(iqr_low, 2), "-", round(iqr_high, 2), "]"
    )
  ) %>%
  dplyr::select(
    CNV_Label, Classification, n,
    Mean_SD_range, Median_IQR
  ) %>%
  dplyr::arrange(match(CNV_Label, desired_order), Classification)

readr::write_tsv(
  cnv_summary_table,
  file.path(
    stats_avg_cnv_dir,
    "FISH_average_SNCA_CNVs_MSA_vs_Control_Average_CNVs_SummaryTable_MeanSDRange_and_MedianIQR.tsv"
  ),
  na = ""
)
