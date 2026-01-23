source(here("R", "00_setup", "01_setup_environment.R"))

format_pval <- function(p) {
  stars <- dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
  
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.0001 ~ paste0("p < 0.0001", stars),
    p < 0.001 ~ paste0("p < 0.001", stars),
    TRUE ~ paste0("p = ", signif(p, 2), stars)
  )
}

make_correlation_plot <- function(df, x_var, y_var, group_var, x_label, y_label) {
  df <- df %>%
    dplyr::filter(!is.na(.data[[x_var]]) & !is.na(.data[[y_var]]))
  if (nrow(df) < 3) stop("Not enough observations")
  
  cor_res <- stats::cor.test(df[[x_var]], df[[y_var]], method = "spearman")
  label <- paste0(
    "Overall: rho = ", round(as.numeric(cor_res$estimate), 2),
    ", ", format_pval(cor_res$p.value)
  )
  
  x_min <- floor(min(df[[x_var]]))
  x_max <- ceiling(max(df[[x_var]]))
  y_min <- floor(min(df[[y_var]]))
  y_max <- ceiling(max(df[[y_var]]))
  annot_y <- y_max - (y_max - y_min) * 0.1
  
  df[[group_var]] <- factor(df[[group_var]], levels = c("OPCA", "SND"))
  
  color_mapping <- c(
    "OPCA" = "aquamarine4",
    "SND"  = "darkorchid4"
  )
  
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[x_var]],
      y = .data[[y_var]],
      color = .data[[group_var]]
    )
  ) +
    ggplot2::geom_point(size = 5, alpha = 0.9) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "black")
  
  for (g in levels(df[[group_var]])) {
    df_g <- df %>% dplyr::filter(.data[[group_var]] == g)
    if (nrow(df_g) >= 2) {
      p <- p + ggplot2::geom_smooth(
        data = df_g,
        ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]]),
        method = "lm",
        se = FALSE,
        linetype = "dashed",
        color = color_mapping[g]
      )
    }
  }
  
  p +
    ggplot2::annotate(
      "text",
      x = (x_min + x_max) / 2,
      y = annot_y,
      label = label,
      size = 6,
      fontface = "italic"
    ) +
    ggplot2::scale_color_manual(values = color_mapping) +
    ggplot2::scale_x_continuous(limits = c(x_min, x_max)) +
    ggplot2::scale_y_continuous(limits = c(y_min, y_max)) +
    ggplot2::labs(x = x_label, y = y_label, color = "Group") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      plot.margin = ggplot2::margin(5, 15, 5, 15)
    )
}

msa_regions <- msa_only_data %>%
  dplyr::filter(
    Classification %in% c("OPCA", "SND"),
    Region %in% c("Cerebellum", "Putamen")
  ) %>%
  tidyr::pivot_wider(
    id_cols = c(Donor_ID, Classification),
    names_from = Region,
    values_from = Percentage_SOX10pos_gain
  ) %>%
  dplyr::mutate(
    most_affected  = dplyr::if_else(Classification == "OPCA", Cerebellum, Putamen),
    least_affected = dplyr::if_else(Classification == "OPCA", Putamen, Cerebellum)
  )

age_onset_pos_df <- avg_cnv_data_msa_only %>%
  dplyr::filter(
    !is.na(Age_of_onset),
    !is.na(Percentage_Average_SOX10pos_gain)
  ) %>%
  dplyr::mutate(group = Classification)

age_onset_neg_df <- avg_cnv_data_msa_only %>%
  dplyr::filter(
    !is.na(Age_of_onset),
    !is.na(Percentage_Average_SOX10neg_gain)
  ) %>%
  dplyr::mutate(group = Classification)

msa_gain_df <- msa_only_data %>%
  dplyr::filter(
    !is.na(Percentage_SOX10neg_gain),
    !is.na(Percentage_SOX10pos_gain)
  ) %>%
  dplyr::mutate(group = Classification)

plot1 <- make_correlation_plot(
  msa_regions,
  x_var = "most_affected",
  y_var = "least_affected",
  group_var = "Classification",
  x_label = expression("Most affected region (% SOX10+" ~ italic(SNCA) ~ " gains)"),
  y_label = expression("Lesser affected region (% SOX10+" ~ italic(SNCA) ~ " gains)")
)



plot2 <- make_correlation_plot(
  age_onset_pos_df,
  x_var = "Age_of_onset",
  y_var = "Percentage_Average_SOX10pos_gain",
  group_var = "group",
  x_label = "Age of onset (years)",
  y_label = expression("Average % SOX10+" ~ italic(SNCA) ~ " gains")
)


plot3 <- make_correlation_plot(
  msa_gain_df,
  x_var = "Percentage_SOX10neg_gain",
  y_var = "Percentage_SOX10pos_gain",
  group_var = "group",
  x_label = expression("% SOX10-" ~ italic(SNCA) ~ " gains"),
  y_label = expression("% SOX10+" ~ italic(SNCA) ~ " gains")
)


plot4 <- make_correlation_plot(
  age_onset_neg_df,
  x_var = "Age_of_onset",
  y_var = "Percentage_Average_SOX10neg_gain",
  group_var = "group",
  x_label = "Age of onset (years)",
  y_label = expression("Average % SOX10-" ~ italic(SNCA) ~ " gains")
)


combined_plot <- (plot1 + plot2) / (plot3 + plot4)

cor_plot_dir <- plots_spearman_dir


cor_plot_dir <- plots_spearman_dir

ggplot2::ggsave(
  file.path(cor_plot_dir, "FISH_Spearman_Correlation_Most_vs_Least_Affected_SOX10pos_SNCA_gain.png"),
  plot1,
  width = 6,
  height = 6,
  dpi = 600
)

ggplot2::ggsave(
  file.path(cor_plot_dir, "FISH_Spearman_Correlation_Avg_SOX10pos_SNCA_gain_vs_Age_of_onset.png"),
  plot2,
  width = 6,
  height = 6,
  dpi = 600
)

ggplot2::ggsave(
  file.path(cor_plot_dir, "FISH_Spearman_Correlation_SOX10neg_vs_SOX10pos_SNCA_gain.png"),
  plot3,
  width = 6,
  height = 6,
  dpi = 600
)

ggplot2::ggsave(
  file.path(cor_plot_dir, "FISH_Spearman_Correlation_Avg_SOX10neg_SNCA_gain_vs_Age_of_onset.png"),
  plot4,
  width = 6,
  height = 6,
  dpi = 600
)

ggplot2::ggsave(
  file.path(cor_plot_dir, "FISH_Spearman_Correlation_SOX10_SNCA_gain_Combined_2x2.png"),
  combined_plot,
  width = 12,
  height = 12,
  dpi = 600
)
