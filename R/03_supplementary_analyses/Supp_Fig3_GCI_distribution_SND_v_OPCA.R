# GCI_summary_from_CNV_data.R
# Summarise + stats + plots for % SOX10+ αSyn inclusions (GCIs) from cnv_data
# Plot shows Cerebellum + Putamen only (summaries still include SN)
# Uses paths defined in 00_setup_paths.R + objects created in 01_setup_environment.R

source(here::here("R", "00_setup", "01_setup_environment.R"))

ggplot2::theme_set(ggplot2::theme_minimal(base_family = "Helvetica"))

iqr_bounds <- function(x) {
  q <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  list(low = q[1], high = q[2])
}

format_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.0001) return("p < 0.0001")
  if (p < 0.001) return("p < 0.001")
  if (p < 0.01) return("p < 0.01")
  paste0("p = ", formatC(p, format = "f", digits = 3))
}

gci_var <- "Percentage_SOX10pos_inclusions"
if (!gci_var %in% colnames(cnv_data)) {
  stop(paste0("Column not found in cnv_data: ", gci_var))
}

dir.create(stats_gci_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_gci_dir, showWarnings = FALSE, recursive = TRUE)

# MSA-only view (OPCA + SND)
gci_data_msa <- cnv_data %>%
  dplyr::filter(Classification %in% c("OPCA", "SND")) %>%
  dplyr::mutate(
    Classification = factor(as.character(Classification), levels = c("OPCA", "SND")),
    Region = factor(as.character(Region), levels = c("Cerebellum", "Putamen", "SN"))
  )

# 1) Summary table by Classification x Region (includes SN)
gci_summary <- gci_data_msa %>%
  dplyr::group_by(Classification, Region) %>%
  dplyr::summarise(
    n = sum(!is.na(.data[[gci_var]])),
    mean = mean(.data[[gci_var]], na.rm = TRUE),
    sd = stats::sd(.data[[gci_var]], na.rm = TRUE),
    min = min(.data[[gci_var]], na.rm = TRUE),
    max = max(.data[[gci_var]], na.rm = TRUE),
    median = stats::median(.data[[gci_var]], na.rm = TRUE),
    iqr_low = iqr_bounds(.data[[gci_var]])$low,
    iqr_high = iqr_bounds(.data[[gci_var]])$high,
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
  dplyr::select(Classification, Region, n, Mean_SD_range, Median_IQR) %>%
  dplyr::arrange(Classification, Region)

readr::write_tsv(
  gci_summary,
  file.path(stats_gci_dir, "GCI_SOX10pos_inclusions_summary_by_region.tsv"),
  na = ""
)

print(gci_summary, n = Inf)

# 2) Within-subtype regional comparisons: Putamen vs Cerebellum (unpaired Wilcoxon)
get_region_wilcox <- function(df, subtype_label) {
  sub <- df %>%
    dplyr::filter(Region %in% c("Putamen", "Cerebellum")) %>%
    tidyr::drop_na(.data[[gci_var]])
  
  if (dplyr::n_distinct(sub$Region) < 2) {
    return(tibble::tibble(
      Classification = subtype_label,
      comparison = "Putamen vs Cerebellum",
      test = "Wilcoxon rank-sum",
      p_value = NA_real_,
      p_text = NA_character_,
      note = "Only one region present after filtering"
    ))
  }
  
  wt <- stats::wilcox.test(stats::as.formula(paste0(gci_var, " ~ Region")), data = sub, exact = FALSE)
  
  tibble::tibble(
    Classification = subtype_label,
    comparison = "Putamen vs Cerebellum",
    test = "Wilcoxon rank-sum",
    p_value = wt$p.value,
    p_text = format_p(wt$p.value),
    note = NA_character_
  )
}

wilcox_region_results <- dplyr::bind_rows(
  get_region_wilcox(gci_data_msa %>% dplyr::filter(Classification == "SND"), "SND"),
  get_region_wilcox(gci_data_msa %>% dplyr::filter(Classification == "OPCA"), "OPCA")
)

readr::write_tsv(
  wilcox_region_results,
  file.path(stats_gci_dir, "GCI_SOX10pos_inclusions_Putamen_vs_Cerebellum_Wilcoxon.tsv"),
  na = ""
)

print(wilcox_region_results, n = Inf)

# 3) Bracket annotations for faceted plot
pval_annotations <- wilcox_region_results %>%
  dplyr::mutate(
    xstart = 1,
    xend = 2,
    y.position = 70,
    label = dplyr::if_else(is.na(p_text), "p = NA", p_text)
  ) %>%
  dplyr::select(Classification, xstart, xend, y.position, label)

# 4) Plot: Cerebellum + Putamen only (faceted by OPCA/SND) with p-value brackets
plot_gci_cb_put <- function(data, title_text, y_lim = c(0, 75)) {
  df <- data %>%
    dplyr::filter(Region %in% c("Cerebellum", "Putamen")) %>%
    tidyr::drop_na(.data[[gci_var]]) %>%
    dplyr::mutate(
      Region = factor(as.character(Region), levels = c("Cerebellum", "Putamen")),
      Classification = factor(as.character(Classification), levels = c("OPCA", "SND"))
    )
  
  ggplot2::ggplot(df, ggplot2::aes(x = Region, y = .data[[gci_var]], fill = Classification)) +
    ggplot2::geom_boxplot(alpha = 0.2, outlier.shape = NA, color = "black", width = 0.5, linewidth = 0.6) +
    ggplot2::geom_jitter(ggplot2::aes(color = Classification), size = 2.5, alpha = 0.8, width = 0.1) +
    ggplot2::facet_wrap(~ Classification, scales = "free_x") +
    ggplot2::coord_cartesian(ylim = y_lim) +
    ggplot2::labs(
      title = title_text,
      x = "",
      y = expression("% SOX10+" ~ alpha * "Syn inclusions (GCIs)")
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, size = 14, colour = "black"),
      axis.text.y = ggplot2::element_text(size = 14, colour = "black"),
      axis.title.y = ggplot2::element_text(size = 16, colour = "black"),
      strip.text = ggplot2::element_text(size = 16, face = "bold", colour = "black"),
      plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5, colour = "black"),
      legend.position = "none",
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    ggplot2::scale_fill_manual(values = c("OPCA" = "aquamarine4", "SND" = "darkorchid4")) +
    ggplot2::scale_color_manual(values = c("OPCA" = "aquamarine4", "SND" = "darkorchid4")) +
    ggpubr::geom_bracket(
      data = pval_annotations,
      ggplot2::aes(xmin = xstart, xmax = xend, y.position = y.position, label = label),
      inherit.aes = FALSE,
      size = 0.5,
      label.size = 4.5,
      label.fontface = "bold",
      color = "black",
      label.color = "black",
      tip.length = 0.01
    )
}

p_gci_cb_put <- plot_gci_cb_put(
  gci_data_msa,
  "% SOX10+ αSyn inclusions (GCIs): Putamen vs Cerebellum"
)

ggplot2::ggsave(
  filename = file.path(plots_gci_dir, "GCI_SOX10pos_inclusions_OPCA_SND_Cerebellum_vs_Putamen.png"),
  plot = p_gci_cb_put,
  width = 9,
  height = 7,
  dpi = 600
)

print(p_gci_cb_put)

message("GCI outputs written to:\n  Plots: ", plots_gci_dir, "\n  Stats: ", stats_gci_dir)
