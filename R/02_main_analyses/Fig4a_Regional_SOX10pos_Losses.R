source(here("R", "00_setup", "01_setup_environment.R"))

plot_dir <- plots_reg_sox10pos_loss
stats_dir <- stats_reg_sox10pos_loss

x_order <- c("OPCA", "SND", "Control")

cnv_clean <- cnv_data %>%
  dplyr::filter(Classification %in% x_order) %>%
  dplyr::mutate(
    Classification = factor(Classification, levels = x_order),
    Region = factor(Region),
    DonorRegion = paste(Donor_ID, Region, sep = "__")
  ) %>%
  dplyr::filter(!is.na(Percentage_SOX10pos_loss))

GLOBAL_Y_MAX <- 35

format_p <- function(p_value) {
  p_value <- as.numeric(p_value)
  if (is.na(p_value)) return("p = NA")
  star <- if (p_value < 0.05) "*" else ""
  if (p_value < 0.0001) {
    paste0("p < 0.0001", star)
  } else if (p_value <= 0.001) {
    paste0("p = 0.001", star)
  } else if (p_value == 1) {
    "p = 1"
  } else if (p_value >= 0.045 & p_value < 0.05) {
    paste0("p = ", formatC(p_value, format = "f", digits = 3), star)
  } else {
    rounded <- round(p_value, 2)
    if (rounded == 0) paste0("p = 0.001", star)
    else paste0("p = ", rounded, star)
  }
}

make_brackets_from_pairs <- function(pair_df, x_levels, y_positions) {
  x_map <- setNames(seq_along(x_levels), x_levels)
  pair_df %>%
    dplyr::mutate(
      xmin = x_map[as.character(group1)],
      xmax = x_map[as.character(group2)],
      y.position = y_positions,
      label = vapply(p_adj, format_p, character(1)),
      group_type = paste0(group1, "_vs_", group2)
    )
}

region_tests_across_groups <- function(df_region) {
  kw <- stats::kruskal.test(Percentage_SOX10pos_loss ~ Classification, data = df_region)
  
  dunn <- dunn.test::dunn.test(
    x = df_region$Percentage_SOX10pos_loss,
    g = df_region$Classification,
    method = "bonferroni",
    kw = FALSE,
    list = TRUE
  )
  
  list(
    kw = tibble::tibble(
      Region = as.character(unique(df_region$Region)),
      test = "Kruskal-Wallis",
      comparison = "Classification (overall)",
      statistic = as.numeric(kw$statistic),
      df = as.numeric(kw$parameter),
      p_raw = as.numeric(kw$p.value),
      p_adj = NA_real_
    ),
    dunn = tibble::tibble(
      Region = as.character(unique(df_region$Region)),
      test = "Dunn (Bonferroni)",
      comparison = dunn$comparisons,
      z = as.numeric(dunn$Z),
      p_raw = as.numeric(dunn$P),
      p_adj = as.numeric(dunn$P.adjusted)
    )
  )
}

parse_dunn_pairs <- function(dunn_tbl) {
  dunn_tbl %>%
    dplyr::mutate(
      group1 = sub(" - .*", "", comparison),
      group2 = sub(".* - ", "", comparison)
    )
}

plot_region_across_groups <- function(df_region, brackets, region_label) {
  ggplot2::ggplot(
    df_region,
    ggplot2::aes(x = Classification, y = Percentage_SOX10pos_loss, fill = Classification)
  ) +
    ggplot2::annotate(
      "rect",
      xmin = 0.5, xmax = length(x_order) + 0.5,
      ymin = GLOBAL_Y_MAX + 0.5, ymax = GLOBAL_Y_MAX + 4,
      fill = "white", color = "black", linewidth = 1
    ) +
    ggplot2::annotate(
      "text",
      x = (length(x_order) + 1) / 2,
      y = GLOBAL_Y_MAX + 2.25,
      label = region_label, size = 8, fontface = "bold"
    ) +
    ggplot2::geom_boxplot(color = "black", outlier.shape = NA, alpha = 0.3, linewidth = 0.85) +
    ggplot2::geom_jitter(ggplot2::aes(color = Classification), width = 0.1, size = 3, alpha = 0.8) +
    ggpubr::geom_bracket(
      data = brackets,
      ggplot2::aes(xmin = xmin, xmax = xmax, y.position = y.position, label = label),
      inherit.aes = FALSE,
      label.size = 5,
      size = 0.7,
      tip.length = 0.01
    ) +
    ggplot2::labs(x = "", y = expression("% SOX10+" ~ italic(SNCA) ~ "Losses")) +
    ggplot2::scale_x_discrete(limits = x_order) +
    ggplot2::scale_fill_manual(values = c(
      "Control" = "deepskyblue3",
      "OPCA" = "aquamarine4",
      "SND" = "darkorchid4"
    )) +
    ggplot2::scale_color_manual(values = c(
      "Control" = "deepskyblue3",
      "OPCA" = "aquamarine4",
      "SND" = "darkorchid4"
    )) +
    ggplot2::scale_y_continuous(
      limits = c(-1, GLOBAL_Y_MAX + 4.5),
      expand = ggplot2::expansion(mult = c(0.05, 0))
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 18),
      axis.text = ggplot2::element_text(size = 16),
      legend.position = "none",
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      plot.background = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
}

regions_to_plot <- c("Cerebellum", "Putamen", "SN")
regions_present <- intersect(as.character(unique(cnv_clean$Region)), regions_to_plot)

all_kw_region <- list()
all_dunn_region <- list()
plots_region <- list()

for (reg in regions_present) {
  df_reg <- cnv_clean %>% dplyr::filter(as.character(Region) == reg)
  
  tst <- region_tests_across_groups(df_reg)
  all_kw_region[[reg]] <- tst$kw
  all_dunn_region[[reg]] <- tst$dunn
  
  dunn_pairs <- parse_dunn_pairs(tst$dunn)
  
  desired <- tibble::tibble(
    group1 = c("OPCA", "OPCA", "SND"),
    group2 = c("SND", "Control", "Control")
  )
  
  dunn_pairs2 <- desired %>%
    dplyr::left_join(dunn_pairs %>% dplyr::select(group1, group2, p_adj), by = c("group1", "group2")) %>%
    dplyr::mutate(
      p_adj = ifelse(
        is.na(p_adj),
        dunn_pairs$p_adj[match(
          paste(group2, group1),
          paste(dunn_pairs$group1, dunn_pairs$group2)
        )],
        p_adj
      )
    )
  
  brackets <- make_brackets_from_pairs(dunn_pairs2, x_order, c(32, 30, 28))
  
  p <- plot_region_across_groups(df_reg, brackets, reg)
  plots_region[[reg]] <- p
  
  ggplot2::ggsave(
    file.path(plot_dir, paste0("SOX10pos_Losses_", reg, "_AcrossGroups.png")),
    p, width = 7, height = 7, dpi = 600
  )
}

combined_plot <- plots_region[["Cerebellum"]] +
  plots_region[["Putamen"]] +
  plots_region[["SN"]] +
  patchwork::plot_layout(ncol = 3)

ggplot2::ggsave(
  file.path(plot_dir, "SOX10pos_Losses_AcrossGroups_Combined.png"),
  combined_plot, width = 14, height = 7, dpi = 600
)

combined_plot

readr::write_tsv(
  dplyr::bind_rows(all_kw_region),
  file.path(stats_dir, "SOX10pos_Losses_AcrossGroupsWithinRegion_Kruskal.tsv")
)

readr::write_tsv(
  dplyr::bind_rows(all_dunn_region),
  file.path(stats_dir, "SOX10pos_Losses_AcrossGroupsWithinRegion_Dunn.tsv")
)

within_group_across_regions_kw_dunn <- function(df, cls) {
  df_cls <- df %>% dplyr::filter(Classification == cls)
  if (dplyr::n_distinct(df_cls$Region) < 2) return(list(kw = tibble::tibble(), dunn = tibble::tibble()))
  
  kw <- stats::kruskal.test(Percentage_SOX10pos_loss ~ Region, data = df_cls)
  
  dunn <- dunn.test::dunn.test(
    x = df_cls$Percentage_SOX10pos_loss,
    g = df_cls$Region,
    method = "bonferroni",
    kw = FALSE,
    list = TRUE
  )
  
  kw_tbl <- tibble::tibble(
    Classification = cls,
    test = "Kruskal-Wallis",
    comparison = "Region (overall)",
    statistic = as.numeric(kw$statistic),
    df = as.numeric(kw$parameter),
    p_raw = as.numeric(kw$p.value),
    p_adj = NA_real_
  )
  
  dunn_tbl <- tibble::tibble(
    Classification = cls,
    test = "Dunn (Bonferroni)",
    comparison = dunn$comparisons,
    z = as.numeric(dunn$Z),
    p_raw = as.numeric(dunn$P),
    p_adj = as.numeric(dunn$P.adjusted)
  )
  
  list(kw = kw_tbl, dunn = dunn_tbl)
}

opca_res <- within_group_across_regions_kw_dunn(cnv_clean, "OPCA")
snd_res <- within_group_across_regions_kw_dunn(cnv_clean, "SND")
ctl_res <- within_group_across_regions_kw_dunn(cnv_clean, "Control")

readr::write_tsv(opca_res$kw, file.path(stats_dir, "SOX10pos_Losses_OPCA_Regions_Kruskal.tsv"))
readr::write_tsv(opca_res$dunn, file.path(stats_dir, "SOX10pos_Losses_OPCA_Regions_Dunn.tsv"))

readr::write_tsv(snd_res$kw, file.path(stats_dir, "SOX10pos_Losses_SND_Regions_Kruskal.tsv"))
readr::write_tsv(snd_res$dunn, file.path(stats_dir, "SOX10pos_Losses_SND_Regions_Dunn.tsv"))

readr::write_tsv(ctl_res$kw, file.path(stats_dir, "SOX10pos_Losses_Control_Regions_Kruskal.tsv"))
readr::write_tsv(ctl_res$dunn, file.path(stats_dir, "SOX10pos_Losses_Control_Regions_Dunn.tsv"))

summary_tbl <- cnv_clean %>%
  dplyr::group_by(Classification, Region) %>%
  dplyr::summarise(
    n_obs = dplyr::n(),
    n_donors = dplyr::n_distinct(Donor_ID),
    n_donor_region = dplyr::n_distinct(DonorRegion),
    mean = mean(Percentage_SOX10pos_loss, na.rm = TRUE),
    sd = stats::sd(Percentage_SOX10pos_loss, na.rm = TRUE),
    median = stats::median(Percentage_SOX10pos_loss, na.rm = TRUE),
    iqr_low = stats::quantile(Percentage_SOX10pos_loss, 0.25, na.rm = TRUE),
    iqr_high = stats::quantile(Percentage_SOX10pos_loss, 0.75, na.rm = TRUE),
    min = min(Percentage_SOX10pos_loss, na.rm = TRUE),
    max = max(Percentage_SOX10pos_loss, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_tsv(
  summary_tbl,
  file.path(stats_dir, "SOX10pos_Losses_ByClassification_ByRegion_Summary.tsv")
)
