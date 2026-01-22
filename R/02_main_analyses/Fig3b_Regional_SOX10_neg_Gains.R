source(here("R", "00_setup", "01_setup_environment.R"))

plot_dir <- plots_reg_sox10neg_gain
stats_dir <- stats_reg_sox10neg_gain

x_order <- c("OPCA", "SND", "Control")

cnv_clean <- cnv_data %>%
  dplyr::filter(Classification %in% x_order) %>%
  dplyr::mutate(
    Classification = factor(Classification, levels = x_order),
    Region = factor(Region),
    DonorRegion = paste(Donor_ID, Region, sep = "__")
  ) %>%
  dplyr::filter(!is.na(Percentage_SOX10neg_gain))

GLOBAL_Y_MAX <- 35

format_p <- function(p_value) {
  p_value <- as.numeric(p_value)
  
  if (is.na(p_value)) {
    return("p = NA")
  }
  
  star <- if (p_value < 0.05) "*" else ""
  
  if (p_value < 0.0001) {
    return(paste0("p < 0.0001", star))
  } else if (p_value <= 0.001) {
    return(paste0("p = 0.001", star))
  } else if (p_value == 1) {
    return("p = 1")
  } else if (p_value >= 0.045 & p_value < 0.05) {
    return(paste0("p = ", formatC(p_value, format = "f", digits = 3), star))
  } else {
    rounded <- round(p_value, 2)
    if (rounded == 0) {
      return(paste0("p = 0.001", star))
    } else {
      return(paste0("p = ", rounded, star))
    }
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
  kw <- stats::kruskal.test(Percentage_SOX10neg_gain ~ Classification, data = df_region)
  
  dunn <- dunn.test::dunn.test(
    x = df_region$Percentage_SOX10neg_gain,
    g = df_region$Classification,
    method = "bonferroni",
    kw = FALSE,
    list = TRUE
  )
  
  dunn_tbl <- tibble::tibble(
    Region = as.character(unique(df_region$Region)),
    test = "Dunn (Bonferroni)",
    comparison = dunn$comparisons,
    z = as.numeric(dunn$Z),
    p_raw = as.numeric(dunn$P),
    p_adj = as.numeric(dunn$P.adjusted)
  )
  
  kw_tbl <- tibble::tibble(
    Region = as.character(unique(df_region$Region)),
    test = "Kruskal-Wallis",
    comparison = "Classification (overall)",
    statistic = as.numeric(kw$statistic),
    df = as.numeric(kw$parameter),
    p_raw = as.numeric(kw$p.value),
    p_adj = NA_real_
  )
  
  list(kw = kw_tbl, dunn = dunn_tbl)
}

parse_dunn_pairs <- function(dunn_tbl) {
  dunn_tbl %>%
    dplyr::mutate(
      group1 = sub(" - .*", "", comparison),
      group2 = sub(".* - ", "", comparison)
    )
}

plot_region_across_groups <- function(df_region, brackets, region_label) {
  y_max <- GLOBAL_Y_MAX
  
  ggplot2::ggplot(
    df_region,
    ggplot2::aes(x = Classification, y = Percentage_SOX10neg_gain, fill = Classification)
  ) +
    ggplot2::annotate(
      "rect",
      xmin = 0.5, xmax = length(x_order) + 0.5,
      ymin = y_max + 0.5, ymax = y_max + 4,
      fill = "white", color = "black", linewidth = 1
    ) +
    ggplot2::annotate(
      "text",
      x = (length(x_order) + 1) / 2,
      y = y_max + 2.25,
      label = region_label,
      size = 8, fontface = "bold"
    ) +
    ggplot2::geom_boxplot(
      color = "black", outlier.shape = NA,
      alpha = 0.3, linewidth = 0.85
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(color = Classification),
      width = 0.1, size = 3, alpha = 0.8
    ) +
    ggpubr::geom_bracket(
      data = brackets,
      ggplot2::aes(xmin = xmin, xmax = xmax, y.position = y.position, label = label),
      inherit.aes = FALSE,
      label.size = 5,
      size = 0.7,
      tip.length = 0.01
    ) +
    ggplot2::labs(x = "", y = expression("% SOX10-" ~ italic(SNCA) ~ "Gains")) +
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
      limits = c(-1, y_max + 4.5),
      expand = ggplot2::expansion(mult = c(0.05, 0))
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18),
      axis.text.x = ggplot2::element_text(size = 16),
      axis.text.y = ggplot2::element_text(size = 16),
      legend.position = "none",
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      plot.background = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      panel.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank(),
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
  
  y_positions <- c(32, 30, 28)
  
  desired <- tibble::tibble(
    group1 = c("OPCA", "OPCA", "SND"),
    group2 = c("SND", "Control", "Control")
  )
  
  dunn_pairs2 <- desired %>%
    dplyr::left_join(
      dunn_pairs %>% dplyr::select(group1, group2, p_adj),
      by = c("group1", "group2")
    ) %>%
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
  
  brackets <- make_brackets_from_pairs(dunn_pairs2, x_order, y_positions)
  
  p <- plot_region_across_groups(df_reg, brackets, reg)
  plots_region[[reg]] <- p
  
  ggplot2::ggsave(
    filename = file.path(plot_dir, paste0("SOX10neg_Gains_", reg, "_AcrossGroups.png")),
    plot = p, width = 7, height = 7, dpi = 600
  )
}

cerebellum_plot <- plots_region[["Cerebellum"]]
putamen_plot <- plots_region[["Putamen"]]
sn_plot <- plots_region[["SN"]]

putamen_plot_clean <- putamen_plot +
  ggplot2::theme(
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  )

sn_plot_clean <- sn_plot +
  ggplot2::theme(
    axis.title.y = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank()
  )

combined_plot <- cerebellum_plot +
  putamen_plot_clean +
  sn_plot_clean +
  patchwork::plot_layout(ncol = 3)

ggplot2::ggsave(
  filename = file.path(plot_dir, "SOX10neg_Gains_AcrossGroups_Combined.png"),
  plot = combined_plot,
  width = 14,
  height = 7,
  dpi = 600
)

combined_plot

kw_region_tbl <- dplyr::bind_rows(all_kw_region)
dunn_region_tbl <- dplyr::bind_rows(all_dunn_region)

readr::write_tsv(
  kw_region_tbl,
  file.path(stats_dir, "SOX10neg_Gains_AcrossGroupsWithinRegion_Kruskal.tsv")
)

readr::write_tsv(
  dunn_region_tbl,
  file.path(stats_dir, "SOX10neg_Gains_AcrossGroupsWithinRegion_Dunn.tsv")
)

within_group_across_regions_opca <- function(df) {
  df_opca <- df %>% dplyr::filter(Classification == "OPCA")
  if (dplyr::n_distinct(df_opca$Region) < 2) return(list(anova = tibble::tibble(), tukey = tibble::tibble()))
  
  fit <- stats::aov(Percentage_SOX10neg_gain ~ Region, data = df_opca)
  anova_sum <- summary(fit)[[1]]
  
  anova_tbl <- tibble::tibble(
    Classification = "OPCA",
    test = "One-way ANOVA",
    term = rownames(anova_sum),
    sum_sq = anova_sum$`Sum Sq`,
    mean_sq = anova_sum$`Mean Sq`,
    f_value = anova_sum$`F value`,
    p_value = anova_sum$`Pr(>F)`
  )
  
  tuk <- stats::TukeyHSD(fit, "Region")
  tuk_df <- as.data.frame(tuk$Region) %>%
    tibble::rownames_to_column("comparison") %>%
    tibble::as_tibble() %>%
    dplyr::rename(
      diff = diff,
      lwr = lwr,
      upr = upr,
      p_adj = `p adj`
    ) %>%
    dplyr::mutate(
      Classification = "OPCA",
      test = "TukeyHSD"
    )
  
  list(anova = anova_tbl, tukey = tuk_df)
}

within_group_across_regions_kw_dunn <- function(df, cls) {
  df_cls <- df %>% dplyr::filter(Classification == cls)
  if (dplyr::n_distinct(df_cls$Region) < 2) return(list(kw = tibble::tibble(), dunn = tibble::tibble()))
  
  kw <- stats::kruskal.test(Percentage_SOX10neg_gain ~ Region, data = df_cls)
  
  dunn <- dunn.test::dunn.test(
    x = df_cls$Percentage_SOX10neg_gain,
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

opca_res <- within_group_across_regions_opca(cnv_clean)
snd_res <- within_group_across_regions_kw_dunn(cnv_clean, "SND")
ctl_res <- within_group_across_regions_kw_dunn(cnv_clean, "Control")

readr::write_tsv(
  opca_res$anova,
  file.path(stats_dir, "SOX10neg_Gains_OPCA_Regions_ANOVA.tsv")
)
readr::write_tsv(
  opca_res$tukey,
  file.path(stats_dir, "SOX10neg_Gains_OPCA_Regions_Tukey.tsv")
)

readr::write_tsv(
  snd_res$kw,
  file.path(stats_dir, "SOX10neg_Gains_SND_Regions_Kruskal.tsv")
)
readr::write_tsv(
  snd_res$dunn,
  file.path(stats_dir, "SOX10neg_Gains_SND_Regions_Dunn.tsv")
)

readr::write_tsv(
  ctl_res$kw,
  file.path(stats_dir, "SOX10neg_Gains_Control_Regions_Kruskal.tsv")
)
readr::write_tsv(
  ctl_res$dunn,
  file.path(stats_dir, "SOX10neg_Gains_Control_Regions_Dunn.tsv")
)

summary_tbl <- cnv_clean %>%
  dplyr::group_by(Classification, Region) %>%
  dplyr::summarise(
    n_obs = dplyr::n(),
    n_donors = dplyr::n_distinct(Donor_ID),
    n_donor_region = dplyr::n_distinct(DonorRegion),
    mean = mean(Percentage_SOX10neg_gain, na.rm = TRUE),
    sd = stats::sd(Percentage_SOX10neg_gain, na.rm = TRUE),
    median = stats::median(Percentage_SOX10neg_gain, na.rm = TRUE),
    iqr_low = stats::quantile(Percentage_SOX10neg_gain, 0.25, na.rm = TRUE),
    iqr_high = stats::quantile(Percentage_SOX10neg_gain, 0.75, na.rm = TRUE),
    min = min(Percentage_SOX10neg_gain, na.rm = TRUE),
    max = max(Percentage_SOX10neg_gain, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_tsv(
  summary_tbl,
  file.path(stats_dir, "SOX10neg_Gains_ByClassification_ByRegion_Summary.tsv")
)
