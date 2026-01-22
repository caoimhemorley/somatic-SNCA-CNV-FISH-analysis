source(here("R", "00_setup", "01_setup_environment.R"))

# 1) Summarize per Donor_ID (per case) within Classification x Region
per_case_cell_count_summary <- cnv_data %>%
  group_by(Classification, Region, Donor_ID) %>%
  summarise(
    Total_cells = if (all(is.na(Total_cells))) NA_real_
    else sum(Total_cells, na.rm = TRUE),
    
    SOX10pos = if (all(is.na(SOX10pos))) NA_real_
    else sum(SOX10pos, na.rm = TRUE),
    
    SOX10neg = if (all(is.na(SOX10neg))) NA_real_
    else sum(SOX10neg, na.rm = TRUE),
    
    .groups = "drop"
  )


readr::write_tsv(
  per_case_cell_count_summary,
  file.path(stats_cellcounts_dir, "FISH_cell_counts_PerCase_SummedCounts_by_Classification_Region_Donor.tsv"),
  na = ""
)

# 2) Summarize by Classification x Region (using per-case totals)
summary_table_cell_count <- per_case_cell_count_summary %>%
  dplyr::group_by(Classification, Region) %>%
  dplyr::summarise(
    Total_cells_sum = sum(Total_cells, na.rm = TRUE),
    SOX10pos_total = sum(SOX10pos, na.rm = TRUE),
    SOX10neg_total = sum(SOX10neg, na.rm = TRUE),
    
    Total_cells_mean = mean(Total_cells, na.rm = TRUE),
    Total_cells_sd = stats::sd(Total_cells, na.rm = TRUE),
    Total_cells_median = stats::median(Total_cells, na.rm = TRUE),
    Total_cells_IQR_low = stats::quantile(Total_cells, 0.25, na.rm = TRUE),
    Total_cells_IQR_high = stats::quantile(Total_cells, 0.75, na.rm = TRUE),
    
    SOX10pos_mean = mean(SOX10pos, na.rm = TRUE),
    SOX10pos_sd = stats::sd(SOX10pos, na.rm = TRUE),
    SOX10pos_median = stats::median(SOX10pos, na.rm = TRUE),
    SOX10pos_IQR_low = stats::quantile(SOX10pos, 0.25, na.rm = TRUE),
    SOX10pos_IQR_high = stats::quantile(SOX10pos, 0.75, na.rm = TRUE),
    
    SOX10neg_mean = mean(SOX10neg, na.rm = TRUE),
    SOX10neg_sd = stats::sd(SOX10neg, na.rm = TRUE),
    SOX10neg_median = stats::median(SOX10neg, na.rm = TRUE),
    SOX10neg_IQR_low = stats::quantile(SOX10neg, 0.25, na.rm = TRUE),
    SOX10neg_IQR_high = stats::quantile(SOX10neg, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_tsv(
  summary_table_cell_count,
  file.path(stats_cellcounts_dir, "FISH_cell_counts_SummaryStats_PerCaseTotals_by_Classification_Region.tsv"),
  na = ""
)

# 3) Kruskal-Wallis tests (per-case totals; across Classification x Region interaction)
kruskal_total <- stats::kruskal.test(
  Total_cells ~ interaction(Classification, Region),
  data = per_case_cell_count_summary
)

kruskal_sox10pos <- stats::kruskal.test(
  SOX10pos ~ interaction(Classification, Region),
  data = per_case_cell_count_summary
)

kruskal_sox10neg <- stats::kruskal.test(
  SOX10neg ~ interaction(Classification, Region),
  data = per_case_cell_count_summary
)

kruskal_results <- tibble::tibble(
  Variable = c("Total_cells", "SOX10pos", "SOX10neg"),
  Test = "Kruskal-Wallis",
  Grouping = "interaction(Classification, Region)",
  Statistic = c(
    as.numeric(kruskal_total$statistic),
    as.numeric(kruskal_sox10pos$statistic),
    as.numeric(kruskal_sox10neg$statistic)
  ),
  df = c(
    as.numeric(kruskal_total$parameter),
    as.numeric(kruskal_sox10pos$parameter),
    as.numeric(kruskal_sox10neg$parameter)
  ),
  p_value = c(
    as.numeric(kruskal_total$p.value),
    as.numeric(kruskal_sox10pos$p.value),
    as.numeric(kruskal_sox10neg$p.value)
  )
)

readr::write_tsv(
  kruskal_results,
  file.path(stats_cellcounts_dir, "FISH_cell_counts_Kruskal_by_Classification_RegionInteraction.tsv"),
  na = ""
)

# 4) Prepare data for plotting (per-case totals, long format)
cellcount_plot_data <- per_case_cell_count_summary %>%
  tidyr::pivot_longer(
    cols = c(Total_cells, SOX10pos),
    names_to = "Variable",
    values_to = "Count"
  ) %>%
  dplyr::mutate(
    Variable = factor(Variable, levels = c("Total_cells", "SOX10pos"))
  )

# 5) Plot
cellcounts_total_sox10_boxplot <- ggplot2::ggplot(
  cellcount_plot_data,
  ggplot2::aes(x = Region, y = Count, fill = Variable)
) +
  ggplot2::geom_boxplot(
    position = ggplot2::position_dodge(0.8),
    width = 0.5,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  ggplot2::geom_jitter(
    ggplot2::aes(color = Variable),
    position = ggplot2::position_jitterdodge(
      jitter.width = 0.1,
      jitter.height = 0,
      dodge.width = 0.8
    ),
    size = 1,
    alpha = 0.7
  ) +
  ggplot2::scale_fill_manual(
    values = c("Total_cells" = "grey", "SOX10pos" = "pink2"),
    labels = c("Total_cells" = "Total Cells", "SOX10pos" = "SOX10+ Cells")
  ) +
  ggplot2::scale_color_manual(
    values = c("Total_cells" = "black", "SOX10pos" = "deeppink3"),
    labels = c("Total_cells" = "Total Cells", "SOX10pos" = "SOX10+ Cells")
  ) +
  ggplot2::labs(
    x = "Region",
    y = "Cell Count",
    title = "Total and SOX10+ Cell Counts by Region and Classification"
  ) +
  ggplot2::facet_wrap(~ Classification, scales = "free_x") +
  ggplot2::coord_cartesian(ylim = c(0, 200)) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 12),
    axis.title = ggplot2::element_text(size = 12),
    plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
    strip.text = ggplot2::element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.text = ggplot2::element_text(size = 12),
    legend.title = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
  )

# 6) Save plot
ggplot2::ggsave(
  filename = file.path(
    plots_cellcounts_dir,
    "FISH_cell_counts_Total_and_SOX10pos_PerCaseTotals_by_Region_FacetedBy_Classification.png"
  ),
  plot = cellcounts_total_sox10_boxplot,
  width = 9.7,
  height = 6,
  dpi = 1200
)
