source(here("R", "00_setup", "01_setup_environment.R"))

ggplot2::theme_set(ggplot2::theme_minimal(base_family = "Helvetica"))

add_consistent_pval <- function(plot_obj, data, y_var, p_val_text, Group_x = c(1, 2)) {
  max_y <- max(data[[y_var]], na.rm = TRUE)
  min_y <- min(data[[y_var]], na.rm = TRUE)
  y_range <- max_y - min_y
  
  if (is.na(y_range) || y_range == 0) y_range <- max_y
  
  y_top <- max_y + y_range * 0.12
  bracket_y <- max_y + y_range * 0.08
  tick_length <- y_range * 0.025
  
  plot_obj +
    ggplot2::expand_limits(y = y_top + tick_length) +
    ggplot2::annotate("text", x = mean(Group_x), y = y_top, label = p_val_text, size = 5) +
    ggplot2::annotate(
      "segment",
      x = Group_x[1], xend = Group_x[2],
      y = bracket_y, yend = bracket_y, linewidth = 0.7
    ) +
    ggplot2::annotate(
      "segment",
      x = Group_x[1], xend = Group_x[1],
      y = bracket_y, yend = bracket_y - tick_length, linewidth = 0.7
    ) +
    ggplot2::annotate(
      "segment",
      x = Group_x[2], xend = Group_x[2],
      y = bracket_y, yend = bracket_y - tick_length, linewidth = 0.7
    )
}

format_p_for_plot <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.0001) return("p < 0.0001*")
  if (p < 0.05) return(paste0("p = ", format(round(p, 3), nsmall = 3), "*"))
  paste0("p = ", format(round(p, 3), nsmall = 3))
}

colors_box <- c(
  "MSA" = scales::alpha("darkseagreen4", 0.5),
  "Control" = scales::alpha("darkmagenta", 0.5)
)

colors_points <- c(
  "MSA" = "darkseagreen4",
  "Control" = "darkmagenta"
)

out_group_dir <- file.path(plots_gh2ax_dir, "Group_Comparisons")
out_paired_dir <- file.path(plots_gh2ax_dir, "Paired_Comparisons")
dir.create(out_group_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_paired_dir, recursive = TRUE, showWarnings = FALSE)

get_summary <- function(data, value_col, Group_col, comparison) {
  data %>%
    dplyr::group_by(.data[[Group_col]]) %>%
    dplyr::summarise(
      n = sum(!is.na(.data[[value_col]])),
      mean = mean(.data[[value_col]], na.rm = TRUE),
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

make_group_plot <- function(df, y_var, title, y_lab) {
  df <- df %>%
    dplyr::filter(group_msa %in% c("MSA", "Control")) %>%
    tidyr::drop_na(.data[[y_var]]) %>%
    dplyr::mutate(group_msa = factor(group_msa, levels = c("Control", "MSA")))
  
  test <- stats::wilcox.test(stats::as.formula(paste0(y_var, " ~ group_msa")), data = df, exact = TRUE)
  p_txt <- format_p_for_plot(test$p.value)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(group_msa, .data[[y_var]], fill = group_msa)) +
    ggplot2::geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.45, color = "black", linewidth = 0.7) +
    ggplot2::geom_jitter(ggplot2::aes(color = group_msa), width = 0.1, size = 3) +
    ggplot2::scale_fill_manual(values = colors_box) +
    ggplot2::scale_color_manual(values = colors_points) +
    ggplot2::labs(title = title, x = NULL, y = y_lab) +
    ggplot2::theme(
      legend.position = "none",
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12)
    )
  
  p <- add_consistent_pval(p, df, y_var, p_txt)
  list(plot = p, p_value = test$p.value)
}

make_paired_asyn_plot <- function(df) {
  df <- df %>%
    dplyr::filter(group_msa == "MSA") %>%
    tidyr::drop_na(Percentage_aSyn_pos_gh2ax_pos, Percentage_aSyn_neg_gh2ax_pos)
  
  test <- stats::wilcox.test(
    df$Percentage_aSyn_pos_gh2ax_pos,
    df$Percentage_aSyn_neg_gh2ax_pos,
    paired = TRUE,
    exact = TRUE
  )
  p_txt <- format_p_for_plot(test$p.value)
  
  long <- df %>%
    dplyr::select(
      Donor,
      Region,
      Percentage_aSyn_pos_gh2ax_pos,
      Percentage_aSyn_neg_gh2ax_pos
    ) %>%
    tidyr::pivot_longer(
      cols = c(Percentage_aSyn_neg_gh2ax_pos, Percentage_aSyn_pos_gh2ax_pos),
      names_to = "aSyn_Status",
      values_to = "percentage_gh2ax_foci"
    ) %>%
    dplyr::mutate(
      aSyn_Status = factor(
        aSyn_Status,
        levels = c("Percentage_aSyn_neg_gh2ax_pos", "Percentage_aSyn_pos_gh2ax_pos"),
        labels = c("aSyn-", "aSyn+")
      )
    )
  
  p <- ggplot2::ggplot(long, ggplot2::aes(aSyn_Status, percentage_gh2ax_foci, fill = aSyn_Status)) +
    ggplot2::geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.45, color = "black", linewidth = 0.7) +
    ggplot2::geom_line(
      ggplot2::aes(group = interaction(Donor, Region)),
      color = "gray40",
      linewidth = 0.6
    ) +
    ggplot2::geom_jitter(width = 0.06, size = 3) +
    ggplot2::scale_fill_manual(values = c(
      "aSyn-" = scales::alpha("darkmagenta", 0.5),
      "aSyn+" = scales::alpha("darkseagreen4", 0.5)
    )) +
    ggplot2::labs(title = "MSA: aSyn+ vs aSyn− (paired)", x = NULL, y = "% cells with ≥1 γH2AX foci") +
    ggplot2::theme(
      legend.position = "none",
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12)
    )
  
  p <- add_consistent_pval(p, long, "percentage_gh2ax_foci", p_txt)
  list(plot = p, p_value = test$p.value)
}

make_most_least_plot <- function(df) {
  df <- df %>%
    dplyr::filter(Group %in% c("OPCA", "SND"), Region %in% c("Cerebellum", "Putamen")) %>%
    dplyr::mutate(
      Affectedness = dplyr::case_when(
        Group == "OPCA" & Region == "Cerebellum" ~ "Most affected",
        Group == "SND" & Region == "Putamen" ~ "Most affected",
        Group == "OPCA" & Region == "Putamen" ~ "Least affected",
        Group == "SND" & Region == "Cerebellum" ~ "Least affected",
        TRUE ~ NA_character_
      ),
      Affectedness = factor(Affectedness, levels = c("Least affected", "Most affected"))
    ) %>%
    tidyr::drop_na(Affectedness, Percentage_all_gh2ax_foci)
  
  test <- stats::wilcox.test(Percentage_all_gh2ax_foci ~ Affectedness, data = df, exact = TRUE)
  p_txt <- format_p_for_plot(test$p.value)
  
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(Affectedness, Percentage_all_gh2ax_foci, fill = Affectedness)
  ) +
    ggplot2::geom_boxplot(
      alpha = 0.3, outlier.shape = NA, width = 0.55,
      color = "black", linewidth = 0.7
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(color = Affectedness),
      width = 0.1,
      size = 3
    ) +
    ggplot2::scale_fill_manual(values = c(
      "Least affected" = scales::alpha("darkmagenta", 0.5),
      "Most affected"  = scales::alpha("darkseagreen4", 0.5)
    )) +
    ggplot2::scale_color_manual(values = c(
      "Least affected" = "darkmagenta",
      "Most affected"  = "darkseagreen4"
    )) +
    ggplot2::labs(
      title = "MSA Regions: Most vs Least affected",
      x = NULL,
      y = "% cells with ≥1 γH2AX foci"
    ) +
    ggplot2::theme(
      legend.position = "none",
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_text(size = 14),
      axis.text = ggplot2::element_text(size = 12)
    )
  
  
  p <- add_consistent_pval(p, df, "Percentage_all_gh2ax_foci", p_txt)
  list(plot = p, p_value = test$p.value)
}

res_all <- make_group_plot(
  gh2ax_data,
  "Percentage_all_gh2ax_foci",
  "All cells: MSA vs Control",
  "% cells with ≥1 γH2AX foci"
)

res_sox10pos <- make_group_plot(
  gh2ax_data,
  "Percentage_SOX10_pos_gh2ax_pos",
  "SOX10+ cells: MSA vs Control",
  "% SOX10+ cells with ≥1 γH2AX foci"
)

res_sox10neg <- make_group_plot(
  gh2ax_data,
  "Percentage_SOX10_neg_gh2ax_pos",
  "SOX10− cells: MSA vs Control",
  "% SOX10− cells with ≥1 γH2AX foci"
)

res_asyn <- make_paired_asyn_plot(gh2ax_data)

res_aff <- make_most_least_plot(gh2ax_data)

ggplot2::ggsave(
  filename = file.path(out_group_dir, "gH2AX_All_cells_MSA_vs_Control.png"),
  plot = res_all$plot,
  width = 5,
  height = 5,
  dpi = 600
)

ggplot2::ggsave(
  filename = file.path(out_group_dir, "gH2AX_SOX10pos_cells_MSA_vs_Control.png"),
  plot = res_sox10pos$plot,
  width = 5,
  height = 5,
  dpi = 600
)

ggplot2::ggsave(
  filename = file.path(out_group_dir, "gH2AX_SOX10neg_cells_MSA_vs_Control.png"),
  plot = res_sox10neg$plot,
  width = 5,
  height = 5,
  dpi = 600
)

ggplot2::ggsave(
  filename = file.path(out_paired_dir, "gH2AX_MSA_aSynpos_vs_aSynneg_paired.png"),
  plot = res_asyn$plot,
  width = 5,
  height = 5,
  dpi = 600
)

ggplot2::ggsave(
  filename = file.path(out_group_dir, "gH2AX_MSA_most_vs_least_affected_regions.png"),
  plot = res_aff$plot,
  width = 5,
  height = 5,
  dpi = 600
)

summary_stats <- dplyr::bind_rows(
  get_summary(
    gh2ax_data %>% dplyr::filter(group_msa %in% c("MSA", "Control")) %>% tidyr::drop_na(Percentage_all_gh2ax_foci),
    "Percentage_all_gh2ax_foci",
    "group_msa",
    "All cells: MSA vs Control"
  ),
  get_summary(
    gh2ax_data %>% dplyr::filter(group_msa %in% c("MSA", "Control")) %>% tidyr::drop_na(Percentage_SOX10_pos_gh2ax_pos),
    "Percentage_SOX10_pos_gh2ax_pos",
    "group_msa",
    "SOX10+ cells: MSA vs Control"
  ),
  get_summary(
    gh2ax_data %>% dplyr::filter(group_msa %in% c("MSA", "Control")) %>% tidyr::drop_na(Percentage_SOX10_neg_gh2ax_pos),
    "Percentage_SOX10_neg_gh2ax_pos",
    "group_msa",
    "SOX10− cells: MSA vs Control"
  ),
  get_summary(
    {
      tmp <- gh2ax_data %>% dplyr::filter(group_msa == "MSA") %>% tidyr::drop_na(Percentage_aSyn_pos_gh2ax_pos, Percentage_aSyn_neg_gh2ax_pos)
      tmp_long <- tmp %>%
        dplyr::select(Donor, Region, Percentage_aSyn_pos_gh2ax_pos, Percentage_aSyn_neg_gh2ax_pos) %>%
        tidyr::pivot_longer(
          cols = c(Percentage_aSyn_neg_gh2ax_pos, Percentage_aSyn_pos_gh2ax_pos),
          names_to = "aSyn_Status",
          values_to = "percentage_gh2ax_foci"
        ) %>%
        dplyr::mutate(
          aSyn_Status = factor(
            aSyn_Status,
            levels = c("Percentage_aSyn_neg_gh2ax_pos", "Percentage_aSyn_pos_gh2ax_pos"),
            labels = c("aSyn-", "aSyn+")
          )
        )
      tmp_long
    },
    "percentage_gh2ax_foci",
    "aSyn_Status",
    "MSA: aSyn+ vs aSyn− (paired)"
  ),
  get_summary(
    gh2ax_data %>%
      dplyr::filter(Group %in% c("OPCA", "SND"), Region %in% c("Cerebellum", "Putamen")) %>%
      dplyr::mutate(
        Affectedness = dplyr::case_when(
          Group == "OPCA" & Region == "Cerebellum" ~ "Most affected",
          Group == "SND" & Region == "Putamen" ~ "Most affected",
          Group == "OPCA" & Region == "Putamen" ~ "Least affected",
          Group == "SND" & Region == "Cerebellum" ~ "Least affected",
          TRUE ~ NA_character_
        )
      ) %>%
      tidyr::drop_na(Affectedness, Percentage_all_gh2ax_foci),
    "Percentage_all_gh2ax_foci",
    "Affectedness",
    "MSA: Most vs Least affected"
  )
) %>%
  dplyr::left_join(
    tibble::tibble(
      comparison = c(
        "All cells: MSA vs Control",
        "SOX10+ cells: MSA vs Control",
        "SOX10− cells: MSA vs Control",
        "MSA: aSyn+ vs aSyn− (paired)",
        "MSA: Most vs Least affected"
      ),
      test = "Wilcoxon",
      paired = c(FALSE, FALSE, FALSE, TRUE, FALSE),
      p_value = c(
        res_all$p_value,
        res_sox10pos$p_value,
        res_sox10neg$p_value,
        res_asyn$p_value,
        res_aff$p_value
      )
    ),
    by = "comparison"
  ) %>%
  dplyr::mutate(
    dplyr::across(c(mean, median, q1, q3, iqr), ~ round(., 2)),
    p_value = signif(p_value, 3)
  ) %>%
  dplyr::mutate(
    dplyr::across(dplyr::where(is.character), ~ gsub("\u2212", "-", .))
  )

readr::write_tsv(
  summary_stats,
  file.path(
    stats_gh2ax_dir,
    "gH2AX_Wilcoxon_SummaryStats_Group_and_Paired_Comparisons.tsv"
  ),
  na = ""
)
