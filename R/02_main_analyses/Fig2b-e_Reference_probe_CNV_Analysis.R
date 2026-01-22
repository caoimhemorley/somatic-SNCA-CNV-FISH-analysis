source(here("R", "00_setup", "01_setup_environment.R"))

gain_vars <- c(
  "Percentage_SOX10pos_Ref_gain", "Percentage_SOX10pos_SNCA_gain",
  "Percentage_SOX10neg_Ref_gain", "Percentage_SOX10neg_SNCA_gain"
)

loss_vars <- c(
  "Percentage_SOX10pos_Ref_loss", "Percentage_SOX10pos_SNCA_loss",
  "Percentage_SOX10neg_Ref_loss", "Percentage_SOX10neg_SNCA_loss"
)

reshape_plot_data <- function(df, vars, event_type) {
  df %>%
    dplyr::mutate(
      Group_recode = dplyr::if_else(Group %in% c("SND", "OPCA"), "MSA", Group),
      Group_recode = factor(Group_recode, levels = c("Control", "MSA"))
    ) %>%
    dplyr::select(Donor_ID, Region, Group_recode, dplyr::all_of(vars)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(vars),
      names_to = "Variable",
      values_to = "Percentage"
    ) %>%
    dplyr::mutate(
      Event = event_type,
      CellType = dplyr::case_when(
        grepl("SOX10pos", Variable) ~ "SOX10+ cells",
        grepl("SOX10neg", Variable) ~ "SOX10- cells",
        TRUE ~ NA_character_
      ),
      CellType = factor(CellType, levels = c("SOX10+ cells", "SOX10- cells")),
      Probe = dplyr::case_when(
        grepl("Ref", Variable) ~ "Reference",
        grepl("SNCA", Variable) ~ "SNCA",
        TRUE ~ NA_character_
      ),
      Probe = factor(Probe, levels = c("Reference", "SNCA"))
    ) %>%
    dplyr::filter(!is.na(CellType), !is.na(Probe))
}

plot_gain_data <- reshape_plot_data(ref_cnv_data, gain_vars, "Gain")
plot_loss_data <- reshape_plot_data(ref_cnv_data, loss_vars, "Loss")

compute_between_pvals <- function(df) {
  df %>%
    dplyr::group_by(CellType, Probe) %>%
    dplyr::summarise(
      p_raw = suppressWarnings(
        stats::wilcox.test(Percentage ~ Group_recode, data = dplyr::pick(dplyr::everything()), exact = TRUE)$p.value
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      p_adj = pmin(p_raw * 4, 1),
      p_label = dplyr::case_when(
        is.na(p_adj) ~ NA_character_,
        p_adj < 0.0001 ~ "p < 0.0001*",
        p_adj < 0.001  ~ "p < 0.001*",
        p_adj < 0.01   ~ "p < 0.01*",
        p_adj < 0.05   ~ paste0("p = ", signif(p_adj, 2), "*"),
        TRUE           ~ paste0("p = ", signif(p_adj, 2))
      )
    )
}

compute_within_pvals <- function(df) {
  df %>%
    dplyr::filter(Probe %in% c("Reference", "SNCA")) %>%
    dplyr::group_by(CellType, Group_recode, Donor_ID, Region) %>%
    dplyr::summarise(
      ref_val = Percentage[Probe == "Reference"][1],
      snca_val = Percentage[Probe == "SNCA"][1],
      .groups = "drop"
    ) %>%
    dplyr::group_by(CellType, Group_recode) %>%
    dplyr::summarise(
      p_raw = suppressWarnings(
        stats::wilcox.test(ref_val, snca_val, paired = TRUE, exact = TRUE)$p.value
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      p_adj = pmin(p_raw * 4, 1),
      Probe = "SNCA vs REF",
      p_label = dplyr::case_when(
        is.na(p_adj) ~ NA_character_,
        p_adj < 0.0001 ~ "p < 0.0001*",
        p_adj < 0.001  ~ "p < 0.001*",
        p_adj < 0.01   ~ "p < 0.01*",
        p_adj < 0.05   ~ paste0("p = ", signif(p_adj, 2), "*"),
        TRUE           ~ paste0("p = ", signif(p_adj, 2))
      )
    )
}

gain_between <- compute_between_pvals(plot_gain_data)
loss_between <- compute_between_pvals(plot_loss_data)
gain_within <- compute_within_pvals(plot_gain_data)
loss_within <- compute_within_pvals(plot_loss_data)

make_brackets <- function(df, between, within, dodge_width = 0.75, global_ymax = NULL) {
  df <- df %>%
    dplyr::mutate(
      Probe = factor(Probe, levels = c("Reference", "SNCA")),
      Group_recode = factor(Group_recode, levels = c("Control", "MSA"))
    )
  
  x_map <- setNames(seq_along(levels(df$Probe)), levels(df$Probe))
  n_groups <- length(levels(df$Group_recode))
  
  ymax_probe <- df %>%
    dplyr::group_by(CellType, Probe) %>%
    dplyr::summarise(ymax = max(Percentage, na.rm = TRUE), .groups = "drop")
  if (!is.null(global_ymax)) ymax_probe$ymax <- global_ymax
  
  between_brackets <- between %>%
    dplyr::left_join(ymax_probe, by = c("CellType", "Probe")) %>%
    dplyr::mutate(
      xmin = x_map[as.character(Probe)] - 0.2,
      xmax = x_map[as.character(Probe)] + 0.2,
      y.position = ymax * 1.05,
      group_type = paste(CellType, Probe, "between")
    )
  
  ymax_cell <- df %>%
    dplyr::group_by(CellType) %>%
    dplyr::summarise(ymax = max(Percentage, na.rm = TRUE), .groups = "drop")
  if (!is.null(global_ymax)) ymax_cell$ymax <- global_ymax
  
  dodge_offsets <- seq(
    -dodge_width / 2 + dodge_width / (2 * n_groups),
    dodge_width / 2 - dodge_width / (2 * n_groups),
    length.out = n_groups
  )
  names(dodge_offsets) <- levels(df$Group_recode)
  
  within_brackets <- within %>%
    dplyr::left_join(ymax_cell, by = "CellType") %>%
    dplyr::mutate(
      xmin = x_map["Reference"] + dodge_offsets[Group_recode],
      xmax = x_map["SNCA"] + dodge_offsets[Group_recode],
      y.position = ymax * dplyr::if_else(Group_recode == "Control", 1.15, 1.25),
      group_type = paste(CellType, Group_recode, "within")
    )
  
  dplyr::bind_rows(between_brackets, within_brackets)
}

gain_global_ymax <- max(plot_gain_data$Percentage, na.rm = TRUE) * 1.05
loss_global_ymax <- max(plot_loss_data$Percentage, na.rm = TRUE) * 1.05

gain_brackets <- make_brackets(plot_gain_data, gain_between, gain_within, global_ymax = gain_global_ymax)
loss_brackets <- make_brackets(plot_loss_data, loss_between, loss_within, global_ymax = loss_global_ymax)

plot_event <- function(df, brackets, event_label) {
  brackets <- brackets %>%
    dplyr::mutate(
      p_label = dplyr::if_else(!is.na(p_adj) & p_adj < 0.05, paste0(p_label, " "), p_label)
    )
  
  ggplot2::ggplot(df, ggplot2::aes(x = Probe, y = Percentage, fill = Group_recode)) +
    ggplot2::geom_boxplot(
      alpha = 0.3, width = 0.5, color = "black", linewidth = 0.7,
      outlier.shape = NA, position = ggplot2::position_dodge(width = 0.75)
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(color = Group_recode),
      position = ggplot2::position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
      size = 2
    ) +
    ggplot2::facet_wrap(~ CellType, nrow = 1) +
    ggpubr::geom_bracket(
      data = brackets,
      ggplot2::aes(
        xmin = xmin, xmax = xmax,
        y.position = y.position,
        label = p_label,
        group = group_type
      ),
      inherit.aes = FALSE,
      size = 0.7,
      tip.length = 0.02,
      label.size = 4.5
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
    ggplot2::scale_fill_manual(name = "Group", values = c("Control" = "deepskyblue4", "MSA" = "slateblue4")) +
    ggplot2::scale_color_manual(name = "Group", values = c("Control" = "deepskyblue4", "MSA" = "slateblue4")) +
    ggplot2::labs(y = paste0("% of cells with ", tolower(event_label)), x = "Probe") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 18, face = "bold", colour = "black"),
      
      axis.title.y = ggplot2::element_text(size = 18, face = "bold", colour = "black"),
      axis.title.x = ggplot2::element_text(size = 18, face = "bold", colour = "black"),
      
      axis.text.x = ggplot2::element_text(size = 16, face = "bold", colour = "black"),
      axis.text.y = ggplot2::element_text(size = 14, colour = "black"),
      
      legend.position = c(0.5, 0.98),
      legend.justification = c("center", "top"),
      legend.direction = "horizontal",
      legend.background = ggplot2::element_rect(
        fill = scales::alpha("white", 0.75), color = NA
      ),
      legend.key.height = grid::unit(0.6, "lines"),
      legend.key.width  = grid::unit(1.2, "lines"),
      legend.text = ggplot2::element_text(size = 14, face = "bold", colour = "black"),
      legend.title = ggplot2::element_text(size = 14, face = "bold", colour = "black"),
      
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1)
    )
}


plot_gain_sox10pos <- plot_event(
  dplyr::filter(plot_gain_data, CellType == "SOX10+ cells"),
  dplyr::filter(gain_brackets, CellType == "SOX10+ cells"),
  "Gains"
)

plot_gain_sox10neg <- plot_event(
  dplyr::filter(plot_gain_data, CellType == "SOX10- cells"),
  dplyr::filter(gain_brackets, CellType == "SOX10- cells"),
  "Gains"
)

plot_loss_sox10pos <- plot_event(
  dplyr::filter(plot_loss_data, CellType == "SOX10+ cells"),
  dplyr::filter(loss_brackets, CellType == "SOX10+ cells"),
  "Losses"
)

plot_loss_sox10neg <- plot_event(
  dplyr::filter(plot_loss_data, CellType == "SOX10- cells"),
  dplyr::filter(loss_brackets, CellType == "SOX10- cells"),
  "Losses"
)

plot_list <- list(
  RefProbe_Gain_SOX10pos = plot_gain_sox10pos,
  RefProbe_Gain_SOX10neg = plot_gain_sox10neg,
  RefProbe_Loss_SOX10pos = plot_loss_sox10pos,
  RefProbe_Loss_SOX10neg = plot_loss_sox10neg
)

for (nm in names(plot_list)) {
  ggplot2::ggsave(
    filename = file.path(
      plots_ref_cnv_dir,
      paste0("FISH_reference_probe_CNVs_", nm, "_boxplot_with_brackets.png")
    ),
    plot = plot_list[[nm]],
    width = 7,
    height = 7,
    dpi = 600
  )
}

summary_table <- function(df) {
  df %>%
    dplyr::filter(!is.na(Percentage)) %>%
    dplyr::group_by(Event, CellType, Probe, Group_recode) %>%
    dplyr::summarise(
      n = dplyr::n_distinct(paste(Donor_ID, Region, sep = "_")),
      mean = mean(Percentage, na.rm = TRUE),
      sd = stats::sd(Percentage, na.rm = TRUE),
      min = min(Percentage, na.rm = TRUE),
      max = max(Percentage, na.rm = TRUE),
      median = stats::median(Percentage, na.rm = TRUE),
      IQR_low = stats::quantile(Percentage, 0.25, na.rm = TRUE),
      IQR_high = stats::quantile(Percentage, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}

summary_ref <- dplyr::bind_rows(
  summary_table(plot_gain_data),
  summary_table(plot_loss_data)
)

readr::write_tsv(
  summary_ref,
  file.path(
    stats_ref_cnv_dir,
    "FISH_reference_probe_CNV_analysis_SummaryStats_Gain_and_Loss_by_CellType_Probe_Group.tsv"
  ),
  na = ""
)

format_tests <- function(df, test_type) {
  
  if (!"Group_recode" %in% names(df)) {
    df <- df %>% dplyr::mutate(Group_recode = NA_character_)
  }
  
  df <- df %>% dplyr::mutate(Test = test_type)
  
  if (test_type == "between") {
    df <- df %>%
      dplyr::mutate(Comparison = paste0("Control vs MSA (", Probe, ")"))
  }
  
  if (test_type == "within") {
    df <- df %>%
      dplyr::mutate(
        Comparison = dplyr::if_else(
          is.na(Group_recode),
          "Reference vs SNCA",
          paste0("Reference vs SNCA (", Group_recode, ")")
        )
      )
  }
  
  df %>%
    dplyr::select(dplyr::any_of(c(
      "Event", "CellType", "Test", "Comparison",
      "p_raw", "p_adj", "p_label"
    )))
}

gain_between$Event <- "Gain"
gain_within$Event  <- "Gain"
loss_between$Event <- "Loss"
loss_within$Event  <- "Loss"

all_tests <- dplyr::bind_rows(
  format_tests(gain_between, "between"),
  format_tests(loss_between, "between"),
  format_tests(gain_within, "within"),
  format_tests(loss_within, "within")
)

readr::write_tsv(
  all_tests,
  file.path(
    stats_ref_cnv_dir,
    "FISH_reference_probe_CNV_analysis_WilcoxonResults_between_Group_and_within_Probe_Gain_and_Loss.tsv"
  ),
  na = ""
)
