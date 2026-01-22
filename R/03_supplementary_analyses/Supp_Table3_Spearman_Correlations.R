# Spearman_Correlations_All.R
# Uses objects + paths created in:
#   R/00_setup/00_setup_paths.R
#   R/00_setup/01_setup_environment.R
# Required objects: cnv_data, avg_cnv_data_msa_only, avg_cnv_data_control_only, stats_spearman_dir

source(here::here("R", "00_setup", "01_setup_environment.R"))

library(dplyr)
library(tidyr)
library(readr)

dir.create(stats_spearman_dir, recursive = TRUE, showWarnings = FALSE)

format_p <- function(p_value) {
  p_value <- as.numeric(p_value)
  if (is.na(p_value)) return(NA_character_)
  
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
    if (rounded == 0) paste0("p = 0.001", star) else paste0("p = ", rounded, star)
  }
}

spearman_row <- function(df, x, y, group_label, var1_label, var2_label) {
  if (!all(c(x, y) %in% colnames(df))) {
    return(tibble(
      Variable_1 = var1_label, Variable_2 = var2_label, Group = group_label,
      n = NA_integer_, rho = NA_real_, p_value = NA_real_, p_display = NA_character_
    ))
  }
  
  d <- df %>%
    dplyr::select(dplyr::all_of(c(x, y))) %>%
    dplyr::filter(!is.na(.data[[x]]), !is.na(.data[[y]]))
  
  if (nrow(d) < 3) {
    return(tibble(
      Variable_1 = var1_label, Variable_2 = var2_label, Group = group_label,
      n = nrow(d), rho = NA_real_, p_value = NA_real_, p_display = NA_character_
    ))
  }
  
  ct <- suppressWarnings(stats::cor.test(d[[x]], d[[y]], method = "spearman", exact = FALSE))
  
  tibble(
    Variable_1 = var1_label,
    Variable_2 = var2_label,
    Group = group_label,
    n = nrow(d),
    rho = unname(ct$estimate),
    p_value = unname(ct$p.value),
    p_display = format_p(unname(ct$p.value))
  )
}

msa_df <- cnv_data %>% dplyr::filter(Classification %in% c("OPCA", "SND"))
control_df <- cnv_data %>% dplyr::filter(Classification == "Control")

make_ma_la <- function(df, value_col) {
  df %>%
    dplyr::filter(Region %in% c("Cerebellum", "Putamen")) %>%
    dplyr::select(
      Donor_ID, Classification, Region,
      !!rlang::sym(value_col),
      Percentage_SOX10pos_inclusions,
      Percentage_SOX10neg_inclusions
    ) %>%
    tidyr::pivot_wider(
      id_cols = c(Donor_ID, Classification),
      names_from = Region,
      values_from = c(
        !!rlang::sym(value_col),
        Percentage_SOX10pos_inclusions,
        Percentage_SOX10neg_inclusions
      ),
      names_sep = "__"
    ) %>%
    dplyr::mutate(
      MA_value = dplyr::if_else(
        Classification == "OPCA",
        .data[[paste0(value_col, "__Cerebellum")]],
        .data[[paste0(value_col, "__Putamen")]]
      ),
      LA_value = dplyr::if_else(
        Classification == "OPCA",
        .data[[paste0(value_col, "__Putamen")]],
        .data[[paste0(value_col, "__Cerebellum")]]
      ),
      MA_incl_pos = dplyr::if_else(
        Classification == "OPCA",
        .data[["Percentage_SOX10pos_inclusions__Cerebellum"]],
        .data[["Percentage_SOX10pos_inclusions__Putamen"]]
      ),
      MA_incl_neg = dplyr::if_else(
        Classification == "OPCA",
        .data[["Percentage_SOX10neg_inclusions__Cerebellum"]],
        .data[["Percentage_SOX10neg_inclusions__Putamen"]]
      )
    )
}

make_control_cer_put <- function(df, value_col) {
  df %>%
    dplyr::filter(Region %in% c("Cerebellum", "Putamen")) %>%
    dplyr::select(Donor_ID, Region, !!rlang::sym(value_col)) %>%
    tidyr::pivot_wider(
      id_cols = Donor_ID,
      names_from = Region,
      values_from = !!rlang::sym(value_col),
      names_sep = "__"
    )
}

msa_ma_la_pos_gain <- make_ma_la(msa_df, "Percentage_SOX10pos_gain")
msa_ma_la_pos_loss <- make_ma_la(msa_df, "Percentage_SOX10pos_loss")
msa_ma_la_neg_gain <- make_ma_la(msa_df, "Percentage_SOX10neg_gain")
msa_ma_la_neg_loss <- make_ma_la(msa_df, "Percentage_SOX10neg_loss")

control_cer_put_pos_gain <- make_control_cer_put(control_df, "Percentage_SOX10pos_gain")
control_cer_put_pos_loss <- make_control_cer_put(control_df, "Percentage_SOX10pos_loss")

out <- list()

out[["pos_gain_all_vs_neg_gain_all_msa"]] <- spearman_row(
  msa_df, "Percentage_SOX10pos_gain", "Percentage_SOX10neg_gain", "MSA",
  "% SOX10+ SNCA Gains (All)", "% SOX10- SNCA Gains (All)"
)

out[["pos_gain_all_vs_neg_gain_all_control"]] <- spearman_row(
  control_df, "Percentage_SOX10pos_gain", "Percentage_SOX10neg_gain", "Control",
  "% SOX10+ SNCA Gains (All)", "% SOX10- SNCA Gains (All)"
)

out[["pos_gain_ma_vs_la_msa"]] <- spearman_row(
  msa_ma_la_pos_gain, "MA_value", "LA_value", "MSA",
  "% SOX10+ SNCA Gains (MA)", "% SOX10+ SNCA Gains (LA)"
)

out[["pos_gain_cer_vs_put_control"]] <- spearman_row(
  control_cer_put_pos_gain, "Cerebellum", "Putamen", "Control",
  "% SOX10+ SNCA Gains (CER)", "% SOX10+ SNCA Gains (PUT)"
)

out[["pos_loss_all_vs_neg_loss_all_msa"]] <- spearman_row(
  msa_df, "Percentage_SOX10pos_loss", "Percentage_SOX10neg_loss", "MSA",
  "% SOX10+ SNCA Losses (All)", "% SOX10- SNCA Losses (All)"
)

out[["pos_loss_all_vs_neg_loss_all_control"]] <- spearman_row(
  control_df, "Percentage_SOX10pos_loss", "Percentage_SOX10neg_loss", "Control",
  "% SOX10+ SNCA Losses (All)", "% SOX10- SNCA Losses (All)"
)

out[["pos_loss_ma_vs_la_msa"]] <- spearman_row(
  msa_ma_la_pos_loss, "MA_value", "LA_value", "MSA",
  "% SOX10+ SNCA Losses (MA)", "% SOX10+ SNCA Losses (LA)"
)

out[["pos_loss_cer_vs_put_control"]] <- spearman_row(
  control_cer_put_pos_loss, "Cerebellum", "Putamen", "Control",
  "% SOX10+ SNCA Losses (CER)", "% SOX10+ SNCA Losses (PUT)"
)

out[["pos_gain_ma_vs_pos_incl_ma_msa"]] <- spearman_row(
  msa_ma_la_pos_gain, "MA_value", "MA_incl_pos", "MSA",
  "% SOX10+ SNCA Gains (MA)", "% SOX10+ alpha-syn inclusions (MA)"
)

out[["neg_gain_ma_vs_neg_incl_ma_msa"]] <- spearman_row(
  msa_ma_la_neg_gain, "MA_value", "MA_incl_neg", "MSA",
  "% SOX10- SNCA Gains (MA)", "% SOX10- alpha-syn inclusions (MA)"
)

out[["pos_loss_ma_vs_pos_incl_ma_msa"]] <- spearman_row(
  msa_ma_la_pos_loss, "MA_value", "MA_incl_pos", "MSA",
  "% SOX10+ SNCA Losses (MA)", "% SOX10+ alpha-syn inclusions (MA)"
)

out[["neg_loss_ma_vs_neg_incl_ma_msa"]] <- spearman_row(
  msa_ma_la_neg_loss, "MA_value", "MA_incl_neg", "MSA",
  "% SOX10- SNCA Losses (MA)", "% SOX10- alpha-syn inclusions (MA)"
)

avg_pairs <- tibble::tibble(
  avg_var = c(
    "Percentage_Average_SOX10pos_gain",
    "Percentage_Average_SOX10neg_gain",
    "Percentage_Average_SOX10pos_loss",
    "Percentage_Average_SOX10neg_loss"
  ),
  avg_label = c(
    "Average % SOX10+ SNCA Gains",
    "Average % SOX10- SNCA Gains",
    "Average % SOX10+ SNCA Losses",
    "Average % SOX10- SNCA Losses"
  )
)

clinical_pairs <- tibble::tibble(
  clin_var = c("Age_of_death", "Age_of_onset", "PMI", "Disease_duration"),
  clin_label = c(
    "Age of death (years)",
    "Age of onset (years)",
    "PMI (hours)",
    "Disease duration (years)"
  )
)

run_avg_block <- function(df, group_label) {
  rows <- list()
  k <- 1
  for (i in seq_len(nrow(avg_pairs))) {
    for (j in seq_len(nrow(clinical_pairs))) {
      rows[[k]] <- spearman_row(
        df,
        avg_pairs$avg_var[i],
        clinical_pairs$clin_var[j],
        group_label,
        avg_pairs$avg_label[i],
        clinical_pairs$clin_label[j]
      )
      k <- k + 1
    }
  }
  dplyr::bind_rows(rows)
}

out[["avg_msa_block"]] <- run_avg_block(avg_cnv_data_msa_only, "MSA")
out[["avg_control_block"]] <- run_avg_block(avg_cnv_data_control_only, "Control")

cor_tbl <- dplyr::bind_rows(out) %>%
  dplyr::arrange(Variable_1, Variable_2, Group)

readr::write_tsv(
  cor_tbl,
  file.path(stats_spearman_dir, "Spearman_Correlations_All.tsv"),
  na = ""
)

print(cor_tbl, n = Inf)
