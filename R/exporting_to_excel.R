# Export_Figure_Input_Data_Exact.R
# Exports "exact values used" for figures into one Excel workbook (multiple tabs)
# Assumes setup_environment.R defines:
# cnv_data, avg_cnv_data, gh2ax_data, aSyn_NS_data, etc.

library(dplyr)
library(tidyr)
library(writexl)
library(here)

source(here("R", "setup_environment.R"))

output_excel <- here("outputs", "Figure_input_data.xlsx")

keep_existing <- function(df, cols) cols[cols %in% colnames(df)]

# Helper: Spearman result for ONE pair (keeps donor IDs for the exact rows used)
spearman_one_pair <- function(df, x, y, id_cols) {
  need_cols <- unique(c(id_cols, x, y))
  need_cols <- keep_existing(df, need_cols)
  
  if (!all(c(x, y) %in% need_cols)) {
    missing <- setdiff(c(x, y), need_cols)
    return(list(
      input = tibble(note = paste0("Missing columns: ", paste(missing, collapse = ", "))),
      result = tibble(x = x, y = y, n = NA_integer_, rho = NA_real_, p_value = NA_real_)
    ))
  }
  
  input_df <- df %>%
    select(all_of(need_cols)) %>%
    drop_na(all_of(c(x, y)))
  
  if (nrow(input_df) < 3 || length(unique(input_df[[x]])) < 2 || length(unique(input_df[[y]])) < 2) {
    return(list(
      input = input_df,
      result = tibble(x = x, y = y, n = nrow(input_df), rho = NA_real_, p_value = NA_real_)
    ))
  }
  
  ct <- suppressWarnings(cor.test(input_df[[x]], input_df[[y]], method = "spearman", exact = FALSE))
  
  list(
    input = input_df,
    result = tibble(
      x = x,
      y = y,
      n = nrow(input_df),
      rho = unname(ct$estimate),
      p_value = ct$p.value
    )
  )
}

# 1) CNV AVERAGED (MSA vs Control) – 4 separate tabs

avg_sheet <- function(df, value_col, sheet_label) {
  if (!value_col %in% colnames(df)) {
    return(tibble(note = paste0("Missing column in avg_cnv_data: ", value_col)))
  }
  
  df %>%
    mutate(
      Classification = ifelse(Classification %in% c("SND", "OPCA"), "MSA", Classification),
      Classification = factor(Classification, levels = c("Control", "MSA"))
    ) %>%
    select(Donor_ID, Classification, all_of(value_col)) %>%
    rename(Percentage = all_of(value_col)) %>%
    mutate(CNV_Label = sheet_label) %>%
    relocate(CNV_Label, .before = Donor_ID)
}

avg_sox10pos_gain <- avg_sheet(avg_cnv_data, "Percentage_Average_SOX10pos_gain", "Avg SOX10+ Gains")
avg_sox10neg_gain <- avg_sheet(avg_cnv_data, "Percentage_Average_SOX10neg_gain", "Avg SOX10- Gains")
avg_sox10pos_loss <- avg_sheet(avg_cnv_data, "Percentage_Average_SOX10pos_loss", "Avg SOX10+ Losses")
avg_sox10neg_loss <- avg_sheet(avg_cnv_data, "Percentage_Average_SOX10neg_loss", "Avg SOX10- Losses")

# 2) CNV REGIONAL – 4 separate tabs (Cer/Put/SN across Control/SND/OPCA)

regional_sheet <- function(df, value_col, sheet_label) {
  if (!value_col %in% colnames(df)) {
    return(tibble(note = paste0("Missing column in cnv_data: ", value_col)))
  }
  
  df %>%
    select(Donor_ID, Classification, Region, all_of(value_col)) %>%
    rename(Percentage = all_of(value_col)) %>%
    mutate(CNV_Label = sheet_label) %>%
    relocate(CNV_Label, .before = Donor_ID)
}

regional_sox10pos_gain <- regional_sheet(cnv_data, "Percentage_SOX10pos_gain", "Regional SOX10+ Gains")
regional_sox10neg_gain <- regional_sheet(cnv_data, "Percentage_SOX10neg_gain", "Regional SOX10- Gains")
regional_sox10pos_loss <- regional_sheet(cnv_data, "Percentage_SOX10pos_loss", "Regional SOX10+ Losses")
regional_sox10neg_loss <- regional_sheet(cnv_data, "Percentage_SOX10neg_loss", "Regional SOX10- Losses")

# 3) GH2AX – exact inputs used by your GH2AX plots

gh2ax_all_cells <- gh2ax_data %>%
  filter(group_msa %in% c("MSA", "Control")) %>%
  select(Donor, Region, group_msa, Percentage_all_gh2ax_foci) %>%
  drop_na(Percentage_all_gh2ax_foci)

gh2ax_sox10_pos <- gh2ax_data %>%
  filter(group_msa %in% c("MSA", "Control")) %>%
  select(Donor, Region, group_msa, Percentage_SOX10_pos_gh2ax_pos) %>%
  drop_na(Percentage_SOX10_pos_gh2ax_pos)

gh2ax_sox10_neg <- gh2ax_data %>%
  filter(group_msa %in% c("MSA", "Control")) %>%
  select(Donor, Region, group_msa, Percentage_SOX10_neg_gh2ax_pos) %>%
  drop_na(Percentage_SOX10_neg_gh2ax_pos)

gh2ax_asyn_paired <- gh2ax_data %>%
  filter(group_msa == "MSA") %>%
  select(Donor, Region, Percentage_aSyn_pos_gh2ax_pos, Percentage_aSyn_neg_gh2ax_pos) %>%
  drop_na()

gh2ax_affected_regions <- gh2ax_data %>%
  filter(Group %in% c("OPCA", "SND"),
         Region %in% c("Cerebellum", "Putamen")) %>%
  mutate(
    Affectedness = case_when(
      Group == "OPCA" & Region == "Cerebellum" ~ "Most affected",
      Group == "SND"  & Region == "Putamen"    ~ "Most affected",
      Group == "OPCA" & Region == "Putamen"    ~ "Least affected",
      Group == "SND"  & Region == "Cerebellum" ~ "Least affected",
      TRUE ~ NA_character_
    )
  ) %>%
  select(Donor, Group, Region, Affectedness, Percentage_all_gh2ax_foci) %>%
  drop_na(Affectedness, Percentage_all_gh2ax_foci)

# 4) aSyn NS vs Sections – paired input (Donor_ID x Region)

asyn_ns_vs_sections <- aSyn_NS_data %>%
  select(Donor_ID, Region, Percent_inclusions_sections, Percent_inclusions_NS) %>%
  pivot_longer(
    cols = c(Percent_inclusions_sections, Percent_inclusions_NS),
    names_to = "Tissue",
    values_to = "Percent_inclusions"
  ) %>%
  mutate(
    Tissue = case_when(
      Tissue == "Percent_inclusions_sections" ~ "Sections",
      Tissue == "Percent_inclusions_NS" ~ "NS",
      TRUE ~ Tissue
    )
  ) %>%
  drop_na(Percent_inclusions)

# 5) GCI distribution – Putamen vs Cerebellum (SND/OPCA only)
# Correct inclusions name per your note: Percentage_SOX10pos_inclusions

gci_put_cb <- cnv_data %>%
  filter(
    Classification %in% c("SND", "OPCA"),
    Region %in% c("Putamen", "Cerebellum")
  ) %>%
  select(Donor_ID, Classification, Region, Percentage_SOX10pos_inclusions) %>%
  drop_na(Percentage_SOX10pos_inclusions)

# 6) Spearman – ONLY the specific figure correlations you listed
# Each correlation gets:
# - an input tab with Donor_ID + the exact x/y used
# - a results tab with rho/p/n

# 6.1) Avg SOX10+ gains vs Age_of_onset (MSA only)
sp1 <- spearman_one_pair(
  df = avg_cnv_data %>% filter(Classification %in% c("SND", "OPCA")),
  x = "Percentage_Average_SOX10pos_gain",
  y = "Age_of_onset",
  id_cols = c("Donor_ID", "Classification")
)

# 6.2) Avg SOX10- gains vs Age_of_onset (MSA only)
# NOTE: if your avg column is missing, this will return a "Missing columns" note sheet.
sp2 <- spearman_one_pair(
  df = avg_cnv_data %>% filter(Classification %in% c("SND", "OPCA")),
  x = "Percentage_Average_SOX10neg_gain",
  y = "Age_of_onset",
  id_cols = c("Donor_ID", "Classification")
)

# 6.3) SOX10+ gains vs SOX10- gains (overall across all regions analysed, MSA only)
# This uses REGIONAL cnv_data (each Donor_ID x Region is an observation),
# keeping Donor_ID and Region so the exact rows are traceable.
sp3 <- spearman_one_pair(
  df = cnv_data %>% filter(Classification %in% c("SND", "OPCA")),
  x = "Percentage_SOX10pos_gain",
  y = "Percentage_SOX10neg_gain",
  id_cols = c("Donor_ID", "Classification", "Region")
)

# 6.4) % SOX10+ inclusions vs % SOX10+ gains for MOST affected region
# Most affected definition: SND Putamen, OPCA Cerebellum
most_affected_cnv <- cnv_data %>%
  filter(
    (Classification == "SND" & Region == "Putamen") |
      (Classification == "OPCA" & Region == "Cerebellum")
  )

sp4 <- spearman_one_pair(
  df = most_affected_cnv,
  x = "Percentage_SOX10pos_inclusions",
  y = "Percentage_SOX10pos_gain",
  id_cols = c("Donor_ID", "Classification", "Region")
)

# WRITE WORKBOOK

write_xlsx(
  list(
    "Avg_SOX10pos_Gains" = avg_sox10pos_gain,
    "Avg_SOX10neg_Gains" = avg_sox10neg_gain,
    "Avg_SOX10pos_Losses" = avg_sox10pos_loss,
    "Avg_SOX10neg_Losses" = avg_sox10neg_loss,
    
    "Regional_SOX10pos_Gains" = regional_sox10pos_gain,
    "Regional_SOX10neg_Gains" = regional_sox10neg_gain,
    "Regional_SOX10pos_Losses" = regional_sox10pos_loss,
    "Regional_SOX10neg_Losses" = regional_sox10neg_loss,
    
    "GH2AX_all_cells" = gh2ax_all_cells,
    "GH2AX_SOX10_pos" = gh2ax_sox10_pos,
    "GH2AX_SOX10_neg" = gh2ax_sox10_neg,
    "GH2AX_aSyn_paired" = gh2ax_asyn_paired,
    "GH2AX_affected_regions" = gh2ax_affected_regions,
    
    "aSyn_NS_vs_Sections" = asyn_ns_vs_sections,
    "GCI_Putamen_vs_Cerebellum" = gci_put_cb,
    
    "Spearman_Input_Avg_SOX10posGain_vs_Onset_MSA" = sp1$input,
    "Spearman_Result_Avg_SOX10posGain_vs_Onset_MSA" = sp1$result,
    
    "Spearman_Input_Avg_SOX10negGain_vs_Onset_MSA" = sp2$input,
    "Spearman_Result_Avg_SOX10negGain_vs_Onset_MSA" = sp2$result,
    
    "Spearman_Input_Regional_SOX10posGain_vs_SOX10negGain_MSA" = sp3$input,
    "Spearman_Result_Regional_SOX10posGain_vs_SOX10negGain_MSA" = sp3$result,
    
    "Spearman_Input_MostAffected_Inclusions_vs_Gains" = sp4$input,
    "Spearman_Result_MostAffected_Inclusions_vs_Gains" = sp4$result
  ),
  path = output_excel
)

message("Excel file written to:")
message(output_excel)
