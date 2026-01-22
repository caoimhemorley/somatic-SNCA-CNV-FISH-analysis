library(tidyverse)
library(here)
library(broom)
library(ggpubr)
library(car)
library(DescTools)
library(dunn.test)
library(moments)
library(FSA)
library(ggcorrplot)
library(viridis)
library(corrplot)
library(reshape2)
library(lme4)
library(janitor)
library(rstatix)
library(cowplot)
library(Hmisc)
library(patchwork)

source(here("R", "00_setup", "00_setup_paths.R"))

clean_column_names <- function(df) {
  nm <- colnames(df)
  nm <- str_replace_all(nm, "\\s+", "_")
  nm <- str_replace_all(nm, "%", "Percentage")
  nm <- str_replace_all(nm, "\\+", "pos")
  nm <- str_replace_all(nm, "\u2212", "-")
  nm <- str_replace_all(nm, "-", "neg")
  colnames(df) <- nm
  df
}

cnv_data <- read_tsv(file.path(data_raw_dir, "SNCA_CNV_Summarised_Percentages.tsv")) %>%
  clean_column_names()

cnv_data <- readr::read_tsv(file.path(data_raw_dir, "SNCA_CNV_Summarised_Percentages.tsv")) %>%
  clean_column_names() %>%
  dplyr::filter(!is.na(Total_cells))


case_info <- read_tsv(file.path(data_raw_dir, "Case_information.tsv")) %>%
  clean_column_names() %>%
  mutate(Sex = toupper(Sex))

ref_cnv_data <- read_tsv(file.path(data_raw_dir, "SNCA_Ref_CNV_Summary.tsv")) %>%
  clean_column_names()

gh2ax_data <- read_tsv(file.path(data_raw_dir, "Summarised_gh2ax_percentages.tsv")) %>%
  clean_column_names() %>%
  mutate(
    group_msa = case_when(
      Group %in% c("OPCA", "SND") ~ "MSA",
      Group == "Control" ~ "Control",
      TRUE ~ NA_character_
    )
  )

aSyn_NS_data <- read_tsv(file.path(data_raw_dir, "aSyn_NS_v_Tissue_sections.tsv")) %>%
  clean_column_names()

region_order <- c("Cerebellum", "Putamen", "SN")
classification_order <- c("Control", "SND", "OPCA")

cnv_data <- cnv_data %>%
  mutate(
    Region = factor(Region, levels = region_order),
    Classification = factor(Classification, levels = classification_order)
  )

avg_cnv_data <- cnv_data %>%
  group_by(Donor_ID, Classification) %>%
  summarise(
    across(starts_with("Percentage_"), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  rename_with(
    ~ paste0("Percentage_Average_", sub("^Percentage_", "", .x)),
    starts_with("Percentage_")
  ) %>%
  left_join(
    case_info %>% select(-Classification),
    by = "Donor_ID"
  )

msa_only_data <- cnv_data %>% filter(Classification != "Control")
control_only_data <- cnv_data %>% filter(Classification == "Control")

avg_cnv_data_msa_only <- avg_cnv_data %>% filter(Classification != "Control")
avg_cnv_data_control_only <- avg_cnv_data %>% filter(Classification == "Control")

putamen_data <- cnv_data %>% filter(Region == "Putamen")
cerebellum_data <- cnv_data %>% filter(Region == "Cerebellum")
sn_data <- cnv_data %>% filter(Region == "SN")

opca_data <- cnv_data %>% filter(Classification == "OPCA")
snd_data <- cnv_data %>% filter(Classification == "SND")
control_data <- cnv_data %>% filter(Classification == "Control")

most_affected_region <- cnv_data %>%
  filter(
    (Classification == "SND" & Region == "Putamen") |
      (Classification == "OPCA" & Region == "Cerebellum")
  )

least_affected_region <- cnv_data %>%
  filter(
    (Classification == "SND" & Region == "Cerebellum") |
      (Classification == "OPCA" & Region == "Putamen")
  )
