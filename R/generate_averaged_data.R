library(dplyr)

avg_cnv_data <- cnv_data %>%
  group_by(Donor_ID, Classification) %>%
  summarise(
    Percentage_Average_SOX10pos_gains =
      mean(Percentage_SOX10pos_gain, na.rm = TRUE),
    
    Percentage_Average_SOX10neg_gains =
      mean(PercentageSOX10neg_gain, na.rm = TRUE),
    
    Percentage_Average_SOX10pos_losses =
      mean(Percentage_SOX10pos_loss, na.rm = TRUE),
    
    Percentage_Average_SOX10neg_losses =
      mean(Percentage_SOX10neg_loss, na.rm = TRUE),
    
    Percentage_Average_SOX10pos_inclusions =
      mean(PercentageSOX10pos_inclusions, na.rm = TRUE),
    
    Percentage_Average_SOX10neg_inclusions =
      mean(Percentage_SOX10neg_inclusions, na.rm = TRUE),
    
    Percentage_Average_Overall_gain =
      mean(Percentage_Overall_gain, na.rm = TRUE),
    
    Percentage_Average_Overall_loss =
      mean(Percentage_Overall_loss, na.rm = TRUE),
    
    .groups = "drop"
  )
avg_cnv_data <- avg_cnv_data %>%
  left_join(case_info_clean, by = c("Donor_ID", "Classification"))
