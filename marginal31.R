### continuation of marginal3.R

summary_table <- all_results %>%
  filter(term %in% c("GD", "G0", "GE")) %>%
  group_by(term, GammaD, D1_Method, Model) %>%
  summarize(
    Mean_log10_p = mean(-log10(p.value)),
    SD_log10_p = sd(-log10(p.value)),
    Mean_p = mean(p.value),
    .groups = 'drop'
  ) %>%
  mutate(
    GammaD_Status = ifelse(GammaD == 0, "0 (Null)", "0.5 (Non-null)"),
    Model = factor(Model, 
                   levels = c("Y1 ~ G", "Y1 ~ Y0 + G", "ΔY ~ G", "ΔY ~ Y0 + G"),
                   labels = c("Y1 ~ G", "Y1 ~ Y0 + G", "ΔY ~ G", "ΔY ~ Y0 + G"))
  ) %>%
  select(Variant = term, GammaD_Status, D1_Method, Model, 
         Mean_log10_p, SD_log10_p, Mean_p) %>%
  arrange(Variant, GammaD_Status, D1_Method, Model)

# Custom formatting function for the table
format_table <- function(df) {
  df %>%
    mutate(
      Mean_log10_p = sprintf("%.3f", Mean_log10_p),
      SD_log10_p = sprintf("%.3f", SD_log10_p),
      Mean_p = format(Mean_p, scientific = TRUE, digits = 3)
    )
}


formatted_table <- summary_table %>% format_table()

knitr::kable(formatted_table, 
             caption = "Summary of Mean p-values Across Simulation Conditions",
             align = c("l", "l", "l", "l", "r", "r", "r"))


write_csv(summary_table, "simulation_results_summary_unrounded.csv")