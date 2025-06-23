library(tidyverse)
library(knitr)
library(xfun)
library(ggplot2)
library(patchwork) 


set.seed(234)

n <- 10000  
noise_sd <- 1  

# fixed across time
G0 <- sample(c(0, 1, 2), n, replace = TRUE, prob = c(0.64, 0.32, 0.04))  # MAF = 0.2
GE <- sample(c(0, 1, 2), n, replace = TRUE, prob = c(0.36, 0.48, 0.16)) # MAF = 0.4
GD <- sample(c(0, 1, 2), n, replace = TRUE, prob = c(0.49, 0.42, 0.09)) # MAF = 0.3

# environmental exposure (fixed across time)
Et <- rnorm(n, mean = 0, sd = 1)


# true effects
beta0 <- 0.3   # Baseline genetic effect
betaE <- 0.2   # Environmental effect
gammaE <- 0.4  # Gene-environment interaction
betaD <- 0.5   # Drug main effect
gammaD <- 0.6  # Pharmacogenetic effect

# Y0 (baseline - no drug effect since D0 = 0)
Y0 <- beta0 * G0 + betaE * Et + gammaE * GE * Et + rnorm(n, 0, noise_sd)


# treatment assignment (D0 = 0 for all and  D1 = 1 if Y0 > threshold)
D0 <- rep(0, n)  # No one treated at baseline
threshold <- quantile(Y0, 0.7)  # Treat top 30% at follow-up
D1 <- ifelse(Y0 > threshold, 1, 0)


# Y1 (follow-up - add drug effects for D1 = 1)
Y1 <- Y0 + betaD * D1 + gammaD * GD * D1 + rnorm(n, 0, noise_sd)
Delta_Y <- Y1 - Y0  # Change from baseline

# Combine into a dataframe
df <- data.frame(ID = 1:n, G0, GE, GD, Et, D0, D1, Y0, Y1, Delta_Y)


######################### Marginal Models ##################################

test_variant <- function(Gx, Gx_name) {
  formulas <- list(
    "Y1 ~ Gx"           = as.formula(paste("Y1 ~", Gx_name)),
    "Y1 ~ Y0 + Gx"      = as.formula(paste("Y1 ~ Y0 +", Gx_name)),
    "Delta_Y ~ Gx"      = as.formula(paste("Delta_Y ~", Gx_name)),
    "Delta_Y ~ Y0 + Gx" = as.formula(paste("Delta_Y ~ Y0 +", Gx_name))
  )
  models <- lapply(formulas, function(f) lm(f, data = df))
  results <- lapply(models, function(m) broom::tidy(m) %>% filter(term != "(Intercept)"))
  return(bind_rows(results, .id = "Model"))
}

# testing all variants
results <- bind_rows(
  mutate(test_variant(G0, "G0"), Variant = "G0"),
  mutate(test_variant(GE, "GE"), Variant = "GE"),
  mutate(test_variant(GD, "GD"), Variant = "GD")
) %>% select(Variant, Model, term, estimate, std.error, p.value)



### Results - Table

results %>%
  mutate(
    # Clean up model names only
    Model = case_when(
      Model == "Y1 ~ Gx" ~ "Y1 ~ G",
      Model == "Y1 ~ Y0 + Gx" ~ "Y1 ~ Y0 + G",
      Model == "Delta_Y ~ Gx" ~ "ΔY ~ G",
      Model == "Delta_Y ~ Y0 + Gx" ~ "ΔY ~ Y0 + G",
      TRUE ~ Model
    )
  ) %>%
  select(
    Variant,
    Model,
    Term = term,
    Estimate = estimate,
    StdError = std.error,
    Pvalue = p.value
  ) %>%
  arrange(Variant, Model) %>%
  as.data.frame()



      
### Results - Graphs


#### Combined


# 1. effect sizes with confidence intervals
effect_plot <- ggplot(results, aes(x = Model, y = estimate, color = Variant)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = estimate - 1.96*std.error, 
                    ymax = estimate + 1.96*std.error),
                width = 0.2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(title = "Effect Sizes Across Models",
       x = "Model", y = "Effect Size (β)", 
       color = "Variant") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 2. -log10(p-values)
pvalue_plot <- ggplot(results, aes(x = Model, y = -log10(p.value), fill = Variant)) +
  geom_col(position = position_dodge(), width = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Statistical Significance",
       x = "Model", y = "-log10(p-value)", 
       fill = "Variant") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# 3. pharmacogenetic effect visualization (GD vs Delta_Y in treated)
pharm_plot <- ggplot(df %>% filter(D1 == 1), 
                     aes(x = factor(GD), y = Delta_Y, fill = factor(GD))) +
  geom_boxplot(show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "red") +
  labs(title = "Pharmacogenetic Effect (Treated Individuals Only)",
       x = "GD Genotype", 
       y = "ΔY (Y1 - Y0)") +
  theme_minimal()


combined_plots <- (effect_plot | pvalue_plot) / pharm_plot +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold"))


print(combined_plots)
ggsave("association_results.png", combined_plots, width = 12, height = 10, dpi = 300)


effect_plot 
pvalue_plot 
pharm_plot
