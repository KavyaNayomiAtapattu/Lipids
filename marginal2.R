library(tidyverse)
library(ggplot2)
library(RColorBrewer)

set.seed(234)

simulate_data <- function(gammaD_value, d1_method = "threshold", threshold_quantile = 0.5) {
  n <- 10000
  noise_sd <- 1
  
  # Genetic variants
  G0 <- sample(c(0, 1, 2), n, replace = TRUE, prob = c(0.64, 0.32, 0.04))  # MAF = 0.2
  GE <- sample(c(0, 1, 2), n, replace = TRUE, prob = c(0.36, 0.48, 0.16))  # MAF = 0.4
  GD <- sample(c(0, 1, 2), n, replace = TRUE, prob = c(0.49, 0.42, 0.09))  # MAF = 0.3
  
  # Environmental exposure
  Et <- rnorm(n, mean = 0, sd = 1)
  
  # True effects
  beta0 <- 0.3   # Baseline genetic effect
  betaE <- 0.2   # Environmental effect
  gammaE <- 0.4  # Gene-environment interaction
  betaD <- 0.5   # Drug main effect
  
  # Baseline measure
  Y0 <- beta0 * G0 + betaE * Et + gammaE * GE * Et + rnorm(n, 0, noise_sd)
  
  # Treatment assignment
  D0 <- rep(0, n)
  if (d1_method == "threshold") {
    threshold <- quantile(Y0, threshold_quantile)
    D1 <- ifelse(Y0 > threshold, 1, 0)
  } else {
    D1 <- rbinom(n, 1, 0.5)  # 50% probability
  }
  
  # Follow-up measure
  Y1 <- beta0 * G0 + betaE * Et + gammaE * GE * Et + betaD * D1 + gammaD_value * GD * D1 + rnorm(n, 0, noise_sd)
  Delta_Y <- Y1 - Y0
  
  data.frame(ID = 1:n, G0, GE, GD, Et, D0, D1, Y0, Y1, Delta_Y)
}

run_models <- function(df, models_to_run = 1:4) {
  model_definitions <- list(
    "Y1 ~ G" = "Y1 ~ Gx",
    "Y1 ~ Y0 + G" = "Y1 ~ Y0 + Gx",
    "ΔY ~ G" = "Delta_Y ~ Gx",
    "ΔY ~ Y0 + G" = "Delta_Y ~ Y0 + Gx"
  )
  
  test_variant <- function(Gx, Gx_name) {
    results <- list()
    
    if (1 %in% models_to_run) {
      m1 <- lm(as.formula(gsub("Gx", Gx_name, model_definitions[[1]])), data = df)
      results[[1]] <- broom::tidy(m1) %>% filter(term != "(Intercept)")
    }
    
    if (2 %in% models_to_run) {
      m2 <- lm(as.formula(gsub("Gx", Gx_name, model_definitions[[2]])), data = df)
      results[[2]] <- broom::tidy(m2) %>% filter(term != "(Intercept)")
    }
    
    if (3 %in% models_to_run) {
      m3 <- lm(as.formula(gsub("Gx", Gx_name, model_definitions[[3]])), data = df)
      results[[3]] <- broom::tidy(m3) %>% filter(term != "(Intercept)")
    }
    
    if (4 %in% models_to_run) {
      m4 <- lm(as.formula(gsub("Gx", Gx_name, model_definitions[[4]])), data = df)
      results[[4]] <- broom::tidy(m4) %>% filter(term != "(Intercept)")
    }
    
    bind_rows(results, .id = "Model") %>% 
      mutate(Model = names(model_definitions)[as.numeric(Model)])
  }
  
  bind_rows(
    mutate(test_variant(G0, "G0"), Variant = "G0"),
    mutate(test_variant(GE, "GE"), Variant = "GE"),
    mutate(test_variant(GD, "GD"), Variant = "GD")
  ) %>% 
    select(Variant, Model, term, estimate, std.error, p.value)
}


run_iterations <- function(n_iter = 1000, gammaD_values = seq(0, 1, by = 0.1), 
                           d1_method = "threshold", models_to_run = 1:4) {
  map_dfr(gammaD_values, function(gammaD) {
    map_dfr(1:n_iter, function(i) {
      df <- simulate_data(gammaD, d1_method)
      results <- run_models(df, models_to_run)
      mutate(results, Iteration = i, GammaD = gammaD)
    })
  })
}




results <- run_iterations(
  n_iter = 1000,
  gammaD_values = c(0, 0.3, 0.6),  # Test these gammaD values
  d1_method = "threshold",         # or "binomial"
  models_to_run = c(1, 2)       # Run models 1, 2
)

# Plot results
ggplot(results %>% filter(Variant == "GD"), 
       aes(x = factor(GammaD), y = -log10(p.value))) +
  geom_boxplot(aes(fill = Model)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "forestgreen") +
  labs(title = "P-values across GammaD values for GD variant",
       x = "GammaD Value",
       y = "-log10(p-value)",
       fill = "Model") +
  theme_minimal() +
  facet_wrap(~Model, scales = "free_y")





ggplot(results %>% filter(Variant == "GD"), 
       aes(x = factor(GammaD), y = -log10(p.value))) +
  geom_boxplot(aes(fill = Model)) +
  scale_fill_brewer(palette = "Set2") +  # Try "Set1", "Set2", "Set3", "Paired", "Dark2"
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "forestgreen") +
  labs(title = "P-values across GammaD values for GD variant",
       x = "GammaD Value",
       y = "-log10(p-value)",
       fill = "Model") +
  theme_minimal() +
  facet_wrap(~Model, scales = "free_y")


