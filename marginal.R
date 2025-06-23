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
D0 <- rep(0, n)  # no one treated at baseline
threshold <- quantile(Y0, 0.5)  # treat top 30% at follow-up
D1 <- ifelse(Y0 > threshold, 1, 0)
#D1 <- rbinom(n,1,0.5)
# Y1 (follow-up - add drug effects for D1 = 1)
Y1 <- beta0 * G0 + betaE * Et + gammaE * GE * Et + betaD * D1+ + gammaD * GD * D1 + rnorm(n, 0, noise_sd)
Delta_Y <- Y1 - Y0  # Change from baseline
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

model1 <- lm(Y1 ~ G0, data = df)
summary(model1)
model2 <- lm(Delta_Y ~ G0, data = df)
summary(model2)