library(tidyverse)
library(ggpmisc)

set.seed(234)
simulate_drug_response <- function(
    n_individuals = 10000,      # Number of individuals
    maf_G0 = 0.2,               # MAF for baseline genetic effect
    maf_GE = 0.2,               # MAF for gene-environment SNP
    maf_GD = 0.2,               # MAF for pharmacogenetic SNP
    beta0 = 0.3,                # Baseline genetic effect
    betaE = 0.3,                # Environmental effect
    gammaE = 0.3,               # Gene-environment interaction effect
    betaD = 0.3,               # Drug effect
    gammaD_range = seq(0.1, 0.5, 0.1), # Range of pharmacogenetic effects to test
    drug_prop = 0.5,            # Proportion receiving drug
    noise_sd = 1,             # Standard deviation of noise
    show_plots = TRUE           # Whether to display plots
) {
  
  # genetic variants
  generate_snp <- function(n, maf) {
    sample(c(0,1,2), n, replace=TRUE, prob=c((1-maf)^2, 2*maf*(1-maf), maf^2)) 
  }
  
  G0 <- generate_snp(n_individuals, maf_G0)
  GE <- generate_snp(n_individuals, maf_GE)
  GD <- generate_snp(n_individuals, maf_GD)
  
  # environmental effect and drug exposure
  Et <- rnorm(n_individuals, mean = 0, sd = 1)
  Dt <- sample(c(0,1), n_individuals, replace = TRUE, prob = c(1-drug_prop, drug_prop))
  
  # genetic variance calculation
  calc_genetic_var <- function(maf, beta) {
    2 * maf * (1 - maf) * beta^2
  }
  
  # second moments calculations
  GE2 <- calc_genetic_var(maf_GE, gammaE) + mean(GE)^2  # E[GE^2]
  Et2 <- 1  # Since Et ~ N(0,1)
  Dt2 <- drug_prop  # For Bernoulli Dt
  
  # initialize results
  results <- data.frame()
  
  for (gammaD in gammaD_range) {
    GD2 <- calc_genetic_var(maf_GD, 1) + mean(GD)^2  # E[GD^2]
    
    # variance of product terms
    var_GE_Et <- GE2 * Et2 - (mean(GE)*mean(Et))^2
    var_GD_Dt <- GD2 * Dt2 - (mean(GD)*mean(Dt))^2
    
    # variance components
    var_components <- c(
      G0 = calc_genetic_var(maf_G0, beta0),
      Et = betaE^2 * var(Et),
      GE.Et = gammaE^2 * var_GE_Et,
      Dt = betaD^2 * var(Dt),
      GD.Dt = gammaD^2 * var_GD_Dt,
      Residuals = noise_sd^2
    )
    
    # % contributions
    var_percent <- 100 * var_components / sum(var_components)
    
    results <- rbind(results, data.frame(
      gammaD = gammaD,
      t(var_percent)
    ))
  }
  
  if (show_plots) {
    # data
    plot_data <- results %>% 
      pivot_longer(-gammaD, names_to = "Component", values_to = "Variance") %>%
      mutate(Component = factor(Component, 
                                levels = c("G0", "Et", "GE.Et", "Dt", "GD.Dt", "Residuals"),
                                labels = c("Baseline Genetic (β0*G0)", 
                                           "Environmental (βE*Et)",
                                           "G×E Interaction (γE*GE*Et)",
                                           "Drug Effect (βD*Dt)",
                                           "Pharmacogenetic (γD*GD*Dt)",
                                           "Residual Variance")))
    
    
    component_colors <- c(
      "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"
    )
    
    # Main plot
    p <- ggplot(plot_data, aes(x = gammaD, y = Variance, fill = Component)) +
      geom_area(position = "stack", alpha = 0.8) +
      geom_line(position = "stack", color = "white", size = 0.3) +
      scale_fill_manual(values = component_colors) +
      labs(
        title = "Variance Explained by Model Components",
        subtitle = sprintf("N = %s | MAF(GD) = %s | Drug Proportion = %s | Noise SD = %s",
                           format(n_individuals, big.mark = ","),
                           maf_GD,
                           drug_prop,
                           noise_sd),
        x = "Pharmacogenetic Effect Size (γD)",
        y = "Percentage of Total Variance",
        fill = "Model Component"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "bottom",
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "gray40"),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold")
      ) +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE))
    
    print(p)
    
    # pharmacogenetic effect
    p_focus <- ggplot(plot_data %>% filter(Component == "Pharmacogenetic (γD*GD*Dt)"), 
                      aes(x = gammaD, y = Variance)) +
      geom_line(color = "#9467bd", size = 1.5) +
      geom_point(color = "#9467bd", size = 3) +
      labs(
        title = "Pharmacogenetic Effect Variance",
        subtitle = "Contribution of GD*Dt interaction to total variance",
        x = "Pharmacogenetic Effect Size (γD)",
        y = "Percentage of Total Variance"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", color = "#9467bd"),
        panel.grid.major = element_line(color = "gray90")
      )
    
    print(p_focus)
  }
  
  return(list(
    variance_components = results,
    parameters = list(
      n_individuals = n_individuals,
      maf = c(G0 = maf_G0, GE = maf_GE, GD = maf_GD),
      effects = c(beta0 = beta0, betaE = betaE, gammaE = gammaE, 
                  betaD = betaD, gammaD_range = gammaD_range),
      proportions = c(drug_prop = drug_prop),
      noise = noise_sd
    )
  ))
  
  # Results table
  variance_table <- results %>%
    rename(
      `γD Value` = gammaD,
      `Baseline Genetic (%)` = G0,
      `Environmental (%)` = Et,
      `G×E Interaction (%)` = GE.Et,
      `Drug Effect (%)` = Dt,
      `Pharmacogenetic (%)` = GD.Dt,
      `Residual (%)` = Residuals
    ) %>%
    mutate(across(where(is.numeric), ~round(., 2)))
  
  # printing the table if plots are shown
  if (show_plots) {
    cat("\nVariance Explained by Each Component (%):\n")
    print(variance_table)
  }
  
  return(list(
    variance_components = variance_table,
    parameters = list(
      n_individuals = n_individuals,
      maf = c(G0 = maf_G0, GE = maf_GE, GD = maf_GD),
      effects = c(beta0 = beta0, betaE = betaE, gammaE = gammaE, 
                  betaD = betaD, gammaD_range = gammaD_range),
      proportions = c(drug_prop = drug_prop),
      noise = noise_sd
    ),
    plots = list(
      main_plot = p,
      focus_plot = p_focus
    )
  ))
}

results <- simulate_drug_response()

# variance table
variance_results <- results$variance_components
print(variance_results)

