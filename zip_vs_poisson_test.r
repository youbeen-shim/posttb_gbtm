# ============================================
# TESTING ZIP vs POISSON ASSUMPTION
# Does using Poisson instead of ZIP lead to erroneous conclusions?
# ============================================

library(tidyverse)
library(lcmm)
library(flexmix)     # Supports ZIP models
library(pscl)        # For zeroinfl function
library(haven)
library(aricode)     # For comparing cluster assignments (ARI, NMI)

# ============================================
# STEP 1: LOAD AND SAMPLE DATA
# ============================================

cat("========================================\n")
cat("Loading and sampling data\n")
cat("========================================\n\n")

# Load full dataset
lcmm_data <- read_sas("/path/to/HOME/lcmm_data.sas7bdat")
# Or: lcmm_data <- read_csv("lcmm_data.csv")

# Sample 8000 individuals for computational efficiency
set.seed(42)
sampled_patients <- lcmm_data %>%
  distinct(INDI_DSCM_NO) %>%
  slice_sample(n = 8000)

analysis_data <- lcmm_data %>%
  semi_join(sampled_patients, by = "INDI_DSCM_NO") %>%
  mutate(
    INDI_DSCM_NO = as.character(INDI_DSCM_NO),
    time = as.numeric(time),
    hospital_visits = as.integer(hospital_visits)
  )

cat(paste0("Sampled ", n_distinct(analysis_data$INDI_DSCM_NO), " patients\n"))
cat(paste0("Total observations: ", nrow(analysis_data), "\n\n"))

# ============================================
# STEP 2: ASSESS ZERO-INFLATION
# ============================================

cat("========================================\n")
cat("Assessing zero-inflation in the data\n")
cat("========================================\n\n")

zero_analysis <- analysis_data %>%
  summarise(
    total_obs = n(),
    n_zeros = sum(hospital_visits == 0),
    pct_zeros_observed = 100 * n_zeros / total_obs,
    mean_visits = mean(hospital_visits),
    var_visits = var(hospital_visits),
    # Expected zeros under Poisson distribution
    pct_zeros_poisson = 100 * exp(-mean_visits),
    # Variance-to-mean ratio (>1 suggests overdispersion)
    dispersion = var_visits / mean_visits
  )

print(zero_analysis)

cat("\nInterpretation:\n")
cat("- If observed zeros >> expected Poisson zeros, data is zero-inflated\n")
cat("- If dispersion >> 1, data is overdispersed\n")
cat("- Both suggest ZIP may be more appropriate than Poisson\n\n")

# Test for zero-inflation using Vuong test
# Fit simple models to demonstrate excess zeros
simple_poisson <- glm(hospital_visits ~ time + I(time^2), 
                      data = analysis_data, 
                      family = poisson)

simple_zip <- zeroinfl(hospital_visits ~ time + I(time^2) | 1, 
                       data = analysis_data,
                       dist = "poisson")

vuong_test <- vuong(simple_zip, simple_poisson)
cat("\nVuong test (ZIP vs Poisson):\n")
print(vuong_test)
cat("\nIf z-statistic > 1.96, ZIP significantly better than Poisson\n\n")

# ============================================
# STEP 3: FIT TRAJECTORY MODELS
# Test 3-5 group solutions for both distributions
# ============================================

cat("========================================\n")
cat("Fitting trajectory models\n")
cat("========================================\n\n")

# We'll test 3-5 groups for computational feasibility
n_groups_to_test <- 3:5

# Storage for results
model_comparison <- tibble()
poisson_models <- list()
zip_models <- list()

# ---- Fit Poisson models using lcmm ----
cat("Fitting Poisson models (lcmm)...\n")

for (k in n_groups_to_test) {
  cat(paste0("  ", k, " groups..."))
  
  set.seed(12345)
  model_pois <- tryCatch({
    hlme(
      fixed = hospital_visits ~ time + I(time^2),
      subject = "INDI_DSCM_NO",
      ng = k,
      data = analysis_data,
      link = "log",
      idiag = TRUE,
      nwg = TRUE,
      maxiter = 500,
      verbose = FALSE
    )
  }, error = function(e) NULL)
  
  if (!is.null(model_pois) && model_pois$conv == 1) {
    poisson_models[[paste0("k", k)]] <- model_pois
    
    model_comparison <- model_comparison %>%
      add_row(
        n_groups = k,
        distribution = "Poisson",
        BIC = model_pois$BIC,
        AIC = model_pois$AIC,
        loglik = model_pois$loglik,
        converged = TRUE
      )
    cat(" converged\n")
  } else {
    cat(" failed\n")
  }
}

# ---- Fit ZIP-like models using flexmix ----
# Note: flexmix handles mixture models differently
# We'll use FLXMRzipoisson for zero-inflated Poisson
cat("\nFitting ZIP models (flexmix)...\n")

# Prepare data for flexmix (needs wide or grouped format)
# We'll fit separate models and compare
for (k in n_groups_to_test) {
  cat(paste0("  ", k, " groups..."))
  
  set.seed(12345)
  model_zip <- tryCatch({
    flexmix(
      hospital_visits ~ time + I(time^2),
      data = analysis_data,
      k = k,
      model = FLXMRzipoisson(),  # Zero-inflated Poisson
      control = list(iter.max = 500, minprior = 0.05)
    )
  }, error = function(e) NULL)
  
  if (!is.null(model_zip)) {
    zip_models[[paste0("k", k)]] <- model_zip
    
    model_comparison <- model_comparison %>%
      add_row(
        n_groups = k,
        distribution = "ZIP",
        BIC = BIC(model_zip),
        AIC = AIC(model_zip),
        loglik = logLik(model_zip),
        converged = TRUE
      )
    cat(" converged\n")
  } else {
    cat(" failed\n")
  }
}

# ============================================
# STEP 4: COMPARE MODEL FIT
# ============================================

cat("\n========================================\n")
cat("Model Fit Comparison\n")
cat("========================================\n\n")

model_comparison_summary <- model_comparison %>%
  arrange(n_groups, distribution) %>%
  group_by(n_groups) %>%
  mutate(
    delta_BIC = BIC - BIC[distribution == "Poisson"],
    delta_AIC = AIC - AIC[distribution == "Poisson"]
  ) %>%
  ungroup()

print(model_comparison_summary)

cat("\nInterpretation:\n")
cat("- Negative delta_BIC/AIC means ZIP fits better than Poisson\n")
cat("- |delta_BIC| > 10 indicates strong evidence for better model\n")
cat("- |delta_BIC| > 6 indicates moderate evidence\n\n")

# Find best model for each distribution
best_poisson <- model_comparison %>%
  filter(distribution == "Poisson") %>%
  slice_min(BIC, n = 1)

best_zip <- model_comparison %>%
  filter(distribution == "ZIP") %>%
  slice_min(BIC, n = 1)

cat("Best Poisson model:", best_poisson$n_groups, "groups, BIC =", 
    round(best_poisson$BIC, 2), "\n")
cat("Best ZIP model:", best_zip$n_groups, "groups, BIC =", 
    round(best_zip$BIC, 2), "\n\n")

# ============================================
# STEP 5: COMPARE GROUP ASSIGNMENTS
# Do Poisson and ZIP assign people to different groups?
# ============================================

cat("========================================\n")
cat("Comparing trajectory group assignments\n")
cat("========================================\n\n")

# Select best model from each for comparison
# Or compare same K (e.g., 4 groups from each)
K_COMPARE <- 4  # Adjust based on your optimal solution

if (paste0("k", K_COMPARE) %in% names(poisson_models) && 
    paste0("k", K_COMPARE) %in% names(zip_models)) {
  
  # Get assignments from Poisson model
  pois_model <- poisson_models[[paste0("k", K_COMPARE)]]
  pois_assignments <- pois_model$pprob %>%
    as_tibble() %>%
    rename(INDI_DSCM_NO = 1, group_poisson = 2) %>%
    mutate(INDI_DSCM_NO = as.character(INDI_DSCM_NO))
  
  # Get assignments from ZIP model
  zip_model <- zip_models[[paste0("k", K_COMPARE)]]
  zip_assignments <- tibble(
    INDI_DSCM_NO = analysis_data %>% 
      distinct(INDI_DSCM_NO) %>% 
      arrange(INDI_DSCM_NO) %>% 
      pull(INDI_DSCM_NO),
    group_zip = clusters(zip_model)
  )
  
  # Merge assignments
  assignment_comparison <- pois_assignments %>%
    inner_join(zip_assignments, by = "INDI_DSCM_NO")
  
  # Calculate agreement metrics
  ari <- ARI(assignment_comparison$group_poisson, 
             assignment_comparison$group_zip)
  nmi <- NMI(assignment_comparison$group_poisson, 
             assignment_comparison$group_zip)
  
  # Crosstab
  crosstab <- table(
    Poisson = assignment_comparison$group_poisson,
    ZIP = assignment_comparison$group_zip
  )
  
  cat(paste0("Comparing ", K_COMPARE, "-group solutions\n\n"))
  cat("Agreement metrics:\n")
  cat(paste0("  Adjusted Rand Index (ARI): ", round(ari, 3), "\n"))
  cat(paste0("  Normalized Mutual Information (NMI): ", round(nmi, 3), "\n\n"))
  cat("  ARI/NMI interpretation:\n")
  cat("    1.0 = perfect agreement\n")
  cat("    0.0 = random agreement\n")
  cat("    <0.5 = poor agreement (different group structures)\n")
  cat("    0.5-0.7 = moderate agreement\n")
  cat("    >0.7 = good agreement\n\n")
  
  cat("Crosstabulation of group assignments:\n")
  print(crosstab)
  cat("\n")
  
  # Percentage agreement
  pct_agreement <- 100 * sum(assignment_comparison$group_poisson == 
                               assignment_comparison$group_zip) / 
                    nrow(assignment_comparison)
  cat(paste0("Exact agreement: ", round(pct_agreement, 1), "%\n\n"))
  
} else {
  cat("Could not compare - models did not converge for K =", K_COMPARE, "\n\n")
}

# ============================================
# STEP 6: COMPARE PREDICTED TRAJECTORIES
# Do the trajectory shapes differ?
# ============================================

cat("========================================\n")
cat("Comparing predicted trajectory shapes\n")
cat("========================================\n\n")

if (paste0("k", K_COMPARE) %in% names(poisson_models) && 
    paste0("k", K_COMPARE) %in% names(zip_models)) {
  
  # Predict trajectories for Poisson model
  pred_data <- expand_grid(
    time = seq(0, 20, by = 0.5),
    trajectory_group = 1:K_COMPARE
  )
  
  pois_pred <- predictY(pois_model, newdata = pred_data, var.time = "time")
  
  pois_traj <- pred_data %>%
    mutate(
      predicted_visits = pois_pred$pred,
      model = "Poisson"
    )
  
  # For ZIP model from flexmix, we need to extract parameters
  # and manually calculate predictions
  zip_params <- parameters(zip_model)
  
  # Calculate ZIP predictions for each component
  zip_traj <- pred_data %>%
    rowwise() %>%
    mutate(
      # Extract coefficients for this component
      intercept = zip_params[1, trajectory_group],
      beta_time = zip_params[2, trajectory_group],
      beta_time2 = zip_params[3, trajectory_group],
      # Linear predictor
      lambda = exp(intercept + beta_time * time + beta_time2 * time^2),
      # ZIP mean (simplified - actual calculation more complex)
      predicted_visits = lambda,
      model = "ZIP"
    ) %>%
    ungroup() %>%
    select(time, trajectory_group, predicted_visits, model)
  
  # Combine trajectories
  combined_traj <- bind_rows(pois_traj, zip_traj)
  
  # Plot comparison
  traj_comparison_plot <- ggplot(combined_traj, 
                                  aes(x = time, y = predicted_visits,
                                      color = model, linetype = model)) +
    geom_line(linewidth = 1) +
    facet_wrap(~ paste("Group", trajectory_group), ncol = 2, scales = "free_y") +
    scale_color_manual(values = c("Poisson" = "blue", "ZIP" = "red")) +
    labs(
      title = paste0("Trajectory Comparison: Poisson vs ZIP (", K_COMPARE, " groups)"),
      x = "Quarter Since Treatment Completion",
      y = "Predicted Hospital Visits",
      color = "Model",
      linetype = "Model"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  print(traj_comparison_plot)
  
  ggsave(
    "trajectory_comparison_poisson_vs_zip.pdf",
    plot = traj_comparison_plot,
    width = 10, height = 8
  )
  
  # Calculate RMSE between trajectories
  traj_rmse <- combined_traj %>%
    select(time, trajectory_group, model, predicted_visits) %>%
    pivot_wider(names_from = model, values_from = predicted_visits) %>%
    group_by(trajectory_group) %>%
    summarise(
      rmse = sqrt(mean((Poisson - ZIP)^2)),
      max_diff = max(abs(Poisson - ZIP)),
      mean_diff = mean(Poisson - ZIP)
    )
  
  cat("Trajectory shape differences (RMSE by group):\n")
  print(traj_rmse)
  cat("\n")
  
} else {
  cat("Could not compare trajectories - models not available\n\n")
}

# ============================================
# STEP 7: CLINICAL IMPLICATIONS
# ============================================

cat("========================================\n")
cat("Summary and Clinical Implications\n")
cat("========================================\n\n")

cat("KEY FINDINGS:\n\n")

cat("1. Zero-inflation assessment:\n")
cat(paste0("   - Observed zeros: ", round(zero_analysis$pct_zeros_observed, 1), "%\n"))
cat(paste0("   - Expected under Poisson: ", round(zero_analysis$pct_zeros_poisson, 1), "%\n"))
cat(paste0("   - Excess zeros: ", 
    round(zero_analysis$pct_zeros_observed - zero_analysis$pct_zeros_poisson, 1), 
    " percentage points\n\n"))

cat("2. Model fit comparison:\n")
if (nrow(model_comparison) > 0) {
  bic_diff <- best_zip$BIC - best_poisson$BIC
  cat(paste0("   - BIC difference (ZIP - Poisson): ", round(bic_diff, 2), "\n"))
  cat("   - Interpretation: ")
  if (abs(bic_diff) < 2) {
    cat("Negligible difference - Poisson adequate\n")
  } else if (bic_diff < -10) {
    cat("ZIP substantially better - Poisson may be inadequate\n")
  } else if (bic_diff < -6) {
    cat("ZIP moderately better - consider using ZIP\n")
  } else if (bic_diff > 10) {
    cat("Poisson substantially better (unusual for zero-inflated data)\n")
  } else {
    cat("Small difference - practical significance unclear\n")
  }
  cat("\n")
}

if (exists("ari")) {
  cat("3. Group assignment agreement:\n")
  cat(paste0("   - ARI: ", round(ari, 3), "\n"))
  cat("   - Interpretation: ")
  if (ari > 0.8) {
    cat("Excellent agreement - Poisson and ZIP identify same groups\n")
  } else if (ari > 0.6) {
    cat("Good agreement - conclusions likely similar\n")
  } else if (ari > 0.4) {
    cat("Moderate agreement - some differences in group structure\n")
  } else {
    cat("Poor agreement - Poisson may identify different trajectory patterns\n")
  }
  cat("\n")
}

if (exists("traj_rmse")) {
  cat("\n4. Trajectory shape differences:\n")
  cat(paste0("   - Average RMSE: ", round(mean(traj_rmse$rmse), 2), " visits\n"))
  cat(paste0("   - Max difference: ", round(max(traj_rmse$max_diff), 2), " visits\n"))
  cat("   - Interpretation: ")
  if (max(traj_rmse$max_diff) < 0.5) {
    cat("Minimal difference in trajectory shapes\n")
  } else if (max(traj_rmse$max_diff) < 2) {
    cat("Small difference - unlikely to affect clinical interpretation\n")
  } else {
    cat("Substantial difference - trajectory shapes may differ meaningfully\n")
  }
  cat("\n")
}

cat("\n\nRECOMMENDATION:\n")
cat("Based on these analyses, using Poisson instead of ZIP is:\n")

# Decision logic
if (exists("bic_diff") && exists("ari")) {
  if (abs(bic_diff) < 6 && ari > 0.7) {
    cat("✓ ACCEPTABLE - Models show similar fit and group structure\n")
    cat("  Continue with Poisson for simplicity\n")
  } else if (bic_diff < -10 || ari < 0.5) {
    cat("✗ PROBLEMATIC - ZIP shows substantially better fit or different groups\n")
    cat("  Consider alternative approaches or acknowledge limitation\n")
  } else {
    cat("⚠ BORDERLINE - Some differences but may not affect main conclusions\n")
    cat("  Proceed with Poisson but report this sensitivity analysis\n")
  }
} else {
  cat("Unable to make recommendation - review individual components above\n")
}

cat("\n========================================\n")
cat("Analysis complete!\n")
cat("========================================\n")