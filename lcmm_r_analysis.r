# ============================================
# GROUP-BASED TRAJECTORY MODELING WITH LCMM
# Hospital visits following TB treatment
# ============================================

library(lcmm)
library(tidyverse)
library(nnet)      # For multinomial logistic regression
library(gtsummary) # For publication-ready tables
library(ggplot2)
library(patchwork)

# ============================================
# LOAD DATA
# ============================================

# Read the SAS dataset - adjust path as needed
library(haven)
lcmm_data <- read_sas("/path/to/HOME/lcmm_data.sas7bdat")

# Or if exported to CSV/Parquet
# lcmm_data <- read_csv("lcmm_data.csv")
# lcmm_data <- arrow::read_parquet("lcmm_data.parquet")

# Ensure proper variable types
lcmm_data <- lcmm_data %>%
  mutate(
    INDI_DSCM_NO = as.character(INDI_DSCM_NO),
    time = as.numeric(time),
    hospital_visits = as.integer(hospital_visits),
    # Categorical variables
    SEX_TYPE = as.factor(SEX_TYPE),
    EPISODE = as.factor(EPISODE),
    GAIBJA_TYPE = as.factor(GAIBJA_TYPE),
    RSLT_FN = as.factor(RSLT_FN),
    KITDIV = as.factor(KITDIV),
    SIDO = as.factor(SIDO),
    INC5 = as.factor(INC5)
  )

# ============================================
# PART 1: MODEL SELECTION (2-10 GROUPS)
# Testing different group solutions systematically
# ============================================

cat("========================================\n")
cat("PART 1: Testing 2-10 group solutions\n")
cat("========================================\n\n")

# Storage for model fit statistics
model_fit_results <- tibble(
  n_groups = integer(),
  BIC = numeric(),
  AIC = numeric(),
  loglik = numeric(),
  converged = logical(),
  n_params = integer()
)

# List to store all models
models_list <- list()

# Test 2-10 groups
for (n_grps in 2:10) {
  
  cat(paste0("\nFitting ", n_grps, "-group solution...\n"))
  
  # Set random seed for reproducibility
  set.seed(12345)
  
  # Fit model with Poisson distribution
  # Using quadratic time trajectories (can adjust polynomial order)
  tryCatch({
    
    model <- hlme(
      fixed = hospital_visits ~ time + I(time^2),  # Quadratic trajectory
      subject = "INDI_DSCM_NO",
      ng = n_grps,
      data = lcmm_data,
      link = "log",           # For Poisson-like count data
      idiag = TRUE,
      nwg = TRUE,
      maxiter = 500,
      verbose = FALSE
    )
    
    # Store model
    models_list[[paste0("grp", n_grps)]] <- model
    
    # Extract fit statistics
    model_fit_results <- model_fit_results %>%
      add_row(
        n_groups = n_grps,
        BIC = model$BIC,
        AIC = model$AIC,
        loglik = model$loglik,
        converged = model$conv == 1,
        n_params = model$N[2]
      )
    
    cat(paste0("  BIC: ", round(model$BIC, 2), 
               " | AIC: ", round(model$AIC, 2), 
               " | Converged: ", model$conv == 1, "\n"))
    
  }, error = function(e) {
    cat(paste0("  Error fitting ", n_grps, "-group model: ", e$message, "\n"))
  })
  
}

# Display model comparison
cat("\n========================================\n")
cat("Model Fit Comparison Table\n")
cat("========================================\n\n")

print(model_fit_results %>%
  mutate(
    delta_BIC = BIC - lag(BIC),
    delta_AIC = AIC - lag(AIC)
  ) %>%
  select(n_groups, BIC, delta_BIC, AIC, loglik, converged))

# Export for supplemental materials
write_csv(model_fit_results, "supplemental_table_model_fit.csv")

# ============================================
# PART 3: FIT FINAL MODEL
# Manually select optimal number of groups
# ============================================

# UPDATE THIS after reviewing model fit results
OPTIMAL_N <- 4

cat("\n========================================\n")
cat(paste0("PART 3: Running final ", OPTIMAL_N, "-group solution\n"))
cat("========================================\n\n")

# Fit final model
set.seed(12345)
final_model <- hlme(
  fixed = hospital_visits ~ time + I(time^2),
  subject = "INDI_DSCM_NO",
  ng = OPTIMAL_N,
  data = lcmm_data,
  link = "log",
  idiag = TRUE,
  nwg = TRUE,
  maxiter = 500,
  verbose = TRUE
)

# Model summary
summary(final_model)

# ============================================
# Calculate posterior probabilities and assign classes
# ============================================

# Get posterior probabilities
posterior_probs <- final_model$pprob %>%
  as_tibble() %>%
  rename(INDI_DSCM_NO = 1, trajectory_group = 2) %>%
  mutate(INDI_DSCM_NO = as.character(INDI_DSCM_NO))

# Calculate average posterior probability by group
avg_pp_by_group <- posterior_probs %>%
  group_by(trajectory_group) %>%
  summarise(
    n = n(),
    pct = 100 * n / nrow(posterior_probs),
    avg_pp = mean(3 + trajectory_group, na.rm = TRUE) # Get max PP column
  )

cat("\nAverage Posterior Probabilities by Group:\n")
print(avg_pp_by_group)

# Check if any groups < 5%
small_groups <- avg_pp_by_group %>%
  filter(pct < 5)

if (nrow(small_groups) > 0) {
  warning("Some groups have < 5% of sample. Consider reducing number of groups.")
}

# ============================================
# PART 4: TRAJECTORY VISUALIZATION
# ============================================

cat("\n========================================\n")
cat("PART 4: Creating trajectory visualizations\n")
cat("========================================\n\n")

# Predict trajectories for each group
pred_data <- expand_grid(
  time = seq(0, 20, by = 0.5),
  trajectory_group = 1:OPTIMAL_N
)

# Get predicted values for each class
predictions <- predictY(
  final_model, 
  newdata = pred_data,
  var.time = "time"
)

# Combine predictions with group info
traj_plot_data <- pred_data %>%
  mutate(
    predicted_visits = predictions$pred
  ) %>%
  left_join(
    avg_pp_by_group %>% select(trajectory_group, n, pct),
    by = "trajectory_group"
  )

# Assign descriptive labels based on visual inspection
# UPDATE THESE after viewing the plot
group_labels <- tibble(
  trajectory_group = 1:OPTIMAL_N,
  label = c(
    "Low-Stable",
    "Moderate-Increasing", 
    "High-Decreasing",
    "Persistently High"
  )[1:OPTIMAL_N]
)

# Add labels to plot data
traj_plot_data <- traj_plot_data %>%
  left_join(group_labels, by = "trajectory_group") %>%
  mutate(
    group_label = paste0(label, " (", round(pct, 1), "%)")
  )

# Create publication-ready trajectory plot
trajectory_plot <- ggplot(traj_plot_data, aes(x = time, y = predicted_visits, 
                                                color = group_label, 
                                                group = trajectory_group)) +
  geom_line(linewidth = 1.2) +
  scale_color_brewer(palette = "Set1", name = "Trajectory Group") +
  scale_x_continuous(
    breaks = seq(0, 20, by = 4),
    limits = c(0, 20)
  ) +
  labs(
    title = "Trajectory Groups of Hospital Visit Patterns",
    subtitle = "Following TB Treatment Completion",
    x = "Quarter Since Treatment Completion",
    y = "Predicted Number of Hospital Visits"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold")
  )

print(trajectory_plot)

# Save high-resolution figure
ggsave(
  "Figure1_Trajectory_Groups.tiff",
  plot = trajectory_plot,
  width = 10, height = 6, units = "in",
  dpi = 300, compression = "lzw"
)

ggsave(
  "Figure1_Trajectory_Groups.pdf",
  plot = trajectory_plot,
  width = 10, height = 6, units = "in"
)

# ============================================
# PART 5: BASELINE CHARACTERISTICS TABLE
# Comparing all trajectory groups
# ============================================

cat("\n========================================\n")
cat("PART 5: Baseline characteristics analysis\n")
cat("========================================\n\n")

# Merge trajectory assignments with baseline data
analysis_data <- lcmm_data %>%
  filter(time == 0) %>%  # Baseline only
  left_join(
    posterior_probs %>% select(INDI_DSCM_NO, trajectory_group),
    by = "INDI_DSCM_NO"
  ) %>%
  left_join(group_labels, by = "trajectory_group") %>%
  mutate(
    trajectory_group = factor(trajectory_group),
    trajectory_label = factor(label)
  )

# Create comprehensive Table 1 using gtsummary
table1 <- analysis_data %>%
  select(
    trajectory_label,
    SEX_TYPE, AGE_NHIS, EPISODE, GAIBJA_TYPE, 
    RSLT_FN, KITDIV, SIDO, INC5,
    hospital_visits
  ) %>%
  tbl_summary(
    by = trajectory_label,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      SEX_TYPE ~ "Sex",
      AGE_NHIS ~ "Age (years)",
      EPISODE ~ "Episode",
      GAIBJA_TYPE ~ "Insurance Type",
      RSLT_FN ~ "Treatment Success",
      KITDIV ~ "Treatment Round",
      SIDO ~ "Province",
      INC5 ~ "Income Quintile",
      hospital_visits ~ "Baseline Hospital Visits"
    ),
    missing = "no"
  ) %>%
  add_overall() %>%
  add_p(
    test = list(
      all_continuous() ~ "aov",  # ANOVA for continuous
      all_categorical() ~ "chisq.test"  # Chi-square for categorical
    )
  ) %>%
  modify_header(label ~ "**Characteristic**") %>%
  modify_caption("**Table 1. Baseline Characteristics by Trajectory Group**") %>%
  bold_labels()

print(table1)

# Save Table 1
table1 %>%
  as_gt() %>%
  gt::gtsave("Table1_Baseline_Characteristics.docx")

# Additional tests for non-normal continuous variables
# Kruskal-Wallis test if needed
if (TRUE) {  # Set to TRUE if you want non-parametric tests
  cat("\nKruskal-Wallis tests for continuous variables:\n")
  
  kw_age <- kruskal.test(AGE_NHIS ~ trajectory_group, data = analysis_data)
  cat(paste0("Age: H = ", round(kw_age$statistic, 2), 
             ", p = ", format.pval(kw_age$p.value, digits = 3), "\n"))
  
  kw_visits <- kruskal.test(hospital_visits ~ trajectory_group, data = analysis_data)
  cat(paste0("Baseline visits: H = ", round(kw_visits$statistic, 2), 
             ", p = ", format.pval(kw_visits$p.value, digits = 3), "\n"))
}

# ============================================
# PART 6: MULTINOMIAL LOGISTIC REGRESSION
# Predictors of trajectory membership
# ============================================

cat("\n========================================\n")
cat("PART 6: Multinomial logistic regression\n")
cat("========================================\n\n")

# Set reference group (e.g., Low-Stable = group 1)
REF_GROUP <- 1
analysis_data <- analysis_data %>%
  mutate(trajectory_group = relevel(trajectory_group, ref = as.character(REF_GROUP)))

# ---- Univariable analyses ----
cat("Running univariable analyses...\n\n")

univar_vars <- c("SEX_TYPE", "AGE_NHIS", "EPISODE", "GAIBJA_TYPE", 
                 "RSLT_FN", "KITDIV", "SIDO", "INC5")

univar_results <- list()

for (var in univar_vars) {
  cat(paste0("  ", var, "\n"))
  
  formula_str <- paste0("trajectory_group ~ ", var)
  model <- multinom(
    as.formula(formula_str),
    data = analysis_data,
    trace = FALSE
  )
  
  univar_results[[var]] <- model
}

# Extract univariable results
univar_summary <- univar_vars %>%
  map_dfr(function(var) {
    model <- univar_results[[var]]
    
    # Get coefficients and standard errors
    coef_table <- summary(model)$coefficients
    se_table <- summary(model)$standard.errors
    
    # Calculate z-values and p-values
    z_vals <- coef_table / se_table
    p_vals <- 2 * (1 - pnorm(abs(z_vals)))
    
    # Convert to odds ratios
    or_table <- exp(coef_table)
    or_lower <- exp(coef_table - 1.96 * se_table)
    or_upper <- exp(coef_table + 1.96 * se_table)
    
    # Format results
    tibble(
      variable = var,
      comparison = rownames(coef_table),
      term = colnames(coef_table),
      OR = as.vector(or_table),
      OR_lower = as.vector(or_lower),
      OR_upper = as.vector(or_upper),
      p_value = as.vector(p_vals)
    )
  })

# Save univariable results
write_csv(univar_summary, "univariable_results.csv")

# ---- Multivariable model ----
cat("\n\nRunning multivariable model...\n")

multivar_model <- multinom(
  trajectory_group ~ SEX_TYPE + AGE_NHIS + EPISODE + GAIBJA_TYPE + 
                     RSLT_FN + KITDIV + SIDO + INC5,
  data = analysis_data,
  trace = FALSE
)

summary(multivar_model)

# Extract multivariable results
coef_multi <- summary(multivar_model)$coefficients
se_multi <- summary(multivar_model)$standard.errors
z_multi <- coef_multi / se_multi
p_multi <- 2 * (1 - pnorm(abs(z_multi)))

# Create formatted results table
multivar_results <- tibble(
  comparison = rep(rownames(coef_multi), each = ncol(coef_multi)),
  term = rep(colnames(coef_multi), nrow(coef_multi)),
  OR = as.vector(exp(coef_multi)),
  OR_lower = as.vector(exp(coef_multi - 1.96 * se_multi)),
  OR_upper = as.vector(exp(coef_multi + 1.96 * se_multi)),
  p_value = as.vector(p_multi)
) %>%
  mutate(
    OR_CI = paste0(
      sprintf("%.2f", OR), " (",
      sprintf("%.2f", OR_lower), "-",
      sprintf("%.2f", OR_upper), ")"
    ),
    p_formatted = case_when(
      p_value < 0.001 ~ "<0.001",
      p_value < 0.01 ~ sprintf("%.3f", p_value),
      TRUE ~ sprintf("%.2f", p_value)
    )
  )

# Display and save multivariable results
cat("\nMultivariable Results (Odds Ratios):\n")
print(multivar_results %>%
  select(comparison, term, OR_CI, p_formatted) %>%
  filter(term != "(Intercept)"))

write_csv(multivar_results, "Table2_Multivariable_Results.csv")

# Create publication table using gtsummary
multivar_table <- analysis_data %>%
  select(trajectory_group, SEX_TYPE, AGE_NHIS, EPISODE, GAIBJA_TYPE,
         RSLT_FN, KITDIV, SIDO, INC5) %>%
  tbl_uvregression(
    method = nnet::multinom,
    y = trajectory_group,
    exponentiate = TRUE,
    hide_n = TRUE
  ) %>%
  modify_caption("**Table 2. Predictors of Trajectory Group Membership (Univariable)**")

print(multivar_table)

# ============================================
# PART 7: SUPPLEMENTAL MATERIALS
# Model diagnostics
# ============================================

cat("\n========================================\n")
cat("PART 7: Supplemental diagnostics\n")
cat("========================================\n\n")

# Group sizes and percentages
group_summary <- posterior_probs %>%
  group_by(trajectory_group) %>%
  summarise(
    n = n(),
    percentage = 100 * n / nrow(posterior_probs)
  ) %>%
  left_join(group_labels, by = "trajectory_group")

cat("\nSupplemental Table S1: Trajectory Group Sizes\n")
print(group_summary)

# Export model information
supplemental_info <- tibble(
  metric = c("Number of groups", "BIC", "AIC", "Log-likelihood", 
             "Number of parameters", "Converged"),
  value = c(
    OPTIMAL_N,
    round(final_model$BIC, 2),
    round(final_model$AIC, 2),
    round(final_model$loglik, 2),
    final_model$N[2],
    final_model$conv == 1
  )
)

cat("\nSupplemental Table S2: Final Model Information\n")
print(supplemental_info)

write_csv(supplemental_info, "supplemental_model_info.csv")

# ============================================
# Session info for reproducibility
# ============================================
cat("\n========================================\n")
cat("Session Information\n")
cat("========================================\n")
sessionInfo()

cat("\n========================================\n")
cat("Analysis complete!\n")
cat("========================================\n")