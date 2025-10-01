library(tidyverse)
library(nonprobsvy)
library(survey)
library(ggplot2)

set.seed(0)

betas <- list(beta_0 = 5, beta_age = 0.15, beta_age2 = -0.002, beta_sex = -0.5, beta_bmi = -0.05, beta_sex_bmi2 = 0.01)
gammas <- list(gamma_0 = -2.5, gamma_1 = -0.05, gamma_2 = -0.1)
probability_of_female <-  0.51
bmi_gamma_distibution <- list(shape=4.5, scale = 2.2, shift = 14)
age_beta_distribution <- list(a=2, b=3, min_age = 18, max_age=80)
avg_bmi <- 23
noise_time_sport_sd <- 3

generation_parameters <- list(
  betas = betas, 
  probability_of_female = probability_of_female, 
  bmi_distribution = bmi_gamma_distibution,
  age_distribution = age_beta_distribution,
  avg_bmi = avg_bmi,
  noise_time_sport_sd = noise_time_sport_sd
)

N <- 10000
n_np <- 1000
n_p <- 2000

generate_data <- function(N, generation_parameters){
  population <- tibble(
    index = 1:N,
    sex = rbinom(N, 1, generation_parameters$probability_of_female),
    bmi = rgamma(N, generation_parameters$bmi_distribution$shape, scale=generation_parameters$bmi_distribution$scale)+generation_parameters$bmi_distribution$shift,
    age = rbeta(N, generation_parameters$age_distribution$a, generation_parameters$age_distribution$b)*(generation_parameters$age_distribution$max_age-generation_parameters$age_distribution$min_age)+generation_parameters$age_distribution$min_age,
    base_time_sport = generation_parameters$betas$beta_0 +
      generation_parameters$betas$beta_age * age +
      generation_parameters$betas$beta_age2 * age^2 +
      generation_parameters$betas$beta_sex * sex +
      generation_parameters$betas$beta_bmi * (bmi - generation_parameters$avg_bmi)^2 +
      generation_parameters$betas$beta_sex_bmi2 * sex * (bmi - generation_parameters$avg_bmi)^2,
    time_sport = pmax(0, base_time_sport + rnorm(N, mean = 0, sd = generation_parameters$noise_time_sport_sd))
  )
  true_mean <- mean(population$time_sport)
  return(list(population=population, true_mean=true_mean))
}

draw_samples <- function(population, gammas, N, n_np, n_p, avg_bmi){

  population <- population |> 
    mutate(
      logit_pi = gammas$gamma_0 + gammas$gamma_1 * age + gammas$gamma_2 * (bmi-avg_bmi)^2,
      pi_np = 1 / (1 + exp(-logit_pi))
    )
  
  
  non_probability_sample <- population  |> 
    slice_sample(n = n_np, weight_by = pi_np) |> 
    select(!c(base_time_sport, logit_pi, pi_np))
  
  naive_mean <- mean(non_probability_sample$time_sport)
  
  probability_sample <- population |>
    anti_join(non_probability_sample, by='index') |> 
    slice_sample(n = n_p) |> 
    select(!c(base_time_sport, time_sport, logit_pi, pi_np)) |> 
    mutate(
      weight = N/n_p,
      population_size = N
    )

  
  probability_sample_svy <- svydesign(
    ids = ~1,
    data = probability_sample,
    weights = ~weight,
    fpc = ~population_size
  )
  
  return(list(
    non_probability_sample=non_probability_sample,
    probability_sample=probability_sample_svy,
    naive_mean=naive_mean
  ))
}

run_estimators <- function(non_probability_sample, probability_sample){

  outcome_formula = time_sport ~ sex + bmi + age
  selection_formula = ~ sex + bmi + age
  target_formula = ~ time_sport
  
  mi_glm <- nonprob(
    data=non_probability_sample, 
    outcome = outcome_formula,
    svydesign = probability_sample,
    method_outcome = 'glm',
    family_outcome = 'gaussian'
    )
  
  mi_glm_var_sel <- nonprob(
    data=non_probability_sample, 
    outcome = outcome_formula,
    svydesign = probability_sample,
    method_outcome = 'glm',
    family_outcome = 'gaussian',
    control_outcome = control_out(penalty = "lasso"), 
    control_inference = control_inf(vars_selection = TRUE)
  )
  
  mi_npar <- nonprob(
    data=non_probability_sample, 
    outcome = outcome_formula,
    svydesign = probability_sample,
    method_outcome = 'npar',
  )
  
  mi_nn_2 <- nonprob(
    data=non_probability_sample, 
    outcome = outcome_formula,
    svydesign = probability_sample,
    method_outcome = 'nn',
    control_outcome = control_out(k=2)
  )
  
  mi_nn_5 <- nonprob(
    data=non_probability_sample, 
    outcome = outcome_formula,
    svydesign = probability_sample,
    method_outcome = 'nn',
    control_outcome = control_out(k=5)
  )
  
  mi_nn_10 <- nonprob(
    data=non_probability_sample, 
    outcome = outcome_formula,
    svydesign = probability_sample,
    method_outcome = 'nn',
    control_outcome = control_out(k=10)
  )
  
  mi_pmm <- nonprob(
    data=non_probability_sample, 
    outcome = outcome_formula,
    svydesign = probability_sample,
    method_outcome = 'pmm',
    family_outcome = 'gaussian'
  )
  
  mi_pmm_2 <- nonprob(
    data=non_probability_sample, 
    outcome = outcome_formula,
    svydesign = probability_sample,
    method_outcome = 'pmm',
    family_outcome = 'gaussian',
    control_outcome = control_out(pmm_match_type = 2)
  )
  
  ipw <- nonprob(
    data=non_probability_sample, 
    selection = selection_formula,
    target = target_formula,
    svydesign = probability_sample,
    method_selection = "logit"
  )
  
  ipw_gee_1 <- nonprob(
    data=non_probability_sample, 
    selection = selection_formula,
    target = target_formula,
    svydesign = probability_sample,
    method_selection = "logit",
    control_selection = control_sel(est_method = "gee", gee_h_fun = 1)
  )
  
  ipw_gee_2 <- nonprob(
    data=non_probability_sample, 
    selection = selection_formula,
    target = target_formula,
    svydesign = probability_sample,
    method_selection = "logit",
    control_selection = control_sel(est_method = "gee", gee_h_fun = 2)
  )
  
  ipw_gee_var_sel <- nonprob(
    data=non_probability_sample, 
    selection = selection_formula,
    target = target_formula,
    svydesign = probability_sample,
    method_selection = "logit",
    control_selection = control_sel(est_method = "gee", gee_h_fun = 1),
    control_inference = control_inf(vars_selection = TRUE)
  )
  
  dr_glm_mle <- nonprob(
    data=non_probability_sample,
    outcome = outcome_formula,
    selection = selection_formula,
    svydesign = probability_sample,
    method_outcome = 'glm',
    family_outcome = 'gaussian',
    method_selection = "logit"
  )
  
  dr_glm_gee <- nonprob(
    data=non_probability_sample,
    outcome = outcome_formula,
    selection = selection_formula,
    svydesign = probability_sample,
    method_outcome = 'glm',
    family_outcome = 'gaussian',
    method_selection = "logit",
    control_selection = control_sel(est_method = "gee", gee_h_fun = 1)
  )
  
  dr_glm_mle_bias_min <- nonprob(
    data=as.data.frame(non_probability_sample),
    selection = selection_formula,
    outcome = outcome_formula,
    svydesign = probability_sample,
    method_outcome = 'glm',
    family_outcome = 'gaussian',
    method_selection = "logit",
    control_inference = control_inf(vars_combine = TRUE, vars_selection = TRUE, bias_correction = TRUE)
  )
  
  dr_glm_gee_bias_min <- nonprob(
    data=as.data.frame(non_probability_sample),
    outcome = outcome_formula,
    selection = selection_formula,
    svydesign = probability_sample,
    method_outcome = 'glm',
    family_outcome = 'gaussian',
    method_selection = "logit",
    control_selection = control_sel(est_method = "gee", gee_h_fun = 1),
    control_inference = control_inf(vars_combine = TRUE, vars_selection = TRUE, bias_correction = TRUE)
  )
  
  methods <- list(
    mi_glm = mi_glm,
    mi_glm_var_sel = mi_glm_var_sel,
    mi_npar = mi_npar,
    mi_nn_2 = mi_nn_2,
    mi_nn_5 = mi_nn_5,
    mi_nn_10 = mi_nn_10,
    mi_pmm = mi_pmm,
    mi_pmm_2 = mi_pmm_2,
    ipw = ipw,
    ipw_gee_1 = ipw_gee_1,
    ipw_gee_2 = ipw_gee_2,
    ipw_gee_var_sel = ipw_gee_var_sel,
    dr_glm_mle = dr_glm_mle,
    dr_glm_gee = dr_glm_gee,
    dr_glm_mle_bias_min = dr_glm_mle_bias_min,
    dr_glm_gee_bias_min = dr_glm_gee_bias_min
  )
  
  return(methods)
}

raw_results <- sapply(methods, function(method){
    c(
      mean = method$output$mean,
      se = method$output$SE,
      lower = method$confidence_interval$lower_bound,
      upper = method$confidence_interval$upper_bound
    )
  }
)
results <- as.data.frame(t(raw_results))


ggplot(data=results, aes(y=row.names(results)))+
  geom_errorbar(aes(xmin = lower, xmax = upper))+
  geom_point(aes(x=mean))+
  geom_vline(aes(xintercept = true_mean, color = "True Mean"), linetype = "dashed") +
  geom_vline(aes(xintercept = naive_mean, color = "Naive Mean"), linetype = "dotted") +
  scale_color_manual(
    name = "Reference Lines",
    values = c("True Mean" = "red", "Naive Mean" = "darkgrey")
  ) +
  #coord_cartesian(xlim = c(5.9, 8.1))+
  labs(
    title = "Comparison of Mean Estimators for Non-Probability Samples",
    subtitle = "Simulation study using the 'nonprobsvy' package",
    x = "Estimated Average Weekly Hours of Sport",
    y = "Estimation Method",
    caption = paste("True population mean (red) =", round(true_mean, 2))
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 10)
  )
