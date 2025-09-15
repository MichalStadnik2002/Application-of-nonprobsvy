library(tidyverse)

N <- 10000
beta_age  <- 0.15
beta_age2 <- -0.002
beta_sex <- -0.5
beta_bmi <- -0.05
beta_sex_bmi2 <- 0.01

population <- tibble(
  index = 1:N,
  sex = rbinom(N, 1, 0.51),
  bmi = rgamma(N, 4.5, scale=2.2)+14,
  age = rbeta(N, 2, 3)*(80-18)+18,
  base_time_sport = 5 +
    beta_age * age +
    beta_age2 * age^2 +
    beta_sex * sex +
    beta_bmi * (bmi - 22)^2 +
    beta_sex_bmi2 * sex * (bmi - 23)^2,
  time_sport = pmax(0, base_time_sport + rnorm(N, mean = 0, sd = 3))
)
true_mean <- mean(population$time_sport)
true_sd <- sd(population$time_sport)

gamma_0 <- -2.5
gamma_1 <- -0.05
gamma_2 <- -0.1

population <- population |> 
  mutate(
    logit_pi = gamma_0 + gamma_1 * age + gamma_2 * (bmi-22)^2,
    pi_np = 1 / (1 + exp(-logit_pi))
  )

n_np <- 1000
non_probability_sample <- population  |> 
  slice_sample(n = n_np, weight_by = pi_np) |> 
  select(!c(base_time_sport, logit_pi, pi_np))

n_p <- 2000
probability_sample <- population |>
  anti_join(non_probability_sample, by='index') |> 
  slice_sample(n = n_p) |> 
  select(!c(base_time_sport, time_sport, logit_pi, pi_np)) |> 
  mutate(
    weight = N/n_p,
    population_size = N
  )

naive_mean <- mean(non_probability_sample$time_sport)
naive_sd <- sd(non_probability_sample$time_sport)

library(nonprobsvy)
library(survey)

probability_sample_svy <- svydesign(
  ids = ~1,
  data = probability_sample,
  weights = ~weight,
  fpc = ~population_size
)

mi_glm <- nonprob(
  data=non_probability_sample, 
  outcome = time_sport ~ sex + bmi + age,
  svydesign = probability_sample_svy,
  method_outcome = 'glm',
  family_outcome = 'gaussian'
  )

mi_glm_var_sel <- nonprob(
  data=non_probability_sample, 
  outcome = time_sport ~ sex + bmi + age,
  svydesign = probability_sample_svy,
  method_outcome = 'glm',
  family_outcome = 'gaussian',
  control_outcome = control_out(penalty = "lasso"), 
  control_inference = control_inf(vars_selection = TRUE)
)

mi_npar <- nonprob(
  data=non_probability_sample, 
  outcome = time_sport ~ sex + bmi + age,
  svydesign = probability_sample_svy,
  method_outcome = 'npar',
)

mi_nn_2 <- nonprob(
  data=non_probability_sample, 
  outcome = time_sport ~ sex + bmi + age,
  svydesign = probability_sample_svy,
  method_outcome = 'nn',
  control_outcome = control_out(k=2)
)

mi_nn_5 <- nonprob(
  data=non_probability_sample, 
  outcome = time_sport ~ sex + bmi + age,
  svydesign = probability_sample_svy,
  method_outcome = 'nn',
  control_outcome = control_out(k=5)
)

mi_nn_10 <- nonprob(
  data=non_probability_sample, 
  outcome = time_sport ~ sex + bmi + age,
  svydesign = probability_sample_svy,
  method_outcome = 'nn',
  control_outcome = control_out(k=10)
)

mi_pmm <- nonprob(
  data=non_probability_sample, 
  outcome = time_sport ~ sex + bmi + age,
  svydesign = probability_sample_svy,
  method_outcome = 'pmm',
  family_outcome = 'gaussian'
)

mi_pmm_2 <- nonprob(
  data=non_probability_sample, 
  outcome = time_sport ~ sex + bmi + age,
  svydesign = probability_sample_svy,
  method_outcome = 'pmm',
  family_outcome = 'gaussian',
  control_outcome = control_out(pmm_match_type = 2)
)

ipw <- nonprob(
  data=non_probability_sample, 
  selection = ~ sex + bmi + age,
  target = ~ time_sport,
  svydesign = probability_sample_svy,
  method_selection = "logit"
)

ipw_gee_1 <- nonprob(
  data=non_probability_sample, 
  selection = ~ sex + bmi + age,
  target = ~ time_sport,
  svydesign = probability_sample_svy,
  method_selection = "logit",
  control_selection = control_sel(est_method = "gee", gee_h_fun = 1)
)

ipw_gee_2 <- nonprob(
  data=non_probability_sample, 
  selection = ~ sex + bmi + age,
  target = ~ time_sport,
  svydesign = probability_sample_svy,
  method_selection = "logit",
  control_selection = control_sel(est_method = "gee", gee_h_fun = 2)
)

ipw_gee_var_sel <- nonprob(
  data=non_probability_sample, 
  selection = ~ sex + bmi + age,
  target = ~ time_sport,
  svydesign = probability_sample_svy,
  method_selection = "logit",
  control_selection = control_sel(est_method = "gee", gee_h_fun = 1),
  control_inference = control_inf(vars_selection = TRUE)
)
