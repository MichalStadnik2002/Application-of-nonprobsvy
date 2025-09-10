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
