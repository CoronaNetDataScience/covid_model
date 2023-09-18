library(cmdstanr)
library(posterior)

real_data <- readRDS("../real_data.rds")
serology <- readRDS("../serology_real.rds")

init_vals <- function() {
  list(
    phi_raw = c(.001, .0001),
    world_infect = as.array(0.1),
    finding = 0.1,
    test_baseline = 0.1,
    mob_alpha_const = rep(0,real_data$G),
    suppress_effect_raw = rep(0, real_data$S),
    lockdown_effect_raw = rep(0, real_data$L),
    mob_effect_raw = rep(0, real_data$G),
    M_sigma = rep(1, real_data$G),
    country_test_raw = rep(0, real_data$num_country),
    country_test_raw2 = rep(0, real_data$num_country),
    country_test_raw3 = rep(0, real_data$num_country),
    mu_test_raw = as.array(0),
    mu_test_raw2 = as.array(0),
    mu_test_raw3 = rep(0, real_data$num_country),
    mu_poly = rep(0, 3),
    sigma_poly = rep(1, 3),
    poly1 = rep(0, real_data$num_country),
    poly2 = rep(0, real_data$num_country),
    poly3 = rep(0, real_data$num_country),
    sero_est = serology$inf_pr,
    pcr_spec = -6,
    sigma_test_raw = as.array(0.01),
    sigma_test_raw2 = as.array(0.01),
    sigma_test_raw3 = as.array(0.01),
    sigma_policy = c(1, 1, 1),
    alpha_infect = -6,
    alpha_test = as.array(-6),
    sigma_fear = 1,
    lockdown_med_raw = matrix(0, nrow = real_data$G,
                              ncol = real_data$L),
    suppress_med_raw = matrix(0, nrow = real_data$G,
                              ncol = real_data$S),
    lockdown_med_raw_fear = rep(0, real_data$L),
    suppress_med_raw_fear = rep(0, real_data$S - 1),
    M_Omega = diag(real_data$G),
    fear_const = 0
  )
}

trans <- TRUE

pan_model_scale <-
  cmdstan_model("../corona_tscs_betab_mix_prior.stan",
                cpp_options = list(stan_threads = TRUE))

us_fit_scale_mod <- pan_model_scale$sample(
  data = real_data,
  threads_per_chain = 4,
  init = init_vals,
  iter_warmup = 1000,
  iter_sampling = 5000,
  max_treedepth = 15,
  chains = 10,
  parallel_chains = 10
)

us_fit_scale_mod$save_object("../data/us_fit_scale_mod.rds")

summaries <- us_fit_scale_mod$summary()

summaries

us_fit_scale <-
  rstan::sflist2stanfit(lapply(us_fit_scale_mod$output_files(),
                               rstan::read_stan_csv))

saveRDS(us_fit_scale, "../data/us_fit_scale.rds")
saveRDS(summaries, "../data/rhat_summaries.rds")