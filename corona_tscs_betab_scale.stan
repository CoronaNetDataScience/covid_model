// Coronavirus tracking model 
// Robert Kubinec and Luiz Carvalho
// New York University Abu Dhabi & Getulio Vargas Foundation
// July 15, 2020
data {
    int time_all;
    int num_country;
    int num_rows;
    int cc[num_rows]; // country counter
    int cases[num_rows];
    int tests[num_rows];
    int S; // number of suppression measures
    int G; // google mobility data (by type of mobility)
    int L; // just lockdown data (for hierarchical predictor)
    int R; // number of seroprevalence essays
    int maxO; // maximum number of post-outbreak time points (could be greater than T due to left-censoring)
    matrix[maxO,3] ortho_time;
    matrix[num_rows,S] suppress; // time-varying suppression measures
    matrix[num_rows,G] mobility; // time-varying mobility measures
    matrix[num_rows,L] lockdown; // hierachical lockdown predictors
    matrix[num_rows,3] count_outbreak;
    vector[num_rows] month_cases;
    int sero_row[R]; // counters for which state/time points have CDC sero surveys
    matrix[R,5] sero; // sero-prevalence datas
    int country_pop[num_rows];
    real phi_scale; // prior on how much change there could be in infection rate over time
}
transformed data {

  matrix[num_rows, S] Q_supp;
  matrix[S, S] R_supp;
  matrix[S, S] R_supp_inverse;
  
  matrix[num_rows, G] Q_mob;
  matrix[G, G] R_mob;
  matrix[G, G] R_mob_inverse;
  
  matrix[num_rows, L] Q_lock;
  matrix[L, L] R_lock;
  matrix[L, L] R_lock_inverse;
  
  // thin and scale the QR decomposition
  Q_supp = qr_Q(suppress)[, 1:S] * sqrt(num_rows - 1);
  R_supp = qr_R(suppress)[1:S, ] / sqrt(num_rows - 1);
  R_supp_inverse = inverse(R_supp);
  
  Q_mob = qr_Q(mobility)[, 1:G] * sqrt(num_rows - 1);
  R_mob = qr_R(mobility)[1:G, ] / sqrt(num_rows - 1);
  R_mob_inverse = inverse(R_mob);
  
  Q_lock = qr_Q(lockdown)[, 1:L] * sqrt(num_rows - 1);
  R_lock = qr_R(lockdown)[1:L, ] / sqrt(num_rows - 1);
  R_lock_inverse = inverse(R_lock);

}
parameters {
  vector[num_country] poly1; // polinomial function of time
  vector[num_country] poly2; // polinomial function of time
  vector[num_country] poly3; // polinomial function of time
  real finding; // difficulty of identifying infected cases 
  //vector<lower=0,upper=1>[R] survey_prop; // variable that stores survey proportions from CDC data
  real world_infect; // infection rate based on number of travelers
  vector[S] suppress_effect_raw; // suppression effect of govt. measures, cannot increase virus transmission rate
  vector[L] lockdown_effect_raw;
  vector[L] suppress_hier_const[G];
  real test_lin_counter;
  real test_baseline;
  real test_lin_counter2;
  real<upper=-4.5> pcr_spec; // anticipated 1 - specificity of RT-PCR tests (taken from literature)
  // constraint equal to upper limit spec of 98 percent, 2% FPR of PCR test results
  vector[3] mu_poly; // hierarchical mean for poly coefficients
  vector[G] mob_effect_raw;
  real test_max_par;
  vector<lower=0>[3] sigma_poly; // varying sigma polys
  vector[G] mob_alpha_const; // mobility hierarchical intercepts
  vector<lower=0>[G] sigma_med;
  vector<lower=0>[num_country] country_test_raw; // unobserved rate at which countries are willing to test vs. number of infected
  vector[num_country] country_test_raw2;
  // we assume that as infection rates increase, more tests will be conducted
  vector[3] alpha; // other intercepts
  vector<lower=0>[2] phi; // shape parameter for infected
  real<lower=0> sigma_test_raw; // estimate of between-state testing heterogeneity
  real<lower=0> sigma_test_raw2;
}
transformed parameters {

  vector[num_rows] prop_infected; // modeled infection rates for domestic transmission
  vector[num_country] poly_nonc1; // non-centered poly parameters
  vector[num_country] poly_nonc2; // non-centered poly parameters
  vector[num_country] poly_nonc3; // non-centered poly parameters
  real<lower=0> sigma_test; 
  // vector[S] suppress_effect_raw;
  // 
  // // add in constrained world infect parameter
  // suppress_effect_raw = append_row(world_infect,suppress_effect_raw_small); 
  
  sigma_test = .1 * sigma_test_raw;
  
  // non-centering of polynomial time trends
  
  poly_nonc1 = mu_poly[1] + sigma_poly[1]*poly1;
  poly_nonc2 = mu_poly[2] + sigma_poly[2]*poly2;
  poly_nonc3 = mu_poly[3] + sigma_poly[3]*poly3;

  // latent infection rate (unobserved), on the logit scale (untransformed)
  prop_infected = alpha[2] + count_outbreak[,1] .* poly_nonc1[cc]  +
                  count_outbreak[,2] .* poly_nonc2[cc] +
                  count_outbreak[,3] .* poly_nonc3[cc] +
                  world_infect*month_cases +
                  Q_supp*suppress_effect_raw +
                  Q_lock*lockdown_effect_raw +
                  Q_mob*mob_effect_raw;

}
model {
  

  poly1 ~ normal(0,1);
  poly2 ~ normal(0,1);
  poly3 ~ normal(0,1);

  
  sigma_poly ~ exponential(.1);
  mu_poly ~ normal(0,10);
  world_infect ~ normal(0,3);
  lockdown_effect_raw ~ normal(0,5);
  alpha ~ normal(0,10); // this can reach extremely low values
  
  phi ~ exponential(phi_scale);
  mob_effect_raw ~ normal(0,5);
  suppress_effect_raw ~ normal(0,5);
  test_max_par ~ normal(0,5);
  test_baseline ~ normal(0,5);
  test_lin_counter ~ normal(0,5);
  test_lin_counter2 ~ normal(0,5);
  
  for(g in 1:G) {
    suppress_hier_const[g] ~ normal(0,5);
  }
  
   mob_alpha_const ~ normal(0,5);
    
  
  finding ~ normal(0,3);
  sigma_test_raw ~ exponential(.1);
  sigma_test_raw2 ~ exponential(.1);
  sigma_med ~ exponential(.1);
  country_test_raw ~ normal(0,sigma_test_raw); // more likely near the middle than the ends
  country_test_raw2 ~ normal(0,sigma_test_raw2);
  
  for(g in 1:G)
    to_vector(mobility[,g]) ~ normal(mob_alpha_const[g] + lockdown*suppress_hier_const[g],sigma_med[g]);
    
  //next model the true infection rate as a function of time since outbreak

    {
    // locations for cases and tests

    vector[num_rows] mix_prop_std = (prop_infected - mean(prop_infected)) ./ sd(prop_infected);
    vector[num_rows] mu_cases = inv_logit(pcr_spec + finding*mix_prop_std);
    vector[num_rows] mu_tests = inv_logit(alpha[1] + test_baseline * mix_prop_std +
    country_test_raw[cc] .* mix_prop_std .* count_outbreak[,1]  +
    test_lin_counter * count_outbreak[,1]);
    //country_test_raw2[cc] .* mix_prop_std .* count_outbreak[,2]  +
    //test_lin_counter2 * count_outbreak[,2]);
    
    
    tests ~ beta_binomial(country_pop,mu_tests*phi[1],(1-mu_tests) * phi[1]);
    cases ~ beta_binomial(tests,mu_cases *phi[2],(1-mu_cases) * phi[2]);
    
    // loop over serology surveys to add informative prior information
    for(r in 1:R) {
          int n = sero_row[r];
          
          real cum_infect = inv_logit(prop_infected[n]);
          
          // Beta prior = uncertainty proportional to sero essay sample size at time point t
          cum_infect ~ beta(.5 + sero[r,1].*sero[r,3],.5 + sero[r,3] - (sero[r,3].*sero[r,1]));
          
          // jacobian adjustment
          
          target += log(cum_infect) + log(1 + cum_infect);     
    }
  
    }
    
}

generated quantities {
  
  // convert QR estimates back to actual numbers
  
  vector[S] suppress_effect; // suppression effect of govt. measures, cannot increase virus transmission rate
  vector[L] lockdown_effect;
  vector[G] mob_effect;
  
  suppress_effect = R_supp_inverse * suppress_effect_raw;
  lockdown_effect = R_lock_inverse * lockdown_effect_raw;
  mob_effect = R_mob_inverse * mob_effect_raw;
  
}

