// Coronavirus tracking model 
// Robert Kubinec and Luiz Carvalho
// New York University Abu Dhabi & Getulio Vargas Foundation
// July 15, 2020
functions {
  
  // if(r_in(3,{1,2,3,4})) will evaluate as 1
  int r_in(int pos,int[] pos_var) {
   
    for (p in 1:(size(pos_var))) {
       if (pos_var[p]==pos) {
       // can return immediately, as soon as find a match
          return 1;
       } 
    }
    return 0;
  }
  
  real partial_sum(int[] y_slice,
                   int start, int end,
                   int[] cases,
                   int[] cases_sero,
                   int[] tests_sero,
                   vector phi,
                   int[] country_pop,
                   int[] country_pop_sero,
                   int num_country,
                   int num_rows,
                   vector mu_poly,
                   vector sigma_poly,
                   vector poly1,
                   vector poly2,
                   vector poly3,
                   vector alpha,
                   matrix count_outbreak,
                   matrix count_outbreak_sero,
                   real world_infect,
                   vector month_cases,
                   vector month_cases_sero,
                   vector suppress_effect_raw,
                   vector lockdown_effect_raw,
                   vector mob_effect_raw,
                   matrix Q_supp,
                   matrix Q_supp_sero,
                   matrix Q_lock,
                   matrix Q_lock_sero,
                   matrix Q_mob,
                   matrix Q_mob_sero,
                   int[] cc,
                   int[] cc_sero,
                   real pcr_spec,
                   real finding,
                   real test_baseline,
                   real test_lin_counter,
                   vector country_test_raw,
                   int S,
                   int G,
                   int L,
                   int R,
                   int[] sero_row,
                   matrix sero,
                   matrix mobility,
                   vector mob_alpha_const,
                   vector[] suppress_hier_const_raw,
                   vector sigma_med,
                   real sigma_test_raw,
                   vector country_test_raw2,
                   real test_lin_counter2,
                   int[,] sero_int,
                   vector sero_est) {
    
    vector[end - start + 1] prop_infected; // modeled infection rates for domestic transmission\
    int this_cc[end - start + 1] = cc[start:end];
    int this_cc2[end - start + 1] = cc[start:end];
    int min_cc = min(this_cc2);
    vector[max(this_cc) - min_cc + 1] poly_nonc1; // non-centered poly parameters
    vector[max(this_cc) - min_cc + 1] poly_nonc2; // non-centered poly parameters
    vector[max(this_cc) - min_cc + 1] poly_nonc3; // non-centered poly parameters
    real log_prob;
    real mu_infect;
    real sd_infect;
    
    vector[end - start + 1] mix_prop_std;
    vector[end - start + 1] mu_cases;
    vector[end - start + 1] mu_tests;
    
    poly_nonc1 = mu_poly[1] + sigma_poly[1]*poly1[min_cc:max(this_cc)];
    poly_nonc2 = mu_poly[2] + sigma_poly[2]*poly2[min_cc:max(this_cc)];
    poly_nonc3 = mu_poly[3] + sigma_poly[3]*poly3[min_cc:max(this_cc)];
    
    // log_prob = normal_lpdf(poly1[min_cc:max(this_cc)]|0,1);
    // log_prob += normal_lpdf(poly2[min_cc:max(this_cc)]|0,1);
    // log_prob += normal_lpdf(poly3[min_cc:max(this_cc)]|0,1);
    
    for(c in 1:num_elements(this_cc2))
        this_cc2[c] = this_cc2[c] - min_cc + 1;
        
    // latent infection rate (unobserved), on the logit scale (untransformed)

        prop_infected = alpha[2] + count_outbreak[start:end,1] .* poly_nonc1[this_cc2]  +
          count_outbreak[start:end,2] .* poly_nonc2[this_cc2] +
          count_outbreak[start:end,3] .* poly_nonc3[this_cc2] +
          world_infect*month_cases[start:end] +
          Q_supp[start:end,1:S]*suppress_effect_raw +
          Q_lock[start:end,1:L]*lockdown_effect_raw +
          Q_mob[start:end,1:G]*mob_effect_raw;
          
    mu_infect = mean(prop_infected);
    sd_infect = sd(prop_infected);
    
    mix_prop_std = (prop_infected - mu_infect) ./ sd_infect;               
    mu_cases = inv_logit(pcr_spec + finding*mix_prop_std);
    
    //log_prob = normal_lpdf(country_test_raw[min_cc:max(this_cc)]|0,sigma_test_raw); // more likely near the middle than the ends
    log_prob = 0;
    
    for(g in 1:G)
      log_prob += normal_lpdf(to_vector(mobility[start:end,g])|mob_alpha_const[g] + Q_lock[start:end,1:L]*suppress_hier_const_raw[g],sigma_med[g]);
    
    
    // mu_tests = inv_logit(alpha[1] + test_baseline * mix_prop_std +
    //                        country_test_raw[this_cc] .* mix_prop_std .* count_outbreak[start:end,1]  +
    //                        test_lin_counter * count_outbreak[start:end,1] +
    //                        test_lin_counter2 * count_outbreak[start:end,2] +
    //                        country_test_raw2[this_cc] .* mix_prop_std .* count_outbreak[start:end,2]);
    mu_tests = inv_logit(alpha[1] + country_test_raw[this_cc] .* mix_prop_std);
    log_prob += beta_binomial_lpmf(y_slice|country_pop[start:end],mu_tests*phi[1],(1-mu_tests) * phi[1]);
    log_prob += beta_binomial_lpmf(cases[start:end]|y_slice,mu_cases*phi[2],(1-mu_cases) * phi[2]);
    
    // loop over serology surveys to add informative prior information
    
      for(r in 1:R) {
        
          real sero_est_std;
          int n = sero_row[r];
          
      if(n>=start && n < end) {
        
        real mu_tests_sero;
        real mu_cases_sero;
        
        log_prob += beta_lpdf(sero_est[r]|.5 + sero[r,1]*sero[r,3],.5 + sero[r,3] - (sero[r,3]*sero[r,1]));  
        
        sero_est_std = (logit(sero_est[r]) - mu_infect)/sd_infect;
        
        mu_tests_sero = inv_logit(alpha[1] + country_test_raw[cc_sero[r]] * sero_est_std);
        mu_cases_sero = inv_logit(pcr_spec + finding*sero_est_std);
        
        log_prob += beta_binomial_lpmf(tests_sero[r]|country_pop_sero[r],mu_tests_sero*phi[1],(1-mu_tests_sero) * phi[1]);
        log_prob += beta_binomial_lpmf(cases_sero[r]|tests_sero[r],mu_cases_sero*phi[2],(1-mu_cases_sero) * phi[2]);
        
        // now do infection model
        
        log_prob += binomial_logit_lpmf(sero_int[r,1]|sero_int[r,2],alpha[2] + count_outbreak_sero[r,1] * poly_nonc1[cc_sero[r] - min_cc + 1]  +
          count_outbreak_sero[r,2] * poly_nonc2[cc_sero[r] - min_cc + 1] +
          count_outbreak_sero[r,3] * poly_nonc3[cc_sero[r] - min_cc + 1] +
          world_infect*month_cases_sero[r] +
          Q_supp_sero[r,1:S]*suppress_effect_raw +
          Q_lock_sero[r,1:L]*lockdown_effect_raw +
          Q_mob_sero[r,1:G]*mob_effect_raw);
      }
    }
          
          // Beta prior = uncertainty proportional to sero essay sample size at time point t
          
          
          // jacobian adjustment
          
          //log_prob += prop_infected[n-start+1] - 2 * log(1 + exp(prop_infected[n-start+1]));    
      
      
    
    return log_prob;
  }
}
data {
  int time_all;
  int num_country;
  int num_rows;
  int num_rows_orig;
  int cc[num_rows]; // country counter
  int cases[num_rows];
  int tests[num_rows];
  int S; // number of suppression measures
  int G; // google mobility data (by type of mobility)
  int L; // just lockdown data (for hierarchical predictor)
  int R; // number of seroprevalence essays
  matrix[num_rows_orig,S] suppress; // time-varying suppression measures
  matrix[num_rows_orig,G] mobility; // time-varying mobility measures
  matrix[num_rows_orig,L] lockdown; // hierachical lockdown predictors
  matrix[num_rows,3] count_outbreak;
  vector[num_rows] month_cases;
  vector[R] month_cases_sero;
  int sero_row[R]; // counters for which state/time points have CDC sero surveys
  matrix[R,5] sero; // sero-prevalence datas
  matrix[R,3] count_outbreak_sero;
  int sero_int[R,2]; // sero integers for binomial model
  int cases_sero[R];
  int tests_sero[R];
  int cc_sero[R];
  int country_pop[num_rows];
  int country_pop_sero[R];
  real phi_scale; // prior on how much change there could be in infection rate over time
}
transformed data {
  
  matrix[num_rows_orig, S] Q_supp_orig;
  matrix[num_rows,S] Q_supp;
  matrix[R,S] Q_supp_sero;
  matrix[S, S] R_supp;
  matrix[S, S] R_supp_inverse;
  
  matrix[num_rows_orig, G] Q_mob_orig;
  matrix[num_rows, G] Q_mob;
  matrix[R, G] Q_mob_sero;
  matrix[G, G] R_mob;
  matrix[G, G] R_mob_inverse;
  
  matrix[num_rows_orig, L] Q_lock_orig;
  matrix[num_rows, L] Q_lock;
  matrix[R, L] Q_lock_sero;
  matrix[L, L] R_lock;
  matrix[L, L] R_lock_inverse;
  
  int i = 0;
  
  // thin and scale the QR decomposition
  Q_supp_orig = qr_Q(suppress)[, 1:S] * sqrt(num_rows - 1);
  R_supp = qr_R(suppress)[1:S, ] / sqrt(num_rows - 1);
  R_supp_inverse = inverse(R_supp);
  
  Q_mob_orig = qr_Q(mobility)[, 1:G] * sqrt(num_rows - 1);
  R_mob = qr_R(mobility)[1:G, ] / sqrt(num_rows - 1);
  R_mob_inverse = inverse(R_mob);
  
  Q_lock_orig = qr_Q(lockdown)[, 1:L] * sqrt(num_rows - 1);
  R_lock = qr_R(lockdown)[1:L, ] / sqrt(num_rows - 1);
  R_lock_inverse = inverse(R_lock);
  
  
  // split data into sero/non-sero data
  for(n in 1:num_rows_orig) {
    if(r_in(n,sero_row)) {
      
      i += 1;
      Q_lock_sero[i,1:L] = Q_lock_orig[n,1:L];
      Q_supp_sero[i,1:S] = Q_supp_orig[n,1:S];
      Q_mob_sero[i,1:G] = Q_lock_orig[n,1:G];
      
    } else {
      Q_lock[n-i,1:L] = Q_lock_orig[n,1:L];
      Q_supp[n-i,1:S] = Q_supp_orig[n,1:S];
      Q_mob[n-i,1:G] = Q_lock_orig[n,1:G];
    }
  }
  
}
parameters {
  vector[num_country] poly1; // polinomial function of time
  vector[num_country] poly2; // polinomial function of time
  vector[num_country] poly3; // polinomial function of time
  vector<lower=0,upper=1>[R] sero_est;
  real finding; // difficulty of identifying infected cases 
  //vector<lower=0,upper=1>[R] survey_prop; // variable that stores survey proportions from CDC data
  real<lower=0> world_infect; // infection rate based on number of travelers
  vector[S] suppress_effect_raw; // suppression effect of govt. measures, cannot increase virus transmission rate
  vector[L] lockdown_effect_raw;
  vector[L] suppress_hier_const_raw[G];
  real<lower=0> test_lin_counter;
  real<lower=0> test_baseline;
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
  vector<upper=0>[num_country] country_test_raw2;
  // we assume that as infection rates increase, more tests will be conducted
  vector[2] alpha; // other intercepts
  vector<lower=0>[2] phi; // shape parameter for infected
  real<lower=0> sigma_test_raw; // estimate of between-state testing heterogeneity
  real<lower=0> sigma_test_raw2;
}
model {
  
  int grainsize = 600;
  
  poly1 ~ normal(0,1);
  poly2 ~ normal(0,1);
  poly3 ~ normal(0,1);
  
  country_test_raw ~ normal(0,sigma_test_raw);
  
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
  
  mob_alpha_const ~ normal(0,5);
  
  
  finding ~ normal(0,3);
  sigma_test_raw ~ exponential(.1);
  sigma_test_raw2 ~ exponential(.1);
  sigma_med ~ exponential(.1);
  country_test_raw2 ~ normal(0,sigma_test_raw2);

target += reduce_sum_static(partial_sum, tests,
                     grainsize,
                     cases, cases_sero,
                     tests_sero,
                     phi,country_pop,
                     country_pop_sero,
                     num_country,
                     num_rows,
                     mu_poly,
                     sigma_poly,
                     poly1,
                     poly2,
                     poly3,
                     alpha,
                     count_outbreak,
                     count_outbreak_sero,
                     world_infect,
                     month_cases,
                     month_cases_sero,
                     suppress_effect_raw,
                     lockdown_effect_raw,
                     mob_effect_raw,
                     Q_supp,
                     Q_supp_sero,
                     Q_lock,
                     Q_lock_sero,
                     Q_mob,
                     Q_mob_sero,
                     cc,
                     cc_sero,
                     pcr_spec,
                     finding,
                     test_baseline,
                     test_lin_counter,
                     country_test_raw,
                     S,
                     G,
                     L,
                     R,
                     sero_row,
                     sero,
                     mobility,
                     mob_alpha_const,
                     suppress_hier_const_raw,
                     sigma_med,
                     sigma_test_raw,
                     country_test_raw2,
                     test_lin_counter2,
                     sero_int,
                     sero_est);

}

generated quantities {
  
  // convert QR estimates back to actual numbers
  
  vector[S] suppress_effect; // suppression effect of govt. measures, cannot increase virus transmission rate
  vector[L] lockdown_effect;
  vector[G] mob_effect;
  vector[L] suppress_hier_const[G];
  
  suppress_effect = R_supp_inverse * suppress_effect_raw;
  lockdown_effect = R_lock_inverse * lockdown_effect_raw;
  mob_effect = R_mob_inverse * mob_effect_raw;
  
  for(g in 1:G)
    suppress_hier_const[g] = R_lock_inverse * suppress_hier_const_raw[g];
  
}

