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
          return p;
       } 
    }
    return 0;
    }
  
  real partial_sum(int[,] y_slice,
                   int start, int end,
                   int[] tests,
                   int[] cases,
                   vector phi,
                   int[] country_pop,
                   int num_country,
                   int num_rows,
                   vector mu_poly,
                   vector sigma_poly,
                   vector poly1,
                   vector poly2,
                   vector poly3,
                   real alpha_test,
                   real alpha_infect,
                   matrix count_outbreak,
                   real world_infect,
                   vector month_cases,
                   vector suppress_effect_raw,
                   vector lockdown_effect_raw,
                   vector mob_effect_raw,
                   matrix Q_supp,
                   matrix Q_supp2,
                   matrix Q_lock,
                   matrix Q_mob,
                   int[] cc,
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
                   vector[] mob_array,
                   vector mob_alpha_const,
                   vector[] lockdown_med_raw,
                   vector[] suppress_med_raw,
                   vector suppress_med_raw_fear,
                   vector lockdown_med_raw_fear,
                   real sigma_test_raw,
                   vector country_test_raw2,
                   real test_lin_counter2,
                   vector sero_est,
                   real test_max_par,
                   vector test_max,
                   matrix lin_counter,
                   real mu_test_raw,
                   real mu_test_raw2,
                   real sigma_test_raw2,
                   matrix M_Sigma,
                   vector fear,
                   real fear_const,
                   real sigma_fear) {
                     
    // big loop over states
    real log_prob = 0;
    for(r in 1:size(y_slice)) {
        
        int s = y_slice[r,1];
        int start2 = y_slice[r,2];
        int end2 = y_slice[r,3];
      
        vector[end2 - start2 + 1] prop_infected; // modeled infection rates for domestic transmission\
        
        int obs = end2 - start2 + 1;
        real poly_nonc1; // non-centered poly parameters
        real poly_nonc2; // non-centered poly parameters
        real poly_nonc3; // non-centered poly parameters
        real country1;
        real country2;
        real mu_infect;
        real sd_infect;
        
        vector[end2 - start2 + 1] mix_prop_std;
        vector[end2 - start2 + 1] mu_cases;
        vector[end2 - start2 + 1] mu_tests;
        vector[G] mu_mob[end2 - start2 + 1];
        
        poly_nonc1 = mu_poly[1] + sigma_poly[1]*poly1[s];
        poly_nonc2 = mu_poly[2] + sigma_poly[2]*poly2[s];
        poly_nonc3 = mu_poly[3] + sigma_poly[3]*poly3[s];
        
        country1 = mu_test_raw + sigma_test_raw*country_test_raw[s];
        country2 = mu_test_raw2 + sigma_test_raw2*country_test_raw2[s];
            
        // latent infection rate (unobserved), on the logit scale (untransformed)
        // constrained to *always* increase
        
        for(i in 1:obs) {
              if(i==1) {
                prop_infected[1] = alpha_infect + count_outbreak[start2,1] * poly_nonc1 +
                    count_outbreak[start2,2] * poly_nonc2 +
                    count_outbreak[start2,3] * poly_nonc3 +
                    world_infect*month_cases[start2] +
                    Q_supp[start2,1:S]*suppress_effect_raw +
                    Q_lock[start2,1:L]*lockdown_effect_raw +
                    Q_mob[start2,1:G]*mob_effect_raw;
              } else {
                prop_infected[i] = exp(alpha_infect + count_outbreak[start2+i-1,1] * poly_nonc1  +
                    count_outbreak[start2+i-1,2] *poly_nonc2 +
                    count_outbreak[start2+i-1,3] * poly_nonc3 +
                    world_infect*month_cases[start2+i-1] +
                    Q_supp[start2+i-1,1:S]*suppress_effect_raw +
                    Q_lock[start2+i-1,1:L]*lockdown_effect_raw +
                    Q_mob[start2+i-1,1:G]*mob_effect_raw) + prop_infected[i - 1];
              }
        }
              
            //need a recursive function to transform to ordered vector
            

              
        mu_infect = mean(prop_infected);
        sd_infect = sd(prop_infected);
        
        mix_prop_std = inv_logit(prop_infected);               
        mu_cases = inv_logit(pcr_spec + finding*mix_prop_std);

        log_prob += normal_lpdf(country_test_raw[s]|0,1); // more likely near the middle than the ends
        log_prob += normal_lpdf(country_test_raw2[s]|0,1); // more likely near the middle than the ends
        
        log_prob += normal_lpdf(poly1[s]|0,1);
        log_prob += normal_lpdf(poly2[s]|0,1);
        log_prob += normal_lpdf(poly3[s]|0,1);
        
        //mobility mediation

        for(g in 1:G) {
          mu_mob[1:obs,g] = to_array_1d(mob_alpha_const[g] +
          Q_lock[start2:end2,1:L]*lockdown_med_raw[g]  +
          Q_supp[start2:end2,1:S]*suppress_med_raw[g]); 
        }
        
          log_prob += multi_normal_cholesky_lpdf(mob_array[start2:end2,1:G]|mu_mob,M_Sigma);
          
        // fear mediation
          log_prob += normal_lpdf(fear[start2:end2]|fear_const + Q_lock[start2:end2,1:L]*lockdown_med_raw_fear +
                                Q_supp2[start2:end2,1:(S-1)]*suppress_med_raw_fear,sigma_fear);

        
    mu_tests = inv_logit(alpha_test + 
                          test_baseline *  mix_prop_std +
                          country1 * mix_prop_std .* lin_counter[start2:end2,1] +
                          test_lin_counter * lin_counter[start2:end2,1]);
        
        // observed data model
        log_prob += beta_binomial_lpmf(tests[start2:end2]|country_pop[start2:end2],mu_tests*phi[1],(1-mu_tests) * phi[1]);
       log_prob += beta_binomial_lpmf(cases[start2:end2]|country_pop[start2:end2],mu_cases*phi[2],(1-mu_cases) * phi[2]);
    
        // loop over serology surveys to add informative prior information
        
        for(n in (start2+1):end2) {
          int q = r_in(n,sero_row);
          real cur_value = inv_logit(prop_infected[n-start2+1]);
          if(q>0) {
            //log_prob += beta_lpdf(prop_infected[n-start2+1]|logit(sero[q,1]),.01);
            log_prob += beta_lpdf(cur_value|.5 + sero[q,1]*sero[q,3],
                                            5 + (sero[q,3] - sero[q,1]*sero[q,3]));
            log_prob += log(prop_infected[n-start2+1] - prop_infected[n-start2]) +
                        log_inv_logit(prop_infected[n-start2+1]) + log1m_inv_logit(prop_infected[n-start2+1]);
          } else {
            // prior implies inflation rate (cases/true infected) is between 2 - 20
            if((n-start2)>60) {
              //log_prob += beta_lpdf(cur_value|.5, 5);
              //log_prob += log(prop_infected[n-start2+1] - prop_infected[n-start2]) +
                //        log_inv_logit(prop_infected[n-start2+1]) + log1m_inv_logit(prop_infected[n-start2+1]);
                // log_prob += lognormal_lpdf(inv_logit(prop_infected[n-start2+1])/((cases[n-start2+1]*1.0)/(country_pop[n-start2+1]*1.0))|2.1,.4);
                // log_prob +=  log(prop_infected[n-start2+1] - prop_infected[n-start2]) - 2 * log(1 + prop_infected[n-start2+1] - prop_infected[n-start2]);
            }
    
          }
        }
        
      }
      return log_prob;
    }
    
    
}
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
  matrix[num_rows,S] suppress; // time-varying suppression measures
  matrix[num_rows,S-1] suppress2; // without COVID poll
  matrix[num_rows,G] mobility; // time-varying mobility measures
  matrix[num_rows,L] lockdown; // hierachical lockdown predictors
  vector[num_rows] fear; // COVID poll
  matrix[num_rows,3] count_outbreak;
  vector[num_rows] month_cases;
  vector[num_rows] test_max;
  matrix[R,5] sero; // sero-prevalence datas
  int sero_row[R];
  int country_pop[num_rows];
  matrix[num_rows,2] lin_counter;
  real phi_scale; // prior on how much change there could be in infection rate over time
  int states[num_country,3];
}
transformed data {
  
  matrix[num_rows,S] Q_supp;
  matrix[S, S] R_supp;
  matrix[S, S] R_supp_inverse;
  
  matrix[num_rows,S-1] Q_supp2;
  matrix[S-1, S-1] R_supp2;
  matrix[S-1, S-1] R_supp_inverse2;
  
  matrix[num_rows, G] Q_mob;
  matrix[G, G] R_mob;
  matrix[G, G] R_mob_inverse;
  
  matrix[num_rows, L] Q_lock;
  matrix[L, L] R_lock;
  matrix[L, L] R_lock_inverse;
  vector[G] mob_array[num_rows];
  
  
  // thin and scale the QR decomposition
  Q_supp = qr_Q(suppress)[, 1:S] * sqrt(num_rows - 1);
  Q_supp2 = qr_Q(suppress2)[, 1:(S-1)] * sqrt(num_rows - 1);
  R_supp = qr_R(suppress)[1:S, ] / sqrt(num_rows - 1);
  R_supp2 = qr_R(suppress2)[1:(S-1), ] / sqrt(num_rows - 1);
  R_supp_inverse = inverse(R_supp);
  R_supp_inverse2 = inverse(R_supp2);
  
  Q_mob = qr_Q(mobility)[, 1:G] * sqrt(num_rows - 1);
  R_mob = qr_R(mobility)[1:G, ] / sqrt(num_rows - 1);
  R_mob_inverse = inverse(R_mob);
  
  Q_lock = qr_Q(lockdown)[, 1:L] * sqrt(num_rows - 1);
  R_lock = qr_R(lockdown)[1:L, ] / sqrt(num_rows - 1);
  R_lock_inverse = inverse(R_lock);
  
  for(g in 1:G) {
    for(n in 1:num_rows) {
      mob_array[n,g] = mobility[n,g];
    }
  }
  
  
}
parameters {
  vector[num_country] poly1; // polinomial function of time
  vector[num_country] poly2; // polinomial function of time
  vector[num_country] poly3; // polinomial function of time
  real mu_test_raw;
  vector<lower=0,upper=1>[R] sero_est;
  real<lower=0> finding; // difficulty of identifying infected cases 
  //vector<lower=0,upper=1>[R] survey_prop; // variable that stores survey proportions from CDC data
  real<lower=0> world_infect; // infection rate based on number of travelers
  vector[S] suppress_effect_raw; // suppression effect of govt. measures, cannot increase virus transmission rate
  vector[L] lockdown_effect_raw;
  vector[L] lockdown_med_raw[G];
  vector[S] suppress_med_raw[G];
  vector[L] lockdown_med_raw_fear;
  vector[S-1] suppress_med_raw_fear;
  real<lower=0> test_lin_counter;
  real<lower=0> test_baseline;
  real test_lin_counter2;
  real mu_test_raw2;
  real pcr_spec; // anticipated 1 - specificity of RT-PCR tests (taken from literature)
  // constraint equal to upper limit spec of 98 percent, 2% FPR of PCR test results
  vector[3] mu_poly; // hierarchical mean for poly coefficients
  vector[G] mob_effect_raw;
  real test_max_par;
  vector<lower=0>[3] sigma_poly; // varying sigma polys
  vector[G] mob_alpha_const; // mobility hierarchical intercepts
  vector[num_country] country_test_raw; // unobserved rate at which countries are willing to test vs. number of infected
  vector[num_country] country_test_raw2;
  real alpha_infect; // other intercepts
  real alpha_test;
  vector<lower=0>[2] phi; // shape parameter for infected
  real<lower=0> sigma_test_raw; // estimate of between-state testing heterogeneity
  real<lower=0> sigma_test_raw2;
  cholesky_factor_corr[G] M_Omega; // these are for the MVN for mobility data
  vector<lower=0>[G] M_sigma;
  real fear_const;
  real sigma_fear;
}
model {
  matrix[G, G] M_Sigma;
  int grainsize = 1;
  
  sigma_poly ~ normal(0,5);
  mu_poly ~ normal(0,10);
  mu_test_raw ~ normal(0,200);
  mu_test_raw2 ~ normal(-5,2);
  world_infect ~ normal(0,3);
  lockdown_effect_raw ~ normal(0,5);
  alpha_infect ~ normal(0,10); // this can reach extremely low values
  alpha_test ~ normal(0,10);

  phi ~ exponential(phi_scale);
  mob_effect_raw ~ normal(0,5);
  suppress_effect_raw ~ normal(0,5);
  test_max_par ~ normal(0,5);
  test_baseline ~ normal(0,100);
  test_lin_counter ~ normal(0,20);
  test_lin_counter2 ~ normal(0,5);
  
  mob_alpha_const ~ normal(0,5);
  pcr_spec ~ normal(0,10);
  
  finding ~ normal(100,50);
  sigma_test_raw ~ exponential(.1);
  sigma_test_raw2 ~ normal(0,.1);
  sigma_fear ~ exponential(.1);
  fear_const ~ normal(0,5);
  sero_est ~ beta(.5 + sero[,1] .* sero[,3],.5 + sero[,3] - (sero[,3] .* sero[,1]));
  
  for(g in 1:G) {
    suppress_med_raw[g] ~ normal(0,5);
    lockdown_med_raw[g] ~ normal(0,5);
  }
  
  suppress_med_raw_fear ~ normal(0,5);
  lockdown_med_raw_fear ~ normal(0,5);
  
  M_Omega ~ lkj_corr_cholesky(4);
  M_sigma ~ cauchy(0, 2.5);
  
  M_Sigma = diag_pre_multiply(M_sigma, M_Omega); // Cholesky decomp for MVN for mobility

target += reduce_sum_static(partial_sum, states,
                     grainsize,
                    tests,
                     cases,
                     phi,country_pop,
                     num_country,
                     num_rows,
                     mu_poly,
                     sigma_poly,
                     poly1,
                     poly2,
                     poly3,
                     alpha_test,
                     alpha_infect,
                     count_outbreak,
                     world_infect,
                     month_cases,
                     suppress_effect_raw,
                     lockdown_effect_raw,
                     mob_effect_raw,
                     Q_supp,
                     Q_supp2,
                     Q_lock,
                     Q_mob,
                     cc,
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
                     mob_array,
                     mob_alpha_const,
                     lockdown_med_raw,
                     suppress_med_raw,
                     suppress_med_raw_fear,
                     lockdown_med_raw_fear,
                     sigma_test_raw,
                     country_test_raw2,
                     test_lin_counter2,
                     sero_est,
                     test_max_par,
                     test_max,
                     lin_counter,
                     mu_test_raw,
                     mu_test_raw2,
                     sigma_test_raw2,
                     M_Sigma,
                     fear,
                     fear_const,
                     sigma_fear);

}
generated quantities {
  
  // convert QR estimates back to actual numbers
  
  vector[S] suppress_effect; // suppression effect of govt. measures, cannot increase virus transmission rate
  vector[L] lockdown_effect;
  vector[G] mob_effect;
  vector[L] lockdown_med[G];
  vector[S] suppress_med[G];
  vector[L] lockdown_med_fear;
  vector[S-1] suppress_med_fear;
  vector[num_rows] prop_infect_out;
  
  suppress_effect = R_supp_inverse * suppress_effect_raw;
  lockdown_effect = R_lock_inverse * lockdown_effect_raw;
  mob_effect = R_mob_inverse * mob_effect_raw;
  
  for(g in 1:G) {
    lockdown_med[g] = R_lock_inverse * lockdown_med_raw[g];
    suppress_med[g] = R_supp_inverse * suppress_med_raw[g];
  }
  
  suppress_med_fear = R_supp_inverse2 * suppress_med_raw_fear;
  lockdown_med_fear = R_lock_inverse * lockdown_med_raw_fear;
    
  
  for(s in 1:num_country) {
    
        real poly_nonc1; // non-centered poly parameters
        real poly_nonc2; // non-centered poly parameters
        real poly_nonc3; // non-centered poly parameters
        
        int start2 = states[s,2];
        int end2 = states[s,3];
        int obs = end2 - start2 + 1;
        
        poly_nonc1 = mu_poly[1] + sigma_poly[1]*poly1[s];
        poly_nonc2 = mu_poly[2] + sigma_poly[2]*poly2[s];
        poly_nonc3 = mu_poly[3] + sigma_poly[3]*poly3[s];
    
    for(i in 1:obs) {
              if(i==1) {
                prop_infect_out[start2] = alpha_infect + count_outbreak[start2,1] * poly_nonc1 +
                    count_outbreak[start2,2] * poly_nonc2  +
                    count_outbreak[start2,3] * poly_nonc3  +
                    world_infect*month_cases[start2] +
                    Q_supp[start2,1:S]*suppress_effect_raw +
                    Q_lock[start2,1:L]*lockdown_effect_raw +
                    Q_mob[start2,1:G]*mob_effect_raw;
              } else {
                prop_infect_out[start2 + i - 1] = exp(alpha_infect + count_outbreak[start2+i-1,1] * poly_nonc1  +
                    count_outbreak[start2+i-1,2] * poly_nonc2 +
                    count_outbreak[start2+i-1,3] * poly_nonc3 +
                    world_infect*month_cases[start2+i-1] +
                    Q_supp[start2+i-1,1:S]*suppress_effect_raw +
                    Q_lock[start2+i-1,1:L]*lockdown_effect_raw +
                    Q_mob[start2+i-1,1:G]*mob_effect_raw) + prop_infect_out[start2 + i - 2];
              }
        }
  }
          
}


