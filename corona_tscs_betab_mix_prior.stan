// Coronavirus tracking model 
// Robert Kubinec and Luiz Carvalho
// New York University Abu Dhabi & Getulio Vargas Foundation
// January 5, 2021
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
                   real[] alpha_test,
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
                   //real test_lin_counter,
                   vector country_test_raw,
                   int S,
                   int G,
                   int L,
                   int R,
                   int RE,
                   int[] sero_row,
                   int[] sero_row_real,
                   int[,] sero,
                   matrix sero_real,
                   matrix mobility,
                   vector[] mob_array,
                   vector mob_alpha_const,
                   vector[] lockdown_med_raw,
                   vector[] suppress_med_raw,
                   vector suppress_med_raw_fear,
                   vector lockdown_med_raw_fear,
                   real[] sigma_test_raw,
                   vector country_test_raw2,
                   vector country_test_raw3,
                   //real test_lin_counter2,
                   //real test_max_par,
                   //vector test_max,
                   matrix lin_counter,
                   real[] mu_test_raw,
                   real[] mu_test_raw2,
                   real[] mu_test_raw3,
                   real[] sigma_test_raw2,
                   real[] sigma_test_raw3,
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
        real country1s;
        real country2s;
        real country3s;
        //real country1f;
        //real country2f;
        //real country3f;
        real mu_infect;
        real sd_infect;
        
        vector[end2 - start2 + 1] prop_success;
        vector[end2 - start2 + 1] prop_fail;
        vector[end2 - start2 + 1] mu_cases;
        vector[end2 - start2 + 1] mu_tests;
        vector[G] mu_mob[end2 - start2 + 1];
        
        poly_nonc1 = mu_poly[1] + sigma_poly[1]*poly1[s];
        poly_nonc2 = mu_poly[2] + sigma_poly[2]*poly2[s];
        poly_nonc3 = mu_poly[3] + sigma_poly[3]*poly3[s];
        
        country1s = mu_test_raw[1] + sigma_test_raw[1]*country_test_raw[s];
        country2s = mu_test_raw2[1] + sigma_test_raw2[1]*country_test_raw2[s];
        country3s = mu_test_raw3[1] + sigma_test_raw3[1]*country_test_raw3[s];
        
        //country1f = mu_test_raw[2] + sigma_test_raw[2]*country_test_raw[2,s];
        //country2f = mu_test_raw2[2] + sigma_test_raw2[2]*country_test_raw2[2,s];
        //country3f = mu_test_raw3[2] + sigma_test_raw3[2]*country_test_raw3[2,s];
        
        //country1 = country_test_raw[1];
        //country2 = country_test_raw2[1];
            
        // latent infection rate (unobserved), on the logit scale (untransformed)
        // constrained to *always* increase
        
        for(i in 1:obs) {
              if(i==1) {
                prop_infected[1] = alpha_infect + 
                count_outbreak[start2,1] * poly_nonc1 +
                  count_outbreak[start2,2] * poly_nonc2 +
                    count_outbreak[start2,3] * poly_nonc3 +
                    world_infect*month_cases[start2] +
                    Q_supp[start2,1:S]*suppress_effect_raw +
                    Q_lock[start2,1:L]*lockdown_effect_raw +
                    Q_mob[start2,1:G]*mob_effect_raw;
              } else {
                prop_infected[i] = exp(alpha_infect + 
                count_outbreak[start2+i-1,1] * poly_nonc1  +
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
        
        prop_success = prop_infected; 
        prop_fail = 1 - inv_logit(prop_infected);
        
        mu_cases = inv_logit(pcr_spec + finding*prop_success);
        //mu_cases[2] = exp(pcr_spec[2] + finding[2]*prop_fail);

        // log_prob += normal_lpdf(to_vector(country_test_raw[1:2,s])|0,1); // more likely near the middle than the ends
        // log_prob += normal_lpdf(to_vector(country_test_raw2[1:2,s])|0,1); // more likely near the middle than the ends
        // log_prob += normal_lpdf(to_vector(country_test_raw3[1:2,s])|0,1);
        
        log_prob += normal_lpdf(country_test_raw[s]|0,1); // more likely near the middle than the ends
        log_prob += normal_lpdf(country_test_raw2[s]|0,1); // more likely near the middle than the ends
        log_prob += normal_lpdf(country_test_raw3[s]|0,1);
        
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

        
          mu_tests = inv_logit(alpha_test[1] + 
                          country1s * lin_counter[start2:end2,1] +
                          country2s * lin_counter[start2:end2,2] +
                          //country3s * lin_counter[start2:end2,3] +
                          test_baseline * prop_success);
        
        // observed data model
        // loop over serology surveys to add informative prior information
        
        for(n in start2:end2) {
          
          int q = r_in(n,sero_row);
          int p = r_in(n, sero_row_real);
          real cur_infect = cases[n]*1.0 / country_pop[n]*1.0;
          
          if(q <= 0 && p <= 0) {
            
            log_prob += beta_binomial_lpmf(cases[n]|country_pop[n],mu_cases[n-start2+1]*phi[1],(1-mu_cases[n-start2+1])*phi[1]);
            log_prob += beta_binomial_lpmf(tests[n]|country_pop[n],mu_tests[n-start2+1]*phi[2],(1-mu_tests[n-start2+1])*phi[2]);

          } else if(p>0) {
            
            // expert survey data
            // beta prior
            // include Jacobian transformation
            
            log_prob += beta_lpdf(inv_logit(prop_infected[n-start2+1])|sero_real[p,1],sero_real[p,2]);
            log_prob += log_inv_logit(prop_infected[n-start2+1]) + log1m_inv_logit(prop_infected[n-start2+1]);
            
            } else if(q > 0) {
            
            // scaling function. we use seroprevalance data to 
            // set a ground truth for the relationship between covariates and
            // infections measured non-parametrically
            
            log_prob += binomial_lpmf(sero[q,1]|sero[q,2],inv_logit(prop_infected[n-start2+1]));
            
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
  int RE; // number of expert surveys
  matrix[num_rows,S] suppress; // time-varying suppression measures
  matrix[num_rows,S-1] suppress2; // without COVID poll
  matrix[num_rows,G] mobility; // time-varying mobility measures
  matrix[num_rows,L] lockdown; // hierachical lockdown predictors
  vector[num_rows] fear; // COVID poll
  matrix[num_rows,3] count_outbreak;
  vector[num_rows] month_cases;
  vector[num_rows] test_max;
  int sero[R,2]; // sero-prevalence datas
  int sero_row[R];
  matrix[RE,2] sero_real; // expert survey datas
  int sero_row_real[RE];
  int country_pop[num_rows];
  matrix[num_rows,3] lin_counter;
  vector[2] phi_scale; // prior on how much change there could be in infection rate over time
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
  real mu_test_raw[1];
  real finding; // difficulty of identifying infected cases 
  real<lower=0> world_infect; // infection rate based on number of travelers
  vector[S] suppress_effect_raw; // suppression effect of govt. measures, cannot increase virus transmission rate
  vector[L] lockdown_effect_raw;
  vector[L] lockdown_med_raw[G];
  vector[S] suppress_med_raw[G];
  vector[L] lockdown_med_raw_fear;
  vector[S-1] suppress_med_raw_fear;
  real test_baseline;
  real mu_test_raw2[1];
  real mu_test_raw3[1];
  real pcr_spec; 
  vector[3] mu_poly; // hierarchical mean for poly coefficients
  vector[G] mob_effect_raw;
  vector<lower=0>[3] sigma_poly; // varying sigma polys
  vector[G] mob_alpha_const; // mobility hierarchical intercepts
  vector[num_country] country_test_raw; // unobserved rate at which countries are willing to test vs. number of infected
  vector[num_country] country_test_raw2;
  vector[num_country] country_test_raw3;
  real alpha_infect; // other intercepts
  real alpha_test[1];
  vector<lower=0>[2] phi_raw; // shape parameter for infected
  real<lower=0> sigma_test_raw[1]; // estimate of between-state testing heterogeneity
  real<lower=0> sigma_test_raw2[1];
  real<lower=0> sigma_test_raw3[1];
  cholesky_factor_corr[G] M_Omega; // these are for the MVN for mobility data
  vector<lower=0>[G] M_sigma;
  real fear_const;
  real<lower=0> sigma_fear;
}
transformed parameters {
  vector[2] phi;
  
  phi = (1 ./ phi_scale) .* phi_raw;
}
model {
  matrix[G, G] M_Sigma;
  int grainsize = 1;
  
  sigma_poly ~ exponential(.1);
  mu_poly ~ normal(0,30);
  mu_test_raw ~ normal(0,20);
  
  mu_test_raw2 ~ normal(0,20);
  mu_test_raw3 ~ normal(0,20);
  world_infect ~ normal(0,3);
  lockdown_effect_raw ~ normal(0,5);
  alpha_infect ~ normal(0,10); // this can reach extremely low values
  alpha_test ~ normal(0,20);

  phi_raw ~ exponential(1);
  mob_effect_raw ~ normal(0,5);
  suppress_effect_raw ~ normal(0,5);
  test_baseline ~ normal(0,20);
  
  mob_alpha_const ~ normal(0,5);
  pcr_spec ~ normal(0,20);
  
  finding ~ normal(0,20);
  sigma_test_raw ~ exponential(1);
  sigma_test_raw2 ~ exponential(1);
  sigma_test_raw3 ~ exponential(1);
  sigma_fear ~ exponential(.1);
  fear_const ~ normal(0,5);
  
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
                     phi,
                     country_pop,
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
                     //test_lin_counter,
                     country_test_raw,
                     S,
                     G,
                     L,
                     R,
                     RE,
                     sero_row,
                     sero_row_real,
                     sero,
                     sero_real,
                     mobility,
                     mob_array,
                     mob_alpha_const,
                     lockdown_med_raw,
                     suppress_med_raw,
                     suppress_med_raw_fear,
                     lockdown_med_raw_fear,
                     sigma_test_raw,
                     country_test_raw2,
                     country_test_raw3,
                     //test_lin_counter2,
                     //test_max_par,
                     //test_max,
                     lin_counter,
                     mu_test_raw,
                     mu_test_raw2,
                     mu_test_raw3,
                     sigma_test_raw2,
                     sigma_test_raw3,
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
  vector[num_rows] cov_out;
  
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
        
        //poly_nonc1 = poly1[1];
        //poly_nonc2 = poly2[2];
        //poly_nonc3 = poly3[3];
    
    for(i in 1:obs) {
              if(i==1) {
                prop_infect_out[start2] = alpha_infect + count_outbreak[start2,1] * poly_nonc1 +
                    count_outbreak[start2,2] * poly_nonc2  +
                    count_outbreak[start2,3] * poly_nonc3  +
                    world_infect*month_cases[start2] +
                    Q_supp[start2,1:S]*suppress_effect_raw +
                    Q_lock[start2,1:L]*lockdown_effect_raw +
                    Q_mob[start2,1:G]*mob_effect_raw;
                    
                 cov_out[start2] = alpha_infect + Q_mob[start2,1:G]*mob_effect_raw;
                    
              } else {
                prop_infect_out[start2 + i - 1] = exp(alpha_infect + count_outbreak[start2+i-1,1] * poly_nonc1  +
                    count_outbreak[start2+i-1,2] * poly_nonc2 +
                    count_outbreak[start2+i-1,3] * poly_nonc3 +
                    world_infect*month_cases[start2+i-1] +
                    Q_supp[start2+i-1,1:S]*suppress_effect_raw +
                    Q_lock[start2+i-1,1:L]*lockdown_effect_raw +
                    Q_mob[start2+i-1,1:G]*mob_effect_raw) + prop_infect_out[start2 + i - 2];
                    
                cov_out[start2 + i - 1] = exp(alpha_infect + Q_mob[start2+i-1,1:G]*mob_effect_raw) + cov_out[start2 + i - 2];
              }
        }
  }
          
}


