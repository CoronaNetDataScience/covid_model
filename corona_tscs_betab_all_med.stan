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
    
  
  int r_in2(int pos,int[] pos_var) {
    
    for (p in 1:(size(pos_var))) {
       int i = 0;
       int prev;
       if (pos_var[p]==pos) {
       // can return immediately, as soon as find a match
          i += 1;
          if(p<size(pos_var)) {
            prev = pos_var[p];
            // remaining elements to check
            if(i==1) {
            // keep going
            } else {
              if(prev==pos_var[p] && p < size(pos_var)) {
                // keep going
                continue;
              } else {
                if(prev==pos_var[p]) {
                  return p;
                } else {
                  // sequence ended last iteration
                  return p - 1;
                }
                
              }
              
            }
          } else {
            // match was at the end of the sequence
            return p;
          }
          
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
                   vector mob_alpha_const,
                   vector[] lockdown_med_raw,
                   vector[] suppress_med_raw,
                   vector sigma_med,
                   real sigma_test_raw,
                   vector country_test_raw2,
                   real test_lin_counter2,
                   vector sero_est,
                   real test_max_par,
                   vector test_max,
                   matrix lin_counter,
                   real mu_test_raw,
                   real mu_test_raw2,
                   real sigma_test_raw2) {
                     
    // big loop over states
    real log_prob = 0;
    for(r in 1:size(y_slice)) {
        
        int s = y_slice[r,1];
        int start2 = y_slice[r,2];
        int end2 = y_slice[r,3];
      
        vector[end2 - start2 + 1] prop_infected; // modeled infection rates for domestic transmission\
        
        int obs = end2 - start2 + 1;
        //int min_cc = min(this_cc2);
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
        
        // increasing constraint
        
        //log_prob += exponential_lpdf(prop_infected[2:obs] - prop_infected[1:(obs-1)]|1);
        
        for(g in 1:G)
          log_prob += normal_lpdf(to_vector(mobility[start2:end2,g])|mob_alpha_const[g] +
          Q_lock[start2:end2,1:L]*lockdown_med_raw[g]  +
          Q_supp[start2:end2,1:S]*suppress_med_raw[g],sigma_med[g]);
        
        
    mu_tests = inv_logit(alpha_test + 
    //test_max_par*test_max[start2:end2] +
                          test_baseline *  mix_prop_std +
                          country1 * mix_prop_std .* lin_counter[start2:end2,1] +
                          test_lin_counter * lin_counter[start2:end2,1]);
                          // country2 * mix_prop_std .* lin_counter[start2:end2,2] +
                          // test_lin_counter2 * lin_counter[start2:end2,2]);
                          //country_test_raw[s] * prop_infected .* lin_counter[start2:end2] +
                          //
                              //test_lin_counter * lin_counter[start2:end2] +
                              //test_max_par * count_outbreak[start2:end2,2] .* prop_infected +

        //mu_tests = inv_logit(alpha[1] + country_test_raw[this_cc] .* mix_prop_std);
        log_prob += beta_binomial_lpmf(tests[start2:end2]|country_pop[start2:end2],mu_tests*phi[1],(1-mu_tests) * phi[1]);
       log_prob += beta_binomial_lpmf(cases[start2:end2]|country_pop[start2:end2],mu_cases*phi[2],(1-mu_cases) * phi[2]);
    
        // loop over serology surveys to add informative prior information
        
        for(n in start2:end2) {
          int q = r_in(n,sero_row);
          if(q>0) {
            //print("beta prior");
            //print(log_prob);
            //log_prob += beta_lpdf(inv_logit(prop_infected[n-start2+1])|.5 + sero[r,1]*sero[r,3],.5 + sero[r,3] - (sero[r,3]*sero[r,1]));  
            //print(log_prob);
            log_prob += normal_lpdf(prop_infected[n-start2+1]|logit(sero[q,1]),.01);
            log_prob += log(prop_infected[n-start2+1] - prop_infected[n-start2]);
            //log_prob += log(prop_infected[n-start2+1] - prop_infected[n-start2]) - 2 * log(1 + prop_infected[n-start2+1] - prop_infected[n-start2]);
          } else {
            // prior implies inflation rate (cases/true infected) is between 2 - 20
            //print("normal prior");
            //print(log_prob);
            if((n-start2)>60) {
                log_prob += lognormal_lpdf(inv_logit(prop_infected[n-start2+1])/((cases[n-start2+1]*1.0)/(country_pop[n-start2+1]*1.0))|2.1,.4);
                log_prob +=  log(prop_infected[n-start2+1] - prop_infected[n-start2]) - 2 * log(1 + prop_infected[n-start2+1] - prop_infected[n-start2]);
                //print(log_prob);
                //log_prob += ;
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
  matrix[num_rows,G] mobility; // time-varying mobility measures
  matrix[num_rows,L] lockdown; // hierachical lockdown predictors
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
  real mu_test_raw;
  vector<lower=0,upper=1>[R] sero_est;
  real<lower=0> finding; // difficulty of identifying infected cases 
  //vector<lower=0,upper=1>[R] survey_prop; // variable that stores survey proportions from CDC data
  real<lower=0> world_infect; // infection rate based on number of travelers
  vector[S] suppress_effect_raw; // suppression effect of govt. measures, cannot increase virus transmission rate
  vector[L] lockdown_effect_raw;
  vector[L] lockdown_med_raw[G];
  vector[S] suppress_med_raw[G];
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
  vector<lower=0>[G] sigma_med;
  vector[num_country] country_test_raw; // unobserved rate at which countries are willing to test vs. number of infected
  vector[num_country] country_test_raw2;
  // we assume that as infection rates increase, more tests will be conducted
  real alpha_infect; // other intercepts
  //real<upper=-5> alpha_test;
  real alpha_test;
  vector<lower=0>[2] phi; // shape parameter for infected
  real<lower=0> sigma_test_raw; // estimate of between-state testing heterogeneity
  real<lower=0> sigma_test_raw2;
}
model {
  
  int grainsize = 1;
  
  sigma_poly ~ normal(0,5);
  mu_poly ~ normal(0,10);
  mu_test_raw ~ normal(0,20);
  mu_test_raw2 ~ normal(-5,2);
  world_infect ~ normal(0,3);
  lockdown_effect_raw ~ normal(0,5);
  alpha_infect ~ normal(0,10); // this can reach extremely low values
  alpha_test ~ normal(0,10);

  phi ~ exponential(phi_scale);
  mob_effect_raw ~ normal(0,5);
  suppress_effect_raw ~ normal(0,5);
  test_max_par ~ normal(0,5);
  test_baseline ~ normal(0,5);
  test_lin_counter ~ normal(0,5);
  test_lin_counter2 ~ normal(0,5);
  
  mob_alpha_const ~ normal(0,5);
  pcr_spec ~ normal(0,10);
  
  finding ~ normal(30,10);
  sigma_test_raw ~ normal(0,5);
  sigma_test_raw2 ~ normal(0,.1);
  sigma_med ~ exponential(.1);
  //country_test_raw2 ~ normal(0,1);
  sero_est ~ beta(.5 + sero[,1] .* sero[,3],.5 + sero[,3] - (sero[,3] .* sero[,1]));

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
                     mob_alpha_const,
                     lockdown_med_raw,
                     suppress_med_raw,
                     sigma_med,
                     sigma_test_raw,
                     country_test_raw2,
                     test_lin_counter2,
                     sero_est,
                     test_max_par,
                     test_max,
                     lin_counter,
                     mu_test_raw,
                     mu_test_raw2,
                     sigma_test_raw2);

}

generated quantities {
  
  // convert QR estimates back to actual numbers
  
  vector[S] suppress_effect; // suppression effect of govt. measures, cannot increase virus transmission rate
  vector[L] lockdown_effect;
  vector[G] mob_effect;
  vector[L] lockdown_med[G];
  vector[S] suppress_med[G];
  vector[num_rows] prop_infect_out;
  // vector[num_rows] out_infected;
  // 
  // vector[num_country] poly_nonc1; // non-centered poly parameters
  // vector[num_country] poly_nonc2; // non-centered poly parameters
  // vector[num_country] poly_nonc3; // non-centered poly parameters
  // 
  // poly_nonc1 = mu_poly[1] + sigma_poly[1]*poly1;
  // poly_nonc2 = mu_poly[2] + sigma_poly[2]*poly2;
  // poly_nonc3 = mu_poly[3] + sigma_poly[3]*poly3;
  
  suppress_effect = R_supp_inverse * suppress_effect_raw;
  lockdown_effect = R_lock_inverse * lockdown_effect_raw;
  mob_effect = R_mob_inverse * mob_effect_raw;
  
  for(g in 1:G) {
    lockdown_med[g] = R_lock_inverse * lockdown_med_raw[g];
    suppress_med[g] = R_supp_inverse * suppress_med_raw[g];
  }
    
  
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
          
  
          
    
  // out_infected = alpha[2] + count_outbreak[,1] .* poly_nonc1[cc]  +
  //         count_outbreak[,2] .* poly_nonc2[cc] +
  //         count_outbreak[,3] .* poly_nonc3[cc] +
  //         world_infect*month_cases +
  //         Q_supp[,1:S]*suppress_effect_raw +
  //         Q_lock[,1:L]*lockdown_effect_raw +
  //         Q_mob[,1:G]*mob_effect_raw;
  // 
}

