//
// This Stan program estimates the parameters of the 
// Beta distribution given an input vector of proportions
// Calculate effective sample size given data we have on proportions
// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] y;
}

parameters {
  real<lower=0,upper=1> mu;
  real<lower=0> kappa;
}

model {
  
  //very weakly informative
  mu ~ normal(0,10);
  kappa ~ exponential(.01);
  y ~ beta_proportion(mu, kappa);

}
generated quantities {
  
  // figure out the denominator of a hypothetical survey
  real alpha;
  real beta;
  real denominator;
  
  alpha = mu*kappa;
  beta = kappa - alpha;
  denominator = alpha + beta;
}
