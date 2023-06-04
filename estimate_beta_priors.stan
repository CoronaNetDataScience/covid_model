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
  real<lower=0> successes;
  real<lower=0> failures;
}

model {
  
  //very weakly informative
  successes ~ exponential(.01);
  failures ~ exponential(.01);
  y ~ beta(successes, failures);

}
generated quantities {
  
  real denominator = successes + failures;
  
}
