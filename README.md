# README

This Github repository contains code and data to reproduce results in the article:

Robert Kubinec, Luiz Max Carvalho, Joan Barceló, Cindy Cheng, Luca Messerschmidt and Matthew Sean Cottrell. "A Bayesian latent variable model for the optimal identification of disease incidence rates given information constraints." *Journal of the Royal Statistical Society Series A: Statistics in Society*. 2024. https://doi.org/10.1093/jrsssa/qnae040 .

A brief description of the files is found below. If you have any questions about the information in the repo, please contact Bob Kubinec at `bobkubinec@gmail.com`.

First, note that the paper relies on a fitted `cmdstanr` model to reproduce results. These model fits are too big to store on Github, but you can access them from this Google drive folder and place them in the `data` sub-folder to reproduce results without fitting models (may take up to a few days):

https://drive.google.com/drive/folders/1hVzD_qL1CnOkTkwI6VH1PEgC1LK44RS1?usp=sharing

## Paper files:

1. `kubinec_model_preprint.Rmd`: This file contains the text and embedded R code to reproduce the figures and tables in the paper.

2. `kubinec_model_SI.Rmd`: This file contains the text and embedded R code to reproduce the supplementary information.

## Code:

1. `corona_tscs_betab_mix_prior_v2.stan` This Stan file contains the code to fit the model described in the paper using Stan (specifically, `cmdstan` accessed via the `cmdstanr` package). See code in the `kubinec_model_preprint.Rmd` file to see how to fit the model from R.

2. `estimate_beta_priors_v2.stan`: This Stan file calculates the uncertainty of the empirical distributions of the estimates of the expert survey about COVID-19 incidence in the early pandemic period.

## Data:

1. `data/combined.rds`: the combined dataset with COVID-19 cases, tests, Census data and expert and serology survey data

2. `nyt_data.rds`, `goog_mobile.rds` and `tests.rds`: New York Times (reported cases), Google mobility data and testing data for the time period described in the paper. Note that the paper code can download these from Github repositories, but these sources may no longer be available.

3. `data/simulation/`: contains masking and Civiqs polls about COVID-19 related fears and behaviors.

4. `data/covid_amp_state_policy_data.xlsx`: contains COVID-AMP state-level policy data as described in the paper

5. `count_pol_covidamp.rds`: Aggregated form of COVID-AMP data to the state level as a count of policies.

6. `cdc_sample_sizes.csv`: CDC serology surveys

7. `data/consensusForecastsDB.csv"`: expert survey of epidemiologists during the early pandemic period

8. `data/rhat_summaries*.rds` Rhat summaries for different models as reported in the paper.



