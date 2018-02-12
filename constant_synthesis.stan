data {
  int num_measurements;
  int<lower=1> measurement_times[num_measurements];
  vector<lower=0>[num_measurements] expression;

  real<lower=0> measurement_sigma_relative;
  real<lower=0> measurement_sigma_absolute;

  real<lower = 0> initial_condition_prior_sigma;
  real<lower = 0> synthesis_over_degradation_prior_sigma;
  real synthesis_over_degradation_prior_mean;
  real degradation_prior_mean;
  real<lower = 0> degradation_prior_sigma;
}


parameters {
  real<lower=0> initial_condition;
  real<lower=0> synthesis_over_degradation;
  real<lower=0> degradation;
}


transformed parameters {
  vector<lower=0>[num_measurements] predicted_expression;

  {
    for(t in 1:num_measurements) {
      predicted_expression[t] = synthesis_over_degradation + exp(-degradation * measurement_times[t]) * (initial_condition - synthesis_over_degradation);
    }
  }
}

model {


    //Observation model
    for(m in 1:num_measurements) {
      real sigma = measurement_sigma_absolute + measurement_sigma_relative * predicted_expression[m];
      expression[m] ~ normal(predicted_expression[m], sigma) T[0,];
    }

    initial_condition ~ normal(0, initial_condition_prior_sigma);
    synthesis_over_degradation ~ normal(synthesis_over_degradation_prior_mean, synthesis_over_degradation_prior_sigma) T[0,];
    degradation ~ lognormal(degradation_prior_mean, degradation_prior_sigma);
}

generated quantities {
  vector<lower=0>[num_measurements] expression_replicates;
  vector[num_measurements] log_likelihood;

  for(m in 1:num_measurements) {
    real sigma = measurement_sigma_absolute + measurement_sigma_relative * predicted_expression[m];
    //Draw from the truncated normal
    real lower_bound = 0;
    real p = normal_cdf(lower_bound, predicted_expression[m], sigma);
    real u = uniform_rng(p, 1);
    expression_replicates[m] = inv_Phi(u) * sigma + predicted_expression[m];

    log_likelihood[m] = normal_lpdf(expression[m]| predicted_expression[m], sigma) - log_diff_exp(1, normal_lcdf(lower_bound | predicted_expression[m], sigma));
  }
}
