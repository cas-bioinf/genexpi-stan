data {
  int num_measurements;
  int num_time;
  int<lower=0,upper=num_time> measurement_times[num_measurements];
  vector<lower=0>[num_measurements] expression;
  real<lower=0> gp_length;
  real<lower=0> gp_variance;

  real<lower=0> measurement_sigma;

  real<lower=0> initial_condition;
  real<lower=0> degradation;
  real<lower=0> sensitivity;

}


transformed data {
    matrix[num_time, num_time] covariance;
    matrix[num_time, num_time] covariance_cholesky;
    real time_real[num_time];
    real basal_transcription = 0;

    for(i in 1:num_time) {
      time_real[i] = i;
    }

    covariance = cov_exp_quad(time_real, gp_variance, gp_length);

    for (I in 1:num_time)
    {
      covariance[I, I] = covariance[I, I] + 1e-12;
    }

    covariance_cholesky = cholesky_decompose(covariance);


}

parameters {
  vector[num_time] regulator_profile_gp;
}


transformed parameters {
  vector[num_time] regulator_profile;
  vector[num_measurements] predicted_expression;

  regulator_profile = covariance_cholesky * regulator_profile_gp;
  {
    //Compute the predicted expression
    real basal_over_degradation = basal_transcription / degradation;
    real degradation_per_step = exp(-degradation);

    real residual;
    vector[num_time] synthesis;

    //synthesis = inv(1 + exp(-regulator_profile));
    synthesis = regulator_profile + 3;

    predicted_expression[1] = initial_condition;

    //Calculating the integral by trapezoid rule in a single pass for all values
    residual = -0.5 * synthesis[1];
    for (measurement in 2:num_measurements)
    {
      for (time in (measurement_times[measurement - 1] + 1):measurement_times[measurement])
      {
        residual = (residual + synthesis[time - 1]) * degradation_per_step;
      }

      { //new block to let me define new vars
        int time = measurement_times[measurement];
        real integral_value = residual + 0.5 * synthesis[time];
        predicted_expression[measurement] = basal_over_degradation + (initial_condition - basal_over_degradation) * exp(-degradation * (time - 1)) + sensitivity * integral_value;
      }

    }

  }
}

model {


    //Observation model
    expression ~ normal(predicted_expression, measurement_sigma);

    //GP prior
    regulator_profile_gp ~ normal(0,1);
}
