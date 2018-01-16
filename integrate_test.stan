data {
  int N;
  vector[N] expression;
  real<lower=0> gp_length;

  real<lower=0> measurement_sigma;
  real<lower=0> degradation;
}



parameters {
  vector[N] regulator_profile_gp;
  real<lower=0> gp_variance;
}


transformed parameters {
  vector[N] regulator_profile;
  vector[N] predicted_expression;

  {
    matrix[N, N] covariance;
    matrix[N, N] covariance_cholesky;
    real time_real[N];

    for(i in 1:N) {
      time_real[i] = i;
    }

    covariance = cov_exp_quad(time_real, gp_variance, gp_length);

    for (I in 1:N)
    {
      covariance[I, I] = covariance[I, I] + 1e-12;
    }

    covariance_cholesky = cholesky_decompose(covariance);

    regulator_profile = covariance_cholesky * regulator_profile_gp;
    predicted_expression[1] = regulator_profile[1];
    for (i in 2:N)
    {
      predicted_expression[i] = predicted_expression[i - 1]* (1 - degradation) + inv(exp(-regulator_profile[i]));
    }

  }
}

model {


    //Observation model
    expression ~ normal(predicted_expression, measurement_sigma);

    //GP prior
    regulator_profile_gp ~ normal(0,1);

    gp_variance ~ normal(0,1);
}
