functions {
  int coeffs_prior_size(int coeffs_prior_given, int num_spline_basis) {
    return coeffs_prior_given ? num_spline_basis : 1;
  }


}

data {
  int num_time;
  int num_measurements;
  int num_targets;
  int<lower=1,upper=num_time> measurement_times[num_measurements];
  int<lower=0,upper=1> regulator_measured;

  vector<lower=0>[regulator_measured ? num_measurements : 0] regulator_expression;
  vector<lower=0>[num_measurements] expression[num_targets];
  int<lower=-1,upper=1> regulation_signs[num_targets];

  int num_spline_basis;
  matrix[num_spline_basis, num_time] spline_basis;  // matrix of B-splines


  int<lower=0,upper=1> coeffs_prior_given;
  vector[coeffs_prior_size(coeffs_prior_given, num_spline_basis)] coeffs_prior_mean;
  cov_matrix[coeffs_prior_size(coeffs_prior_given, num_spline_basis)] coeffs_prior_cov;

  real<lower=0> measurement_sigma_relative;
  real<lower=0> measurement_sigma_absolute;

  real<lower = 0> initial_condition_prior_sigma;
  real<lower = 0> sensitivity_prior_sigma;
  real degradation_prior_mean;
  real<lower = 0> degradation_prior_sigma;
  real<lower = 0> w_prior_sigma;
  real<lower = 0> b_prior_sigma;
  real<lower = 0> intercept_prior_sigma;

  real<lower=0> scale;
}

transformed data {
  matrix[num_spline_basis, num_time] spline_basis_scaled;  // matrix of B-splines
  real time_real[num_time];
  vector[num_targets] regulation_signs_real;
  vector<lower=0>[num_targets] basal_transcription;

  cholesky_factor_cov[coeffs_prior_size(coeffs_prior_given, num_spline_basis)] coeffs_prior_chol_cov = cholesky_decompose(coeffs_prior_cov);

  if(regulator_measured != 0 && coeffs_prior_given != 0) {
    reject("You are not supposed to give coeffs prior if regulator is measured. Use scale parameter instead.");
  }

  for(i in 1:num_time) {
    time_real[i] = i;
  }

  for(t in 1:num_targets) {
    if(regulation_signs[t] == 0) {
      reject("regulation_signs have to be nonzero");
    }
    regulation_signs_real[t] = regulation_signs[t];
  }

  {
    for(m in 2:num_measurements) {
      if(measurement_times[m - 1] >= measurement_times[m]) {
        reject("Measurement times have to be increasing");
      }
    }
  }


  spline_basis_scaled = spline_basis * scale;

  for(t in 1:num_targets) {
    basal_transcription[t] = 0;
  }
}


parameters {
  row_vector[num_spline_basis] coeffs;
  vector<lower=0>[num_targets] initial_condition;
  vector<lower=0>[num_targets] sensitivity;
  vector<lower=0>[num_targets] degradation;
  vector<lower=0>[num_targets] w_raw;
  vector[num_targets] b_raw;
  //real<lower=0> degradation_over_sensitivity;
  real<lower = 0> intercept_raw;
}


transformed parameters {
  vector[num_time] regulator_profile;
  matrix<lower=0>[num_targets,num_time] predicted_expression;
  vector<lower=0>[num_time] predicted_regulator_expression;
  vector[num_targets] b;
  vector[num_targets] w;
  real intercept;

  w = w_raw .* regulation_signs_real;
  if(coeffs_prior_given == 0) {
    regulator_profile = to_vector(coeffs * spline_basis_scaled) ;
  } else {
    row_vector[num_spline_basis] coeffs_transformed = to_row_vector(coeffs_prior_chol_cov * to_vector(coeffs) + coeffs_prior_mean);

    regulator_profile = to_vector(coeffs_transformed * spline_basis_scaled);
  }
  {
    real min_intercept = -min(regulator_profile);
    if(regulator_measured == 0) {
      intercept = min_intercept;
    } else {
      intercept = min_intercept + intercept_raw;
    }
    predicted_regulator_expression = regulator_profile + intercept;
  }
  b = b_raw + (intercept * w);
  for(t in 1:num_targets)
  {
    real basal_over_degradation = basal_transcription[t] / degradation[t];
    real degradation_per_step = exp(-degradation[t]);

    real residual;
    vector[num_time] synthesis;

    synthesis = inv(1 + exp(-regulator_profile * w[t] - b_raw[t]));

    predicted_expression[t,1] = initial_condition[t];

    //Calculating the integral by trapezoid rule in a single pass for all values
    residual = -0.5 * synthesis[1];
    for (measurement in 2:num_measurements)
    {
      for (time in (measurement_times[measurement - 1] + 1):measurement_times[measurement])
      {
        residual = (residual + synthesis[time - 1]) * degradation_per_step;

        { //new block to let me define new vars
          real integral_value = residual + 0.5 * synthesis[time];
          predicted_expression[t,time] = basal_over_degradation + (initial_condition[t] - basal_over_degradation) * exp(-degradation[t] * (time - 1)) + sensitivity[t] * integral_value;
        }

      }
    }
  }
}

model {


    //Observation model
    for(m in 1:num_measurements) {
      vector[num_targets] sigma = measurement_sigma_absolute + measurement_sigma_relative * predicted_expression[,measurement_times[m]];
      for(t in 1:num_targets) {
        expression[t, m] ~ normal(predicted_expression[t,measurement_times[m]], sigma[t]) T[0,];
      }

      if(regulator_measured != 0) {
        real sigma_regulator = measurement_sigma_absolute + measurement_sigma_relative * predicted_regulator_expression[measurement_times[m]];
        regulator_expression[m] ~ normal(predicted_regulator_expression[measurement_times[m]], sigma_regulator)  T[0,];
      }
      //log(expression[m]) ~ normal(log(predicted_expression[measurement_times[m]]), measurement_sigma);
    }

    initial_condition ~ normal(0, initial_condition_prior_sigma);
    coeffs ~ normal(0,1); //coeffs are rescaled by the scale parameter
    //scale ~ normal(0,scale_prior_sigma);
    intercept_raw ~ normal(0,intercept_prior_sigma);
    sensitivity ~ normal(0, sensitivity_prior_sigma);
    //degradation_over_sensitivity ~ lognormal(0,1);
    degradation ~ lognormal(degradation_prior_mean, degradation_prior_sigma);
    w_raw ~ normal(0,w_prior_sigma);
    b_raw ~ normal(0,b_prior_sigma);
}

generated quantities {
  vector<lower=0>[num_measurements] expression_replicates[num_targets];
  vector[num_measurements] log_likelihood[num_targets];

  for(m in 1:num_measurements) {
    vector[num_targets] sigma = measurement_sigma_absolute + measurement_sigma_relative * predicted_expression[,measurement_times[m]];
    for(t in 1:num_targets) {
      //Draw from the truncated normal
      real lower_bound = 0;
      real p = normal_cdf(lower_bound, predicted_expression[t,measurement_times[m]], sigma[t]);
      real u = uniform_rng(p, 1);
      expression_replicates[t,m] = inv_Phi(u) * sigma[t] + predicted_expression[t,measurement_times[m]];

      log_likelihood[t,m] = normal_lpdf(expression[t,m]| predicted_expression[t,measurement_times[m]], sigma[t]) - log_diff_exp(1, normal_lcdf(lower_bound | predicted_expression[t,measurement_times[m]], sigma[t]));
      /*
      while (1) {
        expression_replicates[t,m] = normal_rng(predicted_expression[t,measurement_times[m]], sigma[t]);
        if(expression_replicates[t,m] >= 0) {
          break;
        }
      }*/
    }
  }
}
