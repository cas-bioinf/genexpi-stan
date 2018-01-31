//Spline code from https://github.com/milkha/Splines_in_Stan/blob/master/b_spline.stan

functions {
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    // INPUTS:
    //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    ind:        the index of the b_spline
    //    order:      the order of the b-spline
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
             (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) /
                 (ext_knots[ind+order] - ext_knots[ind+1]);
      // Calculating the B-spline recursively as linear interpolation of two lower-order splines
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) +
                 w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }

  int count_num_basis(int num_knots, int spline_degree) {
    return num_knots + spline_degree - 1;
  }

  int regulator_expression_size(int regulator_measured, int num_measurements) {
    if (regulator_measured == 0) {
      return 0;
    } else {
      return num_measurements;
    }
  }

  int coeffs_prior_size(int coeffs_prior_given, int num_knots, int spline_degree) {
    if(coeffs_prior_given == 0) {
      return 1;
    } else {
      return count_num_basis(num_knots, spline_degree);
    }
  }
}

data {
  int num_time;
  int num_measurements;
  int num_targets;
  int<lower=1,upper=num_time> measurement_times[num_measurements];
  int<lower=0,upper=1> regulator_measured;

  vector<lower=0>[regulator_expression_size(regulator_measured, num_measurements)] regulator_expression;
  vector<lower=0>[num_measurements] expression[num_targets];
  int<lower=-1,upper=1> regulation_signs[num_targets];

  int num_knots;            // num of knots
  vector[num_knots] knots;  // the sequence of knots
  int spline_degree;        // the degree of spline (is equal to order - 1)


  int<lower=0,upper=1> coeffs_prior_given;
  vector[coeffs_prior_size(coeffs_prior_given, num_knots, spline_degree)] coeffs_prior_mean;
  cov_matrix[coeffs_prior_size(coeffs_prior_given, num_knots, spline_degree)] coeffs_prior_cov;

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
  int num_basis = count_num_basis(num_knots, spline_degree); // total number of B-splines
  matrix[num_basis, num_time] B;  // matrix of B-splines
  vector[spline_degree + num_knots] ext_knots_temp;
  vector[2*spline_degree + num_knots] ext_knots; // set of extended knots
  real time_real[num_time];
  vector[num_targets] regulation_signs_real;
  vector<lower=0>[num_targets] basal_transcription;

  cholesky_factor_cov[coeffs_prior_size(coeffs_prior_given, num_knots, spline_degree)] coeffs_prior_chol_cov = cholesky_decompose(coeffs_prior_cov);

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

  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
  for (ind in 1:num_basis) {
    B[ind,:] = to_row_vector(build_b_spline(time_real, to_array_1d(ext_knots), ind, spline_degree + 1));
  }
  B[num_basis, num_time] = 1;

  B = B * scale;

  for(t in 1:num_targets) {
    basal_transcription[t] = 0;
  }
}


parameters {
  row_vector[num_basis] coeffs;
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
    regulator_profile = to_vector(coeffs * B) ;
  } else {
    row_vector[num_basis] coeffs_transformed = to_row_vector(coeffs_prior_chol_cov * to_vector(coeffs) + coeffs_prior_mean);

    regulator_profile = to_vector(coeffs_transformed * B);
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
