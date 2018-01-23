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
}

data {
  int num_time;
  int num_measurements;
  vector<lower=0>[num_measurements] expression;
  int<lower=0,upper=num_time> measurement_times[num_measurements];

  int num_knots;            // num of knots
  vector[num_knots] knots;  // the sequence of knots
  int spline_degree;        // the degree of spline (is equal to order - 1)

  real<lower=0> measurement_sigma;

  real<lower = 0> scale_prior_sigma;
  real<lower = 0> sensitivity_prior_sigma;
  real<lower = 0> degradation_prior_sigma;

  real<lower=0> scale;
}

transformed data {
  int num_basis = num_knots + spline_degree - 1; // total number of B-splines
  matrix[num_basis, num_time] B;  // matrix of B-splines
  vector[spline_degree + num_knots] ext_knots_temp;
  vector[2*spline_degree + num_knots] ext_knots; // set of extended knots
  real time_real[num_time];
  real basal_transcription = 0;

  for(i in 1:num_time) {
    time_real[i] = i;
  }


  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
  for (ind in 1:num_basis) {
    B[ind,:] = to_row_vector(build_b_spline(time_real, to_array_1d(ext_knots), ind, spline_degree + 1));
  }
  B[num_basis, num_time] = 1;
}


parameters {
  row_vector[num_basis] coeffs;
  real<lower=0> initial_condition;
  real<lower=0> sensitivity;
  //real<lower=0> degradation_over_sensitivity;
  real intercept;
  real<lower=0> degradation;



  //real<lower=0> spline_variance;
}


transformed parameters {
  vector[num_time] regulator_profile;
  vector[num_time] predicted_expression;
  //real<lower=0> degradation = degradation_over_sensitivity * sensitivity;

  regulator_profile = to_vector(coeffs*B) * scale + intercept;
  {
    real basal_over_degradation = basal_transcription / degradation;
    real degradation_per_step = exp(-degradation);

    real residual;
    vector[num_time] synthesis;

    synthesis = inv(1 + exp(-regulator_profile));

    predicted_expression[1] = initial_condition;

    //Calculating the integral by trapezoid rule in a single pass for all values
    residual = -0.5 * synthesis[1];
    for (measurement in 2:num_measurements)
    {
      for (time in (measurement_times[measurement - 1] + 1):measurement_times[measurement])
      {
        residual = (residual + synthesis[time - 1]) * degradation_per_step;

        { //new block to let me define new vars
          real integral_value = residual + 0.5 * synthesis[time];
          predicted_expression[time] = basal_over_degradation + (initial_condition - basal_over_degradation) * exp(-degradation * (time - 1)) + sensitivity * integral_value;
        }

      }
    }

/*
    predicted_expression[1] = initial_condition;
    for (i in 2:num_time)
    {
      predicted_expression[i] = predicted_expression[i - 1]* (1 - degradation) + sensitivity/( 1 + exp(-regulator_profile[i]));
    }
*/
  }
}

model {


    //Observation model
    for(m in 1:num_measurements) {
      log(expression[m]) ~ normal(log(predicted_expression[measurement_times[m]]), measurement_sigma);
    }

    initial_condition ~ normal(0,1);
    coeffs ~ normal(0,1);
    //scale ~ normal(0,scale_prior_sigma);
    intercept ~ normal(0,1);
    sensitivity ~ normal(0, sensitivity_prior_sigma);
    //degradation_over_sensitivity ~ lognormal(0,1);
    degradation ~ normal(0, degradation_prior_sigma);
}
