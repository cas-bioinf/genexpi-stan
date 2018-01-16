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
  int N;
  vector<lower=0>[N] expression;

  int num_knots;            // num of knots
  vector[num_knots] knots;  // the sequence of knots
  int spline_degree;        // the degree of spline (is equal to order - 1)

  real<lower=0> measurement_sigma;

  real<lower = 0> scale_prior_sigma;
  real<lower = 0> sensitivity_prior_sigma;
  real intercept;

  real<lower=0> scale;
  real min_scale;

}

transformed data {
  int num_basis = num_knots + spline_degree - 1; // total number of B-splines
  matrix[num_basis, N] B;  // matrix of B-splines
  vector[spline_degree + num_knots] ext_knots_temp;
  vector[2*spline_degree + num_knots] ext_knots; // set of extended knots
  real time_real[N];

  for(i in 1:N) {
    time_real[i] = i;
  }


  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
  for (ind in 1:num_basis) {
    B[ind,:] = to_row_vector(build_b_spline(time_real, to_array_1d(ext_knots), ind, spline_degree + 1));
  }
  B[num_basis, N] = 1;
}


parameters {
  row_vector[num_basis] coeffs;
  real<lower=0> initial_condition;
  real<lower=0> sensitivity;
  real<lower=0> degradation_over_sensitivity;
  //real<lower=0> degradation;



  //real<lower=0> spline_variance;
}


transformed parameters {
  vector[N] regulator_profile;
  vector[N] predicted_expression;
  real<lower=0> degradation = degradation_over_sensitivity * sensitivity;

  {

    regulator_profile = to_vector(coeffs*B) * (scale + min_scale) + intercept;
    predicted_expression[1] = initial_condition;
    for (i in 2:N)
    {
      predicted_expression[i] = predicted_expression[i - 1]* (1 - degradation) + sensitivity/( 1 + exp(-regulator_profile[i]));
    }

  }
}

model {


    //Observation model
    log(expression) ~ normal(log(predicted_expression), measurement_sigma);

    initial_condition ~ normal(0,1);
    coeffs ~ normal(0,1);
    scale ~ normal(0,scale_prior_sigma);
    intercept ~ normal(0,1);
    sensitivity ~ normal(0, sensitivity_prior_sigma);
    degradation_over_sensitivity ~ lognormal(0,1);
    //degradation ~ normal(0,1);
}
