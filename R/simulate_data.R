target_ODE <- function(t, state, parameters, protein_transform = identity)
{
  with(as.list(c(state, parameters)), {
    regulatoryInput = bias + weight * protein_transform(protein(t));
    #dX = basal_transcription + (sensitivity/(1 + exp(-regulatoryInput))) - degradation * x;
    dX = basal_transcription + sensitivity * (regulatoryInput + 3) - degradation * x;
    list(dX)
  })
}

gp_covariance <- function(distance, gp_scale, gp_length, periodic = FALSE, period = 1) {
  if(periodic) {
    return((gp_scale^2) * exp(-2*(sin(3.141592 * abs(distance) / period)^2) / (gp_length ^ 2)))
  } else {
    return((gp_scale^2) * exp( (-0.5 / (gp_length ^ 2)) * (distance ^ 2) ));
  }

}


generate_random_profile <- function(time, scale, length, mean_func = 0, periodic = FALSE, period = 1, positive_transform = TRUE) {
  # Construct the squared exponential covariance matrix
  cov_matrix = array(0,c(length(time), length(time)));
  maxTime = max(time)
  minTime = min(time)
  for(i in 1:length(time))  {
    cov_matrix[i,i] = scale^2 + 0.00000001; #Adding jitter to ensure positive semi-definiteness
    if (i < length(time)) {
      for(j in (i+1):length(time)) {
        distance = (time[i] - time[j])
        covariance = gp_covariance(distance, scale, length, periodic, period)
        cov_matrix[i,j] = covariance
        cov_matrix[j,i] = covariance
      }
    }
  }
  # Draw from the multinormal distribution using the cholesky decomposition of the cov. matrix
  chol_cov = t(chol(cov_matrix));
  raw_profile = chol_cov %*% rnorm(length(time)) + mean_func;
  # Transform to strictly positive values
  if(positive_transform) {
    positive_profile = log1p(exp(raw_profile))
    return(t(positive_profile))
  } else {
    return(t(raw_profile))
  }
}

plot_random_profiles <- function(n, time, scale, length, true_time = time, true_profile = NULL, mean_func = 0, periodic = FALSE, period = 1, positive_transform = TRUE) {
  profiles = array(0, c(n, length(time)));
  for(i in 1:n) {
    profiles[i,] = generate_random_profile(time, scale, length, mean_func = mean_func, periodic = periodic, period = period, positive_transform = positive_transform);
  }


  result = ggmatplot(time, t(profiles))
  if(!is.null(true_profile)) {
    result = result + geom_point(data = data.frame(x = true_time, y = true_profile), aes(x=x, y=y))
  }
  return(result)
}

simulate_spline <- function(num_time, num_knots, measurement_times, measurement_sigma_absolute, measurement_sigma_relative, integrate_ode45 = TRUE) {
  time <- 1:num_time;
  #spline_variance <- rnorm(1,0,1)

  #TODO: initial condition

  spline_degree = 3
  knots <- seq(from = 1, to = num_time, length.out = num_knots)
  spline_basis <- bs(1:num_time, knots = knots, degree = spline_degree)

  num_coeff <- dim(spline_basis)[2] - 1
  if(any(spline_basis[,num_coeff + 1] != 0)) {
    stop("Broken assumption of zero column in bs output")
  }
  spline_basis = spline_basis[,1:num_coeff]

  coeffs <- rnorm(num_coeff, 0, 1)
  intercept <- rnorm(1, 0, 1)

  scale_prior_sigma = 5
  #scale <- abs(rnorm(1, 0, scale_prior_sigma))
  #min_scale <- 0.5
  scale <- 5

  regulator_profile <- array((spline_basis %*% coeffs)  * scale + intercept, num_time)

  sensitivity_prior_sigma <-  1
  sensitivity <- abs(rnorm(1, 0, sensitivity_prior_sigma))

  #degradation_over_sensitivity <- rlnorm(1,0,1)
  #degradation <- degradation_over_sensitivity * sensitivity;
  degradation_prior_sigma <- 0.3
  degradation <- abs(rnorm(1,0, degradation_prior_sigma))



  expression_true <- array(-1, num_time)

  initial_condition_prior_sigma <- 1
  initial_condition <- abs(rnorm(1,0,initial_condition_prior_sigma))

  if(integrate_ode45){
    params <- c(degradation = degradation, bias = 0, sensitivity = sensitivity, weight = 1, basal_transcription = 0, protein = approxfun(time, regulator_profile, rule=2));
    expression_true <-  ode( y = c(x = initial_condition), times = time, func = target_ODE, parms = params, method = "ode45")[,"x"];
  } else {
    expression_true <- numerical_integration(0, degradation, initial_condition, sensitivity, weight = 1, bias = 0, regulator_profile, num_time)
  }

  #expression_observed <- exp(rnorm(length(measurement_times), log(expression_true[measurement_times]), measurement_sigma))
  expression_observed <- rtruncnorm(length(measurement_times), mean = expression_true[measurement_times], sd =  measurement_sigma_absolute + measurement_sigma_relative * expression_true[measurement_times], a = 0)

  #Handle regulator
  w_prior_sigma <- 1
  b_prior_sigma <- 5
  inv_w <- rlnorm(1,0,w_prior_sigma)

  # if (w > 0) {
  #   max_b <- min(expression_true)
  #   b <- rtruncnorm(1, mean = 0, sd = b_prior_sigma, b = max_b)
  # } else {
  #   min_b <- max(expression_true)
  #   b <- rtruncnorm(1, mean = 0, sd = b_prior_sigma, a = min_b)
  # }
  min_b <- -min(inv_w * regulator_profile)
  b_raw <- abs(rnorm(1,0,b_prior_sigma))
  regulator_expression_true <- regulator_profile * inv_w + min_b + b_raw
  w <- 1/inv_w;
  b <- -(min_b + b_raw) / inv_w;
  regulator_expression_observed <- rtruncnorm(length(measurement_times), mean = regulator_expression_true[measurement_times], sd =  measurement_sigma_absolute + measurement_sigma_relative * regulator_expression_true[measurement_times], a = 0)

  return(list(
    true = list (
      coeffs = coeffs,
      regulator_profile = regulator_profile,
      regulator_expression = regulator_expression_true,
      w = w,
      b = b,
      initial_condition = initial_condition,
      sensitivity = sensitivity,
      expression = expression_true,
      intercept = intercept,
      degradation = degradation
    ), observed = list(
      num_time = num_time,
      num_measurements = length(measurement_times),
      measurement_times = measurement_times,
      num_knots = num_knots,
      knots = knots,
      spline_degree = spline_degree,
      expression = expression_observed,
      regulator_expression = regulator_expression_observed,
      measurement_sigma_absolute = measurement_sigma_absolute,
      measurement_sigma_relative = measurement_sigma_relative,
      initial_condition_prior_sigma = initial_condition_prior_sigma,
      sensitivity_prior_sigma = sensitivity_prior_sigma,
      degradation_prior_sigma = degradation_prior_sigma,
      w_prior_sigma = w_prior_sigma,
      b_prior_sigma = b_prior_sigma,
      scale = scale
    )
  ))
}


simulate_gp <- function(num_time = 101, measurement_times = seq(1, num_time, by = 10))
{
  gp_length <- 3
  gp_variance <- 1

  time <-  1:num_time;

  regulator_profile <- generate_random_profile(time, gp_variance, gp_length, positive_transform = FALSE) %>%
    array(num_time) #Transform to 1D array



  initial_condition <-  exp(rnorm(1, -0.5,1));
  basal_transcription <-  0;#exp(rnorm(numTargets, -0.8,0.2));
  degradation <-  exp(rnorm(1, 2,1));
  sensitivity <-  exp(rnorm(1, 2,1));

  #measurement_sigma <- abs(rcauchy(1,0,0.1)) + 0.00001
  measurement_sigma <- 0.00001


  expression_true <- array(0, c(num_time));
  expression_observed = array(0, c(num_time));

  while((sensitivity) / degradation < 0.5) {
    degradation = exp(rnorm(1, 2,1));
    sensitivity = exp(rnorm(1, 2,1));
  }

  #params <- c(degradation = degradation, bias = 0, sensitivity = sensitivity, weight = 1, basal_transcription = basal_transcription, protein = approxfun(time, regulator_profile, rule=2));
  #expression_true <-  ode( y = c(x = initial_condition), times = time, func = target_ODE, parms = params, method = "ode45")[,"x"];

  expression_true <- numerical_integration(basal_transcription, degradation, initial_condition, sensitivity, 1, 0, regulator_profile, num_time)


  expression_observed <- expression_true[measurement_times] + rnorm(length(measurement_times),0,measurement_sigma)
  expression_observed[expression_observed < 0] <- 0


  data <- list(
    true = list(
      expression_full = expression_true,
      expression = expression_true[measurement_times],
      regulator_profile = regulator_profile
    ),
    observed = list(
      num_measurements = length(measurement_times),
      num_time = num_time,
      measurement_times = measurement_times,
      expression = expression_observed,
      gp_length = gp_length,
      gp_variance = gp_variance,
      measurement_sigma = measurement_sigma,
      initial_condition = initial_condition,
      degradation = degradation,
      sensitivity = sensitivity
    )
  )
}

numerical_integration <- function(basal_transcription, degradation, initial_condition, sensitivity, weight, bias, protein, num_time){

  numerical_result = numeric(num_time);
  numerical_result[1] = initial_condition
  integrationTimeStep = 1

  basal_over_degradation = basal_transcription / degradation;

  regulation_input = bias + weight * protein;
  synthesis = integrationTimeStep * ( 1 / (1 + exp(-regulation_input))) #Don't forget to change the line at the end

  residual = -0.5 * synthesis[1];
  degradationPerStep = exp(-degradation * integrationTimeStep)


  for(time in 2:num_time){
    integral = 0;
    integration_end <- time
    for(previousIndex in 1:integration_end)
    {
      if(previousIndex == 1 || previousIndex == integration_end)
      {
        h = 0.5;
      }
      else
      {
        h = 1;
      }
      degradationCoeff = exp(-degradation * (integration_end - previousIndex) * integrationTimeStep)
      integral = integral + h * synthesis[previousIndex] * degradationCoeff;
    }
    numerical_result[time] = basal_over_degradation + (initial_condition - basal_over_degradation) * exp(-degradation * (integration_end)) + sensitivity * integral
  }
  return(numerical_result)
}

