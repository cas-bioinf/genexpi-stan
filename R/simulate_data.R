target_ODE <- function(t, state, parameters, protein_transform = identity)
{
  with(as.list(c(state, parameters)), {
    regulatoryInput = bias + weight * protein_transform(protein(t));
    dX = basal_transcription + (sensitivity/(1 + exp(-regulatoryInput))) - degradation * x;
    list(dX)
  })
}

constant_synthesis_ODE <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX = synthesis - degradation * x;
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


simulate_multiple_targets_spline <- function(num_targets, num_time, num_knots, measurement_times, measurement_sigma_absolute, measurement_sigma_relative, regulator_measured = TRUE, integrate_ode45 = TRUE) {
  time <- 1:num_time;

  spline_degree = 3
  knots <- seq(from = 1, to = num_time, length.out = num_knots)
  spline_basis <- bs(1:num_time, knots = knots, degree = spline_degree)

  num_coeff <- dim(spline_basis)[2] - 1
  if(any(spline_basis[,num_coeff + 1] != 0)) {
    stop("Broken assumption of zero column in bs output")
  }
  spline_basis = spline_basis[,1:num_coeff]

  scale <- 5

  #Rejection sampling to get interesting profile
  n_rejections <- 0
  repeat {
    coeffs <- rnorm(num_coeff, 0, 1)
    regulator_profile <- array((spline_basis %*% coeffs)  * scale, num_time)
    if(sd(regulator_profile) > 1
       && sd(regulator_profile[1:floor(num_time / 2)]) > 1
       && sd(regulator_profile[ceiling(num_time / 2):num_time]) > 1
       && mean(diff(regulator_profile) > 1e-2) > 0.2 #Avoid monotonicity
       && mean(diff(regulator_profile) < 1e-2) > 0.2
    ) {
      break;
    }
    n_rejections <- n_rejections + 1
  }
  cat(n_rejections, " rejections regulator\n")

  min_intercept <-  -min(regulator_profile)
  if(regulator_measured) {
    intercept_prior_sigma <- 1#scale
    intercept <- min_intercept + abs(rnorm(1, 0, intercept_prior_sigma))
  } else {
    intercept <- min_intercept
    intercept_prior_sigma <- 0.01
  }




  expression_true <- array(-1, c(num_targets,num_time))
  expression_observed <- array(-1, c(num_targets,length(measurement_times)))

  sensitivity_prior_sigma <-  1
  degradation_prior_mean <- -3
  degradation_prior_sigma <- 1
  b_prior_sigma <- 2
  initial_condition_prior_sigma <- 1
  w_prior_sigma <- 2

  b_raw <- array(-Inf, num_targets)
  b <- array(-Inf, num_targets)
  w <- array(-Inf, num_targets)
  initial_condition <- array(-Inf, num_targets)
  degradation <- array(-Inf, num_targets)
  sensitivity <- array(-Inf, num_targets)


  n_rejections <- 0
  for(t in 1:num_targets) {
    #Rejection sampling to have interesting profiles
    repeat {
      sensitivity[t] <- abs(rnorm(1, 0, sensitivity_prior_sigma))

      #degradation_over_sensitivity[t] <- rlnorm(1,0,1)
      #degradation[t] <- degradation_over_sensitivity[t] * sensitivity[t];
      #degradation[t] <- abs(rnorm(1, 0, degradation_prior_sigma))
      degradation[t] <- rlnorm(1, degradation_prior_mean, degradation_prior_sigma)

      w[t] <- rnorm(1, 0, w_prior_sigma)
      b_raw[t] <- rnorm(1, 0, b_prior_sigma)
      b[t] <- b_raw[t] + intercept * w[t]
      initial_condition[t] <- abs(rnorm(1, 0,initial_condition_prior_sigma))

      if(integrate_ode45){
        params <- c(degradation = degradation[t], bias = b_raw[t], sensitivity = sensitivity[t], weight = w[t], basal_transcription = 0, protein = approxfun(time, regulator_profile, rule=2));
        expression_true[t,] <-  ode( y = c(x = initial_condition[t]), times = time, func = target_ODE, parms = params, method = "ode45")[,"x"];
      } else {
        expression_true[t,] <- numerical_integration(0, degradation[t], initial_condition[t], sensitivity[t], weight = w[t], bias = b_raw[t], regulator_profile, num_time)
      }

      if(all(expression_true[t,] > 0)
         && sd(expression_true[t,]) > 1
         #&& any(expression_true[t,] > 2 * measurement_sigma_absolute)
         && mean(diff(expression_true[t,]) > 1e-2) > 0.2 #Avoid monotonicity
         && mean(diff(expression_true[t,]) < 1e-2) > 0.2
      ) {
        break;
      }
      n_rejections <- n_rejections + 1
    }
    expression_observed[t,] <- rtruncnorm(length(measurement_times), mean = expression_true[t,measurement_times], sd =  measurement_sigma_absolute + measurement_sigma_relative * expression_true[t,measurement_times], a = 0)
  }

  cat(n_rejections, " rejections targets\n")

  regulator_expression_true <- regulator_profile + intercept
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
      num_targets = num_targets,
      regulator_measured = if (regulator_measured) { 1 } else { 0 },

      coeffs_prior_given = 0, #No coeffs prior. Sadly have to give at least size 1
      coeffs_prior_mean = array(0.0,1),
      coeffs_prior_cov = array(1, c(1,1)),

      measurement_times = measurement_times,

      num_spline_basis = num_coeff,
      spline_basis = t(spline_basis),

      expression = expression_observed,
      regulation_signs = sign(w),
      regulator_expression = if (regulator_measured) { regulator_expression_observed } else { numeric(0) },
      measurement_sigma_absolute = measurement_sigma_absolute,
      measurement_sigma_relative = measurement_sigma_relative,
      intercept_prior_sigma = intercept_prior_sigma,
      initial_condition_prior_sigma = initial_condition_prior_sigma,
      sensitivity_prior_sigma = sensitivity_prior_sigma,
      degradation_prior_mean = degradation_prior_mean,
      degradation_prior_sigma = degradation_prior_sigma,
      w_prior_sigma = w_prior_sigma,
      b_prior_sigma = b_prior_sigma,
      scale = scale
    )
  ))
}


simulate_constant_synthesis <- function(measurement_times, measurement_sigma_absolute, measurement_sigma_relative, integrate_ode45 = TRUE) {



  expression_true <- array(-1, length(measurement_times))
  expression_observed <- array(-1, length(measurement_times))

  synthesis_over_degradation_prior_mean <- 3
  synthesis_over_degradation_prior_sigma <-  3
  degradation_prior_mean <- -3
  degradation_prior_sigma <- 1
  initial_condition_prior_sigma <- 1

  n_rejections <- 0
  #Rejection sampling to have interesting profiles
  repeat {
    synthesis_over_degradation <- rnorm(1,synthesis_over_degradation_prior_mean,synthesis_over_degradation_prior_sigma)
    degradation <- rlnorm(1, degradation_prior_mean, degradation_prior_sigma)
    synthesis = synthesis_over_degradation * degradation

    initial_condition <- abs(rnorm(1, 0,initial_condition_prior_sigma))

    if(integrate_ode45){
      params <- c(degradation = degradation, synthesis = synthesis);
      expression_true <-  ode( y = c(x = initial_condition), times = measurement_times, func = constant_synthesis_ODE, parms = params, method = "ode45")[,"x"];
    } else {
      expression_true <- numerical_integration(synthesis, degradation, initial_condition, 0, 0, 0, 0, max[measurement_times])[measurement_times]
    }

    if(all(expression_true > 0))
    {
      break;
    }
    n_rejections <- n_rejections + 1
  }
  expression_observed <- rtruncnorm(length(measurement_times), mean = expression_true, sd =  measurement_sigma_absolute + measurement_sigma_relative * expression_true, a = 0)

  cat(n_rejections, " rejections targets\n")

  return(list(
    true = list (
      initial_condition = initial_condition,
      synthesis_over_degradation = synthesis_over_degradation,
      degradation = degradation,
      expression = expression_true
    ), observed = list(
      num_measurements = length(measurement_times),
      measurement_times = measurement_times,
      expression = expression_observed,
      measurement_sigma_absolute = measurement_sigma_absolute,
      measurement_sigma_relative = measurement_sigma_relative,
      initial_condition_prior_sigma = initial_condition_prior_sigma,
      synthesis_over_degradation_prior_mean = synthesis_over_degradation_prior_mean,
      synthesis_over_degradation_prior_sigma = synthesis_over_degradation_prior_sigma,
      degradation_prior_mean = degradation_prior_mean,
      degradation_prior_sigma = degradation_prior_sigma
    )
  ))
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

