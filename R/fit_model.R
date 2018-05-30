regulated_init_fun <- function(data){
  function(chain_id = 0) {
    list(
      initial_condition = array(data$expression[1,], data$num_targets)
    )
  }
}


regulated_control <- list(adapt_delta = 0.95)

fit_regulated <- function(data, model = stan_model(file = here::here('stan', 'regulated.stan')), ...) {
  rstan::sampling(model, data = data,
                  init = regulated_init_fun(data), control = regulated_control,
                  ...)
}

fit_regulated_multi <- function(data, output.dir,
                                       model = stan_model(file = here::here('stan', 'regulated.stan')),
                                       chains = 4, ids_to_compute = 1:length(data),
                                       ...) {
  #init_per_item <- lapply(data, function(x) { lapply(1:chains, regulated_init_fun(x)) })
  init_per_item <- lapply(data, regulated_init_fun)
  sampling_multi(model = model, data = data, output.dir = output.dir,
                 init_per_item = init_per_item, control = regulated_control,
                 chains = chains, ids_to_compute = ids_to_compute,
                 ...)
}


get_waic_regulated <- function(fit, target = 1) {
  samples_log_lik <- rstan::extract(fit, "log_likelihood")$log_likelihood
  lik_data <- samples_log_lik[,,target]

  waic(lik_data)$estimates["waic","Estimate"]
}

get_waic_csynth <- function(fit) {
  samples_log_lik <- rstan::extract(fit, "log_likelihood")$log_likelihood
  lik_data <- samples_log_lik

  waic(lik_data)$estimates["waic","Estimate"]
}


get_waic_multi <- function(results, waic_fun, cores = parallel::detectCores()) {
  cl <- parallel::makeCluster(cores, useXDR = FALSE)
  on.exit(parallel::stopCluster(cl))


  extract_waic <- function(id) {
    fit <- sampling_multi_read_fit(results, id)
    waic_fun(fit)
  }

  dependencies <- c("rstan","Rcpp","genexpiStan","loo")
  .paths <- unique(c(.libPaths(), sapply(dependencies, FUN = function(d) {
    dirname(system.file(package = d))
  })))
  .paths <- .paths[.paths != ""]
  parallel::clusterExport(cl, varlist = ".paths", envir = environment())
  parallel::clusterEvalQ(cl, expr = .libPaths(.paths))
  parallel::clusterEvalQ(cl, expr =
                           suppressPackageStartupMessages(require(rstan, quietly = TRUE)))
  parallel::clusterEvalQ(cl, expr =
                           suppressPackageStartupMessages(require(loo, quietly = TRUE)))

  parallel::clusterExport(cl, varlist = c("results"), envir = environment())

  parallel::parSapplyLB(cl, X = 1:length(results), FUN = extract_waic)
}

get_waic_regulated_multi <- function(results, target = 1, cores = parallel::detectCores()) {

  waic_fun <- function(fit) {
    get_waic_regulated(fit, target)
  }

  get_waic_multi(results, waic_fun, cores)
}

get_waic_csynth_multi <- function(results, cores = parallel::detectCores()) {

  get_waic_multi(results, get_waic_csynth, cores)
}
