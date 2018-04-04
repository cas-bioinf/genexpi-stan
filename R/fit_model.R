multiple_targets_init_fun <- function(data){
  function(chain_id = 0) {
    list(
      initial_condition = array(data$expression[1,], data$num_targets)
    )
  }
}


multiple_targets_control <- list(adapt_delta = 0.95)

fit_multiple_targets <- function(data, model = stan_model(file = 'multiple_targets_splines.stan'), ...) {
  rstan::sampling(model, data = data,
                  init = multiple_targets_init_fun(data), control = multiple_targets_control,
                  ...)
}

fit_multiple_targets_multi <- function(data, output.dir,
                                       model = stan_model(file = 'multiple_targets_splines.stan'),
                                       chains = 4,
                                       ...) {
  #init_per_item <- lapply(data, function(x) { lapply(1:chains, multiple_targets_init_fun(x)) })
  init_per_item <- lapply(data, multiple_targets_init_fun)
  sampling_multi(model = model, data = data, output.dir = output.dir,
                 init_per_item = init_per_item, control = multiple_targets_control,
                 chains = chains,
                 ...)
}


get_waic_multiple_targets <- function(fit, target = 1) {
  samples_log_lik <- rstan::extract(fit, "log_likelihood")$log_likelihood
  lik_data <- samples_log_lik[,,target]

  waic(lik_data)$waic
}

get_waic_multiple_targets_multi <- function(results, target = 1, cores = parallel::detectCores()) {
  cl <- parallel::makeCluster(cores, useXDR = FALSE)
  on.exit(parallel::stopCluster(cl))


  waic_fun <- function(id) {
    fit <- sampling_multi_read_fit(results, id)
    get_waic_multiple_targets(fit, target)
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

  parallel::clusterExport(cl, varlist = c("results", "target"), envir = environment())

  parallel::parSapplyLB(cl, X = 1:length(results), FUN = waic_fun)
}
