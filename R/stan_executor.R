sampling_multi_filename <- function(output.dir, id, chain) {
  paste0(output.dir, "/stan_out_",id,"_c",chain,".csv")
}

sampling_multi <- function(model, data, output.dir, chains = 4, cores = parallel::detectCores(), ...) {
  cl <- parallel::makeCluster(cores, useXDR = FALSE)
  on.exit(parallel::stopCluster(cl))

#  objects <- ls()
  fit_fun <- function(i) {
    id <- floor( (i - 1) / chains) + 1
    chain_id <- ((i - 1) %% chains) + 1
    .dotlist$chain_id <- chain_id

    out_file <- sampling_multi_filename(output.dir, id, chain_id)
    if(is.list(.dotlist$init)) .dotlist$init <- .dotlist$init[chain_id]
    .dotlist$sample_file <- out_file
    .dotlist$data <- data[[id]]
    out <- do.call(rstan::sampling, args = .dotlist)
    return(out_file)
  }

#  .dotlist <- c(sapply(objects, simplify = FALSE, FUN = get,
#                       envir = environment()), list(...))
  .dotlist <- c(list(model, chains = 1L, cores = 0L), list(...))


  dependencies <- c("rstan", "Rcpp", "genexpiStan")
  .paths <- unique(c(.libPaths(), sapply(dependencies, FUN = function(d) {
    dirname(system.file(package = d))
  })))
  .paths <- .paths[.paths != ""]
  parallel::clusterExport(cl, varlist = ".paths", envir = environment())
  parallel::clusterEvalQ(cl, expr = .libPaths(.paths))
  parallel::clusterEvalQ(cl, expr =
                           suppressPackageStartupMessages(require(rstan, quietly = TRUE)))

  parallel::clusterExport(cl, varlist = ".dotlist", envir = environment())

  parallel::clusterExport(cl, varlist = c("data", "output.dir"), envir = environment())
  items <- 1:(chains * length(data))
  results_flat <-  parallel::parSapplyLB(cl, X = items, FUN = fit_fun)

  ids <- floor( (items - 1) / chains) + 1
  results <- list()
  for(i in 1:length(data)) {
    results[[i]] <- results_flat[ids == i]
  }
  results
}