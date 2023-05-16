
# Adapted from: bridgesampling::bridge_sampler
# Hacky fix for simplex inaccuracies
bridge_sampler_bsvar <- function(samples = NULL, stanfit_model = samples,
                                 repetitions = 1, method = "normal", cores = 1,
                                 use_neff = TRUE, maxiter = 1000, silent = FALSE,
                                 verbose = FALSE,
                                 ...) {
  # cores > 1 only for unix:
  if (!(.Platform$OS.type == "unix") & (cores != 1)) {
    warning("cores > 1 only possible on Unix/MacOs. Uses 'core = 1' instead.", call. = FALSE)
    cores <- 1L
  }

  # get simplex parameter names
  model_text <- strsplit(stanfit_model@stanmodel@model_code, "[;]")[[1]]
  model_text <- model_text[grep("simplex[", model_text, fixed = TRUE)]
  simplex_par_names <- unname(sapply(model_text, function(x) rev(strsplit(x, " ", fixed = TRUE)[[1]])[1]))

  # convert samples into matrix and fix inaccurate simplex values
  if (!requireNamespace("rstan")) stop("package rstan required")
  if(is.null(samples)) samples <- stanfit_model
  ex <- rstan::extract(samples, permuted = FALSE)
  skeleton <- bridgesampling:::.create_skeleton(samples@sim$pars_oi,
                                                samples@par_dims[samples@sim$pars_oi])
  upars <- apply(ex, 1:2, FUN = function(theta) {
    x <- bridgesampling:::.rstan_relist(theta, skeleton)
    for(i in 1:length(x)) if(names(x)[i] %in% simplex_par_names) x[[i]] <- t(apply(x[[i]], 1, function(row) row / sum(row)))
    rstan::unconstrain_pars(stanfit_model, x)
  })

  if (length(dim(upars)) == 2) { # for one parameter models
    dim(upars) <- c(1, dim(upars))
  }

  nr <- dim(upars)[2]
  samples4fit_index <- seq_len(nr) %in% seq_len(round(nr/2)) # split samples in two parts
  samples_4_fit <- apply(upars[,samples4fit_index,,drop=FALSE], 1, rbind)

  samples_4_iter_stan <- upars[,!samples4fit_index,,drop=FALSE]
  samples_4_iter_tmp <- vector("list", dim(upars)[3])
  for (i in seq_along(samples_4_iter_tmp)) {
    samples_4_iter_tmp[[i]] <- coda::as.mcmc(t(samples_4_iter_stan[,,i]))
  }
  samples_4_iter_tmp <- coda::as.mcmc.list(samples_4_iter_tmp)

  if (use_neff) {
    neff <- tryCatch(median(coda::effectiveSize(samples_4_iter_tmp)), error = function(e) {
      warning("effective sample size cannot be calculated, has been replaced by number of samples.", call. = FALSE)
      return(NULL)
    })
  } else {
    neff <- NULL
  }

  samples_4_iter <- apply(samples_4_iter_stan, 1, rbind)

  parameters <- paste0("x", (seq_len(dim(upars)[1])))

  transTypes <- rep("unbounded", length(parameters))
  names(transTypes) <- parameters

  # prepare lb and ub
  lb <- rep(-Inf, length(parameters))
  ub <- rep(Inf, length(parameters))
  names(lb) <- names(ub) <- parameters

  colnames(samples_4_iter) <- paste0("trans_", parameters)
  colnames(samples_4_fit) <- paste0("trans_", parameters)

  # run bridge sampling (slightly modified)
  fun_to_call <- getFromNamespace(paste0(".bridge.sampler.", method), "bridgesampling")
  if (cores == 1) {
    bridge_output <- do.call(what = fun_to_call,
                             args = list(samples_4_fit = samples_4_fit,
                                         samples_4_iter = samples_4_iter,
                                         neff = neff,
                                         log_posterior = bridgesampling:::.stan_log_posterior,
                                         data = list(stanfit = stanfit_model),
                                         lb = lb, ub = ub,
                                         param_types = rep("real", ncol(samples_4_fit)),
                                         transTypes = transTypes,
                                         repetitions = repetitions, cores = cores,
                                         packages = "rstan", maxiter = maxiter, silent = silent,
                                         verbose = verbose,
                                         r0 = 0.5, tol1 = 1e-10, tol2 = 1e-4))
  } else {
    bridge_output <- do.call(what = fun_to_call,
                             args = list(samples_4_fit = samples_4_fit,
                                         samples_4_iter = samples_4_iter,
                                         neff = neff,
                                         log_posterior = bridgesampling:::.stan_log_posterior,
                                         data = list(stanfit = stanfit_model),
                                         lb = lb, ub = ub,
                                         param_types = rep("real", ncol(samples_4_fit)),
                                         transTypes = transTypes,
                                         repetitions = repetitions, varlist = "stanfit",
                                         envir = sys.frame(sys.nframe()),
                                         cores = cores, packages = "rstan", maxiter = maxiter,
                                         silent = silent, verbose = verbose,
                                         r0 = 0.5, tol1 = 1e-10, tol2 = 1e-4))
  }

  return(bridge_output)

}

