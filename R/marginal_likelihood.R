
#' @export
marginal_likelihood <- function(fit, ...) {
  if(rstan::stan_version() < 2.26) stop("Unfortunately, 'rstan::stan_version() >= 2.26' required for 'marginal_likelihood' to work.\nSee: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started")
  is_cmdstanr_model <- length(fit@stanmodel@model_code) == 0
  if(is_cmdstanr_model) {
    standata <- attributes(fit)$standata
    standata$parallel_likelihood <- 0
    if(!exists("mod_rstan")) {
      cat("First time calling 'marginal_likelihood' in a session (frustrating, I know)\nan 'rstan' model is built, which might take a moment or two...\n")
      mod_rstan <<- rstan::stan_model(system.file("bsvar.stan", package = "bsvar"))
    }
    cat("You can safely ignore the possible 'the number of chains ... sampling not done' message below.\n")
    mod_rstan_empty_sample <- rstan::sampling(mod_rstan, chains = 0, data = standata)
    out <- bridge_sampler_bsvar(samples = fit, stanfit_model = mod_rstan_empty_sample, res = res, ...)
  } else {
    out <- bridge_sampler_bsvar(stanfit_model = fit, res = res, ...)
  }
  out
}





