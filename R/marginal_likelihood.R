
#' @export
marginal_likelihood <- function(fit, ...) {
  if(rstan::stan_version() < "2.26") {
    cat("It seems that your rstan installation is not up to date ('marginal_likelihood' requires rstan 2.26.0 or higher).\n")
    cat("Would you like me to try and install the latest version of rstan for you?")
    ans <- readline(prompt = "(y/n): ")
    if(ans == "y") {
      try(
        install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
      )
    }
  }
  if(rstan::stan_version() < "2.26") {
    stop("Ugh, 'rstan version >= 2.26' still not found. For more instructions on how to get your hands on the latest version of rstan, see: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started")
  }
  is_cmdstanr_model <- length(fit@stanmodel@model_code) == 0
  if(is_cmdstanr_model) {
    standata <- attributes(fit)$standata
    standata$parallel_likelihood <- 0
    if(!exists("mod_rstan")) {
      cat("First time calling 'marginal_likelihood' in a session an 'rstan' model is built,\n(frustrating, I know) which might take a moment or two...\n")
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





