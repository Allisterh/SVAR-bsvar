
#' Print a concise summary.
#'
#' Wrapper for \code{\link[rstan]{print}}, that prints by default a summary
#' of only some of the parameters from the \code{stanfit} object outputted
#' by \code{\link{bsvar}}.
#'
#' @export
prnt <- function(x, pars = NULL,
                 probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                 digits_summary = 2, include = TRUE, ...) {
  if(is.null(pars)) {
    pars_bool <- fit@sim$pars_oi %in% c("lambda", "log_p", "log_q", "constant", "hyper_shrinkage",
                                        "hyper_lags", "hyper_ownlags", "hyper_soc", "hyper_dio",
                                        "hyper_B", "garch_param", "B", "A", "lp__")
    pars <- fit@sim$pars_oi[which(pars_bool)]
  }
  print(x, pars, probs, digits_summary, include, ...)
}

