

#' Extract samples from a fitted model
#'
#' Wrapper for \code{\link[rstan]{extract}}, but with an additional argument
#' (\code{apply_restriction}).
#'
#' @param object An instance of class \code{stanfit} as outputted by \code{\link{bsvar}}.
#' @param pars A character vector of parameter names.If not specified, all parameters
#' and other quantities are used.
#' @param apply_restriction Logical. Defaults to \code{TRUE} in which case the posterior
#' sample is restricted according to a binary vector given by \code{attributes(object)$restriction},
#' if exists (typically it does not), such that ones stand for posterior draws to be dropped from
#' the restricted sample. This can be used for imposing sign restrictions a posteriori.
#' @seealso \code{\link[rstan]{extract}}
#' @export
extract <- function(object, pars, apply_restriction = TRUE) {
  out <- rstan::extract(object, pars)
  res <- attributes(object)$restriction
  drop <- which(res == 1)
  if(apply_restriction && !is.null(res)) {
    out <- lapply(out, function(x) if(length(dim(x)) == 3) x[-drop,,] else if(length(dim(x)) == 2) x[-drop,] else x[-drop])
  }
  out
}

#' Print a concise summary
#'
#' Wrapper for \code{\link[rstan]{print.stanfit}}, that prints by default a summary
#' of only the most important parameters in the \code{stanfit} object outputted
#' by \code{\link{bsvar}}. By default \code{use_cache = FALSE}, which makes \code{prnt} typically
#' much faster than \code{\link[rstan]{print.stanfit}} if there is no cached data (right after estimation),
#' but much slower if there is cached data (if \code{\link[rstan]{print.stanfit}} has already been called).
#' @param x An instance of class \code{stanfit} as outputted by \code{\link{bsvar}}.
#' @param pars A character vector of parameter names. Defaults to a predefined set of the
#' most important parameters. If \code{include = FALSE}, then the specified parameters are excluded
#' from the printed summary.
#' @param probs A numeric vector of quantiles of interest.
#' @param digits_summary The number of significant digits to use when printing the summary, defaulting to 2.
#' @param include Logical. Indicates whether to include or exclude the parameters named by the pars argument.
#' @param ... Additional arguments passed to the summary method for stanfit objects.
#' @seealso \code{\link[rstan]{print.stanfit}}
#' @export
prnt <- function(x, pars = NULL,
                 probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
                 digits_summary = 2, include = TRUE,
                 use_cache = FALSE, ...) {
  if(is.null(pars)) {
    pars_bool <- x@sim$pars_oi %in% c("lambda", "log_p", "log_q", "constant", "hyper_shrinkage",
                                        "hyper_lags", "hyper_ownlags", "hyper_soc", "hyper_dio",
                                        "hyper_B", "garch_param", "B", "A", "lp__")
    pars <- x@sim$pars_oi[which(pars_bool)]
  }
  print(x, pars, probs, digits_summary, include, ...)
}

#' \loadmathjax{}
#' \mjseqn{\hat{R}} convergence diagnostic for Markov Chains
#'
#' A convenience function that extracts \mjseqn{\hat{R}} convergence diagnostic values
#' from the summary of the \code{stanfit} object outputted by \code{\link{bsvar}}
#' and returns them as a named vector.
#'
#' @param object An instance of class \code{stanfit} as outputted by \code{\link{bsvar}}.
#' @param pars A character vector of parameter names. Defaults to a predefined set of the
#' most important parameters.
#' @seealso \code{\link[rstan]{Rhat}}
#' @export
rhat <- function(object, pars = NULL) {
  if(is.null(pars)) {
    all_pars <- object@sim$pars_oi
    pars_bool <- all_pars %in% c("lambda", "log_p", "log_q", "constant", "hyper_shrinkage",
                                 "hyper_lags", "hyper_ownlags", "hyper_soc", "hyper_dio",
                                 "hyper_B", "garch_param", "B", "A", "lp__")
    pars <- all_pars[which(pars_bool)]
  }
  s <- rstan::summary(object, pars = pars, use_cache = FALSE)
  s$summary[,"Rhat"]
}

