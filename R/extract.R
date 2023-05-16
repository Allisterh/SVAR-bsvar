
#' @export
extract <- function(object, pars, permuted = TRUE, inc_warmup = FALSE,
                    include = TRUE, apply_restriction = TRUE) {
  out <- rstan::extract(object, pars)
  res <- attributes(object)$restriction
  if(apply_restriction && !is.null(res)) {
    out <- lapply(out, function(x) if(length(dim(x)) == 3) x[-res,,] else if(length(dim(x)) == 2) x[-res,] else x[-res])
  }
  out
}
