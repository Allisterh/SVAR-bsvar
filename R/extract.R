
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
