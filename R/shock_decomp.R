
#' \code{shock_decomp} computes ... for models
#' estimated with \code{bsvar}.
#'
#' @param obj An object of class `stanfit` as returned by `bsvar()`.
#' @return ...
#'
#' @export
shock_decomp <- function(obj,
                         sub_N = NULL,
                         apply_restriction = TRUE) {

  fit <- obj
  decomp <- list()
  decomp$residuals <- extract(fit, pars = "residuals", apply_restriction = apply_restriction)[[1]]
  decomp$shocks <- extract(fit, pars = "shocks", apply_restriction = apply_restriction)[[1]]
  decomp$volatility <- extract(fit, pars = "volatility", apply_restriction = apply_restriction)[[1]]
  include_garch_groups <- attributes(fit)$standata$include_garch_groups
  if(include_garch_groups) {
    garch_group_mat <- attributes(fit)$standata$garch_group_mat
    vol_strecth <- array(NA, dim = dim(decomp$shocks))
    for(j in 1:dim(vol_strecth)[3]) vol_strecth[,,j] <- decomp$volatility[,,which(garch_group_mat[j,] == 1)]
    decomp$volatility <- vol_strecth
  }
  decomp$shocks_raw <- decomp$shocks * decomp$volatility

  include_garch <- attributes(fit)$standata$include_garch
  if(include_garch == TRUE) {
    decomp$residuals <- decomp$residuals[,-1,]
    decomp$shocks <- decomp$shocks[,-1,]
    decomp$volatility <- decomp$volatility[,-1,]
    decomp$shocks_raw <- decomp$shocks_raw[,-1,]
  }
  if("ts" %in% attributes(attributes(fit)$standata$y_raw)$class) {
    attributes(decomp)$tsp_ <- attributes(attributes(fit)$standata$y_raw)$tsp
    if(include_garch) attributes(decomp)$tsp_[1] <- attributes(decomp)$tsp_[1] + (1 / attributes(decomp)$tsp_[3])
    lags <- attributes(fit)$standata$lags
    if(lags > 0) attributes(decomp)$tsp_[1] <- attributes(decomp)$tsp_[1] + (lags / attributes(decomp)$tsp_[3])
  }
  invisible(decomp)
}

shock_decomp_pick <- function(decomp_obj,
                              show_date,
                              start_date,
                              freq) {

  row_indices <- ts(1:dim(decomp_obj$shocks)[2],
                    start = start_date,
                    frequency = freq)
  row_index <- window(row_indices,
                      start = show_date,
                      end = show_date)

  E <- decomp_obj$shocks[,row_index,]
  U <- decomp_obj$residuals[,row_index,]
  ret <- list("E" = E,
              "U" = U)
  ret
}

narrative_sign_probs <- function(decomp_obj,
                                 start_date,
                                 freq,
                                 dates,
                                 signs) {
  N <- dim(decomp_obj$shocks)[1]
  m <- dim(decomp_obj$shocks)[3]
  probs <- matrix(NA, nrow = length(dates), ncol = m)
  total_probs <- rep(NA, m)
  for(i in 1:m) {
    for(j in 1:length(dates)) {
      picked <- shock_decomp_pick(decomp_obj, show_date = dates[[j]], start_date = start_date, freq = freq)
      if(signs[j] == 1) agree_dummy <- picked$E[,i] > 0
      if(signs[j] == -1) agree_dummy <- picked$E[,i] < 0
      if(j == 1) {
        agree_dummies <- agree_dummy
      } else {
        agree_dummies <- cbind(agree_dummies, agree_dummy)
      }
    }
    probs[,i] <- apply(agree_dummies, 2, mean)
    total_probs[i] <- mean(apply(agree_dummies, 1, function(x) ifelse(sum(x) == length(x), TRUE, FALSE)))
  }
  ret <- rbind(probs, total_probs)
  rownames(ret) <- c(1:length(dates), "Total")
  ret
}







