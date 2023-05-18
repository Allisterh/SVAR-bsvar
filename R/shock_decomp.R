
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

#' @export
plot_volatility <- function(decomp, y_labels = NULL, probs = NULL, control = list()) {
  if(is.null(control$col)) control$col <- "tomato"
  if(is.null(probs)) probs <- c(0.05, 0.16, 0.84, 0.95)
  probs <- c(probs[1] - 0.01, probs, probs[length(probs)] + 0.01)
  if(min(probs) < 0 || max(probs) > 1) stop("'probs' must have values in between 0.01 and 0.99.")
  if(!is.null(decomp$tsp_)) {
    first_period <- decomp$tsp_[1]; freq <- decomp$tsp_[3]
  } else {
    first_period <- 1; freq <- 1
  }
  med_volatility <- ts(apply(decomp$volatility, c(2,3), median), start = first_period, frequency = freq)
  y_labels <- paste0("Volatility of Shock ", 1:ncol(med_volatility))
  if(!is.null(control$mfrow)) {
    rows <- control$mfrow[1]
    cols <- control$mfrow[2]
  } else {
    cols <- rows <- ceiling(sqrt(length(y_labels)))
    while(cols * rows - cols >= length(y_labels)) rows <- rows - 1
  }
  par(mfrow = c(rows, cols))
  if(!is.null(control$mfcol)) {
    rows <- control$mfcol[1]
    cols <- control$mfcol[2]
    par(mfcol = c(rows, cols))
  }
  if(is.null(control$mar)) control$mar <- c(3, 4.5, 1, 1)
  par(mar = control$mar)
  for(i in 1:ncol(med_volatility)) {
    sub_volatility <- apply(decomp$volatility[,,i], 2, function(x) quantile(x, probs = probs))
    plot(med_volatility[,i], main = "", ylab = y_labels[i], lwd = 1, ylim = c(0, max(sub_volatility) + 0.1))
    grid()
    fanplot::fan(data = sub_volatility, data.type = "values", probs = probs,
                 start = first_period, frequency = freq,
                 fan.col = colorRampPalette(c(control$col, "white")),
                 rlab = NULL, ln = NULL)
    lines(med_volatility[,i], lwd = 1, lty = 2)
  }
  par(mfrow = c(1,1))
}





