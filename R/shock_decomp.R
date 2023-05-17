
#' \code{shock_decomp} computes ... for models
#' estimated with \code{bsvar}.
#'
#' @param obj An object of class `stanfit` as returned by `bsvar()`.
#' @return ...
#'
#' @export
shock_decomp <- function(obj,
                         sub_N = NULL,
                         apply_restriction = TRUE,
                         progress_bar = TRUE) {

  fit <- obj
  lags <- attributes(fit)$standata$lags
  if(lags > 0) A <- extract(fit, pars = "A", apply_restriction = apply_restriction)[[1]]
  if(attributes(fit)$standata$B_inverse) {
    B <- extract(fit, pars = "B_inv", apply_restriction = apply_restriction)[[1]]
  } else {
    B <- extract(fit, pars = "B", apply_restriction = apply_restriction)[[1]]
  }
  N <- dim(B)[1]
  M <- dim(B)[2]
  if(!is.null(sub_N)) {
    if(sub_N < N) {
      sub_i <- sample.int(N, sub_N, replace = TRUE)
      if(lags > 0) A <- A[sub_i,,]
      B <- B[sub_i,,]
      N <- dim(B)[1]
    }
  }
  yy <- attributes(fit)$standata$y
  xx <- attributes(fit)$standata$x
  if(attributes(fit)$standata$include_constant == TRUE) {
    const <- extract(fit, pars = "constant", apply_restriction = apply_restriction)[[1]]
  } else {
    const <- matrix(0, N, M)
  }
  decomp <- list()
  decomp$residuals <- array(NA, dim = c(N, nrow(yy), ncol(yy)))
  decomp$shocks <- array(NA, dim = c(N, nrow(yy), ncol(yy)))
  include_garch <- attributes(fit)$standata$include_garch
  garch_group_mat <- attributes(fit)$standata$garch_group_mat
  eta_form <- attributes(fit)$standata$garch_eta_form
  if(include_garch == TRUE) {
    fix_moments <- attributes(fit)$standata$fix_moments
    garch_c <- extract(fit, pars = "garch_c", apply_restriction = apply_restriction)[[1]]
    garch_C <- extract(fit, pars = "garch_C", apply_restriction = apply_restriction)[[1]]
    garch_D <- extract(fit, pars = "garch_D", apply_restriction = apply_restriction)[[1]]
    decomp$volatility <- array(NA, dim = c(N, nrow(yy), ncol(garch_c)))
    decomp$volatility[,1,] = 1
    decomp$shocks_raw <- array(NA, dim = c(N, nrow(yy), ncol(yy)))
  }

  # R implementation:
  if(progress_bar) pb <- txtProgressBar(0, N, style = 3)
  for(i in 1:N) {
    Bi <- B[i,,]
    if(lags > 0) Ai <- A[i,,]
    consti <- const[i,]

    # Residuals
    if(lags > 0) u <- yy - (consti + xx %*% Ai) else u <- yy - consti
    decomp$residuals[i,,] <- u

    # Structural shocks
    e <- u %*% t(Bi)
    decomp$shocks[i,,] <- e
    if(include_garch == TRUE) decomp$shocks_raw[i,,] <- e

    # Volatility
    if(include_garch == TRUE) {
      garch_ci <- garch_c[i,]
      garch_Ci <- garch_C[i,]
      garch_Di <- garch_D[i,,]
      if(length(garch_group_mat) == 0) {
        for(j in 2:nrow(yy)) {
          if(eta_form) {
            decomp$volatility[i,j,] <- sqrt( garch_ci + garch_Ci * decomp$volatility[i,j-1,]^2 + garch_Di %*% e[j-1,]^2 )
            decomp$shocks_raw[i,j,] <- e[j,] / decomp$volatility[i,j,]
          } else {
            decomp$volatility[i,j,] <- sqrt( garch_ci + garch_Ci * decomp$volatility[i,j-1,]^2 + garch_Di %*% (e[j-1,] * decomp$volatility[i,j-1,])^2 )
            decomp$shocks_raw[i,j,] <- e[j,] / decomp$volatility[i,j,]
          }
        }
      } else {
        for(j in 2:nrow(yy)) {
          mean_last_e_squared <- rep(NA, ncol(garch_group_mat))
          if(eta_form) {
            for(h in 1:ncol(garch_group_mat)) mean_last_e_squared[h] <- mean(e[j-1, which(garch_group_mat[,h] == 1)]^2)
            decomp$volatility[i,j,] <- sqrt( garch_ci + garch_Ci * decomp$volatility[i,j-1,]^2 + garch_Di %*% mean_last_e_squared)
            for(h in 1:ncol(garch_group_mat)) decomp$shocks_raw[i,j,which(garch_group_mat[,h] == 1)] <- e[j,which(garch_group_mat[,h] == 1)] / decomp$volatility[i,j,h]
          } else {
            for(h in 1:ncol(garch_group_mat)) mean_last_e_squared[h] <- mean((e[j-1, which(garch_group_mat[,h] == 1)] * decomp$volatility[i,j-1,h])^2)
            decomp$volatility[i,j,] <- sqrt( garch_ci + garch_Ci * decomp$volatility[i,j-1,]^2 + garch_Di %*% mean_last_e_squared)
            for(h in 1:ncol(garch_group_mat)) decomp$shocks_raw[i,j,which(garch_group_mat[,h] == 1)] <- e[j,which(garch_group_mat[,h] == 1)] / decomp$volatility[i,j,h]
          }
        }
      }

    }
    if(progress_bar) setTxtProgressBar(pb, i)
  }
  if(progress_bar) close(pb)

  # Rcpp implementation:
  # ...

  if(include_garch == TRUE) {
    decomp$residuals <- decomp$residuals[,-1,]
    decomp$shocks <- decomp$shocks[,-1,]
    decomp$shocks_raw <- decomp$shocks_raw[,-1,]
    decomp$volatility <- decomp$volatility[,-1,]
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
plot_shock_volatility <- function(decomp, first_period = 0, freq = 1, y_labels = NULL, probs = NULL, fit = NULL, controls = list()) {
  if(!is.null(fit)) {
    y <- attributes(fit)$original_y
    if("ts" %in% class(y)) {
      lags <- attributes(fit)$standata$lags
      freq <- attributes(y)$tsp[3]
      first_period <- (attributes(y)$tsp[1] + 1 / freq) + (lags / freq)
    }
  }
  if(is.null(controls$col)) controls$col <- "tomato"
  if(is.null(probs)) probs <- c(0.05, 0.16, 0.84, 0.95)
  probs <- c(probs[1] - 0.01, probs, probs[length(probs)] + 0.01)
  if(min(probs) < 0 || max(probs) > 1) stop("'probs' must have values in between 0.01 and 0.99.")
  med_volatility <- ts(apply(decomp$volatility, c(2,3), median), start = first_period, frequency = freq)
  y_labels <- paste0("Volatility of shock ", 1:ncol(med_volatility))
  if(!is.null(controls$mfrow)) {
    rows <- controls$mfrow[1]
    cols <- controls$mfrow[2]
  } else {
    cols <- rows <- ceiling(sqrt(length(y_labels)))
    while(cols * rows - cols >= length(y_labels)) rows <- rows - 1
  }
  par(mfrow = c(rows, cols))
  if(!is.null(controls$mfcol)) {
    rows <- controls$mfcol[1]
    cols <- controls$mfcol[2]
    par(mfcol = c(rows, cols))
  }
  if(is.null(controls$mar)) controls$mar <- c(3, 4.5, 1, 1)
  par(mar = controls$mar)
  for(i in 1:ncol(med_volatility)) {
    sub_volatility <- apply(decomp$volatility[,,i], 2, function(x) quantile(x, probs = probs))
    plot(med_volatility[,i], main = "", ylab = y_labels[i], lwd = 1, ylim = c(0, max(sub_volatility) + 0.1))
    grid()
    fanplot::fan(data = sub_volatility, data.type = "values", probs = probs,
                 start = first_period, frequency = freq,
                 fan.col = colorRampPalette(c(controls$col, "white")),
                 rlab = NULL, ln = NULL)
    lines(med_volatility[,i], lwd = 1, lty = 2)
  }
  par(mfrow = c(1,1))
}

plot_shock_volatility <- function(decomp, first_period = 0, freq = 1, y_labels = NULL, probs = NULL, fit = NULL, controls = list()) {
  if(!is.null(fit)) {
    y <- attributes(fit)$original_y
    if("ts" %in% class(y)) {
      lags <- attributes(fit)$standata$lags
      freq <- attributes(y)$tsp[3]
      first_period <- (attributes(y)$tsp[1] + 1 / freq) + (lags / freq)
    }
  }
  if(is.null(controls$col)) controls$col <- "tomato"
  if(is.null(probs)) probs <- c(0.05, 0.16, 0.84, 0.95)
  probs <- c(probs[1] - 0.01, probs, probs[length(probs)] + 0.01)
  if(min(probs) < 0 || max(probs) > 1) stop("'probs' must have values in between 0.01 and 0.99.")
  med_volatility <- ts(apply(decomp$volatility, c(2,3), median), start = first_period, frequency = freq)
  y_labels <- paste0("Volatility of shock ", 1:ncol(med_volatility))
  if(!is.null(controls$mfrow)) {
    rows <- controls$mfrow[1]
    cols <- controls$mfrow[2]
  } else {
    cols <- rows <- ceiling(sqrt(length(y_labels)))
    while(cols * rows - cols >= length(y_labels)) rows <- rows - 1
  }
  par(mfrow = c(rows, cols))
  if(!is.null(controls$mfcol)) {
    rows <- controls$mfcol[1]
    cols <- controls$mfcol[2]
    par(mfcol = c(rows, cols))
  }
  if(is.null(controls$mar)) controls$mar <- c(3, 4.5, 1, 1)
  par(mar = controls$mar)
  for(i in 1:ncol(med_volatility)) {
    sub_volatility <- apply(decomp$volatility[,,i], 2, function(x) quantile(x, probs = probs))
    plot(med_volatility[,i], main = "", ylab = y_labels[i], lwd = 1, ylim = c(0, max(sub_volatility) + 0.1))
    grid()
    fanplot::fan(data = sub_volatility, data.type = "values", probs = probs,
                 start = first_period, frequency = freq,
                 fan.col = colorRampPalette(c(controls$col, "white")),
                 rlab = NULL, ln = NULL)
    lines(med_volatility[,i], lwd = 1, lty = 2)
  }
  par(mfrow = c(1,1))
}






