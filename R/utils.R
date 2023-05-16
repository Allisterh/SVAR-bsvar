
build_xy <- function(y, lags, constant = FALSE) {
  if(lags == 0) {
    yy <- y
    xx <- as.array(matrix(0, nrow = nrow(y), ncol = 0))
    if(constant == TRUE) xx <- matrix(rep(1, nrow(yy)), ncol = 1)
  } else {
    n <- ncol(y)
    for(i in 1:lags) {
      if(i == 1) {
        temp <- rbind(rep(NA, n), y[-(nrow(y)),])
        xx <- temp
      } else {
        temp <- rbind(rep(NA, n), temp[-(nrow(y)),])
        xx <- cbind(xx, temp)
      }
    }
    if(constant == TRUE) xx <- cbind(rep(1, nrow(xx)), xx)
    xx <- xx[-c(1:lags),]
    yy <- y[-c(1:lags),]
  }
  ret <- list("xx" = xx, "yy" = yy)
  ret
}

stackA <- function(A) {
  m <- nrow(A)
  lags <- ncol(A)/m
  eye <- diag(m*lags-m)
  A <- rbind(A, cbind(eye, matrix(0, ncol = m, nrow= nrow(eye))))
  A
}

minnesota_sds <- function(M, lags, hypers, reduced_scales = rep(1, M)) {
  sd_mat <- matrix(NA, nrow = M*lags, ncol = M)
  for(h in 1:lags) {
    for(i in 1:M) {
      for(j in 1:M) {
        if (i == j) {
          if (h == 1) {
            sd_mat[(h-1)*M + i, j] <- hypers[1]
          } else {
            sd_mat[(h-1)*M + i, j] <- hypers[1] / (h^hypers[2])
          }
        } else {
          sd_mat[(h-1)*M + i, j] <- (hypers[1] * hypers[3] * reduced_scales[j]) / ((h^hypers[2]) * reduced_scales[i])
        }
      }
    }
  }
  sd_mat
}

missing_data_prior_template <- function(y, full_inf = FALSE) {
  missing_data_prior <- list()
  missing_data_prior$mu <- missing_data_prior$sigma <- missing_data_prior$df <- y
  missing_data_prior$mu[!is.na(y)] <- missing_data_prior$sigma[!is.na(y)] <- missing_data_prior$df[!is.na(y)] <- Inf
  if(full_inf) missing_data_prior$mu[is.na(y)] <- missing_data_prior$sigma[is.na(y)] <- missing_data_prior$df[is.na(y)] <- Inf
  missing_data_prior
}


