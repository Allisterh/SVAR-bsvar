
#' @export
fevd <- function(obj,
                 horizon,
                 sub_N = NULL,
                 print_means = TRUE,
                 controls = list(),
                 apply_restriction = TRUE,
                 R_progress_bar = TRUE) {

  if(class(obj) == "stanfit") {
    fit <- obj
    A <- extract(fit, pars = "A", apply_restriction = apply_restriction)[[1]]
    if(attributes(fit)$standata$B_inverse) {
      B <- extract(fit, pars = "B_inv", apply_restriction = apply_restriction)[[1]]
    } else {
      B <- extract(fit, pars = "B", apply_restriction = apply_restriction)[[1]]
    }
    N <- dim(A)[1]
    M <- dim(B)[2]
    if(!is.null(sub_N)) {
      if(sub_N < N) {
        sub_i <- sample.int(N, sub_N, replace = TRUE)
        if(attributes(fit)$standata$lags > 0) A <- A[sub_i,,]
        B <- B[sub_i,,]
        N <- dim(B)[1]
      }
    }

    # Check for variance/scale things here...
    if(attributes(obj)$standata$fix_moments != 2) stop("FEVD can only be computed if 'fix_moments' == 2.")

  } else {
    stop("'obj' needs to be 'stanfit' object.")
  }

  # R implementation:
  D_array <- array(NA, dim = c(M, M, horizon + 1, N))
  if(R_progress_bar) pb <- txtProgressBar(max = N, style = 3)
  for(sample_index in 1:N) {
    Bi <- B[sample_index,,]
    Ai <- A[sample_index,,]
    AAi <- stackA(t(Ai))
    for(h in 1:(horizon+1)) {
      if(h == 1) {
        zero <- Bi
        D_array[,,h,sample_index] <- zero
      } else {
        D_array[,,h,sample_index] <- (expm::`%^%`(AAi, (h-1))[1:M, 1:M] %*% zero)
      }
    }
    if(R_progress_bar) setTxtProgressBar(pb, sample_index)
  }
  if(R_progress_bar) close(pb)

  dim(D_array)
  D_array[,,1,1]
  D_array[,,2,1]
  D_array[,,3,1]

  fevd_list <- list()
  for(i in 1:M) fevd_list[[i]] <- array(NA, dim = c(horizon + 1, M, N))
  for(sample_index in 1:N) {
    fevd_matrix <- matrix(0, ncol = M, nrow = M)
    for(h in 1:(horizon + 1)) {
      fevd_matrix <- fevd_matrix + D_array[,,h,sample_index]^2
      for(i in 1:M) fevd_list[[i]][h,,sample_index] <- fevd_matrix[i,] / sum(fevd_matrix[i,])
    }
  }

  if(print_means) {
    cat("\n")
    cat("----- FEVD Posterior Means ----- \n")
    cat("\n")
    for(i in 1:M) {
      cat(paste0(colnames(attributes(obj)$original_y)[i], ": \n"))
      cat("\n")
      fevd_post_mean <- apply(fevd_list[[i]], c(1,2), mean)
      colnames(fevd_post_mean) <- paste0("Shock ", 1:M)
      print(fevd_post_mean)
      cat("\n")
    }
  }

  invisible(fevd_list)
}

collect_fevd_means <- function(fevd_object,
                               horizons = c(0, 4, 8, 12, 16, 20),
                               varnames = NULL,
                               shocknames = NULL,
                               digits = 3) {
  m <- length(fevd_object)
  if(is.null(varnames)) varnames <- paste0("Variable ", 1:m)
  if(is.null(shocknames)) shocknames <- paste0("Shock ", 1:m)
  m <- length(fevd_object)
  ret <- list()
  for(i in 1:m) {
    mat <- matrix(nrow = length(horizons), ncol = m)
    rownames(mat) <- as.character(horizons)
    colnames(mat) <- shocknames
    mat[,] <- apply(fevd_object[[i]], c(1, 2), mean)[horizons + 1,]
    ret[[i]] <- round(mat, digits)
  }
  names(ret) <- varnames
  ret
}
