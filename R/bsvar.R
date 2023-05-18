
#' @export
bsvar <- function(y,
                  lags = 1,
                  include_constant = TRUE,
                  dist = "st",
                  dist_control = list("fix_moments" = 2,
                                      "p_q_mins" = NULL),
                  B_inverse = TRUE,
                  B_pos_diag = TRUE,
                  B_restrictions = NULL,
                  A_restrictions = NULL,
                  include_minnesota = TRUE,
                  minnesota_control = list("minnesota_means" = 1,
                                           "hyper_free" = c(1, 0, 0),
                                           "hyper_fixed_values" = c(0.2, 1, 0.5),
                                           "hyper_shrinkage_prior" = c(0, 1, 0.001),
                                           "soc" = FALSE, "dio" = FALSE),
                  include_garch = FALSE,
                  garch_control = list("garch_prior" = 10/12,
                                       "garch_dep" = FALSE,
                                       "garch_shrinkage" = 1,
                                       "garch_groups" = NULL,
                                       "garch_eta_form" = FALSE),
                  vol_breaks = NULL,
                  other_priors = list("B_prior" = c(2^2, NA),
                                      "prior_elasticity" = NULL,
                                      "constant_prior" = c(0, Inf),
                                      "lambda_prior" = c(4, 4),
                                      "p_prior" = c(1, 1),
                                      "q_prior" = c(2, 1),
                                      "p_q_prior_shift" = NULL,
                                      "hyper_B_prior" = c(0, 1),
                                      "hyper_B_soft_prior" = c(4, 4),
                                      "vol_breaks_prior" = 2),
                  chains = 4,
                  parallel_chains = ifelse(chains <= 4, chains, 4),
                  threads_per_chain = 1,
                  max_treedepth = 15,
                  control = list(),
                  ...) {

  if(cmdstanr::cmdstan_version() >= "2.32.1" || cmdstanr::cmdstan_version() < "2.29.0") {
    installed_cmdstans <- list.files(cmdstanr::cmdstan_default_install_path())
    suitable_cmdstans <- paste0("cmdstan-2.", c("32.0", "31.0", "30.1", "30.0", "29.2", "29.1", "29.0"))
    chosen_cmdstan <- NULL
    for(version in suitable_cmdstans) {
      if(version %in% installed_cmdstans) {
        chosen_cmdstan <- version
        break
      }
    }
    if(is.null(chosen_cmdstan)) {
      cat("A suitable cmdstan installation was not found (for now >=2.29.0 but < 2.32.1 required).\n")
      cat("Would you like me to try and download one for you?\n\n")
      cat("(This might fail if you are using Rstudio! If this indeed fails, then you should exit Rstudio, open R without Rstudio, and run the same command there. ")
      cat("Given you succeed there, you may restart Rstudio and continue there. Also, this may fill your console with all kinds of red and black text, but don't panic! It's all very normal!)\n")
      ans <- readline(prompt = "(y/n): ")
      if(ans == "y") {
        try({
          cmdstanr::install_cmdstan(version = "2.29.2") # bsvar has been most thoroughly tested with 2.29.2
          cmdstanr::set_cmdstan_path(file.path(cmdstanr::cmdstan_default_install_path(), chosen_cmdstan))
        })
      }
    } else {
      cmdstanr::set_cmdstan_path(file.path(cmdstanr::cmdstan_default_install_path(), chosen_cmdstan))
    }
    if(cmdstanr::cmdstan_version() >= "2.32.1" || cmdstanr::cmdstan_version() < "2.29.0") {
      stop("Unfortunately I need to stop you right here, as I cannot find a suitable cmdstan installation, or perhaps I have failed you in some other way. In any case, see https://github.com/jetroant/bsvar for more instructions.")
    }
  }
  stan_file <- system.file("bsvar.stan", package = "bsvar")
  mod <- cmdstanr::cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

  con <- list("transform" = TRUE,
              "no_init_ols" = NULL,
              "init_optim" = TRUE,
              "init_hyperpar_optim" = FALSE,
              "optim_iter" = 100000,
              "init_fun" = NULL,
              "metric" = "diag_e",
              "init_inv_metric" = "auto",
              "skip_sampling" = FALSE,
              "verbose" = TRUE,
              "recompile" = FALSE)
  con[(names(con) <- names(control))] <- control
  list2env(con, envir = environment())
  include_vol_breaks <- !is.null(vol_breaks)
  if(include_garch + include_vol_breaks == 2) stop("'include_garch = TRUE' only allowed when 'vol_breaks' is NULL.")
  if(include_garch == 1) {
   if(dist_control$fix_moments < 2) stop("'include_garch = TRUE' only allowed if 'dist_control$fix_moments = 2'.")
   if(!is.null(garch_control$garch_groups) && !is.null(garch_control$garch_dep)) {
     if(garch_control$garch_dep) stop("Either 'garch_control$garch_dep' needs to be FALSE or 'garch_control$garch_groups' needs to be NULL.")
   }
  }
  if(is.null(no_init_ols)) no_init_ols <- include_minnesota
  N <- nrow(y) - lags
  if(N < 0) N <- 0
  M <- ncol(y)

  #######################################
  ### Interpolation of quarterly data ###
  #######################################

  if(!is.null(attributes(y)$quarterly)) {
    if(nrow(y) %% 3 != 0 || sum(is.na(y[,attributes(y)$quarterly])) > 0) {
      stop("For now quarterly data to be interpolated needs to be balanced, ",
           "\ni.e no missing values and number of rows divisible by three.")
    }
    quarterly_binary <- as.numeric(1:M %in% attributes(y)$quarterly)
    number_of_quarterly <- sum(quarterly_binary)
    init_monthly <- matrix(NA, ncol = 3, nrow = nrow(y) * number_of_quarterly / 3)
    quarterly_sums <- rep(NA, nrow(y) * number_of_quarterly / 3)
    yq <- y[,attributes(y)$quarterly]
    if(is.null(nrow(yq))) yq <- matrix(yq, ncol = 1)
    for(i in 1:length(quarterly_sums)) {
      sub <- yq[(i * 3 - 2):(i * 3)]
      quarterly_sums[i] <- sum(exp(sub)) # Log-levels assumed
      init_monthly[i,] <- sub / sum(sub)
    }
  } else {
    number_of_quarterly <- 0
    quarterly_binary <- rep(0, M)
    init_monthly <- matrix(0, ncol = 3, nrow = 0)
    quarterly_sums <- rep(0, 0)
  }

  ##################################
  ### Distribution of the shocks ###
  ##################################

  if(dist == "sgt") {
    dist_control$sgt_fixed_values = c(0, log(2), log(5))
    dist_control$sgt_free <- c(1, 1, 1)
  } else if(dist == "st") {
    dist_control$sgt_fixed_values = c(0, log(2), log(5))
    dist_control$sgt_free <- c(1, 0, 1)
  } else if(dist == "t") {
    dist_control$sgt_fixed_values = c(0, log(2), log(5))
    dist_control$sgt_free <- c(0, 0, 1)
  } else if(dist == "snormal") {
    dist_control$sgt_fixed_values = c(0, log(2), Inf)
    dist_control$sgt_free <- c(1, 0, 0)
  } else if(dist == "normal") {
    dist_control$sgt_fixed_values = c(0, log(2), Inf)
    dist_control$sgt_free <- c(0, 0, 0)
  } else if(dist == "custom") {
    if(is.null(dist_control$sgt_fixed_values) || is.null(dist_control$sgt_free)) stop("'dist_control$sgt_fixed_values' and 'dist_control$sgt_free' need to be specified if 'dist' is set to 'custom'.")
  } else {
    stop("Allowed values for 'dist' are: 'sgt', 'st' (default, recommended), 't', 'snormal', 'normal' (not recommended), and 'custom'.")
  }
  if(is.null(dist_control$fix_moments)) dist_control$fix_moments <- 2
  if(!(dist_control$fix_moments %in% c(0, 1, 2))) stop("Only allowed values for 'dist_control$fix_moments' are 0, 1 and 2.")
  if(is.null(dist_control$p_q_mins)) {
    if(dist_control$sgt_free[2] == 0 && dist_control$sgt_free[3] == 1) {
      if(dist_control$fix_moments == 2) dist_control$p_q_mins <- c(exp(dist_control$sgt_fixed_values[2]), 2/exp(dist_control$sgt_fixed_values[2]))
      if(dist_control$fix_moments == 1) dist_control$p_q_mins <- c(exp(dist_control$sgt_fixed_values[2]), 1/exp(dist_control$sgt_fixed_values[2]))
      if(dist_control$fix_moments == 0) dist_control$p_q_mins <- c(exp(dist_control$sgt_fixed_values[2]), 0)
    } else if(dist_control$sgt_free[2] == 1 && dist_control$sgt_free[3] == 0) {
      if(dist_control$fix_moments == 2) dist_control$p_q_mins <- c(2/exp(dist_control$sgt_fixed_values[3]), exp(dist_control$sgt_fixed_values[3]))
      if(dist_control$fix_moments == 1) dist_control$p_q_mins <- c(1/exp(dist_control$sgt_fixed_values[3]), exp(dist_control$sgt_fixed_values[3]))
      if(dist_control$fix_moments == 0) dist_control$p_q_mins <- c(0, exp(dist_control$sgt_fixed_values[3]))
    } else {
      if(dist_control$fix_moments == 2) dist_control$p_q_mins <- c(1, 2)
      if(dist_control$fix_moments == 1) dist_control$p_q_mins <- c(1, 1)
      if(dist_control$fix_moments == 0) dist_control$p_q_mins <- c(0, 0)
    }
  }

  #####################################################
  ### Restrictions and priors on parameter matrix B ###
  #####################################################

  if(B_pos_diag) {
    if(!is.null(B_restrictions)) {
      warning("User provided 'B_restrictions' overrides 'B_pos_diag == TRUE'.")
      B_pos_diag <- FALSE
    } else {
      B_restrictions = diag(ncol(y))
      B_restrictions[B_restrictions == 0] <- NA
    }
  }
  if(is.null(B_restrictions)) {
    B_restrictions <- matrix(NA, nrow = M, ncol = M)
  } else {
    not_NA_Br <- B_restrictions[!is.na(B_restrictions)]
    if(length(not_NA_Br) > 0) {
      if(nrow(B_restrictions) != M || ncol(B_restrictions) != M) stop("'B_restrictions' must be a 'ncol(y)' times 'ncol(y)' square matrix.")
      if(sum((not_NA_Br == 0) + (abs(not_NA_Br))) != length(not_NA_Br)) stop("Allowed values in 'B_restrictions' are NA, 0, 1 and -1.")
      if(length(which(B_restrictions == 0)) > 0) {
        B_det_test <- B_restrictions; B_det_test[is.na(B_det_test)] <- 1;
        if(sum(diag(B_det_test) == 0) > 0) stop("Zero restrictions on the diagonal of B not allowed.")
      }
    }
  }
  Bzr <- B_restrictions
  Bzr[abs(Bzr) == 1] <- NA
  Bzr[Bzr == 0] <- 1
  Bzr[is.na(Bzr)] <- 0
  Bzr <- c(Bzr)
  B_zero_res <- cbind(Bzr, cumsum(Bzr))
  Bsr <- B_restrictions
  Bsr[is.na(Bsr)] <- 0
  B_pos_res <- c(Bsr)
  B_pos_res[B_pos_res == -1] <- 0
  B_pos_res <- cbind(B_pos_res, cumsum(B_pos_res))
  B_neg_res <- c(Bsr) * (-1)
  B_neg_res[B_neg_res == -1] <- 0
  B_neg_res <- cbind(B_neg_res, cumsum(B_neg_res))
  if(!is.null(other_priors$B_prior)) B_prior <- other_priors$B_prior else B_prior <- c(2^2, NA)
  if(class(B_prior) == "list") {
    include_hyper_B <- FALSE
    include_B_prior <- TRUE
  } else if(class(B_prior) == "numeric" && length(B_prior) == 2) {
    if(is.na(B_prior[1])) stop("'is.na(B_prior[1]) == TRUE' not allowed.")
    if(is.na(B_prior[2])) include_hyper_B <- 1L else include_hyper_B <- 0L
    if(!include_hyper_B && B_prior[1] == Inf && B_prior[2] == Inf) include_B_prior <- FALSE else include_B_prior <- TRUE
  } else {
    stop("'B_prior' misspecified. Must either be a list or a numeric vector of length two.")
  }
  if(!include_B_prior) {
    B_prior_mean <- c(diag(M))
    B_prior_cov <- diag(M^2)
  } else {
    if(class(B_prior) == "list") {
      if(length(B_prior) != 2) stop("'B_prior' misspecified, length(B_prior) != 2.")
      if(sum(names(B_prior) == c("mean", "cov")) != 2) stop("'B_prior' misspecified, names(B_prior) != c('mean', 'cov').")
      if(length(B_prior$mean) != M^2) stop("'B_prior' misspecified, length(B_prior$mean) != M^2.")
      if(sum(dim(B_prior$cov) == M^2) != 2) stop("'B_prior' misspecified, dim(B_prior$cov) != c(M^2, M^2).")
      if(sum(B_prior$cov == Inf) > 0) stop("'B_prior' misspecified, diagonal of prior covariance matrix allows for only finite values.")
      if(sum(diag(B_prior$cov) <= 0) > 0) stop("'B_prior' misspecified, diagonal of prior covariance matrix allows for only positive values.")
      B_prior_mean <- B_prior$mean
      B_prior_cov <- B_prior$cov
    } else if(class(B_prior) == "numeric") {
      if(!include_hyper_B) {
        if(sum(B_prior[2] <= 0 || B_prior[1] <= 0) > 0) stop("'B_prior' misspecified, diagonal of prior covariance matrix allows for only positive values.")
      }
      B_prior_mean <- c(diag(M))
      if(sum(B_pos_res[,1]) > 0) B_prior_mean[which(B_pos_res[,1] == 1)] <- 1
      if(sum(B_neg_res[,1]) > 0) B_prior_mean[which(B_neg_res[,1] == 1)] <- -1
      B_prior_cov <- diag(M^2)
      diag(B_prior_cov)[diag(M) == 1] <- B_prior[1]
      diag(B_prior_cov)[diag(M) == 0] <- B_prior[2]
      B_prior_cov[is.na(B_prior_cov)] <- -99 # -99 == NA (Stan does not allow for NAs)
    } else {
      stop("'B_prior' misspecified. Must be either a list (prior mean vector and prior covariance matrix) \n or a numeric vector of length two (prior variance of diagonal and off-diagonal elements assuming diagonal prior covariance and prior mean equal to identity matrix).")
    }
  }

  ######################################################
  ### Zero restrictions on autoregressive parameters ###
  ######################################################

  if(is.null(A_restrictions)) {
    A_restrictions <- matrix(NA, nrow = M * lags, ncol = M)
  } else {
    not_NA_Ar <- A_restrictions[!is.na(A_restrictions)]
    if(length(not_NA_Ar) > 0) {
      if(nrow(A_restrictions) != M * lags || ncol(A_restrictions) != M) stop("'A_restrictions' must be a 'ncol(y) * lags' times 'ncol(y)' square matrix.")
      if(sum(not_NA_Ar == 0) != length(not_NA_Ar)) stop("Allowed values in 'A_restrictions' are NA and 0.")
    }
  }
  Azr <- A_restrictions
  Azr[Azr == 0] <- 1
  Azr[is.na(Azr)] <- 0
  Azr <- c(Azr)
  A_zero_res <- cbind(Azr, cumsum(Azr))

  #######################
  ### Minnesota prior ###
  #######################

  # Control arguments
  if(is.null(minnesota_control$minnesota_means)) minnesota_means <- rep(1, M) else minnesota_means <- minnesota_control$minnesota_means
  if(length(minnesota_means) == 1) minnesota_means <- rep(minnesota_means, M)
  if(is.null(minnesota_control$hyper_free)) hyper_free <- c(1, 0, 0) else hyper_free <- minnesota_control$hyper_free
  if(is.null(minnesota_control$hyper_fixed_values)) hyper_fixed_values <- c(0.2, 1, 0.5) else hyper_fixed_values <- minnesota_control$hyper_fixed_values
  if(is.null(minnesota_control$hyper_shrinkage_prior)) hyper_shrinkage_prior <- c(0, 1, 0.001) else hyper_shrinkage_prior <- minnesota_control$hyper_shrinkage_prior
  if(is.null(minnesota_control$hyper_lags_prior)) hyper_lags_prior <- c(0.5^2, 0.5) else hyper_lags_prior <- minnesota_control$hyper_lags_prior # Log-normal prior: c(0.5^2, 0.5) -> mode == 1
  if(is.null(minnesota_control$hyper_ownlags_prior)) hyper_ownlags_prior <- c(4, 4) else hyper_ownlags_prior <- minnesota_control$hyper_ownlags_prior # Beta prior: c(x, x) -> mode == 0.5
  if(is.null(minnesota_control$parallel)) minnesota_parallel <- TRUE else minnesota_parallel <- minnesota_control$parallel
  if(is.null(minnesota_control$no_cross_shrinkage)) minnesota_control$no_cross_shrinkage <- c(1)
  no_cross_shrinkage_vec <- rep(0, M); no_cross_shrinkage_vec[minnesota_control$no_cross_shrinkage] <- 1
  if(is.null(minnesota_control$soc)) minnesota_control$soc <- Inf
  if(is.null(minnesota_control$dio)) minnesota_control$dio <- Inf
  if(is.na(minnesota_control$soc)) minnesota_control$soc <- -99
  if(is.na(minnesota_control$dio)) minnesota_control$dio <- -99
  if(minnesota_control$soc == FALSE) minnesota_control$soc <- Inf
  if(minnesota_control$dio == FALSE) minnesota_control$dio <- Inf
  if(minnesota_control$soc == TRUE) minnesota_control$soc <- -99
  if(minnesota_control$dio == TRUE) minnesota_control$dio <- -99
  if(is.null(minnesota_control$soc_dio_dep)) minnesota_control$soc_dio_dep <- c(0, 0)
  if(is.null(minnesota_control$soc_prior)) if(minnesota_control$soc_dio_dep[1] == 1) minnesota_control$soc_prior <- c(log(5) + 1, 1) else minnesota_control$soc_prior <- c(1, 1)
  if(is.null(minnesota_control$dio_prior)) if(minnesota_control$soc_dio_dep[2] == 1) minnesota_control$dio_prior <- c(log(5) + 1, 1) else minnesota_control$dio_prior <- c(1, 1)
  if(lags == 0) include_minnesota <- FALSE
  if(include_minnesota == FALSE) {
    hyper_free <- c(0, 0, 0)
    minnesota_control$soc <- Inf
    minnesota_control$dio <- Inf
  }
  hyper_shrinkage_lim <- c(hyper_shrinkage_prior[3], Inf)
  hyper_shrinkage_prior <- hyper_shrinkage_prior[1:2]

  # Initial values for the autoregressive parameters according to the Minnesota prior
  if(lags > 0) {
    A0 <- diag(M)
    diag(A0) <- minnesota_means
    if(lags > 1) {
      for(i in 2:lags) {
        A0 <- rbind(A0, matrix(0, ncol = M, nrow = M))
      }
    }
  } else {
    A0 <- as.array(matrix(0, ncol = M, nrow = 0))
  }

  ####################
  ### Other priors ###
  ####################

  if(is.null(other_priors$constant_prior)) other_priors$constant_prior <- c(0, Inf)
  if(is.null(other_priors$lambda_prior)) other_priors$lambda_prior <- c(4, 4) # Beta prior on (lambda + 1 / 2): c(x, x) -> mode == 0
  if(is.null(other_priors$p_prior)) other_priors$p_prior <- c(1, 1) # Log-normal prior: c(1, 1) -> mode == 1 + p_q_mins[1] (e.q. for default sgt 1 + 1 == 2)
  if(is.null(other_priors$q_prior)) other_priors$q_prior <- c(2, 1) # Log-normal prior: c(2, 1) -> mode ~ 2.7 (dof ~ 5.4), mean ~ 12 (dof ~ 24)
  if(is.null(other_priors$p_q_prior_shift)) other_priors$p_q_prior_shift <- dist_control$p_q_mins
  if(is.null(other_priors$hyper_B_prior)) other_priors$hyper_B_prior <- c(0, 1) # Log-normal prior: c(0, 1) -> mode ~ 0.37
  if(is.null(other_priors$hyper_B_soft_prior)) other_priors$hyper_B_soft_prior <- c(4, 4) # Symmetric Beta prior
  if(is.null(other_priors$vol_breaks_prior)) other_priors$vol_breaks_prior <- 2 # Dirichlet prior, either uniform (scalar), or not (vector)
  prior_elast_binary <- rep(0, M * M)
  if(!is.null(other_priors$prior_elasticity)) {
    prior_elast_num <- length(other_priors$prior_elasticity)
    if(prior_elast_num < 1 || class(other_priors$prior_elasticity) != "list") stop("'other_priors$prior_elasticity' must either be NULL or a list of positive length.")
    prior_elast_mat <- matrix(NA, ncol = 5, nrow = prior_elast_num)
    for(i in 1:prior_elast_num) {
      if(length(other_priors$prior_elasticity[[i]]) != 5) stop("Each element of 'other_priors$prior_elasticity' must be a vector of length 5: (numerator, denumerator, mu, sigma > 0, df > 0).")
      prior_elast_mat[i,] <- other_priors$prior_elasticity[[i]]
    }
    if(length(unique(prior_elast_mat[,1:2])) != prior_elast_num * 2) stop("Elements of B (or B_inv) may appear only once in 'other_priors$prior_elasticity'.")
    if(max(prior_elast_mat[,1:2]) > M * M || min(prior_elast_mat[,1:2]) < 1) stop("'other_priors$prior_elasticity' out of bounds.")
    if(min(prior_elast_mat[,4:5]) <= 0) stop("Each element of 'other_priors$prior_elasticity' must be a vector of length 5: (numerator, denumerator, mu, sigma > 0, df > 0).")
    prior_elast_binary[prior_elast_mat[,1]] <- 1
  } else {
    prior_elast_num <- 0
    prior_elast_mat <- matrix(0, ncol = 5, nrow = prior_elast_num)
  }

  ##########################
  ### Missing data (1/2) ###
  ##########################

  if(length(which(is.na(y))) > 0) include_missing_data <- TRUE else include_missing_data <- FALSE
  if(N > 0) if(sum(is.na(y[1,])) > 0) stop("Time series are not allowed to begin with missing values.")

  ###########################
  ### Data transformation ###
  ###########################

  data_transformation <- list("mean" = rep(0, M), "scale" = rep(1, M))
  if(transform == TRUE && N > 0) {
    data_transformation$mean <- apply(y, 2, mean, na.rm = TRUE)
    for(i in 1:ncol(y)) y[,i] <- y[,i] - data_transformation$mean[i]
    xy1 <- build_xy(y, lags)
    if(include_missing_data == TRUE) {
      missing_rows <- which(is.na(apply(cbind(xy1$yy, xy1$xx), 1, sum)))
      xy1$yy <- xy1$yy[-missing_rows,]
      xy1$xx <- xy1$xx[-missing_rows,]
    }
    if(lags > 0) {
      if(ncol(xy1$xx) < nrow(xy1$xx) && no_init_ols == FALSE) {
        A1 <- chol2inv(chol(crossprod(xy1$xx))) %*% t(xy1$xx) %*% xy1$yy
        U1 <- xy1$yy - xy1$xx %*% A1
      } else {
        U1 <- xy1$yy - xy1$xx %*% A0
      }
    } else {
      U1 <- xy1$yy
    }
    S1 <- crossprod(U1) / nrow(U1)
    data_transformation$scale <- sqrt(diag(S1))
    for(i in 1:ncol(y)) y[,i] <- y[,i] / data_transformation$scale[i]
    if(!is.null(attributes(y)$quarterly)) {
      yq <- y[,attributes(y)$quarterly]
      if(is.null(nrow(yq))) yq <- matrix(yq, ncol = 1)
      for(i in 1:length(quarterly_sums)) {
        sub <- yq[(i * 3 - 2):(i * 3)]
        quarterly_sums[i] <- sum(exp(sub)) # Log-levels assumed
      }
    }
  }
  xy <- xy2 <- build_xy(y, lags)

  #################################################
  ### Initial values for the parameter matrix B ###
  #################################################

  if(include_missing_data == TRUE) {
    missing_rows <- which(is.na(apply(cbind(xy2$yy, xy2$xx), 1, sum)))
    xy2$yy <- xy2$yy[-missing_rows,]
    xy2$xx <- xy2$xx[-missing_rows,]
  }
  if(lags > 0) {
    if(ncol(xy2$xx) < nrow(xy2$xx) && no_init_ols == FALSE) {
      A2 <- chol2inv(chol(crossprod(xy2$xx))) %*% t(xy2$xx) %*% xy2$yy
      U2 <- xy2$yy - xy2$xx %*% A2
    } else {
      U2 <- xy2$yy - xy2$xx %*% A0
    }
  } else {
    U2 <- xy2$yy
  }
  S2 <- crossprod(U2) / nrow(U2)
  if(N > 0) B0 <- expm::sqrtm(S2) else B0 <- diag(M)
  if(is.complex(B0)) B0 <- Re(B0)
  B0_inv <- solve(B0)
  E0 <- U2 %*% t(B0_inv)
  if(B_inverse == TRUE) B0 <- B0_inv

  ############################################################
  ### Initial values for the shock distribution parameters ###
  ############################################################

  sgt0 <- matrix(NA, ncol = M, nrow = 3)
  par0 <- c(dist_control$sgt_fixed_values[dist_control$sgt_free == 1])
  if(length(par0) > 0) {
    mean_cent <- as.logical(ifelse(dist_control$fix_moments > 0, 1, 0))
    var_adj <- as.logical(ifelse(dist_control$fix_moments > 1, 1, 0))
    sgt_obj <- function(par, E_col, dist_control, mean_cent, var_adj) {
      sgt_par <- c(dist_control$sgt_fixed_values)
      sgt_par[which(dist_control$sgt_free == 1)] <- par
      out <- -sum(sgt::dsgt(E_col, mu = 0, sigma = 1, lambda = sgt_par[1], p = exp(sgt_par[2]), q = exp(sgt_par[3]),
                            mean.cent = mean_cent, var.adj = var_adj, log = T))
      if(is.nan(out)) out <- 10000000
      out
    }
    for(i in 1:M) {
      if(N == 0) {
        sgt0[, i] <- dist_control$sgt_fixed_values
        next
      }
      sgt0[which(dist_control$sgt_free == 1), i] <- stats::optim(par0,
                                                                 sgt_obj,
                                                                 method = "L-BFGS-B",
                                                                 lower = c(-0.9, log(dist_control$p_q_mins[1] + 0.1), log(dist_control$p_q_mins[2] + 0.1))[which(dist_control$sgt_free == 1)],
                                                                 upper = c(0.9, 5, 5)[which(dist_control$sgt_free == 1)],
                                                                 E_col = E0[,i],
                                                                 dist_control = dist_control,
                                                                 mean_cent = mean_cent,
                                                                 var_adj = var_adj)$par
    }
  }

  ##########################
  ### Missing data (2/2) ###
  ##########################

  missing_data_location <- matrix(NA, ncol = 2, nrow = sum(is.na(y)))
  missing_data_raw_init <- as.array(rep(0, sum(is.na(y))))
  missing_data_count <- 0
  if(include_missing_data == TRUE) {
    for(col in 1:ncol(y)) {
      for(row in 1:nrow(y)) {
        if(is.na(y[row, col])) {
          missing_data_count <- missing_data_count + 1
          missing_data_location[missing_data_count,] <- c(row, col)
        }
      }
    }
  } else {
    missing_data_raw_init <- as.array(rep(0, 0))
    missing_data_location <- matrix(NA, ncol = 2, nrow = sum(is.na(y)))
  }

  ################################
  ### GARCH / Shock volatility ###
  ################################

  if(is.null(garch_control$garch_prior)) garch_control$garch_prior <- 10/12
  if(is.null(garch_control$garch_dep)) garch_control$garch_dep <- FALSE
  if(is.null(garch_control$garch_shrinkage)) garch_control$garch_shrinkage <- 1
  if(is.null(garch_control$garch_groups)) include_garch_groups <- FALSE else include_garch_groups <- TRUE
  if(is.null(garch_control$garch_eta_form)) garch_control$garch_eta_form <- FALSE

  if(garch_control$garch_prior <= 0 || garch_control$garch_prior >= 1) {
    stop("'garch_control$garch_prior' not valid (must have a value between 0 and 1).")
  }
  if(garch_control$garch_shrinkage <= 0) stop("'garch_control$garch_shrinkage' > 0 required.")
  if(include_garch) {

    if(!garch_control$garch_dep) {

      # Independent shock volatility
      if(!include_garch_groups) {
        garch_init <- matrix(NA, nrow = M, ncol = 3)
        garch_init[,2] <- garch_control$garch_prior
        garch_init[,-2] <- (1 - garch_control$garch_prior) / (ncol(garch_init) - 1)

      # Grouped shock volatility (partly dependent)
      } else {
        if(length(garch_control$garch_groups) == 0) stop("'length(garch_control$garch_groups) == 0' not allowed.")
        if(1 %in% sapply(garch_control$garch_groups, length)) stop("Elements of length one in 'garch_control$garch_groups' not allowed.\n  If you need groups of size one, just remove those groups from 'garch_control$garch_groups'.")
        garch_group_mat <- diag(M)
        for(i in 1:length(garch_control$garch_groups)) {
          garch_group_mat <- cbind(garch_group_mat, apply(garch_group_mat[,garch_control$garch_groups[[i]]], 1, sum))
        }
        garch_group_mat <- garch_group_mat[,-unlist(garch_control$garch_groups)]
        if(is.null(ncol(garch_group_mat))) garch_group_mat <- matrix(garch_group_mat, ncol = 1)
        garch_group_num <- ncol(garch_group_mat)
        garch_init <- matrix(NA, nrow = garch_group_num, ncol = 3)
        garch_init[,2] <- garch_control$garch_prior
        garch_init[,-2] <- (1 - garch_control$garch_prior) / (ncol(garch_init) - 1)
      }
      garch_control$garch_prior <- c(1, 2 / (1 - garch_control$garch_prior) - 2, 1)

    # Dependent shock volatility
    } else {
      garch_init <- matrix(NA, nrow = M, ncol = M + 2)
      garch_init[,2] <- garch_control$garch_prior
      garch_init[,-2] <- (1 - garch_control$garch_prior) / (ncol(garch_init) - 1)
      garch_control$garch_prior <- c(1, (M + 1) / (1 - garch_control$garch_prior) - (M + 1), rep(1, M))
    }

  } else {
    garch_init <- matrix(NA, nrow = 0, ncol = ifelse(garch_control$garch_dep == 0, 3, M + 2))
    garch_control$garch_prior <- rep(0, 0)
  }
  if(!include_garch_groups) {
    garch_group_mat <- matrix(0, nrow = 0, ncol = 0)
    garch_group_num <- ncol(garch_group_mat)
  }
  garch_control$garch_prior <- garch_control$garch_prior / garch_control$garch_shrinkage

  #########################
  ### Volatility breaks ###
  #########################

  if(include_vol_breaks == TRUE) {
    if(!("ts" %in% class(y))) stop("'y' must be a time series object if 'vol_breaks' is not NULL.")
    if(length(vol_breaks) < 1) stop("'vol_breaks' misspecified (length(vol_breaks) < 1).")
    vol_breaks_binary <- rep(0, N)
    tsi <- ts(1:N, start = attributes(y)$tsp[1], frequency = attributes(y)$tsp[3])
    for(i in 1:length(vol_breaks)) {
      break_index <- as.numeric(window(tsi, start = vol_breaks[[i]])[1]) - lags
      if(break_index < 3) stop("'vol_breaks' misspecified. You have probably set the first the break too early. \n Note that the first break must be at least (lags + 2) periods after the first period in 'y'.")
      vol_breaks_binary[break_index] <- 1
    }
  } else {
    vol_breaks_binary <- rep(0, N)
  }
  if(length(other_priors$vol_breaks_prior) != sum(vol_breaks_binary) + 1 && length(other_priors$vol_breaks_prior) != 1) stop("'other_priors$vol_breaks_prior' misspecified, needs to be a scalar or of length sum(vol_breaks) + 1.")
  if(sum(other_priors$vol_breaks_prior <= 0) > 0) stop("'other_priors$vol_breaks_prior' misspecified, needs to be positive.")
  if(length(other_priors$vol_breaks_prior) == 1) other_priors$vol_breaks_prior <- rep(other_priors$vol_breaks_prior, sum(vol_breaks_binary) + 1)
  init_relative_vol <- matrix(1, nrow = ncol(y) * include_vol_breaks, ncol = sum(vol_breaks_binary) + 1) / (sum(vol_breaks_binary) + 1)

  ##############################
  ### Collect initial values ###
  ##############################

  inits <- list("A" = A0, "B" = B0, "sgt" = sgt0, "constant" = rep(0, M),
                "garch_param" = garch_init,
                "relative_vol_raw" = init_relative_vol,
                "missing_data_raw" = missing_data_raw_init,
                "monthly_raw" = init_monthly)
  if(hyper_free[1] != 0) inits$hyper_shrinkage <- hyper_fixed_values[1]
  if(hyper_free[2] != 0) inits$hyper_lags <- hyper_fixed_values[2]
  if(hyper_free[3] != 0) inits$hyper_ownlags <- hyper_fixed_values[3]
  if(minnesota_control$soc < 0) if(minnesota_control$soc_dio_dep[1] == 1) inits$hyper_soc <- log(5) else inits$hyper_soc <- log(1)
  if(minnesota_control$dio < 0) if(minnesota_control$soc_dio_dep[2] == 1) inits$hyper_dio <- log(5) else inits$hyper_dio <- log(1)
  if(include_hyper_B) inits$hyper_B <- exp(other_priors$hyper_B_prior[1] - other_priors$hyper_B_prior[2]^2) # Mode of log-normal default prior

  ########################
  ### Collect standata ###
  ########################

  y_raw <- y # Only relevant if 'include_missing_data == TRUE'
  if(include_missing_data == TRUE) {
    y_raw[which(is.na(y_raw))] <- Inf # Here: Inf is a placeholder for NA (Stan does not allow for NAs)
    xy$yy[which(is.na(xy$yy))] <- Inf
    xy$xx[which(is.na(xy$xx))] <- Inf
  }
  B_prior_cov_diagonal <- ifelse(sum(B_prior_cov[c(which(lower.tri(B_prior_cov)), which(upper.tri(B_prior_cov)))]) == 0, 1, 0)
  if(B_prior_cov_diagonal == 0) stop("Non-diagonal B prior covariance matrix not yet supported.")
  standata <- list(y = xy$yy,
                   x = xy$xx,
                   y_raw = y_raw,
                   lags = lags,
                   N = N,
                   M = M,
                   include_constant = include_constant,
                   B_inverse = B_inverse,

                   # Distribution of the shocks
                   fix_moments = dist_control$fix_moments,
                   sgt_len_vec = dist_control$sgt_free * M,
                   sgt_fixed = dist_control$sgt_fixed_values,
                   p_q_mins = dist_control$p_q_mins,
                   pq_len = M * !(dist_control$sgt_free[2] == 0 && dist_control$sgt_free[3] == 0),

                   # Restrictions and priors on B
                   B_zero_res_num = B_zero_res[length(B_zero_res)],
                   B_pos_res_num = B_pos_res[length(B_pos_res)],
                   B_neg_res_num = B_neg_res[length(B_neg_res)],
                   B_zero_res = as.array(B_zero_res),
                   B_pos_res = as.array(B_pos_res),
                   B_neg_res = as.array(B_neg_res),
                   include_B_prior = include_B_prior,
                   B_prior_mean = as.array(B_prior_mean),
                   B_prior_cov = as.array(B_prior_cov),
                   B_prior_cov_diagonal = B_prior_cov_diagonal,

                   # Restrictions on A and the Minnesota prior
                   A_zero_res_num = A_zero_res[length(A_zero_res)],
                   A_zero_res = as.array(A_zero_res),
                   include_minnesota = include_minnesota,
                   hyper_len_vec = as.array(hyper_free),
                   hyper_fixed = hyper_fixed_values,
                   minnesota_means = as.array(minnesota_means),
                   hyper_shrinkage_prior = hyper_shrinkage_prior,
                   hyper_lags_prior = hyper_lags_prior,
                   hyper_ownlags_prior = hyper_ownlags_prior,
                   hyper_shrinkage_lim = hyper_shrinkage_lim,
                   include_soc = ifelse(minnesota_control$soc == Inf, 0, 1),
                   include_dio = ifelse(minnesota_control$dio == Inf, 0, 1),
                   soc_free = ifelse(minnesota_control$soc < 0, 1, 0),
                   dio_free = ifelse(minnesota_control$dio < 0, 1, 0),
                   soc = minnesota_control$soc,
                   dio = minnesota_control$dio,
                   soc_prior = minnesota_control$soc_prior,
                   dio_prior = minnesota_control$dio_prior,
                   soc_dio_dep = minnesota_control$soc_dio_dep,

                   # GARCH / Shock volatility
                   include_garch = include_garch,
                   garch_dependence = garch_control$garch_dep,
                   garch_prior = garch_control$garch_prior,
                   include_garch_groups = include_garch_groups,
                   garch_group_num = garch_group_num,
                   garch_group_mat = as.array(garch_group_mat),
                   garch_eta_form = garch_control$garch_eta_form,

                   # Volatility breaks
                   include_vol_breaks = include_vol_breaks,
                   vol_breaks_num = sum(vol_breaks_binary),
                   vol_breaks = vol_breaks_binary,
                   vol_breaks_prior = as.array(other_priors$vol_breaks_prior),

                   # Other priors
                   constant_prior = other_priors$constant_prior,
                   lambda_prior = other_priors$lambda_prior,
                   p_prior = other_priors$p_prior,
                   q_prior = other_priors$q_prior,
                   p_q_prior_shift = other_priors$p_q_prior_shift,
                   include_hyper_B = include_hyper_B,
                   hyper_B_prior = other_priors$hyper_B_prior,
                   prior_elast_binary = prior_elast_binary,
                   prior_elast_num = prior_elast_num,
                   prior_elast_ind = matrix(prior_elast_mat[,1:2], nrow = prior_elast_num, ncol = 2),
                   prior_elast_par = matrix(prior_elast_mat[,3:5], nrow = prior_elast_num, ncol = 3),

                   # Missing data
                   include_missing_data = include_missing_data,
                   missing_data_count = missing_data_count,
                   missing_data_location = missing_data_location,

                   # Interpolation of quarterly data
                   number_of_quarterly = number_of_quarterly,
                   quarterly_binary = quarterly_binary,
                   quarterly_sums = quarterly_sums,

                   # Other stuff
                   parallel_likelihood = ifelse(threads_per_chain == 1, 0, 1),
                   minnesota_parallel = minnesota_parallel,
                   data_transformation_mean = data_transformation$mean,
                   data_transformation_scale = data_transformation$scale,
                   hyperbolic_transformation = 0,
                   no_cross_shrinkage_vec = no_cross_shrinkage_vec
  )

  ###########################################################
  ### Function that provides Stan with the initial values ###
  ###########################################################

  if(is.null(init_fun)) {
    init_fun <- function() {
      list2env(standata, envir = environment())
      list2env(inits, envir = environment())
      if(lags > 0) {
        if(A_zero_res[length(A_zero_res)] > 0) {
          for(i in 1:nrow(A_zero_res)) if(A_zero_res[i, 1] == 1) A[i] <- NA
        }
      }
      if(B_zero_res[length(B_zero_res)] > 0) {
        for(i in 1:nrow(B_zero_res)) if(B_zero_res[i, 1] == 1) B[i] <- NA
      }
      B_pos <- rep(0.01, 0)
      if(B_pos_res[length(B_pos_res)] > 0) {
        for(i in 1:nrow(B_pos_res)) {
          if(B_pos_res[i, 1] == 1) {
            if(B0[i] <= 0) B_pos <- c(B_pos, 0.01) else B_pos <- c(B_pos, B[i])
            B[i] <- NA
          }
        }
      }
      B_neg <- rep(-0.01, 0)
      if(B_neg_res[length(B_neg_res)] > 0) {
        for(i in 1:nrow(B_neg_res)) {
          if(B_neg_res[i, 1] == 1) {
            if(B[i] >= 0) B_neg <- c(B_neg, -0.01) else B_neg <- c(B_neg, B[i])
            B[i] <- NA
          }
        }
      }
      if(dist_control$sgt_free[1] == 1) lambda0 <- sgt0[1,] else lambda0 <- rep(0, 0)
      if(dist_control$sgt_free[2] == 1) log_p0 <- sgt0[2,] else log_p0 <- rep(0, 0)
      if(dist_control$sgt_free[3] == 1) log_q0 <- sgt0[3,] else log_q0 <- rep(0, 0)
      if(include_constant == FALSE) constant <- rep(0, 0)
      hypers <- c(hyper_fixed, 1)
      if(hyper_len_vec[1] != 0) hypers[1] <- inits$hyper_shrinkage
      if(hyper_len_vec[2] != 0) hypers[2] <- inits$hyper_lags
      if(hyper_len_vec[3] != 0) hypers[3] <- inits$hyper_ownlags
      hyper_soc <- soc
      hyper_dio <- dio
      if(soc_free == 1) hyper_soc <- inits$hyper_soc
      if(dio_free == 1) hyper_dio <- inits$hyper_dio
      if(include_hyper_B) hypers[4] <- inits$hyper_B
      list(lambda = as.array(lambda0),
           log_p = as.array(log_p0),
           log_q = as.array(log_q0),
           constant = as.array(constant),
           B_par = c(B)[!is.na(B)],
           B_pos = as.array(B_pos),
           B_neg = as.array(B_neg),
           hyper_shrinkage = as.array(rep(hypers[1], hyper_len_vec[1])),
           hyper_lags = as.array(rep(hypers[2], hyper_len_vec[2])),
           hyper_ownlags = as.array(rep(hypers[3], hyper_len_vec[3])),
           hyper_soc = as.array(rep(hyper_soc, soc_free)),
           hyper_dio = as.array(rep(hyper_dio, dio_free)),
           hyper_B = as.array(rep(hypers[4], include_hyper_B)),
           garch_param = as.array(inits$garch_param),
           relative_vol_raw = as.array(inits$relative_vol_raw),
           missing_data_raw = as.array(inits$missing_data_raw),
           A_par = c(A)[!is.na(A)],
           monthly_raw = as.array(inits$monthly_raw))
    }
  } else {
    if(class(init_fun) == "list") {
      user_inits <- init_fun
      init_fun <- function() {
        list2env(user_inits, envir = environment())
        user_inits
      }
    }
  }

  #################################################
  ### Initial point estimation / MLE estimation ###
  #################################################

  if(init_optim == TRUE) {
    cat("Initial value optimization... \n")
    standata_optim <- standata
    if(init_hyperpar_optim == FALSE) { # NB: 'hyper_B_soft' not included here for now
      standata_optim$hyper_len_vec[] <- 0
      standata_optim$include_soc <- 0
      standata_optim$include_dio <- 0
      standata_optim$include_hyper_B <- 0
      if(sum(standata_optim$B_prior_cov == -99) != 0) standata_optim$B_prior_cov[standata_optim$B_prior_cov == -99] <- inits$hyper_B
    }
    opt <- mod$optimize(data = standata_optim,
                        init = init_fun,
                        threads = threads_per_chain,
                        iter = optim_iter,
                        refresh = 5000)
    opt_par <- opt$mle()
    inits$A <- matrix(unname(opt_par[grep("A[", names(opt_par), fixed = TRUE)]), ncol = M)
    if(length(grep("hyper_B[", names(opt_par), fixed = TRUE)) > 0) inits$B <- matrix(unname(opt_par[grep("B[", names(opt_par), fixed = TRUE)[-1]]), ncol = M)
    if(length(grep("hyper_B[", names(opt_par), fixed = TRUE)) == 0) inits$B <- matrix(unname(opt_par[grep("B[", names(opt_par), fixed = TRUE)]), ncol = M)
    if(dist_control$sgt_free[1] == 1) inits$sgt[1,] <- opt_par[grep("lambda", names(opt_par))]
    if(dist_control$sgt_free[2] == 1) inits$sgt[2,] <- opt_par[grep("log_p", names(opt_par))]
    if(dist_control$sgt_free[3] == 1) inits$sgt[3,] <- opt_par[grep("log_q", names(opt_par))]
    if(include_constant == TRUE) inits$constant <- unname(opt_par[grep("constant", names(opt_par))])
    if(include_garch == TRUE) {
      inits$garch_param <- matrix(unname(opt_par[grep("garch_param", names(opt_par))]), nrow = nrow(inits$garch_param))
      for(i in 1:nrow(inits$garch_param)) inits$garch_param[i,] <- inits$garch_param[i,] / sum(inits$garch_param[i,])
    }
    if(include_vol_breaks == TRUE) {
      inits$relative_vol_raw[,] <- unname(opt_par[grep("relative_vol_raw", names(opt_par))])
      for(i in 1:nrow(inits$relative_vol_raw)) inits$relative_vol_raw[i,] <- inits$relative_vol_raw[i,] / sum(inits$relative_vol_raw[i,])
    }
    if(include_missing_data == TRUE) inits$missing_data <- unname(opt_par[grep("missing_data_raw", names(opt_par))])
    if(init_hyperpar_optim == TRUE) {
      if(standata$hyper_len_vec[1] != 0) inits$hyper_shrinkage <- unname(opt_par[grep("hyper_shrinkage", names(opt_par))])
      if(standata$hyper_len_vec[2] != 0) inits$hyper_lags <- unname(opt_par[grep("hyper_lags", names(opt_par))])
      if(standata$hyper_len_vec[3] != 0) inits$hyper_ownlags <- unname(opt_par[grep("hyper_ownlags", names(opt_par))])
      if(soc_free == 1) inits$hyper_soc <- unname(opt_par[grep("hyper_soc", names(opt_par))])
      if(dio_free == 1) inits$hyper_dio <- unname(opt_par[grep("hyper_dio", names(opt_par))])
      if(standata$include_hyper_B != 0) inits$hyper_B <- unname(opt_par[grep("hyper_B[", names(opt_par), fixed = TRUE)])
      if(standata$B_soft_res_num > 0) inits$hyper_B_soft <- unname(opt_par[grep("hyper_B_soft[", names(opt_par), fixed = TRUE)])
    }
  }

  #############################################
  ### Initial guess for posterior variances ###
  #############################################

  if(init_inv_metric[1] == "auto") {
    inits_val <- init_fun()
    inits_var <- c()
    param_names <- names(mod$variables()$parameters)
    for(i in 1:length(param_names)) {
      param_name <- param_names[i]
      param_inits <- unlist(inits_val[param_name])
      if(length(param_inits) != 0) {
        if(param_name == "A") { # !!! (NB:  'A_restrictions')
          inits_var <- c(inits_var, c(minnesota_sds(M, lags, hyper_fixed_values))^2)
        } else if(param_name %in% c("garch_param", "relative_vol_raw")) {
          if(param_name == "garch_param" && include_garch_groups) {
            inits_var <- c(inits_var, rep(0.2^2, length(param_inits) -  garch_group_num))
          } else {
            inits_var <- c(inits_var, rep(0.2^2, length(param_inits) -  M))
          }
        } else if(param_name == "monthly_raw") {
          inits_var <- c(inits_var, rep(0.2^2, length(param_inits) * (2 / 3)))
        } else {
          inits_var <- c(inits_var, rep(0.2^2, length(param_inits)))
        }
      }
    }
    init_inv_metric <- inits_var
    if(metric == "dense_e") init_inv_metric <- diag(init_inv_metric)
  }

  ##########################
  ### Posterior sampling ###
  ##########################

  if(!skip_sampling) {
    cmdfit <- mod$sample(
      data = standata,
      init = init_fun,
      chains = chains,
      parallel_chains = parallel_chains,
      threads_per_chain = threads_per_chain,
      inv_metric = init_inv_metric,
      metric = metric,
      save_warmup = TRUE,
      ...
    )
    stanfit <- tryCatch(
      {
        rstan::read_stan_csv(cmdfit$output_files())
      },
      error = function(cond) {
        cmdfit$save_output_files()
        warning("There was a problem with the construction of the rstan object. \n Output files have been saved in the current directory.")
        message("Original error message:")
        message(cond)
        return(NA)
      }
    )
  } else {
    cmdfit <- NA
    stanfit <- NA
  }

  attributes(stanfit)$cmdfit <- cmdfit
  attributes(stanfit)$standata <- standata
  attributes(stanfit)$data_transformation <- data_transformation
  attributes(stanfit)$original_y <- y
  attributes(stanfit)$init_fun_out <- init_fun()
  if(init_optim == TRUE) attributes(stanfit)$init_optim <- opt
  stanfit
}
