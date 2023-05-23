
#' @export
dist_plot <- function(fit,
                      quantiles = c(0.68, 0.90),
                      shock_names = NULL,
                      skew_line = 0,
                      pq_line = 5,
                      skew_lim = NULL,
                      pq_lim = NULL,
                      col = "darkblue",
                      side_by_side = TRUE,
                      mar = c(5, 5, 2, 2),
                      override_graph = FALSE,
                      apply_restriction = TRUE) {

  if(!override_graph) par(mar = mar)
  M <- attributes(fit)$standata$M
  sgt_len_vec <- attributes(fit)$standata$sgt_len_vec
  if(is.null(shock_names)) shock_names <- paste("Shock", 1:M)
  if(sum(sgt_len_vec) == 0) {
    cat("dist_plot: Nothing to plot.\n")
    return(NULL)
  }
  if(length(quantiles) != 2) stop("length(quantiles) != 2.")
  if(min(quantiles) < 0 || max(quantiles) > 1) stop("Values in 'quantiles' must be between 0 and 1.")
  quants <- c( (1 - max(quantiles)) / 2, (1 - min(quantiles)) / 2, 1 - (1 - min(quantiles)) / 2, 1 - (1 - max(quantiles)) / 2 )

  # Get posterior draws
  fit_draws = as.data.frame(fit)
  res <- attributes(fit)$restriction
  if(apply_restriction && !is.null(res)) {
    N <- nrow(fit_draws) - sum(res)
    if(sgt_len_vec[1] > 0) fit_draws[1:N, grep("lambda", colnames(fit_draws))] <- extract(fit, apply_restriction = apply_restriction)$lambda
    if(sgt_len_vec[2] > 0 || sgt_len_vec[3] > 0) fit_draws[1:N, grep("pq", colnames(fit_draws))] <- extract(fit, apply_restriction = apply_restriction)$pq
    fit_draws <- fit_draws[1:N,]
  }

  # Fix ggplot shock ordering for large models (hacky...)
  if(M >= 10 && sgt_len_vec[1] > 0) {
    fill_with_these <- rep(c(9, 99, 999, 9999, 99999), each = 10)
    correct_these <- grep("lambda", colnames(fit_draws))[-c(1:9)]
    colnames(fit_draws)[correct_these] <- paste0(substr(colnames(fit_draws)[correct_these], 1, 7),
                                                 fill_with_these[1:length(correct_these)],
                                                 substr(colnames(fit_draws)[correct_these], 8, 10))
  }
  if(M >= 10 && (sgt_len_vec[2] > 0 || sgt_len_vec[3] > 0)) {
    fill_with_these <- rep(c(9, 99, 999, 9999, 99999), each = 10)
    correct_these <- grep("pq", colnames(fit_draws))[-c(1:9)]
    colnames(fit_draws)[correct_these] <- paste0(substr(colnames(fit_draws)[correct_these], 1, 3),
                                                 fill_with_these[1:length(correct_these)],
                                                 substr(colnames(fit_draws)[correct_these], 4, 6))
  }

  # Generate prior samples
  if(sgt_len_vec[1] > 0) {
    lambda_prior <- attributes(fit)$standata$lambda_prior
    prior_sample_lambda <- 2 * rbeta(nrow(fit_draws), lambda_prior[1], lambda_prior[2]) - 1
  }
  if(sgt_len_vec[2] > 0 || sgt_len_vec[3] > 0) {
    q_prior <- attributes(fit)$standata$q_prior; q_fixed <- exp(attributes(fit)$standata$sgt_fixed[2])
    p_prior <- attributes(fit)$standata$p_prior; p_fixed <- exp(attributes(fit)$standata$sgt_fixed[3])
    p_q_mins <- attributes(fit)$standata$p_q_mins
    if(sgt_len_vec[2] > 0 && sgt_len_vec[3] > 0) {
      prior_sample_pq <- (exp(rnorm(nrow(fit_draws), p_prior[1], p_prior[2])) + p_q_mins[1]) * (exp(rnorm(nrow(fit_draws), q_prior[1], q_prior[2])) + p_q_mins[2])
    }
    if(sgt_len_vec[2] > 0 && sgt_len_vec[3] == 0) {
      prior_sample_pq <- (exp(rnorm(nrow(fit_draws), p_prior[1], p_prior[2])) + p_q_mins[1]) * q_fixed
    }
    if(sgt_len_vec[2] == 0 && sgt_len_vec[3] > 0) {
      prior_sample_pq <- (exp(rnorm(nrow(fit_draws), q_prior[1], q_prior[2])) + p_q_mins[2]) * p_fixed
    }
  }

  if(side_by_side && !override_graph) par(mfrow = c(1, 2))
  if(!side_by_side && !override_graph) par(mfrow = c(2, 1))

  # Skewness parameters
  if(length(grep("lambda", names(fit))) > 0) {
    fit_lambda <- fit_draws[, grep("lambda[", colnames(fit_draws), fixed = TRUE)]
    fit_lambda$zzz_prior <- prior_sample_lambda
    if(is.null(skew_lim)) skew_lim <- c(-0.6, 0.6)
    plot(0, 0, type = "n", ylim = c(0, M + 1), xlim = skew_lim,
         xlab = "Skewness parameter", ylab = "", yaxt = "n", bty = "l")
    axis(2, at = seq(0, M, by = 1) + 0.5, labels = c(shock_names, "Prior"), las = 1)
    for(i in 1:ncol(fit_lambda)) {
      q <- quantile(fit_lambda[,i], quants)
      rect(q[2], i - 0.75, q[3], i - 0.25, border = NA, col = col)
      lines(cbind(c(q[1], q[4]), c(i - 0.5, i - 0.5)), lwd = 4, col = col)
    }
    abline(v = skew_line, lwd = 1.5, lty = 2)
  }

  if(length(grep("pq", names(fit))) > 0) {
    fit_pq <- fit_draws[, grep("pq[", colnames(fit_draws), fixed = TRUE)]
    fit_pq$zzz_prior <- prior_sample_pq
    if(is.null(pq_lim)) pq_lim <- c(p_q_mins[1] * p_q_mins[2], 30)
    plot(0, 0, type = "n", ylim = c(0, M + 1), xlim = pq_lim,
         xlab = "Degree-of-freedom parameter", ylab = "", yaxt = "n", xaxt = "n", bty = "l")
    axis(2, at = seq(0, M, by = 1) + 0.5, labels = c(shock_names, "Prior"), las = 1)
    axis(1, at = c(2, 10, 30, pq_line), labels = c(2, 10, 30, pq_line), las = 1)
    for(i in 1:ncol(fit_pq)) {
      q <- quantile(fit_pq[,i], quants)
      rect(q[2], i - 0.75, q[3], i - 0.25, border = NA, col = col)
      lines(cbind(c(q[1], q[4]), c(i - 0.5, i - 0.5)), lwd = 4, col = col)
    }
    abline(v = pq_line, lwd = 1.5, lty = 2)
  }

}

