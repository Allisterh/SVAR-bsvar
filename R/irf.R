
#' \code{irf} computes and plots impulse response function for one shock.
#'
#' @export
irf <- function(obj,
                horizon = 40,
                sub_N = NULL,
                shock = 1,
                variables = "all",
                cumulative = c(),
                shock_size = 1,
                fixed_impact = 0,
                transform = NULL,
                controls = list(),
                apply_restriction = TRUE,
                R_progress_bar = TRUE) {

  if(is.null(shock)) stop("Please, specify with respect to which shock to compute the IRFs.")
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
    if(variables[1] == "all") variables <- c(1:M)

    dt <- attributes(fit)$data_transformation
    if(is.null(transform)) transform <- unname(!is.na(dt$mean[1]))
    if(transform == TRUE && !is.na(dt$mean[1])) {
      if(is.null(controls$varscale)) controls$varscale <- rep(1, M)
      controls$varscale <- controls$varscale * dt$scale
    }

    # R implementation:
    e <- rep(0, M)
    e[shock] <- shock_size
    irf_array <- array(NA, dim = c(M, horizon + 1, N))
    if(R_progress_bar) pb <- txtProgressBar(max = N, style = 3)
    for(sample_index in 1:N) {
      Bi <- B[sample_index,,]
      Ai <- A[sample_index,,]
      AAi <- stackA(t(Ai)) # utils.R
      if(fixed_impact > 0) {
        e[shock] <- shock_size / Bi[fixed_impact, shock] / controls$varscale[shock]
      } else {
        e[shock] <- shock_size / controls$varscale[shock]
      }
      for(h in 1:(horizon + 1)) {
        if(h == 1) {
          zero <- Bi %*% e
          zero_long <- c(zero, rep(0, (ncol(AAi) - length(zero))))
          irf_array[,h,sample_index] <- zero
          AAi_power <- AAi
        } else {
          irf_array[,h,sample_index] <- (AAi_power %*% zero_long)[1:length(zero)]
          if(h < (horizon + 1)) AAi_power <- AAi_power %*% AAi
        }
      }
      if(R_progress_bar) setTxtProgressBar(pb, sample_index)
    }
    if(R_progress_bar) close(pb)
    if(length(cumulative) > 0) {
      for(i in 1:length(cumulative)) {
        for(j in 1:N) {
          irf_array[cumulative[i],,j] <- cumsum(irf_array[cumulative[i],,j])
        }
      }
    }

    if(!is.null(controls$varscale)) {
      for(j in 1:M) {
        irf_array[j,,] <- irf_array[j,,] * controls$varscale[j]
      }
    }

  } else if(class(obj) == "bsvar_irf") {
    irf_array <- obj
    shock <- attributes(irf_array)$shock
    if(variables[1] == "all") variables <- c(1:dim(irf_array)[1])
  } else {
    stop("Unknown obj argument.")
  }

  if(is.null(controls$varnames)) {
    if(class(obj) == "stanfit") {
      if(!is.null(colnames(attributes(obj)$standata$y))) {
        if(length(colnames(attributes(obj)$standata$y)) == nrow(irf_array)) {
          controls$varnames <- colnames(attributes(obj)$standata$y)
        }
      }
    }
  }
  if(!is.null(controls$noplot)) if(controls$noplot == TRUE) {
    class(irf_array) <- "bsvar_irf"
    attributes(irf_array)$shock <- shock
    return(invisible(irf_array))
  }
  if(is.null(controls$col)) controls$col <- "tomato"
  if(is.null(controls$median)) controls$median <- FALSE
  if(is.null(controls$prob)) controls$prob <- c(0.05, 0.16, 0.84, 0.95)
  if(is.null(controls$override_graph)) if(is.null(controls$mar)) controls$mar <- c(2, 2.5, 2, 1)
  if(is.null(controls$varnames)) {
    varnames <- paste0("Variable ", 1:length(variables))
  } else {
    varnames <- controls$varnames
  }
  if(is.null(controls$shockname)) controls$shockname <- paste0("Shock ", shock)
  if(is.null(controls$ylab)) controls$ylab <- ""
  if(is.null(controls$line)) controls$line <- 3
  if(is.null(controls$override_graph)) if(is.null(controls$mfrow)) controls$mfrow <- c(ceiling(sqrt(length(variables))), ceiling(sqrt(length(variables))))
  if(is.null(controls$shock_rows)) controls$shock_rows <- TRUE

  probs <- controls$prob
  if(max(probs) > 0.99 | min(probs) < 0.01) stop("Values in 'controls$probs' need to be between 0.01 and 0.99.")
  probs <- c(probs[1] - 0.01, probs, probs[length(probs)] + 0.01)
  if(is.null(controls$override_graph)) par(mar = controls$mar)
  if(is.null(controls$override_graph)) par(mfrow = controls$mfrow)
  if(!is.null(controls$main)) main_text <- controls$main
  for(variable in variables) {
    if(is.null(controls$main)) main_text <- paste0(controls$shockname, " on ", varnames[variable]) else main_text <- controls$main
    if(is.null(controls$ylab)) ylab_text <- "" else ylab_text <- controls$ylab
    sub_irfs <- t(irf_array[variable,,])
    median_sub_irfs <- ts(apply(sub_irfs, 2, median), start = 0)
    quant <- function(column) quantile(column, probs = probs)
    quantiles_sub_irfs <- apply(sub_irfs, 2, quant)
    if(is.null(controls$ylims)) {
      ylim <- c(min(quantiles_sub_irfs), max(quantiles_sub_irfs))
    } else {
      if(controls$shock_rows) {
        ylim <- c(controls$ylims[[1]][variable, shock], controls$ylims[[2]][variable, shock])
      } else {
        ylim <- c(controls$ylims[[1]][shock, variable], controls$ylims[[2]][shock, variable])
      }
    }
    if(length(ylab_text) > 1) {
      if(controls$shock_rows) ylab_index <- shock else ylab_index <- variable
    } else {
      ylab_index <- 1
    }
    if(length(main_text) > 1) {
      if(controls$shock_rows) main_index <- shock else main_index <- variable
    } else {
      main_index <- 1
    }
    plot(unname(median_sub_irfs), lwd = 2, lty = 2, col = controls$col, yaxt = "n", xlab = "", ylab = "",
         main = main_text[main_index], ylim = ylim)
    axis(side = 2, las = 2, mgp = c(3, 0.75, 0), ylab = ylab_text[ylab_index])
    title(ylab = ylab_text[ylab_index], line = controls$line)
    grid()
    fanplot::fan(data = quantiles_sub_irfs, data.type = "values", probs = probs,
                 start = 0, fan.col = colorRampPalette(c(controls$col, "white")),
                 rlab = NULL, ln = NULL)
    if(controls$median) lines(median_sub_irfs, lty = 2)
    abline(h = 0, lwd = 2, lty = 2)
  }
  if(is.null(controls$override_graph)) par(mfrow = c(1,1))

  class(irf_array) <- "bsvar_irf"
  attributes(irf_array)$shock <- shock
  invisible(irf_array)
}

#' \code{irfs} computes and plots impulse response functions for multiple shocks.
#'
#' @export
irfs <- function(obj,
                 horizon = 40,
                 sub_N = NULL,
                 shocks = "all",
                 variables = "all",
                 cumulative = c(),
                 shock_sizes = 1,
                 fixed_impacts = 0,
                 transform = NULL,
                 controls = list(),
                 apply_restriction = TRUE,
                 R_progress_bar = TRUE) {

  if(is.null(shocks)) stop("Please, specify with respect to which shocks to compute the IRFs.")
  if(shocks[1] == "all" && class(obj) == "stanfit") shocks <- 1:attributes(obj)$standata$M
  if(variables[1] == "all" && class(obj) == "stanfit") variables <- 1:attributes(obj)$standata$M
  if(shocks[1] == "all" && class(obj) == "bsvar_irfs") shocks <- 1:length(obj)
  if(variables[1] == "all" && class(obj) == "bsvar_irfs") variables <- 1:dim(obj[[1]])[1]
  if(class(obj) == "bsvar_irfs") shocks <- attributes(obj)$shocks
  override_graph <- controls$override_graph
  shock_rows <- controls$shock_rows
  if(is.null(override_graph)) override_graph <- FALSE
  if(is.null(shock_rows)) shock_rows <- TRUE
  controls$override_graph <- TRUE
  if(is.null(controls$mfrow)) controls$mfrow <- c(length(shocks), length(variables))
  if(is.null(controls$mar)) controls$mar <- c(2, 2.5, 2, 1)
  irf_arrays <- list()
  if(!override_graph) par(mar = controls$mar)
  if(!override_graph) if(shock_rows) par(mfrow = controls$mfrow) else par(mfcol = controls$mfrow)
  if(length(shock_sizes) == 1) shock_sizes <- rep(shock_sizes, length(shocks))
  if(length(fixed_impacts) == 1) fixed_impacts <- rep(fixed_impacts, length(shocks))
  for(i in 1:length(shocks)) {
    if(!is.null(controls$varnames)) controls$varname <- controls$varnames[shocks[i]]
    if(!is.null(controls$shocknames)) controls$shockname <- controls$shocknames[shocks[i]]
    if(!is.null(controls$ylabs)) controls$ylab <- controls$ylabs[[shocks[i]]]
    if(!is.null(controls$mains)) controls$main <- controls$mains[[shocks[i]]]
    cat(paste0("Computing IRFs... (", i, "/", length(shocks), ") \n"))
    if(class(obj) == "stanfit") irf_arrays[[i]] <- irf(obj, horizon, sub_N, shocks[i], variables, cumulative, shock_sizes[i], fixed_impacts[i], transform, controls, apply_restriction, R_progress_bar)
    if(class(obj) == "bsvar_irfs") irf_arrays[[i]] <- irf(obj[i], horizon, sub_N, shocks[i], variables, cumulative, shock_sizes[i], fixed_impacts[i], transform, controls, apply_restriction, R_progress_bar)
  }
  if(!override_graph) par(mfrow = c(1, 1))
  class(irf_arrays) <- "bsvar_irfs"
  attributes(irf_arrays)$shocks <- shocks
  invisible(irf_arrays)
}


