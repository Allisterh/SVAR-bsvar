
#' @export
plot_volatility <- function(obj,
                            y_labels = NULL,
                            probs = NULL,
                            apply_restriction = TRUE,
                            first_period = NA,
                            freq = NA,
                            control = list()) {

  if("stanfit" %in% class(obj)) {
    decomp <- shock_decomp(obj, apply_restriction = apply_restriction)
  } else {
    decomp <- obj
  }
  if(is.null(control$col)) control$col <- "tomato"
  if(is.null(probs)) probs <- c(0.05, 0.16, 0.84, 0.95)
  probs <- c(probs[1] - 0.01, probs, probs[length(probs)] + 0.01)
  if(min(probs) < 0 || max(probs) > 1) stop("'probs' must have values in between 0.01 and 0.99.")
  if(!is.null(decomp$tsp_)) {
    first_period <- decomp$tsp_[1]; freq <- decomp$tsp_[3]
  } else {
    if(is.na(first_period[1])) first_period <- 1
    if(is.na(freq)) freq <- 1
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
