fit_philip <- function(data, time, rate, na.rm = TRUE) {

  # ---------------------------
  # Tidy evaluation of columns
  # ---------------------------
  time  <- data[[deparse(substitute(time))]]
  infil <- data[[deparse(substitute(rate))]]

  # ---------------------------
  # Handle NA values
  # ---------------------------
  if (na.rm) {
    ok <- stats::complete.cases(time, infil)
    time  <- time[ok]
    infil <- infil[ok]
  }

  # ---------------------------
  # Basic checks
  # ---------------------------
  if (!is.numeric(time) || !is.numeric(infil))
    stop("time and rate must be numeric.")

  if (length(time) < 5)
    stop("At least 5 observations are required.")

  # ---------------------------
  # Ensure increasing time
  # ---------------------------
  ord <- order(time)
  time  <- time[ord]
  infil <- infil[ord]

  # Avoid zero time (division by sqrt)
  if (any(time <= 0))
    stop("All time values must be > 0 for the Philip model.")

  # ---------------------------
  # Philip infiltration equation
  # f(t) = fc + 0.5 * S * t^(-1/2)
  # ---------------------------
  philip_fun <- function(t, fc, S) {
    fc + 0.5 * S * t^(-0.5)
  }

  # ---------------------------
  # Smart starting values
  # ---------------------------

  # steady infiltration from tail
  fc_start <- mean(tail(infil, max(3, floor(length(infil) * 0.1))))

  # sorptivity via correct linearization:
  # (f - fc) ~ 1/sqrt(t)
  z <- 1 / sqrt(time)
  y <- infil - fc_start

  valid <- y > 0
  if (sum(valid) >= 2) {
    slope <- coef(stats::lm(y[valid] ~ z[valid]))[2]
    S_start <- 2 * slope
  } else {
    S_start <- diff(range(infil))  # safe fallback
  }

  S_start <- max(S_start, 1e-6)

  # ---------------------------
  # Nonlinear fit (Levenberg–Marquardt)
  # ---------------------------
  fit <- minpack.lm::nlsLM(
    infil ~ philip_fun(time, fc, S),
    start = list(fc = fc_start, S = S_start),
    lower = c(fc = 0, S = 0),
    control = minpack.lm::nls.lm.control(maxiter = 200)
  )

  # ---------------------------
  # Predictions & residuals
  # ---------------------------
  pred  <- stats::predict(fit)
  resid <- infil - pred

  # ---------------------------
  # Goodness-of-fit statistics
  # ---------------------------
  rmse <- sqrt(mean(resid^2))
  mae  <- mean(abs(resid))
  nrmse <- rmse / (max(infil) - min(infil)) * 100
  pbias <- 100 * sum(resid) / sum(infil)

  r2  <- 1 - sum(resid^2) / sum((infil - mean(infil))^2)
  nse <- 1 - sum(resid^2) / sum((infil - mean(infil))^2)

  aic <- stats::AIC(fit)

  # ---------------------------
  # Structured tidy output
  # ---------------------------
  list(
    model = "Philip",

    parameters = data.frame(
      fc = coef(fit)[["fc"]],
      S  = coef(fit)[["S"]]
    ),

    fitted = data.frame(
      time      = time,
      observed  = infil,
      predicted = pred,
      residual  = resid
    ),

    perfomance = data.frame(
      RMSE  = rmse,
      MAE   = mae,
      NRMSE = nrmse,
      PBIAS = pbias,
      R2    = r2,
      NSE   = nse,
      AIC   = aic
    ),

    nls_obj = fit
  )
}
