#' Fit the Philip Infiltration Model
#'
#' Fits the two-parameter Philip infiltration model to observed soil water
#' infiltration experimental data using nonlinear least squares
#' (Levenberg–Marquardt algorithm) via \code{minpack.lm::nlsLM()}. The
#' function provides automatic starting values, parameter estimation, model
#' diagnostics, and commonly used goodness-of-fit statistics.
#'
#' The Philip model describes infiltration rate as a function of sorptivity
#' and steady-state infiltration rate:
#' \deqn{f(t) = f_c + \frac{1}{2} S t^{-1/2}}
#'
#' where \eqn{f_c} is the steady-state infiltration rate and
#' \eqn{S} is sorptivity.
#'
#' @param data A data frame containing infiltration observations.
#' @param time Unquoted column name representing cumulative time
#'   (must be numeric and strictly > 0).
#' @param rate Unquoted column name representing infiltration rate
#'   (numeric).
#' @param na.rm Logical; if `TRUE` (default), incomplete cases are removed
#'   before fitting.
#'
#' @returns A list with class \code{"philip_model"} containing:
#' \itemize{
#'   \item \code{model} — Model name.
#'   \item \code{parameters} — Estimated parameters (\eqn{f_c}, \eqn{S}).
#'   \item \code{fitted} — Data frame with observed, predicted, and residual values.
#'   \item \code{performance} — Data frame of goodness-of-fit
#'       statistics (RMSE, MAE, NRMSE, PBIAS, R², NSE, AIC).
#'   \item \code{nls_obj} — The fitted \code{nlsLM} model object.
#' }
#'
#' @details
#' Starting values are estimated internally:
#' \itemize{
#'   \item \eqn{f_c} is approximated from the tail of the observed series.
#'   \item \eqn{S} is estimated via linear regression of
#'   \eqn{f(t) - f_c} against \eqn{1/\sqrt{t}}.
#' }
#'
#' Time values must be strictly positive because the model contains
#' a \eqn{t^{-1/2}} term.
#'
#' The function computes standard hydrological performance metrics
#' including Nash–Sutcliffe Efficiency (NSE), Percent Bias (PBIAS),
#' and Normalized Root Mean Square Error (NRMSE).
#'
#' @references
#' Philip, J. R. (1957). The theory of infiltration: 1. The infiltration equation
#' and its solution. *Soil Science*, 83(5), 345–358.
#'
#' @export
#'
#' @examples
#' # Example using synthetic data
#' set.seed(123)
#' time <- seq(5, 60, by = 5)
#' rate <- 5 + 0.5 * 12 * time^(-0.5) + rnorm(length(time), 0, 0.3)
#'
#' data <- data.frame(time = time, rate = rate)
#'
#' model <- fit_philip(data, time = time, rate = rate)
#'
#' print(model$parameters)
#' print(model$performance)

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

    performance = data.frame(
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
