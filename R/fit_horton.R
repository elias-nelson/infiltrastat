#' Fit the Horton Infiltration Model
#'
#' Fits the Horton (1940) exponential decay model to observed soil
#' infiltration rate data using nonlinear least squares with bounded
#' optimization via \code{minpack.lm::nlsLM}. The function provides
#' automatic starting values, parameter estimation, model diagnostics,
#' and commonly used goodness-of-fit statistics.
#'
#' The Horton model is defined as:
#' \deqn{f(t) = f_c + (f_0 - f_c) e^{-k t}}
#' where \eqn{f_0} is the initial infiltration rate,
#' \eqn{f_c} is the equilibrium (final) infiltration rate,
#' and \eqn{k} is the decay constant.
#'
#' @param data Optional data frame containing the variables specified
#'   in \code{time} and \code{rate}. If supplied, unquoted column names
#'   should be provided for \code{time} and \code{rate}.
#' @param time Numeric vector of time values (must be increasing),
#'   or unquoted column name if \code{data} is provided.
#' @param rate Numeric vector of observed infiltration rates,
#'   or unquoted column name if \code{data} is provided.
#' @param na.rm Logical; if \code{TRUE} (default), missing values are
#'   removed before model fitting.
#'
#' @returns A structured list with class \code{"infiltra_model"}
#'   containing:
#'   \itemize{
#'     \item \code{model} – Character string identifying the model.
#'     \item \code{parameters} – Data frame of estimated parameters
#'       (\eqn{f_c}, \eqn{f_0}, \eqn{k}).
#'     \item \code{fitted} – Data frame with time, observed values,
#'       predicted values, and residuals.
#'     \item \code{performance} – Data frame of goodness-of-fit
#'       statistics (RMSE, MAE, NRMSE, PBIAS, R², NSE, AIC).
#'     \item \code{nls_obj} – The fitted \code{nlsLM} model object.
#'   }
#'
#' @details
#' Starting parameter values are automatically estimated from the data
#' using heuristic rules to improve convergence stability. The decay
#' constant is initialized via linearization of the exponential form
#' when possible. Parameter bounds are constrained to non-negative values.
#'
#' The function computes standard hydrological performance metrics
#' including Nash–Sutcliffe Efficiency (NSE), Percent Bias (PBIAS),
#' and Normalized Root Mean Square Error (NRMSE).
#'
#' @export
#'
#' @importFrom minpack.lm nlsLM nls.lm.control
#' @importFrom stats predict AIC coef lm
#' @importFrom utils tail
#'
#' @examples
#' # Example using synthetic data
#' set.seed(123)
#' time <- seq(0, 60, by = 5)
#' rate <- 5 + (20 - 5) * exp(-0.08 * time) + rnorm(length(time), 0, 0.5)
#'
#' model <- fit_horton(time = time, rate = rate)
#'
#' model$parameters
#' model$performance

fit_horton <- function(data = NULL, time, rate, na.rm = TRUE) {

  # ---------------------------
  # Extract from data if provided
  # ---------------------------
  if (!is.null(data)) {
    time  <- data[[deparse(substitute(time))]]
    infil <- data[[deparse(substitute(rate))]]
  } else {
    infil <- rate
  }

  # ---------------------------
  # Input checks
  # ---------------------------
  if (na.rm) {
    ok <- stats::complete.cases(time, infil)
    time  <- time[ok]
    infil <- infil[ok]
  }

  if (!is.numeric(time) || !is.numeric(infil))
    stop("time and rate must be numeric.")

  if (length(time) < 5)
    stop("At least 5 observations are required.")

  # ---------------------------
  # Ensure time is increasing
  # ---------------------------
  ord <- order(time)
  time  <- time[ord]
  infil <- infil[ord]

  # ---------------------------
  # Horton equation
  # ---------------------------
  horton_fun <- function(t, fc, f0, k) {
    fc + (f0 - fc) * exp(-k * t)
  }

  # ---------------------------
  # Smart starting values
  # ---------------------------
  f0_start <- infil[1]
  fc_start <- mean(tail(infil, max(3, floor(length(infil) * 0.1))))

  eps <- 1e-6
  valid <- infil > fc_start + eps

  if (sum(valid) >= 2) {
    y_lin <- log(infil[valid] - fc_start)
    k_start <- -coef(stats::lm(y_lin ~ time[valid]))[2]
  } else {
    k_start <- 0.01
  }

  k_start <- max(k_start, 1e-4)

  # ---------------------------
  # Nonlinear fit with bounds
  # ---------------------------
  fit <- minpack.lm::nlsLM(
    infil ~ horton_fun(time, fc, f0, k),
    start = list(fc = fc_start, f0 = f0_start, k = k_start),
    lower = c(fc = 0, f0 = 0, k = 0),
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
  rmse  <- sqrt(mean(resid^2))
  mae   <- mean(abs(resid))
  nrmse <- rmse / (max(infil) - min(infil)) * 100
  pbias <- 100 * sum(resid) / sum(infil)

  r2  <- 1 - sum(resid^2) / sum((infil - mean(infil))^2)
  nse <- r2

  aic <- stats::AIC(fit)

  # ---------------------------
  # Structured tidy output
  # ---------------------------
  list(
    model = "Horton",

    parameters = data.frame(
      fc = coef(fit)[["fc"]],
      f0 = coef(fit)[["f0"]],
      k  = coef(fit)[["k"]]
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
