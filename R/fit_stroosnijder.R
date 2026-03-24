#' Fit the Stroosnijder Infiltration Model
#'
#' Fits the Stroosnijder (1976) infiltration model to observed soil water
#' infiltration data using nonlinear least squares with bounded optimization
#' via \code{minpack.lm::nlsLM}. The function provides automatic starting values,
#' parameter estimation, model diagnostics, and commonly used goodness-of-fit statistics.
#'
#' The Stroosnijder model is defined as:
#' \deqn{f(t) = \frac{S}{2 \sqrt{t + \tau}} + K}
#' where \eqn{S} is sorptivity, \eqn{K} is the saturated hydraulic conductivity,
#' and \eqn{\tau} is an empirical time-shift parameter that stabilizes early-time behavior.
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
#'     \item \code{params} – Data frame of estimated parameters
#'       (\eqn{S}, \eqn{K}, \eqn{\tau}).
#'     \item \code{fitted} – Data frame with time, observed values,
#'       predicted values, and residuals.
#'     \item \code{smooth} – Data frame of fine-grid predictions for
#'       plotting (predicted infiltration rates and cumulative infiltration).
#'     \item \code{stats} – Data frame of goodness-of-fit
#'       statistics (RMSE, MAE, NRMSE, PBIAS, R², NSE, KGE, AIC).
#'     \item \code{nls_obj} – The fitted \code{nlsLM} model object.
#'   }
#'
#' @details
#' Starting parameter values are automatically estimated from the data
#' using heuristic rules to improve convergence stability. Sorptivity (\eqn{S})
#' is initialized relative to the range of observed infiltration rates,
#' conductivity (\eqn{K}) from late-time observations, and the time-shift
#' parameter (\eqn{\tau}) is initialized to a small positive value.
#' Parameter bounds are constrained to non-negative values.
#'
#' The function computes standard hydrological performance metrics
#' including Nash–Sutcliffe Efficiency (NSE), Percent Bias (PBIAS),
#' Kling–Gupta Efficiency (KGE), and Normalized Root Mean Square Error (NRMSE).
#'
#' @references
#' Stroosnijder, L. (1976). Infiltration measurements on a sandy soil
#' with a rainfall simulator. *Netherlands Journal of Agricultural Science*, 24, 55–67.
#'
#' Brutsaert, W. (1977). The infiltration process in water resources modeling.
#' *Water Resources Research*, 13(3), 637–644.
#'
#' Philip, J. R. (1957). The theory of infiltration: 4. Sorptivity and
#' algebraic infiltration equations. *Soil Science*, 84(3), 257–264.
#'
#' @export
#'
#' @importFrom minpack.lm nlsLM nls.lm.control
#' @importFrom stats predict AIC coef cor
#' @importFrom utils tail
#' @importFrom stats approx
#'
#' @examples
#' # Example using synthetic data
#' set.seed(123)
#' time <- seq(5, 60, by = 5)
#'
#' # Simulated infiltration curve with Stroosnijder form
#' S <- 30; K <- 5; tau <- 2
#' rate <- S / (2 * sqrt(time + tau)) + K
#' rate <- rate + rnorm(length(rate), 0, 0.5)  # add noise
#'
#' data <- data.frame(time = time, rate = rate)
#'
#' model <- fit_stroosnijder(data, time = time, rate = rate)
#'
#' model$params
#' model$stats

fit_stroosnijder <- function(data = NULL, time, rate, na.rm = TRUE) {

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
  # Handle missing values
  # ---------------------------
  if (na.rm) {
    ok <- stats::complete.cases(time, infil)
    time  <- time[ok]
    infil <- infil[ok]
  }

  # ---------------------------
  # Checks
  # ---------------------------
  if (!is.numeric(time) || !is.numeric(infil))
    stop("time and rate must be numeric.")
  if (length(time) < 5)
    stop("At least 5 observations are required.")
  if (any(time <= 0))
    stop("time must be > 0.")

  # ---------------------------
  # Ensure increasing time
  # ---------------------------
  ord <- order(time)
  time  <- time[ord]
  infil <- infil[ord]

  # ---------------------------
  # Stroosnijder equation
  # f(t) = S / (2 * sqrt(t + tau)) + K
  # ---------------------------
  stroos_fun <- function(t, S, K, tau) {
    S / (2 * sqrt(t + tau)) + K
  }

  # ---------------------------
  # Starting values
  # ---------------------------
  K_start   <- mean(tail(infil, max(3, floor(length(infil) * 0.1))))
  S_start   <- diff(range(infil)) * sqrt(mean(time))
  tau_start <- 0.1

  K_start   <- max(K_start, 1e-6)
  S_start   <- max(S_start, 1e-6)
  tau_start <- max(tau_start, 1e-6)

  # ---------------------------
  # Nonlinear fit
  # ---------------------------
  fit <- minpack.lm::nlsLM(
    infil ~ stroos_fun(time, S, K, tau),
    start = list(S = S_start, K = K_start, tau = tau_start),
    lower = c(S = 0, K = 0, tau = 0),
    control = minpack.lm::nls.lm.control(maxiter = 500)
  )

  # ---------------------------
  # Predictions & residuals
  # ---------------------------
  pred  <- stats::predict(fit)
  resid <- infil - pred

  # ---------------------------
  # Parameters
  # ---------------------------
  S_fit   <- coef(fit)[["S"]]
  K_fit   <- coef(fit)[["K"]]
  tau_fit <- coef(fit)[["tau"]]

  # ---------------------------
  # Smooth predictions on fine grid
  # ---------------------------
  new_time <- seq(min(time), max(time), length.out = 200)
  new_rate <- stroos_fun(new_time, S_fit, K_fit, tau_fit)
  new_cum  <- S_fit * sqrt(new_time + tau_fit) + K_fit * new_time   # cumulative infiltration

  smooth_pred <- data.frame(
    time = new_time,
    predicted_rate = new_rate,
    predicted_cum  = new_cum
  )

  # ---------------------------
  # Performance metrics
  # ---------------------------
  rmse  <- sqrt(mean(resid^2))
  mae   <- mean(abs(resid))
  nrmse <- rmse / (max(infil) - min(infil)) * 100
  pbias <- 100 * sum(resid) / sum(infil)

  nse <- 1 - sum(resid^2) / sum((infil - mean(infil))^2)
  r2  <- cor(infil, pred)^2

  r     <- cor(infil, pred)
  alpha <- sd(pred) / sd(infil)
  beta  <- mean(pred) / mean(infil)
  kge   <- 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)

  aic <- stats::AIC(fit)

  # ---------------------------
  # Structured tidy output
  # ---------------------------
  list(
    model = "Stroosnijder",
    params = data.frame(S = S_fit, K = K_fit, tau = tau_fit),
    fitted = data.frame(
      time      = time,
      observed  = infil,
      predicted = pred,
      residual  = resid
    ),
    smooth = smooth_pred,   # fine-grid predictions
    stats = data.frame(
      RMSE  = rmse,
      MAE   = mae,
      NRMSE = nrmse,
      PBIAS = pbias,
      R2    = r2,
      NSE   = nse,
      KGE   = kge,
      AIC   = aic
    ),
    nls_obj = fit
  )
}
