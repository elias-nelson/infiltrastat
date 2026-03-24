#' Fit the Holtan Infiltration Model
#'
#' Fits the Holtan (1961) semi-empirical infiltration model to observed
#' soil water infiltration data using nonlinear least squares with bounded
#' optimization via \code{minpack.lm::nlsLM}. The function provides
#' automatic starting values, parameter estimation, model diagnostics,
#' and commonly used goodness-of-fit statistics.
#'
#' The Holtan model is defined as:
#' \deqn{f(t) = f_c + (f_0 - f_c) e^{-k t} + m P(t)}
#' where \eqn{f_0} is the initial infiltration rate,
#' \eqn{f_c} is the equilibrium (final) infiltration rate,
#' \eqn{k} is the decay constant,
#' \eqn{m} is a soil moisture storage coefficient,
#' and \eqn{P(t)} is rainfall intensity (constant or time-varying).
#'
#' @param data Optional data frame containing the variables specified
#'   in \code{time}, \code{rate}, and optionally \code{rain}.
#'   If supplied, unquoted column names should be provided.
#' @param time Numeric vector of time values (must be increasing),
#'   or unquoted column name if \code{data} is provided.
#' @param rate Numeric vector of observed infiltration rates,
#'   or unquoted column name if \code{data} is provided.
#' @param rain Optional numeric vector or column name representing
#'   rainfall intensity (mm/hr or cm/hr). If not supplied, defaults to zero.
#' @param na.rm Logical; if \code{TRUE} (default), missing values are
#'   removed before model fitting.
#'
#' @returns A structured list with class \code{"infiltra_model"}
#'   containing:
#'   \itemize{
#'     \item \code{model} – Character string identifying the model.
#'     \item \code{params} – Data frame of estimated parameters
#'       (\eqn{f_c}, \eqn{f_0}, \eqn{k}, \eqn{m}).
#'     \item \code{fitted} – Data frame with time, observed values,
#'       predicted values, and residuals.
#'     \item \code{smooth} – Data frame of fine-grid predictions for
#'       plotting (predicted infiltration rates).
#'     \item \code{stats} – Data frame of goodness-of-fit
#'       statistics (RMSE, MAE, NRMSE, PBIAS, R², NSE, KGE, AIC).
#'     \item \code{nls_obj} – The fitted \code{nlsLM} model object.
#'   }
#'
#' @details
#' Starting parameter values are automatically estimated from the data
#' using heuristic rules to improve convergence stability. The equilibrium
#' infiltration rate \eqn{f_c} is initialized from late-time observations,
#' while \eqn{f_0} is taken from early-time observations. The rainfall
#' intensity term \eqn{P(t)} can be constant or time-varying. Parameter
#' bounds are constrained to non-negative values.
#'
#' The function computes standard hydrological performance metrics
#' including Nash–Sutcliffe Efficiency (NSE), Percent Bias (PBIAS),
#' Kling–Gupta Efficiency (KGE), and Normalized Root Mean Square Error (NRMSE).
#'
#' @references
#' Holtan, H. N. (1961). A concept for infiltration estimates in watershed engineering.
#' *USDA Agricultural Research Service Publication*, ARS 41-51.
#'
#' Brutsaert, W. (1977). The infiltration process in water resources modeling.
#' *Water Resources Research*, 13(3), 637–644.
#'
#' USDA Soil Conservation Service (1972). National Engineering Handbook,
#' Section 4: Hydrology. U.S. Department of Agriculture, Washington, D.C.
#'
#' @export
#'
#' @importFrom minpack.lm nlsLM nls.lm.control
#' @importFrom stats predict AIC coef cor
#' @importFrom utils tail
#'
#' @examples
#' # Example using synthetic data
#' set.seed(123)
#' time <- seq(5, 60, by = 5)
#' rain <- rep(10, length(time))  # constant rainfall intensity
#'
#' # Simulated infiltration curve with Holtan form
#' f0 <- 40; fc <- 5; k <- 0.02; m <- 0.1
#' rate <- fc + (f0 - fc) * exp(-k * time) + m * rain
#' rate <- rate + rnorm(length(rate), 0, 0.5)  # add noise
#'
#' data <- data.frame(time = time, rate = rate, rain = rain)
#'
#' model <- fit_holtan(data, time = time, rate = rate, rain = rain)
#'
#' model$params
#' model$stats
fit_holtan <- function(data = NULL, time, rate, rain = NULL, na.rm = TRUE) {

  # ---------------------------
  # Extract from data if provided
  # ---------------------------
  if (!is.null(data)) {
    time  <- data[[deparse(substitute(time))]]
    infil <- data[[deparse(substitute(rate))]]
    if (!is.null(substitute(rain))) {
      rain <- data[[deparse(substitute(rain))]]
    }
  } else {
    infil <- rate
  }

  # ---------------------------
  # Handle missing values
  # ---------------------------
  if (na.rm) {
    ok <- stats::complete.cases(time, infil)
    if (!is.null(rain)) ok <- ok & stats::complete.cases(rain)
    time  <- time[ok]
    infil <- infil[ok]
    if (!is.null(rain)) rain <- rain[ok]
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
  if (!is.null(rain)) rain <- rain[ord]

  # ---------------------------
  # Holtan equation
  # f(t) = fc + (f0 - fc) * exp(-k * t) + m * P(t)
  # ---------------------------
  holtan_fun <- function(t, fc, f0, k, m, P) {
    fc + (f0 - fc) * exp(-k * t) + m * P
  }

  # ---------------------------
  # Starting values
  # ---------------------------
  f0_start <- infil[1]
  fc_start <- mean(tail(infil, max(3, floor(length(infil) * 0.1))))
  k_start  <- 0.01
  m_start  <- ifelse(is.null(rain), 0, 0.1)

  # ---------------------------
  # Nonlinear fit
  # ---------------------------
  if (is.null(rain)) rain <- rep(0, length(time))

  fit <- minpack.lm::nlsLM(
    infil ~ holtan_fun(time, fc, f0, k, m, rain),
    start = list(fc = fc_start, f0 = f0_start, k = k_start, m = m_start),
    lower = c(fc = 0, f0 = 0, k = 0, m = 0),
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
  fc_fit <- coef(fit)[["fc"]]
  f0_fit <- coef(fit)[["f0"]]
  k_fit  <- coef(fit)[["k"]]
  m_fit  <- coef(fit)[["m"]]

  # ---------------------------
  # Smooth predictions on fine grid
  # ---------------------------
  new_time <- seq(min(time), max(time), length.out = 200)
  new_rain <- if (length(unique(rain)) == 1) rep(mean(rain), length(new_time)) else approx(time, rain, new_time)$y
  new_rate <- holtan_fun(new_time, fc_fit, f0_fit, k_fit, m_fit, new_rain)

  smooth_pred <- data.frame(
    time = new_time,
    predicted_rate = new_rate
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
    model = "Holtan",
    params = data.frame(fc = fc_fit, f0 = f0_fit, k = k_fit, m = m_fit),
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
