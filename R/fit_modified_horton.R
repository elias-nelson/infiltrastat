#' Fit the Modified Horton Infiltration Model
#'
#' Fits the Modified Horton exponential decay model to observed soil water
#' infiltration experimental data using nonlinear least squares with bounded
#' optimization via \code{minpack.lm::nlsLM}. The function provides
#' automatic starting values, parameter estimation, model diagnostics,
#' and commonly used goodness-of-fit statistics.
#'
#' The Modified Horton model is defined as:
#' \deqn{f(t) = f_c + (f_0 - f_c) e^{-k t^n}}
#' where \eqn{f_0} is the initial infiltration rate,
#' \eqn{f_c} is the equilibrium (final) infiltration rate,
#' \eqn{k} is the decay constant,
#' and \eqn{n} is an empirical exponent that generalizes the rate of decay.
#' When \eqn{n = 1}, the model reduces to the classical Horton equation.
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
#'       (\eqn{f_c}, \eqn{f_0}, \eqn{k}, \eqn{n}).
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
#' using heuristic rules to improve convergence stability. The decay
#' constant is initialized via linearization of the exponential form
#' when possible. The exponent \eqn{n} is initialized to 1, corresponding
#' to the classical Horton model. Parameter bounds are constrained to
#' non-negative values.
#'
#' The function computes standard hydrological performance metrics
#' including Nash–Sutcliffe Efficiency (NSE), Percent Bias (PBIAS),
#' Kling–Gupta Efficiency (KGE), and Normalized Root Mean Square Error (NRMSE).
#'
#' @references
#' Horton, R. E. (1940). An approach toward a physical interpretation of
#' infiltration capacity. *Soil Science Society of America Proceedings*, 5, 399–417.
#'
#' Brutsaert, W. (1977). The infiltration process in water resources modeling.
#' *Water Resources Research*, 13(3), 637–644.
#'
#' Singh, V. P., & Yu, F. X. (1990). Derivation of infiltration equations
#' using the concept of variable infiltration capacity.
#' *Water Resources Bulletin*, 26(2), 253–262.
#'
#' @export
#'
#' @importFrom minpack.lm nlsLM nls.lm.control
#' @importFrom stats predict AIC coef lm cor
#' @importFrom utils tail
#'
#' @examples
#' # Example using synthetic data
#' set.seed(123)
#' time <- seq(5, 60, by = 5)
#'
#' # Simulated infiltration curve with modified Horton form
#' f0 <- 50; fc <- 5; k <- 0.02; n <- 1.2
#' rate <- fc + (f0 - fc) * exp(-k * time^n)
#' rate <- rate + rnorm(length(rate), 0, 0.5)  # add noise
#'
#' data <- data.frame(time = time, rate = rate)
#'
#' model <- fit_modified_horton(data, time = time, rate = rate)
#'
#' model$params
#' model$stats
fit_modified_horton <- function(data = NULL, time, rate, na.rm = TRUE) {

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
  # Modified Horton equation
  # f(t) = fc + (f0 - fc) * exp(-k * t^n)
  # ---------------------------
  mod_horton_fun <- function(t, fc, f0, k, n) {
    fc + (f0 - fc) * exp(-k * t^n)
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
  n_start <- 1   # default Horton decay exponent

  # ---------------------------
  # Nonlinear fit
  # ---------------------------
  fit <- minpack.lm::nlsLM(
    infil ~ mod_horton_fun(time, fc, f0, k, n),
    start = list(fc = fc_start, f0 = f0_start, k = k_start, n = n_start),
    lower = c(fc = 0, f0 = 0, k = 0, n = 0),
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
  n_fit  <- coef(fit)[["n"]]

  # ---------------------------
  # Smooth predictions on fine grid
  # ---------------------------
  new_time <- seq(min(time), max(time), length.out = 200)
  new_rate <- mod_horton_fun(new_time, fc_fit, f0_fit, k_fit, n_fit)

  smooth_pred <- data.frame(
    time = new_time,
    predicted_rate = new_rate
  )

  # ---------------------------
  # Goodness-of-fit statistics
  # ---------------------------
  rmse  <- sqrt(mean(resid^2))
  mae   <- mean(abs(resid))
  nrmse <- rmse / (max(infil) - min(infil)) * 100
  pbias <- 100 * sum(resid) / sum(infil)

  # NSE (Nash & Sutcliffe, 1970)
  nse <- 1 - sum(resid^2) / sum((infil - mean(infil))^2)

  # R² as squared correlation
  r2  <- cor(infil, pred)^2

  # KGE (Gupta et al. 2009)
  r     <- cor(infil, pred)
  alpha <- sd(pred) / sd(infil)
  beta  <- mean(pred) / mean(infil)
  kge   <- 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)

  # AIC
  aic <- stats::AIC(fit)

  # ---------------------------
  # Structured tidy output
  # ---------------------------
  list(
    model = "Modified Horton",
    params = data.frame(fc = fc_fit, f0 = f0_fit, k = k_fit, n = n_fit),
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
