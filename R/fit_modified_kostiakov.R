#' Fit the Modified Kostiakov Infiltration Model
#'
#' Fits the Modified Kostiakov infiltration model to observed soil water
#' infiltration experimental data using nonlinear least squares with bounded
#' optimization via \code{minpack.lm::nlsLM}. The function provides
#' automatic starting values, parameter estimation, model diagnostics,
#' and commonly used goodness-of-fit statistics.
#'
#' The Modified Kostiakov model is defined as:
#' \deqn{f(t) = k a t^{a-1} + f_c}
#' where \eqn{k} and \eqn{a} are empirical constants,
#' and \eqn{f_c} is the equilibrium (final) infiltration rate.
#' This modification ensures that infiltration asymptotically approaches
#' a non-zero steady-state rate, unlike the original Kostiakov model.
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
#'       (\eqn{k}, \eqn{a}, \eqn{f_c}).
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
#' using heuristic rules to improve convergence stability. The equilibrium
#' infiltration rate \eqn{f_c} is initialized from late-time observations,
#' while \eqn{k} and \eqn{a} are estimated via log-linear regression.
#' Parameter bounds are constrained to non-negative values.
#'
#' The function computes standard hydrological performance metrics
#' including Nash–Sutcliffe Efficiency (NSE), Percent Bias (PBIAS),
#' Kling–Gupta Efficiency (KGE), and Normalized Root Mean Square Error (NRMSE).
#'
#' @references
#' Kostiakov, A. N. (1932). The dynamics of the coefficient of water
#' percolation in soils and the necessity for studying it from a dynamic
#' point of view for the purpose of amelioration. *Society of Soil Science*, Moscow.
#'
#' Lewis, M. R. (1949). The infiltration equation and its application.
#' *Transactions of the American Geophysical Union*, 30(3), 459–470.
#'
#' USDA Soil Conservation Service (1956). National Engineering Handbook,
#' Section 4: Hydrology. U.S. Department of Agriculture, Washington, D.C.
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
#' # Simulated infiltration curve with Modified Kostiakov form
#' k <- 20; a <- 0.7; fc <- 5
#' rate <- k * a * time^(a - 1) + fc
#' rate <- rate + rnorm(length(rate), 0, 0.5)  # add noise
#'
#' data <- data.frame(time = time, rate = rate)
#'
#' model <- fit_modified_kostiakov(data, time = time, rate = rate)
#'
#' model$params
#' model$stats
fit_modified_kostiakov <- function(data = NULL, time, rate, na.rm = TRUE) {

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
  # Modified Kostiakov equation
  # f(t) = k * a * t^(a - 1) + fc
  # ---------------------------
  mk_fun <- function(t, k, a, fc) {
    k * a * t^(a - 1) + fc
  }

  # ---------------------------
  # Smart starting values
  # ---------------------------
  fc_start <- mean(tail(infil, max(3, floor(length(infil) * 0.2))))

  # Try log-log regression for a and k
  log_t <- log(time)
  log_f <- log(pmax(infil - fc_start, 1e-6))
  lm_fit <- try(stats::lm(log_f ~ log_t), silent = TRUE)

  if (!inherits(lm_fit, "try-error")) {
    a_start <- max(min(coef(lm_fit)[2] + 1, 1), 0.1)  # keep between 0.1 and 1
    k_start <- max(exp(coef(lm_fit)[1]) / a_start, 1e-6)
  } else {
    # Fallback if regression fails
    a_start <- 0.5
    k_start <- mean(infil) / mean(time)
  }

  fc_start <- max(fc_start, 1e-6)

  # ---------------------------
  # Nonlinear fit
  # ---------------------------
  fit <- minpack.lm::nlsLM(
    infil ~ mk_fun(time, k, a, fc),
    start = list(k = k_start, a = a_start, fc = fc_start),
    lower = c(k = 0, a = 0, fc = 0),
    control = minpack.lm::nls.lm.control(maxiter = 200)
  )

  # ---------------------------
  # Predictions & residuals
  # ---------------------------
  pred  <- stats::predict(fit)
  resid <- infil - pred

  # ---------------------------
  # Parameters
  # ---------------------------
  k_fit  <- coef(fit)[["k"]]
  a_fit  <- coef(fit)[["a"]]
  fc_fit <- coef(fit)[["fc"]]

  # ---------------------------
  # Smooth predictions on fine grid
  # ---------------------------
  new_time <- seq(min(time), max(time), length.out = 200)
  new_rate <- mk_fun(new_time, k_fit, a_fit, fc_fit)
  new_cum  <- k_fit * new_time^a_fit + fc_fit * new_time

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
    model = "Modified Kostiakov",
    params = data.frame(k = k_fit, a = a_fit, fc = fc_fit),
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
