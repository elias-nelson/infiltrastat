#' Fit the Smith-Parlange Infiltration Model
#'
#' Fits the Smith-Parlange infiltration model to observed
#' infiltration rate data using nonlinear least squares with the
#' Levenberg-Marquardt algorithm.
#'
#' The Smith-Parlange infiltration rate function is:
#'
#' \deqn{f(t) = f_c + \frac{(f_0 - f_c)}{(1 + k \cdot t)^n}}
#'
#' where:
#' \itemize{
#'   \item \eqn{f(t)} - instantaneous infiltration rate,
#'   \item \eqn{f_0} - initial infiltration rate,
#'   \item \eqn{f_c} - steady infiltration rate,
#'   \item \eqn{k} - decay constant,
#'   \item \eqn{n} - shape parameter,
#'   \item \eqn{t} - time.
#' }
#'
#' The cumulative infiltration function is:
#'
#' \deqn{F(t) = f_c \cdot t + \frac{f_0 - f_c}{k (1-n)} \left[1 - (1 + k t)^{1-n}\right]}
#'
#' for \eqn{n \neq 1}. If \eqn{n = 1}, cumulative infiltration is
#' approximated numerically by integrating the fitted rate curve.
#'
#' @param data A data frame containing infiltration observations.
#' @param time Unquoted column name representing cumulative time
#'   (numeric, non-negative; typically minutes).
#' @param rate Unquoted column name representing observed infiltration rate
#'   (mm/hr or cm/hr).
#' @param na.rm Logical. If \code{TRUE} (default), incomplete observations
#'   are removed before fitting.
#'
#' @returns A list of class \code{"smith_parlange_model"} with components:
#' \itemize{
#'   \item \code{model} - Model name.
#'   \item \code{params} - Estimated parameters (\eqn{f_0}, \eqn{f_c}, \eqn{k}, \eqn{n}).
#'   \item \code{fitted} - Data frame with:
#'     \itemize{
#'       \item \code{time} - elapsed time,
#'       \item \code{observed} - measured infiltration rate,
#'       \item \code{predicted} - fitted infiltration rate,
#'       \item \code{residual} - residuals (observed - fitted rate),
#'       \item \code{cumul_fit} - fitted cumulative infiltration.
#'     }
#'   \item \code{stats} - Goodness-of-fit statistics (RMSE, MAE, NRMSE,
#'   PBIAS, R², NSE, AIC).
#'   \item \code{nls_obj} - The fitted \code{nlsLM} model object.
#' }
#'
#' @details
#' The Smith-Parlange model generalizes Green-Ampt and Horton forms by
#' introducing a flexible decay function with shape parameter \eqn{n}.
#' It is particularly useful for soils with variable ponding and
#' heterogeneity. Cumulative infiltration is computed analytically
#' when \eqn{n \neq 1}, and numerically otherwise.
#'
#' Starting values are estimated internally: \eqn{f_0} from the maximum
#' observed rate, \eqn{f_c} from late-time observations, and default
#' values for \eqn{k} and \eqn{n}.
#'
#' @references
#' Smith, R. E., & Parlange, J. Y. (1978). A parameter-efficient
#' infiltration model. *Water Resources Research*, 14(3), 533-538.
#'
#' @export
#' @importFrom stats predict
#'
#' @examples
#' # Example dataset
#' infilt_data <- data.frame(
#'   Time = c(5, 10, 15, 20, 30, 45, 60),
#'   Rate = c(35, 25, 18, 14, 11, 9, 8)  # mm/hr
#' )
#'
#' # Fit the Smith-Parlange model
#' sp_fit <- fit_smith_parlange(data = infilt_data, time = Time, rate = Rate)
#'
#' # Inspect estimated parameters
#' sp_fit$params
#'
#' # Plot observed vs fitted infiltration rate
#' with(sp_fit$fitted, {
#'   plot(time, observed, pch = 16, col = "blue",
#'        xlab = "Time (minutes)", ylab = "Infiltration rate (mm/hr)")
#'   lines(time, predicted, col = "red", lwd = 2)
#'   legend("topright", legend = c("Observed", "Fitted"),
#'          col = c("blue", "red"), pch = c(16, NA), lty = c(NA, 1))
#' })
fit_smith_parlange <- function(data = NULL, time, rate, na.rm = TRUE) {

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
  if (any(time < 0))
    stop("time must be >= 0.")

  # ---------------------------
  # Ensure increasing time
  # ---------------------------
  ord <- order(time)
  time  <- time[ord]
  infil <- infil[ord]

  # ---------------------------
  # Smith-Parlange function
  # f(t) = fc + (f0 - fc) / ((1 + k * t)^n)
  # ---------------------------
  sp_fun <- function(t, f0, fc, k, n) {
    fc + (f0 - fc) / ((1 + k * t)^n)
  }

  # ---------------------------
  # Starting values
  # ---------------------------
  f0_start <- max(infil)
  fc_start <- mean(tail(infil, max(3, floor(length(infil) * 0.1))))
  k_start  <- 0.01
  n_start  <- 1

  # ---------------------------
  # Nonlinear fit
  # ---------------------------
  fit <- minpack.lm::nlsLM(
    infil ~ sp_fun(time, f0, fc, k, n),
    start = list(f0 = f0_start, fc = fc_start, k = k_start, n = n_start),
    lower = c(f0 = 0, fc = 0, k = 0, n = 0),
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
  f0_fit <- coef(fit)[["f0"]]
  fc_fit <- coef(fit)[["fc"]]
  k_fit  <- coef(fit)[["k"]]
  n_fit  <- coef(fit)[["n"]]

  # ---------------------------
  # Smooth predictions on fine grid
  # ---------------------------
  new_time <- seq(min(time), max(time), length.out = 200)
  new_rate <- sp_fun(new_time, f0_fit, fc_fit, k_fit, n_fit)

  # Cumulative infiltration (analytical if n is not equal to 1, else numerical)
  if (abs(n_fit - 1) > 1e-6) {
    new_cum <- fc_fit * new_time + (f0_fit - fc_fit) / (k_fit * (1 - n_fit)) *
      (1 - (1 + k_fit * new_time)^(1 - n_fit))
  } else {
    dt <- diff(c(0, new_time))
    new_cum <- cumsum(new_rate * dt)
  }

  smooth_pred <- data.frame(
    time = new_time,
    predicted_rate = new_rate,
    predicted_cum  = new_cum
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
    model = "Smith-Parlange",
    params = data.frame(f0 = f0_fit, fc = fc_fit, k = k_fit, n = n_fit),
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
