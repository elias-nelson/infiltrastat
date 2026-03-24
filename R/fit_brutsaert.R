#' Fit the Brutsaert Infiltration Model
#'
#' This function fits the Brutsaert infiltration model to observed infiltration
#' rate data using nonlinear least squares optimization. The Brutsaert model
#' expresses infiltration rate as:
#'
#' \deqn{f(t) = \frac{S}{2 \sqrt{\pi t}} + K}
#'
#' where \eqn{S} is the sorptivity parameter and \eqn{K} is the saturated
#' hydraulic conductivity. The function estimates these parameters from field
#' data and returns fitted values, residuals, smooth predictions, and
#' performance statistics.
#'
#' @param data Optional data frame containing time and infiltration rate columns.
#' @param time Numeric vector or column name in `data` representing time (minutes).
#' @param rate Numeric vector or column name in `data` representing infiltration rate.
#' @param na.rm Logical; if TRUE, missing values are removed before fitting.
#'
#' @details
#' The fitting procedure uses the Levenberg–Marquardt algorithm implemented in
#' \code{minpack.lm::nlsLM}. Starting values for \eqn{S} and \eqn{K} are
#' automatically generated from early and late portions of the infiltration curve.
#'
#' The output includes:
#' \itemize{
#'   \item Estimated parameters (\eqn{S}, \eqn{K})
#'   \item Fitted vs. observed infiltration rates
#'   \item Smooth predictions on a fine time grid
#'   \item Performance statistics (RMSE, MAE, NSE, R², KGE, AIC)
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{model}{Model name ("Brutsaert")}
#'   \item{params}{Data frame of fitted parameters}
#'   \item{fitted}{Data frame of observed, predicted, and residual values}
#'   \item{smooth}{Data frame of smooth predictions (rate and cumulative infiltration)}
#'   \item{stats}{Data frame of performance metrics}
#'   \item{nls_obj}{Underlying \code{nlsLM} object}
#' }
#'
#' @references
#' Brutsaert, W. (1977). Vertical infiltration in dry soil. \emph{Water Resources Research}, 13(2), 363–368.
#'
#' @examples
#' # Example dataset
#' time <- c(5, 10, 15, 20, 30, 45, 60)
#' rate <- c(30, 20, 15, 12, 9, 7, 6)  # mm/hr
#' data <- data.frame(time = time, rate = rate)
#'
#' # Fit Brutsaert model
#' model <- fit_brutsaert(data, time = time, rate = rate)
#' model$params
#' plot(data$time, data$rate, pch = 16, xlab = "Time (min)", ylab = "Infiltration rate (mm/hr)")
#' lines(model$smooth$time, model$smooth$predicted_rate, col = "blue")
#'
#' @export
fit_brutsaert <- function(data = NULL, time, rate, na.rm = TRUE) {

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
  # Brutsaert equation
  # f(t) = S / (2 * sqrt(pi * t)) + K
  # ---------------------------
  brut_fun <- function(t, S, K) {
    S / (2 * sqrt(pi * t)) + K
  }

  # ---------------------------
  # Starting values
  # ---------------------------
  K_start <- mean(tail(infil, max(3, floor(length(infil) * 0.1))))
  S_start <- diff(range(infil)) * sqrt(mean(time))

  K_start <- max(K_start, 1e-6)
  S_start <- max(S_start, 1e-6)

  # ---------------------------
  # Nonlinear fit
  # ---------------------------
  fit <- minpack.lm::nlsLM(
    infil ~ brut_fun(time, S, K),
    start = list(S = S_start, K = K_start),
    lower = c(S = 0, K = 0),
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
  S_fit <- coef(fit)[["S"]]
  K_fit <- coef(fit)[["K"]]

  # ---------------------------
  # Smooth predictions on fine grid
  # ---------------------------
  new_time <- seq(min(time), max(time), length.out = 200)
  new_rate <- brut_fun(new_time, S_fit, K_fit)
  new_cum  <- S_fit * sqrt(new_time) + K_fit * new_time   # cumulative infiltration

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
    model = "Brutsaert",
    params = data.frame(S = S_fit, K = K_fit),
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
