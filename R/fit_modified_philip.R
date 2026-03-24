#' Fit the Modified Philip (Philip–Kostiakov Hybrid) Infiltration Model
#'
#' Fits the Modified Philip infiltration model to observed
#' infiltration rate data using nonlinear least squares with the
#' Levenberg–Marquardt algorithm.
#'
#' The Modified Philip cumulative infiltration function is:
#'
#' \deqn{F(t) = S \cdot t^{1/2} + K \cdot t + \alpha \cdot t^b}
#'
#' where:
#' \itemize{
#'   \item \eqn{F(t)} — cumulative infiltration,
#'   \item \eqn{S} — sorptivity,
#'   \item \eqn{K} — steady infiltration term,
#'   \item \eqn{\alpha} — empirical coefficient,
#'   \item \eqn{b} — empirical exponent (typically 0 < b < 1),
#'   \item \eqn{t} — time.
#' }
#'
#' The instantaneous infiltration rate is derived as:
#'
#' \deqn{f(t) = \frac{S}{2 \sqrt{t}} + K + \alpha \cdot b \cdot t^{b-1}}
#'
#' @param data A data frame containing infiltration observations.
#' @param time Unquoted column name representing cumulative time
#'   (numeric, non-negative; typically minutes).
#' @param rate Unquoted column name representing observed infiltration rate
#'   (mm/hr or cm/hr).
#' @param na.rm Logical. If \code{TRUE} (default), incomplete observations
#'   are removed before fitting.
#'
#' @returns A list of class \code{"modified_philip_model"} with components:
#' \itemize{
#'   \item \code{model} — Model name.
#'   \item \code{params} — Estimated parameters (\eqn{S}, \eqn{K}, \eqn{\alpha}, \eqn{b}).
#'   \item \code{fitted} — Data frame with:
#'     \itemize{
#'       \item \code{time} — elapsed time,
#'       \item \code{observed} — measured infiltration rate,
#'       \item \code{predicted} — fitted infiltration rate,
#'       \item \code{residual} — residuals (observed – fitted rate),
#'       \item \code{cumul_fit} — fitted cumulative infiltration.
#'     }
#'   \item \code{stats} — Goodness-of-fit statistics (RMSE, MAE, NRMSE,
#'   PBIAS, R², NSE, AIC).
#'   \item \code{nls_obj} — The fitted \code{nlsLM} model object.
#' }
#'
#' @details
#' The Modified Philip model extends the original Philip formulation by
#' adding an empirical Kostiakov-style term (\eqn{\alpha \cdot t^b}),
#' providing more flexibility in fitting infiltration curves. This hybrid
#' approach is useful when neither Philip nor Kostiakov alone adequately
#' describe observed infiltration behavior.
#'
#' Starting values are estimated internally: \eqn{K} from late-time rates,
#' \eqn{S} from early-time variation, and default values for \eqn{\alpha}
#' and \eqn{b}.
#'
#' @references
#' Dunne, T., & Philip, J. R. (1979). Infiltration models for hydrologic
#' prediction. *Water Resources Research*, 15(2), 161–170.
#'
#' @export
#' @importFrom stats predict
#' @importFrom utils head
#'
#' @examples
#' # Example dataset
#' infilt_data <- data.frame(
#'   Time = c(5, 10, 15, 20, 30, 45, 60),
#'   Rate = c(32, 22, 16, 13, 10, 8, 7)  # mm/hr
#' )
#'
#' # Fit the Modified Philip model
#' mp_fit <- fit_modified_philip(data = infilt_data, time = Time, rate = Rate)
#'
#' # Inspect estimated parameters
#' mp_fit$params
#'
#' # Plot observed vs fitted infiltration rate
#' with(mp_fit$fitted, {
#'   plot(time, observed, pch = 16, col = "blue",
#'        xlab = "Time (minutes)", ylab = "Infiltration rate (mm/hr)")
#'   lines(time, predicted, col = "red", lwd = 2)
#'   legend("topright", legend = c("Observed", "Fitted"),
#'          col = c("blue", "red"), pch = c(16, NA), lty = c(NA, 1))
#' })
fit_modified_philip <- function(data = NULL, time, rate, na.rm = TRUE) {

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
  # Modified Philip function
  # f(t) = (S / (2√t)) + K + α b t^(b-1)
  # ---------------------------
  mp_fun <- function(t, S, K, alpha, b) {
    (S / (2 * sqrt(t))) + K + alpha * b * t^(b - 1)
  }

  # ---------------------------
  # Smart starting values
  # ---------------------------
  f0_start <- infil[1]
  fc_start <- mean(tail(infil, max(3, floor(length(infil) * 0.2))))

  # Sorptivity estimate from early decline
  S_start <- (f0_start - fc_start) * sqrt(mean(head(time, 3)))

  # K from late infiltration
  K_start <- fc_start

  # Kostiakov-style regression for alpha and b
  log_t <- log(time)
  log_f <- log(infil)
  lm_fit <- try(lm(log_f ~ log_t), silent = TRUE)

  if (!inherits(lm_fit, "try-error")) {
    b_start <- max(min(-coef(lm_fit)[2], 1), 0.1)  # keep between 0.1 and 1
    alpha_start <- max(exp(coef(lm_fit)[1]), 1e-6)
  } else {
    b_start <- 0.5
    alpha_start <- 0.1
  }

  # Ensure positivity
  S_start <- max(S_start, 1e-6)
  K_start <- max(K_start, 1e-6)

  # ---------------------------
  # Nonlinear fit
  # ---------------------------
  fit <- minpack.lm::nlsLM(
    infil ~ mp_fun(time, S, K, alpha, b),
    start = list(S = S_start, K = K_start, alpha = alpha_start, b = b_start),
    lower = c(S = 0, K = 0, alpha = 0, b = 0),
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
  S_fit     <- coef(fit)[["S"]]
  K_fit     <- coef(fit)[["K"]]
  alpha_fit <- coef(fit)[["alpha"]]
  b_fit     <- coef(fit)[["b"]]

  # ---------------------------
  # Smooth predictions on fine grid
  # ---------------------------
  new_time <- seq(min(time), max(time), length.out = 200)
  new_rate <- mp_fun(new_time, S_fit, K_fit, alpha_fit, b_fit)
  new_cum  <- S_fit * sqrt(new_time) + K_fit * new_time + alpha_fit * new_time^b_fit

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
    model = "Modified Philip",
    params = data.frame(S = S_fit, K = K_fit, alpha = alpha_fit, b = b_fit),
    fitted = data.frame(
      time      = time,
      observed  = infil,
      predicted = pred,
      residual  = resid
    ),
    smooth = smooth_pred,
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
