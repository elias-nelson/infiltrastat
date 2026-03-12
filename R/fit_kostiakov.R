#' Fit the Kostiakov Infiltration Model
#'
#' Fits the empirical Kostiakov infiltration model to observed
#' infiltration rate data using nonlinear least squares with the
#' Levenberg–Marquardt algorithm.
#'
#' The Kostiakov cumulative infiltration function is:
#'
#' \deqn{F(t) = k \cdot t^a}
#'
#' where:
#' \itemize{
#'   \item \eqn{F(t)} — cumulative infiltration,
#'   \item \eqn{k} — empirical scaling coefficient,
#'   \item \eqn{a} — empirical exponent (typically 0 < a < 1),
#'   \item \eqn{t} — time.
#' }
#'
#' The instantaneous infiltration rate is derived as:
#'
#' \deqn{f(t) = k \cdot a \cdot t^{a-1}}
#'
#' @param data A data frame containing infiltration observations.
#' @param time Unquoted column name representing cumulative time
#'   (numeric, non-negative; typically minutes).
#' @param rate Unquoted column name representing observed infiltration rate
#'   (mm/hr or cm/hr).
#' @param na.rm Logical. If \code{TRUE} (default), incomplete observations
#'   are removed before fitting.
#'
#' @returns A list of class \code{"kostiakov_model"} with components:
#' \itemize{
#'   \item \code{model} — Model name.
#'   \item \code{params} — Estimated parameters (\eqn{k}, \eqn{a}).
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
#' The Kostiakov model is empirical and widely used for infiltration studies.
#' It assumes cumulative infiltration follows a power-law function of time.
#' Instantaneous infiltration rate is obtained by differentiation.
#' Starting values are estimated via log–log regression of observed rate
#' against time. At least five observations are required.
#'
#' @references
#' Kostiakov, A. N. (1932). On the dynamics of the coefficient of water-percolation
#' in soils and on the necessity for studying it from a dynamic point of view
#' for purposes of amelioration. *Transactions of the Sixth Commission, International
#' Society of Soil Science*, Part A, 17–21.
#'
#' @export
#' @importFrom stats lm predict
#' @examples
#' # Example dataset
#' infilt_data <- data.frame(
#'   Time = c(5, 10, 15, 20, 30, 45, 60),
#'   Rate = c(25, 18, 14, 12, 9, 7, 6)  # mm/hr
#' )
#'
#' # Fit the Kostiakov model
#' kostiakov_fit <- fit_kostiakov(data = infilt_data, time = Time, rate = Rate)
#'
#' # Inspect estimated parameters
#' kostiakov_fit$params
#'
#' # Plot observed vs fitted infiltration rate
#' with(kostiakov_fit$fitted, {
#'   plot(time, observed, pch = 16, col = "blue",
#'        xlab = "Time (minutes)", ylab = "Infiltration rate (mm/hr)")
#'   lines(time, predicted, col = "red", lwd = 2)
#'   legend("topright", legend = c("Observed", "Fitted"),
#'          col = c("blue", "red"), pch = c(16, NA), lty = c(NA, 1))
#' })
#'
#' # Plot fitted cumulative infiltration curve
#' with(kostiakov_fit$fitted, {
#'   plot(time, cumul_fit, type = "l", col = "red", lwd = 2,
#'        xlab = "Time (minutes)", ylab = "Cumulative infiltration (mm)")
#' })
fit_kostiakov <- function(data, time, rate, na.rm = TRUE) {

  # ---------------------------
  # Extract columns
  # ---------------------------
  time  <- data[[deparse(substitute(time))]]
  infil <- data[[deparse(substitute(rate))]]

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
  # Kostiakov functions
  # ---------------------------
  kostiakov_fun <- function(t, k, a) {
    k * a * t^(a - 1)   # instantaneous rate
  }

  kostiakov_cum <- function(t, k, a) {
    k * t^a            # cumulative infiltration
  }

  # ---------------------------
  # Starting values (log-log regression)
  # ---------------------------
  log_t <- log(time)
  log_f <- log(infil)
  lm_fit <- stats::lm(log_f ~ log_t)
  a_start <- coef(lm_fit)[2] + 1
  k_start <- exp(coef(lm_fit)[1]) / a_start

  k_start <- max(k_start, 1e-6)
  a_start <- max(min(a_start, 1), 1e-6)

  # ---------------------------
  # Nonlinear fit
  # ---------------------------
  fit <- minpack.lm::nlsLM(
    infil ~ kostiakov_fun(time, k, a),
    start = list(k = k_start, a = a_start),
    lower = c(k = 0, a = 0),
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
  k_fit <- coef(fit)[["k"]]
  a_fit <- coef(fit)[["a"]]
  fitted_cum <- kostiakov_cum(time, k_fit, a_fit)

  # ---------------------------
  # Smooth predictions on fine grid
  # ---------------------------
  new_time <- seq(min(time), max(time), length.out = 200)
  new_rate <- kostiakov_fun(new_time, k_fit, a_fit)
  new_cum  <- kostiakov_cum(new_time, k_fit, a_fit)

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
    model = "Kostiakov",
    params = data.frame(k = k_fit, a = a_fit),
    fitted = data.frame(
      time      = time,
      observed  = infil,
      predicted = pred,
      residual  = resid,
      cumul_fit = fitted_cum
    ),
    smooth = smooth_pred,   # NEW: fine-grid predictions
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
