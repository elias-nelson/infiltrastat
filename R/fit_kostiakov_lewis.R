#' Fit the Kostiakov-Lewis Infiltration Model
#'
#' Fits the empirical Kostiakov-Lewis infiltration model to observed
#' infiltration rate data using nonlinear least squares with the
#' Levenberg-Marquardt algorithm.
#'
#' The Kostiakov-Lewis cumulative infiltration function is:
#'
#' \deqn{F(t) = k \cdot t^a + f_c \cdot t}
#'
#' where:
#' \itemize{
#'   \item \eqn{F(t)} - cumulative infiltration,
#'   \item \eqn{k} - empirical scaling coefficient,
#'   \item \eqn{a} - empirical exponent (typically 0 < a < 1),
#'   \item \eqn{f_c} - steady infiltration rate,
#'   \item \eqn{t} - time.
#' }
#'
#' The instantaneous infiltration rate is derived as:
#'
#' \deqn{f(t) = k \cdot a \cdot t^{a-1} + f_c}
#'
#' @param data A data frame containing infiltration observations.
#' @param time Unquoted column name representing cumulative time
#'   (numeric, non-negative; typically minutes).
#' @param rate Unquoted column name representing observed infiltration rate
#'   (mm/hr or cm/hr).
#' @param na.rm Logical. If \code{TRUE} (default), incomplete observations
#'   are removed before fitting.
#'
#' @returns A list of class \code{"kostiakov_lewis_model"} with components:
#' \itemize{
#'   \item \code{model} - Model name.
#'   \item \code{params} - Estimated parameters (\eqn{k}, \eqn{a}, \eqn{f_c}).
#'   \item \code{fitted} - Data frame with:
#'     \itemize{
#'       \item \code{time} - elapsed time,
#'       \item \code{observed} - measured infiltration rate,
#'       \item \code{predicted} - fitted infiltration rate,
#'       \item \code{residual} - residuals (observed – fitted rate),
#'       \item \code{cumul_fit} - fitted cumulative infiltration.
#'     }
#'   \item \code{stats} - Goodness-of-fit statistics (RMSE, MAE, NRMSE,
#'   PBIAS, R², NSE, AIC).
#'   \item \code{nls_obj} - The fitted \code{nlsLM} model object.
#' }
#'
#' @details
#' The Kostiakov-Lewis model extends the original Kostiakov formulation by
#' adding a steady infiltration term (\eqn{f_c}), allowing the rate to
#' approach a constant value at long times. Starting values are estimated
#' internally: \eqn{f_c} from late-time observed rates, and \eqn{k}, \eqn{a}
#' via log-log regression. At least five observations are required.
#'
#' @export
#' @importFrom stats cor sd
#'
#' @examples
#' # Example using synthetic data
#' set.seed(123)
#' time <- seq(5, 60, by = 5)
#' k <- 20; a <- 0.7; fc <- 5
#' rate <- k * a * time^(a - 1) + fc
#' rate <- rate + rnorm(length(rate), 0, 0.5)  # add noise
#'
#' data <- data.frame(Time = time, Rate = rate)
#'
#' model <- fit_kostiakov_lewis(data, time = Time, rate = Rate)
#' model$params
#' model$stats
fit_kostiakov_lewis <- function(data = NULL, time, rate, na.rm = TRUE) {

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

  if (any(time <= 0))
    stop("time must be > 0.")

  # ---------------------------
  # Ensure time is increasing
  # ---------------------------
  ord <- order(time)
  time  <- time[ord]
  infil <- infil[ord]

  # ---------------------------
  # Kostiakov-Lewis equation
  # f(t) = k * a * t^(a - 1) + fc
  # ---------------------------
  kl_fun <- function(t, k, a, fc) {
    k * a * t^(a - 1) + fc
  }

  # ---------------------------
  # Smart starting values
  # ---------------------------
  # steady infiltration from late-time values
  fc_start <- mean(tail(infil, max(3, floor(length(infil) * 0.1))))
  fc_start <- max(min(fc_start, min(infil)), 1e-6)

  # log-linear regression for k and a
  log_t <- log(time)
  log_f <- log(pmax(infil - fc_start, 1e-3))  # avoid negatives
  lm_fit <- stats::lm(log_f ~ log_t)

  a_start <- coef(lm_fit)[2] + 1
  a_start <- max(min(a_start, 0.9), 0.1)      # keep in (0.1, 0.9)

  k_start <- exp(coef(lm_fit)[1]) / a_start
  k_start <- max(k_start, 1e-6)

  # ---------------------------
  # Nonlinear fit
  # ---------------------------
  fit <- minpack.lm::nlsLM(
    infil ~ kl_fun(time, k, a, fc),
    start = list(k = k_start, a = a_start, fc = fc_start),
    lower = c(k = 0, a = 0, fc = 0),
    upper = c(k = Inf, a = 2, fc = Inf),   # allow a up to 2
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
  k_fit  <- coef(fit)[["k"]]
  a_fit  <- coef(fit)[["a"]]
  fc_fit <- coef(fit)[["fc"]]

  # ---------------------------
  # Smooth predictions on fine grid
  # ---------------------------
  new_time <- seq(min(time), max(time), length.out = 200)
  new_rate <- kl_fun(new_time, k_fit, a_fit, fc_fit)
  new_cum  <- k_fit * new_time^a_fit + fc_fit * new_time

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
    model = "Kostiakov-Lewis",
    params = data.frame(k = k_fit, a = a_fit, fc = fc_fit),
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
