#' Fit the Modified Green–Ampt Infiltration Model
#'
#' Fits the Modified Green–Ampt infiltration model to observed
#' infiltration data using nonlinear least squares with the
#' Levenberg–Marquardt algorithm. This variant introduces an
#' empirical adjustment parameter (\eqn{\alpha}) to improve
#' flexibility under non-ideal soil conditions.
#' This implementation follows the Modified Green–Ampt approach
#' as described by Mein & Larson (1973) and Rawls et al. (1983).
#'
#' The Modified Green–Ampt cumulative infiltration equation is defined implicitly as:
#'
#' \deqn{F(t) - K t - \alpha N_s \log\left(1 + \frac{F(t)}{N_s}\right) = 0}
#'
#' where:
#' \itemize{
#'   \item \eqn{F(t)} is cumulative infiltration,
#'   \item \eqn{K} is the saturated hydraulic conductivity,
#'   \item \eqn{N_s} is the effective suction–moisture deficit term,
#'   \item \eqn{\alpha} is an empirical adjustment parameter,
#'   \item \eqn{t} is time.
#' }
#'
#' The instantaneous infiltration rate is derived from the fitted parameters as:
#'
#' \deqn{f(t) = K \left(1 + \frac{\alpha N_s}{F(t)}\right)}
#'
#' @param data A data frame containing infiltration observations.
#' @param time Unquoted column name representing cumulative time
#'   (numeric, non-negative; typically minutes).
#' @param rate Unquoted column name representing observed infiltration rate
#'   (mm/hr or cm/hr).
#' @param na.rm Logical. If \code{TRUE} (default), incomplete observations
#'   are removed before fitting.
#'
#' @returns A list of class \code{"modified_greenampt_model"} with components:
#' \itemize{
#'   \item \code{model} — Model name.
#'   \item \code{params} — Estimated parameters (\eqn{K}, \eqn{N_s}, \eqn{\alpha}).
#'   \item \code{fitted} — Data frame with:
#'     \itemize{
#'       \item \code{time_min}, \code{time_hr} — elapsed time,
#'       \item \code{observed_rate} — measured infiltration rate,
#'       \item \code{predicted_rate} — rate curve implied by fitted parameters,
#'       \item \code{cumul_obs} — observed cumulative infiltration,
#'       \item \code{predicted_cum} — fitted cumulative infiltration,
#'       \item \code{residual_cum} — residuals (observed – fitted cumulative).
#'     }
#'   \item \code{smooth} — Data frame of fine-grid predictions for plotting
#'     (predicted cumulative infiltration and infiltration rate).
#'   \item \code{stats} — Goodness-of-fit statistics (RMSE, MAE, NRMSE,
#'   PBIAS, R², NSE, KGE, AIC).
#'   \item \code{nls_obj} — The fitted \code{nlsLM} model object.
#' }
#'
#' @details
#' The model solves the implicit Modified Green–Ampt equation for cumulative
#' infiltration at each time step using \code{uniroot()}. The fitted rate curve
#' is then derived from the estimated parameters. Starting values are estimated
#' internally: \eqn{K} from late-time observed rates, \eqn{N_s} relative to mean
#' cumulative infiltration, and \eqn{\alpha} initialized to 1. At least five
#' observations are required.
#'
#' @references
#' Green, W. H., & Ampt, G. A. (1911). Studies on soil physics:
#' The flow of air and water through soils. *Journal of Agricultural Science*, 4, 1–24.
#'
#' Rawls, W. J., Brakensiek, D. L., & Miller, N. (1983).
#' Green-Ampt infiltration parameters from soil data.
#' *Journal of Hydraulic Engineering*, 109(1), 62–70.
#'
#' Mein, R. G., & Larson, C. L. (1973).
#' Modeling infiltration during a steady rain.
#' *Water Resources Research*, 9(2), 384–394.
#'
#' @export
#' @importFrom stats uniroot predict AIC coef cor
#' @importFrom minpack.lm nlsLM nls.lm.control
#' @importFrom utils tail
#'
#' @examples
#' # Example using synthetic data
#' set.seed(123)
#' time <- seq(5, 60, by = 5)
#'
#' # Simulated infiltration curve
#' rate <- 10 * (1 + 50 / (time + 10))
#' rate <- rate + rnorm(length(rate), 0, 0.5)
#'
#' data <- data.frame(time = time, rate = rate)
#'
#' model <- fit_modified_greenampt(data, time = time, rate = rate)
#'
#' model$params
#' model$stats
fit_modified_greenampt <- function(data = NULL, time, rate, na.rm = TRUE) {

  # ---------------------------
  # Extract data columns
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
  # Convert time to hours (assuming input is minutes)
  # ---------------------------
  time_hr <- time / 60

  # ---------------------------
  # Compute incremental dt and cumulative infiltration
  # ---------------------------
  dt_hr <- diff(c(0, time_hr))
  cumul_obs <- cumsum(infil * dt_hr)

  # ---------------------------
  # Modified Green-Ampt functions
  # Adds empirical adjustment parameter 'alpha'
  # ---------------------------
  MGA_F <- function(t, Ks, PsiDT, alpha) {
    f_root <- function(F) F - Ks * t - alpha * PsiDT * log(1 + F / PsiDT)
    upper_bound <- max(cumul_obs, na.rm = TRUE) * 10
    uniroot(f_root, interval = c(0, upper_bound), tol = 1e-8)$root
  }

  MGA_rate <- function(t, Ks, PsiDT, alpha) {
    F_t <- vapply(t, MGA_F, FUN.VALUE = numeric(1), Ks = Ks, PsiDT = PsiDT, alpha = alpha)
    Ks * (1 + (alpha * PsiDT) / F_t)
  }

  # ---------------------------
  # Starting values
  # ---------------------------
  Ks_start    <- mean(tail(infil, max(3, floor(length(infil) * 0.1))))
  PsiDT_start <- mean(cumul_obs) * 0.5
  alpha_start <- 1

  Ks_start    <- max(Ks_start, 1e-6)
  PsiDT_start <- max(PsiDT_start, 1e-6)
  alpha_start <- max(alpha_start, 1e-6)

  # ---------------------------
  # Nonlinear fit via nlsLM
  # ---------------------------
  fit <- minpack.lm::nlsLM(
    cumul_obs ~ vapply(time_hr, MGA_F, FUN.VALUE = numeric(1),
                       Ks = Ks, PsiDT = PsiDT, alpha = alpha),
    start = list(Ks = Ks_start, PsiDT = PsiDT_start, alpha = alpha_start),
    lower = c(Ks = 0, PsiDT = 0, alpha = 0),
    control = minpack.lm::nls.lm.control(maxiter = 500)
  )

  # ---------------------------
  # Predictions and residuals
  # ---------------------------
  pred <- stats::predict(fit)
  resid <- cumul_obs - pred

  # ---------------------------
  # Derive fitted rate curve from parameters
  # ---------------------------
  Ks_fit    <- coef(fit)[["Ks"]]
  PsiDT_fit <- coef(fit)[["PsiDT"]]
  alpha_fit <- coef(fit)[["alpha"]]
  fitted_rate <- MGA_rate(time_hr, Ks_fit, PsiDT_fit, alpha_fit)

  # ---------------------------
  # Smooth predictions on fine grid
  # ---------------------------
  new_time_hr <- seq(min(time_hr), max(time_hr), length.out = 200)
  new_cum <- vapply(new_time_hr, MGA_F, FUN.VALUE = numeric(1),
                    Ks = Ks_fit, PsiDT = PsiDT_fit, alpha = alpha_fit)
  new_rate <- MGA_rate(new_time_hr, Ks_fit, PsiDT_fit, alpha_fit)

  smooth_pred <- data.frame(
    time_hr = new_time_hr,
    time_min = new_time_hr * 60,
    predicted_cum = new_cum,
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
    model = "Modified Green-Ampt",
    params = data.frame(Ks = Ks_fit, PsiDT = PsiDT_fit, alpha = alpha_fit),
    fitted = data.frame(
      time_min       = time,
      time_hr        = time_hr,
      observed_rate  = infil,
      predicted_rate = fitted_rate,
      cumul_obs      = cumul_obs,
      predicted_cum  = pred,
      residual_cum   = resid
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
