#' Fit the Green-Ampt Infiltration Model
#'
#' Fits the physically based Green–Ampt infiltration rate model to observed
#' infiltration data using nonlinear least squares with the
#' Levenberg–Marquardt algorithm.
#'
#' The Green–Ampt cumulative infiltration equation is defined implicitly as:
#'
#' \deqn{F - K t - N_s \log\left(1 + \frac{F}{N_s}\right) = 0}
#'
#' where:
#' \itemize{
#'   \item \eqn{F} is cumulative infiltration,
#'   \item \eqn{K} is the hydraulic conductivity parameter,
#'   \item \eqn{N_s} is the effective capillary suction term,
#'   \item \eqn{t} is time.
#' }
#'
#' The instantaneous infiltration rate is computed as:
#'
#' \deqn{f(t) = K \left(1 + \frac{N_s}{F(t)}\right)}
#'
#' where \eqn{F(t)} is obtained numerically via root-finding.
#'
#' @param data A data frame containing infiltration observations.
#' @param time Unquoted column name representing cumulative time
#'   (must be numeric and >= 0).
#' @param rate Unquoted column name representing observed infiltration rate.
#' @param na.rm Logical. If \code{TRUE} (default), incomplete observations
#'   are removed before fitting.
#'
#' @returns A list with class \code{"greenampt_model"} containing:
#' \itemize{
#'   \item \code{model} — Model name.
#'   \item \code{params} — Estimated parameters (\eqn{K}, \eqn{N_s}).
#'   \item \code{fitted} — Data frame with observed, predicted,
#'   residual, and cumulative observed infiltration.
#'   \item \code{stats} — Goodness-of-fit statistics (RMSE, MAE, NRMSE,
#'   PBIAS, R², NSE, AIC).
#'   \item \code{nls_obj} — The fitted \code{nlsLM} model object.
#' }
#'
#' @details
#' The Green–Ampt model is physically based and requires solving an
#' implicit nonlinear equation for cumulative infiltration at each time step.
#' This implementation computes \eqn{F(t)} numerically using
#' \code{uniroot()} within each iteration of the nonlinear regression.
#'
#' Starting values are estimated internally:
#' \itemize{
#'   \item \eqn{K} is approximated from the late-time observed infiltration rate.
#'   \item \eqn{N_s} is initialized relative to mean cumulative infiltration.
#' }
#'
#' A minimum of five observations is required for fitting.
#'
#' @references
#' Green, W. H., & Ampt, G. A. (1911). Studies on soil physics:
#' The flow of air and water through soils. *Journal of Agricultural Science*, 4, 1–24.
#'
#' @export
#'
#' @importFrom stats uniroot
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
#' model <- fit_greenampt(data, time = time, rate = rate)
#'
#' model$params
#' model$stats
fit_greenampt <- function(data, time, rate, na.rm = TRUE) {

  # ---------------------------
  # Extract data columns
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
  # Model definition
  # Implicit F(t) solved by root-finding
  # ---------------------------
  GA_F <- function(t, Ks, PsiDT) {
    f_root <- function(F) F - Ks * t - PsiDT * log(1 + F / PsiDT)
    upper_bound <- max(cumul_obs, na.rm = TRUE) * 10
    uniroot(f_root, interval = c(0, upper_bound), tol = 1e-8)$root
  }

  GA_rate <- function(t, Ks, PsiDT) {
    F_t <- vapply(t, GA_F, FUN.VALUE = numeric(1), Ks = Ks, PsiDT = PsiDT)
    Ks * (1 + PsiDT / F_t)
  }

  # ---------------------------
  # Starting values
  # ---------------------------
  Ks_start    <- mean(tail(infil, max(3, floor(length(infil) * 0.1))))
  PsiDT_start <- mean(cumul_obs) * 0.5

  Ks_start    <- max(Ks_start, 1e-6)
  PsiDT_start <- max(PsiDT_start, 1e-6)

  # ---------------------------
  # Nonlinear fit via nlsLM (fit cumulative infiltration)
  # ---------------------------
  fit <- minpack.lm::nlsLM(
    cumul_obs ~ vapply(time_hr, GA_F, FUN.VALUE = numeric(1), Ks = Ks, PsiDT = PsiDT),
    start = list(Ks = Ks_start, PsiDT = PsiDT_start),
    lower = c(Ks = 0, PsiDT = 0),
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
  fitted_rate <- GA_rate(time_hr, Ks_fit, PsiDT_fit)

  # ---------------------------
  # Goodness-of-fit metrics
  # ---------------------------
  rmse  <- sqrt(mean(resid^2))
  mae   <- mean(abs(resid))
  nrmse <- rmse / (max(cumul_obs) - min(cumul_obs)) * 100
  pbias <- 100 * sum(resid) / sum(cumul_obs)

  r2  <- 1 - sum(resid^2) / sum((cumul_obs - mean(cumul_obs))^2)
  nse <- 1 - sum(resid^2) / sum((cumul_obs - mean(cumul_obs))^2)
  aic <- stats::AIC(fit)

  # ---------------------------
  # Structured tidy output
  # ---------------------------
  list(
    model = "Green-Ampt",

    params = data.frame(
      Ks     = Ks_fit,
      PsiDT  = PsiDT_fit
    ),

    fitted = data.frame(
      time_min       = time,
      time_hr        = time_hr,
      observed_rate  = infil,
      predicted_rate = fitted_rate,
      cumul_obs      = cumul_obs,
      predicted_cum  = pred,
      residual_cum   = resid
    ),

    stats = data.frame(
      RMSE  = rmse,
      MAE   = mae,
      NRMSE = nrmse,
      PBIAS = pbias,
      R2    = r2,
      NSE   = nse,
      AIC   = aic
    ),

    nls_obj = fit
  )
}
