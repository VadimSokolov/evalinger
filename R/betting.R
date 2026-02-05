#' Betting Martingale E-Process for Two-Arm Binary Trials
#'
#' Construct an e-process from paired binary outcomes using the betting
#' martingale \eqn{E_n = \prod_{i=1}^n (1 + \lambda_i D_i)}, where
#' \eqn{D_i = X_i^T - X_i^C}.
#'
#' @param x_T Integer vector of treatment arm binary outcomes (0/1).
#' @param x_C Integer vector of control arm binary outcomes (0/1), same length as \code{x_T}.
#' @param lambda Either a single fixed betting fraction in (0, 1), or a numeric
#'   vector of the same length as \code{x_T} giving predictable (pre-determined)
#'   betting fractions. Default is NULL, which uses the GROW-optimal fraction
#'   for the design alternative specified by \code{p_T_design} and \code{p_C_design}.
#' @param p_T_design Design alternative for treatment arm (used to compute GROW
#'   lambda when \code{lambda} is NULL). Default 0.45.
#' @param p_C_design Design alternative for control arm (used to compute GROW
#'   lambda when \code{lambda} is NULL). Default 0.30.
#' @param alpha Significance level for the rejection threshold (default 0.025).
#'
#' @return An object of class \code{"eprocess"} with components:
#' \describe{
#'   \item{log_evalue}{Numeric vector of cumulative log e-values at each observation.}
#'   \item{evalue}{Numeric vector of e-values at each observation.}
#'   \item{lambda}{The betting fraction(s) used.}
#'   \item{n}{Number of observations (per arm).}
#'   \item{alpha}{Significance level.}
#'   \item{threshold}{Log threshold \eqn{\log(1/\alpha)}.}
#'   \item{rejected}{Logical: did the e-process ever cross the threshold?}
#'   \item{rejection_time}{Integer: first observation at which threshold was crossed, or NA.}
#'   \item{D}{The pairwise differences \eqn{D_i}.}
#'   \item{endpoint}{Character: "binary".}
#' }
#'
#' @details
#' Under the null hypothesis \eqn{H_0: p_T = p_C}, the pairwise difference
#' \eqn{D_i = X_i^T - X_i^C} has mean zero regardless of the common response
#' rate. The product \eqn{E_n = \prod (1 + \lambda_i D_i)} is therefore a
#' nonnegative martingale (hence an e-process) under \eqn{H_0}, valid for
#' the composite null without requiring estimation of the nuisance parameter.
#'
#' The test rejects when \eqn{E_n \ge 1/\alpha}. By Ville's inequality,
#' \eqn{P_{H_0}(\sup_n E_n \ge 1/\alpha) \le \alpha}.
#'
#' @examples
#' set.seed(42)
#' x_T <- rbinom(200, 1, 0.45)
#' x_C <- rbinom(200, 1, 0.30)
#' ep <- eprocess_binary(x_T, x_C, lambda = 0.31)
#' summary(ep)
#' plot(ep)
#'
#' @export
eprocess_binary <- function(x_T, x_C, lambda = NULL,
                            p_T_design = 0.45, p_C_design = 0.30,
                            alpha = 0.025) {
  # --- Input validation ---
  stopifnot(is.numeric(x_T), is.numeric(x_C))
  stopifnot(length(x_T) == length(x_C))
  stopifnot(all(x_T %in% c(0L, 1L, 0, 1)))
  stopifnot(all(x_C %in% c(0L, 1L, 0, 1)))
  n <- length(x_T)
  stopifnot(n >= 1)

  # --- Determine lambda ---
  if (is.null(lambda)) {
    lambda <- grow_lambda(p_T_design, p_C_design)
  }
  if (length(lambda) == 1) {
    stopifnot(lambda > 0, lambda < 1)
    lambda_vec <- rep(lambda, n)
  } else {
    stopifnot(length(lambda) == n)
    stopifnot(all(lambda > 0 & lambda < 1))
    lambda_vec <- lambda
  }

  # --- Compute e-process ---
  D <- as.numeric(x_T) - as.numeric(x_C)
  log_increments <- log(pmax(1 + lambda_vec * D, .Machine$double.xmin))
  log_evalue <- cumsum(log_increments)
  evalue <- exp(log_evalue)

  # --- Threshold ---
  threshold <- log(1 / alpha)
  crossed <- log_evalue >= threshold
  if (any(crossed)) {
    rejection_time <- which(crossed)[1]
    rejected <- TRUE
  } else {
    rejection_time <- NA_integer_
    rejected <- FALSE
  }

  structure(
    list(
      log_evalue = log_evalue,
      evalue = evalue,
      lambda = if (length(unique(lambda_vec)) == 1) lambda_vec[1] else lambda_vec,
      n = n,
      alpha = alpha,
      threshold = threshold,
      rejected = rejected,
      rejection_time = rejection_time,
      D = D,
      endpoint = "binary"
    ),
    class = "eprocess"
  )
}
