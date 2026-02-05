#' Hybrid E-Process + Group Sequential Monitor
#'
#' Run both an e-process and a group sequential boundary on the same trial data
#' simultaneously, reporting decisions from both approaches at each look.
#' This enables a direct comparison and provides dual evidence for DSMBs.
#'
#' @param x_T Integer vector of treatment arm binary outcomes (0/1).
#' @param x_C Integer vector of control arm binary outcomes (0/1).
#' @param look_times Integer vector of interim analysis times (observations per arm).
#' @param lambda Betting fraction for the e-process. If NULL, uses GROW-optimal.
#' @param alpha Significance level (default 0.025).
#' @param p_T_design Design alternative for treatment arm (default 0.45).
#' @param p_C_design Design alternative for control arm (default 0.30).
#' @param gs_type Type of group sequential boundary: \code{"obf"} (O'Brien-Fleming,
#'   default) or \code{"pocock"}.
#' @param gs_c Optional calibrated OBF constant \eqn{c} for an OBF-style boundary
#'   \eqn{z \ge c/\sqrt{t}}. If provided and \code{gs_type = "obf"}, the boundary
#'   uses this constant instead of \code{qnorm(1 - alpha)}. This is useful for
#'   matching simulation-calibrated group sequential boundaries.
#'
#' @return A data frame with one row per interim look and columns:
#' \describe{
#'   \item{look}{Look number (integer).}
#'   \item{n}{Cumulative observations per arm.}
#'   \item{info_frac}{Information fraction.}
#'   \item{delta_hat}{Estimated treatment effect.}
#'   \item{z_stat}{Z-statistic.}
#'   \item{gs_boundary}{Group sequential boundary.}
#'   \item{gs_reject}{Logical: GS rejection at this look.}
#'   \item{log_evalue}{Cumulative log e-value.}
#'   \item{evalue}{E-value.}
#'   \item{e_reject}{Logical: e-process rejection at this look.}
#'   \item{av_pvalue}{Always-valid p-value at this look.}
#' }
#'
#' @examples
#' set.seed(42)
#' x_T <- rbinom(200, 1, 0.45)
#' x_C <- rbinom(200, 1, 0.30)
#' hm <- hybrid_monitor(x_T, x_C, look_times = c(50, 100, 150, 200))
#' print(hm)
#'
#' @export
hybrid_monitor <- function(x_T, x_C, look_times, lambda = NULL,
                           alpha = 0.025,
                           p_T_design = 0.45, p_C_design = 0.30,
                           gs_type = "obf",
                           gs_c = NULL) {
  stopifnot(length(x_T) == length(x_C))
  n_total <- length(x_T)
  Nmax <- max(look_times)
  stopifnot(Nmax <= n_total)

  if (is.null(lambda)) {
    lambda <- grow_lambda(p_T_design, p_C_design)
  }

  # Compute full e-process
  D <- as.numeric(x_T) - as.numeric(x_C)
  log_e <- cumsum(log(pmax(1 + lambda * D, .Machine$double.xmin)))
  threshold <- log(1 / alpha)

  # GS boundaries
  info_frac <- look_times / Nmax
  z_alpha <- if (is.null(gs_c)) stats::qnorm(1 - alpha) else gs_c

  if (gs_type == "obf") {
    gs_bound <- z_alpha / sqrt(info_frac)
  } else if (gs_type == "pocock") {
    gs_bound <- rep(z_alpha * sqrt(log(1 + (exp(1) - 1) * info_frac[length(info_frac)])),
                    length(look_times))
    # Simplified Pocock-like constant boundary
    gs_bound <- rep(z_alpha + 0.5 * log(length(look_times)), length(look_times))
  } else {
    stop("gs_type must be 'obf' or 'pocock'")
  }

  # Build results
  rows <- vector("list", length(look_times))
  for (k in seq_along(look_times)) {
    t <- look_times[k]
    p_T_hat <- sum(x_T[1:t]) / t
    p_C_hat <- sum(x_C[1:t]) / t
    delta_hat <- p_T_hat - p_C_hat
    se <- sqrt(p_T_hat * (1 - p_T_hat) / t + p_C_hat * (1 - p_C_hat) / t)
    if (se < 1e-12) se <- 1e-6
    z <- delta_hat / se

    log_e_k <- log_e[t]
    e_k <- exp(log_e_k)
    max_log_e_k <- max(log_e[1:t])

    rows[[k]] <- data.frame(
      look = k,
      n = t,
      info_frac = info_frac[k],
      delta_hat = delta_hat,
      z_stat = z,
      gs_boundary = gs_bound[k],
      gs_reject = z >= gs_bound[k],
      log_evalue = log_e_k,
      evalue = e_k,
      e_reject = max_log_e_k >= threshold,
      av_pvalue = min(1, 1 / exp(max_log_e_k)),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

#' Convert gsDesign Object to E-Process Design
#'
#' Extract the interim analysis schedule from a \pkg{gsDesign} object and
#' compute the corresponding e-process design parameters.
#'
#' @param gs_obj A \code{gsDesign} object from the \pkg{gsDesign} package.
#' @param p_C Control arm response probability.
#' @param p_T Treatment arm response probability.
#'
#' @return A list with components:
#' \describe{
#'   \item{look_times}{Integer vector of look times.}
#'   \item{lambda}{GROW-optimal betting fraction.}
#'   \item{alpha}{Alpha from the gsDesign object.}
#'   \item{gs_upper}{Upper boundaries from the gsDesign.}
#'   \item{gs_lower}{Lower boundaries from the gsDesign (if present).}
#' }
#'
#' @examples
#' \dontrun{
#' library(gsDesign)
#' gs <- gsDesign(k = 4, test.type = 1, alpha = 0.025, beta = 0.20,
#'                sfu = sfLDOF)
#' as_eprocess_design(gs, p_C = 0.30, p_T = 0.45)
#' }
#'
#' @export
as_eprocess_design <- function(gs_obj, p_C, p_T) {
  if (!requireNamespace("gsDesign", quietly = TRUE)) {
    stop("Package 'gsDesign' is required but not installed.")
  }
  stopifnot(inherits(gs_obj, "gsDesign"))

  # Extract info fractions and compute look times
  info_frac <- gs_obj$timing
  n_max <- ceiling(gs_obj$n.I[length(gs_obj$n.I)])
  look_times <- round(info_frac * n_max)

  lambda <- grow_lambda(p_T, p_C)

  list(
    look_times = as.integer(look_times),
    lambda = lambda,
    alpha = gs_obj$alpha,
    gs_upper = gs_obj$upper$bound,
    gs_lower = if (!is.null(gs_obj$lower)) gs_obj$lower$bound else NULL
  )
}
