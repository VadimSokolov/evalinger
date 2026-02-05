#' Create a Stateful E-Value Monitor
#'
#' Initialize a monitoring object that can be incrementally updated with
#' new patient data. Suitable for real-time DSMB monitoring.
#'
#' @param alpha Significance level (default 0.025).
#' @param lambda Betting fraction (default 0.31).
#' @param delta_min Minimum clinically important difference for futility
#'   assessment (default NULL, no futility monitoring).
#'
#' @return An object of class \code{"emonitor"}.
#'
#' @examples
#' mon <- emonitor(alpha = 0.025, lambda = 0.31)
#' # Update with first batch of 10 patients per arm
#' set.seed(1)
#' mon <- update(mon, x_T = rbinom(10, 1, 0.45), x_C = rbinom(10, 1, 0.30))
#' print(mon)
#'
#' @export
emonitor <- function(alpha = 0.025, lambda = 0.31, delta_min = NULL) {
  stopifnot(alpha > 0, alpha < 1, lambda > 0, lambda < 1)
  structure(
    list(
      alpha = alpha,
      lambda = lambda,
      delta_min = delta_min,
      n = 0L,
      log_evalue = 0,
      max_log_evalue = 0,
      av_pvalue = 1,
      rejected = FALSE,
      rejection_time = NA_integer_,
      history_log_e = numeric(0),
      history_D = numeric(0),
      cs = NULL,
      cumsum_T = 0L,
      cumsum_C = 0L
    ),
    class = "emonitor"
  )
}

#' Update an E-Value Monitor with New Data
#'
#' Add a new batch of patient outcomes to an existing e-value monitor and
#' update all statistics (e-value, always-valid p-value, confidence sequence).
#'
#' @param object An \code{"emonitor"} object.
#' @param x_T Integer vector of new treatment arm outcomes (0/1).
#' @param x_C Integer vector of new control arm outcomes (0/1), same length
#'   as \code{x_T}.
#' @param ... Additional arguments (ignored).
#'
#' @return Updated \code{"emonitor"} object.
#'
#' @examples
#' mon <- emonitor(alpha = 0.025, lambda = 0.31)
#' set.seed(1)
#' mon <- update(mon, x_T = rbinom(50, 1, 0.45), x_C = rbinom(50, 1, 0.30))
#' print(mon)
#'
#' @export
update.emonitor <- function(object, x_T, x_C, ...) {
  stopifnot(length(x_T) == length(x_C))
  stopifnot(all(x_T %in% c(0, 1)), all(x_C %in% c(0, 1)))

  D_new <- as.numeric(x_T) - as.numeric(x_C)
  log_inc <- log(pmax(1 + object$lambda * D_new, .Machine$double.xmin))

  new_log_e <- object$log_evalue + cumsum(log_inc)
  final_log_e <- new_log_e[length(new_log_e)]

  object$history_log_e <- c(object$history_log_e, new_log_e)
  object$history_D <- c(object$history_D, D_new)
  object$n <- object$n + length(x_T)
  object$log_evalue <- final_log_e
  object$max_log_evalue <- max(object$max_log_evalue, max(new_log_e))
  object$av_pvalue <- min(1, 1 / exp(object$max_log_evalue))
  object$cumsum_T <- object$cumsum_T + sum(x_T)
  object$cumsum_C <- object$cumsum_C + sum(x_C)

  threshold <- log(1 / object$alpha)
  if (!object$rejected && any(new_log_e >= threshold)) {
    object$rejected <- TRUE
    crossed_idx <- which(new_log_e >= threshold)[1]
    object$rejection_time <- object$n - length(x_T) + crossed_idx
  }

  # Simple confidence sequence via normal approximation
  if (object$n >= 10) {
    phat_T <- object$cumsum_T / object$n
    phat_C <- object$cumsum_C / object$n
    # Wald-type interval width inflated by log factor for CS validity
    se <- sqrt((phat_T * (1 - phat_T) + phat_C * (1 - phat_C)) / object$n)
    # Time-uniform width: use sqrt(2 * log(2/alpha) / n) as boundary
    cs_width <- sqrt(2 * log(2 / object$alpha) * (phat_T * (1 - phat_T) +
                     phat_C * (1 - phat_C)) / object$n)
    delta_hat <- phat_T - phat_C
    object$cs <- c(delta_hat - cs_width, delta_hat + cs_width)
  }

  object
}
