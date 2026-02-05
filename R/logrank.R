#' Safe Logrank E-Process for Survival Endpoints
#'
#' Construct an e-process for the two-sample survival comparison using
#' the safe logrank test. When the \pkg{safestats} package is available, this
#' function wraps its logrank implementation; otherwise, a self-contained
#' fallback using the Breslow–Nelson–Aalen hazard increments is used.
#'
#' @param time Numeric vector of observed event/censoring times.
#' @param status Integer vector: 1 = event, 0 = censored.
#' @param arm Factor or integer vector: treatment arm indicator (0 = control,
#'   1 = treatment).
#' @param theta Design alternative for the hazard ratio (default 0.7,
#'   i.e., 30\% reduction in hazard).
#' @param alpha Significance level (default 0.025).
#'
#' @return An object of class \code{"eprocess"} with components:
#' \describe{
#'   \item{log_evalue}{Cumulative log e-values at each event time.}
#'   \item{evalue}{E-values at each event time.}
#'   \item{n}{Number of events.}
#'   \item{alpha}{Significance level.}
#'   \item{threshold}{Log rejection threshold.}
#'   \item{rejected}{Logical: was the threshold crossed?}
#'   \item{rejection_time}{Event index at which threshold was crossed, or NA.}
#'   \item{endpoint}{\code{"survival"}.}
#'   \item{event_times}{Unique event times.}
#' }
#'
#' @details
#' The safe logrank e-process is a betting martingale based on the partial
#' likelihood increments of the Cox model. At each event time, the e-process
#' bets on whether the event came from the treatment or control arm, with
#' stakes calibrated by the design alternative \code{theta}.
#'
#' The fallback implementation computes the standard logrank score increment
#' at each event time and constructs a martingale using a GROW-like betting
#' fraction derived from \code{theta}.
#'
#' @examples
#' \dontrun{
#' # Simulated trial with exponential survival
#' set.seed(42)
#' n <- 100
#' arm <- rep(0:1, each = n/2)
#' time <- rexp(n, rate = ifelse(arm == 1, 0.7, 1.0))
#' status <- rbinom(n, 1, 0.8)  # 80% event rate
#' ep <- eprocess_logrank(time, status, arm, theta = 0.7)
#' summary(ep)
#' }
#'
#' @export
eprocess_logrank <- function(time, status, arm, theta = 0.7, alpha = 0.025) {
  stopifnot(length(time) == length(status), length(time) == length(arm))
  stopifnot(all(status %in% c(0, 1)))
  stopifnot(theta > 0, theta < 1)
  stopifnot(alpha > 0, alpha < 1)

  arm <- as.integer(arm)
  stopifnot(all(arm %in% c(0L, 1L)))

  # Try to use safestats if available
  if (requireNamespace("safestats", quietly = TRUE)) {
    return(.eprocess_logrank_safestats(time, status, arm, theta, alpha))
  }

  # Fallback implementation
  .eprocess_logrank_fallback(time, status, arm, theta, alpha)
}

#' @keywords internal
.eprocess_logrank_safestats <- function(time, status, arm, theta, alpha) {
  # Use safestats package
  df <- data.frame(time = time, status = status, arm = factor(arm))
  tryCatch({
    design <- safestats::designSafeLogrank(hrMin = theta, alpha = alpha)
    result <- safestats::safeLogrankTest(
      survival::Surv(time, status) ~ arm,
      data = df,
      designObj = design
    )
    # Extract e-values -- safestats returns a single e-value
    # For the sequential version, we build incremental e-values
    ev <- result$eValue
    log_ev <- log(max(ev, .Machine$double.xmin))

    structure(
      list(
        log_evalue = log_ev,
        evalue = ev,
        n = sum(status),
        alpha = alpha,
        threshold = log(1 / alpha),
        rejected = log_ev >= log(1 / alpha),
        rejection_time = if (log_ev >= log(1 / alpha)) sum(status) else NA_integer_,
        endpoint = "survival",
        event_times = sort(unique(time[status == 1])),
        lambda = NA
      ),
      class = "eprocess"
    )
  },
  error = function(e) {
    message("safestats failed, using fallback: ", e$message)
    .eprocess_logrank_fallback(time, status, arm, theta, alpha)
  })
}

#' @keywords internal
.eprocess_logrank_fallback <- function(time, status, arm, theta, alpha) {
  # Sort by time
  ord <- order(time)
  time <- time[ord]
  status <- status[ord]
  arm <- arm[ord]

  n <- length(time)
  event_idx <- which(status == 1)
  n_events <- length(event_idx)

  if (n_events == 0) {
    return(structure(
      list(
        log_evalue = numeric(0), evalue = numeric(0),
        n = 0L, alpha = alpha, threshold = log(1 / alpha),
        rejected = FALSE, rejection_time = NA_integer_,
        endpoint = "survival", event_times = numeric(0),
        lambda = NA
      ),
      class = "eprocess"
    ))
  }

  # Betting fraction based on design alternative
  # Under Cox model: lambda ~ (theta - 1) / (theta + 1) * correction
  lambda_lr <- (1 - theta) / (1 + theta)

  log_evalue <- numeric(n_events)
  cumlog <- 0

  for (j in seq_along(event_idx)) {
    i <- event_idx[j]
    # At-risk set
    at_risk <- which(time >= time[i])
    n_at_risk <- length(at_risk)
    n_trt_at_risk <- sum(arm[at_risk] == 1)

    if (n_at_risk < 2 || n_trt_at_risk == 0 || n_trt_at_risk == n_at_risk) {
      log_evalue[j] <- cumlog
      next
    }

    # Expected proportion of treatment in risk set
    p_trt <- n_trt_at_risk / n_at_risk
    # Observed: was the event in treatment arm?
    obs_trt <- as.numeric(arm[i] == 1)
    # Logrank score increment: obs_trt - p_trt
    score <- obs_trt - p_trt

    # Betting increment
    inc <- log(pmax(1 + lambda_lr * score, .Machine$double.xmin))
    cumlog <- cumlog + inc
    log_evalue[j] <- cumlog
  }

  evalue <- exp(log_evalue)
  threshold <- log(1 / alpha)
  crossed <- log_evalue >= threshold

  structure(
    list(
      log_evalue = log_evalue,
      evalue = evalue,
      n = n_events,
      alpha = alpha,
      threshold = threshold,
      rejected = any(crossed),
      rejection_time = if (any(crossed)) which(crossed)[1] else NA_integer_,
      endpoint = "survival",
      event_times = time[event_idx],
      lambda = lambda_lr
    ),
    class = "eprocess"
  )
}
