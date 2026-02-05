#' E-BH Procedure for Multiple Testing with E-Values
#'
#' Applies the e-Benjamini-Hochberg procedure to control the false discovery
#' rate (FDR) when testing multiple hypotheses using e-values. Unlike the
#' classical BH procedure with p-values, the e-BH procedure remains valid
#' under arbitrary dependence between the e-values.
#'
#' @param e_values Numeric vector of e-values (nonnegative).
#' @param alpha Target FDR level (default 0.05).
#'
#' @return A list with components:
#' \describe{
#'   \item{rejected}{Logical vector: which hypotheses are rejected.}
#'   \item{n_rejected}{Number of rejections.}
#'   \item{threshold}{The e-value threshold used.}
#'   \item{adjusted_e}{Adjusted e-values for reference.}
#'   \item{alpha}{Target FDR level.}
#' }
#'
#' @details
#' The e-BH procedure (Wang and Ramdas, 2022) sorts the e-values in decreasing
#' order and finds the largest \eqn{k} such that \eqn{e_{(k)} \ge K / (k \alpha)},
#' where \eqn{K} is the total number of hypotheses. It then rejects the
#' hypotheses corresponding to the \eqn{k} largest e-values.
#'
#' The key advantage over p-value-based BH is that e-BH controls FDR under
#' arbitrary dependence, without requiring independence or PRDS assumptions.
#'
#' @references
#' Wang, R. and Ramdas, A. (2022). False discovery rate control with e-values.
#' \emph{Journal of the Royal Statistical Society: Series B}, 84(3), 822-852.
#'
#' @examples
#' # 10 hypotheses, 3 with signal
#' set.seed(42)
#' e_vals <- c(rexp(3, rate = 0.01), rexp(7, rate = 1))  # 3 large, 7 small
#' result <- ebh(e_vals, alpha = 0.05)
#' result$n_rejected
#'
#' @export
ebh <- function(e_values, alpha = 0.05) {
  stopifnot(is.numeric(e_values), all(e_values >= 0))
  stopifnot(alpha > 0, alpha < 1)

  K <- length(e_values)
  if (K == 0) {
    return(list(rejected = logical(0), n_rejected = 0L,
                threshold = Inf, adjusted_e = numeric(0), alpha = alpha))
  }

  # Sort e-values in decreasing order
  ord <- order(e_values, decreasing = TRUE)
  e_sorted <- e_values[ord]

  # Find largest k such that e_{(k)} >= K / (k * alpha)
  thresholds <- K / (seq_len(K) * alpha)
  eligible <- which(e_sorted >= thresholds)

  if (length(eligible) == 0) {
    return(list(
      rejected = rep(FALSE, K),
      n_rejected = 0L,
      threshold = K / alpha,
      adjusted_e = e_values,
      alpha = alpha
    ))
  }

  k_star <- max(eligible)
  # Reject the top k_star hypotheses
  rejected <- rep(FALSE, K)
  rejected[ord[seq_len(k_star)]] <- TRUE

  list(
    rejected = rejected,
    n_rejected = k_star,
    threshold = K / (k_star * alpha),
    adjusted_e = e_values,
    alpha = alpha
  )
}

#' Platform Trial Multi-Arm Monitoring
#'
#' Set up and run e-value monitoring for a platform trial with multiple
#' treatment arms sharing a common control. Each arm is monitored
#' independently with its own e-process, while multiplicity is handled
#' via the e-BH procedure applied at each interim look.
#'
#' @param K Number of treatment arms.
#' @param look_times Integer vector of interim analysis times.
#' @param x_C Integer vector of shared control arm outcomes (0/1) of length
#'   at least \code{max(look_times)}.
#' @param x_T_list A list of \code{K} integer vectors, each of length at least
#'   \code{max(look_times)}, giving the outcomes for each treatment arm.
#' @param lambda_vec Numeric vector of length \code{K} with betting fractions
#'   per arm. If NULL, all arms use the same default (0.31).
#' @param alpha Overall Type I error rate (default 0.025).
#' @param fdr_alpha FDR level for the e-BH procedure at each look (default 0.05).
#'
#' @return A list with components:
#' \describe{
#'   \item{arm_results}{A list of \code{K} lists, each with \code{log_evalue},
#'     \code{rejected}, \code{rejection_time}.}
#'   \item{ebh_results}{A list of e-BH results at each look time.}
#'   \item{summary}{A data frame summarizing which arms were rejected and when.}
#' }
#'
#' @examples
#' set.seed(42)
#' K <- 3
#' n <- 200
#' x_C <- rbinom(n, 1, 0.30)
#' x_T_list <- list(
#'   rbinom(n, 1, 0.45),  # effective
#'   rbinom(n, 1, 0.30),  # null
#'   rbinom(n, 1, 0.40)   # moderate effect
#' )
#' pm <- platform_monitor(K = K, look_times = c(50, 100, 150, 200),
#'                        x_C = x_C, x_T_list = x_T_list)
#' pm$summary
#'
#' @export
platform_monitor <- function(K, look_times, x_C, x_T_list,
                             lambda_vec = NULL, alpha = 0.025,
                             fdr_alpha = 0.05) {
  stopifnot(length(x_T_list) == K)
  Nmax <- max(look_times)
  stopifnot(length(x_C) >= Nmax)
  for (k in seq_len(K)) {
    stopifnot(length(x_T_list[[k]]) >= Nmax)
  }

  if (is.null(lambda_vec)) {
    lambda_vec <- rep(0.31, K)
  }
  stopifnot(length(lambda_vec) == K)

  threshold <- log(1 / alpha)

  # Compute e-processes for each arm
  arm_results <- vector("list", K)
  for (k in seq_len(K)) {
    D <- as.numeric(x_T_list[[k]]) - as.numeric(x_C)
    log_e <- cumsum(log(pmax(1 + lambda_vec[k] * D[1:Nmax], .Machine$double.xmin)))

    crossed <- log_e >= threshold
    arm_results[[k]] <- list(
      log_evalue = log_e,
      evalue = exp(log_e),
      rejected = any(crossed),
      rejection_time = if (any(crossed)) which(crossed)[1] else NA_integer_
    )
  }

  # Apply e-BH at each look
  ebh_results <- vector("list", length(look_times))
  for (j in seq_along(look_times)) {
    t <- look_times[j]
    e_at_look <- sapply(arm_results, function(a) a$evalue[t])
    ebh_results[[j]] <- ebh(e_at_look, alpha = fdr_alpha)
  }

  # Summary
  summary_df <- data.frame(
    arm = seq_len(K),
    individual_rejected = sapply(arm_results, function(a) a$rejected),
    individual_reject_time = sapply(arm_results, function(a) {
      if (is.na(a$rejection_time)) NA_integer_ else a$rejection_time
    }),
    ebh_rejected_final = ebh_results[[length(look_times)]]$rejected,
    stringsAsFactors = FALSE
  )

  list(
    arm_results = arm_results,
    ebh_results = ebh_results,
    summary = summary_df
  )
}
