#' evalinger: E-Values for Adaptive Clinical Trial Monitoring
#'
#' Implements e-value and e-process methodology for interim monitoring
#' of adaptive clinical trials. Provides betting-martingale e-process
#' construction for two-arm binary trials, safe logrank tests for survival
#' endpoints, always-valid confidence sequences, futility monitoring,
#' design calibration (GROW optimal betting fractions, expected stopping
#' times), and simulation tools for comparing e-value monitoring with group
#' sequential and Bayesian approaches.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{eprocess_binary}}}{Construct e-process for binary endpoints}
#'   \item{\code{\link{eprocess_logrank}}}{Safe logrank e-process for survival endpoints}
#'   \item{\code{\link{grow_lambda}}}{GROW-optimal betting fraction}
#'   \item{\code{\link{edesign_binary}}}{Design calibration}
#'   \item{\code{\link{emonitor}}}{Stateful incremental monitor}
#'   \item{\code{\link{confseq_binary}}}{Always-valid confidence sequences}
#'   \item{\code{\link{simulate_comparison}}}{Compare monitoring methods}
#'   \item{\code{\link{hybrid_monitor}}}{Side-by-side e-process + GS monitoring}
#'   \item{\code{\link{platform_monitor}}}{Multi-arm platform trial monitoring}
#'   \item{\code{\link{ebh}}}{E-BH procedure for FDR control}
#' }
#'
#' @docType package
#' @name evalinger-package
#' @keywords internal
#'
#' @importFrom graphics abline barplot legend lines par plot polygon
#' @importFrom stats integrate
"_PACKAGE"

## Suppress R CMD check NOTEs for ggplot2 .data pronoun usage.
## These variables appear in aes() calls inside code guarded by
## requireNamespace("ggplot2", quietly = TRUE).
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".data"
  ))
}
