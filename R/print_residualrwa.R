#' Print residualrwa objects
#'
#' @param x An object of class \code{"residualrwa"}.
#' @param ... Currently ignored.
#'
#' @return Print summarized information about the \code{residualrwa} object.
#'
#' @examples
#'
#' n <- 100
#' X1 <- rnorm(n)
#' X2 <- rnorm(n)
#' X3 <- rnorm(n)
#' Y <- X2^3 + 10 * X1 * X2
#' data <- as.data.frame(cbind(Y, X1, X2, X3))
#'
#' ex <- residualrwa(
#'   response = "Y",
#'   control = NULL,
#'   fixed = NULL,
#'   free = c("X1", "X2", "X3"),
#'   data = data,
#'   include_interactions = TRUE,
#'   boot_ci = TRUE,
#'   n_boot = 5
#' )
#'
#' print(ex)
#'
#' @export

print.residualrwa <- function(x, ...) {
  cat("\nCall:\n", deparse(x[["call"]]), "\n", sep = "")
  cat("\nNumber of variables:", ncol(x[["data"]]) - 1, "\n")
  cat("\nNumber of observations:", nrow(x[["data"]]), "\n")
  cat("\nFree variables:", x[["variables"]][["free"]], "\n")
  cat("\nFixed variables:", x[["variables"]][["fixed"]], "\n")
  cat("\nControl variables:", x[["variables"]][["control"]], "\n")
  cat("\nInteraction variables:", x[["variables"]][["interaction"]], "\n")
  print(as.data.frame(x$summary))
}
