#' Plot for residualrwa objects
#'
#' @param x An object of class \code{"residualrwa"}.
#' @param boot_ci A boolean. Its value determines if the function draws or not
#'   the boostrap confidence. If misssing, then the function uses the value in
#'   \code{x$boot_ci}.
#' @param font_size A numeric value. Overall font size. It uses the same
#'   mechanism as \code{\link[cowplot]{theme_cowplot}}.
#' @param ... Currently ignored.
#'
#' @return Returns a ggplot object
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
#' plot(ex)
#'
#' @export

plot.residualrwa <- function(x, boot_ci, font_size = 12, ...) {
  ## Pacify R checks
  variable <- weight <- ci_low <- ci_up <- NULL

  df <- x$data_frame
  df$variable <- stats::reorder(df$variable, df$weight)

  if (missing(boot_ci)) {
    boot_ci <- x$boot
  } else {
    if (inherits(boot_ci, "logical")) {
      stop("boot_ci must be logical.")
    }
  }

  p <- ggplot2::ggplot(data = df) +
    ggplot2::geom_point(
      ggplot2::aes(x = variable, y = weight),
      color = "red", fill = "red", size = 3, shape = 23
    )

  if (boot_ci == TRUE) {
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(
      x = variable,
      ymin = ci_low,
      ymax = ci_up
    ))
  } else {
    p <- p + ggplot2::geom_segment(ggplot2::aes(
      x = variable,
      xend = variable,
      y = 0,
      yend = weight
    ))
  }

  p <- p + ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::coord_flip() +
    cowplot::theme_minimal_grid(font_size = font_size)


  return(p)
}
