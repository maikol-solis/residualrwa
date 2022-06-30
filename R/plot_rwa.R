
# Plot for residualrwa objects
#'
#' @param x an object of class "residualrwa"
#' @param title an optional title for the plot
#'
#' @return Returns a ggplot object
#' @export
#'
#' @examples
plot_rwa <- function(x, boot_ci, font_size = 12) {
  ## Pacify R checks
  variable <- weight <- ci_low <- ci_up <- NULL

  df <- x$data_frame
  df$variable <- reorder(df$variable, df$weight)

  if (missing(boot_ci)) {
    boot_ci <- x$boot
  } else {
    if (class(boot_ci) != "logical") {
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
