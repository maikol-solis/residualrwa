
# Plot for residualrwa objects
#'
#' @param x an object of class "residualrwa"
#' @param title an optional title for the plot
#'
#' @return Returns a ggplot object
#' @export
#'
#' @examples
plot_rwa <- function(x, title = "") {
  df <- x$data_frame
  df$variable <- reorder(df$variable, df$weight)

  if (x$boot_ci == TRUE) {
    df_boot <- x$data_frame_boot
    df_boot$variable <- reorder(df_boot$variable, df_boot$weight)
    ggplot2::ggplot(data = df_boot) +
      ggplot2::geom_boxplot(ggplot2::aes(x = variable, y = weight)) +
      ggplot2::geom_point(
        data = df,
        ggplot2::aes(x = variable, y = weight),
        color = "red", fill = "red", size = 2, shape = 23
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::scale_y_continuous(labels = scales::percent) +
      ggplot2::ggtitle(title) +
      ggplot2::coord_flip()
  } else {
    ggplot2::ggplot(data = df) +
      ggplot2::geom_segment(ggplot2::aes(
        x = variable,
        xend = variable,
        y = 0,
        yend = weight
      )) +
      ggplot2::geom_point(
        ggplot2::aes(x = variable, y = weight),
        color = "red", fill = "red", size = 2, shape = 23
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::scale_y_continuous(labels = scales::percent) +
      ggplot2::ggtitle(title) +
      ggplot2::coord_flip()
  }
}
