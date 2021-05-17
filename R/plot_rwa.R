
# Plot for modelrwa objects
#'
#' @param x an object of class "modelrwa"
#' @param title an optional title for the plot
#'
#' @return Returns a ggplot object
#' @export
#'
#' @examples
#'
plot_rwa <- function(x, title = "") {

  df <- x$data_frame
  df$Variable <- reorder(df$Variable, df$Weight)

  p2 <- ggplot2::ggplot(data = df) +
    ggplot2::geom_bar(ggplot2::aes(Variable, Weight), stat = "identity") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::ggtitle(title) +
    ggplot2::coord_flip()


  return(p2)

}
