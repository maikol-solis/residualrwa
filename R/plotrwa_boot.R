# Plot for bootmodelrwa objects
#'
#' @param x an object of class "modelrwa"
#' @param title an optional title for the plot
#'
#' @return Returns a ggplot object
#' @export
#'
#' @examples
#'




plotbootrwa <- function(x,title=""){

  x$Variable <- reorder(x$Variable, x$Weight)
  p2 <-
    ggplot2::ggplot(x) +
    geom_boxplot(ggplot2::aes(Variable, Weight)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::ggtitle(title) +
    ggplot2::coord_flip()

  return(p2)

}
