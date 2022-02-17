#' #' Detect influential variables in model with Relative Weight Analysis
#'
#' @param response.name name for independent variable
#' @param control character vector with names for control variables
#' @param fixed character vector with names for fixed variables
#' @param free character vector with names for free variables
#' @param data object data.frame with data
#' @param family type of regression (see \code{\link[stats]{family}})
#' @param include.interactions A boolean indicating if the model should calculate all the pairwise interactions between variables. It uses the names in the parameters \code{fixed} and \code{variables}.
#' @param alpha level of significance use if \code{method = "p"}.
#' @param method method used to keep interactions if \code{include.interactions=TRUE}. See \code{\link[rms]{fastbw}}.
#' @param name.control,name.fixed,name.variables,name.interactions Names used to label the summary tables
#' @param verbose if \code{TRUE}, \code{residualrwa} shows all the stepwise process. Defaults to \code{FALSE}.
#' @return A modelrwa object
#' @export
#'
#' @examples
#'
#' X1 <- rnorm(1000)
#' X2 <- rnorm(1000)
#' X3 <- rnorm(1000)
#'
#' Y <- plogis(X2^3+10*X1*X2)>0.5
#'
#' data <- as.data.frame(cbind(Y,X1,X2,X3))
#'
#' ex1 <- modelrwa(
#'   response.name = "Y" ,
#'   control = NULL,
#'   fixed = NULL,
#'   variables = c("X1","X2", "X3"),
#'   data = data,
#'   family = binomial,
#'   include.interactions = TRUE
#' )
#'
#'



residualrwa_boot <- function(response.name,
                             control = NULL,
                             fixed = NULL,
                             free,
                             data,
                             family = stats::gaussian(),
                             include.interactions = FALSE,
                             alpha = 0.1,
                             method = c("aic", "p"),
                             name.control = "Control",
                             name.fixed = "Fixed",
                             name.free = "Free",
                             name.interactions = "Interactions",
                             verbose = FALSE,
                             nboot = 100) {
  rwaBoots <- lapply(X = 1:nboot,
                     function(i) {
                       message(paste0("Boot sample #", i))
                       data_boot <- data[sample(nrow(data), nrow(data), replace = T), ]

                       exRWA <- residualrwa(
                         response.name,
                         control,
                         fixed,
                         free,
                         data = data_boot,
                         family,
                         include.interactions,
                         alpha,
                         method,
                         name.control,
                         name.fixed,
                         name.free,
                         name.interactions,
                         verbose
                       )

                       return(data.frame(exRWA$data_frame, nboot = rep(i, nrow(exRWA$data_frame))))

                     })

  result <- dplyr::bind_rows(rwaBoots)

  result <- tidyr::complete(data = result, Variable, nboot)

  result$Weight <- ifelse(is.na(result$Weight), 0, result$Weight)

  #Completo en caso de que no aparezca en la selecci'on

  return(result)

}
