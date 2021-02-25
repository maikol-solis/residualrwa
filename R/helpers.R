
#' Relative Weight Analysis for logistic regression
#'
#' @param X matrix or data.frame with inputs
#' @param Y numeric vector or data.frame with output.
#' @param data data.frame
#'
#' @return a list with the R^2, adjusted R^2, and relative weights for the X with respect to Y.
#'
#' @keywords internal
RWA <- function(X, Y, data) {
  X.svd <- svd(X)
  Q <- X.svd$v
  P <- X.svd$u
  Z <- P %*% t(Q)

  Z.stand <- scale(Z)

  Lambda <-
    qr.solve(t(Z.stand) %*% Z.stand) %*% t(Z.stand) %*% as.matrix(X) # Obtaining Lambda from equation 7 from Johnson (2000) pg 8

  logrfit <- stats::glm(Y ~ Z.stand, family = stats::binomial)
  unstCoefs <- stats::coef(logrfit)
  b <- unstCoefs[2:length(unstCoefs)]

  LpredY <-
    stats::predict(logrfit, newdata = data, type = "response")
  lYhat <- log(LpredY / (1 - LpredY)) # Creating logit-Y-hat
  stdlYhat <- stats::sd(lYhat) # Getting stdev of logit-Y-hat
  getting.Rsq <- stats::lm(LpredY ~ Y) # Getting R-sq
  Rsq <- summary(getting.Rsq)$r.squared
  AdjRsq <- summary(getting.Rsq)$adj.r.squared
  beta <-
    b * ((sqrt(Rsq)) / stdlYhat) # Computing standardized logistic regression coefficients
  epsilon <- Lambda ^ 2 %*% beta ^ 2
  sum.epsilons <- sum(epsilon)
  PropWeights <- (epsilon / sum.epsilons)
  names(PropWeights) <- colnames(X)
  return(
    list(
      AdjRsq = AdjRsq,
      Rsq = Rsq,
      sum.epsilons = sum.epsilons,
      epsilon = epsilon,
      PropWeights = PropWeights
    )
  )
}


extract_column_names <- function(data, type_variable) {
  unlist(sapply(colnames(data), function(column.name)
    na.omit(
      stringr::str_match(type_variable, column.name)
    )))
}


# R2M <- function(x, formula, data) {
#   LM <- stats::update(x, formula, data = data)
#   response.name <- all.vars(formula)[1]
#   L0 <-  stats::update(x, formula(paste0(response.name, "~", 1)),
#                        data = data)
#   return(1 - stats::logLik(LM) / stats::logLik(L0))
# }


