


#' Relative Weight Analysis for logistic regression
#'
#' @param X matrix or data.frame with inputs
#' @param Y numeric vector or data.frame with output.
#' @param data data.frame
#'
#' @return a list with the R^2, adjusted R^2, and relative weights for the X with respect to Y.
#'
#' @keywords internal
RWA <- function(X, Y, data, family) {
  X.svd <- svd(X)
  Q <- X.svd$v
  P <- X.svd$u
  Z <- P %*% t(Q)

  Z.stand <- scale(Z)

  # Obtaining Lambda from equation 7 from Johnson (2000) pg 8
  Lambda <- qr.solve(t(Z.stand) %*% Z.stand) %*% t(Z.stand) %*% as.matrix(X)

  fit <- stats::glm(Y ~ Z.stand, family = family)
  unstCoefs <- stats::coef(fit)
  predY <- stats::predict(fit, newdata = data, type = "response")
  Yhat <- fit$fitted.values # Creating Y-hat

  getting.Rsq <- stats::lm(predY ~ Y) # Getting R^2
  Rsq <- summary(getting.Rsq)$r.squared
  AdjRsq <- summary(getting.Rsq)$adj.r.squared


  if (family$family == "binomial") {
    b <- unstCoefs[2:length(unstCoefs)]
    # Getting stdev of logit-Y-hat
    stdYhat <- stats::sd(Yhat)
    # Computing standardized logistic regression coefficients
    beta <- b * ((sqrt(Rsq)) / stdYhat)
  } else if (family$family == "gaussian") {
    beta <- unstCoefs[2:length(unstCoefs)]
  }



  epsilon <- Lambda^2 %*% beta^2
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
  unlist(sapply(colnames(data), function(column.name) {
    na.omit(
      stringr::str_match(type_variable, column.name)
    )
  }))
}


extract_X_RWA <- function(model, interactions, family) {
  # Define Design matrix
  XDesign <- model$x
  idx <- model$assign

  # Case when there are interactions
  if (interactions) {
    col.interactions <-
      which(model$Design$assume == "interaction")
    # min.interaction <- min(col.interactions)
    # design.interaction <- model$Design$interactions


    for (k in col.interactions) {
      xx <- stringr::str_split(colnames(model$x)[idx[[k]] - 1], " \\* |:")

      if (length(xx) == 1) {
        ll <- lm(model$x[, idx[[k]] - 1] ~ -1 + model$x[, xx[[1]]])
        XDesign[, idx[[k]] - 1] <- resid(ll)
      } else {
        for (i in seq_along(xx)) {
          ll <- lm(model$x[, idx[[k]] - 1][, i] ~ -1 + model$x[, xx[[i]]])
          XDesign[, idx[[k]] - 1][, i] <- resid(ll)
        }
      }
    }
  }

  # Define the new fitting from the Design matrix
  if (family$family == "binomial") {
    ff <- rms::lrm.fit(x = XDesign, y = model$y, tol = 1e-9)
  } else if (family$family == "gaussian") {
    ff <- rms::ols(model$y ~ XDesign, tol = 1e-9)
  }

  # ff <- stats::glm.fit(x = cbind(1,XDesign), y = model$y, family = family)
  cc <- ff$coefficients[-1]
  X <- NULL


  # Consolidate a new design matrix grouping the splines terms
  for (j in seq_along(model$assign)) {
    nn <- names(idx[j])
    m <- as.matrix(XDesign[, idx[[j]] - 1])
    sc <- cc[idx[[j]] - 1]
    Xtmp <- m %*% sc
    colnames(Xtmp) <- nn
    X <- cbind(X, Xtmp)
  }

  return(X)
}


# R2M <- function(x, formula, data) {
#   LM <- stats::update(x, formula, data = data)
#   response.name <- all.vars(formula)[1]
#   L0 <-  stats::update(x, formula(paste0(response.name, "~", 1)),
#                        data = data)
#   return(1 - stats::logLik(LM) / stats::logLik(L0))
# }
