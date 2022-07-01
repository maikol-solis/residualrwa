## Estimate relative weights
estimate_rwa <- function(x, y, data, family) {
  x_svd <- svd(x)
  q_matrix <- x_svd$v
  p_matrix <- x_svd$u
  z_matrix <- p_matrix %*% t(q_matrix)

  z_standard <- scale(z_matrix)

  # Obtaining Lambda from equation 7 from Johnson (2000) pg 8
  lambda <- qr.solve(t(z_standard) %*% z_standard) %*%
    t(z_standard) %*% as.matrix(x)

  fit <- stats::glm(y ~ z_standard, family = family)
  unstandardized_coefs <- stats::coef(fit)
  prediction_y <- stats::predict(fit, newdata = data, type = "response")
  y_hat <- fit$fitted.values # Creating Y-hat

  getting_rsq <- stats::lm(prediction_y ~ y) # Getting R^2
  r_sq <- summary(getting_rsq)$r.squared
  adj_r_sq <- summary(getting_rsq)$adj.r.squared


  if (family$family == "binomial") {
    b <- unstandardized_coefs[2:length(unstandardized_coefs)]
    # Getting stdev of logit-Y-hat
    std_y_hat <- stats::sd(y_hat)
    # Computing standardized logistic regression coefficients
    beta <- b * ((sqrt(r_sq)) / std_y_hat)
  } else if (family$family == "gaussian") {
    beta <- unstandardized_coefs[2:length(unstandardized_coefs)]
  }

  epsilon <- lambda^2 %*% beta^2
  sum_epsilons <- sum(epsilon)
  prop_weights <- (epsilon / sum_epsilons)
  names(prop_weights) <- colnames(x)

  return(list(
    adj_r_sq = adj_r_sq,
    r_sq = r_sq,
    sum_epsilons = sum_epsilons,
    epsilon = epsilon,
    prop_weights = prop_weights
  ))
}

## Extract columns names from data
extract_column_names <- function(data, type_variable) {
  unlist(sapply(
    X = colnames(data),
    FUN = function(column_name) {
      stats::na.omit(stringr::str_match(type_variable, column_name))
    }
  ))
}

## Consolidate spline columns into one factor
consolidate_design_matrix <- function(model, interactions, family) {
  # Define Design matrix
  x_design <- model$x
  idx <- model$assign

  # Case when there are interactions
  if (interactions) {
    cols_interactions <- which(model$Design$assume == "interaction")


    for (k in cols_interactions) {
      xx <- stringr::str_split(colnames(model$x)[idx[[k]] - 1], " \\* |:")

      if (length(xx) == 1) {
        ll <- stats::lm(model$x[, idx[[k]] - 1] ~ -1 + model$x[, xx[[1]]])
        x_design[, idx[[k]] - 1] <- stats::resid(ll)
      } else {
        for (i in seq_along(xx)) {
          ll <- stats::lm(
            model$x[, idx[[k]] - 1][, i] ~ -1 + model$x[, xx[[i]]]
          )
          x_design[, idx[[k]] - 1][, i] <- stats::resid(ll)
        }
      }
    }
  }

  # Define the new fitting from the Design matrix
  if (family$family == "binomial") {
    ff <- rms::lrm.fit(x = x_design, y = model$y, tol = 1e-9)
  } else if (family$family == "gaussian") {
    ff <- rms::ols(model$y ~ x_design, tol = 1e-9)
  }
  cc <- ff$coefficients[-1]
  x <- NULL


  # Consolidate a new design matrix grouping the splines terms
  for (j in seq_along(model$assign)) {
    nn <- names(idx[j])
    m <- as.matrix(x_design[, idx[[j]] - 1])
    sc <- cc[idx[[j]] - 1]
    x_tmp <- m %*% sc
    colnames(x_tmp) <- nn
    x <- cbind(x, x_tmp)
  }

  return(x)
}


# R2M <- function(x, formula, data) {
#   LM <- stats::update(x, formula, data = data)
#   response.name <- all.vars(formula)[1]
#   L0 <-  stats::update(x, formula(paste0(response.name, "~", 1)),
#                        data = data)
#   return(1 - stats::logLik(LM) / stats::logLik(L0))
# }
