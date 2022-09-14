#' @title Nonlinear relative weight analysis with residualization
#'
#' @description Method to detecting influential variables in nonlinear model
#'   with residualized relative weight analysis.
#'

#'
#' @param response,free,fixed,control Character or character vectors with names
#'   for each parameter. Used in \code{residualrwa} to determine how to treat
#'   each variable. \code{response} is the model output. \code{free} are
#'   variables that can be included or removed (even in interactions) in all the
#'   models without restrictions. \code{fixed} are persistent variables across
#'   all the models, however they can be included/removed on interactions.
#'   \code{control} they only act as main effects and never are used for
#'   interactions.
#' @param data A \code{data.frame}.
#' @param family A string with default \code{"gaussian"}. It specifies the model
#'   link function. Accepted values are \code{"gaussian"} or \code{"binomial"}.
#' @param include_interactions A boolean with default \code{FALSE}. Determine if
#'   the model should calculate all the pairwise interactions between variables.
#'   It uses the names in the parameters \code{free} and \code{fixed}.
#' @param name_free,name_fixed,name_control,name_interactions A string for type
#'   of variable with defaults \code{"Free"}, \code{"Fixed"}, \code{"Control"}
#'   and \code{"Interaction"}. Names used to label the summary tables
#' @param boot_ci A boolean with default \code{FALSE}. Determine if a bootstrap
#'   procedure should be used to estimate confidence intervals for the weights.
#'   Defaults to \code{FALSE}.
#' @param n_boot A numeric with default \code{100}. Number of bootstrap samples
#'   used to estimate the confidence intervals.
#' @param mc_cores A numeric with default \code{1}. Number of cores used when
#'   performing the bootstrap samples.
#'
#' @return A residualrwa object with this structure:
#'
#' \describe{
#'
#' \item{summary}{a \code{data.frame} with the consolidated sum of relative
#' weights according to each type of component; free, fixed, control and
#' interactions.}
#'
#' \item{data_frame}{a \code{data.frame} with the individual relative weights
#' for each variable. The \code{data.frame} has columns: \code{variable},
#' \code{weight} and \code{type}. If \code{boot_ci = TRUE}, two optional columns
#' \code{ci_low} and \code{ci_up} indicating the lower and upper confidence
#' intervals for each variable respectively.}
#'
#' \item{model}{a \code{\link[rms]{Glm}} object with the final model used to
#' estimate the relative weights.}
#'
#' \item{data}{original \code{data.frame}.}
#'
#' \item{variables}{a list with character vectors \code{free}, \code{fixed},
#' \code{control} and \code{interactions} with the final variables used in the
#' model.}
#'
#' \item{boot_ci}{a boolean indicating if the bootstrap procedure was used.}
#'
#' }
#'
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
#' @export
residualrwa <- function(response,
                        control = NULL,
                        fixed = NULL,
                        free,
                        data,
                        family = "gaussian",
                        include_interactions = FALSE,
                        name_control = "Control",
                        name_fixed = "Fixed",
                        name_free = "Free",
                        name_interactions = "Interactions",
                        boot_ci = FALSE,
                        n_boot = 100,
                        mc_cores = 1) {
  out <- list(call = match.call())
  ## Setting NULL variables to pacify R checks.
  variable <- weight <- NULL

  if (!is.data.frame(data)) {
    stop("The parameter 'data' must be a data.frame")
  }

  if (!is.character(family)) {
    stop("Parameter family must be a string: 'gaussian' or 'binomial'.")
  }

  if (length(family) > 1) {
    stop("Parameter family must be of length 1.")
  }

  if (!(family %in% c("gaussian", "binomial"))) {
    stop(paste0("Parameter family = ", family, ", is not recognized."))
  }

  if (is.character(control) && length(control) == 0) {
    stop("control parameter must be a character vector or NULL")
  }

  if (is.character(fixed) && length(fixed) == 0) {
    stop("fixed parameter must be a character vector or NULL")
  }

  control_reorder <- extract_column_names(data, control)
  fixed_reorder <- extract_column_names(data, fixed)
  free_reorder <- extract_column_names(data, free)

  data <- data[, c(response, control_reorder, fixed_reorder, free_reorder)]


  out_residualrwa <- estimate_residualrwa(
    response = response,
    control = control,
    fixed = fixed,
    free = free,
    data = data,
    family = family,
    include_interactions = include_interactions,
    name_control = name_control,
    name_fixed = name_fixed,
    name_free = name_free,
    name_interactions = name_interactions
  )

  out <- append(out, out_residualrwa)

  if (boot_ci) {
    cat("\nEstimating bootstrap confidence intervals.")

    ## For debuggin purposes
    # rwa_boot <- lapply(

    rwa_boot <- pbmcapply::pbmclapply(
      mc.cores = mc_cores,
      X = 1:n_boot,
      FUN = function(i) {
        data_boot <- data[sample(nrow(data), nrow(data), replace = TRUE), ]

        run_rwa <- estimate_residualrwa(
          response = response,
          control = control,
          fixed = fixed,
          free = free,
          data = data_boot,
          family = family,
          include_interactions = include_interactions,
          name_control = name_control,
          name_fixed = name_fixed,
          name_free = name_free,
          name_interactions = name_interactions
        )

        df_out <- data.frame(
          cbind(run_rwa$data_frame,
            n_boot = rep(i, nrow(run_rwa$data_frame))
          )
        )
        return(df_out)
      }
    )
    cat("\n")

    out_boot <- dplyr::bind_rows(rwa_boot)
    out_boot <- tidyr::complete(data = out_boot, variable, n_boot)
    out_boot$weight <- ifelse(is.na(out_boot$weight), 0, out_boot$weight)

    out_boot_ci <- dplyr::group_by(.data = out_boot, variable)
    out_boot_ci <- dplyr::summarise(
      .data = out_boot_ci,
      ci_low = stats::quantile(weight, 0.025),
      ci_up = stats::quantile(weight, 0.975)
    )

    out$data_frame <- suppressMessages(
      dplyr::left_join(out$data_frame, out_boot_ci)
    )

    out <- append(out, list(
      data_frame_boot = out_boot,
      boot_ci = TRUE
    ))
  } else {
    out <- append(out, list(boot_ci = FALSE))
  }

  class(out) <- "residualrwa"
  return(out)
}

estimate_residualrwa <- function(response,
                                 control,
                                 fixed,
                                 free,
                                 data,
                                 family,
                                 include_interactions,
                                 name_control,
                                 name_fixed,
                                 name_free,
                                 name_interactions) {
  ## Setting NULL variables to pacify R checks.
  type <- weight <- NULL

  ## Core function starts here.
  control <- stringr::str_remove(control, "\\s")
  fixed <- stringr::str_remove(fixed, "\\s")
  free <- stringr::str_remove(free, "\\s")

  formula_base_model <- stats::formula(
    paste0(
      response,
      "~ 1 ",
      ifelse(length(control) != 0, paste0("+", paste0(
        control,
        collapse = "+"
      )), ""),
      ifelse(length(fixed) != 0,
        paste0("+", paste0(fixed, collapse = "+")), ""
      ),
      ifelse(length(free) != 0,
        paste0("+", paste0(free, collapse = "+")), ""
      )
    )
  )

  base_model <- rms::Glm(
    formula = formula_base_model,
    data = data,
    model = TRUE,
    maxit = 1000,
    family = family,
    x = TRUE,
    y = TRUE
  )

  interaction_model <- include_interactions_fn(
    formula = formula_base_model,
    response = response,
    data = data,
    control = control,
    fixed = fixed,
    free = free,
    include_interactions = include_interactions,
    family = family
  )

  names_in_model <- stats::terms(base_model$formula)
  names_in_model <- attr(names_in_model, "term.labels")

  model <- interaction_model$final_model

  x <- consolidate_design_matrix(
    model = model,
    interactions = include_interactions,
    family = family
  )
  x <- scale(x)
  y <- data[, response]

  base_names <- colnames(interaction_model$final_model$x)
  base_names <- stringr::str_replace(base_names, "\\[1\\]", "")
  base_names <- stringr::str_replace(base_names, "\\=.*", "")

  rwa_values <- estimate_rwa(x, y, data, family = family)

  columns_names <- colnames(x)
  columns_names <- stringr::str_remove(columns_names, "\\s")
  columns_names <- stringr::str_replace(columns_names, "\\[1\\]", "")
  columns_names <- stringr::str_replace(columns_names, "\\=.*", "")

  df_rwa_summary <- data.frame(
    variable = columns_names,
    weight = rwa_values$prop_weights
  )

  rownames(df_rwa_summary) <- columns_names

  cols_with_control <- columns_names[columns_names %in% control]

  cols_with_fixed <- columns_names[columns_names %in% fixed]

  cols_with_free <- columns_names[columns_names %in% free]

  cols_with_interactions <- columns_names[
    !(columns_names %in% control) &
      !(columns_names %in% fixed) &
      !(columns_names %in% free)
  ]

  df_rwa_summary[cols_with_control, "type"] <- name_control
  df_rwa_summary[cols_with_fixed, "type"] <- name_fixed
  df_rwa_summary[cols_with_free, "type"] <- name_free
  df_rwa_summary[cols_with_interactions, "type"] <- name_interactions

  resume <- df_rwa_summary
  resume <- dplyr::group_by(resume, type)
  resume <- dplyr::summarise(resume, total_rsq = sum(weight))

  out <- list(
    summary = resume,
    data_frame = df_rwa_summary,
    model = interaction_model$final_model,
    data = data,
    variables = list(
      free = df_rwa_summary[cols_with_free, "variable"],
      fixed = df_rwa_summary[cols_with_fixed, "variable"],
      control = df_rwa_summary[cols_with_control, "variable"],
      interactions = df_rwa_summary[cols_with_interactions, "variable"]
    )
  )
  return(out)
}

include_interactions_fn <- function(formula,
                                    response,
                                    data,
                                    control,
                                    fixed,
                                    free,
                                    include_interactions,
                                    family) {
  if (include_interactions) {
    combinations <- utils::combn(c(fixed, free), 2)

    idxcombi <- sapply(
      X = seq_along(combinations[1, ]),
      FUN = function(i) {
        all(
          stringr::str_detect(combinations[, i], "\\A1\\z", negate = TRUE)
        )
      }
    )

    combinations <- combinations[, idxcombi]

    if (is.matrix(combinations)) {
      free_interactions <- apply(combinations, 2, paste0, collapse = "%ia%")
    } else {
      free_interactions <- paste0(combinations, collapse = "%ia%")
    }

    idxia <- stringr::str_detect(free_interactions, "rcs")

    free_interactions <- c(
      free_interactions[idxia],
      stringr::str_replace(free_interactions[!idxia], "%ia%", "*")
    )

    frm_full <- stats::formula(
      paste0(
        response,
        "~",
        paste0(c(control, fixed), collapse = "+"),
        ifelse(length(control) == 0 & length(fixed) == 0, "", "+"),
        paste0(free, collapse = "+"),
        "+",
        paste0(free_interactions, collapse = "+")
      )
    )
  } else {
    frm_full <- stats::formula(
      paste0(
        response,
        "~ 1 ",
        ifelse(length(control) != 0,
          paste0("+", paste0(control, collapse = "+")), ""
        ),
        ifelse(length(fixed) != 0,
          paste0("+", paste0(fixed, collapse = "+")), ""
        ),
        ifelse(length(free) != 0, paste0("+", paste0(free, collapse = "+")), "")
      )
    )
  }

  if (family == "gaussian") {
    full_model <- rms::ols(
      formula = frm_full,
      data = data,
      model = TRUE,
      x = TRUE,
      y = TRUE
    )
  } else if (family == "binomial") {
    full_model <- rms::lrm(
      formula = frm_full,
      data = data,
      model = TRUE,
      x = TRUE,
      y = TRUE
    )
  }

  names_factors <- stringr::str_remove(
    names(full_model$model)[-1], "\\s"
  )

  idx_control <- which(names_factors == control)
  idx_fixed <- which(names_factors == fixed)

  force_control <- unlist(full_model$assign[idx_control])
  force_fixed <- unlist(full_model$assign[idx_fixed])

  stepwise_model <- rms::fastbw(
    fit = full_model,
    type = "individual",
    force = c(force_control, force_fixed),
    k.aic = log(nrow(data))
  )

  selected_vars <- names_factors[stepwise_model$factors.kept]

  selected_main_vars <- stringr::str_split(
    selected_vars, "( ?%ia% ?|\\*)",
    simplify = TRUE
  )

  selected_main_vars <- as.character(selected_main_vars)
  selected_main_vars <- stringr::str_remove(selected_main_vars, "\\s")
  selected_main_vars <- selected_main_vars[nchar(selected_main_vars) != 0]

  selected_vars <- unique(c(selected_main_vars, selected_vars))

  frm_selected <- stats::formula(
    paste0(
      response, "~",
      paste0(c(control, fixed), collapse = "+"),
      ifelse(length(control) == 0 & length(fixed) == 0, "", "+"),
      ifelse(length(selected_vars) == 0, "1", ""),
      paste0(selected_vars, collapse = "+")
    )
  )

  final_model <- rms::Glm(
    formula = frm_selected,
    data = data,
    family = family,
    model = TRUE,
    x = TRUE,
    y = TRUE
  )

  idx <- stringr::str_detect(
    attr(stats::terms(final_model), "term.labels"), "%ia%|\\*|:"
  )

  interactions <- attr(stats::terms(final_model), "term.labels")[idx]
  interactions <- stringr::str_remove(interactions, "\\s")

  return(
    list(
      control = control,
      fixed = fixed,
      free = unique(selected_main_vars),
      interactions = interactions,
      final_model = final_model
    )
  )
}
