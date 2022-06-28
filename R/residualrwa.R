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
#' Y <- plogis(X2^3 + 10 * X1 * X2) > 0.5
#'
#' data <- as.data.frame(cbind(Y, X1, X2, X3))
#'
#' ex1 <- modelrwa(
#'   response.name = "Y",
#'   control = NULL,
#'   fixed = NULL,
#'   variables = c("X1", "X2", "X3"),
#'   data = data,
#'   family = binomial,
#'   include.interactions = TRUE
#' )
#'
residualrwa <- function(response_name,
                        control = NULL,
                        fixed = NULL,
                        free,
                        data,
                        family = stats::gaussian(),
                        include_interactions = FALSE,
                        alpha = 0.1,
                        method = c("aic", "p"),
                        name_control = "Control",
                        name_fixed = "Fixed",
                        name_free = "Free",
                        name_interactions = "Interactions",
                        verbose = FALSE) {
  if (!is.data.frame(data)) {
    stop("The parameter 'data' must be a data.frame")
  } else {
    # Reorder the dataframe

    control_reorder <- extract_column_names(data, control)
    fixed_reorder <- extract_column_names(data, fixed)
    free_reorder <- extract_column_names(data, free)
  }

  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (is.null(family$family) |
    !(family$family %in% c("gaussian", "binomial"))) {
    print(family)
    stop("'family' not recognized for residualrwa")
  }


  data <-
    data[, c(
      response_name,
      control_reorder,
      fixed_reorder,
      free_reorder
    )]

  if (is.character(control) & length(control) == 0) {
    stop("control parameter must be a character vector or NULL")
  }

  if (is.character(fixed) & length(fixed) == 0) {
    stop("fixed parameter must be a character vector or NULL")
  }

  # if (is.null(control))
  #   control <- 1
  #
  # if (is.null(fixed))
  #   fixed <- 1

  control <- stringr::str_remove(control, "\\s")
  fixed <- stringr::str_remove(fixed, "\\s")
  free <- stringr::str_remove(free, "\\s")

  # pos.control <- which(colnames(data) %in% control_reorder )
  # pos.fixed <- which(colnames(data) %in% fixed_reorder )
  # pos.free <- which(colnames(data) %in% free_reorder )

  formula_base_model <-
    stats::formula(paste0(
      response_name,
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
    ))

  base_model <-
    rms::Glm(
      formula = formula_base_model,
      data = data,
      model = TRUE,
      maxit = 1000,
      family = family,
      x = TRUE,
      y = TRUE
    )



  interaction_model <- include_interactions(
    formula = formula_base_model,
    data = data,
    control = control,
    fixed = fixed,
    free = free,
    # pos.fixed = pos.fixed,
    # pos.free = pos.free,
    # pos.control = pos.control,
    response.name = response_name,
    include.interactions = include_interactions,
    alpha = alpha,
    method = method,
    family = family,
    verbose = verbose
  )

  names_in_model <- stats::terms(base_model$formula)
  names_in_model <- attr(names_in_model, "term.labels")




  # Old residualization implementation
  # if (include.interactions &
  #     length(interaction_model$interactions) != 0) {
  #   interaction_names <- interaction_model$interactions
  #   individual_var_names <-
  #     stringr::str_split(interaction_names , "( ?%ia% ?|\\*|:)", simplify = TRUE)
  #
  #   interaction_names <-
  #     stringr::str_replace(interaction_names, ":", "*")
  #
  #
  #   residualized_var_list <- NULL
  #
  #   for (i in seq_along(interaction_names)) {
  #     w <-
  #       stats::as.formula(
  #         paste0(
  #           interaction_names[i],
  #           "~",
  #           "-1",
  #           sep = "+",
  #           individual_var_names[i, 1],
  #           sep = "+",
  #           individual_var_names[i, 2]
  #         )
  #       )
  #
  #
  #
  #     residualized_model <- stats::lm(w, data)
  #
  #     residualized_variable <-  stats::resid(residualized_model)
  #
  #
  #     if (interactions) {
  #       col.interactions <-
  #         which(model$Design$assume == "interaction")
  #       # min.interaction <- min(col.interactions)
  #       # design.interaction <- model$Design$interactions
  #
  #
  #
  #
  #       for (k in col.interactions) {
  #         xx <- str_split(colnames(model$x)[idx[[k]] - 1], " \\* |:")
  #
  #         if (length(xx) == 1) {
  #           ll <- lm(model$x[, idx[[k]] - 1] ~ -1 + model$x[, xx[[1]]])
  #           XDesign[, idx[[k]] - 1] <- resid(ll)
  #
  #         } else{
  #           for (i in seq_along(xx)) {
  #             ll <- lm(model$x[, idx[[k]] - 1][, i] ~ -1 + model$x[, xx[[i]]])
  #             XDesign[, idx[[k]] - 1][, i] <- resid(ll)
  #           }
  #         }
  #
  #
  #       }
  #
  #
  #     }
  #
  #
  #     # New implementation of residualization
  #     # if (!is.matrix(residualized_variable)) {
  #     #   residualized_variable <- unclass(residualized_variable)
  #     #   residualized_variable <- data.frame(residualized_variable)
  #     #   colnames(residualized_variable) <-
  #     #     colnames(residualized_model$model)[1]
  #     # } else {
  #     #   # residualized_variable <- rowSums(residualized_variable)
  #     #   residualized_variable <-
  #     #     as.data.frame.matrix(residualized_variable)
  #     # }
  #
  #
  #     residualized_var_list[[interaction_names[i]]] <-
  #       residualized_variable
  #   }
  #   names(residualized_var_list) <- NULL
  #
  #
  #
  #   if (is.null(residualized_var_list)) {
  #     X <- as.data.frame(base_model$x)
  #   } else {
  #     names <- colnames(interaction_model$final_model$x)
  #     is_main_effect <- !stringr::str_detect(names, "\\*")
  #
  #     X <-
  #       cbind(as.data.frame(interaction_model$final_model$x[, is_main_effect]),
  #             residualized_var_list)
  #   }
  #
  #   Y <- data[, response.name]
  # } else {
  #   X <- interaction_model$final_model$x
  #   Y <- data[, response.name]
  #   is_main_effect <- !logical(ncol(X))
  # }


  model <- interaction_model$final_model

  X <- extract_X_RWA(
    model = model,
    interactions = include_interactions,
    family = family
  )
  X <- scale(X)
  Y <- data[, response_name]

  base_names <- colnames(interaction_model$final_model$x)
  base_names <- stringr::str_replace(base_names, "\\[1\\]", "")
  base_names <- stringr::str_replace(base_names, "\\=.*", "")



  rwa_values <- RWA(X, Y, data, family = family)

  names_in_model_with_interactions <- colnames(X)
  names_in_model_with_interactions <- stringr::str_remove(names_in_model_with_interactions, "\\s")
  names_in_model_with_interactions <-
    stringr::str_replace(names_in_model_with_interactions, "\\[1\\]", "")
  names_in_model_with_interactions <-
    stringr::str_replace(names_in_model_with_interactions, "\\=.*", "")

  df_rwa_summary <-
    data.frame(
      Variable = names_in_model_with_interactions,
      Weight = rwa_values$PropWeights
    )
  rownames(df_rwa_summary) <- names_in_model_with_interactions

  #
  #     control_columns_names <-
  #       as.character(sapply(control, function(x)
  #         attr(eval(parse(
  #           text = x
  #         )), "colnames")))
  #
  #     fixed_columns_names <-
  #       as.character(sapply(fixed, function(x)
  #         attr(eval(parse(
  #           text = x
  #         )), "colnames")))
  #
  #     free_columns_names <-
  #       as.character(sapply(free, function(x)
  #         attr(eval(parse(
  #           text = x
  #         )), "colnames")))


  cols.with.control <-
    names_in_model_with_interactions[names_in_model_with_interactions %in% control]

  cols.with.fixed <-
    names_in_model_with_interactions[names_in_model_with_interactions %in% fixed]

  cols.with.free <-
    names_in_model_with_interactions[names_in_model_with_interactions %in% free]

  cols.with.interactions <-
    names_in_model_with_interactions[!(names_in_model_with_interactions %in% control) &
      !(names_in_model_with_interactions %in% fixed) &
      !(names_in_model_with_interactions %in% free)]


  df_rwa_summary[cols.with.control, "Type"] <-
    name_control

  df_rwa_summary[cols.with.fixed, "Type"] <- name_fixed

  df_rwa_summary[cols.with.free, "Type"] <-
    name_free

  df_rwa_summary[cols.with.interactions, "Type"] <-
    name_interactions


  resume <- df_rwa_summary
  resume <- dplyr::group_by(resume, Type)
  resume <- dplyr::summarise(resume, Effects_sum = sum(Weight))


  out <- list(
    summary = resume,
    data_frame = df_rwa_summary,
    model = interaction_model$final_model,
    rwa_model = rwa_values
  )

  class(out) <- "residualrwa"

  return(out)
}





#' Include interactions in a model
#'
#' @param formula  an object of class formula with main effects.
#' @param data  a data.frame.
#' @param response.name output name of the data.
# @param pos.control,pos.fixed,pos.free position of the control, fixed and variables inputs in the data.frame.
#' @param include.interactions a boolean indicating if create or no a pairwise set of interactions.
#' @param family type of regression (see \code{\link[stats]{family}})
#' @param alpha level of significance use if \code{method = "p"}.
#' @param method method used to keep interactions if \code{include.interactions=TRUE}. See \code{\link[rms]{fastbw}}.
#'
#' @return a list with the control, fixed, free and interactions used in the model. Also it returns the the fitted model.
#'
#' @keywords internal
include_interactions <- function(formula = NULL,
                                 response.name = NULL,
                                 data = NULL,
                                 control = control,
                                 fixed = fixed,
                                 free = free,
                                 include.interactions = FALSE,
                                 family,
                                 alpha = 0.01,
                                 method = c("aic", "p"),
                                 verbose = FALSE) {
  # if (is.null(pos.fixed)) {
  #   fixed <- NULL
  #
  #   free <-
  #     attr(stats::terms(formula), "term.labels")
  # } else {
  #   # It doesnt assume any order
  #   idxcontrol <-
  #     attr(stats::terms(formula), "term.labels") %in% colnames(data)[pos.control]
  #
  #   idxfixed <-
  #     attr(stats::terms(formula), "term.labels") %in% colnames(data)[pos.fixed]
  #
  #
  #   idxvariable <-
  #     attr(stats::terms(formula), "term.labels") %in% colnames(data)[pos.free]
  #
  #   control <-
  #     attr(stats::terms(formula), "term.labels")[idxcontrol]
  #   # if (length(control) == 0) {
  #   #   control <- 1
  #   # }
  #
  #   fixed <-
  #     attr(stats::terms(formula), "term.labels")[idxfixed]
  #   # if (length(fixed) == 0) {
  #   #   fixed <- 1
  #   # }
  #
  #   free <-
  #     attr(stats::terms(formula), "term.labels")[idxvariable]
  # }

  if (include.interactions) {
    combinations <- utils::combn(c(fixed, free), 2)

    idxcombi <-
      sapply(seq_along(combinations[1, ]), function(i) {
        all(
          stringr::str_detect(combinations[, i], "\\A1\\z", negate = TRUE)
        )
      })

    combinations <- combinations[, idxcombi]

    # free_interactions <-
    #   apply(combinations, 2, paste0, collapse = "*")

    free_interactions <-
      apply(combinations, 2, paste0, collapse = "%ia%")

    idxia <-
      stringr::str_detect(free_interactions, "rcs")

    free_interactions <- c(
      free_interactions[idxia],
      stringr::str_replace(free_interactions[!idxia], "%ia%", "*")
    )


    frm_base <-
      stats::formula(paste0(
        response.name,
        "~ 1 ",
        ifelse(length(control) != 0, paste0("+", paste0(
          control,
          collapse = "+"
        )), ""),
        ifelse(length(fixed) != 0, paste0("+", paste0(
          fixed,
          collapse = "+"
        )), ""),
        ifelse(length(free) != 0, paste0("+", paste0(free, collapse = "+")), "")
      ))



    frm_full <-
      stats::formula(paste0(
        response.name,
        "~",
        paste0(c(control, fixed), collapse = "+"),
        ifelse(length(control) == 0 &
          length(fixed) == 0, "", "+"),
        paste0(free, collapse = "+"),
        "+",
        paste0(free_interactions, collapse = "+")
      ))
  } else {
    frm_base <-
      stats::formula(paste0(
        response.name,
        "~ 1 ",
        ifelse(length(control) != 0, paste0("+", paste0(
          control,
          collapse = "+"
        )), ""),
        ifelse(length(fixed) != 0, paste0("+", paste0(
          fixed,
          collapse = "+"
        )), ""),
        ifelse(length(free) != 0, paste0("+", paste0(free, collapse = "+")), "")
      ))

    frm_full <- frm_base
  }

  base_model <-
    stats::glm(
      formula = frm_base,
      data = data,
      family = family,
      model = TRUE,
      x = TRUE,
      y = TRUE
    )



  # full_model <-
  #   rms::Glm(
  #     formula = frm_full,
  #     data = data,
  #     family = family,
  #     model = TRUE,
  #     x = TRUE,
  #     y = TRUE,
  #     maxit = 500
  #   )

  # force_control <-
  #   which(colnames(full_model$x) %in% extract_column_names(data, control))
  #
  # force_fixed <-
  #   which(colnames(full_model$x) %in% extract_column_names(data, fixed))


  # force_control <- extract_column_names(data, control)
  # force_fixed <- extract_column_names(data, fixed)

  if (length(control) == 0 & length(fixed) == 0) {
    frm_lower <- formula("~1")
  } else {
    frm_lower <-
      formula(paste0("~", paste0(control, fixed, collapse = "+")))
  }

  stepwise_model <- MASS::stepAIC(base_model,
    scope = list(
      lower = frm_lower,
      upper = frm_full
    ), trace = verbose
  )

  selected_vars <- attr(stepwise_model$terms, "term.labels")


  # stepwise_model <-
  #   rms::fastbw(
  #     full_model,
  #     rule = method,
  #     type = "individual",
  #     sls = alpha,
  #     force = c(force_control, force_fixed)
  #   )
  #

  # selected_vars <-
  #   attr(terms(full_model$formula), "term.label")[stepwise_model$factors.kept]

  selected_main_vars <- stringr::str_split(selected_vars, "( ?%ia% ?|\\*)", simplify = TRUE)
  selected_main_vars <- as.character(selected_main_vars)
  selected_main_vars <- stringr::str_remove(selected_main_vars, "\\s")
  selected_main_vars <- selected_main_vars[nchar(selected_main_vars) != 0]


  selected_vars <- unique(c(selected_main_vars, selected_vars))




  # selected_vars <-  stringr::str_replace(selected_vars, ":", "*")

  frm_selected <-
    stats::formula(paste0(
      response.name,
      "~",
      paste0(c(control, fixed), collapse = "+"),
      ifelse(length(control) == 0 & length(fixed) == 0, "", "+"),
      ifelse(length(selected_vars) == 0, "1", ""),
      paste0(selected_vars, collapse = "+")
    ))


  final_model <- rms::Glm(
    formula = frm_selected,
    data = data,
    family = family,
    model = TRUE,
    x = TRUE,
    y = TRUE
  )

  # terms_labels <-
  #   stringr::str_remove(attr(stats::terms(final_model), "term.labels"), "\\s")
  #
  # is_variable_interaction <-
  #   !(terms_labels %in% control |
  #       terms_labels %in% fixed |
  #       terms_labels %in% free)

  # interactions <-
  #   as.character(attr(stats::terms(final_model)[is_variable_interaction], "term.labels"))

  idx <- stringr::str_detect(
    attr(
      terms(final_model),
      "term.labels"
    ),
    "%ia%|\\*|:"
  )
  interactions <- attr(terms(final_model), "term.labels")[idx]

  # interactions <-
  #   stringr::str_replace(interactions, ":", "*")

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
