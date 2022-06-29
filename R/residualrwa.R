#' #' Detect influential variables in model with Relative Weight Analysis
#'
#' @param response.name name for independent variable
#' @param control character vector with names for control variables
#' @param fixed character vector with names for fixed variables
#' @param free character vector with names for free variables
#' @param data object data.frame with data
#' @param family type of regression (see \code{\link[stats]{family}})
#' @param include.interactions A boolean indicating if the model should calculate all the pairwise interactions between variables. It uses the names in the parameters \code{fixed} and \code{variables}.
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
                        name_control = "Control",
                        name_fixed = "Fixed",
                        name_free = "Free",
                        name_interactions = "Interactions",
                        boot_ci = FALSE,
                        n_boot = 100,
                        mc_cores = 1,
                        verbose = FALSE) {
  if (!is.data.frame(data)) {
    stop("The parameter 'data' must be a data.frame")
  } else {
    # Reorder the dataframe

    control_reorder <- extract_column_names(data, control)
    fixed_reorder <- extract_column_names(data, fixed)
    free_reorder <- extract_column_names(data, free)

    data <- data[, c(
      response_name,
      control_reorder,
      fixed_reorder,
      free_reorder
    )]
  }

  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (is.null(family$family) ||
    !(family$family %in% c("gaussian", "binomial"))) {
    print(family)
    stop("'family' not recognized for residualrwa")
  }


  if (is.character(control) && length(control) == 0) {
    stop("control parameter must be a character vector or NULL")
  }

  if (is.character(fixed) && length(fixed) == 0) {
    stop("fixed parameter must be a character vector or NULL")
  }


  out_residualrwa <- estimate_residualrwa(
    response_name,
    control,
    fixed,
    free,
    data,
    family,
    include_interactions,
    name_control,
    name_fixed,
    name_free,
    name_interactions,
    verbose
  )

  out <- out_residualrwa

  if (boot_ci) {
    message("\n")
    rwa_boot <- parallel::mclapply(
      mc.cores = mc_cores,
      X = 1:n_boot,
      function(i) {
        message(paste0("Boot sample #", i))
        data_boot <- data[sample(nrow(data), nrow(data), replace = TRUE), ]

        run_rwa <- estimate_residualrwa(
          response_name,
          control,
          fixed,
          free,
          data = data_boot,
          family,
          include_interactions,
          name_control,
          name_fixed,
          name_free,
          name_interactions,
          verbose
        )

        df_out <- data.frame(
          cbind(run_rwa$data_frame,
            n_boot = rep(i, nrow(run_rwa$data_frame))
          )
        )
        return(df_out)
      }
    )

    out_boot <- dplyr::bind_rows(rwa_boot)

    out_boot <- tidyr::complete(data = out_boot, variable, n_boot)

    out_boot$weight <- ifelse(is.na(out_boot$weight), 0, out_boot$weight)

    out <- append(out, list(data_frame_boot = out_boot, boot_ci = TRUE))
  } else {
    out <- append(out, list(boot_ci = FALSE))
  }

  class(out) <- "residualrwa"

  return(out)
}


estimate_residualrwa <- function(response_name,
                                 control,
                                 fixed,
                                 free,
                                 data,
                                 family,
                                 include_interactions,
                                 name_control,
                                 name_fixed,
                                 name_free,
                                 name_interactions,
                                 verbose) {
  control <- stringr::str_remove(control, "\\s")
  fixed <- stringr::str_remove(fixed, "\\s")
  free <- stringr::str_remove(free, "\\s")


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

  interaction_model <- include_interactions_fn(
    formula = formula_base_model,
    response_name = response_name,
    data = data,
    control = control,
    fixed = fixed,
    free = free,
    include_interactions = include_interactions,
    family = family,
    verbose = verbose
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
  y <- data[, response_name]

  base_names <- colnames(interaction_model$final_model$x)
  base_names <- stringr::str_replace(base_names, "\\[1\\]", "")
  base_names <- stringr::str_replace(base_names, "\\=.*", "")

  rwa_values <- estimate_rwa(x, y, data, family = family)

  columns_names <- colnames(x)
  columns_names <- stringr::str_remove(columns_names, "\\s")
  columns_names <- stringr::str_replace(columns_names, "\\[1\\]", "")
  columns_names <- stringr::str_replace(columns_names, "\\=.*", "")

  df_rwa_summary <-
    data.frame(
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
  resume <- dplyr::summarise(resume, Effects_sum = sum(weight))

  out <- list(
    summary = resume,
    data_frame = df_rwa_summary,
    model = interaction_model$final_model,
    rwa_model = rwa_values
  )
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
#'
#' @return a list with the control, fixed, free and interactions used in the model. Also it returns the the fitted model.
#'
#' @keywords internal
include_interactions_fn <- function(formula = NULL,
                                    response_name = NULL,
                                    data = NULL,
                                    control = control,
                                    fixed = fixed,
                                    free = free,
                                    include_interactions = FALSE,
                                    family,
                                    verbose = FALSE) {
  if (include_interactions) {
    combinations <- utils::combn(c(fixed, free), 2)

    idxcombi <-
      sapply(seq_along(combinations[1, ]), function(i) {
        all(
          stringr::str_detect(combinations[, i], "\\A1\\z", negate = TRUE)
        )
      })

    combinations <- combinations[, idxcombi]

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
        response_name,
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
        response_name,
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
        response_name,
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

  if (length(control) == 0 && length(fixed) == 0) {
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

  selected_main_vars <- stringr::str_split(selected_vars,
    "( ?%ia% ?|\\*)",
    simplify = TRUE
  )
  selected_main_vars <- as.character(selected_main_vars)
  selected_main_vars <- stringr::str_remove(selected_main_vars, "\\s")
  selected_main_vars <- selected_main_vars[nchar(selected_main_vars) != 0]


  selected_vars <- unique(c(selected_main_vars, selected_vars))




  frm_selected <-
    stats::formula(paste0(
      response_name,
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


  idx <- stringr::str_detect(
    attr(
      terms(final_model),
      "term.labels"
    ),
    "%ia%|\\*|:"
  )
  interactions <- attr(terms(final_model), "term.labels")[idx]


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
