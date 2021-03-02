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
#'
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
#'
#'

modelrwa <-
  function(response.name,
           control = NULL,
           fixed = NULL,
           free,
           data,
           family = stats::gaussian(),
           include.interactions = FALSE,
           alpha = 0.1,
           method = c("aic", "p"),
           name.control = "Control variables",
           name.fixed = "Fixed variables",
           name.free = "Free variables",
           name.interactions = "Interactions") {
    if (!is.data.frame(data))
      stop("The parameter 'data' must be a data.frame")

    # Reorder the dataframe

    control_reorder <-   extract_column_names(data, control)
    fixed_reorder <-   extract_column_names(data, fixed)
    free_reorder <-   extract_column_names(data, free)

    data <-
      data[, c(response.name,
               control_reorder,
               fixed_reorder,
               free_reorder)]

    if (is.character(control) & length(control) == 0)
    {
      stop("control parameter must be a character vector or NULL")
    }

    if (is.character(fixed) & length(fixed) == 0)
    {
      stop("fixed parameter must be a character vector or NULL")
    }

    if (is.null(control))
      control <- 1

    if (is.null(fixed))
      fixed <- 1

    control <- stringr::str_remove(control, "\\s")
    fixed <- stringr::str_remove(fixed, "\\s")
    free <- stringr::str_remove(free, "\\s")

    # pos.control <- which(colnames(data) %in% control_reorder )
    # pos.fixed <- which(colnames(data) %in% fixed_reorder )
    # pos.free <- which(colnames(data) %in% free_reorder )

    formula_base_model <-
      stats::formula(paste0(
        response.name,
        "~",
        paste0(control, collapse = "+"),
        "+",
        paste0(fixed, collapse = "+"),
        "+",
        paste0(free, collapse = "+")
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
      response.name = response.name,
      include.interactions = include.interactions,
      alpha = alpha,
      method = method,
      family = family
    )

    names_in_model <- stats::terms(base_model$formula)
    names_in_model <- attr(names_in_model, "term.labels")

    if (include.interactions &
        length(interaction_model$interactions) != 0) {
      interaction_names <- interaction_model$interactions
      individual_var_names <-
        stringr::str_split(interaction_names , "( ?%ia% ?|\\*|:)", simplify = TRUE)

      interaction_names <-
        stringr::str_replace(interaction_names, ":", "*")


      residualized_var_list <- NULL

      for (i in seq_along(interaction_names)) {
        w <-
          stats::as.formula(
            paste0(
              interaction_names[i],
              "~",
              "-1",
              sep = "+",
              individual_var_names[i, 1],
              sep = "+",
              individual_var_names[i, 2]
            )
          )

        residualized_model <- stats::lm(w, data)

        residualized_variable <-  stats::resid(residualized_model)



        if (!is.matrix(residualized_variable)) {
          residualized_variable <- unclass(residualized_variable)
          residualized_variable <- data.frame(residualized_variable)
          colnames(residualized_variable) <-
            colnames(residualized_model$model)[1]
        } else {
          # residualized_variable <- rowSums(residualized_variable)
          residualized_variable <-
            as.data.frame.matrix(residualized_variable)
        }


        residualized_var_list[[interaction_names[i]]] <-
          residualized_variable
      }
      names(residualized_var_list) <- NULL



      if (is.null(residualized_var_list)) {
        X <- as.data.frame(base_model$x)
      } else {
        is_main_effect <-
          !(attr(
            stats::terms(interaction_model$final_model$formula),
            "order"
          ) == 2)
        X <-
          cbind(as.data.frame(interaction_model$final_model$x[, is_main_effect]),
                residualized_var_list)
      }

      Y <- data[, response.name]
    } else {
      X <- interaction_model$final_model$x
      Y <- data[, response.name]
    }


    base_names <- colnames(interaction_model$final_model$x)
    base_names <- stringr::str_replace(base_names, "\\[1\\]", "")
    base_names <- stringr::str_replace(base_names, "\\=.*", "")



    (rwa_values <- RWA(X, Y, data))

    names_in_model_with_interactions <- colnames(X)
    names_in_model_with_interactions <-
      stringr::str_replace(names_in_model_with_interactions, "\\[1\\]", "")
    names_in_model_with_interactions <-
      stringr::str_replace(names_in_model_with_interactions, "\\=.*", "")

    df_rwa_summary <-
      data.frame(x = names_in_model_with_interactions,
                 y = rwa_values$PropWeights)
    rownames(df_rwa_summary) <- names_in_model_with_interactions



    cols.with.control <-
      names_in_model_with_interactions[names_in_model_with_interactions %in% control]

    cols.with.fixed <-
      names_in_model_with_interactions[names_in_model_with_interactions %in% fixed]

    cols.with.free <-
      names_in_model_with_interactions[(names_in_model_with_interactions %in% free)]

    cols.with.interactions <-
      names_in_model_with_interactions[!(names_in_model_with_interactions %in% cols.with.control) &
                                         !(names_in_model_with_interactions %in% cols.with.fixed) &
                                         !(names_in_model_with_interactions %in% cols.with.free)]


    df_rwa_summary[cols.with.control, "Variable_Type"] <-
      name.control

    df_rwa_summary[cols.with.fixed, "Variable_Type"] <- name.fixed

    df_rwa_summary[cols.with.free, "Variable_Type"] <-
      name.free

    df_rwa_summary[cols.with.interactions, "Variable_Type"] <-
      name.interactions


    resume <- df_rwa_summary
    resume <- dplyr::group_by(resume, "Variable_Type")
    resume <- dplyr::summarise(resume, "Effects_sum" = sum(y))


    out <- list(
      summary = resume,
      data_frame = df_rwa_summary,
      model = interaction_model$final_model,
      rwa_model = rwa_values
    )

    class(out) <- "modelrwa"

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
include_interactions <-
  function(formula = NULL,
           response.name = NULL,
           data = NULL,
           control = control,
           fixed = fixed,
           free = free,
           include.interactions = FALSE,
           family,
           alpha = 0.01,
           method = c("aic", "p")) {
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
        sapply(seq_along(combinations[1,]),  function(i)
          all(
            stringr::str_detect(combinations[, i], "\\A1\\z", negate = TRUE)
          ))

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
          "~",
          paste0(c(control, fixed), collapse = "+"),
          ifelse(length(control) == 0 &
                   length(fixed) == 0, "", "+"),
          paste0(free, collapse = "+")
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
          "~",
          paste0(c(control, fixed), collapse = "+"),
          ifelse(length(control) == 0 &
                   length(fixed) == 0, "", "+"),
          paste0(free, collapse = "+")
        ))

      frm_full <- frm_base
    }

    base_model <-
      rms::Glm(
        formula = frm_base,
        data = data,
        family = family,
        model = TRUE,
        x = TRUE,
        y = TRUE
      )



    full_model <-
      rms::Glm(
        formula = frm_full,
        data = data,
        family = family,
        model = TRUE,
        x = TRUE,
        y = TRUE,
        maxit = 500
      )

    force_control <-
      which(colnames(full_model$x) %in% extract_column_names(data, control))

    force_fixed <-
      which(colnames(full_model$x) %in% extract_column_names(data, fixed))


    stepwise_model <-
      rms::fastbw(
        full_model,
        rule = method,
        type = "individual",
        sls = alpha,
        force = c(force_control, force_fixed)
      )

    selected_main_vars <-
      stringr::str_split(stepwise_model$names.kept , "( ?%ia% ?|\\*)", simplify = TRUE)
    selected_main_vars  <- as.character(selected_main_vars)
    selected_main_vars <-
      stringr::str_remove(selected_main_vars, "\\s")
    selected_main_vars <-
      selected_main_vars[nchar(selected_main_vars) != 0]


    selected_vars <-
      unique(c(selected_main_vars, stepwise_model$names.kept))




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

    idx <- attr(terms(final_model), "order") == 2
    interactions <- attr(terms(final_model), "term.labels")[idx]

    # interactions <-
    #   stringr::str_replace(interactions, ":", "*")

    interactions <- stringr::str_remove(interactions, "\\s")

    return(
      list(
        control = control,
        fixed = fixed,
        free = free,
        interactions = interactions,
        final_model = final_model
      )
    )
  }
