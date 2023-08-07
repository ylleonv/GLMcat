#' Fitting models for categorical responses
#' @description Estimate generalized linear models implemented under the unified
#' specification ( ratio,cdf,Z) where \code{ratio} represents the ratio of probabilities
#' (reference, cumulative, adjacent, or sequential), \code{cdf} the cumulative distribution function
#' for the linkage, and Z the design matrix which must be specified through the \code{parallel}
#' and the \code{threshold} arguments.
#' @title Generalized linear models for categorical responses
#' @rdname glmcat
#' @name glmcat
#' @param formula formula a symbolic description of the model to be fit. An expression of the form `y ~ predictors` is interpreted as a specification that the response `y` is modeled by a linear predictor specified by `predictors`.
#' @param ratio a string indicating the ratio (equivalently to the family) options are: reference, adjacent, cumulative and sequential.  It is mandatory for the user to specify the desired ratio option as there is no default value.
#' @param cdf The inverse distribution function to be used as part of the link function.
#'   - If the distribution has no parameters to specify, then it should be entered as a string indicating the name, e.g., `cdf = "normal"`. The default value is `cdf = "logistic"`.
#'   - If there are parameters to specify, then a list must be entered. For example, for Student's distribution: `cdf = list("student", df=2)`. For the non-central distribution of Student: `cdf = list("noncentralt", df=2, mu=1)`.
#' @param categories_order a character vector indicating the incremental order of the categories, e.g., `c("a", "b", "c")` for `a < b < c`. Alphabetical order is assumed by default. Order is relevant for adjacent, cumulative, and sequential ratio.
#' @param ref_category a string indicating the reference category. This option is suitable for models with reference ratio.
#' @param parallel a character vector indicating the name of the variables with a parallel effect. If a variable is categorical, specify the name and the level of the variable as a string, e.g., `"namelevel"`.
#' @param data a dataframe object in R, with the dependent variable as a factor.
#' @param threshold a restriction to impose on the thresholds. Options are: `standard`, `equidistant`, or `symmetric`. This is valid only for the cumulative ratio.
#' @param control a list of control parameters for the estimation algorithm.
#'   - `maxit`: The maximum number of iterations for the Fisher scoring algorithm.
#'   - `epsilon`: A double to change the convergence criterion of GLMcat models.
#'   - `beta_init`: An appropriately sized vector for the initial iteration of the algorithm.
#' @param normalization the quantile to use for the normalization of the estimated coefficients when the logistic distribution is used as the base cumulative distribution function.
#' @param na.action an argument to handle missing data. Available options are `na.omit`, `na.fail`, and `na.exclude`. It does not include the `na.pass` option.
#' @param find_nu a logical argument to indicate whether the user intends to utilize the Student CDF and seeks an optimization algorithm to identify an optimal degrees of freedom setting for the model.
#' @param ... additional arguments.
#'
#' @details This function fits generalized linear models for categorical responses using the unified specification framework introduced by Peyhardi, Trottier, and Guédon (2015).
#'
#' @export
#' @references
#'  Peyhardi J, Trottier C, Guédon Y (2015). “A new specification of generalized linear models
#'  for categorical responses.” \emph{Biometrika}, 102(4), 889–906. doi:10.1093/biomet/asv042.
#' @examples
#' data(DisturbedDreams)
#' ref_log_com <- glmcat(formula = Level ~ Age, data = DisturbedDreams,
#'     ref_category = "Very.severe",
#'     cdf = "logistic", ratio = "reference")
#' @seealso
#' \code{\link{summary.glmcat}}
#' @keywords generalized linear model, categorical variables
glmcat <-
  function(formula,
           data,
           ratio = c("reference", "cumulative", "sequential","adjacent"),
           # cdf = "logistic",
           cdf = list(),
           parallel = NA,
           categories_order = NA,
           ref_category = NA,
           threshold = c("standard", "symmetric", "equidistant"),
           control = list(),
           normalization = 1,
           na.action = "na.omit",
           find_nu = FALSE,
           # doFit = TRUE, na.action,
           # contrasts, model = TRUE,
           ...)
  {
    # Check if the ratio argument is missing
    if(missing(ratio)){
      stop("The ratio was not specified and is required. Please specify one among the options: cumulative, sequential, adjacent, or reference.")
    }

    # check_categorical <- is.factor(model.frame(formula = formula, data)[, 1])
    # if (!check_categorical) {
    #   warning("The response variable is not defined as a categorical variable.")
    # }

    # Check if the response variable is defined as a categorical variable
    check_ordered <- is.factor(model.frame(formula = formula, data)[,1])
    if ( check_ordered == FALSE ) { warning( "The response variable is not defined as a categorical variable" ) }

    # Check if the response variable is defined as an ordered variable
    check_ordered <- is.ordered(model.frame(formula = formula, data)[,1])
    if ( check_ordered == FALSE & ratio != "reference") { warning( "The response variable is not defined as an ordered variable. Recall that the the reference ratio is appropiate for nominal responses, while for ordinal responses the ratios to use are cumulative, sequential or adjacent." ) }
    if ( check_ordered == TRUE & ratio == "reference") { warning( "The response variable is defined as an ordered variable. Recall that the the reference ratio is appropiate for nominal responses, while for ordinal responses the ratios to use are cumulative, sequential or adjacent." ) }

    # Set the default value for cdf
    if(length(cdf)==0){cdf[[1]] = "logistic"}

    # Match the value of cdf to a valid option
    cdf[[1]] <- match.arg(cdf[[1]], c("logistic", "normal", "gumbel", "gompertz", "cauchy",
                                      "laplace", "student", "noncentralt"))
    ratio <- match.arg(ratio)

    # Set the default value for threshold
    threshold <- match.arg(threshold)
    contrasts <- NULL
    control <- do.call(control_glmcat, c(control, list(...)))

    # Default for reference ratio should be the complete design
    if(ratio == "reference" && is.na(parallel)){
      parallel = F
    }

    # Set the na.action based on the specified option
    na.action <- match.arg(na.action, c("na.omit", "na.fail", "na.exclude"))
    if(na.action == "na.omit"){
      data <- na.omit(data)
    }else if (na.action == "na.fail"){
      data <- na.fail(data)
    }else if (na.action == "na.exclude"){
      data <- na.exclude(data)
    }

    cdf_sel <- cdf

    # if(find_nu == TRUE) {
    #   # Estimate the models with Student link where ν = 1 and ν = 8
    #   cdf_1 <- list("student", 1)
    #   cdf_8 <- list("student", 8)
    #
    #   model_1 <- .GLMcat(formula = formula, data = data, ratio = ratio, cdf = cdf_1, parallel = parallel, categories_order = categories_order,
    #                      ref_category = ref_category, threshold = threshold , control = control, normalization = normalization)
    #
    #   model_8 <- .GLMcat(formula = formula, data = data, ratio = ratio, cdf = cdf_8, parallel = parallel, categories_order = categories_order,
    #                      ref_category = ref_category, threshold = threshold , control = control, normalization = normalization)
    #
    #
    #   # Check if lν=8 > lν=1
    #   if (model_8[["LogLikelihood"]] > model_1[["LogLikelihood"]]) {
    #     # Estimate the log-likelihood lp of a binary model with the probit link
    #     model_p <- .GLMcat(formula = formula, data = data, ratio = ratio, cdf = "normal", parallel = parallel, categories_order = categories_order,
    #                        ref_category = ref_category, threshold = threshold , control = control, normalization = normalization)
    #
    #     # Check if lp > lν=8
    #     if (model_p[["LogLikelihood"]] > model_8[["LogLikelihood"]]) {
    #       # Use the probit link
    #       cdf_sel <- "normal"
    #     } else {
    #       # Use the logit link
    #       cdf_sel <- "logistic"
    #     }
    #   } else {
    #     # Use optimize() to find the best ν ∈ (0.25, 1) of the Student CDF
    #     optimize_likelihood <- function(nu) {
    #       model_nu <- .GLMcat(formula = formula, data = data, ratio = ratio, cdf = list("student", nu), parallel = parallel, categories_order = categories_order,
    #                           ref_category = ref_category, threshold = threshold , control = control, normalization = normalization)
    #
    #       model_nu[["LogLikelihood"]]
    #     }
    #
    #     opt_result <- optimize(optimize_likelihood, interval = c(0.25, 1), maximum = TRUE)
    #     best_nu <- opt_result$maximum
    #
    #     # Use the best ν found
    #     cdf_sel <- list("student", best_nu)
    #   }
    # }

    # Use optimize() to find the best ν ∈ (0.25, 1) of the Student CDF
    optimize_likelihood <- function(nu) {
      model_nu <- .GLMcat(formula = formula, data = data, ratio = ratio, cdf = list("student", nu), parallel = parallel, categories_order = categories_order,
                          ref_category = ref_category, threshold = threshold , control = control, normalization = normalization)

      model_nu[["LogLikelihood"]]
    }

    if(find_nu == TRUE){
      opt_result <- optimize(optimize_likelihood, interval = c(0.25, 8), maximum = TRUE)
      best_nu <- opt_result$maximum
      # Use the best ν found
      cdf_sel <- list("student", best_nu)
    }

    cdf <- cdf_sel

    # cdf_sel <- cdf_sel1
    # Call the GLMcat C++ function
    fit_old <- .GLMcat(formula = formula, data = data, ratio = ratio, cdf = cdf_sel, parallel = parallel, categories_order = categories_order,
                       ref_category = ref_category, threshold = threshold , control = control, normalization = normalization)



    # Store the model frame and data in the fit_old object
    fit_old[["model"]] <- model.frame(formula = formula, data)
    fit_old[["data"]] <- data

    # Generate the table summary

    fit_old$table_summary <- table_summary(fit_old)

    fit_old <- fit_old[sort(names(fit_old))]
    class(fit_old) <- "glmcat"

    return(fit_old)
  }

#' Family of models for Discrete Choice
#' @description Family of models for Discrete Choice. Fits discrete choice models which require data in long form.
#' For each individual (or decision maker), there are multiple observations (rows),
#' one for each of the alternatives the individual could have chosen.
#' A group of observations of the same individual is a "case".
#' It is important to note that each case represents a single statistical observation
#' although it comprises multiple observations.
#'
#' @title Discrete Choice Models
#' @rdname discrete_cm
#' @name discrete_cm
#'
#' @param formula a symbolic description of the model to be fit.
#'   An expression of the form y ~ predictors is interpreted as a specification
#'   that the response y is modeled by a linear predictor specified symbolically by model.
#'   A particularity for the formula is that for the case-specific variables,
#'   the user can define a specific effect for a category (in the parameter `alternative_specific`).
#' @param case_id a string with the name of the column that identifies each case.
#' @param alternatives a string with the name of the column that identifies
#'   the vector of alternatives the individual could have chosen.
#' @param reference a string indicating the reference category.
#' @param alternative_specific a character vector with the name of the explanatory variables
#'   that are different for each case, these are the alternative-specific variables.
#'   By default, the case-specific variables are the explanatory variables
#'   that are not identified here but are part of the formula.
#' @param data a dataframe (in long format) object in R, with the dependent variable as a factor.
#' @param cdf a parameter specifying the inverse distribution function to be used as part of the link function.
#'   If the distribution has no parameters to specify, it should be entered as a string indicating the name.
#'   The default value is 'logistic'. If there are parameters to specify, a list must be entered.
#'   For example, for Student's distribution, it would be `list("student", df=2)`.
#'   For the non-central distribution of Student, it would be `list("noncentralt", df=2, mu=1)`.
#' @param intercept if set to "conditional", the design will be equivalent to the conditional logit model.
#' @param normalization the quantile to use for the normalization of the estimated coefficients
#'   where the logistic distribution is used as the base cumulative distribution function.
#' @param na.action an argument to handle missing data.
#'   Available options are na.omit, na.fail, and na.exclude.
#'   It comes from the stats library and does not include the na.pass option.
#' @param find_nu a logical argument to indicate whether the user intends to utilize the Student CDF and seeks an optimization algorithm to identify an optimal degrees of freedom setting for the model.
#' @param control a list specifying additional control parameters.
#'   - `maxit`: the maximum number of iterations for the Fisher scoring algorithm.
#'   - `epsilon`: a double value to fix the epsilon value.
#'   - `beta_init`: an appropriately sized vector for the initial iteration of the algorithm.
#'
#' @examples
#' library(GLMcat)
#' data(TravelChoice)
#'
#' discrete_cm(formula = choice ~ hinc + gc + invt,
#'             case_id = "indv", alternatives = "mode", reference = "air",
#'             data = TravelChoice,
#'             cdf = "logistic")
#'
#' #' Model with alternative specific effects for gc and invt:
#' discrete_cm(formula = choice ~ hinc + gc + invt,
#'             case_id = "indv", alternatives = "mode", reference = "air",
#'             data = TravelChoice, alternative_specific = c("gc", "invt"),
#'             cdf = "logistic")
#'
#' @note For these models, it is not allowed to exclude the intercept.
#' @export
discrete_cm <-
  function(
    formula,
    case_id,
    alternatives,
    reference,
    alternative_specific = NA,
    data ,
    cdf = list(),
    intercept = "standard",
    normalization = 1,
    control = list(),
    na.action = "na.omit",
    find_nu = FALSE){

    # check_categorical <- is.factor(model.frame(formula = formula, data)[, 1])
    # if (!check_categorical) {
    #   warning("The response variable is not defined as a categorical variable.")
    # }


    # check_ordered <- is.factor(model.frame(formula = formula, data)[,1])
    # if ( check_ordered == F ) { warning( "The response variable is not defined as a categorical variable" ) }

    # Set the na.action based on the specified option
    na.action <- match.arg(na.action, c("na.omit", "na.fail", "na.exclude"))
    if(na.action == "na.omit"){
      data <- na.omit(data)
    }else if (na.action == "na.fail"){
      data <- na.fail(data)
    }else if (na.action == "na.exclude"){
      data <- na.exclude(data)
    }

    if (normalization <= 0) {
      stop("Error: 'normalization' must be a positive integer.")
    }


    # Set the default value for cdf
    if(length(cdf)==0){cdf[[1]] = "logistic"}
    # Match the value of cdf to a valid option
    cdf[[1]] <- match.arg(cdf[[1]], c("logistic", "normal", "gumbel", "gompertz", "cauchy",
                                      "laplace", "student", "noncentralt"))

    # Match the value of intercept to a valid option
    intercept <- match.arg(intercept, c("standard","conditional"))
    control <- do.call(control_glmcat, c(control, list()))



    cdf_sel <- cdf

    # if(find_nu == TRUE) {
    #   # Estimate the models with Student link where ν = 1 and ν = 8
    #   cdf_1 <- list("student", 1)
    #   cdf_8 <- list("student", 8)
    #   model_1 <- .Discrete_CM(formula = formula,
    #                           case_id = case_id,
    #                           alternatives = alternatives,
    #                           reference = reference,
    #                           alternative_specific = alternative_specific,
    #                           data = data,
    #                           cdf = cdf_1,
    #                           intercept = intercept,
    #                           normalization = normalization,
    #                           control = control)
    #   model_8 <- .Discrete_CM(formula = formula,
    #                           case_id = case_id,
    #                           alternatives = alternatives,
    #                           reference = reference,
    #                           alternative_specific = alternative_specific,
    #                           data = data,
    #                           cdf = cdf_8,
    #                           intercept = intercept,
    #                           normalization = normalization,
    #                           control = control)
    #
    #   # Check if lν=8 > lν=1
    #   if (model_8[["LogLikelihood"]] > model_1[["LogLikelihood"]]) {
    #     # Estimate the log-likelihood lp of a binary model with the probit link
    #     model_p <- .Discrete_CM(formula = formula,
    #                             case_id = case_id,
    #                             alternatives = alternatives,
    #                             reference = reference,
    #                             alternative_specific = alternative_specific,
    #                             data = data,
    #                             cdf = "normal",
    #                             intercept = intercept,
    #                             normalization = normalization,
    #                             control = control)
    #     # Check if lp > lν=8
    #     if (model_p[["LogLikelihood"]] > model_8[["LogLikelihood"]]) {
    #       # Use the probit link
    #       cdf_sel <- "normal"
    #     } else {
    #       # Use the logit link
    #       cdf_sel <- "logistic"
    #     }
    #   } else {
    #     # Use optimize() to find the best ν ∈ (0.25, 1) of the Student CDF
    #     optimize_likelihood <- function(nu) {
    #       model_nu <- .Discrete_CM(formula = formula,
    #                                case_id = case_id,
    #                                alternatives = alternatives,
    #                                reference = reference,
    #                                alternative_specific = alternative_specific,
    #                                data = data,
    #                                cdf = list("student", nu),
    #                                intercept = intercept,
    #                                normalization = normalization,
    #                                control = control)
    #
    #       model_nu[["LogLikelihood"]]
    #     }
    #
    #     opt_result <- optimize(optimize_likelihood, interval = c(0.25, 8), maximum = TRUE)
    #     best_nu <- opt_result$maximum
    #
    #     # Use the best ν found
    #     cdf_sel <- list("student", best_nu)
    #   }
    #   # print(cdf_sel)
    # }

    #     # Use optimize() to find the best ν ∈ (0.25, 1) of the Student CDF
    optimize_likelihood <- function(nu) {
      model_nu <- .Discrete_CM(formula = formula,
                               case_id = case_id,
                               alternatives = alternatives,
                               reference = reference,
                               alternative_specific = alternative_specific,
                               data = data,
                               cdf = list("student", nu),
                               intercept = intercept,
                               normalization = normalization,
                               control = control)

      model_nu[["LogLikelihood"]]
    }

    if(find_nu == TRUE){
      opt_result <- optimize(optimize_likelihood, interval = c(0.25, 8), maximum = TRUE)
      best_nu <- opt_result$maximum
      # Use the best ν found
      cdf_sel <- list("student", best_nu)
    }

    cdf <- cdf_sel

    # if(is.null(alternative_specific)){alternative_specific <- NA}
    # Run the Discrete Choice Modeling C++ function
    fit_old <- .Discrete_CM(formula = formula,
                            case_id = case_id,
                            alternatives = alternatives,
                            reference = reference,
                            alternative_specific = alternative_specific,
                            data = data,
                            cdf = cdf,
                            intercept = intercept,
                            normalization = normalization,
                            control = control)

    formula1 <- paste(format(fit_old$formula),"+",fit_old$arguments$case_id,"+",
                      fit_old$arguments$alternatives,sep = "")

    fit_old$formula <- formula1
    fit_old[["model"]] <- model.frame(formula = gsub("\\[.*\\]","",as.character(formula1)),
                                      data)

    fit_old$table_summary <- table_summary(fit_old)

    fit_old <- fit_old[sort(names(fit_old))]
    class(fit_old) <- "glmcat"

    return(fit_old)

  }
