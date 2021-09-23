#' Fitting models for categorical responses
#' @description Estimate generalized linear models implemented under the unified
#' specification ( ratio,cdf,Z\) where \code{ratio} represents the ratio of probabilities
#' (reference, cumulative, adjacent, or sequential), \code{cdf} the cumulative distribution function
#' for the linkage, and Z the design matrix which must be specified through the \code{parallel}
#' and the \code{threshold} arguments.
#' @title glmcat
#' @rdname glmcat
#' @name glmcat
#' @param formula a symbolic description of the model to be fit. An expression of the form y ~ predictors is interpreted as a specification that the response y is modelled by a linear predictor specified symbolically by model.
#' @param ratio a string indicating the ratio (equivantly to the family) options are: reference, adjacent, cumulative and sequential. Default value is reference.
#' @param cdf
#' \describe{
#' The inverse distribution function to be used as part of the link function.
#' If the distribution has no parameters to specify then it should be entered as a
#' string indicating the name, e.g., \code{cdf = "normal"}, the default value is \code{cdf = "logistic"}.
#' If there are parameters to specify then a list must be entered,
#' so far this would only be the case for Student's distribution which would be
#' \code{list("student", df=2)},
#' and for the non-central distribution of student, \code{list("noncentralt", df=2, mu=1)},
#' }
#' @param categories_order a character vector indicating the incremental order of the categories: c("a", "b", "c"); a<b<c. Alphabetical order is assumed by default. Order is relevant for adjacent, cumulative and sequential ratio.
#' @param ref_category a string indicating the reference category. Proper option for models with reference ratio.
#' @param parallel a character vector indicating the name of the variables with a parallel effect. If variable is categorical, specify the name and the level of the variable as a string "namelevel".
#' @param data a dataframe object in R, with the dependent variable as factor.
#' @param threshold restriction to impose on the thresholds, options are: standard, equidistant or symmetric (Valid only for the cumulative ratio).
#' @param control
#' \describe{
#' \item{\code{maxit}:}{the maximum number of iterations for the Fisher scoring algorithm.}
#' \item{\code{epsilon}:}{a double to change update the convergence criterion of GLMcat models.}
#' \item{\code{beta_init}:}{an appropiate sized vector for the initial iteration of the algorithm.}
#' }
#' @param normalization the quantile to use for the normalization of the estimated coefficients where the logistic distribution is used as the base cumulative distribution function.
#' @param ... additional arguments.
#' @export
#' @references
#' \insertRef{peyhardi_new}
#' @examples
#' data(DisturbedDreams)
#' ref_log_com <- glmcat(formula = Level ~ Age, data = DisturbedDreams,
#'     ref_category = "Very.severe",
#'     cdf = "logistic", ratio = "reference")
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
           # doFit = TRUE, na.action,
           # contrasts, model = TRUE,
           ...)
  {
    if(length(cdf)==0){cdf[[1]] = "logistic"}
    cdf[[1]] <- match.arg(cdf[[1]], c("logistic", "normal", "gumbel", "gompertz", "cauchy",
                               "laplace", "student", "noncentralt"))
    ratio <- match.arg(ratio)
    threshold <- match.arg(threshold)
    contrasts <- NULL
    control <- do.call(control_glmcat, c(control, list(...)))

    # Default for reference ratio should be the complete design
    if(ratio == "reference" && is.na(parallel)){
      parallel = F
    }

    fit_old <- .GLMcat(formula = formula, data = data, ratio = ratio, cdf = cdf, parallel = parallel, categories_order = categories_order,
                      ref_category = ref_category, threshold = threshold , control = control, normalization = normalization)

    fit_old[["model"]] <- model.frame(formula = formula,data)

    fit_old$table_summary <- table_summary(fit_old)

    fit_old <- fit_old[sort(names(fit_old))]
    class(fit_old) <- "glmcat"

    return(fit_old)
  }

#' Family of models for Discrete Choice
#' @description Discrete choice model: Requires data in long form.
#' For each individual (or decision maker), there are multiple observations (rows),
#' one for each of the alternatives the individual could have chosen.
#' We call the group of observations for an individual a “case”.
#' Each case represents a single statistical observation although it comprises
#' multiple observations.
#' @title Discrete_CM
#' @rdname discrete_cm
#' @name discrete_cm
#' @param formula a symbolic description of the model to be fit. An expression of the form y ~ predictors is interpreted as a specification that the response y is modelled by a linear predictor specified symbolically by model. A particularity for the formula is that for the case-specific variables, the user can define a specific effect for a category.
#' @param case_id a string with the name of the column that identifies each case.
#' @param alternatives a string with the name of the column that identifies the vector of alternatives the individual could have chosen.
#' @param reference a string indicating the reference category
#' @param alternative_specific a character vector with the name of the explanatory variables that are different for each case, these are the alternative specific variables. By default, the case specific variables are the explanatory variables that are not identify in here, but that are part of the formula.
#' @param data a dataframe (in a long format) object in R, with the dependent variable as factor.
#' @param cdf
#' \describe{
#' The inverse distribution function to be used as part of the link function.
#' If the distribution has no parameters to specify then it should be entered as a
#' string indicating the name, e.g., \code{cdf = "normal"}, the default value is \code{cdf = "logistic"}.
#' If there are parameters to specify then a list must be entered,
#' so far this would only be the case for Student's distribution which would be
#' \code{list("student", df=2)},
#' and for the non-central distribution of student, \code{list("noncentralt", df=2, mu=1)},
#' }
#' @param intercept if "conditional" then the design will be equivalent to the conditional logit model
#' @param normalization the quantile to use for the normalization of the estimated coefficients where the logistic distribution is used as the base cumulative distribution function.
#' @param control
#' \describe{
#' \item{\code{maxit}:}{the maximum number of iterations for the Fisher scoring algorithm.}
#' \item{\code{epsilon}:}{a double with to fix the epsilon value}
#' \item{\code{beta_init}:}{an appropiate sized vector for the initial iteration of the algorithm}}
#' @examples
#' library(GLMcat)
#' data(TravelChoice)
#' discrete_cm(formula = choice ~ hinc + gc + invt,
#' case_id = "indv",alternatives = "mode", reference = "air",
#' data = TravelChoice,  alternative_specific = c("gc", "invt"),
#' cdf = "logistic")
#' @note For these models it is not allowed to exclude the intercept.
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
    control = list() ){

    if(length(cdf)==0){cdf[[1]] = "logistic"}
    cdf[[1]] <- match.arg(cdf[[1]], c("logistic", "normal", "gumbel", "gompertz", "cauchy",
                                      "laplace", "student", "noncentralt"))

    intercept <- match.arg(intercept, c("standard","conditional"))
    control <- do.call(control_glmcat, c(control, list()))


    # if(is.null(alternative_specific)){alternative_specific <- NA}

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
