#' Variance-Covariance Matrix for a Fitted glmcat Model Object
#' @description Returns the variance-covariance matrix of the main parameters of a fitted \code{glmcat} model object.
#' @rdname glmcat2
#' @param formula a symbolic description of the model to be fit. An expression of the form y ~ predictors is interpreted as a specification that the response y is modelled by a linear predictor specified symbolically by model.
#' @param ratio a string indicating the F cdf, options are: reference, adjacent, cumulative and sequential. Default value is reference.
#' @param cdf
#' \describe{
  #' \item{\code{cdf}:}{a string indicating the F cdf, options are: logistic, normal, cauchy, student (any df), noncentralt, gompertz, gumbel and laplace.}
  #' \item{\code{df}:}{an integer with the degrees of freedom of the 'cdf'}
  #' \item{\code{mu}:}{an integer with the mu parameter of the 'cdf'}
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
#' @export
glmcat2 <-
  function(formula,
           data,
           ratio = c("reference", "cumulative", "sequential","adjacent"),
           cdf = c("logistic", "normal", "gumbel", "gompertz", "cauchy", "laplace"),
           parallel = NA,
           categories_order = NA,
           ref_category = NA,
           threshold = c("standard", "symmetric", "equidistant"),
           control = list(),
           normalization = 1,
           doFit = TRUE, na.action,
           contrasts,
           model = TRUE, ...)
  {
    cdf <- match.arg(cdf)
    ratio <- match.arg(ratio)
    threshold <- match.arg(threshold)
    contrasts <- NULL
    control <- do.call(control.glmcat, c(control, list(...)))

    fit_old <- GLMcat(formula = formula, data = data, ratio = ratio, cdf = cdf, parallel = parallel, categories_order = categories_order,
                       ref_category = ref_category, threshold = threshold , control = control, normalization = normalization)
    # return(list(formula = formula, data = data, ratio = ratio, cdf = cdf, parallel = parallel, categories_order = categories_order,
    #             ref_category = ref_category, threshold = threshold , control = control, normalization = normalization))

    fit_old[["model"]] <- model.frame(formula = formula,data)

    fit_old$table_summary <- table_summary(fit_old)

    fit_old <- fit_old[sort(names(fit_old))]
    class(fit_old) <- "glmcat"

    return(fit_old)
  }
