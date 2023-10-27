#' Print method for a fitted \code{glmcat} model object
#' @description \code{print} method for a fitted \code{glmcat} model object.
#' @param x an object of class \code{glmcat}.
#' @param ... additional arguments.
#' @rdname print
#' @examples
#' model <- glmcat(formula = Level ~ Age, data = DisturbedDreams,
#'                 ref_category = "Very.severe", ratio = "cumulative")
#' print(model)
#' @exportS3Method
print.glmcat <- function(x, ...) {
  cat("\nFormula:\n")
  print(x$formula)
  print(x$table_summary)
  cat("\nCoefficients:\n")
  print(coef(x, with_baseline = FALSE))
  ll <- logLik(x)
  cat("\nLog-Likelihood:\n ", ll, " (df = ", attr(ll, "df"), ")", sep = "")
  cat("\n\n")
  invisible(x)
}

#' Plot method for a fitted \code{glmcat} model object
#' @description \code{plot} of the log-likelihood profile for a fitted \code{glmcat} model object.
#' @param x an object of class \code{glmcat}.
#' @param ... additional arguments.
#' @rdname plot
#' @exportS3Method
plot.glmcat <- function(x, ...) {
  log_iter <- x$LogLikIter
  plot(log_iter[-1],main = "Log-likelihood profile", xlab = "Iteration", ylab = "Log-likelihood")
  lines(log_iter[-1])
}


#' Variance-Covariance Matrix for a fitted \code{glmcat} model object
#' @description Returns the variance-covariance matrix of the main parameters of a fitted \code{glmcat} model object.
#' @param object an object of class \code{glmcat}.
#' @param ... additional arguments.
#' @rdname vcov
#' @method vcov glmcat
#' @usage \method{vcov}{glmcat}(object,...)
#' @exportS3Method
vcov.glmcat <- function(object,...) {
  colnames(object$cov_beta) <- rownames(object$cov_beta) <- rownames(object$coefficients)
  return(object$cov_beta)
}

#' Terms of a fitted \code{glmcat} model object
#' @description Returns the terms of a fitted \code{glmcat} model object.
#' @param x an object of class \code{glmcat}.
#' @param ... additional arguments.
#' @rdname terms
#' @method terms glmcat
#' @usage \method{terms}{glmcat}(x, ...)
#' @exportS3Method
terms.glmcat <- function(x,...) {
  return(terms(x$formula))
}

#' Predict method for a a fitted \code{glmcat} model object
#' @description Obtains predictions of a fitted \code{glmcat} model object.
#' @param object a fitted object of class \code{glmcat}.
#' @param newdata optionally, a data frame in which to look for the variables involved in the model. If omitted, the fitted linear predictors are used.
# #' @param se.fit should standard errors of the predictions be provided? Not applicable and ignored when \code{type = "class"}.
# #' @param interval should confidence intervals for the predictions be provided?  Not applicable and ignored when \code{type = "class"}.
# #' @param level the confidence level.
#' @param type the type of prediction required.
#' The default is \code{"prob"} which gives the probabilities, the other option is
#' \code{"linear.predictor"} which gives predictions on the scale of the linear predictor.
# #' @param na.action function determining what should be done with missing values in \code{newdata}. The default is to predict \code{NA}.
#' @param ... further arguments.
#' The default is \code{"prob"} which gives the probabilities, the other option is
#' \code{"linear.predictor"} which gives predictions on the scale of the linear predictor.
#' @rdname predict
#' @method predict glmcat
#' @usage \method{predict}{glmcat}(object, newdata, type, ...)
#' @exportS3Method
predict.glmcat <- function(object,
                           newdata,
                           # se.fit = FALSE,
                           # interval = FALSE,
                           # level = 0.95,
                           type = c("prob", "linear.predictor"),
                           # na.action = na.pass,
                           ...){
  if (missing(type)) { type <- "prob" }
  type <- match.arg(type, c("prob", "linear.predictor"))
  if (missing(newdata)) {
    # if (object$Function == "DiscreteCM"){
    #   object1 <- object
    #   formula1 <- paste(format(object1$formula),"+",object1$arguments$case_id,"+",
    #                     object1$arguments$alternatives,sep = "")
    #   object1$formula <- formula1
    #   newdata <- model.frame(object1)
    #   print(object1$formula)
    #
    # }else{
    newdata <- model.frame(object)
    # }
  }

  newdata[with(attributes(terms(object)), as.character(variables[response+1]))] <- object$categories_order[1]

  return(.predict_glmcat(model_object = object, data = newdata, type = type))
}

#' Confidence intervals for parameters of a fitted \code{glmcat} model object
#' @description Computes confidence intervals from a fitted \code{glmcat} model object for all the parameters.
#' @param object an fitted object of class \code{glmcat}.
#' @param parm a numeric or character vector indicating which regression coefficients should be displayed
#' @param level the confidence level.
#' @param ... other parameters.
#' @rdname confint.glmcat
#' @method confint glmcat
#' @usage \method{confint}{glmcat}(object, parm, level, ...)
#' @exportS3Method
confint.glmcat <-
  function(object, parm = NULL, level = 0.95, ...)
  {
    stopifnot(is.numeric(level) && length(level) == 1 && level > 0 && level < 1)
    lev <- (1 - level)/2
    lev <- c(lev, 1 - lev)
    pct <- paste(format(100 * lev, trim = TRUE, scientific = FALSE,  digits = 3), "%")
    fac <- qnorm(lev)
    coefs <- coef(object)
    parnames <- rownames(coefs)
    if(is.character(parm))
      parm <- match(parm, parnames, nomatch = 0)
    if(is.null(parm)){
      parm = seq_len(length(coefs))
    }
    if(!all(parm %in% seq_len(length(coefs))))
      stop("invalid 'parm' argument")
    stopifnot(length(parm) > 0)

    ses <- coef(summary(object))[, 2]
    ci <- array(NA, dim = c(length(coefs), 2L), dimnames = list(names(coefs), pct))
    ci[] <- cbind(coefs,coefs) + ses %o% fac
    rownames(ci) <- rownames(coefs)

    ci <- ci[parm,]
    return(ci)
  }


#' Summary method for a fitted \code{glmcat} model object
#' @description Summary method for a fitted `glmcat` model object.
#' @param object an fitted object of class `glmcat`.
#' @param normalized if `TRUE`, the summary method yields the normalized coefficients.
#' @param correlation if `TRUE`, prints the correlation matrix.
#' @param ... additional arguments affecting the summary produced.
#' @rdname summary
#' @method summary glmcat
#' @exportS3Method
#' @examples
#' mod1 <- discrete_cm(formula = choice ~ hinc + gc + invt,
#'                     case_id = "indv", alternatives = "mode", reference = "air",
#'                     data = TravelChoice,  alternative_specific = c("gc", "invt"),
#'                     cdf = "normal", normalization = 0.8)
#' summary(mod1, normalized = TRUE)
summary.glmcat <- function(object, normalized = FALSE, correlation = FALSE,...) {
  vcov <- object$cov_beta
  coefs <- matrix(NA, length(object$coefficients), 4,
                  dimnames = list(names(object$coefficients),
                                  c("Estimate", "Std. Error", "z value", "Pr(>|z|)")))

  coefs[, 1] <- object$coefficients
  coefs[, 2] <- sd <- sqrt(diag(vcov))
  # Check if normalized coefficients are requested
  if(normalized){
    cat("Normalized coefficients with s0 = ",object$normalization_s0, "\n")
    coefs <- matrix(NA, length(object$coefficients*object$normalization_s0), 4,
                    dimnames = list(names(object$coefficients*object$normalization_s0),
                                    c("Estimate", "Std. Error", "z value", "Pr(>|z|)")))
    coefs[, 1] <- object$coefficients*object$normalization_s0
    coefs[, 2] <- sd <- sqrt(diag(vcov))*object$normalization_s0
  }

  if(!all(is.finite(vcov))) {
    ## warning("Variance-covariance matrix of the parameters is not defined")
    coefs[, 2:4] <- NA
    if(correlation) warning("Correlation matrix is unavailable")
  }
  else {
    # alias <- unlist(object$aliased)

    ## Cond is Inf if Hessian contains NaNs:
    object$cond.H <-
      if(any(is.na(object$Hessian))) Inf
    else with(eigen(object$Hessian, symmetric=TRUE, only.values = TRUE),
              abs(max(values) / min(values)))
    coefs[, 3] <- coefs[, 1]/coefs[, 2]
    coefs[, 4] <- 2 * pnorm(abs(coefs[, 3]),
                            lower.tail=FALSE)
    if(correlation)
      object$correlation <- cov2cor(vcov)
  }

  # coefs[, 1] <- object$coefficients

  rownames(coefs) <- rownames(object$coefficients)

  object$coefficients <- coefs
  class(object) <- "summary.glmcat"
  object
}

#' Model coefficients of a fitted \code{glmcat} model object
#' @description Returns the coefficient estimates of the fitted \code{glmcat} model object.
#' @param object an fitted object of class \code{glmcat}.
#' @param na.rm TRUE for NA coefficients to be removed, default is FALSE.
#' @rdname coef
#' @param ... additional arguments affecting the \code{coef} method.
#' @exportS3Method
coef.glmcat <- function(object, na.rm = FALSE, ...) {
  if (na.rm) {
    coefs <- object$coefficients
    coefs[!is.na(coefs)]
  }
  else {
    object$coefficients
  }
}

#' Number of observations of a fitted \code{glmcat} model object
#' @description Extract the number of observations of the fitted \code{glmcat} model object.
#' @param object an fitted object of class \code{glmcat}.
#' @param ... additional arguments affecting the \code{nobs} method.
#' @rdname nobs
#' @method nobs glmcat
#' @exportS3Method
nobs.glmcat <- function(object,...) {
  return(object$nobs_glmcat)
}

#' Log-likelihood of a fitted \code{glmcat} model object
#' @description Extract Log-likelihood of a fitted \code{glmcat} model object.
#' @rdname logLik
#' @param object an fitted object of class \code{glmcat}.
#' @param ... additional arguments affecting the loglik.
#' @method logLik glmcat
#' @exportS3Method
logLik.glmcat <- function(object,...) {
  structure(object$LogLikelihood,
            df = object$df, nobs_glmcat = object$nobs_glmcat,
            class = "logLik"
  )
}

#' Extract AIC from a fitted \code{glmcat} model object
#' @description Method to compute the (generalized) Akaike An Information Criterion for a fitted object of class \code{glmcat}.
#' @rdname extractAIC
#' @param fit an fitted object of class \code{glmcat}.
#' @param ... further arguments (currently unused in base R).
#' @method extractAIC glmcat
#' @examples
#' model <- glmcat(formula = Level ~ Age, data = DisturbedDreams,
#'                 ref_category = "Very.severe", ratio = "cumulative")
#' extractAIC(model)
#' @exportS3Method
extractAIC.glmcat <- function(fit, ...) {
  scale = 0
  k = 2
  edf <- fit$df
  c(edf, -2*fit$LogLikelihood + k * edf)
}

#' Control parameters for \code{glmcat} models
#' @description Set control parameters for \code{glmcat} models.
#' @rdname control_glmcat
#' @param maxit the maximum number of the Fisher's Scoring Algorithm iterations. Defaults to 25.
#' @param epsilon a double to change update the convergence criterion of GLMcat models.
#' @param beta_init an appropriate sized vector for the initial iteration of the algorithm.
#' @export
control_glmcat <- function(maxit = 25, epsilon = 1e-06, beta_init = NA) {
  return(list("maxit" = maxit, "epsilon" = epsilon, "beta_init" = beta_init))
  # return(maxit)
}

# Computes the lower tail probability of the chi-square distribution
# with non-negative degrees of freedom, replacing non-positive values
# of degrees of freedom with NA.
# Parameters:
# - q: quantiles
# - df: degrees of freedom
# - ...: additional arguments passed to the pchisq function
safe_pchisq <- function(q, df, ...) {
  df[df <= 0] <- NA
  pchisq(q = q, df = df, ...)
}

# Performs single term deletions from a model object and returns the resulting analysis of deviance table.
# Parameters:
# - object: a model object
# - scope: a character vector specifying the terms to be dropped from the model
# - data: a data frame or matrix containing the data
# - scale: a numeric value specifying the scale parameter
# - test: a character string specifying the test type ("none" or "Chisq")
# - trace: a logical value indicating whether to display trace information
# - ...: additional arguments
drop2 <- function(object, scope, data, scale = 0, test=c("none", "Chisq"),
                  trace = FALSE,  ...)
{

  tl <- attr(terms(object), "term.labels")
  if(missing(scope)) scope <- drop.scope(object)
  else {
    if(!is.character(scope))
      scope <- attr(terms(update.formula(object$formula, scope)), "term.labels")
    if(!all(match(scope, tl, 0L) > 0L))
      stop("scope is not a subset of term labels")
  }
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1L, ncol = 2L,
                dimnames =  list(c("<none>", scope), c("df", "AIC")))
  ans[1, ] <- AIC(object)
  n0 <- nobs(object, use.fallback = TRUE)
  env <- environment(formula(object))
  for(i in seq_len(ns)) {
    tt <- scope[i]
    if(trace > 1) {
      cat("trying -", tt, "\n", sep = "")
      flush.console()
    }

    nfit <- update(object$formula, as.formula(paste("~ . -", tt)),
                   evaluate = FALSE)

    fun_in = object$Function
    if(fun_in == "GLMcat"){
      object$parallel <- object$parallel[object$parallel != tt]
      if(length(object$parallel) == 0) {object$parallel <- NA}
      nfit <- .GLMcat(nfit, data, object$ratio, object$cdf, object$parallel,
                     object$categories_order, object$ref_category,
                     object$threshold, control_glmcat(object$control$maxit, object$control$epsilon, object$control$beta_init),
                     object$normalization_s0)
    }else{
      object$arguments$alternative_specific <- object$arguments$alternative_specific[object$arguments$alternative_specific != tt]
      if(length(object$arguments$alternative_specific) == 0) {object$arguments$alternative_specific <- NA}
      nfit <- .Discrete_CM(formula = nfit,
                          data = data,
                          cdf = object$cdf,
                          case_id = object$arguments$case_id,
                          alternatives = object$arguments$alternatives,
                          alternative_specific = object$arguments$alternative_specific,
                          # object$categories_order,
                          intercept = object$arguments$intercept,
                          reference = object$arguments$reference,
                          control_glmcat(object$control$maxit,
                                         object$control$epsilon, object$control$beta_init),
                          normalization = object$normalization_s0)

    }
    ans[i+1, ] <- AIC(nfit)
    nnew <- nobs(nfit, use.fallback = TRUE)
    if(all(is.finite(c(n0, nnew))) && nnew != n0)
      stop("number of rows in use has changed: remove missing values?")
  }
  dfs <- ans[1L , 1L] - ans[, 1L]
  dfs[1L] <- NA
  aod <- data.frame(Df = dfs, AIC = ans[,2])
  test <- match.arg(test)
  if(test == "Chisq") {
    dev <- ans[, 2L] - 2*ans[, 1L]
    dev <- dev - dev[1L] ; dev[1L] <- NA
    nas <- !is.na(dev)
    P <- dev
    P[nas] <- safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
    aod[, c("LRT", "Pr(>Chi)")] <- list(dev, P)
  }
  head <- c("Single term deletions", "\nModel:", deparse(formula(object)),
            if(scale > 0) paste("\nscale: ", format(scale), "\n"))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}

# Performs single term additions to a model object and returns the resulting analysis of deviance table.
# Parameters:
# - object: a model object
# - scope: a character vector specifying the terms to be added to the model
add2 <- function(object, scope, data, scale = 0, test=c("none", "Chisq"),
                 trace = FALSE, ...)
{
  if(missing(scope) || is.null(scope)) stop("no terms in scope")
  if(!is.character(scope))
    scope <- add.scope(object, update.formula(object, scope))
  if(!length(scope))
    stop("no terms in scope for adding to object")
  #     newform <- update.formula(object,
  #                               paste(". ~ . +", paste(scope, collapse="+")))
  #     data <- model.frame(update(object, newform)) # remove NAs
  #     object <- update(object, data = data)
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1L, ncol = 2L,
                dimnames = list(c("<none>", scope), c("df", "AIC")))
  ans[1L,  ] <- AIC(object)
  n0 <- nobs(object, use.fallback = TRUE)
  # env <- environment(formula(object))
  for(i in seq_len(ns)) {
    tt <- scope[i]
    if(trace > 1) {
      cat("trying +", tt, "\n", sep = "")
      flush.console()
    }
    # nfit <- update(object, as.formula(paste("~ . +", tt)),
    #                evaluate = FALSE)
    # nfit <- eval(nfit, envir=env) # was  eval.parent(nfit)

    nfit <- update(object$formula, as.formula(paste("~ . +", tt)),
                   evaluate = FALSE)

    # object$parallel <- object$parallel[object$parallel != tt]
    fun_in = object$Function
    if(fun_in == "GLMcat"){
      # object$parallel <- object$parallel[object$parallel != tt]
      # if(length(object$parallel) == 0) {object$parallel <- NA}
      nfit <- .GLMcat(nfit, data, object$ratio, object$cdf, object$parallel,
                     object$categories_order, object$ref_category,
                     object$threshold, control_glmcat(object$control$maxit, object$control$epsilon, object$control$beta_init),
                     object$normalization_s0)
    }else{
      # object$arguments$alternative_specific <- object$arguments$alternative_specific[object$arguments$alternative_specific != tt]
      # if(length(object$arguments$alternative_specific) == 0) {object$arguments$alternative_specific <- NA}
      nfit <- .Discrete_CM(formula = nfit,
                          data = data,
                          cdf = object$cdf,
                          alternative_specific = object$arguments$alternative_specific,
                          # object$categories_order,
                          case_id = object$arguments$case_id,
                          alternatives = object$arguments$alternatives,
                          reference = object$arguments$reference,
                          intercept = object$arguments$intercept,
                          control_glmcat(object$control$maxit,
                                         object$control$epsilon, object$control$beta_init),
                          normalization = object$normalization_s0)
    }

    ans[i+1L, ] <- AIC(nfit)
    nnew <- nobs(nfit, use.fallback = TRUE)
    if(all(is.finite(c(n0, nnew))) && nnew != n0)
      stop("number of rows in use has changed: remove missing values?")
  }
  dfs <- ans[, 1L] - ans[1L, 1L]
  dfs[1L] <- NA
  aod <- data.frame(Df = dfs, AIC = ans[, 2L])
  test <- match.arg(test)
  if(test == "Chisq") {
    dev <- ans[, 2L] - 2*ans[, 1L]
    dev <- dev[1L] - dev; dev[1L] <- NA
    nas <- !is.na(dev)
    P <- dev
    P[nas] <- safe_pchisq(dev[nas], dfs[nas], lower.tail=FALSE)
    aod[, c("LRT", "Pr(>Chi)")] <- list(dev, P)
  }
  head <- c("Single term additions", "\nModel:", deparse(formula(object)),
            if(scale > 0) paste("\nscale: ", format(scale), "\n"))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}

#' Stepwise for a \code{glmcat} model object
#' @description Stepwise for a \code{glmcat} model object based on the AIC.
#' @param object an fitted object of class \code{glmcat}.
#' @param scope defines the range of models examined in the stepwise search (same as in the step function of the stats package). This should be either a single formula, or a list containing components upper and lower, both formulae.
#' @param direction the mode of the stepwise search.
#' @param trace to print the process information.
#' @param steps the maximum number of steps.
#' @rdname step
#' @method step glmcat
#' @usage \method{step}{glmcat}(object, scope, direction, trace, steps)
#' @exportS3Method
step.glmcat <- function (object,
                         scope,
                         direction = c("both", "backward", "forward"),
                         trace = 1, steps = 1000)
{
  data <- object$data
  mydeviance <- function(x, ...) {
    dev <- deviance(x)
    if (!is.null(dev))
      dev
    else AIC(x)
  }
  cut.string <- function(string) {
    if (length(string) > 1L)
      string[-1L] <- paste0("\n", string[-1L])
    string
  }
  # re.arrange <- function(keep) {
  #   namr <- names(k1 <- keep[[1L]])
  #   namc <- names(keep)
  #   nc <- length(keep)
  #   nr <- length(k1)
  #   array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr,
  #                                                          namc))
  # }
  step.results <- function(models, fit, object, usingCp = FALSE) {
    change <- sapply(models, "[[", "change")
    rd <- sapply(models, "[[", "deviance")
    dd <- c(NA, abs(diff(rd)))
    rdf <- sapply(models, "[[", "df.resid")
    ddf <- c(NA, diff(rdf))
    AIC <- sapply(models, "[[", "AIC")
    heading <- c("Stepwise Model Path \nAnalysis of Deviance Table",
                 "\nInitial Model:", deparse(formula(object)), "\nFinal Model:",
                 deparse(formula(fit)), "\n")
    aod <- data.frame(Step = I(change), Df = ddf, Deviance = dd,
                      `Resid. Df` = rdf, `Resid. Dev` = rd, AIC = AIC,
                      check.names = FALSE)
    if (usingCp) {
      cn <- colnames(aod)
      cn[cn == "AIC"] <- "Cp"
      colnames(aod) <- cn
    }
    attr(aod, "heading") <- heading
    fit$anova <- aod
    fit
  }

  object$terms <-   terms(formula(object$formula), data = data)

  Terms <- terms(object)
  object$call$formula <- object$formula <- Terms
  md <- missing(direction)
  direction <- match.arg(direction)
  backward <- direction == "both" | direction == "backward"
  forward <- direction == "both" | direction == "forward"
  if (missing(scope)) {
    fdrop <- numeric()
    fadd <- attr(Terms, "factors")
    if (md)
      forward <- FALSE
  }
  else {
    if (is.list(scope)) {
      fdrop <- if (!is.null(fdrop <- scope$lower))
        attr(terms(update.formula(object, fdrop)), "factors")
      else numeric()
      fadd <- if (!is.null(fadd <- scope$upper))
        attr(terms(update.formula(object, fadd)), "factors")
    }
    else {
      fadd <- if (!is.null(fadd <- scope))
        attr(terms(update.formula(object, scope)), "factors")
      fdrop <- numeric()
    }
  }
  models <- vector("list", steps)
  # if (!is.null(keep))
  #   keep.list <- vector("list", steps)
  n <- nobs(object, use.fallback = TRUE)
  fit <- object
  bAIC <- AIC(fit)
  edf <- length(attr(fit$terms,"term.labels"))
  # bAIC <- bAIC[2L]
  if (is.na(bAIC))
    stop("AIC is not defined for this model, so 'step' cannot proceed")
  if (bAIC == -Inf)
    stop("AIC is -infinity for this model, so 'step' cannot proceed")
  nm <- 1
  if (trace) {
    cat("Start:  AIC=", format(round(bAIC, 2)), "\n", cut.string(deparse(formula(fit))),
        "\n\n", sep = "")
    flush.console()
  }
  models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n -
                         edf, change = "", AIC = bAIC)
  # if (!is.null(keep))
  #   keep.list[[nm]] <- keep(fit, bAIC)
  usingCp <- FALSE
  while (steps > 0) {
    steps <- steps - 1
    AIC <- bAIC
    ffac <- attr(Terms, "factors")
    scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
    aod <- NULL
    change <- NULL
    if (backward && length(scope$drop)) {
      aod <- drop2(fit, scope$drop, data, scale = 0, trace = trace)
      rn <- row.names(aod)
      # print(aod)
      row.names(aod) <- c(rn[1L], paste("-", rn[-1L]))
      if (any(aod$Df == 0, na.rm = TRUE)) {
        zdf <- aod$Df == 0 & !is.na(aod$Df)
        change <- rev(rownames(aod)[zdf])[1L]
      }
    }
    if (is.null(change)) {
      if (forward && length(scope$add)) {
        aodf <- add2(fit, scope$add, data, scale = 0,
                     trace = trace)
        rn <- row.names(aodf)
        row.names(aodf) <- c(rn[1L], paste("+", rn[-1L]))
        aod <- if (is.null(aod))
          aodf
        else rbind(aod, aodf[-1, , drop = FALSE])
      }
      attr(aod, "heading") <- NULL
      nzdf <- if (!is.null(aod$Df))
        aod$Df != 0 | is.na(aod$Df)
      aod <- aod[nzdf, ]
      if (is.null(aod) || ncol(aod) == 0)
        break
      nc <- match(c("Cp", "AIC"), names(aod))
      nc <- nc[!is.na(nc)][1L]
      o <- order(aod[, nc])
      if (trace)
        print(aod[o, ])
      if (o[1L] == 1)
        break
      change <- rownames(aod)[o[1L]]
    }
    usingCp <- match("Cp", names(aod), 0L) > 0L
    # fit <- update(fit, paste("~ .", change), evaluate = FALSE)
    # fit <- eval.parent(fit)

    form1 <- update(fit$formula, as.formula(paste("~ .", change)),
                    evaluate = FALSE)
    object$parallel <- object$parallel[object$parallel != str_trim(sub("-","", change))]
    if(length(object$parallel) == 0) {object$parallel <- NA}

    # fit <- GLMcat(form1, data, object$ratio, object$cdf, object$parallel,
    #               object$categories_order, object$ref_category,
    #               object$threshold, control_glmcat(object$control$maxit, object$control$epsilon, object$control$beta_init),
    #               object$normalization_s0)

    fun_in = object$Function
    if(fun_in == "GLMcat"){
      fit <- .GLMcat(form1, data, object$ratio, object$cdf, object$parallel,
                    object$categories_order, object$ref_category,
                    object$threshold, control_glmcat(object$control$maxit, object$control$epsilon, object$control$beta_init),
                    object$normalization_s0)
    }else{
      # object$arguments$alternative_specific <- object$arguments$alternative_specific[object$arguments$alternative_specific != tt]
      # if(length(object$arguments$alternative_specific) == 0) {object$arguments$alternative_specific <- NA}
      fit <- .Discrete_CM(formula = form1,
                         data = data,
                         cdf = object$cdf,
                         alternative_specific = object$arguments$alternative_specific,
                         # object$categories_order,
                         case_id = object$arguments$case_id,
                         alternatives = object$arguments$alternatives,
                         reference = object$arguments$reference,
                         intercept = object$arguments$intercept,
                         control_glmcat(object$control$maxit,
                                        object$control$epsilon, object$control$beta_init),
                         normalization = object$normalization_s0)
    }

    fit$terms <-   terms(formula(fit$formula), data = fit$data)

    nnew <- nobs(fit, use.fallback = TRUE)
    if (all(is.finite(c(n, nnew))) && nnew != n)
      stop("number of rows in use has changed: remove missing values?")
    Terms <- terms(fit)
    bAIC <- AIC(fit)
    edf <- length(attr(fit$terms,"term.labels"))
    if (trace) {
      cat("\nStep:  AIC=", format(round(bAIC, 2)), "\n",
          cut.string(deparse(formula(fit))), "\n\n", sep = "")
      flush.console()
    }
    if (bAIC >= AIC + 1e-07)
      break
    nm <- nm + 1
    models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n -
                           edf, change = change, AIC = bAIC)
    # if (!is.null(keep))
    #   keep.list[[nm]] <- keep(fit, bAIC)
  }
  # if (!is.null(keep))
  #   fit$keep <- re.arrange(keep.list[seq(nm)])
  step.results(models = models[seq(nm)], fit, object, usingCp)
}

