#' Summary of models
#' @description \code{summary} method for GLMcat objects.
#' @param object a GLMcat model
#' @param ... additional arguments affecting the summary produced.
#' @rdname summary
#' @export
summary.glmcat <- function(object, ...) {
  coef <- object$coefficients
  se <- object$stderr
  s0 <- object$normalization_s0

  tval <- coef / se

  object$coefficients <- cbind(
    "Estimate" = coef,
    "Std. Error" = se,
    "z value" = tval,
    "Pr(>|z|)" = 2 * pnorm(-abs(tval))
  )

  colnames(object$coefficients) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  sum_ma <- object$coefficients
  printCoefmat(object$coefficients, P.values = TRUE, has.Pvalue = TRUE, ...)

  if(s0 !=1 ){
    print("Normalized coefficients")
    object$coefficients <- cbind(
      "Estimate" = coef * s0,
      "Std. Error" = se * s0,
      "z value" = tval,
      "Pr(>|z|)" = 2 * pnorm(-abs(tval))
    )
    colnames(object$coefficients) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    printCoefmat(object$coefficients, P.values = TRUE, has.Pvalue = TRUE, ...)
  }

}

#' Model coefficients
#' @description Extract model coefficients from a glmcat object.
#' @param object a GLMcat model.
#' @param na.rm TRUE for NA coefficients to be removed, default is FALSE.
#' @param ...	other arguments.
#' @rdname coef
#' @export
#' @examples
#' data(DisturbedDreams)
#' mod1 <- GLMcat(
#'   formula = Level ~ Age,
#'   ref_category = "Very.severe",
#'   data = DisturbedDreams, cdf = "logistic"
#' )
#' coef(mod1)
coef.glmcat <- function(object, na.rm = FALSE, ...) {
  if (na.rm) {
    coefs <- object$coefficients
    coefs[!is.na(coefs)]
  }
  else {
    object$coefficients
  }
}

#' Number of observations in a glmcat model
#' @description Extract the number of observations from a GLMcat model.
#' @param object a GLMcat model.
#' @param ...	other arguments.
#' @rdname nobs_glmcat
#' @export
#' @examples
#' data(DisturbedDreams)
#' mod1 <- GLMcat(
#'   formula = Level ~ Age,
#'   categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
#'   data = DisturbedDreams, cdf = "logistic"
#' )
#' nobs_glmcat(mod1)
nobs_glmcat <- function(object, ...) {
  return(object$nobs_glmcat)
}

#' LogLikelihood glmcat models
#' @description Extract LogLikelihood for GLMcat models.
#' @rdname logLik
#' @param object a GLMcat model.
#' @param ...	other arguments.
#' @export
#' @examples
#' data(DisturbedDreams)
#' mod1 <- GLMcat(
#'   formula = Level ~ Age,
#'   categories_order = c("Not.severe", "Severe.1", "Severe.2", "Very.severe"),
#'   data = DisturbedDreams, cdf = "logistic"
#' )
#' logLik(mod1)
logLik.glmcat <- function(object, ...) {
  structure(object$LogLikelihood,
            df = object$df, nobs_glmcat = object$nobs_glmcat,
            class = "logLik"
  )
}


#' control glmcat models
#' @description Set control parameters for GLMcat models.
#' @rdname control
#' @param maxit the maximum number of the Fisher's Scorng Algorithm iterations. Defaults to 25.
#' @param epsilon a double to change update the convergence criterion of GLMcat models.
#' @param beta_init an appropiate sized vector for the initial iteration of the algorithm.
#' @export
control.glmcat <- function(maxit = 25, epsilon = 1e-06, beta_init = NA) {
  return(list("maxit" = maxit, "epsilon" = epsilon, "beta_init" = beta_init))
  # return(maxit)
}

#' student glmcat models
#' @description Student distribution
#' @rdname student
#' @param df degrees_freedom
#' @export
student.glmcat <- function(df = 7) {
  return(list("cdf" = "student", "df" = df))
}

#' Noncentral t cdf for glmcat models
#' @description Noncentral t cdf
#' @rdname noncentralt
#' @param df degrees_freedom
#' @param mu non centrality parameter
#' @export
noncentralt.glmcat <- function(df = 7, mu = 0) {
  return(list("cdf" = "noncentralt", "df" = df, "mu" = mu))
}

safe_pchisq <- function(q, df, ...) {
  df[df <= 0] <- NA
  pchisq(q = q, df = df, ...)
}

drop2 <- function(object, scope, scale = 0, test=c("none", "Chisq"),
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
  n0 <- nobs_glmcat(object, use.fallback = TRUE)
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
      nfit <- GLMcat(nfit, object$data, object$ratio, object$cdf, object$parallel,
                     object$categories_order, object$ref_category,
                     object$threshold, control.glmcat(object$control$maxit, object$control$epsilon, object$control$beta_init),
                     object$normalization_s0)
    }else{
      object$arguments$alternative_specific <- object$arguments$alternative_specific[object$arguments$alternative_specific != tt]
      if(length(object$arguments$alternative_specific) == 0) {object$arguments$alternative_specific <- NA}
      nfit <- Discrete_CM(formula = nfit,
                          data = object$data,
                          cdf = object$cdf,
                          case_id = object$arguments$case_id,
                          alternatives = object$arguments$alternatives,
                          alternative_specific = object$arguments$alternative_specific,
                          # object$categories_order,
                          intercept = object$arguments$intercept,
                          reference = object$arguments$reference,
                          control.glmcat(object$control$maxit,
                                         object$control$epsilon, object$control$beta_init),
                          normalization = object$normalization_s0)

    }
    ans[i+1, ] <- AIC(nfit)
    nnew <- nobs_glmcat(nfit, use.fallback = TRUE)
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

add2 <- function(object, scope, scale = 0, test=c("none", "Chisq"),
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
  n0 <- nobs_glmcat(object, use.fallback = TRUE)
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
      nfit <- GLMcat(nfit, object$data, object$ratio, object$cdf, object$parallel,
                     object$categories_order, object$ref_category,
                     object$threshold, control.glmcat(object$control$maxit, object$control$epsilon, object$control$beta_init),
                     object$normalization_s0)
    }else{
      # object$arguments$alternative_specific <- object$arguments$alternative_specific[object$arguments$alternative_specific != tt]
      # if(length(object$arguments$alternative_specific) == 0) {object$arguments$alternative_specific <- NA}
      nfit <- Discrete_CM(formula = nfit,
                          data = object$data,
                          cdf = object$cdf,
                          alternative_specific = object$arguments$alternative_specific,
                          # object$categories_order,
                          case_id = object$arguments$case_id,
                          alternatives = object$arguments$alternatives,
                          reference = object$arguments$reference,
                          intercept = object$arguments$intercept,
                          control.glmcat(object$control$maxit,
                                         object$control$epsilon, object$control$beta_init),
                          normalization = object$normalization_s0)
    }

    ans[i+1L, ] <- AIC(nfit)
    nnew <- nobs_glmcat(nfit, use.fallback = TRUE)
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

#' Stepwise for glmcat models
#' @description Stepwise based on the AIC
#' @rdname step_glmcat
#' @param object a GLMcat model.
#' @param scope defines the range of models examined in the stepwise search (same as in the step function of the stats package). This should be either a single formula, or a list containing components upper and lower, both formulae.
#' @param direction the mode of the stepwise search.
#' @param trace to print the process information.
#' @param steps the maximum number of steps.
#' @export
step_glmcat <- function (object, scope,
                         direction = c("both", "backward",
                                       "forward"),
                         trace = 1, steps = 1000)
{
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

  object$terms <-   terms(formula(object$formula), data = object$data)

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
  n <- nobs_glmcat(object, use.fallback = TRUE)
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
      aod <- drop2(fit, scope$drop, scale = 0, trace = trace)
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
        aodf <- add2(fit, scope$add, scale = 0,
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

    # fit <- GLMcat(form1, object$data, object$ratio, object$cdf, object$parallel,
    #               object$categories_order, object$ref_category,
    #               object$threshold, control.glmcat(object$control$maxit, object$control$epsilon, object$control$beta_init),
    #               object$normalization_s0)

    fun_in = object$Function
    if(fun_in == "GLMcat"){
      fit <- GLMcat(form1, object$data, object$ratio, object$cdf, object$parallel,
                    object$categories_order, object$ref_category,
                    object$threshold, control.glmcat(object$control$maxit, object$control$epsilon, object$control$beta_init),
                    object$normalization_s0)
    }else{
      # object$arguments$alternative_specific <- object$arguments$alternative_specific[object$arguments$alternative_specific != tt]
      # if(length(object$arguments$alternative_specific) == 0) {object$arguments$alternative_specific <- NA}
      fit <- Discrete_CM(formula = form1,
                         data = object$data,
                         cdf = object$cdf,
                         alternative_specific = object$arguments$alternative_specific,
                         # object$categories_order,
                         case_id = object$arguments$case_id,
                         alternatives = object$arguments$alternatives,
                         reference = object$arguments$reference,
                         intercept = object$arguments$intercept,
                         control.glmcat(object$control$maxit,
                                        object$control$epsilon, object$control$beta_init),
                         normalization = object$normalization_s0)
    }

    fit$terms <-   terms(formula(fit$formula), data = fit$data)

    nnew <- nobs_glmcat(fit, use.fallback = TRUE)
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

