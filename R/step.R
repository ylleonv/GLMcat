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
#' @param trace to print the process information.XLASL
#' @param steps the maximum number of steps.
#' @rdname step
#' @method step glmcat
#' @usage \method{step}{glmcat}(object, scope, direction, trace, steps)
#' @exportS3Method
step.glmcat <-
  function (object, scope, direction = c("both", "backward", "forward"), trace = 1, steps = 1000)
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
    results <- step.results(models = models[seq(nm)], fit, object, usingCp)

    class(results) <- c("step.glmcat")
    return(results)
  }
