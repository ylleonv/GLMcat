#' Variance-Covariance Matrix for a Fitted glmcat Model Object
#' @description Returns the variance-covariance matrix of the main parameters of a fitted \code{glmcat} model object.
#' @rdname anova
#' @method anova glmcat
#' @usage \method{anova}{glmcat}(object, object2, type, ...)
#' @export
anova.glmcat <-
  function(object, object2, type, ...)
  {
    mlist <- c(list(object),list(object2))
    AIC <- sapply(mlist, AIC)
    logLiks <- sapply(mlist, logLik)
    statistic <- c(NA, 2 * diff(sapply(mlist, logLik)))
    no.par <- sapply(mlist, function(x) x$`df of the model`)
    df <- c(NA, diff(no.par))
    pval <- c(NA, pchisq(statistic[-1], df[-1], lower.tail = FALSE))
    pval[!is.na(df) & df == 0] <- NA
    tab <- data.frame(no.par, AIC, logLiks, statistic, df, pval)
    tab.names <- c("no.par", "AIC", "logLik", "LR.stat", "df",
                   "Pr(>Chisq)")
    # mnames <- sapply(as.list(mc), Deparse)[-1]
    colnames(tab) <- tab.names
    # rownames(tab) <- rownames(models) <- mnames[ord]
    # colnames(models) <- models.names
    models <- sapply(mlist, function(x) x$formula)
    attr(tab, "models") <- models
    attr(tab, "heading") <- "Likelihood ratio tests of cumulative link models:\n"
    class(tab) <- c("anova.clm", "data.frame")
    return(tab)
  }
