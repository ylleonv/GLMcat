table_summary <- function(object, ...) {
  # object <- mod_pol
  names <- c("ratio"," link", "nobs", "niter", "logLik")
  info <- data.frame("ratio" = object$ratio,
               "cdf" = object$cdf[1],
               "nobs" = object$nobs,
               "niter" = paste("(", object$iteration, ")", sep=""),
               "logLik" = object$LogLikelihood,
               row.names = "Model info:"
    )
  return(info)
}

#' Summarising \code{glmcat} Model Fits
#' @description \code{print.summary} method for GLMcat objects.
#' @param object an object of class \code{"glmcat"}.
#' @rdname print.summary
#' @method print summary.glmcat
#' @export
print.summary.glmcat <-
  function(x, digits = max(3, getOption("digits") - 3),
           signif.stars = getOption("show.signif.stars"), ...)
  {
    # x <- glmcat_wine_3
    # cat("formula:", x$formula, fill=TRUE)
    print(x$formula)
    print(x$table_summary)

    printCoefmat(x$coefficients[, , drop=FALSE],
                 digits=digits, signif.stars=signif.stars,
                 has.Pvalue=TRUE, ...)

    if(!is.null(correl <- x$correlation)) {
      cat("\nCorrelation of Coefficients:\n")
      ll <- lower.tri(correl)
      correl[ll] <- format(round(correl[ll], digits))
      correl[!ll] <- ""
      print(correl[-1, -ncol(correl)], quote = FALSE, ...)
    }
    return(invisible(x))
  }
