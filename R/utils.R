table_summary <- function(object, ...) {
  # object <- mod_pol
  names <- c("ratio"," link", "nobs", "niter")
  info <- data.frame("ratio" = object$ratio,
               "cdf" = object$cdf[1],
               "nobs" = object$nobs,
               "niter" = paste("(", object$iteration, ")", sep=""),
               row.names = "Model info:"
    )
  return(info)
}

