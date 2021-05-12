## usethis namespace: start
#' @import Rcpp
# #' @export GLMadj
# #' @export GLMcum
# #' @export GLMseq
# #' @export GLMref
#' @export GLMcat
#' @export Discrete_CM
#' @export summary.glmcat
#' @export coef.glmcat
#' @export nobs_glmcat
#' @export logLik.glmcat
#' @export student.glmcat
#' @export noncentralt.glmcat
#' @export control.glmcat
#' @useDynLib GLMcat, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom utils flush.console
#' @importFrom stats pchisq add1 model.frame model.matrix as.formula pnorm printCoefmat add.scope deviance drop.scope AIC factor.scope formula terms update update.formula
## usethis namespace: end

loadModule("GLMcatmodule", TRUE)
loadModule("Discrete_CMmodule", TRUE)
