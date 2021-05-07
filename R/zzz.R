## usethis namespace: start
#' @import Rcpp
# #' @export GLMadj
# #' @export GLMcum
# #' @export GLMseq
# #' @export GLMref
#' @export GLMcat
#' @export Discrete_CM
# #' @export predict_glmcat_Response
#' @export summary.glmcat
#' @export coef.glmcat
#' @export nobs.glmcat
#' @export logLik.glmcat
#' @export normalization.glmcat
#' @export student.glmcat
#' @export noncentralt.glmcat
#' @export control.glmcat
#' @useDynLib GLMcat, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom utils flush.console
#' @importFrom stringr str_trim
#' @importFrom stats pchisq add1 model.frame model.matrix as.formula pnorm printCoefmat add.scope deviance drop.scope AIC factor.scope formula terms update update.formula
## usethis namespace: end

loadModule("GLMcatmodule", TRUE)
loadModule("Discrete_CMmodule", TRUE)
