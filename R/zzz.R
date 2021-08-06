## usethis namespace: start
#' @import Rcpp
# #' @export GLMadj
# #' @export GLMcum
# #' @export GLMseq
# #' @export GLMref
#' @export GLMcat
#' @export Discrete_CM
#' @export print.glmcat
#' @export vcov.glmcat
#' @export summary.glmcat
#' @export coef.glmcat
#' @export logLik.glmcat
#' @export student_glmcat
#' @export noncentralt_glmcat
#' @export glmcat_control
#' @useDynLib GLMcat, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom utils flush.console
#' @importFrom stats pchisq add1 model.frame model.matrix as.formula pnorm printCoefmat add.scope deviance drop.scope AIC factor.scope formula terms update update.formula
## usethis namespace: end

loadModule("GLMcatmodule", TRUE)
loadModule("Discrete_CMmodule", TRUE)
