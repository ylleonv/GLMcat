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
#' @export nobs_glmcat
#' @export logLik.glmcat
#' @export normalization.glmcat
#' @export student.glmcat
#' @export noncentralt.glmcat
#' @export control.glmcat
# #' @export ReferenceF
#' @useDynLib GLMcat, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom utils flush.console
#' @import MASS add
#' @import stringr str_trim
#' @importFrom stats model.frame model.matrix as.formula pnorm printCoefmat add.scope deviance drop.scope AIC factor.scope formula terms update update.formula
## usethis namespace: end

loadModule("GLMcatmodule", TRUE)
loadModule("discretemodule", TRUE)
# loadModule("discretemodule", TRUE)
# loadModule("cumulativemodule", TRUE)
# loadModule("exportmod", TRUE)
# loadModule("sequentialmodule", TRUE)
# loadModule("adjacentmodule", TRUE)
# loadModule("referencemodule", TRUE) # predict_glmcatION
