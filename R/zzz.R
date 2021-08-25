## usethis namespace: start
#' @import Rcpp
# #' @export GLMref
# #' @export Discrete_CM
# #' @export print.glmcat
# #' @export vcov.glmcat
# #' @export summary.glmcat
# #' @export coef.glmcat
# #' @export logLik.glmcat
# #' @export student_glmcat
# #' @export noncentralt_glmcat
# #' @export glmcat_control
#' @useDynLib GLMcat, .registration = TRUE
#' @importFrom Rcpp sourceCpp evalCpp
#' @importFrom utils flush.console
#' @importFrom stringr str_trim
#' @importFrom stats na.pass coef cov2cor logLik qnorm vcov step nobs pchisq add1 model.frame model.matrix as.formula pnorm printCoefmat add.scope deviance drop.scope AIC factor.scope formula terms update update.formula
## usethis namespace: end

# Rcpp::loadModule("GLMcatmodule", TRUE)
# Rcpp::loadModule("Discrete_CMmodule", TRUE)
