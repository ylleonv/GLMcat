\name{step.glmcat}
\alias{step.glmcat}
\title{Stepwise for a \code{glmcat} model object}
\usage{
\method{step}{glmcat}(object, scope, scale, direction, trace, keep, steps, k, ...)
}
\arguments{
\item{object}{an fitted object of class \code{glmcat}.}
\item{scope}{defines the range of models examined in the stepwise search (same as in the step function of the stats package). This should be either a single formula, or a list containing components upper and lower, both formulae.}
\item{scale}{the scaling parameter (if applicable).}
\item{direction}{the mode of the stepwise search.}
\item{trace}{to print the process information.}
\item{keep}{a logical value indicating whether to keep the models from all steps.}
\item{steps}{the maximum number of steps.}
\item{k}{additional arguments (if needed).}
\item{\dots}{additional arguments passed to the function.}
}
\description{
Stepwise for a \code{glmcat} model object based on the AIC.
}
