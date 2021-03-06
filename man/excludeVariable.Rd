% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preProcessing.R
\name{excludeVariable}
\alias{excludeVariable}
\title{Determine whether to exclude a variable}
\usage{
excludeVariable(x, pmaxcut, sdmincut, k)
}
\arguments{
\item{x}{A variable}

\item{pmaxcut}{For discrete variables or those with few unique values
(set by the argument \code{k}), what is the highest proportion of
non-missing scores that can be on any one value before the variable
is excluded?}

\item{sdmincut}{For continuous variabbles, what is the smallest
standard deviation of non missing values allowed before the variable
is excluded?}

\item{k}{An integer indicating the number of unique values requires for
numeric or integer variables to be treated as discrete and the
\code{pmaxcut} threshold to be used instead of the standard
deviation threshold. Note that factor, character, and logical
class variables are always treated as discrete regardless of
how many unique values they have.}
}
\value{
A single logical value, \code{TRUE} or \code{FALSE}
  indicating whether to exclude the variable or not.
}
\description{
This function attempts to determine whether the variation or
non missing information in a variable is so low it should
be excluded.
}
\details{
Note that options for \dQuote{soc} can be set to control the default
arguments. Currently, the defaults are null but internally in functions,
if no options are set, they fall back to:
\code{options(soc = list(pmaxcut = .95, sdmincut = .01, k = 5))}
All or none of these options can be specified by the user to globally
control the behavior of the functions.#'
}
\examples{
excludeVariable(mtcars$mpg)
excludeVariable(mtcars$cyl)
excludeVariable(factor(mtcars$cyl))
excludeVariable(c(rep(1, 20), 2))
excludeVariable(c(rep(1, 20), 2), k = 1L)
options(soc = list(k = 1L))
excludeVariable(c(rep(1, 20), 2))
options(soc = NULL)
excludeVariable(c(rep(1, 20), 2))
excludeVariable(c(1, 1, 1, 1, 1, 1))
excludeVariable(1:10)
}
