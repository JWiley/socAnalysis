% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preProcessing.R
\name{varType}
\alias{varType}
\title{Check the type of a variable}
\usage{
varType(x)
}
\arguments{
\item{x}{A variable}
}
\value{
A character string indicating the type of the variable,
  either \dQuote{c} for continuous,
  \dQuote{p} for ordered categorical,
  \dQuote{d} for unordered categorical,
  \dQuote{cannot be classified} for variables that do not fit in.
}
\description{
Check the type of a variable
}
\examples{
varType(mtcars$mpg)
varType(mtcars$cyl)
varType(factor(mtcars$cyl))
varType("test")
varType(factor("b"))
varType(factor(c("a", "b"), ordered = TRUE))
varType(TRUE)
}
