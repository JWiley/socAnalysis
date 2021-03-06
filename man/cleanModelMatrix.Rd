% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preProcessing.R
\name{cleanModelMatrix}
\alias{cleanModelMatrix}
\title{Function to create a model matrix for PCA excluding variables as needed}
\usage{
cleanModelMatrix(dat, seed = 1234)
}
\arguments{
\item{dat}{A data frame to process}

\item{seed}{The random seed to use for imputation, if needed.
Note that if data are imputed, a single imputation is done using
predictive mean matching.}
}
\value{
A list with the model matrix in the first position,
  the raw data with variables ultimately included (imputed once if needed),
  and the variable names in the second position.
}
\description{
Excludes variables with too little variability, computes the design matrix,
and then again excludes any variables or dummy codes with too low variability.
See \code{excludeVariable} for details.
}
\examples{
tmpd <- data.frame(
  mpg = 1:30,
  cyl = factor(rep(letters[1:3], times = c(20, 9, 1))))
cleanModelMatrix(tmpd)

tmpd <- data.frame(
  mpg = c(1:15, NA, 17:30),
  cyl = factor(rep(letters[1:3], times = c(15, 10, 5))))
cleanModelMatrix(tmpd)
}
