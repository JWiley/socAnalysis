#' Check the type of a variable
#'
#' @param x A variable
#' @return A character string indicating the type of the variable,
#'   either \dQuote{c} for continuous,
#'   \dQuote{p} for ordered categorical,
#'   \dQuote{d} for unordered categorical,
#'   \dquote{cannot be classified} for variables that do not fit in.
#' @export
#' @examples
#' varType(mtcars$mpg)
#' varType(mtcars$cyl)
#' varType(factor(mtcars$cyl))
#' varType("test")
#' varType(factor("b"))
#' varType(factor(c("a", "b"), ordered = TRUE))
#' varType(TRUE)
varType <- function(x) {
  x <- na.omit(x)
  if ((length(unique(x)) > 5) && (is.integer(x) || is.numeric(x))) {
    "c"
  } else if (is.ordered(x)) {
    "p"
  } else if ((is.factor(x) && length(levels(x)) == 2) || length(unique(x)) == 2) {
    "d"
  } else {
    ("cannot be classified")
  }
}

#' Determine whether to exclude a variable
#'
#' This function attempts to determine whether the variation or
#' non missing information in a variable is so low it should
#' be excluded.
#'
#' Note that options for \dQuote{soc} can be set to control the default
#' arguments. Currently, the defaults are null but internally in functions,
#' if no options are set, they fall back to:
#' \code{options(soc = list(pmaxcut = .95, sdmincut = .01, k = 5))}
#' All or none of these options can be specified by the user to globally
#' control the behavior of the functions.#'
#'
#' @param x A variable
#' @param pmaxcut For discrete variables or those with few unique values
#'   (set by the argument \code{k}), what is the highest proportion of
#'   non-missing scores that can be on any one value before the variable
#'   is excluded?
#' @param sdmincut For continuous variabbles, what is the smallest
#'   standard deviation of non missing values allowed before the variable
#'   is excluded?
#' @param k An integer indicating the number of unique values requires for
#'   numeric or integer variables to be treated as discrete and the
#'   \code{pmaxcut} threshold to be used instead of the standard
#'   deviation threshold. Note that factor, character, and logical
#'   class variables are always treated as discrete regardless of
#'   how many unique values they have.
#' @export
#' @return A single logical value, \code{TRUE} or \code{FALSE}
#'   indicating whether to exclude the variable or not.
#' @examples
#' excludeVariable(mtcars$mpg)
#' excludeVariable(mtcars$cyl)
#' excludeVariable(factor(mtcars$cyl))
#' excludeVariable(c(rep(1, 20), 2))
#' excludeVariable(c(rep(1, 20), 2), k = 1L)
#' options(soc = list(k = 1L))
#' excludeVariable(c(rep(1, 20), 2))
#' options(soc = NULL)
#' excludeVariable(c(rep(1, 20), 2))
#' excludeVariable(c(1, 1, 1, 1, 1, 1))
#' excludeVariable(1:10)
excludeVariable <- function(x, pmaxcut, sdmincut, k) {
  x <- na.omit(x)
  if (missing(pmaxcut)) {
    pmaxcut <- getOption("soc")$pmaxcut
    ## back up default
    if (is.null(pmaxcut)) {
      pmaxcut <- .95
    }
  }
  if (missing(sdmincut)) {
    sdmincut <- getOption("soc")$sdmincut
    ## back up default
    if (is.null(sdmincut)) {
      sdmincut <- .01
    }
  }
  if (missing(k)) {
    k <- getOption("soc")$k
    ## back up default
    if (is.null(k)) {
      k <- 5L
    }
  }

  if (!length(x)) {
    out <- TRUE
  } else if (length(unique(x)) == 1) {
    out <- TRUE
  } else if (is.factor(x) || is.character(x) || is.logical(x) || (length(unique(x)) < k)) {
    if (any(prop.table(table(x)) > pmaxcut)) {
      out <- TRUE
    } else {
      out <- FALSE
    }
  } else if (is.numeric(x) || is.integer(x)) {
    if (sd(x) < .01) {
      out <- TRUE
    } else {
      out <- FALSE
    }
  }
  return(out)
}

#' Function to create a model matrix for PCA excluding variables as needed
#'
#' Excludes variables with too little variability, computes the design matrix,
#' and then again excludes any variables or dummy codes with too low variability.
#' See \code{excludeVariable} for details.
#'
#' @param dat A data frame to process
#' @param seed The random seed to use for imputation, if needed.
#'   Note that if data are imputed, a single imputation is done using
#'   predictive mean matching.
#' @return A list with the model matrix in the first position,
#'   the raw data with variables ultimately included (imputed once if needed),
#'   and the variable names in the second position.
#' @export
#' @examples
#' tmpd <- data.frame(
#'   mpg = 1:30,
#'   cyl = factor(rep(letters[1:3], times = c(20, 9, 1))))
#' cleanModelMatrix(tmpd)
#'
#' tmpd <- data.frame(
#'   mpg = c(1:15, NA, 17:30),
#'   cyl = factor(rep(letters[1:3], times = c(15, 10, 5))))
#' cleanModelMatrix(tmpd)
cleanModelMatrix <- function(dat, seed = 1234) {
  dat <- as.data.frame(dat)
  stopifnot(is.data.frame(dat))

  if (anyNA(dat)) {
    tempData <- mice(dat, m = 1,
                     maxit = 100,
                     method = 'pmm',
                     seed = seed)
    dat <- complete(tempData, action = 1L)
  }

  exclude1 <- apply(dat, 2, excludeVariable)
  dat <- dat[, !exclude1]
  vnames <- names(dat)

  f <- as.formula(sprintf("~ %s", paste(vnames, collapse = " + ")))
  x <- model.matrix(f, data = dat)
  exclude2 <- apply(x, 2, excludeVariable)
  vnames <- attr(terms(f), "term.labels")[unique(attr(x, "assign")[!exclude2])]
  assignment <- attr(x, "assign")
  x <- x[, !exclude2]
  return(list(
    X = x,
    Data = dat[, vnames],
    Variables = vnames,
    Assignment = assignment[!exclude2]))
}

#' Dimension reduction on data using PCAmethods
#'
#' This function uses PCA to reduce dimensionality by calculating
#' dimension scores and by listing the top variables for each
#' dimension.
#'
#' @param dat A data frame to process.
#' @param method A character string indicating the type of PCA: either \dQuote{svd}
#'   or \dQuote{nlpca}.
#' @param cumR2 A numeric value indicating the cumulative R2 value desired, used to
#'   determine the number of dimensions to extract or use.
#' @param targetK An integer indicating the target number of variables desired.
#'   When extracting raw variables, tries to (very) roughly match this target.
#' @export
#' @return A list.
#' @examples
#' dimensionReduction(
#'  dat = mtcars,
#'  method = "svd",
#'  cumR2 = .95,
#'  targetK = 6)
#' dimensionReduction(
#'  dat = mtcars,
#'  method = "nlpca",
#'  cumR2 = .95,
#'  targetK = 6)
dimensionReduction <- function(dat, method = c("svd", "nlpca"), cumR2 = .95, targetK = 20L) {
  ## source("https://bioconductor.org/biocLite.R")
  ## biocLite("pcaMethods")
  x <- cleanModelMatrix(dat)
  res <- pca(object = x$X,
             method = method,
             nPcs = ncol(x$X),
             center=TRUE,
             scale = "uv")
  dimensions <- min(which(res@R2cum  >= cumR2))
  dimensions <- 1L:dimensions
  cumR2 <- as.numeric(res@R2cum[max(dimensions)])

  if (method == "svd") {
  usetopK <- targetK %/% max(dimensions)
  topVars <- apply(abs(res@loadings[, dimensions]), 2, function(y) {
    ry <- sort(rank(-y))
    names(ry)[1:min(c(usetopK, length(ry)))]
  })
  topVars <- unique(as.vector(unlist(topVars)))
  topRealVars <- x$Variables[
    unique(x$Assignment[colnames(x$X) %in% topVars])]
  } else {
    topVars <- colnames(x$X)
    topRealVars <- x$Variables
  }

  out <- list(
    PCAScores = res@scores[, dimensions],
    Variables = topRealVars,
    ModelMatrix = x$X[, topVars],
    Data = x$Data[, topRealVars],
    cumR2 = cumR2,
    D = length(dimensions))

  return(out)
}
