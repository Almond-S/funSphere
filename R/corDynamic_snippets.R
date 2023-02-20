
#' Title
#'
#' @param object
#' @param ...
#'
#' @import nlme
#' @export
#'
corMatrix.corDynamic <- function(object, ...) {
  if(attr(object, "initialized")) {
    stop("corDynamic is only providing infrastructure for dynamic correlation structures.
         No corDynamic-subclass specified which would implement particular correlation.")
  } else {
    return(corMatrix(attr(object, "working_correlation")))
  }
}

#' Title
#'
#' @param object
#' @param ...
#'
#' @import nlme
#' @export
corFactor.corDynamic <- function(object, ...) {
  if(attr(object, "initialized")) {
    stop("corDynamic is only providing infrastructure for dynamic correlation structures.
         No corDynamic-subclass specified which would implement particular correlation.")
  } else {
    return(corFactor(attr(object, "working_correlation")))
  }
}

#' Title
#'
#' @param object
#' @param ...
#'
#' @import nlme
#' @export
coef.corDynamic <- function(object, ...) {
  if(attr(object, "initialized")) {
    stop("corDynamic is only providing infrastructure for dynamic correlation structures.
         No corDynamic-subclass specified which would implement particular correlation.")
  } else {
    return(corFactor(attr(object, "working_correlation")))
  }
}

#' @export
#' @import nlme
#' @rdname nlme::recalc
recalc.corDynamic <- function(object, conLin, ...) {
  if(attr(object, "initialized")) {
    stop("corDynamic is only providing infrastructure for dynamic correlation structures.
         No corDynamic-subclass specified which would implement particular correlation.")
  } else {
    return(recalc(attr(object, "working_correlation")))
  }
}
