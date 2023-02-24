
#' Dynamic correlation structure
#'
#' @param value a \code{corStruct} object indicating the working correlation
#' structure assumed initially when computing the dynamic correlation structure
#' subclass requires pre-fitting the model. The default \code{NULL} indicates working independence.
#' @param formula a formula
#' @param fixed an optional logical value indicating whether the coefficients
#' should be allowed to vary in the optimization, or kept fixed at their initial value.
#' Defaults to \code{FALSE}, in which case the coefficients are allowed to vary.
#'
#' @import nlme
#' @export
#'
#' @examples
corDynamic <- function(value = 0, form = ~1, fixed = FALSE,
                       working_correlation = NULL, working_control = list()) {
  # Prepare 'sleeping' dynamic covariance structure
  attr(value, "formula") <- form
  attr(value, "fixed") <- fixed
  class(value) <- c(setdiff(class(value), c("corDynamic", "corStruct")),
                    c("corDynamic", "corStruct"))

  if(length(form) == 3) { # => two-sided formula
    # data for initialization of corDynamic should be given on left side
    return(value)
  } else { # initialize with working correlation first
    # ... and initially choose working correlation
    if(is.null(working_correlation)) stop("Working independence not implemented, yet.") else {
      stopifnot(is.list(working_control))
      if(!is.null(working_control$form))
        stop("Supplied formula is also taken for working correlation and no other formula can be specified.")
      working_control$form <- form
      working_correlation <- do.call(working_correlation, working_control)
      attr(working_correlation, "dynamic") <- value
      class(working_correlation) <- c("corDynamic_init", class(working_correlation))
      # append(class(working_correlation), "corDynamic_init",
      #                                    after = length(class(working_correlation))-1)
      return(working_correlation)
    }
  }
}


#' Initialization of a dynamic correlation structure
#'
#' Correlations structures of class \code{corDynamic} are not intended for
#' direct use but provide infrastructure for subclasses which require access to
#' objects which are not usually passed to it.
#' In particular, the initialization method provides the environment of the
#' \code{lme.formula} call executing it and to model residuals computed from
#' the last call.
#'
#' @export
#' @import nlme
#'
Initialize.corDynamic_init <- function(object, data, ...) {
  class(object) <- class(object)[-1]
  object <- Initialize(object, data, ...)
  class(object) <- c("corDynamic_init", class(object))

  # need to update formula for use in gamm()
  attr(attr(object, "dynamic"), "formula") <- formula(object)

  # Initialize also dynamic component
  attr(object, "dynamic") <- Initialize(attr(object, "dynamic"), data, ...)

  object
}

#' @export
#' @import nlme
#'
Initialize.corDynamic <- function(object, data, ...) {
  # which_corDyn <- which(class(object) == "corDynamic")
  # class(object) <- class(object)[-which_corDyn]
  # object <- Initialize(object, data, ...)
  # class(object) <- append(class(object), "corDynamic", which_corDyn-1)
  object <- NextMethod()

  # get first residuals (or other required data)
  if(length(formula(object)) == 3)
    attr(object, "residuals") <- getResponse(object, data = data)

  # catch lme environment
  lme_env <- parent.frame(3)
  er <- is.null(lme_env$.Method)
  if(!er) er <- lme_env$.Method != "lme.formula"
  if(er) {
    lme_env <- parent.frame(4)
    er <- is.null(lme_env$.Method)
    if(!er) er <- lme_env$.Method != "lme.formula"
  }
  if(er) {
    warning("corDynamic seems not to be executed from within lme.formula().
                 No environment and functionality for dynamic updates available.")
    return(object)
  }

  # prepare dynamic infrastructure
  fitted_ <- function(level = lme_env$Q) {
    # copied from nlme:::lme.formula
    ## fitted.values and residuals (in original order)
    Fitted <- matrix(fitted(lme_env$lmeSt, level = level,
                            conLin = if (lme_env$decomp) oldConLin else
                              attr(lme_env$lmeSt, "conLin")), ncol = length(level))[
                                lme_env$revOrder, , drop = FALSE]
    rownames(Fitted) <- lme_env$origOrder
    Fitted
  }

  residuals_ <- function(level = lme_env$Q) {
    Fitted <- fitted_(level)
    Resid <- lme_env$y[lme_env$revOrder] - Fitted
    Resid
  }

  attr(object, "lme_env") <- lme_env
  attr(object, "get_fitted") <- fitted_
  attr(object, "get_residuals") <- residuals_

  object
}


#' @export
#' @import nlme
#' @rdname nlme::needUpdate
needUpdate.corDynamic_init <- function(object) {
  f <- parent.frame(3)
  if(is.null(f$.Generic))
    return(TRUE)
  if(f$.Generic == "Initialize")
    return(FALSE)
  TRUE
}

#' @export
#' @import nlme
#' @rdname nlme::needUpdate
needUpdate.corDynamic <- function(object) needUpdate.corDynamic_init(object)

#' @export
#' @import nlme
update.corDynamic_init <- function(object, data) {
  new_object <- attr(object, "dynamic")
  attr(object, "dynamic") <- NULL
  attr(new_object, "working_correlation") <- object
  # ensure update iteration
  attr(new_object, "lme_env")$oldPars[] <- Inf
  # store residuals already (as they will be useful for basically any corDynamic)
  attr(new_object, "residuals") <- attr(new_object, "get_residuals")()

  coef(new_object) <- coef(new_object)
  new_object
}


#' @export
#' @import nlme
update.corDynamic <- function(object, data) {
  # store residuals already (as they will be useful for basically any corDynamic)
  attr(object, "residuals") <- attr(object, "get_residuals")()

  coef(object) <- coef(object)
  object
}


# #' @export
# #' @import nlme
# corFactor.corDynamic_init <- function(object, ...) {
#   class(object) <- setdiff(class(object), "corDynamic_init")
#   corFactor(object, ...)
# }

#' @export
#' @import nlme
#'
getCovariate.corDynamic <-
  getCovariate.corDynamic_init <- function(object, form = formula(object), data, ...) {

  # Copy of getCovariate.corStruct replacing getCovariate.data.frame by get_all_vars

  if (!missing(form)) {
    form <- formula(object)
    warning("cannot change 'form'")
  }

  if (is.null(covar <- attr(object, "covariate"))) {
    if (missing(data)) {
      stop("need data to calculate covariate of \"corStruct\" object")
    }
    covForm <- getCovariateFormula(form)
    # for covariance smoothing add also response to formula
    grps <- if (!is.null(getGroupsFormula(form)))
      getGroups(object, data = data)
    if (length(all.vars(covForm)) > 0) {
      if (is.null(grps)) {
        covar <- get_all_vars(formula = covForm, data = data)
      }
      else {
        if (all(all.vars(covForm) == sapply(splitFormula(covForm,
                                                         "+"), function(el) deparse(el[[2]])))) {
          covar <- split(get_all_vars(data = data, formula = covForm),
                         grps)
        }
        else {
          covar <- lapply(split(data, grps), get_all_vars,
                          formula = covForm)
        }
      }
    }
    else {
      if (is.null(grps)) {
        covar <- 1:nrow(data)
      }
      else {
        covar <- lapply(split(grps, grps), function(x) seq_along(x))
      }
    }
    if (!is.null(grps)) {
      covar <- as.list(covar)
    }
  }

  covar
}

# #' @export
# #' @import nlme
# #'
# getCovariate.corDynamic <- function(object, form = formula(object), data, ...) {
#   getCovariate.corDynamic_init(object = object, form = form, data = data, ...)
# }


#' @export
#' @import nlme
#'
getResponse <- function(object, form = formula(object), data, ...) {
  UseMethod("getResponse")
}

#' @export
#' @import nlme
#'
getResponse.corDynamic <- function(object, form = formula(object), data, ...) {

  # Copy of getCovariate.corStruct replacing getCovariateFormula by getResponseFormula

  if (!missing(form)) {
    form <- formula(object)
    warning("cannot change 'form'")
  }

  if (is.null(covar <- attr(object, "residuals"))) {
    if (missing(data)) {
      stop("need data to extract from \"corStruct\" object")
    }
    covForm <- getResponseFormula(form)
    # for covariance smoothing add also response to formula
    grps <- if (!is.null(getGroupsFormula(form)))
      getGroups(object, data = data)
    if (length(all.vars(covForm)) > 0) {
      if (is.null(grps)) {
        covar <- getCovariate(form = covForm, object = data)
      }
      else {
        if (all(all.vars(covForm) == sapply(splitFormula(covForm,
                                                         "+"), function(el) deparse(el[[2]])))) {
          covar <- split(getCovariate(object = data, form = covForm),
                         grps)
        }
        else {
          covar <- lapply(split(data, grps), getCovariate,
                          form = covForm)
        }
      }
    }
    else {
      if (is.null(grps)) {
        covar <- 1:nrow(data)
      }
      else {
        covar <- lapply(split(grps, grps), function(x) seq_along(x))
      }
    }
    if (!is.null(grps)) {
      covar <- as.list(covar)
    }
  }

  covar
}

