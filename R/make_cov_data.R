
#' Generate dataset for covariance estimation of irregular functional data
#'
#' @param data a data.frame
#' @param response_var character string indicating the name of the response in \code{data}.
#' @param id_var character string indicating the name of the id variable in \code{data}
#'               identifying the functional observations.
#' @param sep character string seperating argument/covariate names and their number 1 and 2.
#'
#' @return data.frame with data for covariance estimation
#' @export
#'
make_cov_data <- function(data, response_var, id_var = NULL, sep = "") {
  stopifnot(is.data.frame(data))
  covs <- setdiff(names(data), c(response_var, id_var))
  data <- split(data, if(is.null(id_var)) 0 else data[[id_var]])
  data <- lapply(data, function(x) {
    combs <- combn(1:nrow(x), 2)
          d <- data.frame(
           yy = x[[response_var]][combs[1,]] * x[[response_var]][combs[2,]]
           )
         names(d) <- c(
           paste(response_var, response_var, sep = sep)
         )
         for(i in 1:2)
          d[paste(covs, i, sep = sep)] <- x[combs[i,], covs]
         if(!is.null(id_var))
           d[[id_var]] <- rep(x[[id_var]][1], ncol(combs))
    d
  })
  do.call(rbind, data)
}
