
#' Title
#'
#' @param data
#' @param response_var
#' @param id_var
#' @param var_sep
#'
#' @return
#' @export
#'
#' @examples
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
