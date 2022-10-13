
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
make_cov_data <- function(data, response_var, id_var = NULL, sep = "", diag = FALSE) {
  stopifnot(is.data.frame(data))
  covs <- setdiff(names(data), c(response_var, id_var))
  data <- split(data, id_var)
  data <- lapply(data, function(x) {
    combs <- combn(1:nrow(x), 2)
    with(x, {
         d <- data.frame(
           yy = x[[response_var]][combs[1,]] * x[[response_var]][combs[2,]],
           id1 = x[[id_var]][combs[1,]],
           id2 = x[[id_var]][combs[2,]]
           )
         names(d) <- c(
           paste(response_var, response_var, sep = sep),
           paste(id_var, 1:2, sep = sep)
         )
         for(i in 1:2)
          d[paste(covs, i, sep = sep)] <- d[combs[i,], covs]
      }
    )
  })
  do
}
