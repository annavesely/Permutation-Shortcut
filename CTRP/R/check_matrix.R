#' @title Condition check: matrix of finite numbers
#'
#' @description Internal function. It checks whether an object is a matrix of finite numbers.
#'
#' @param x an R object
#'
#' @return The function returns \code{TRUE} if \code{x} is a matrix and all its elements
#' are finite numbers.
#'
#' @author
#'
#' @examples
#' M <- matrix(c(0, 1, 2, 3), nrow=2)
#' check_matrix(M)
#' @export



check_matrix <- function(x){
  out <- is.matrix(x) & all(is.finite(x))
  return(out)
}
