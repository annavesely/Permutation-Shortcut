#' @title Condition check: matrix of indices
#'
#' @description Internal function. It checks whether a matrix is a well-defined
#' matrix of indices with given dimensions.
#'
#' @param x matrix of finite numbers
#' @param R number of rows
#' @param C number of columns
#'
#' @return The function returns \code{TRUE} if \code{x} has \code{R} rows
#' and \code{C} columns, and each row contains an unordered sequence of
#' integers from 1 to \code{C}.
#'
#' @author
#'
#' @examples
#' M <- matrix(c(1, 3, 2, 2, 3, 1, 2, 1, 3, 2, 3, 1), ncol=3, byrow=TRUE)
#' check_indices(M, R=4, C=3)
#' M <- matrix(c(9, 3, 2, 2, 3, 1, 2, 1, 3, 2, 3, 1), ncol=3, byrow=TRUE)
#' check_indices(M, R=4, C=3)
#' @export


check_indices <- function(x, R, C){
  check_row <- function(r){
    cond <- length(setdiff(r, (1:C))) == 0
    return(cond)
  }
  out <- nrow(x)==R & ncol(x)==C & all(apply(x, 1, check_row))
  return(out)
}
