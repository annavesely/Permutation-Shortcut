#' @title Condition check: subset under test
#'
#' @description Internal function. It checks whether an object is a well-defined
#' vector of indices to be tested.
#'
#' @param x an R object
#' @param N maximum index
#'
#' @return The function returns \code{TRUE} if \code{x} is a vector, and all its elements
#' are integers greater or equal to 1 and smaller or equal to \code{N}.
#'
#' @author
#'
#' @examples
#' S <- c(5, 7, 1, 9)
#' check_hp(S, 20)
#' check_hp(S, 8)
#' S <- c(1.5, 7, 1, 9)
#' check_hp(S, 20)
#' @export



check_hp <- function(x, N){
  out <- is.vector(x) & all(floor(x)==x) & min(x)>=1 & max(x)<=N
  return(out)
}
