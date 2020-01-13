#' @title Sign of the sample quantile
#'
#' @description Internal function. It checks whether zero is greater than the sample quaantile
#' of a vector.
#'
#' @param x numerical vector
#' @param alpha significance level
#'
#' @return The function returns \code{TRUE} if zero is greater than the
#' (1-\code{alpha})-quantile of \code{x}.
#'
#' @author
#'
#' @examples
#' Y <- rnorm(100, mean=-2, sd=1)
#' hist(Y)
#' above_quantile(Y, 0.05)
#'
#' Y <- rnorm(100, mean=0, sd=1)
#' hist(Y)
#' above_quantile(Y, 0.05)
#' @export



above_quantile <- function(x, alpha=0.05){
  c <- quantile(x, 1-alpha, names=FALSE)
  out <- sign(c)==-1
  return(out)
}



# TRUE if zero is greater than the sample quantile of a vector X

