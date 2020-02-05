#' @title Extremes check
#'
#' @description Internal function. When testing a subset of indices by closed testing,
#' the function checks the bounds for the lowest and highest superset sizes
#' (where only the indices of the subset and all the indices are considered, respectively).
#'
#' @param Ds numerical vector, sum of the centered test statistics of the indices
#' under testing (for each permutation)
#' @param Dfull numerical vector, sum of the centered test statistics of all the indices
#' (for each permutation)
#' @param alpha significance level
#'
#' @return The function returns a list with the following objects:
#' \code{condition} (\code{FALSE} if a non-rejection has been found) and
#' \code{order_incr} (\code{TRUE} if the rejection for the lowest size is the least extreme,
#' and thus the other sizes should be explored in increasing order).
#'
#' @author
#'
#' @examples
#' Dc_sum <- c(0.00, 2.41, -0.12, 4.69, 0.68, 3.19, 3.88, 3.20, 0.85, 3.18)
#' Ds <- c(0, 0.14, 2.47, -1.35, -2.52, -0.90, -1.94, 0.05, 0.42, -0.78)
#' Dfull <- Ds + Dc_sum
#' extremes(Ds, Dfull)
#'
#' Ds <- c(0, -4.86, -2.53, -6.35, -7.52, -5.90, -6.94, -4.95, -4.58, -5.78)
#' Dfull <- Ds + Dc_sum
#' extremes(Ds, Dfull)
#' @export


extremes <- function(Ds, Dfull, alpha=0.05){
  c0 <- quantile(Ds, 1-alpha, names=FALSE)
  cfull <- quantile(Dfull, 1-alpha, names=FALSE)

  # stop (not rejected) if 0 <= c0 or 0 <= cfull
  if(sign(c0)>-1 | sign(cfull)>-1){
    cond <- FALSE
    ord <- NULL
  }
  else{
    cond <- TRUE
    ord <- (cfull <= c0)
  }
  out <- list("condition"=cond, "order_incr"=ord)
  return(out)
}
