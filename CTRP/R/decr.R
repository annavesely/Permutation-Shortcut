#' @title Decreasing check
#'
#' @description Internal function. When testing a subset of indices by closed testing,
#' the function checks the bounds for all the possible superset sizes, excluding
#' the lowest and the highest. The sizes are examined in decreasing order.
#'
#' @param Dfull numerical vector, sum of the centered test statistics of all the indices
#' (for each permutation)
#' @param Dc matrix of the centered test statistics of the indices not under testing,
#' where each row corresponds to a permutation
#' @param m number of indices not under testing
#' @param alpha significance level
#'
#' @return The function returns a list with the following objects:
#' \code{condition} (\code{FALSE} if a non-rejection has been found) and
#' \code{indecisive} (vector of sizes with indecisive outcome).
#'
#' @author
#'
#' @examples
#' x1 <- c(0, 1.38, 1.67, 1.16, -0.51, -0.03, 1.03, 0.15, 0.35, 0.12)
#' x2 <- c(0.00, 1.03, -1.49, 1.97, -1.12, 1.56, 2.12, 2.62, -0.18, 0.83)
#' x3 <- c(0, -0.79, -1.98, 0.11, 0.87, -0.22, -1.69, -0.36, -1.70, 1.36)
#' x4 <- c(0, 0.79, 1.68, 1.45, 1.44, 1.88, 2.42, 0.79, 2.38, 0.87)
#' Dc <- matrix(c(x1, x2, x3, x4), ncol=4)
#' Ds <- c(0, 0.14, 2.47, -1.35, -2.52, -0.90, -1.94, 0.05, 0.42, -0.78)
#' decr(Ds, Dc, m=4)
#'
#' Ds <- c(0, -4.86, -2.53, -6.35, -7.52, -5.90, -6.94, -4.95, -4.58, -5.78)
#' decr(Ds, Dc, m=4)
#' @export



decr <- function(Dfull, Dc, m=ncol(Dc), alpha=0.05){
  L <- Dfull
  U <- Dfull
  cond <- TRUE
  ind <- 1:(m-1)

  for(v in rev(ind)){
    L <- L - Dc[,m-v]
    # stop (not rejected) if 0 <= quantile(L)
    if(!above_quantile(L, alpha)){
      cond <- FALSE
      break
    }

    # v is not indecisive if 0 > quantile(U)
    U <- U - Dc[,v+1]
    if(above_quantile(U, alpha)){
      ind[v] <- NA
    }
  }

  if(!cond | all(is.na(ind))){ind <- NULL}
  else{ind <- na.omit(ind)}
  out <- list("condition"=cond, "indecisive"=ind)
  return(out)
}
