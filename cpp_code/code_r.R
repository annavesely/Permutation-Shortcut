setwd('C:/Users/coshk/Desktop/Perm_project/cpp_code')
require(Rcpp)
Rcpp::sourceCpp('code_c.cpp')



# Given a matrix G of test statistics, the function returns
# D, matrix of centered test statistics where the rows are sorted so that
# the first one in ascending order
# R, matrix where each row is sorted in descending order
# I, matrix of the indices corresponding to the elements of D

ctrp_set <- function(G, higher=TRUE){
  f <- ncol(G)
  B <- nrow(G)
  
  if(!higher) G <- -1 * G
  
  # ordering according to the first row
  I_incr <- order(G[1,], decreasing=F)
  D <- G[,I_incr]
  
  # centered test statistics
  D <- sweep(D, 2, D[1,])
  
  # centered test statistics where each row is in ascending order
  o <- t(apply(D, 1, order, decreasing=T))
  R <- t(sapply(seq(B), function(x) D[x, o[x,]]))
  
  # matrix of indices in R
  I <- matrix(rep(I_incr,B), ncol=f, byrow=T)
  I <- t(sapply(seq(B), function(x) I[x, o[x,]]))
  I[1,] <- rev(I[1,]) # since they are all 0, we need to reverse the order
  
  out <- list("D"=D, "R"=R,"I"=I)
  return(out)
}




# Internal function
# Given a set S of indices, and the matrices D, R and I given by ctrp_set,
# m=ncol(D)-length(S), B=nrow(D) and s=length(S), it splits D, R and I and returns:
# ds, vector of the test statistics for S
# matrix D with only the indices not in S
# matrix R with only the indices not in S
# matrices Dsum and Rsum of the cumulative sums of ds with D and R

gen_sub <- function(S, D, R, I, f=ncol(D), m=ncol(D)-length(S), B=nrow(D), s=length(S)){
  
  if(m==0){
    ds <- rowSums(D)
    out <- list("ds"=ds, "D"=NULL, "R"=NULL, "I"=NULL)
    return(out)
  }
  
  i <- match(S,I[1,])
  i_D <- f+1-i
  Dc <- D[,-i_D]
  
  if(s==1){
    ds <- D[,i_D]
  }else{
    ds <- rowSums(D[,i_D])
  }
  
  # if m=1, then Rc=Dc and Ic are vectors
  if(m==1){
    Ic <- rep(I[1,-i], B)
    Rc <- Dc
  }else{
    Ic <- matrix(rep(NA, B*m), ncol=m)
    Rc <- Ic
    Ic[1,] <- I[1,-i]
    Rc[1,] <- R[1,-i]
    for(x in (2:B)){
      i <- match(S, I[x,])
      Ic[x,] <- I[x,-i]
      Rc[x,] <- R[x,-i]
    }
  }
  Dsum <- t(apply(cbind(ds, Dc), 1, cumsum))
  Rsum <- t(apply(cbind(ds, Rc), 1, cumsum))
  out <- list("ds"=ds, "D"=as.matrix(Dc), "R"=as.matrix(Rc), "I"=as.matrix(Ic),
              "Dsum"=Dsum, "Rsum"=Rsum)
  return(out)
}




# Given a vector of indices S, the matrices D, R, I as given by ctrp_set,
# the significance level alpha, the maximum number n_max of BAB iterations
# the condition from_low (T if the BAB starts from the lowest statistic)
# and the condition first_rem (T if the BAB starts by removing the statistic),
# the function tests S. It returns non_rej (T if there is a non-rejection,
# F if S is rejected, and NULL if the maximum number of iterations is reached before
# a decisive output) and BAB, the number of iterations.

ctrp_test <- function(S, D, R, I, alpha=0.05, n_max=10000, from_low=T, first_rem=T){
  f <- ncol(D)
  B <- nrow(D)
  s <- length(S)
  m <- f-s
  k <- ceiling((1-alpha)*B)
  
  # if S=F:
  if(m==0){
    # lower and upper bounds (equal)
    low <- Q(rowSums(D), k, B)
    out <- list("non_rej"=!low, "BAB"=0)
    return(out)
  }
  
  g <- gen_sub(S, D, R, I, f, m, B, s)
  out <- cpp_test(g$D, g$R, g$I, g$Dsum, g$Rsum, k, m, B, from_low, first_rem, n_max)
  return(out)
}



