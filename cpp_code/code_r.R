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




# ------------------------------------------------------------#
# TRUE DISCOVERY PROPORTION



# Internal function
# Given a set S of indices, and the matrices D, R and I given by ctrp_set,
# m=ncol(D)-length(S), B=nrow(D) and s=length(S), it splits D, R and I and returns:
# ds, vector of the test statistics for S
# matrix Ds with only the indices in S
# matrix D with only the indices not in S
# matrix Rs with only the indices in S
# matrix R with only the indices not in S
# matrices Dsum and Rsum of the cumulative sums of ds with D and R

# FINISH WITH CASES SUCH AS s=1 (as in gen_sub)
# IT WILL SUBSTITUTE gen_sub
# MAYBE WRITE IN C?

gen_sub2 <- function(S, D, R, I, f=ncol(D), m=ncol(D)-length(S), B=nrow(D), s=length(S)){
  
  Sord <- I[1,][I[1,]%in%S]
  i <- match(Sord, I[1,])
  iD <- rev(f+1-i)
  
  #i_D <- rev(i_D)
  Ds <- D[,iD]
  Dc <- D[,-iD]
  
  # if m=1, then Rc=Dc and Ic are vectors
  if(m==1){
    Ic <- rep(I[1,-i], B)
    Rc <- Dc
  }else{
    Ic <- matrix(rep(0, B*m), ncol=m)
    Rc <- Ic
    Is <- matrix(rep(0, B*s), ncol=s)
    Rs <- Is
    Ic[1,] <- I[1,-i]
    Is[1,] <- I[1,i]
    for(x in (2:B)){
      Sord <- I[x,][I[x,]%in%S]
      i <- match(Sord, I[x,])
      Ic[x,] <- I[x,-i]
      Rc[x,] <- R[x,-i]
      Is[x,] <- I[x,i]
      Rs[x,] <- R[x,i]
    }
  }
  ds <- rowSums(Ds)
  Dsum <- t(apply(cbind(ds, Dc), 1, cumsum))
  Rsum <- t(apply(cbind(ds, Rc), 1, cumsum))
  out <- list("ds"=ds, "D"=as.matrix(Dc), "R"=as.matrix(Rc), "I"=as.matrix(Ic),
              "Dsum"=Dsum, "Rsum"=Rsum, "Ds"=as.matrix(Ds), "Rs"=as.matrix(Rs),
              "Is"=as.matrix(Is))
}




# ------------------------------------------------------------#
# PLOTS FOR THE BOUNDS



# Given a vector X and a value k, it returns the critical value
# i.e. the k-th statistic (when sorted in increasing order)
Qu <- function(X, k){
  Xord <- sort(X, na.last = NA, decreasing=F, method="quick")
  return(Xord[k])
}




# Given the matrices Dsum and Rsum as returned by gen_sub,
# the index k and m=ncol(Dsum)-1 (with m>0),
# it returns the vectors of the lower and upper critical values
# for each superset size, from 0 to m
write_bounds <- function(Dsum, Rsum, k, m=ncol(Dsum)-1){
  low <- rep(NA,m+1)
  low[1] <- Qu(Dsum[,1], k)
  low[m+1] <- Qu(Dsum[,m+1], k)
  up <- low
  
  v <- 1
  while(v < m){
    v <- v+1
    low[v] <- Qu(Dsum[,v], k)
    up[v] <- Qu(Rsum[,v], k)
  }
  out <- list("low"=low, "up"=up)
  return(out)
}




# Given a vector of indices S, the matrices D, R, I as given by ctrp_set
# and the significance level alpha,
# it returns the vectors of the lower and upper critical values
# for each superset size, from 0 to m
ctrp_bounds <- function(S, D, R, I, alpha=0.05){
  f <- ncol(D)
  B <- nrow(D)
  s <- length(S)
  m <- f-s
  k <- ceiling((1-alpha)*B)
  
  # if S=F:
  if(m==0){
    # lower and upper bounds (equal)
    low <- Qu(rowSums(D), k)
    out <- list("low"=low, "up"=low)
    return(out)
  }
  
  g <- gen_sub(S, D, R, I, f, m, B, s)
  out <- write_bounds(g$Dsum, g$Rsum, k, m)
  return(out)
}




# Given the vectors of the lower and upper bounds
# and a boolean add_point (TRUE if points corresponding to observations
# should be added to the plot), it plots the lower bound in blue,
# the upper bound in red and the value zero in black
# A correction in applied for vectors having values <= -1e+155

require(ggplot2)
ctrp_plot <- function(low, up, add_point=FALSE){
  
  a <- min(0,low)
  b <- max(0,up)
  M <- max(abs(a), abs(b))
  div <- 10^(floor(log(M,10)))
  
  low <- low/div
  up <- up/div
  df <- data.frame(size=seq(length(low))-1, low=low, up=up, obs=0)
  
  a <- floor(a/div)
  b <- ceiling(b/div)
  c <- (a+b)/2
  
  gp <-  ggplot(df, aes()) +
    geom_line(aes(size, obs), linetype = "dashed", size=1) +
    geom_line(aes(size, low), colour="blue", size=1) +
    geom_line(aes(size, up), colour="red", size=1) +
    ylab("") + xlab("additional indices") +
    scale_y_continuous(limits=c(a,b), breaks=c(a,c,0,b), labels=c(a*div, c*div, 0, b*div))
  
  if(add_point){
    gp <- gp + 
      geom_point(aes(size, low), colour="blue", size=2) +
      geom_point(aes(size, up), colour="red", size=2)
  }
  return(gp)
}




# Given the matrices Dsum, Rsum, Ds and Rs as returned by gen_sub,
# the index k, s=ncol(Ds) and m=ncol(Dsum)-1,
# it returns the s x (m+1) matrices L and U
# where the z-th row is the vector of the lower/upper critical values
# for each superset size v (z=0,...,s and v=0,...,m)

# FINISH FOR m=0 and s=1

tdp_bounds <- function(Dsum, Rsum, Ds, Rs, k, s=ncol(Ds), m=ncol(Dsum)-1){
  # row z: bounds for v varying from 0 to m
  L <- matrix(rep(NA, s*(m+1)), ncol=m+1)
  U <- L
  
  Dsum_new <- Dsum
  Rsum_new <- Rsum
  w <- write_bounds(Dsum_new, Rsum_new, k, m)
  L[1,] <- w$low
  U[1,] <- w$up
  
  
  for(z in(1:(s-1))){
    Dsum_new <- Dsum_new - g$Ds[,s-z+1]
    Rsum_new <- Rsum_new - g$Rs[,s-z+1]
    w <- write_bounds(Dsum_new, Rsum_new, k, m)
    L[z+1,] <- w$low
    U[z+1,] <- w$up
  }
  out <- list("L"=L, "U"=U)
  return(out)
}




