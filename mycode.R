# Given a matrix of test statistics, the function returns
# M, matrix of centered test statistics where the rows are sorted so that
# the first one in ascending order
# D, matrix where each row is sorted in descending order
# I, matrix of the indices corresponding to the elements of D

ctrp_set <- function(G){
  f <- ncol(G)
  B <- nrow(G)
  
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
  I[1,] <- rev(I[1,])
  
  out <- list("D"=D, "R"=R,"I"=I)
  return(out)
}




# Internal function
# Given a vector X and a value k, it returns the critical value
# i.e. the k-th statistic (when sorted in increasing order)

Q <- function(X, k=ceiling(0.95*length(X))){
  Xord <- sort(X, na.last = NA, decreasing=F, method="quick")
  return(Xord[k])
}




# Internal function
# Given a vector X and a value k, it computes the lower critical value L
# and returns TRUE if L<0 (no non-rejection found)

Lcond <- function(X, k=ceiling(0.95*length(X))){
  L <- Q(X, k)
  c <- sign(L)==-1
  return(c)
}




# Internal function
# Given a vector X and a value k, it computes the upper critical value U
# and returns TRUE if U<=0 (no certain rejection found)

Ucond <- function(X, k=ceiling(0.95*length(X))){
  U <- Q(X, k)
  c <- sign(U)>-1
  return(c)
}




# Internal function
# Given a set S of indices, and the matrices D, R and I given by ctrp_set,
# it splits D, R and I, returning:
# ds, vector of the test statistic for S
# matrix D with only the indices not in S
# matrix R with only the indices not in S


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
  
  out <- list("ds"=ds, "D"=Dc, "R"=Rc, "I"=Ic)
  return(out)
}



# Internal function.
# Given ds and the complementary matrices D and R given by gen_sub, finds
# it finds the first column R[,j] such that all the el. are negative/null
# It checks the bounds up to j-1, and then the following bounds
# until one upper bound becomes negative.

compute_bounds <- function(indecisive, R, Dsum, Rsum, k=ceiling(0.95*nrow(R)), m=ncol(R)){
  
  H <- length(indecisive)
  up <- rep(F, H) # keeps track of indecisive sizes
  
  h <- 1
  v <- indecisive[1]
  low <- Lcond(Dsum[,v+1], k)
  up[1] <- Ucond(Rsum[,v+1], k)
  
  # first column in R with no positive element
  # (if m=1, then R is a vector)
  if(m==1){
    j <- Find(function(x) sign(max(R))<1 , 1)
  }else{
    j <- Find(function(x) sign(max(R[,x]))<1 , seq(m))
  }
  
  # if there is no such column
  if(length(j)==0){
    j <- m+1
  }
  
  # bounds for v < j (before the upper c.v. decreases)
  while(low & v<j-1 & h<H){
    h <- h+1
    v <- indecisive[h]
    low <- Lcond(Dsum[,v+1], k)
    up[h] <- Ucond(Rsum[,v+1], k)
  }
  
  # bounds for v >= j (until the upper c.v. becomes negative)
  while(low & up[h] & h<H){
    h <- h+1
    v <- indecisive[h]
    low <- Lcond(Dsum[,v+1], k)
    up[h] <- Ucond(Rsum[,v+1], k)
  }

  out <- list("non_rej"=!low, "indecisive"=indecisive[up])
  return(out)
}





# Internal function.
# Given ds and the complementary matrices D and R given by gen_sub, finds
# it finds the first column R[,j] such that all the el. are negative/null
# It computes the upper bounds up to j-1, and then the following bounds
# until one becomes negative.

compute_bounds_fast <- function(indecisive, R, Dsum, Rsum, k=ceiling(0.95*nrow(R)), m=ncol(R)){
  
  H <- length(indecisive)
  up <- rep(F, H) # keeps track of indecisive sizes
  
  h <- 1
  v <- indecisive[1]
  up[1] <- Ucond(Rsum[,v+1], k)
  
  # first column in R with no positive element
  j <- Find(function(x) sign(max(R[,x]))<1 , seq(m))
  
  # if there is no such column
  if(length(j)==0){
    j <- m+1
  }
  
  # upper bounds for v < j (before the upper c.v. decreases)
  while(v<j-1 & h<H){
    h <- h+1
    v <- indecisive[h]
    up[h] <- Ucond(Rsum[,v+1], k)
  }
  
  # bounds for v >= j (until the upper c.v. becomes negative)
  while(up[h] & h<H){
    h <- h+1
    v <- indecisive[h]
    up[h] <- Ucond(Rsum[,v+1], k)
  }
  
  out <- indecisive[up]
  return(out)
}




# Given D, R, I as given by ctrp_set,
# a vector of indices S, the significance level alpha
# and the maximum number n_max of BAB iterations (0 if no BAB),
# the function tests S and returns:
# non_rej, TRUE if a non-rejection has been found,
# indecisive, vector of indecisive sizes (not null if the BAB iterations stop
# before a decisive outcome)
# BAB, the number of BAB iterations

ctrp_test <- function(S, D, R, I, alpha=0.05, n_max=10000){
  f <- ncol(D)
  B <- nrow(D)
  s <- length(S)
  m <- f-s
  k <- ceiling((1-alpha)*B)
  
  # if S=F:
  if(m==0){
    # lower and upper bounds (equal)
    low <- Lcond(rowSums(D), k)
    out <- list("non_rej"=!low, "BAB"=0)
    return(out)
  }

  g <- gen_sub(S, D, R, I, f, m, B, s)
  
  Dsum <- t(apply(cbind(g$ds, g$D), 1, cumsum))
  Rsum <- t(apply(cbind(g$ds, g$R), 1, cumsum))
  indecisive <- (0:m)

  cb <- compute_bounds(indecisive, g$R, Dsum, Rsum, k, m)
  
  # if there is a non-rejection or there are no indecisives
  if(cb$non_rej | length(cb$indecisive)==0){
    out <- list("non_rej"=cb$non_rej, "BAB"=0)
    return(out)
  }
  
  #out <- ctrp_bab(cb$indecisive, g$D, g$R, g$I, Dsum, Rsum, k, m, B, n_max)
  
  #out <- list("non_rej"=cb$non_rej, "indecisive"=cb$indecisive)
  
  out <- list("non_rej"=cb$non_rej, "indecisive"=cb$indecisive, "ds"=g$ds, "D"=g$D, "R"=g$R, "I"=g$I, "Dsum"=Dsum, "Rsum"=Rsum)
 return(out)
}





