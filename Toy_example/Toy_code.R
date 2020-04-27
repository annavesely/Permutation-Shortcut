
myseed <- 33
set.seed(myseed)



# Given eta, the function generates y according to the logistic model
# (+ vectorized function)
gen_y0 <- function(eta){
  set.seed(myseed)
  p <- exp(eta)/(1+exp(eta))
  y <- rbinom(1, 1, prob=p)
  return(y)
}

gen_y <- Vectorize(gen_y0)




# Given the number of observations n,
# the number of predictors f,
# and the number of permutations,
# the intercept beta0 and the vector of coefficients beta (with length f),
# the function simulates the covariate matrix X from a standard normal distrnution.
# Then it computes the matrix of the global test statistics
# (obtained by randomly permuting y)

gt <- function(n, f, B, beta0, beta){
  set.seed(myseed)
  X <- matrix(rnorm(n*f,0,1), ncol=f)
  eta <- beta0 + X %*% beta
  y <- gen_y(eta)
  H <- matrix(rep(1/n, n*n), ncol=n)

  g_mat <- matrix(rep(NA,B*f), ncol=f)
  g_mat[1,] <- (t(y) %*% (diag(n)-H) %*% X)^2
  for(b in(2:B)){
    y_perm <- sample(y, n, replace=FALSE)
    g_mat[b,] <- (t(y_perm) %*% (diag(n)-H) %*% X)^2
  }
  return(g_mat)
}




# -------------------------------------------------------------------------------- #


# TOY EXAMPLE

# matrix of global test statistics
G <- gt(n=20, f=5, B=10, beta0=0, beta=c(20,10,5,0,0))
G <- G[,c(1,4,2,5,3)]
c <- ctrp_set(G)

# set under testing
S <- c(5)
te <- ctrp_test(S, c$D, c$R, c$I, 0.20, 20, F, T)
te


# SAME, COLUMNS NOT SORTED

# matrix of global test statistics
G <- gt(n=20, f=5, B=10, beta0=0, beta=c(20,10,5,0,0))
c <- ctrp_set(G)

# set under testing
S <- c(3)
te <- ctrp_test(S, c$D, c$R, c$I, 0.20, 20, from_low=T, first_rem=T)
te




# TDP EXAMPLE
f <- 10
B <- 10
G <- gt(n=20, f=10, B=10, beta0=0, beta=c(30,20,10,5,5,5,5,5,5,5))
G <- G[,c(10,8,6,5,7,9,2,4,1,3)]
c <- ctrp_set(G)

S <- c(7,8,9,10)
s <- length(S)
m <- f-s
te <- ctrp_test(S, c$D, c$R, c$I, 0.20, 20, F, T)
te


#G <- G[,c(4,6,5,7,1,3,2)]
#G <- G[,c(3,2,5,7,1,4,6)]
#G <- G[,c(1,5,7,3,4,6,2)]

#S <- c(5,6,7)



gen_sub2 <- function(S, D, R, I, f=ncol(D), m=ncol(D)-length(S), B=nrow(D), s=length(S)){
  
  i <- match(S,I[1,])
  i_D <- f+1-i
  i_D <- rev(i_D)
  Ds <- D[,i_D]
  Dc <- D[,-i_D]
  
  # if m=1, then Rc=Dc and Ic are vectors
  if(m==1){
    Ic <- rep(I[1,-i], B)
    Rc <- Dc
  }else{
    Ic <- matrix(rep(NA, B*m), ncol=m)
    Rc <- Ic
    Is <- matrix(rep(NA, B*s), ncol=s)
    Rs <- Is
    Ic[1,] <- I[1,-i]
    Rc[1,] <- R[1,-i]
    Is[1,] <- I[1,i]
    Rs[1,] <- R[1,i]
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


g <- gen_sub2(S, c$D, c$R, c$I, f, m, B, s)







alpha <- 0.2
k <- ceiling(alpha * B)



# Given a vector X and a value k, it returns the critical value
# i.e. the k-th statistic (when sorted in increasing order)
Qu <- function(X, k){
  Xord <- sort(X, na.last = NA, decreasing=F, method="quick")
  return(Xord[k])
}


tdp_bounds <- function(Dsum, Rsum, Ds, Rs, s=ncol(Rs), m=ncol(Dsum)-1, k){
  
  z <- 0
  cond <- TRUE
  Dsum_z <- Dsum
  Rsum_z <- Rsum
  
  low <- rep(NA,m+1)
  up <- low
  
  while(cond & z<s){
    z <- z+1
    Dsum_z <- Dsum_z - Ds[,s-z+1]
    Rsum_z <- Rsum_z - Rs[,s-z+1]
    low[1] <- Qu(Dsum_new[,1], k)
    up[1] <- low[1]
    low[m+1] <- Qu(Dsum_new[,m+1], k)
    up[m+1] <- low[m+1]
    
    for(v in(2:m)){
      low[v] <- Qu(Dsum_new[,v], k)
      up[v] <- Qu(Rsum_new[,v], k)
    }
  }
  
  return(list("low"=low, "up"=up))
}



# z=1
Dsum_new <- g$Dsum - g$Ds[,s]
Rsum_new <- g$Rsum - g$Rs[,s]


low <- rep(NA,m+1)
low[1] <- Qu(Dsum_new[,1], k)
up <- low

for(v in(2:(m+1))){
  low[v] <- Qu(Dsum_new[,v], k)
  up[v] <- Qu(Rsum_new[,v], k)
}

round(up,2)
round(low,2)



# z=2
Dsum_new <- Dsum_new - g$Ds[,s-1]
Rsum_new <- Rsum_new - g$Rs[,s-1]


low <- rep(NA,m+1)
low[1] <- Qu(Dsum_new[,1], k)
up <- low

for(v in(2:(m+1))){
  low[v] <- Qu(Dsum_new[,v], k)
  up[v] <- Qu(Rsum_new[,v], k)
}


round(up,2)
round(low,2)






### -------------------------------------- ###


mycheck <- function(S, D, R, I, alpha, n_max=20){
  
  lr <- ctrp_test(S, D, R, I, alpha, n_max, from_low=T, first_rem=T) # remove low
  lk <- ctrp_test(S, D, R, I, alpha, n_max, from_low=T, first_rem=F) # keep low
  hr <- ctrp_test(S, D, R, I, alpha, n_max, from_low=F, first_rem=T) # remove high
  hk <- ctrp_test(S, D, R, I, alpha, n_max, from_low=F, first_rem=F) # keep high
  
  cond <- (lr$non_rej == lk$non_rej & lk$non_rej == hr$non_rej & hr$non_rej == hk$non_rej)
  babs <- c(lr$BAB, lk$BAB, hr$BAB, hk$BAB)
  out <- list("non_rej"=lr$non_rej, "BABs"=babs, "cond"=cond)
  return(out)
}


G <- gt(n=20, f=5, B=10, beta0=0, beta=c(20,10,5,0,0))
c <- ctrp_set(G)

mycheck(3, c$D, c$R, c$I, 0.20, 20) # 2,4,3,3
mycheck(1, c$D, c$R, c$I, 0.20, 20) # all 0
mycheck(c(3,4), c$D, c$R, c$I, 0.20, 20) # all 0
mycheck(c(1,4), c$D, c$R, c$I, 0.20, 20) # all 0
mycheck(c(1,2,3,4), c$D, c$R, c$I, 0.20, 20) # all 0
mycheck(c(1,2,3,4,5), c$D, c$R, c$I, 0.20, 20) # all 0

G <- gt(n=20, f=10, B=400, beta0=0, beta=c(20,0,5,0,0,0,0,0,0,0))
G <- G[,c(9,4,5,3,10,6,2,8,7,1)]
c <- ctrp_set(G)
mycheck(3, c$D, c$R, c$I, 0.20, 20)
mycheck(c(2,10), c$D, c$R, c$I, 0.20, 50) # 13,13,4,4
mycheck(c(3,4), c$D, c$R, c$I, 0.20, 20)
mycheck(c(1,4), c$D, c$R, c$I, 0.20, 20)
mycheck(c(1,2,3,4), c$D, c$R, c$I, 0.20, 20)
mycheck(c(1,2,3,4,5), c$D, c$R, c$I, 0.20, 20)

mycheck(10, c$D, c$R, c$I, 0.20, 50) #34,34,14,14
mycheck(c(3,10), c$D, c$R, c$I, 0.20, 50) # 16,16,6,6
#G <- gt(n=20, f=5, B=100, beta0=0, beta=c(100,1,0,0,0))

