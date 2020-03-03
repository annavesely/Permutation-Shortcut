
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
S <- c(1:5)
te <- ctrp_test(S, c$D, c$R, c$I, 0.20, 20)
te


# SAME, COLUMNS NOT SORTED

# matrix of global test statistics
G <- gt(n=20, f=5, B=10, beta0=0, beta=c(20,10,5,0,0))
c <- ctrp_set(G)

# set under testing
S <- c(3)
te <- ctrp_test(S, c$D, c$R, c$I, 0.20, 20, from_low=T, first_rem=T)
te





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

