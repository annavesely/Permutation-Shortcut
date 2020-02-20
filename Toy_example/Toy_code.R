
myseed <- 33
set.seed(myseed)



# Given eta, the function generates y according to the logistic model
# (+ vectorized function)
gen_y0 <- function(eta){
  set.seed(myseed)
  p <- exp(eta)/(1+exp(eta))
  y <- rbinom(n=1, size=1, prob=p)
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
c <- ctrp_set(G)

D <- c$D
R <- c$R
I <- c$I



#G <- gt(n=20, f=10, B=400, beta0=0, beta=c(20,0,5,0,0,0,0,0,0,0))
#c <- ctrp_set(G)

# set under testing
S <- c(3)
te <- ctrp_test(S, c$D, c$R, c$I, 0.20)
te

ds <- te$ds
Rsum_0 <- te$Rsum
Dsum_0 <- te$Dsum
I_0 <- te$I
D_0 <- te$D
R_0 <- te$R
I_0 <- te$I
indecisive_0 <- te$indecisive
alpha <- 0.2
k <- 8
m <- 4
f <- 5
s <- 1
B <- 10
