
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
aB <- alpha*B


# CASE 1
f <- function(D, R, I){
  m <- ncol(D)
  B <- nrow(D)
  A <- D[,m]
  D_new <- D[,1:(m-1)]
  r <- I[1,1]
  I_new <- matrix(rep(NA, B*(m-1)), ncol=(m-1))
  R_new <- I_new
  
  for(x in seq(B)){
    i <- match(r, I[x,])
    I_new[x,] <- I[x,-i]
    R_new[x,] <- R[x,-i]
  }
  out <- list("D"=D_new, "R"=R_new, "I"=I_new, "A"=A)
  return(out)
}


# step 1: -1
P <- f(D_0, R_0, I_0)
Dsum <- t(apply(cbind(ds,P$D),1,cumsum))
Rsum <- t(apply(cbind(ds,P$R),1,cumsum))
c(Q(Rsum[,2],k), Q(Rsum[,3],k)) # upper bounds
c(Q(Dsum[,2],k), Q(Dsum[,3],k)) # lower bounds
A1 <- P$A
# step 2: -1,-4
P <- f(P$D, P$R, P$I)
Dsum <- t(apply(cbind(ds,P$D),1,cumsum))
Rsum <- t(apply(cbind(ds,P$R),1,cumsum))
c(Q(Rsum[,2],k), Q(Rsum[,3],k)) # upper bounds
c(Q(Dsum[,2],k), Q(Dsum[,3],k)) # lower bounds
# step 3: -1,+4
Dsum <- t(apply(cbind(ds+P$A,P$D),1,cumsum))
Rsum <- t(apply(cbind(ds+P$A,P$R),1,cumsum))
c(Q(Rsum[,1],k), Q(Rsum[,2],k)) # upper bounds
c(Q(Dsum[,1],k), Q(Dsum[,2],k)) # lower bounds
X <- Dsum[,1]
a <- (aB - length(X[X > 0]))/length(X[X == 0])
a # 0 - > not reject



# CASE 2
f <- function(D, R, I){
  m <- ncol(D)
  B <- nrow(D)
  A <- D[,1]
  D_new <- D[,2:m]
  r <- I[1,m]
  I_new <- matrix(rep(NA, B*(m-1)), ncol=(m-1))
  R_new <- I_new
  
  for(x in seq(B)){
    i <- match(r, I[x,])
    I_new[x,] <- I[x,-i]
    R_new[x,] <- R[x,-i]
  }
  out <- list("D"=D_new, "R"=R_new, "I"=I_new, "A"=A)
  return(out)
}


# step 1: -5
P <- f(D_0, R_0, I_0)
Dsum <- t(apply(cbind(ds,P$D),1,cumsum))
Rsum <- t(apply(cbind(ds,P$R),1,cumsum))
c(Q(Rsum[,2],k), Q(Rsum[,3],k)) # upper bounds
c(Q(Dsum[,2],k), Q(Dsum[,3],k)) # lower bounds
# indecisive: 1
A1 <- P$A
# step 2: -5,-2
P <- f(P$D, P$R, P$I)
Dsum <- t(apply(cbind(ds,P$D),1,cumsum))
Rsum <- t(apply(cbind(ds,P$R),1,cumsum))
Q(Rsum[,2],k) # upper bounds
Q(Dsum[,2],k) # lower bounds
X <- Dsum[,2]
a <- (aB - length(X[X > 0]))/length(X[X == 0])
a #0 - > not reject



# CASE 3
# step 1: +5
P <- f(D_0, R_0, I_0)
Dsum <- t(apply(cbind(ds + P$A,P$D),1,cumsum))
Rsum <- t(apply(cbind(ds + P$A,P$R),1,cumsum))
c(Q(Rsum[,1],k), Q(Rsum[,2],k)) # upper bounds
c(Q(Dsum[,1],k), Q(Dsum[,2],k)) # lower bounds
# indecisive: 2
A1 <- P$A
# step 2: +5,+2
P <- f(P$D, P$R, P$I)
Dsum <- t(apply(cbind(ds + P$A + A1,P$D),1,cumsum))
Rsum <- t(apply(cbind(ds + P$A + A1,P$R),1,cumsum))
Q(Rsum[,1],k) # upper bounds
Q(Dsum[,1],k) # lower bounds
# step 3: +5,-2
Dsum <- t(apply(cbind(ds + A1,P$D),1,cumsum))
Rsum <- t(apply(cbind(ds + A1,P$R),1,cumsum))
Q(Rsum[,2],k) # upper bounds
Q(Dsum[,2],k) # lower bounds
# from above: +5...
# from above: +5...



