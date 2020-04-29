
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
bo <- ctrp_bounds(S, c$D, c$R, c$I, 0.20)
ctrp_plot(bo$low, bo$up, TRUE)

S <- c(4,5)
s <- length(S)
f <- 5
m <- f-s
B <- 10
k <- 8
g <- gen_sub2(S, c$D, c$R, c$I, f, m, B, s)
td <- tdp_bounds(g$Dsum, g$Rsum, g$Ds, g$Rs, k, s, m)
ctrp_plot(td$L[1,], td$U[1,], TRUE)









# TDP EXAMPLE
G0 <- gt(n=20, f=10, B=10, beta0=0, beta=c(30,30,20,10,10,5,5,5,5,5))
G0 <- G0[,c(10,1,4,3,7,9,2,6,8,5)]
G <- G0[,c(1,3,4,6,7,9,10)]

f <- 7
B <- 10
alpha <- 0.2
k <- ceiling((1-alpha) * B)

c <- ctrp_set(G)
S <- c(5,6,7)
s <- length(S)
m <- f-s
te <- ctrp_test(S, c$D, c$R, c$I, 0.20, 20, F, T)
te
bo <- ctrp_bounds(S, c$D, c$R, c$I, 0.20)
ctrp_plot(bo$low, bo$up, TRUE)

g <- gen_sub2(S, c$D, c$R, c$I, f, m, B, s)
td <- tdp_bounds(g$Dsum, g$Rsum, g$Ds, g$Rs, k, s, m)

ctrp_plot(td$L[2,], td$U[2,], TRUE) # z=1, BAB needed
ctrp_plot(td$L[3,], td$U[3,], TRUE) # z=2, non-rejection


round(td$U[2,],2)
round(td$L[2,],2)


#BAB - removal of highest not in S

new_node <- function(dfixed, D0, R0, I0, m0){
  m <- m0 - 1
  D <- D0[,-m0]
  I <- matrix(rep(0, B*m), ncol=m)
  R <- I
  e <- I0[1,1]
  I[1,] <- I0[1,-1]
  for(x in (2:B)){
    i <- match(e, I0[x,])
    I[x,] <- I0[x,-i]
    R[x,] <- R0[x,-i]
  }
  Dsum <- t(apply(cbind(dfixed, D), 1, cumsum))
  Rsum <- t(apply(cbind(dfixed, R), 1, cumsum))
  out <- list("dfixed"=dfixed, "D"=D, "R"=R, "I"=I, "Dsum"=Dsum, "Rsum"=Rsum)
  return(out)
}

# check (2,3) -> rejection
h <- 1
N1 <- new_node(g$ds, g$D, g$R, g$I, m)
td <- tdp_bounds(N1$Dsum, N1$Rsum, g$Ds, g$Rs, k, s, m-1)
ctrp_plot(td$L[2,], td$U[2,], TRUE)
round(td$U[2,],2)
round(td$L[2,],2)

# dx: check (1,2) -> (1,2)
N1$Dsum <- N1$Dsum + g$D[,m-h+1]
N1$Rsum <- N1$Rsum + g$D[,m-h+1]
N1$dfixed <- N1$dfixed + g$D[,m-h+1]
td <- tdp_bounds(N1$Dsum, N1$Rsum, g$Ds, g$Rs, k, s, m-h)
ctrp_plot(td$L[2,], td$U[2,], TRUE)
round(td$U[2,],2)
round(td$L[2,],2)

# check (1,2) -> rejection
h <- 2
N2 <- new_node(N1$dfixed, N1$D, N1$R, N1$I, m-h+1)
td <- tdp_bounds(N2$Dsum, N2$Rsum, g$Ds, g$Rs, k, s, m-h)
ctrp_plot(td$L[2,], td$U[2,], TRUE)
round(td$U[2,],2)
round(td$L[2,],2)


# dx: check (1,2) -> non-rejection
N2$Dsum <- N2$Dsum + g$D[,m-h+1]
N2$Rsum <- N2$Rsum + g$D[,m-h+1]
N2$dfixed <- N2$dfixed + g$D[,m-h+1]
td <- tdp_bounds(N2$Dsum, N2$Rsum, g$Ds, g$Rs, k, s, m-h)
ctrp_plot(td$L[2,], td$U[2,], TRUE)
round(td$U[2,],2)
round(td$L[2,],2)




#BAB - removal of highest in S

new_node2 <- function(dfixed, Ds0, Rs0, Is0, Dsum, Rsum, s0){
  s <- s0 - 1
  Ds <- Ds0[,-s0]
  Is <- matrix(rep(0, B*s), ncol=s)
  Rs <- Is
  e <- Is0[1,1]
  Is[1,] <- Is0[1,-1]
  for(x in (2:B)){
    i <- match(e, Is0[x,])
    Is[x,] <- Is0[x,-i]
    Rs[x,] <- Rs0[x,-i]
  }
  dfixed <- dfixed - Ds0[,s0]
  Dsum <- Dsum - Ds0[,s0]
  Rsum <- Rsum - Ds0[,s0]
  out <- list("dfixed"=dfixed, "Ds"=Ds, "Rs"=Rs, "Is"=Is, "Dsum"=Dsum, "Rsum"=Rsum)
  return(out)
}


# check (2,3) -> rejection
h <- 1
N1 <- new_node2(g$ds, g$Ds, g$Rs, g$Is, g$Dsum, g$Rsum, s)
td <- tdp_bounds(N1$Dsum, N1$Rsum, N1$Ds, N1$Rs, k, s-1, m)
ctrp_plot(td$L[1,], td$U[1,], TRUE)
round(td$U[1,],2)
round(td$L[1,],2)








