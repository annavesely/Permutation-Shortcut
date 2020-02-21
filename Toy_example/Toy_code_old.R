
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




# Given a matrix G of test statistics, the function returns
# M, matrix of centered test statistics where the first row is sorted in descending order
# D, matrix D where each row is sorted in descending order
# I, matrix of the indices corresponding to the elements of D

ctrp_set <- function(G){
  f <- ncol(G)
  B <- nrow(G)
  
  # ordering according to the first row
  Im <- order(G[1,], decreasing=T)
  M <- G[,Im]
  
  # centered test statistics
  M <- sweep(M, 2, M[1,])
  
  # centered test statistics where each row is in ascending order
  o <- t(apply(M, 1, order, decreasing=T))
  D <- t(sapply(seq(B), function(x) M[x, o[x,]]))
  I <- matrix(rep(Im,B), ncol=f, byrow=T)
  I <- t(sapply(seq(B), function(x) I[x, o[x,]]))
  
  out <- list("M"=M, "D"=D ,"I"=I)
  return(out)
}




# Given a vector X and a significance level, the function returns TRUE
# if 0 is greater than the quantile
aq <- function(X, alpha=0.05){
  k <- ceiling((1-alpha)*length(X))
  Xord <- sort(X, na.last = NA, decreasing=F, method="quick")
  c <- Xord[k]
  #out <- sign(c)<0
  #return(out)
  return(c)
}




# Internal function, used in ctrp_set
# Given a set S of indices, it splits M, D and I, returning:
# Ms, vector of the test statistic for S
# Mc, matrix M with only the indices not in S
# Dc, matrix D with only the indices not in S

gen_sub <- function(S, M, D, I){
  B <- nrow(M)
  s <- length(S)
  
  #matrix/vector of indices of S in I
  ind <- t(sapply(seq(B), function(x) match(S, I[x,])))
  
   # if s=1, then ind is a vector, otherwise it is a matrix
  if(s==1){
    Ms <- M[,ind[1]]
    Mc <- M[,-ind[1]]
    Dc <- t(sapply(seq(B), function(x) D[x, -ind[x]]))
    Ic <- t(sapply(seq(B), function(x) I[x, -ind[x]]))
  }else{
    Ms <- rowSums(M[,ind[1,]])
    Mc <- M[,-ind[1,]]
    Dc <- t(sapply(seq(B), function(x) D[x, -ind[x,]]))
    Ic <- t(sapply(seq(B), function(x) I[x, -ind[x,]]))
  }
  
  out <- list("Ms"=Ms, "Mc"=Mc, "Dc"=Dc, "Ic"=Ic)
  return(out)
}









# Given a vector of indices S, the function returns the vectors of the lower
# and upper bounds.

ctrp_test <- function(S, M, D, I, alpha=0.05){
  f <- ncol(M)
  B <- nrow(M)
  s <- length(S)
  m <- f-s
  
  # if s=f, we are testing the whole set of indices
  if(s==f){
    # lower and upper bounds (equal)
    low <- aq(rowSums(M), alpha)
    up <- low
  }else{
    # lower and upper bounds for each possible superset size
    low <- rep(NA, m+1)
    up <- low
    
    g <- gen_sub(S, M, D, I)

    # v=0 additional columns (only the set S)
    L <- g$Ms
    U <- L
    low[1] <- aq(L, alpha)
    up[1] <- low[1]
    
    # v=1 additional column
    v <- 1
    
    # if m=1, then M is a vector
    # and the only possible additional size is v=1
    if(m==1){
      L <- L + g$Mc
      U <- U + g$Dc
      low[2] <- aq(L, alpha)
      up[2] <- aq(U, alpha)
    }else{
      while(v<m+1){
        L <- L + g$Mc[,m-v+1] # L + the v-th smallest
        U <- U + g$Dc[,v] # U + the v-th highest
        low[v+1] <- aq(L, alpha)
        up[v+1] <- aq(U, alpha)
        v <- v+1
      }
    }
  }
  
  out <- list("low"=low, "up"=up, "Ms"=g$Ms, "Mc"=g$Mc, "Dc"=g$Dc, "Ic"=g$Ic)
  return(out)
}




# Given the vectors of the lower and upper bounds
# the function plots the lower bound in blue,
# the upper bound in red,
# and the value zero in black

require(ggplot2)

ctrp_plot <- function(low, up){
  df <- data.frame(size=seq(length(low))-1, low=low, up=up, obs=0)
  ggplot(df, aes()) +
    geom_point(aes(size, low), colour="blue", size=2) +
    geom_line(aes(size, low), colour="blue") +
    geom_point(aes(size, up), colour="red", size=2) +
    geom_line(aes(size, up), colour="red") +
    geom_point(aes(size, obs), size=2) +
    ylab("") + xlab("v")
}





# -------------------------------------------------------------------------------- #

# TOY EXAMPLE

# matrix of global test statistics
initial_data <- gt(n=20, f=5, B=10, beta0=0, beta=c(20,10,5,0,0))
c <- ctrp_set(initial_data)


# set under testing
S <- c(3)
te <- ctrp_test(S, c$M, c$D, c$I, 0.20)
ctrp_plot(te$low, te$up)

te$up
te$low

# -------------------------------------------------------------------------------- #
alpha <- 0.2
B <- 10
m <- 4
Ms <- te$Ms
Mc <- t(sapply(seq(B), function(x) te$Mc[x,rev(seq(4))]))
Dc <- te$Dc
Ic <- te$Ic
indecisive <- c(1,2)

p <- prova(from_top=F, indecisive, Mc, Dc, Ic, B, alpha)
# rem: 1, keep:2

pr <- prova(from_top=T, p$indecisive_r, p$Mc, p$Dc, p$Ic, B, alpha)
# STOP

prova <- function(from_top, indecisive, Mc, Dc, Ic, B, alpha){
  if(from_top){
    i <- 1
  }else{
    i <- ncol(Mc)
  }
  Q <- Ic[1,i]
  ind <- t(sapply(seq(B), function(x) match(Q, Ic[x,])))
  ind_m <- ncol(Mc)+1-ind[1]
  
  Mc2 <- Mc[,-ind_m]
  Dc2 <- t(sapply(seq(B), function(x) Dc[x, -ind[x]]))
  Ic2 <- t(sapply(seq(B), function(x) Ic[x, -ind[x]]))
  
  if(nrow(Dc2)==1){
    Dsum2 <- t(apply(cbind(Ms,t(Dc2)), 1, cumsum))
    Msum2 <- t(apply(cbind(Ms,t(Mc2)), 1, cumsum))
  }else{
    Dsum2 <- t(apply(cbind(Ms,Dc2), 1, cumsum))
    Msum2 <- t(apply(cbind(Ms,Mc2), 1, cumsum))
  }

  # lower bounds
  
  n_ind <- length(indecisive)
  low_r <- rep(NA, n_ind)
  low_k <- low_r
  up_r <- low_r
  up_k <- low_r
  
  h <- 1
  cond_low <- T
  while(cond_low & h<=length(indecisive)){
    v <- indecisive[h]
    low_r[h] <- quantile(Msum2[,v+1], 1-alpha, names=F)
    low_k[h] <- quantile(Msum2[,v]+Mc[,ind_m], 1-alpha, names=F)
    cond_low <- (sign(low_r[h])==-1) & (sign(low_k[h])==-1)
    h <- h+1
  }
  
  if(!cond_low){
    out <- list("non_rej"=T, "indecisive_r"=NULL, "indecisive_k"=NULL,
                "Mc"=NULL, "Dc"=NULL, "Ic"=NULL)
    return(out)
  }
  
  up_r <- sapply((1:n_ind), function(x)
    quantile(Dsum2[,indecisive[x]+1], 1-alpha, names=F))
  ind_r <- sapply((1:n_ind), function(x)
    sign(up_r[x]>-1))
  
  if(all(!ind_r)){
    indecisive_r <- NULL
  }else{
    indecisive_r <- indecisive[ind_r]
  }
  
  up_k <- sapply((1:n_ind), function(x)
    quantile(Dsum2[,indecisive[x]]+Mc[,ind_m], 1-alpha, names=F))
  ind_k <- sapply((1:n_ind), function(x)
    sign(up_k[x]>-1))
  
  if(all(!ind_k)){
    indecisive_k <- NULL
  }else{
    indecisive_k <- indecisive[ind_k]
  }
  out <- list("non_rej"=F, "indecisive_r"=indecisive_r, "indecisive_k"=indecisive_k,
              "low_r"=low_r, "up_r"=up_r, "low_k"=low_k, "up_k"=up_k,
              "Mc"=Mc2, "Dc"=Dc2, "Ic"=Ic2)
  return(out)
}




# -------------------------------------------------------------------------------- #


