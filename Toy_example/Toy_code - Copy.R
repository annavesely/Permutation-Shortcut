
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




# Given a matrix of test statistics, the function returns
# M, matrix of centered test statistics where the first row is sorted in descending order
# D, matrix D where each row is sorted in descending order
# I, matrix of the indices corresponding to the elements of D

ctrp_set <- function(M){
  f <- ncol(M)
  B <- nrow(M)
  
  # ordering according to the first row
  Im <- order(M[1,], decreasing=T)
  M <- M[,Im]
  
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
  c <- quantile(X, 1-alpha, names=F)
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
    ylab("") + xlab("additional size")
}





# -------------------------------------------------------------------------------- #

# TOY EXAMPLE

# matrix of global test statistics
M <- gt(n=20, f=5, B=10, beta0=0, beta=c(20,10,5,0,0))
c <- ctrp_set(M)


# set under testing
S <- c(3)
te <- ctrp_test(S, c$M, c$D, c$I, 0.20)
ctrp_plot(te$low, te$up)




# -------------------------------------------------------------------------------- #


# somme cumulate
#t(apply(M, 1, cumsum))
#ma devo invertire l'ordine
# perché sommo le colonne dal fondo
# Ma ha senso salvare un'altra matrice
# o è meno peggio rifare le somme tutte le volte?

Mnew <- cbind(te$Ms, t(apply(te$Mc,1,rev)))
Dnew <- cbind(te$Ms, te$Dc)

Msum <- t(apply(Mnew, 1, cumsum))
Dsum <- t(apply(Dnew, 1, cumsum))
round(Msum,2)
round(Dsum,2)





Mnew <- t(apply(te$Mc,1,rev))
Dnew <- te$Dc

Msum <- t(apply(Mnew, 1, cumsum))
Dsum <- t(apply(Dnew, 1, cumsum))
round(Msum,2)
round(Dsum,2)

nq <- function(v, sum_matrix, Ms, alpha=0.05){
  out <- quantile(Ms + sum_matrix[,v], 1-alpha, names=F)
  return(out)
}

sapply(seq(4), nq, Msum, te$Ms, 0.20)
sapply(seq(4), nq, Dsum, te$Ms, 0.20)


round(Dsum-Msum,2)


aq(Mnew[,1], 0.8)
round(apply(Msum, 2, aq, alpha=0.8),2)
round(apply(Dsum, 2, aq, alpha=0.8),2)




xprova <- runif(10, min=-10, max=10)
yprova <- runif(10, min=0, max=10)
zprova <- xprova-yprova
qx <- quantile(xprova, 0.80, names=F)
qz <- quantile(zprova, 0.80, names=F)
out <- qz<=qx
out


# -------------------------------------------------------------------------------- #

# sizes v=1 and v=2 are indecisive





newsub <- gen_sub(1, te$Mc, te$Dc, te$I)
Dnew <- newsub$Dc
Inew <- newsub$Ic
Mfixed <- newsub$Ms

# REMOVE 1
# low = same as before
# up = Ms + sum of the v highest elements of Dnew

low_rem <- c(NA, NA)
up_rem <- c(NA, NA)

# v=1
L <- te$Ms + te$Mc[,4] # same as before
low_rem[1] <- aq(L, 0.20)

U <- te$Ms + Dnew[,1]
up_rem[1] <- aq(U, 0.20)

# v=2
L <- L +te$Mc[,3] # same as before
low_rem[2] <- aq(L, 0.20)

U <- U + Dnew[,2]
up_rem[2] <- aq(U, 0.20)


# KEEP ONE
# low = Ms + Mfixed + sum of last (v-1) elements of M
# up = Ms + Mfixed + sum of the (v-1) highest elements of Dew

low_keep <- c(NA, NA)
up_keep <- c(NA, NA)

# v=1
L <- te$Ms + te$Mc[,1] # U = L
low_keep[1] <- aq(L, 0.20)


