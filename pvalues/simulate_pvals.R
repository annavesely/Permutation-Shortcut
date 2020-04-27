#require(mvtnorm)
#myseed <- 33
#set.seed(myseed)

require(pracma)


# Generate X
# Given the number f of variables, the number f1 of false null hp
# the level rho of equi-correlation between pairs of variables
# the signal-to-noise ratio SNR and the number n of observations,
# the function returns a nxf matrix X where each column contains the n observations
# of a variable. The first ceiling(f*perc) variables have mean SNR/sqrt(m),
# and the remaining have mean 0.
# Notice: rho in (-1/(f-1), 1)

#X_generate <- function(f, f1, rho, SNR, n){
#  set.seed(myseed)
#  Sigma <- (matrix(rep(rho, f*f), nrow=f) + diag(1-rho, f))
#  epsilon <- rmvnorm(n, rep(0,f), Sigma)
#  mu <- c(rep(SNR*f/f1, f1), rep(0, f - f1))
#  #q <- pnorm(SNR*f1/f, lower.tail=FALSE)
#  #mu <- c(rep(SNR*f/f1, f1), rep(0, f - f1))
#  #mu <- c(rep(q, f1), rep(0, f - f1))
#  X <- t(apply(epsilon, 1, function(x) x+mu))
#  return(X)
#}

#z1 <- matrix(rnorm(k), k, n)
#z0 <- matrix(rnorm(n*k), k, n)
#z <- sqrt(1-rho)*z0 + sqrt(rho)*z1



# Generate X
# Given the number f of variables, the number f1 of false null hp
# the levels rho1 and rho2 of equi-correlation between pairs of variables
# within the groups of the false and true null hp
# (assuming 0 correlation between the two groups),
# the signal-to-noise ratio SNR and the number n of observations,
# and a condition decr
# the function returns a nxf matrix X where each column contains the n observations
# The first f1 variables have mean q (where q=SNR*f/f1 if decr=T, q=SNR otherwise)
# and correlation rho1,
# and the remaining have mean 0 and correlation rho2
X_generate <- function(f, f1, rho1, rho2, SNR, n, decr=TRUE){
  
  # if the null hp are all all true or all false
  if(f1==f | f1==0){
    if(f1==f){
      rho <- rho1
      mu <- SNR
    }else{
      rho <- rho2
      mu <- 0
    }
    shared <- matrix(rep(rnorm(n),f), ncol=f) # shared variation
    own <- matrix(rnorm(f*n), ncol=f) # own variation
    X <- mu + sqrt(1-rho)*own + sqrt(rho)*shared
    return(X)
  }
  
  # active variables
  q <- ifelse(decr, SNR*f/f1, SNR)
  shared1 <- matrix(rep(rnorm(n),f1), ncol=f1) # shared variation
  own1 <- matrix(rnorm(f1*n), ncol=f1) # own variation
  X1 <- q + sqrt(1-rho1)*own1 + sqrt(rho1)*shared1
  
  # inactive variables
  f2 <- f-f1
  shared2 <- matrix(rep(rnorm(n),f2), ncol=f2) # shared variation
  own2 <- matrix(rnorm(f2*n), ncol=f2) # own variation
  X2 <- sqrt(1-rho2)*own2 + sqrt(rho2)*shared2
  
  X <- cbind(X1,X2)
  return(X)
}




# Given a vector V of length n, the function returns the p-value for
# the one-sample t-test for the null hp that the mean is zero
p_compute <- function(V, n=length(V)){
  sd <- sqrt(var(V))
  if(sd==0) sd <- .Machine$double.xmin
  tstat <- mean(V)/(sd/sqrt(n)) # t-statistic to test whether the mean is 0
  #p <- 2 * pt(abs(tstat), df=n-1, lower.tail=FALSE)
  # p <- 2 * pt(abs(tstat), df=n-1, lower.tail=FALSE)
  p <- pt(tstat, df=n-1, lower.tail=FALSE)
  return(p)
}




# Given a nxf matrix X, n=nrow(X), f=ncol(X) and the number B of permutations
# the function returns a Bxf matrix G where
# G[i,j] is the p-value for the t-test on the j-th transformation of X[,i]
# The t-test is a one-sample t-test for the null hp that the mean is zero
# The transformation is sign-flipping

p_matrix <- function(X, n=nrow(X), f=ncol(X), B){
  G = matrix(rep(NA, B*f), ncol=f)
  G[1,] = apply(X, 2, p_compute)
  
  # signs
  z <- matrix(sample(c(-1,1), n*(B-1), replace=TRUE), ncol=(B-1))
  Xnew <- matrix(rep(NA, f*n), ncol=f)

  for(b in(2:B)){
    Xnew <- apply(X, 2, function(x) z[,b-1]*x) # sign-flipping of X columns
    G[b,] = apply(Xnew, 2, p_compute) # corresponding vector of p-values
  }
  return(G)
}



# Given an indexr and a matrix A of p-values, the function returns
# log(A) if r=0, and A^r otherwise
r_p_matrix <- function(r, A, B=nrow(A), f=ncol(A)){
  
  if(r>0){
    G <- A^r
  }
  
  else if(r==0){
    G <- log(A)
  }
  
  else{
    a <- min(A)
    if(is.finite(f*a^r)){
      G <- A^r
    }else{
      lambda <- ceiling(nthroot(f/.Machine$double.xmax, -r)/a)
      G <- (lambda*A)^r
    }
  }
  return(G)
}



# -------------------------------------------------------------- #



# Given the number of variables f,
# the size s of the subset S, and the number of false null hp in S
# it generates a vector of indices of length s
# having the first s1 and the last s2 indices of the full model
S_generate <- function(f, s, s1){
  s2 <- s-s1
  if(s2==0){
    S <- (1:s1)
  }else if (s1==0){
    S <- ((f-s2+1):f)
  }else{
    S <- c((1:s1),((f-s2+1):f))
  } 
  return(S)
}




# Conversion of non-rejection to integer:
# NULL -> 0
# TRUE -> 1
# FALSE -> -1
nr_num <- function(nr){
  ifelse(length(nr)==0, 0, ifelse(nr, 1, -1))
}




# Test for the BAB
# Given a subset S, the matrices D, R, I as given by ctrp_set,
# a significance level alpha and the maximum number of iterations n_max,
# it returns non_rej (numeric) and the different number of BAB iterations
# for all four scenarios
check_bab <- function(S, D, R, I, alpha, n_max){
  C1 <- ctrp_test(S, D, R, I, alpha, n_max, T, T) # remove lowest
  C2 <- ctrp_test(S, D, R, I, alpha, n_max, T, F) # keep lowest
  C3 <- ctrp_test(S, D, R, I, alpha, n_max, F, T) # remove highest
  C4 <- ctrp_test(S, D, R, I, alpha, n_max, F, F) # keep highest
  
  # check on the outcomes (error if they are non-null and different)
  v <- c(nr_num(C1$non_rej), nr_num(C2$non_rej), nr_num(C3$non_rej), nr_num(C4$non_rej))
  v <- v[v != 0] # non-null outcomes
  L <- length(v)
  if(L == 0){
    n_tot <- 0
  }else{
    # if there are some rejections AND some non-rejections
    if(length(v[v==1]) != 0 & length(v[v==-1]) != 0) stop("Error BAB")
    n_tot <- v[1]
  }
  out <- list("non_rej"=n_tot, "b1"=C1$BAB, "b2"=C2$BAB, "b3"=C3$BAB, "b4"=C4$BAB)
  return(out)
}



# Given the vector of indices r and R=length(r),
# the number f of variables and the number f1 of false null hp,
# the levels rho1 and rho2 of equi-correlation between active and inactive pairs of variables
# the signal-to-noise ratio SNR, the number n of observations,
# the number B of transformations (sign-flipping),
# the subset S, the significance level alpha
# and the maximum number of BAB iterations,
# the function returns a vector V of length R having with
# V[r]=1 if the index r leads to a rejection, and 0 otherwise
sim <- function(r, R, f, f1, rho1, rho2, SNR, n, B, S, alpha, n_max){
  V <- rep(0,R)
  # generate a new matrix of p-values
  X <- X_generate(f, f1, rho1, rho2, SNR, n)
  G_init <- p_matrix(X, n, f, B)
  
  # analyze the mtrix with all the r values
  for(j in (1:R)){
    G <- r_p_matrix(r[j], G_init)
    if(!all(is.finite(G))) stop("Error: infinite")
    c <- ctrp_set(G, higher=(r[j]<0))
    cb <- check_bab(S, c$D, c$R, c$I, alpha, n_max)
    if(cb$non_rej == -1) V[j] <- 1 # if there is a rejection
  }
  return(V)
}



# Given the vector of indices r, the number f of variables,
# the levels rho1 and rho2 of equi-correlation between pairs of active and inactive variables
# the signal-to-noise ratio SNR, the number n of observations,
# the number B of transformations (sign-flipping),
# the percentage s_size of variables contained in S,
# the percentage s_star of false null hp within S,
# the vector of percentages o_star of false null hp outside S,
# the significance level alpha, the maximum number of BAB iterations
# and the number W of repeated simulations
# the function returns the matrix M having M[i,j] = percentage
# of rejections for o_star[i] and r[j] (over all W simulations)
bab_pvals <- function(r=c(-100,-10,-2,-1,-0.5,-0.1,0,0.1,0.5,1,2,10,100),
                      f=100, rho1=0, rho2=0, SNR=5, n=20, B=100, s_size=0.2,
                      s_star=0.5, o_star=c(0,0.1,0.2,0.5,0.8,1),
                      alpha=0.20, n_max=10000, W){
  
  s <- ceiling(f * s_size) # size of S
  s1 <- ceiling(s * s_star) # false null hp inside S
  S <- S_generate(f, s, s1)
  
  R <- length(r)
  O <- length(o_star)
  M <- matrix(rep(0,R*O), ncol=R) # matrix for non-rejections
  
  for(i in (1:O)){
    o1 <- ceiling((f-s) * o_star[i]) # false null hp outside S
    f1 <- s1 + o1 # total number of false null hp
    
    w <- 0
    while(w<W){
      M[i,] <- M[i,] + sim(r, R, f, f1, rho1, rho2, SNR, n, B, S, alpha, n_max)
      w <- w+1
    }
  }
  return(M*100/W)
}


#rho1: 0, 0.99
# SNR: 0.001 (la situa bella è quando rho1=0.99 e rho2=0)
# 0,0.1,0.5,0.9,1

bab_pvals(r=c(-100,-10,-2,-1,-0.5,-0.1,0,0.1,0.5,1,2,10,100),
          f=50, rho1=0, rho2=0, SNR=0.1, n=20, B=50, s_size=0.2,
          s_star=0, o_star=c(0,0.1,0.5,0.9,1),
          alpha=0.20, n_max=10000, W=100)


bab_pvals(r=c(1),
          f=100, rho1=0, rho2=0, SNR=100, n=20, B=100, s_size=0.2,
          s_star=0, o_star=c(1),
          alpha=0.20, n_max=10000, W=100)

f <- 100
rho1 <- 0
rho2 <- 0
SNR <- 100
n <- 20
B <- 100
s_size <- 0.2
s_star <- 0
o_star <- 1
alpha <- 0.20
n_max <- 10000


# Given the vector of indices r, the number f of variables,
# the levels rho1 and rho2 of equi-correlation between pairs of active and inactive variables
# the signal-to-noise ratio SNR, the number n of observations,
# the number B of transformations (sign-flipping),
# the percentage s_size of variables contained in S,
# the percentage s_star of false null hp within S,
# the vector of percentages o_star of false null hp outside S,
# the significance level alpha and the maximum number of BAB iterations
# the function returns all th cases which needed the BAB, with the
# corresponding outcome (rejection/non-rejection/null) and the number
# of iterations needed by the four scenarios (remove/keep lowest/highest statistic)
bab_iter <- function(r=c(-100,-10,-2,-1,-0.5,-0.1,0,0.1,0.5,1,2,10,100),
                      f=100, rho1=0, rho2=0, SNR=5, n=20, B=100, s_size=0.2,
                      s_star=0.5, o_star=c(0,0.1,0.2,0.5,0.8,1),
                      alpha=0.20, n_max=10000){
  
  s <- ceiling(f * s_size) # size of S
  s1 <- ceiling(s * s_star) # false null hp inside S
  S <- S_generate(f, s, s1)
  
  R <- length(r)
  O <- length(o_star)
  M <- matrix(rep(NA,R*O*11), ncol=11) # matrix for BAB iterations
  colnames(M) <- c("r", "s_star", "o_star", "s_size", "rho1", "rho2", "non_rej", "RL", "KL", "RH", "KH")
  z <- 0
  
  for(i in (1:O)){
    o1 <- ceiling((f-s) * o_star[i]) # false null hp outside S
    f1 <- s1 + o1
    X <- X_generate(f, f1, rho1, rho2, SNR, n)
    G_init <- p_matrix(X, n, f, B)
    
    for(j in (1:R)){
      G <- r_p_matrix(r[j], G_init)
      c <- ctrp_set(G, higher=(r[j]<0))
      cb <- check_bab(S, c$D, c$R, c$I, alpha, n_max)
      
      if(cb$b1>0){
        z <- z+1
        M[z,] <- c(r[j], s_star*100, o_star[i]*100, s_size*100, rho1, rho2, cb$non_rej, cb$b1, cb$b2, cb$b3, cb$b4)
      }
    }
  }
  return(M[(1:z),])
}

# -------------------------------------------------------------- #


# Given a vector X and a value k, it returns the critical value
# i.e. the k-th statistic (when sorted in increasing order)
Qu <- function(X, k){
  Xord <- sort(X, na.last = NA, decreasing=F, method="quick")
  return(Xord[k])
}






# Given a vector of indices S, the matrices D, R, I as given by ctrp_set,
# the significance level alpha, the function returns the vectors of
# the upper and lower bounds for each superset size
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
    up <- low
    out <- list("up"=up, "low"=low)
    return(out)
  }
  
  g <- gen_sub(S, D, R, I, f, m, B, s)
  ind <- (0:m)
  low <- rep(NA,m+1)
  up <- low
  
  for(i in(1:(m+1))){
    low[i] <- Qu(g$Dsum[,i], k)
    up[i] <- Qu(g$Rsum[,i], k)
  }
  
  nr <- ctrp_test(S, D, R, I, alpha, 0, F, F)$non_rej
  
  out <- list("up"=up, "low"=low, "non_rej"=nr)
  return(out)
}



myseed <- 33
set.seed(myseed)

# Given an index r, the number f of variables,
# the levels rho1 and rho2 of equi-correlation between pairs of active and inactive variables
# the signal-to-noise ratio SNR, the number n of observations,
# the number B of transformations (sign-flipping),
# the percentage s_size of variables contained in S,
# the percentage s_star of false null hp within S,
# the percentage o_star of false null hp outside S,
# and the significance level alpha
# the function returns the vectors of the lower and upper bounds
bounds_pvals <- function(r=0, f=100, rho1=0, rho2=0, SNR=5, n=20, B=100, s_size=0.2,
                         s_star=0.5, o_star=0.8, alpha=0.20){
  set.seed(myseed)
  s <- ceiling(f * s_size) # size of S
  s1 <- ceiling(s * s_star) # false null hp inside S
  S <- S_generate(f, s, s1)
  
  o1 <- ceiling((f-s) * o_star) # false null hp outside S
  f1 <- s1 + o1
  X <- X_generate(f, f1, rho1, rho2, SNR, n)
  G_init <- p_matrix(X, n, f, B)
  G <- r_p_matrix(r, G_init)
  c <- ctrp_set(G, higher=(r<0))
  
  out <- ctrp_bounds(S, c$D, c$R, c$I, alpha)
  return(out)
}




# Given the vectors of the lower and upper bounds
# the function plots the lower bound in blue, the upper bound in red,
# and the value zero in black
require(ggplot2)

ctrp_plot <- function(low, up){
  df <- data.frame(size=seq(length(low))-1, low=low, up=up, obs=0)
  ggplot(df, aes()) +
    geom_line(aes(size, obs), linetype = "dashed", size=1) +
    geom_line(aes(size, up), colour="red", size=1) +
    geom_line(aes(size, low), colour="blue", size=1) +
    ylab("") + xlab("additional indices")
}



# repeat for all r values
bp <- bounds_pvals(r=-100, f=50, rho1=0.99, rho2=0, SNR=0.1, n=20, B=50, s_size=0.2, s_star=0.5,
                o_star=0.1, alpha=0.20)
ctrp_plot(bp$low, bp$up)


# r=-10
bp$low <- bp$low/1e+222
bp$up <- bp$up/1e+222
df <- data.frame(size=seq(length(bp$low))-1, low=bp$low, up=bp$up, obs=0)
ggplot(df, aes()) +
  geom_line(aes(size, obs), linetype = "dashed", size=1) +
  geom_line(aes(size, up), colour="red", size=1) +
  geom_line(aes(size, low), colour="blue", size=1) +
  ylab("") + xlab("additional indices") +
  #scale_y_continuous(limits=c(-3.5e+06,1))
  scale_y_continuous(limits=c(-2.1,0), breaks=c(-2,-1,0), labels=c(-2e+222, -1e+222, 0))

# -------------------------------------------------------------- #

