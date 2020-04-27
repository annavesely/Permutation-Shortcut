setwd('C:/Users/coshk/Desktop/Perm_project/cpp_code')
myseed <- 33
set.seed(33)


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




# Beta: given the number of variables f,
# the percentage perc of active variables,
# a vector beta_active of non-zero beta values
# and a vector beta_inactive of zeros
# it generates a vector of length f with ceiling(f * perc) active values
generate_beta <- function(f, perc, beta_active, beta_inactive){
  n1 <- ceiling(f * perc) # active
  n2 <- f - n1 # inactive
  
  if(n1==0){
    beta <- beta_inactive[1:n2]
  }else if(n2==0){
    beta <- beta_active[1:n1]
  }else{
    beta <- c(beta_active[1:n1], beta_inactive[1:n2])
  }
  return(beta)
}




# S: given the number of variables f,
# the percentage s_size of variables contained in the subset S,
# and the percentage s_active of active variables within S,
# it generates a vector of indices of length ceiling(f * s_size)
# with ceiling(s * s_active) active indices
generate_S <- function(f, s_size, s_active){
  s <- ceiling(f * s_size) # size
  n1 <- ceiling(s * s_active) # active
  n2 <- s - n1 # inactive
  full <- (1:f)
  if(n1==0){
    S <- full[(f-n2+1):f]
  }else if(n2==0){
    S <- full[1:n1]
  }else{
    S <- c(full[1:n1], full[(f-n2+1):f])
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




# Test for the main body (without BAB)
# Given a subset S, the matrices D, R, I as given by ctrp_set,
# and a significance level alpha, it returns cont=T if if the BAB is needed
# It returns error if R and Rcpp give different outputs
ctrp_test0 <- function(S, D, R, I, alpha){
  f <- ncol(D)
  B <- nrow(D)
  s <- length(S)
  m <- f-s
  k <- ceiling((1-alpha)*B)
  
  # if S=F:
  if(m==0){
    # lower and upper bounds (equal)
    low_R <- R_Q(rowSums(D), k, B)
    low <- Q(rowSums(D), k, B)
    if(low_R != low) stop("R and C give different outputs")
    out <- list("cont"=F, "non_rej"=!low, "BAB"=0)
    return(out)
  }
  
  g <- gen_sub(S, D, R, I, f, m, B, s)
  ind <- (0:m)
  cb_R <- R_compute_bounds(ind, g$R, g$Dsum, g$Rsum, k, m, B, both=T)
  cb = compute_bounds(ind, g$Dsum, g$Rsum, g$R, g$I, m, k, B, TRUE)
  if(cb_R$non_rej != cb$non_rej) stop("Error1")
  if(!cb$non_rej & !all(cb_R$ind == cb$ind)) stop("Error2")
  
  # if there is a non-rejection or there are no indecisives
  if(cb$non_rej | length(cb$ind)==0){
    out <- list("cont"=F, "non_rej"=cb$non_rej, "BAB"=0)
    return(out)
  }
  out <- list("cont"=T, "non_rej"=cb$non_rej, "ind"=cb$ind, "k"=k, "m"=m, "B"=B)
  return(out)
}




# Test for the BAB
# Given a subset S, the matrices D, R, I as given by ctrp_set,
# a significance level alpha and the maximum number of iterations n_max,
# it returns error if R and Rcpp give different outputs or
# use a different number of BAB iterations
# for all four scenarios
check_bab <- function(S, D, R, I, alpha, n_max){
  # 1=lr, 2=lk, 3=hr, 4=hk
  # 1: remove lowest (from_low=T, first_rem=T)
  R1 <- R_ctrp_test(S, D, R, I, alpha, n_max, T, T)
  C1 <- ctrp_test(S, D, R, I, alpha, n_max, T, T)
  if(nr_num(R1$non_rej) != nr_num(C1$non_rej) | R1$BAB != C1$BAB) stop("Error BAB1")
  
  # 2: keep lowest (from_low=T, first_rem=F)
  R2 <- R_ctrp_test(S, D, R, I, alpha, n_max, T, F)
  C2 <- ctrp_test(S, D, R, I, alpha, n_max, T, F)
  if(nr_num(R2$non_rej) != nr_num(C2$non_rej) | R2$BAB != C2$BAB) stop("Error BAB2")
  
  # 3: remove highest (from_low=F, first_rem=T)
  R3 <- R_ctrp_test(S, D, R, I, alpha, n_max, F, T)
  C3 <- ctrp_test(S, D, R, I, alpha, n_max, F, T)
  if(nr_num(R3$non_rej) != nr_num(C3$non_rej) | R3$BAB != C3$BAB) stop("Error BAB3")
  
  
  # 4: keep highest (from_low=F, first_rem=F)
  R4 <- R_ctrp_test(S, D, R, I, alpha, n_max, F, F)
  C4 <- ctrp_test(S, D, R, I, alpha, n_max, F, F)
  if(nr_num(R4$non_rej) != nr_num(C4$non_rej) | R4$BAB != C4$BAB) stop("Error BAB4")
}




# Test
# Given the intervals for the number of variables f (f_int)
# the percentage perc of active variables (perc_int),
# the number of permutations B (B_int),
# the percentage s_size of variables contained in the subset S (s_size_int),
# the percentage s_active of active variables within S (s_active_int),
# the significance level alpha (alpha_int),
# and given the mean m and standard deviation for the active beta values
# and the maximum number of iterations,
# it examines all combinations of values and
# returns error if R and Rcpp give different outputs or
# use a different number of BAB iterations
check_RC <- function(f_int=c(10,50,100), perc_int=c(0, 0.01,0.1,0.2,0.5,0.8),
                     B_int=c(10,50,100), s_size_int=c(0.01,0.1,0.2,0.5,0.8, 1),
                     s_active_int=perc_int, alpha_int=c(0.05,0.2),
                     m=10, sd=5, n_max=10000){
  
  set.seed(myseed)
  fmax <- max(f_int)
  beta_active <- rnorm(fmax, m, sd)
  beta_inactive <- rep(0, fmax)
  
  for(f in f_int){
    
    for(perc in perc_int){
      beta <- generate_beta(f, perc, beta_active, beta_inactive)
      
      for(B in B_int){
        G <- gt(f, f, B, 0, beta)
        c <- ctrp_set(G)
        
        for(s_size in s_size_int){
          
          for(s_active in s_active_int[s_active_int <= perc/s_size & s_active_int >= 1 - (1 - perc)/s_size]){
            S <- generate_S(f, s_size, s_active)
            
            for(alpha in alpha_int){
              te <- ctrp_test0(S, c$D, c$R, c$I, alpha)
              if(te$cont) check_bab(S, c$D, c$R, c$I, alpha, n_max)
            }
          }
        }
      }
    }
  }
  return(TRUE)
}







# Time
# Given the number of variables f
# the interval for the percentage perc of active variables (perc_int),
# the number of permutations B,
# the interval for the percentage s_size of variables contained in the subset S (s_size_int),
# the interval for the percentage s_active of active variables within S (s_active_int),
# the significance level alpha,
# and given the mean m and standard deviation for the active beta values
# and the maximum number of iterations,
# it computes the times that R and Rcpp take to examine each scenario
# and returns the two mean times (in milliseconds)
get_time <- function(f, perc_int=c(0, 0.01,0.1,0.2,0.5,0.8),
                     B, s_size_int=c(0.01,0.1,0.2,0.5,0.8, 1),
                     s_active_int=perc_int, alpha,
                     m=10, sd=5, n_max=10000){
  
  set.seed(myseed)
  L <- length(perc_int) * length(s_size_int) * length(s_active_int)
  M <- matrix(rep(NA, 2*L), ncol=2)
  colnames(M) <- c("tR", "tC")
  
  beta_active <- rnorm(f, m, sd)
  beta_inactive <- rep(0, f)
  i <- 0
  
  for(perc in perc_int){
    beta <- generate_beta(f, perc, beta_active, beta_inactive)
    G <- gt(f, f, B, 0, beta)
    c <- ctrp_set(G)
    
    for(s_size in s_size_int){
      
      for(s_active in s_active_int[s_active_int <= perc/s_size & s_active_int >= 1 - (1 - perc)/s_size]){
        S <- generate_S(f, s_size, s_active)
        #tR <- system.time(R_ctrp_test(S, c$D, c$R, c$I, alpha, n_max, F, T))[3]
        tC <- system.time(ctrp_test(S, c$D, c$R, c$I, alpha, n_max, F, T))[3]
        i <- i+1
        #M[i,] <- c(tR, tC)
        M[i,] <- c(0, tC)
      }
    }
  }
  M <- M[(1:i),]
  out <- list("R_mean"=round(mean(M[,1]),2), "C_mean"=round(mean(M[,2]),2))
  return(out)
}




# Tests:
perc_int <- c(0,0.01,0.1,0.2,0.5,0.8)
B_int <- c(10,50,100)
s_active_int <- c(0,0.01,0.1,0.2,0.5,0.8)
s_size_int <- c(0.01,0.1,0.2,0.5,0.8)
alpha_int <- c(0.05,0.2)
m=10
sd=5
n_max=10000

check_RC(f_int=c(10), perc_int, B_int, s_size_int, s_active_int, alpha_int, m, sd, n_max)
check_RC(f_int=c(50), perc_int, B_int, s_size_int, s_active_int, alpha_int, m, sd, n_max)
check_RC(f_int=c(100), perc_int, B_int, s_size_int, s_active_int, alpha_int, m, sd, n_max)


# Time:
alpha <- 0.05
get_time(f=100, perc_int, B=100, s_size_int, s_active_int, alpha, m, sd, n_max)
get_time(f=100, perc_int, B=500, s_size_int, s_active_int, alpha, m, sd, n_max)
get_time(f=500, perc_int, B=100, s_size_int, s_active_int, alpha, m, sd, n_max)
get_time(f=500, perc_int, B=500, s_size_int, s_active_int, alpha, m, sd, n_max)






s_active <- 1
perc <- 0.01
s_size <- 0.01
f <- 100
B <- 100
alpha <- 0.2
m <- 10
sd <- 5
n_max <- 10000

beta_active <- rnorm(f, m, sd)
beta_inactive <- rep(0, f)
beta <- generate_beta(f, perc, beta_active, beta_inactive)
G <- gt(f, f, B, 0, beta)
c <- ctrp_set(G)
S <- generate_S(f, s_size, s_active)


Rprof()
prova2 <- ctrp_test(S, c$D, c$R, c$I, alpha, n_max, F, T)
Rprof(NULL)
summaryRprof()















