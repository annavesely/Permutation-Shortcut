
### CASE 2: KEEPING OF THE LOWEST STATISTIC ###


# Internal function - Branch and Bound
# Used when keeping the lowest statistic
# Given D_0 (initial matrix) and
# d_kept, R, I, Dsum, Rsum at the previous step with level lev-1
# (where the level is the number of indices that have already been considered)
# the function considers the index of the lev-th lowest statistic
# and removes the corresponding elements from R, I,
# and adds it to d_kept, Dsum and Rsum

us1_2 <- function(lev, d_kept, D_0, R, I, Dsum, Rsum, m=ncol(D_0), B=nrow(D_0)){
  
  # number of columns of I_new
  #(ncol(Rsum_new)=ncol(I)=N+1, ncol(Rsum_prev)=N+2)
  N <- m-lev
  
  # index to be considered
  r <- I[1,N+1]
  
  d_kept <- d_kept + D_0[,lev]
  I_new <- matrix(rep(NA, B*N), ncol=N)
  R_new <- I_new
  Rsum_new <- matrix(rep(NA, B*(N+1)), ncol=N+1)
  Dsum_new <- Dsum[,2:(N+2)]
  
  for(x in seq(B)){
    i <- match(r, I[x,])
    I_new[x,] <- I[x,-i]
    R_new[x,] <- R[x,-i]
    Rsum_new[x,] <- Rsum[x,-i]
    if(i>1){
      Rsum_new[x, 1:(i-1)] <- Rsum_new[x, 1:(i-1)] + D_0[x,lev]
    }
  }
  
  out <- list("R"=R_new, "I"=I_new, "Dsum"=Dsum_new, "Rsum"=Rsum_new, "d_kept"=d_kept)
  return(out)
}




# Internal function - Branch and Bound
# Used when removing the lowest statistic at the same level of the previous step
# Given D_0 (initial matrix) and
# d_kept, Dsum, Rsum at the previous step with the same level lev
# (where the level is the number of indices that have already been considered)
# the function considers the index of the lev-th highest statistic
# and removes the corresponding elements from Dsum, Rsum and d_kept
# (R and I do not vary)

us2_2 <- function(lev, d_kept, D_0, R, I, Dsum, Rsum, m=ncol(D_0)){
  
  #index to be considered
  N <- m-lev
  A <- D_0[,lev]
  
  d_kept <- d_kept - A
  Dsum_new <- Dsum - A
  Rsum_new <- Rsum - A
  
  out <- list("R"=R, "I"=I, "Dsum"=Dsum_new, "Rsum"=Rsum_new, "d_kept"=d_kept)
  return(out)
}




# Internal function - Branch and Bound
# Used when removing the lowest statistics at a lower level than the previous step
# Given I_0, R_0 and Dsum_0 (initial matrices)
# the test statistic d_kept at level lev
# (where the level is the number of indices that have already been considered)
# the function considers the indices of the lev lowest statistics,
# removes the corresponding elements from d_kept and Dsum, and
# removes the corresponding elements from I and R
# and finally computes the cumulative sums Rsum

us3_2 <- function(lev, d_kept, D_0, R_0, I_0, Dsum_0, m=ncol(R_0), B=nrow(R_0)){
  
  # new index to be considered
  N <- m-lev
  
  Dsum_new <- Dsum_0[,(lev+1):(m+1)] - Dsum_0[,lev+1] + Dsum_0[,1] + d_kept
  
  #vector of the indices to be considered
  r <- I_0[1,(N+1):m]
  
  I_new <- matrix(rep(NA, B*N), ncol=N)
  R_new <- I_new
  
  for(x in seq(B)){
    i <- match(r, I_0[x,])
    I_new[x,] <- I_0[x,-i]
    R_new[x,] <- R_0[x,-i]
  }
  
  Rsum_new <- t(apply(cbind(Dsum_0[,1] + d_kept, R_new), 1, cumsum))

  out <- list("R"=R_new, "I"=I_new, "Dsum"=Dsum_new, "Rsum"=Rsum_new, "d_kept"=d_kept)
  return(out)
}



# Internal function - Branch and Bound
# Given the initial matrices D_0, R_0, I_0, Dsum_0 and Rsum_0
# as defined in ctrp_test
# the function partitions the total space according to the highest observed statistic
# Firstly, it removes as many indices as it can until a certain rejection is found
# Then it explores the node where the last (lev-th) statistic is kept instead than removed
# If it is certainly rejected, it proceeds by adding another statistic ((lev-1)-th)
# It continues by removing other indices, iterating the previous steps.
# The function returns non_rej (TRUE if S is not rejected, and NULL if the algorithm
# makes n_max steps without a decisive outcome)
# and BAB, the number of steps made.

ctrp_bab2 <- function(indecisive_0, D_0, R_0, I_0, Dsum_0, Rsum_0,
                     k=ceiling(0.95*nrow(R)), m=ncol(D_0), B=nrow(D_0), n_max=10000){
  
  list_lev <- list(0)
  list_kept <- list(0)
  list_ind <- list(indecisive_0)
  
  BAB <- 0 # steps
  lev <- 0
  non_rej <- F
  cond <- T
  
  us <- list("R"=R_0, "I"=I_0, "Dsum"=Dsum_0, "Rsum"=Rsum_0, "d_kept"=0)
  
  while(length(list_lev)>0 & BAB<n_max){
    
    # we keep keeping the lowest statistic until we can close a branch
    while(cond & BAB<n_max){
      BAB <- BAB + 1
      lev <- lev + 1
      us <- us1_2(lev, us$d_kept, D_0, us$R, us$I, us$Dsum, us$Rsum, m, B)
      
      indecisive <- compute_bounds_fast(tail(list_ind,1)[[1]]-1, us$R, us$Dsum, us$Rsum, k, m-lev)
      cond <- length(indecisive)>0
      
      if(cond){
        list_lev <- append(list_lev, lev)
        list_kept <- append(list_kept, list(us$d_kept))
        list_ind <- append(list_ind, list(indecisive))
      }
    }
    
    # then we explore the branch right next to it
    BAB <- BAB + 1
    us <- us2_2(lev, us$d_kept, D_0, us$R, us$I, us$Dsum, us$Rsum, m)
    cb <- compute_bounds(tail(list_ind,1)[[1]], us$R, us$Dsum, us$Rsum, k, m-lev)
    if(cb$non_rej){
      out <- list("non_rej"=T, "BAB"=BAB)
      return(out)
    }
    
    # remove previous element from list, since now that the right branch has been
    # generated, it has been completely analyzed
    L <- length(list_lev)
    list_lev <- list_lev[-L]
    list_kept <- list_kept[-L]
    list_ind <- list_ind[-L]
    L <- L-1
    
    cond <- length(cb$indecisive)>0
    
    # if the set is not indecisive, we take the next element from the list, and remove
    # the lev-th statistic
    while(!cond & L>0 & BAB<n_max){
      BAB <- BAB + 1
      lev <- tail(list_lev,1)[[1]] + 1
      us <- us3_2(lev, tail(list_kept,1)[[1]], D_0, R_0, I_0, Dsum_0, m, B)
      cb <- compute_bounds(tail(list_ind,1)[[1]], us$R, us$Dsum, us$Rsum, k, m-lev)
      
      if(cb$non_rej){
        out <- list("non_rej"=T, "BAB"=BAB)
        return(out)
      }
      
      # remove previous element from list, since now that the right branch has been
      # generated, it has been completely analyzed
      L <- length(list_lev)
      list_lev <- list_lev[-L]
      list_kept <- list_kept[-L]
      list_ind <- list_ind[-L]
      L <- L-1
      
      cond <- length(cb$indecisive)>0
    }
    
    # if one indecisive set has been found, it is added
    # (otherwise it means that the list has become empty without finding such a set)
    if(cond){
      list_lev <- append(list_lev, lev)
      list_kept <- append(list_kept, list(us$d_kept))
      list_ind <- append(list_ind, list(cb$indecisive))
    }
  }
  
  # if there are some non-explored nodes (i.e. the algorithm stopped because BAB>n_max)
  if(length(list_lev)>0){
    out <- list("non_rej"=NULL, "BAB"=BAB)
  }else{
    out <- list("non_rej"=F, "BAB"=BAB)
  }
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

ctrp_test2 <- function(S, D, R, I, alpha=0.05, n_max=10000){
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
  
  out <- ctrp_bab2(cb$indecisive, g$D, g$R, g$I, Dsum, Rsum, k, m, B, n_max)
  return(out)
}






