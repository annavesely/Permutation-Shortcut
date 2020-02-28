# ERROR CHECK



# Internal function - Branch and Bound
# Used when keeping the highest statistics at a lower level than the previous step
# Given I_0, R_0 and Dsum_0 (initial matrices)
# the test statistic d_kept at level lev
# (where the level is the number of indices that have already been considered)
# the function considers the indices of the lev highest statistics,
# adds the corresponding elements to d_kept and Dsum, and
# removes the corresponding elements from I and R
# and finally computes the cumulative sums Rsum

us_check <- function(lev, d_kept, rem=T, D_0, R_0, I_0, m=ncol(R_0), B=nrow(R_0)){
  
  # new index to be considered
  N <- m-lev
  
  if(!rem){
    A <- D_0[,N+1]
    d_kept <- d_kept + A
  }
  
  
  #vector of the indices to be considered
  r <- I_0[1,1:lev]
  
  I_new <- matrix(rep(NA, B*N), ncol=N)
  R_new <- I_new
  D_new <- I_new
  
  for(x in seq(B)){
    i <- match(r, I_0[x,])
    iD <- N-i+2
    I_new[x,] <- I_0[x,-i]
    R_new[x,] <- R_0[x,-i]
    D_new[x,] <- D_0[x,-iD]
  }
  
  first_col <- Dsum_0[,1] + d_kept
  
  Rsum_new <- t(apply(first_col, R_new), 1, cumsum)
  
  out <- list("R"=R_new, "Dsum"=Dsum_new, "Rsum"=Rsum_new, "d_kept"=d_kept)
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

ctrp_bab_check <- function(indecisive_0, D_0, R_0, I_0, Dsum_0, Rsum_0,
                      k=ceiling(0.95*nrow(R)), m=ncol(D_0), B=nrow(D_0), n_max=10000){
  
  list_lev <- list(0)
  list_kept <- list(0)
  list_ind <- list(indecisive_0)
  
  BAB <- 0 # steps
  lev <- 0
  non_rej <- F
  cond <- T
  
  us <- list("R"=R_0, "Dsum"=Dsum_0, "Rsum"=Rsum_0, "d_kept"=0)
  
  while(length(list_lev)>0 & BAB<n_max){
    
    # we keep removing the highest statistic until we can close a branch
    while(cond & BAB<n_max){
      BAB <- BAB + 1
      lev <- lev + 1
      us <- us_check(lev, us$d_kept, rem=T, D_0, R_0, I_0, m, B)
      indecisive <- compute_bounds_fast(tail(list_ind,1)[[1]], us$R, us$Dsum, us$Rsum, k, m-lev)
      cond <- length(indecisive)>0
      
      if(cond){
        list_lev <- append(list_lev, lev)
        list_kept <- append(list_kept, list(us$d_kept))
        list_ind <- append(list_ind, list(indecisive))
      }
    }
    
    # then we explore the branch right next to it
    BAB <- BAB + 1
    us <- us_check(lev, us$d_kept, rem=F, D_0, R_0, I_0, m, B)
    cb <- compute_bounds(tail(list_ind,1)[[1]]-1, us$R, us$Dsum, us$Rsum, k, m-lev)
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
    
    # if the set is not indecisive, we take the next element from the list, and add
    # the lev-th statistic
    while(!cond & L>0 & BAB<n_max){
      BAB <- BAB + 1
      lev <- tail(list_lev,1)[[1]] + 1
      us <- us_check(lev, tail(list_kept,1)[[1]], rem=F, D_0, R_0, I_0, m, B)
      cb <- compute_bounds(tail(list_ind,1)[[1]]-1, us$R, us$Dsum, us$Rsum, k, m-lev)
      
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

ctrp_test_check <- function(S, D, R, I, alpha=0.05, n_max=10000){
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
  
  out <- ctrp_bab_check(cb$indecisive, g$D, g$R, g$I, Dsum, Rsum, k, m, B, n_max)
  return(out)
}


