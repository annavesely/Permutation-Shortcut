# Given a matrix G of test statistics, the function returns
# D, matrix of centered test statistics where the rows are sorted so that
# the first one in ascending order
# R, matrix where each row is sorted in descending order
# I, matrix of the indices corresponding to the elements of D

R_ctrp_set <- function(G){
  f <- ncol(G)
  B <- nrow(G)
  
  # ordering according to the first row
  I_incr <- order(G[1,], decreasing=F)
  D <- G[,I_incr]
  
  # centered test statistics
  D <- sweep(D, 2, D[1,])
  
  # centered test statistics where each row is in ascending order
  o <- t(apply(D, 1, order, decreasing=T))
  R <- t(sapply(seq(B), function(x) D[x, o[x,]]))
  
  # matrix of indices in R
  I <- matrix(rep(I_incr,B), ncol=f, byrow=T)
  I <- t(sapply(seq(B), function(x) I[x, o[x,]]))
  I[1,] <- rev(I[1,])
  
  out <- list("D"=D, "R"=R,"I"=I)
  return(out)
}




# Internal function
# Given a vector X and a value k, it returns the critical value
# i.e. the k-th statistic (when sorted in increasing order)

#Q <- function(X, k=ceiling(0.95*length(X))){
  #Xord <- sort(X, na.last = NA, decreasing=F, method="quick")
  #return(Xord[k])
#}




# Internal function
# Given a vector X and a value k, it computes the lower critical value L
# and returns TRUE if L<0 (no non-rejection found)

#Lcond <- function(X, k=ceiling(0.95*length(X))){
  #L <- Q(X, k)
  #c <- sign(L)==-1
  #return(c)
#}



# Internal function
# Given a vector X and a value k, it computes the upper critical value U
# and returns TRUE if U>=0 (no certain rejection found)

#Ucond <- function(X, k=ceiling(0.95*length(X))){
  #U <- Q(X, k)
  #c <- sign(U)>-1
  #return(c)
#}




# Internal function
# Given a vector X, a value k, and B=length(X),
# it returns TRUE if the critical value
# (k-th statistic when sorted in increasing order)
# is negative
# Notice that L<0 when no non-rejection has been found
# and U<0 when a certain rejection has been found

R_Q <- function(X, k=ceiling(0.95*length(X)), B=length(X)){
  
  t <- B-k+1 # threshold: min. number of non-negative el. to have q>= 0
  n <- 0 # number of non-negative el.
  
  for(i in seq(B)){
    if(sign(X[i]) > -1){
      n <- n+1
      if(n >= t){ # if n reaches the threshold, then q>=0
        return(F)
      }
    }else{
      if(i-n >= k){ # if the negative el. found are already k or more, then q<0
        return(T)
      }
    }
  }
}




# Internal function
# Given a set S of indices, and the matrices D, R and I given by R_ctrp_set,
# it splits D, R and I, returning:
# ds, vector of the test statistic for S
# matrix D with only the indices not in S
# matrix R with only the indices not in S
# matrices Dsum and Rsum of the cumulative sums of ds with D and R

R_gen_sub <- function(S, D, R, I, f=ncol(D), m=ncol(D)-length(S), B=nrow(D), s=length(S)){
  
  if(m==0){
    ds <- rowSums(D)
    out <- list("ds"=ds, "D"=NULL, "R"=NULL, "I"=NULL)
    return(out)
  }
  
  i <- match(S,I[1,])
  i_D <- f+1-i
  Dc <- D[,-i_D]
  
  if(s==1){
    ds <- D[,i_D]
  }else{
    ds <- rowSums(D[,i_D])
  }
  
  # if m=1, then Rc=Dc and Ic are vectors
  if(m==1){
    Ic <- rep(I[1,-i], B)
    Rc <- Dc
  }else{
    Ic <- matrix(rep(NA, B*m), ncol=m)
    Rc <- Ic
    Ic[1,] <- I[1,-i]
    Rc[1,] <- R[1,-i]
    for(x in (2:B)){
      i <- match(S, I[x,])
      Ic[x,] <- I[x,-i]
      Rc[x,] <- R[x,-i]
    }
  }
  Dsum <- t(apply(cbind(ds, Dc), 1, cumsum))
  Rsum <- t(apply(cbind(ds, Rc), 1, cumsum))
  out <- list("ds"=ds, "D"=Dc, "R"=Rc, "I"=Ic, "Dsum"=Dsum, "Rsum"=Rsum)
  return(out)
}




# Internal function.
# Given the vector ind of indecisive sizes,
# the matrices R, Dsum and Rsum,
# k and m=ncol(R) and B=nrow(R),
# it checks both bounds. It returns non_rej (T if there is a non-rejection)
# and the new vector of indecisive sizes.
# In particular, it finds the first column R[,j] such that all the el. are negative/null
# It checks the bounds up to j-1, and then the following bounds
# until one upper bound becomes negative.

R_bounds_both <- function(ind, R, Dsum, Rsum, k=ceiling(0.95*nrow(R)), m=ncol(R), B=nrow(R)){
  
  H <- length(ind)
  up <- rep(F, H) # keeps track of indecisive sizes
  h <- 1
  v <- ind[1]
  low <- R_Q(Dsum[,v+1], k, B)
  if(!low){
    out <- list("non_rej"=T, "ind"=ind)
    return(out)
  }
  
  # if v=0, then up=F
  if(v > 0){
    up[1] <- !R_Q(Rsum[,v+1], k, B)
  }
  
  # first column in R with no positive element
  # (if m=1, then R is a vector)
  if(m==1){
    j <- ifelse(sign(max(R))<1, 1, 2)
  }else{
    j <- Find(function(x) sign(max(R[,x]))<1 , seq(m))
    # if there is no such column
    if(length(j)==0){
      j <- m+1
    }
  }
  
  # bounds for v < j (before the upper c.v. decreases)
  while(low & v<j-1 & h<H){
    h <- h+1
    v <- ind[h]
    low <- R_Q(Dsum[,v+1], k, B)
    up[h] <- !R_Q(Rsum[,v+1], k, B)
  }
  
  # bounds for v >= j (until the upper c.v. becomes negative)
  while(low & up[h] & h<H){
    h <- h+1
    v <- ind[h]
    low <- R_Q(Dsum[,v+1], k, B)
    up[h] <- !R_Q(Rsum[,v+1], k, B)
  }
  out <- list("non_rej"=!low, "ind"=ind[up])
  return(out)
}





# Internal function.
# Given the vector ind of indecisive sizes,
# the matrices R, Dsum and Rsum,
# k and m=ncol(R) and B=nrow(R),
# it checks the upper bound. It returns non_rej=F
# and the new vector of indecisive sizes.
# In particular, it finds the first column R[,j] such that all the el. are negative/null
# It checks the bounds up to j-1, and then the following bounds
# until one upper bound becomes negative.

R_bounds_upper <- function(ind, R, Dsum, Rsum, k=ceiling(0.95*nrow(R)), m=ncol(R), B=nrow(R)){
  
  H <- length(ind)
  up <- rep(F, H) # keeps track of indecisive sizes
  h <- 1
  v <- ind[1]
  up[1] <- !R_Q(Rsum[,v+1], k, B)
  
  # first column in R with no positive element
  # (if m=1, then R is a vector)
  if(m==1){
    j <- ifelse(sign(max(R))<1, 1, 2)
  }else{
    j <- Find(function(x) sign(max(R[,x]))<1 , seq(m))
    # if there is no such column
    if(length(j)==0){
      j <- m+1
    }
  }
  
  # upper bounds for v < j (before the upper c.v. decreases)
  while(v<j-1 & h<H){
    h <- h+1
    v <- ind[h]
    up[h] <- !R_Q(Rsum[,v+1], k, B)
  }
  
  # bounds for v >= j (until the upper c.v. becomes negative)
  while(up[h] & h<H){
    h <- h+1
    v <- ind[h]
    up[h] <- !R_Q(Rsum[,v+1], k, B)
  }
  out <- list("non_rej"=F, "ind"=ind[up])
  return(out)
}




# Internal function.
# Given the vector ind of indecisive sizes,
# the matrices R, Dsum and Rsum,
# k and m=ncol(R) and B=nrow(R)
# and the condition both (T if both bounds must be checked, F if ony the upper bound)
# it checks the bounds for the indecisive sizes.

R_compute_bounds <- function(ind, R, Dsum, Rsum, k=ceiling(0.95*nrow(R)), m=ncol(R), B=nrow(R), both=T){
  if(both){
    out <- R_bounds_both(ind, R, Dsum, Rsum, k, m, B)
  }else{
    out <- R_bounds_upper(ind, R, Dsum, Rsum, k, m, B)
  }
  return(out)
}




# Given D, R, I as given by R_ctrp_set,
# a vector of indices S, the significance level alpha
# the maximum number n_max of BAB iterations (0 if no BAB),
# the condition from_low (T if the BAB starts from the lowest statistic)
# and the condition first_rem (T if the BAB starts by removing the statistic),
# the function tests S. It returns non_rej (T if S is not rejected, F if it is rejected,
# and NULL if the number of steps reached the maximum before an indecisive outcome),
# and BAB, the number of iterations.

R_ctrp_test <- function(S, D, R, I, alpha=0.05, n_max=10000, from_low=T, first_rem=T){
  f <- ncol(D)
  B <- nrow(D)
  s <- length(S)
  m <- f-s
  k <- ceiling((1-alpha)*B)
  
  # if S=F:
  if(m==0){
    # lower and upper bounds (equal)
    low <- R_Q(rowSums(D), k, B)
    out <- list("non_rej"=!low, "BAB"=0)
    return(out)
  }
  
  g <- R_gen_sub(S, D, R, I, f, m, B, s)
  ind <- (0:m)
  cb <- R_compute_bounds(ind, g$R, g$Dsum, g$Rsum, k, m, B, both=T)
  
  # if there is a non-rejection or there are no indecisives
  if(cb$non_rej | length(cb$ind)==0){
    out <- list("non_rej"=cb$non_rej, "BAB"=0)
    return(out)
  }
  out <- R_ctrp_bab(cb$ind, g$D, g$R, g$I, g$Dsum, g$Rsum, k, m, B, n_max, from_low, first_rem)
  #out <- list("non_rej"=cb$non_rej, "ind"=cb$ind, "D"=g$D, "R"=g$R, "I"=g$I, "Dsum"=g$Dsum, "Rsum"=g$Rsum, "k"=k, "m"=m, "B"=B)
  
  #out <- R_generate_nodes(cb$ind, g$Dsum, g$Rsum, g$R, g$I, g$D, m-1, m, B, from_low, first_rem)

  return(out)
}



### ------------------------------------------------------------- ###




# Internal function - Branch and Bound
# Given a node (ind, Dsum, Rsum, R, I) and the dimenzion m=ncol(R)
# the initial matrix D_0 and m_0=ncol(D_0)
# it computes the two subspaces 1 (remove) and 2 (keep)
# when considering the lowest statistic
# In subspace 2, the indecisive sizes decrease of 1 unit
# and, since only the upper bound will be studied,
# zero is removed from the indecisive sizes.
# If zero is the only indecisive size, node 2 is not computed

R_gen_nodes_low <- function(ind, Dsum, Rsum, R, I, D_0, m=ncol(R), m_0=ncol(D_0), B=nrow(D_0)){
  h <- I[1,m+1] # index
  iD <- m_0 - m # index in D_0
  
  # Special case:
  # if the only indecisive size is 1, then ind_2=0 and we remove node 2
  if(length(ind)==1 & ind[1]==1){
    out <- list("ind_1"=ind, # 1 = remove case
                "Dsum_1"=Dsum[,2:(m+2)] - D_0[,iD],
                "Rsum_1"=matrix(rep(NA, B*(m+1)), ncol=m+1),
                "ind_2"=NULL, # 2 = keep case
                "Dsum_2"=NULL,
                "Rsum_2"=NULL,
                "R"=matrix(rep(NA, B*m), ncol=m),
                "I"=matrix(rep(NA, B*m), ncol=m))
    
    for(x in seq(B)){
      i <- match(h, I[x,])
      out$I[x,] <- I[x,-i]
      out$R[x,] <- R[x,-i]
      out$Rsum_1[x,] <- Rsum[x,-i]
      if(i<m+2){
        out$Rsum_1[x, i:(m+1)] <- out$Rsum_1[x, i:(m+1)] - D_0[x,iD]
      }
    }
    return(out)
  }
  
  ind_2 <- ind - 1
  if(ind_2[1]==0){ # if zero is an indecisive size, we remove it
    ind_2 <- ind_2[-1]
  }
  out <- list("ind_1"=ind, # 1 = remove case
              "Dsum_1"=Dsum[,2:(m+2)] - D_0[,iD],
              "Rsum_1"=matrix(rep(NA, B*(m+1)), ncol=m+1),
              "ind_2"=ind_2, # 2 = keep case
              "Dsum_2"=Dsum[,2:(m+2)],
              "Rsum_2"=matrix(rep(NA, B*(m+1)), ncol=m+1),
              "R"=matrix(rep(NA, B*m), ncol=m),
              "I"=matrix(rep(NA, B*m), ncol=m))
  
  for(x in seq(B)){
    i <- match(h, I[x,])
    out$I[x,] <- I[x,-i]
    out$R[x,] <- R[x,-i]
    out$Rsum_1[x,] <- Rsum[x,-i]
    out$Rsum_2[x,] <- out$Rsum_1[x,]
    if(i<m+2){
      out$Rsum_1[x, i:(m+1)] <- out$Rsum_1[x, i:(m+1)] - D_0[x,iD]
    }
    if(i>1){
      out$Rsum_2[x, 1:(i-1)] <- out$Rsum_2[x, 1:(i-1)] + D_0[x,iD]
    }
  }
  return(out)
}




# Internal function - Branch and Bound
# Given a node (ind, Dsum, Rsum, R, I) and the dimenzion m=ncol(R)
# the initial matrix D_0 and m_0=ncol(D_0)
# it computes the two subspaces 1 (remove) and 2 (keep)
# when considering the highest statistic
# In subspace 2, the indecisive sizes decrease of 1 unit

R_gen_nodes_high <- function(ind, Dsum, Rsum, R, I, D_0, m=ncol(R), m_0=ncol(D_0), B=nrow(D_0)){
  h <- I[1,1] # index
  iD <- m+1 # index in D_0
  
  out <- list("ind_1"=ind, # 1 = remove case
              "Dsum_1"=Dsum[,1:(m+1)],
              "Rsum_1"=matrix(rep(NA, B*(m+1)), ncol=m+1),
              "ind_2"=ind-1, # 2 = keep case
              "Dsum_2"=Dsum[,1:(m+1)] + D_0[,iD],
              "Rsum_2"=matrix(rep(NA, B*(m+1)), ncol=m+1),
              "R"=matrix(rep(NA, B*m), ncol=m),
              "I"=matrix(rep(NA, B*m), ncol=m))
  
  for(x in seq(B)){
    i <- match(h, I[x,])
    out$I[x,] <- I[x,-i]
    out$R[x,] <- R[x,-i]
    out$Rsum_1[x,] <- Rsum[x,-i]
    out$Rsum_2[x,] <- out$Rsum_1[x,]
    if(i<m+2){
      out$Rsum_1[x, i:(m+1)] <- out$Rsum_1[x, i:(m+1)] - D_0[x,iD]
    }
    if(i>1){
      out$Rsum_2[x, 1:(i-1)] <- out$Rsum_2[x, 1:(i-1)] + D_0[x,iD]
    }
  }
  return(out)
}




# Internal function - Branch and Bound
# Given a node (ind, Dsum, Rsum, R, I) and the dimenzion m=ncol(R)
# the initial matrix D_0 and m_0=ncol(D_0)
# the settings from_low and first_rem,
# it computes the two subspaces 1 and 2
# (1 corresponds to the subspace that will be explored first)

R_generate_nodes <- function(ind, Dsum, Rsum, R, I, D_0, m=ncol(R), m_0=ncol(D_0),
                           B=nrow(D_0), from_low=T, first_rem=T){
  if(from_low){
    out <- R_gen_nodes_low(ind, Dsum, Rsum, R, I, D_0, m, m_0, B)
  }else{
    out <- R_gen_nodes_high(ind, Dsum, Rsum, R, I, D_0, m, m_0, B)
  }
  
  # if we start by keeping the selected statistic, then we exchange the names:
  # 1 = keep case, 2 = remove case
  if(!first_rem){
    names(out)[1:6] <- c("ind_2", "Dsum_2", "Rsum_2", "ind_1", "Dsum_1", "Rsum_1")
  }
  return(out)
}




R_ctrp_bab <- function(ind_0, D_0, R_0, I_0, Dsum_0, Rsum_0,
                     k=ceiling(0.95*nrow(R)), m_0=ncol(D_0), B=nrow(D_0), n_max=10000,
                     from_low=T, first_rem=T){
  
  # when from_low=first_rem (either T or F), in the first loop we compute both bounds
  # and in the second only the upper
  both_first <- (from_low == first_rem)

  nodes <- list(list("ind"=ind_0, "Dsum"=Dsum_0, "Rsum"=Rsum_0, "R"=R_0, "I"=I_0))
  L <- length(nodes)
  BAB <- 0 # steps
  m <- m_0
  
  while(BAB<n_max){
    # we take the last element added to the list, generate the two branches
    m <- m-1
    g <- R_generate_nodes(nodes[[L]]$ind, nodes[[L]]$Dsum,
                        nodes[[L]]$Rsum, nodes[[L]]$R,
                        nodes[[L]]$I, D_0, m, m_0, B, from_low, first_rem)
    cond <- length(g$ind_1)>0
    
    # then we remove the original node and add the node 2 (if not null)
    if(length(g$ind_2)>0){
      nodes[[L]] <- list("ind"=g$ind_2, "R"=g$R, "I"=g$I, "Dsum"=g$Dsum_2, "Rsum"=g$Rsum_2)
    }else{
      nodes <- nodes[-L]
      L <- L-1
    }
    
    # we evaluate node 1, and iterate until a node 1 can be closed
    while(cond & BAB<n_max){
      BAB <- BAB + 1
      cb <- R_compute_bounds(g$ind_1, g$R, g$Dsum_1, g$Rsum_1, k, m, B, both_first)
      if(cb$non_rej){
        out <- list("non_rej"=T, "BAB"=BAB)
        return(out)
      }
      
      cond <- length(cb$ind)>0
      if(cond){
        m <- m-1
        # we generate new nodes, after updating the indecisive sizes (g$ind_1 becomes cb$ind)
        g <- R_generate_nodes(cb$ind, g$Dsum_1, g$Rsum_1, g$R, g$I, D_0, m, m_0, B, from_low, first_rem)
        cond <- length(g$ind_1)>0
        
        if(length(g$ind_2)>0){
          nodes <- append(nodes, list(list("ind"=g$ind_2, "R"=g$R, "I"=g$I, "Dsum"=g$Dsum_2, "Rsum"=g$Rsum_2)))
          L <- L+1
        }
      }
    }
    
    # if the list becomes empty before a non-rejection is found, we reject
    if(L==0){
      out <- list("non_rej"=F, "BAB"=BAB)
      return(out)
    }
    
    while(!cond & BAB<n_max){
      BAB <- BAB + 1
      
      if(is.matrix(nodes[[L]]$R)){
        m <- ncol(nodes[[L]]$R)
      }else{
        m <- 1
      }
      cb <- R_compute_bounds(nodes[[L]]$ind, nodes[[L]]$R,
                           nodes[[L]]$Dsum, nodes[[L]]$Rsum, k, m, B,
                           !both_first)
      if(cb$non_rej){
        out <- list("non_rej"=T, "BAB"=BAB)
        return(out)
      }
      
      cond <- length(cb$ind)>0
      if(cond){
        nodes[[L]]$ind <- cb$ind # update of indecisive sizes after evaluation
      }else{
        nodes <- nodes[-L] # the node is removed from the list
        L <- L-1
        
        # if the list becomes empty before a non-rejection is found, we reject
        if(L==0){
          out <- list("non_rej"=F, "BAB"=BAB)
          return(out)
        }
      }
    }
  }
  out <- list("non_rej"=NULL, "BAB"=BAB)
}
