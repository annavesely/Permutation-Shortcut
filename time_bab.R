
setwd("C:/Users/coshk/Desktop/Perm_project")

myseed <- 33
set.seed(33)


#given f and perc
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




# given s_size and s_active
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




get_time0 <- function(f_int=c(10,100,500,1000), perc_int=c(0, 0.01,0.1,0.2,0.5,0.8),
                     B_int=c(10,100,500,1000), s_size_int=c(0.01,0.1,0.2,0.5,0.8),
                     s_active_int=perc_int, alpha_int=c(0.05,0.2),
                     m=10, sd=5, n=20){
  
  set.seed(myseed)
  L <- length(f_int) * length(perc_int) * length(B_int) * length(s_size_int) *
    length(s_active_int) * length(alpha_int)
  M <- matrix(rep(NA, 8*L), ncol=8)
  colnames(M) <- c("f", "perc", "B", "s_size", "s_active", "alpha", "BAB", "time")
  
  fmax <- max(f_int)
  beta_active <- rnorm(fmax, m1, sd1)
  beta_inactive <- rep(0, fmax)
  
  i <- 0
  
  for(f in f_int){
    
    for(perc in perc_int){
      beta <- generate_beta(f, perc, beta_active, beta_inactive)
      
      for(B in B_int){
        G <- gt(n, f, B, 0, beta)
        c <- ctrp_set(G)
        
        for(s_size in s_size_int){
          
          for(s_active in s_active_int[s_active_int <= perc/s_size]){
            S <- generate_S(f, s_size, s_active)
            
            for(alpha in alpha_int){
              i <- i+1
              t <- system.time(te <- ctrp_test1(S, c$D, c$R, c$I, alpha))[3]
              M[i,] <- c(f, perc, B, s_size, s_active, alpha, te$BAB, t)
            }
          }
        }
      }
    }
  }
  M <- M[(1:i),]
  M <- M[M[,7]>0,]
  out <- list("M"=M, "beta_active"=beta_active, "beta_inactive"=beta_inactive, "n"=n)
  return(out)
}



t1 <- get_time0()
t1

write.table(t1, file="time1.txt", sep="", row.names=F, col.names=F)

my_data <- read.table("time1.txt", sep ="", header=F, dec=".")


get_time <- function(M, beta_active, beta_inactive, n=20){
  set.seed(myseed)
  w <- nrows(M)
  C <- rep(NA, w)
  M <- cbind(M, C, C, C, C)
  colnames(M) <- c("f", "perc", "B", "s_size", "s_active", "alpha", "BAB1", "time1",
                   "BAB2", "time2", "BAB3", "time3")
  
  for(j in(1:w)){
    f <- M[j,1]
    perc <- M[j,2]
    B <- M[j,3]
    s_size <- M[j,4]
    s_active <- M[j,5]
    alpha <- M[j,6]
    
    beta <- generate_beta(f, perc, beta_active, beta_inactive)
    G <- gt(n, f, B, 0, beta)
    c <- ctrp_set(G)
    S <- generate_S(f, s_size, s_active)
    t2 <- system.time(te2 <- ctrp_test2(S, c$D, c$R, c$I, alpha))[3]
    t3 <- system.time(te3 <- ctrp_test3(S, c$D, c$R, c$I, alpha))[3]
    M[j,(9:12)] <- c(te2$BAB, t2, te3$BAB, t3)
  }
  
  return(M)
}










