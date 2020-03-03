
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


f <- 1000
perc <- 0.1
B <- 10
s_size <- 0.2
s_active <- 0.5
alpha <- 0.2
m <- 10
sd <- 5
n_max <- 200
n <- f

get_time0 <- function(f_int=c(10,100,500,1000), perc_int=c(0, 0.01,0.1,0.2,0.5,0.8),
                     B_int=c(10,100,500,1000), s_size_int=c(0.01,0.1,0.2,0.5,0.8),
                     s_active_int=perc_int, alpha_int=c(0.05,0.2),
                     m=10, sd=5, n_max=1000){
  
  set.seed(myseed)
  L <- length(f_int) * length(perc_int) * length(B_int) * length(s_size_int) *
    length(s_active_int) * length(alpha_int)
  M <- matrix(rep(NA, 9*L), ncol=9)
  colnames(M) <- c("f", "perc", "B", "s_size", "s_active", "alpha", "non_rej", "BAB", "time")
  
  fmax <- max(f_int)
  beta_active <- rnorm(fmax, m, sd)
  beta_inactive <- rep(0, fmax)
  
  i <- 0
  
  for(f in f_int){
    
    for(perc in perc_int){
      beta <- generate_beta(f, perc, beta_active, beta_inactive)
      
      for(B in B_int){
        G <- gt(f, f, B, 0, beta)
        c <- ctrp_set(G)
        
        for(s_size in s_size_int){
          
          for(s_active in s_active_int[s_active_int <= perc/s_size]){
            S <- generate_S(f, s_size, s_active)
            
            for(alpha in alpha_int){
              i <- i+1
              t <- system.time(te <- ctrp_test(S, c$D, c$R, c$I, alpha, n_max, from_low=T, first_rem=T))[3]
              nr <- ifelse(length(te$non_rej)==0, 0, ifelse(te$non_rej, 1, -1))
              M[i,] <- c(f, perc, B, s_size, s_active, alpha, nr, te$BAB, t)
            }
          }
        }
      }
    }
  }
  M <- M[(1:i),]
  M <- M[M[,8]>0,]
  out <- list("M"=M, "beta_active"=beta_active, "beta_inactive"=beta_inactive)
  return(out)
}


mytime <- get_time0(f_int=c(10), perc_int=c(0, 0.01,0.1,0.2,0.5,0.8),
                    B_int=c(100), s_size_int=c(0.01,0.1,0.2,0.5,0.8),
                    s_active_int=c(0, 0.01,0.1,0.2,0.5,0.8), alpha_int=c(0.05,0.2),
                    m=10, sd=5, n_max=1000)

# f=10, B=10 -> nessuno
# f=10, B=100 -> 4
get_time(mytime$M, mytime$beta_active, mytime$beta_inactive)



#write.table(mytime$M, file="time1.txt", sep="", row.names=F, col.names=F)
#my_data <- read.table("time1.txt", sep ="", header=F, dec=".")


get_time <- function(M, beta_active, beta_inactive){
  set.seed(myseed)
  w <- nrow(M)
  C <- rep(NA, w)
  M <- cbind(M, C, C, C, C, C, C)
  # 1=lr, 2=lk, 3=hr, 4=hk
  colnames(M) <- c("f", "perc", "B", "s_size", "s_active", "alpha", "flag", "BAB1", "time1",
                   "BAB2", "time2", "BAB3", "time3", "BAB4", "time4")
  
  for(j in(1:w)){
    f <- M[j,1]
    perc <- M[j,2]
    B <- M[j,3]
    s_size <- M[j,4]
    s_active <- M[j,5]
    alpha <- M[j,6]
    nr1 <- M[j,7]
    
    beta <- generate_beta(f, perc, beta_active, beta_inactive)
    G <- gt(f, f, B, 0, beta)
    c <- ctrp_set(G)
    S <- generate_S(f, s_size, s_active)
    t2 <- system.time(te2 <- ctrp_test(S, c$D, c$R, c$I, alpha, n_max, from_low=T, first_rem=F))[3]
    t3 <- system.time(te3 <- ctrp_test(S, c$D, c$R, c$I, alpha, n_max, from_low=F, first_rem=T))[3]
    t4 <- system.time(te4 <- ctrp_test(S, c$D, c$R, c$I, alpha, n_max, from_low=F, first_rem=F))[3]
    nr2 <- ifelse(length(te2$non_rej)==0, 0, ifelse(te2$non_rej, 1, -1))
    nr3 <- ifelse(length(te3$non_rej)==0, 0, ifelse(te3$non_rej, 1, -1))
    nr4 <- ifelse(length(te4$non_rej)==0, 0, ifelse(te4$non_rej, 1, -1))

    M[j,(10:15)] <- c(te2$BAB, t2, te3$BAB, t3, te4$BAB, t4)
    M[j,7] <- !(nr1 == nr2 & nr2 == nr3 & nr3 == nr4)
    # F if there are no problems
    
    
  }
  
  return(M)
}



f <- 100
perc <- 0.20
B <- 10
s_size <- 0.2
s_active <- 0.8
alpha <- 0.2


ctrp_test(S, c$D, c$R, c$I, alpha, n_max, from_low=F, first_rem=F)




