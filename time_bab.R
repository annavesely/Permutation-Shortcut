
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



# NULL: 0
# TRUE: 1
# FALSE: -1
nr_num <- function(nr){
  ifelse(length(nr)==0, 0, ifelse(nr, 1, -1))
}




ctrp_test0 <- function(S, D, R, I, alpha=0.05, n_max=10000){
  f <- ncol(D)
  B <- nrow(D)
  s <- length(S)
  m <- f-s
  k <- ceiling((1-alpha)*B)
  
  # if S=F:
  if(m==0){
    # lower and upper bounds (equal)
    low <- Q(rowSums(D), k, B)
    out <- list("cont"=F, "non_rej"=!low, "BAB"=0)
    return(out)
  }
  
  g <- gen_sub(S, D, R, I, f, m, B, s)
  ind <- (0:m)
  cb <- compute_bounds(ind, g$R, g$Dsum, g$Rsum, k, m, B, both=T)
  
  # if there is a non-rejection or there are no indecisives
  if(cb$non_rej | length(cb$ind)==0){
    out <- list("cont"=F, "non_rej"=cb$non_rej, "BAB"=0)
    return(out)
  }
  out <- list("cont"=T, "non_rej"=cb$non_rej, "ind"=cb$ind, "D"=g$D, "R"=g$R, "I"=g$I,
              "Dsum"=g$Dsum, "Rsum"=g$Rsum, "k"=k, "m"=m, "B"=B)
  return(out)
}



get_time0 <- function(te, n_max){
  # 1=lr, 2=lk, 3=hr, 4=hk
  # 1: remove lowest (from_low=T, first_rem=T)
  t1 <- system.time(b1 <- ctrp_bab(te$ind, te$D, te$R, te$I, te$Dsum, te$Rsum, te$k, te$m, te$B, n_max, T, T))[3]
  nr1 <- nr_num(b1$non_rej)
  
  # 2: keep lowest (from_low=T, first_rem=F)
  t2 <- system.time(b2 <- ctrp_bab(te$ind, te$D, te$R, te$I, te$Dsum, te$Rsum, te$k, te$m, te$B, n_max, T, F))[3]
  nr2 <- nr_num(b2$non_rej)
  
  # 3: remove highest (from_low=F, first_rem=T)
  t3 <- system.time(b3 <- ctrp_bab(te$ind, te$D, te$R, te$I, te$Dsum, te$Rsum, te$k, te$m, te$B, n_max, F, T))[3]
  nr3 <- nr_num(b3$non_rej)
  
  # 4: keep highest (from_low=F, first_rem=F)
  t4 <- system.time(b4 <- ctrp_bab(te$ind, te$D, te$R, te$I, te$Dsum, te$Rsum, te$k, te$m, te$B, n_max, F, F))[3]
  nr4 <- nr_num(b4$non_rej)
  
  # T (1) if the outcome differs
  flag0 <- !(nr1==nr2 & nr2==nr3 & nr3==nr4)
  flag <- ifelse(flag0, 1, 0)
  
  c(t1, t2, t3, t4, b1$BAB, b2$BAB, b3$BAB, b4$BAB, nr1, nr2, nr3, nr4, flag)
}





get_time <- function(f_int=c(10,50,100), perc_int=c(0, 0.01,0.1,0.2,0.5,0.8),
                     B_int=c(10,50,100), s_size_int=c(0.01,0.1,0.2,0.5,0.8, 1),
                     s_active_int=perc_int, alpha_int=c(0.05,0.2),
                     m=10, sd=5, n_max=10000){
  
  set.seed(myseed)
  L <- length(f_int) * length(perc_int) * length(B_int) * length(s_size_int) *
    length(s_active_int) * length(alpha_int)
  M <- matrix(rep(NA, 19*L), ncol=19)
  colnames(M) <- c("s_active", "perc", "s_size", "f", "B", "alpha",
                  "t1", "t2", "t3", "t4", "BAB1", "BAB2", "BAB3", "BAB4",
                  "nr1", "nr2", "nr3", "nr4", "flag")
  # 1=lr, 2=lk, 3=hr, 4=hk
  
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
          
          for(s_active in s_active_int[s_active_int <= perc/s_size &
                                       s_active_int >= 1 - (1 - perc)/s_size]){
            S <- generate_S(f, s_size, s_active)
            
            for(alpha in alpha_int){
              te <- ctrp_test0(S, c$D, c$R, c$I, alpha, n_max)
              
              if(te$cont){
                i <- i+1
                M[i,] <- c(s_active, perc, s_size, f, B, alpha, get_time0(te, n_max))
              }
            }
          }
        }
      }
    }
  }
  M <- M[(1:i),]
  return(M)
}




f_int <- c(100)
perc_int <- c(0.2)
B_int <- c(100)
s_size_int <- c(0.2)
alpha_int <- c(0.05)
m <- 10
sd <- 5
n_max <- 10000
s_active_int <- c(1)





mytime <- get_time(f_int, perc_int, B_int, s_size_int, s_active_int, alpha_int, m, sd, n_max)
mytime[15:19] # flag
cbind(mytime[,1:6], mytime[,15], mytime[,11:14]) # table
mytime[,7:10] # times

c(mytime[1:6], mytime[15], mytime[11:14])




### -------------------------------------------------------------- ###



get_time2 <- function(f_int=c(10,50,100), perc_int=c(0, 0.01,0.1,0.2,0.5,0.8),
                     B_int=c(10,50,100), s_size_int=c(0.01,0.1,0.2,0.5,0.8, 1),
                     s_active_int=perc_int, alpha_int=c(0.05,0.2),
                     m=10, sd=5, n_max=10000){
  
  set.seed(myseed)
  L <- length(f_int) * length(perc_int) * length(B_int) * length(s_size_int) *
    length(s_active_int) * length(alpha_int)
  M <- matrix(rep(NA, 9*L), ncol=9)
  colnames(M) <- c("s_active", "perc", "s_size", "f", "B", "alpha", "nr",
                   "BAB", "t")
  # 1=lr, 2=lk, 3=hr, 4=hk
  
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
          
          for(s_active in s_active_int[s_active_int <= perc/s_size &
                                       s_active_int >= 1 - (1 - perc)/s_size]){
            S <- generate_S(f, s_size, s_active)
            
            for(alpha in alpha_int){
              i <- i+1
              t <- system.time(te <- ctrp_test(S, c$D, c$R, c$I, alpha, n_max, from_low=F, first_rem=T))[3]
              M[i,] <- c(s_active, perc, s_size, f, B, alpha, nr_num(te$non_rej), te$BAB, t)
            }
          } 
        }
      }
    }
  }
  M <- M[(1:i),]
  return(M)
}



w <- function(x){
  out <- round(x*1000,1)
  return(out)
}




f_int <- c(250)
alpha_int <- c(0.05)


perc_int <- c(0, 0.01,0.1,0.2,0.5,0.8)
s_size_int <- c(0.01,0.1,0.2,0.5,0.8, 1)
m <- 10
sd <- 5
n_max <- 10000
s_active_int <- c(0, 0.01,0.1,0.2,0.5,0.8)


mytime1 <- get_time2(f_int, perc_int, B_int=c(10), s_size_int, s_active_int, alpha_int, m, sd, n_max)
mytime2 <- get_time2(f_int, perc_int, B_int=c(50), s_size_int, s_active_int, alpha_int, m, sd, n_max)
mytime3 <- get_time2(f_int, perc_int, B_int=c(100), s_size_int, s_active_int, alpha_int, m, sd, n_max)


w(mean(mytime1[,9]))
w(mean(mytime2[,9]))
w(mean(mytime3[,9]))

w(mean(c(mytime1[,9], mytime2[,9], mytime3[,9])))








# f= 10, B=10,100,500
# f=100, B=10, 100

# da fare:
# f= 10, B=50
# f=50, B=10,50,100
# f= 100, B=50




get_time(f_int=c(100), perc_int=c(1),
                     B_int=c(10), s_size_int=c(0.5),
                     s_active_int=c(1), alpha_int=c(0.05),
                     m=10, sd=5, n_max=100000)




f <- 100
B <- 10
perc <- 100/100
s_size <- 50/100
s_active <- 100/100
alpha <- 0.20
m <- 10
sd <- 5
n_max <- 100000

beta_active <- rnorm(f, m, sd)
beta_inactive <- rep(0, f)
beta <- generate_beta(f, perc, beta_active, beta_inactive)
G <- gt(f, f, B, 0, beta)
c <- ctrp_set(G)
S <- generate_S(f, s_size, s_active)

te <- ctrp_test(S, c$D, c$R, c$I, alpha, n_max, from_low=T, first_rem=T)

b1 <- ctrp_bab(te$ind, te$D, te$R, te$I, te$Dsum, te$Rsum, te$k, te$m, te$B, n_max, T, T)
b1
b2 <- ctrp_bab(te$ind, te$D, te$R, te$I, te$Dsum, te$Rsum, te$k, te$m, te$B, n_max, T, F)
b2
b3 <- ctrp_bab(te$ind, te$D, te$R, te$I, te$Dsum, te$Rsum, te$k, te$m, te$B, n_max, F, T)
b3
b4 <- ctrp_bab(te$ind, te$D, te$R, te$I, te$Dsum, te$Rsum, te$k, te$m, te$B, n_max, F, F)
b4






c(mytime[2:6], mytime[15], mytime[11:14]) # table


# f=10, B=10 -> nessuno
# f=10, B=100 -> 4
get_time(mytime$M, mytime$beta_active, mytime$beta_inactive)



#write.table(mytime$M, file="time1.txt", sep="", row.names=F, col.names=F)
#my_data <- read.table("time1.txt", sep ="", header=F, dec=".")






