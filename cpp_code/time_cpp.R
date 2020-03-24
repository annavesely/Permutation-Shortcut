
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


setwd('C:/Users/coshk/Desktop/Perm_project/cpp_code')

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
