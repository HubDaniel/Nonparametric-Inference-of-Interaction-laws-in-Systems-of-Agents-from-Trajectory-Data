## get_Riemann_based_time_modifier


get_Riemann_based_time_modifier <- function(N, d, time_vec, Riemann){
  L <- length(time_vec)
  time_mod <- time_vec[2:L]-time_vec[1:(L-1)]
  if(Riemann == 0){
    time_mod <- rep(1,L)
  }else if(Riemann == 1){
    time_mod <- c(time_mod,0)
  }else if(Riemann == 2){
    time_mod <- c(0,time_mod)
  }else if(Riemann ==3){
    time_mod <- (c(time_mod,0)+c(0,time_mod))/2
  }
  time_mod <- sqrt(time_mod)
  time_mod <- rep(time_mod,each=N*d)
  
  return(time_mod)
}