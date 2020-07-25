# generateICs

generateICs <- function(M, sys_info){
  
  y_inits <- matrix(0,ncol=M,nrow=sys_info$d * sys_info$N)
  for (m in 1:M){
    y_init <- sys_info$mu0()
    y_inits[,m] <- y_init
  }
  return(y_inits)
}



#test
#generateICs(2,sys_info)
