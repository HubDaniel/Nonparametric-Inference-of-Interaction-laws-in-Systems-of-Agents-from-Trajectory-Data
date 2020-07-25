## true kernel function phi(.)
# input: r: vector of pair wise distance
# output: phi(r) 

kernel_function <- function(r,epsilon = 10,sigma = 1){
  result <- rep(0,length(r))     # result holder
  ind <- r!=0
  rinv <- sigma/r[ind]
  result[ind] <- 24 * epsilon/(sigma^2) * (-2*(rinv^14) + (rinv^8))
  
  ind <- r==0
  result[ind] <- -Inf
  
  return(result)
}



