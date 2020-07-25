# state_norm

# input: state_vec: vector of length d*N, this is the diff of vector containing the location of N agents at time t_l

# output: scalor value of 
#                         ||X||_2^2 = \sum_{k = 1}^K 1/N_k \sum_{i = 1}^N_k||x_i||_2^2
# at time t_l


state_norm <- function(state_vec, sys_info, N_ks, C_ks){
  x <- matrix(state_vec,nrow = sys_info$d,ncol = sys_info$N)
  the_sum <- colSums(x^2)  # square component wise and sum over the rows, vector of length N
  state_norm <- rep(0,sys_info$K)
  for(k in 1:sys_info$K){
    state_norm[k] <- sum(the_sum[as.logical(C_ks[k,])])/N_ks[k]
    # ||X||_2^2 = \sum_{k = 1}^K 1/N_k \sum_{i = 1}^N_k||x_i||_2^2
  }
  state_norm <- sqrt(sum(state_norm))
  return(state_norm)
}