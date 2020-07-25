# traj_norm

# input: two individual traj, each is a matrix of dim dN*L


traj_norm <- function(traj_1, traj_2,sys_info){
  
  N_ks <- rep(0,sys_info$K)
  C_ks <- matrix(0,nrow=sys_info$K, ncol=sys_info$N)
  
  for(k in 1:sys_info$K){
    agents_Ck <- sys_info$type_info==k   # true and false
    C_ks[k,] <- agents_Ck                # matrix of K*N, each row contains true or false (0 or 1)
    N_ks[k] <- sum(agents_Ck)
  }
  
  L <- dim(traj_1)[2]
  if(dim(traj_2)[2]!=L) print("traj_2 must have the same length in time!!!!!")
  
  state_diff <- rep(0,L)
  
  for(l in 1:L){
    state_diff[l] <- state_norm(traj_1[,l] - traj_2[,l], sys_info, N_ks, C_ks)
  }
  
  traj_dist <- max(state_diff)
  
  return(traj_dist)
  
}