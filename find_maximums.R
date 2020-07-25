## find_maximums

find_maximums <- function(traj, sys_info){
  L <- ncol(traj)
  rs_each_type   <- rep(0, sys_info$K)
  one_block <- sys_info$N * sys_info$d
  agent_each_type <- matrix(0, nrow = sys_info$K, ncol = sys_info$N)

  
  for(l in 1:L){
    x_at_tl <- traj[,l]
    x_at_tl <- matrix(x_at_tl,ncol = sys_info$N,nrow = sys_info$d)
    rs_at_tl <- colSums(x_at_tl^2)^(0.5)
    
    
    for(k in 1:sys_info$K){
      if(l==1){
        agent_Ck              <- sys_info$type_info == k;
        agent_each_type[k,]   <- agent_Ck
      }else{
        agent_Ck              <- agent_each_type[k,]
      }
      rs_each_type[k]         <- max(rs_at_tl[as.logical(agent_Ck)])
    }
    
    
    if(l==1){
      max_rs <- matrix(rep(rs_each_type,sys_info$K),ncol=sys_info$K) + 
        matrix(rep(rs_each_type,each=sys_info$K),ncol=sys_info$K)
    }else{
      a_max <- matrix(rep(rs_each_type,sys_info$K),ncol=sys_info$K) + 
        matrix(rep(rs_each_type,each=sys_info$K),ncol=sys_info$K)
      max_rs <- pmax(max_rs, a_max)
    }
  }
  # it is a K*K matrix where the diagonal contains 2 times the || of the max of each type, 
  # off-diagonal is max|type1| +  max|type2|
  
  return(max_rs) 
}


# test
# traj <- matrix(1:60,ncol = 3)
# traj
# find_maximums(traj,sys_info)

