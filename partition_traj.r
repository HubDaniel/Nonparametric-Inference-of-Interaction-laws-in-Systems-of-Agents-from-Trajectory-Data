## partition_traj
# input: one single trajectory
# output1: max and min pairwise distance between agents on this trajectory
# output2: pairwise distance matrix of dimension N(N-1)/2 * L. each column is pairwise distance vector at time t_l

partition_traj <- function(traj, sys_info){
  
  print("start part traj")
  
  L <- dim(traj)[2]
  all_max_rs <- array(0,dim=c(sys_info$K,sys_info$K,L))
  all_min_rs <- array(0,dim=c(sys_info$K,sys_info$K,L))
  num_agents_each_type <- rep(0,sys_info$K)
  type_ind <- vector("list",sys_info$K)
  pdist_x <- vector("list",sys_info$K)
  for(i in 1:length(pdist_x)){
    pdist_x[[i]] <- vector("list",sys_info$K)
  }
  block_size <- sys_info$d * sys_info$N
  
 
  
  for ( l in 1:L){
    

    x_at_t <- traj[1:block_size,l]
    x_at_t <- matrix(x_at_t,nrow=sys_info$d,ncol=sys_info$N)
    pdist_x_at_t <- sqdist_mod(x_at_t)     # square root is taken in the function
    
    # print("dim(pdist_x_at_t)")
    # print(dim(pdist_x_at_t))
    
    
    for(k_1 in 1:sys_info$K){
      
      
      # print("k_1")
      # print(k_1)

      
      if(l==1 && k_1==1){
        agents_Ck1                 <- which(sys_info$type_info == k_1)
        type_ind[[k_1]]            <- agents_Ck1
        num_agents_Ck1             <- length(agents_Ck1)
        num_agents_each_type[k_1]  <- num_agents_Ck1
      }else{
        agents_Ck1                 <- type_ind[[k_1]]
        num_agents_Ck1             <- num_agents_each_type[k_1]
      }
      
      # print("agents_Ck1")
      # print(agents_Ck1)
      
      
      for(k_2 in 1:sys_info$K){
        
        # print("k_2")
        # print(k_2)
        
        if(l==1){
          if(k_2==1){
            agents_Ck2                 <- type_ind[[k_2]]
            num_agents_Ck2             <- num_agents_each_type[k_2]
          }else{
            agents_Ck2                 <- which(sys_info$type_info == k_2)
            type_ind[[k_2]]            <- agents_Ck2
            num_agents_Ck2             <- length(agents_Ck2)
            num_agents_each_type[k_2]  <- num_agents_Ck2
          }
        }else{
          agents_Ck2 <- type_ind[[k_2]]
          num_agents_Ck2 <- num_agents_each_type[k_2]
        }
        
        # print("agents_Ck2")
        # print(agents_Ck2)
        
        
        if(l==1){
          if(k_1==k_2){
            if(num_agents_Ck1 > 1){
              pdist_x[[k_1]][[k_2]] <- matrix(0,nrow = num_agents_Ck1 * (num_agents_Ck1 - 1)/2, ncol = L)
            }
          }else{
            pdist_x[[k_1]][[k_2]] <- matrix(0,nrow = num_agents_Ck1 * num_agents_Ck2, ncol = L)
          }
        }
        
        
        # print(paste0("pdist_x",k_1,k_2))
        # print(dim(pdist_x[[k_1]][[k_2]]))
        
        
        pdist_x_Ck1_Ck2 <- pdist_x_at_t[agents_Ck1,agents_Ck2]
      
        # print(paste0("pdist_x_ck1_ck2 ",k_1,k_2))
        # print((pdist_x_Ck1_Ck2))
        
        
        if(k_1 == k_2){
          if(num_agents_Ck1 > 1){
            

            pdist_x[[k_1]][[k_2]][,l] <- pdist_x_Ck1_Ck2[lower.tri(pdist_x_Ck1_Ck2)]
          }
        }else{
          pdist_x[[k_1]][[k_2]][,l] <- as.vector(pdist_x_Ck1_Ck2)
        }
        
        if(!is.null(pdist_x[[k_1]][[k_2]])){
          all_max_rs[k_1, k_2, l] <- max(pdist_x[[k_1]][[k_2]][, l])
          all_min_rs[k_1, k_2, l] <- min(pdist_x[[k_1]][[k_2]][, l])
        }
      }
    }
  }
  
  
  max_r <- apply(all_max_rs,c(1,2),max)
  max_r[max_r==0] <- 1
  min_r <- apply(all_min_rs,c(1,2),min)
  
  
  
  
  result <- list()
  result[["max_r"]] <- max_r
  result[["min_r"]] <- min_r
  result[["pdist_x"]] <- pdist_x
  
  print("finish part traj")
  
  return(result)
  
  
}