## partition_x_and_v

partition_x_and_v <- function(x,v,sys_info, learn_info){
  
  print("start part on x and v")
  L <- dim(x)[2]
  N <- sys_info$N
  d <- sys_info$d
  K <- sys_info$K
  type_info <- sys_info$type_info
  energy_pdist <- vector("list",K)
  energy_pdiff <- vector("list",K)
  
  for(i in 1:length(energy_pdist)){
    energy_pdist[[i]] <- vector("list",K)
    energy_pdiff[[i]] <- vector("list",K)
  }
  
  
  Psi_row_ind <- vector("list",K)
  num_agents_each_class <- rep(0,K)
  agents_class_indicator <- vector("list",K)
  row_ind_Phi_all_class <- vector("list",K)
  
  
  
  for (l in 1:L){
    
    #print(paste0("l=",l))
    
    x_at_t <- x[,l]
    x_at_t <- matrix(x_at_t,nrow = d, ncol = N)
    x_at_tt <- t(x_at_t)
    x_pdist <- pdist(x_at_tt)                 # pairwise distance in matrix form (symmetric), each row is one obs.
    x_pdiff <- find_pair_diff(x_at_t)

    
    for(k_1 in 1:K){
      #print(paste0("k_1=",k_1))
      if(l==1 && k_1==1){
        agents_Ck1                                    <- which(type_info == k_1)     
        agents_class_indicator[[k_1]]                 <- agents_Ck1     
        num_agents_Ck1                                <- length(agents_Ck1)  
        num_agents_each_class[k_1]                    <- num_agents_Ck1
      }else{
        agents_Ck1                                    <- agents_class_indicator[[k_1]]
        num_agents_Ck1                                <- num_agents_each_class[k_1]
      }
      
      
      
      #print("agents_Ck1")
      #print(agents_Ck1)
      
      
      
      row_ind_in_pdist <- 1:num_agents_Ck1 + (l-1)*num_agents_Ck1
      row_ind_in_pdiff <- 1:(num_agents_Ck1*d) + (l-1)*num_agents_Ck1*d
      
      if(l==1){
        row_ind_Phi <- rep(1:d,num_agents_Ck1)
        row_ind_Phi <- row_ind_Phi + rep((agents_Ck1-1)*d,each=d)
        row_ind_Phi_all_class[[k_1]] <- row_ind_Phi
      }else{
        row_ind_Phi <- row_ind_Phi_all_class[[k_1]]
      }
      
      #print("row_ind_Phi")
      #print(row_ind_Phi)
      
      
      if(l==1){
        Psi_row_ind[[k_1]] <- rep(0,L * num_agents_Ck1 * d)
      }
      
      Psi_row_ind[[k_1]][row_ind_in_pdiff] <- row_ind_Phi + (l - 1) * N * d
      
      
      
      #print(head(Psi_row_ind[[k_1]],50))
      
      
      
      
      
      for(k_2 in 1:K){
        #print(paste0("k_2=",k_2))
        if(l==1){
          if(k_2==1){
            agents_Ck2                                <- agents_class_indicator[[k_2]]
            num_agents_Ck2                            <- num_agents_each_class[k_2]
          }else{
            agents_Ck2                                <- which(type_info == k_2)
            agents_class_indicator[[k_2]]             <- agents_Ck2
            num_agents_Ck2                            <- length(agents_Ck2)
            num_agents_each_class[k_2]                <- num_agents_Ck2
          }
        }else{
          agents_Ck2                                  <- agents_class_indicator[[k_2]]
          num_agents_Ck2                              <- num_agents_each_class[k_2]
        }
        
        if(l==1){
          energy_pdist[[k_1]][[k_2]] <- matrix(0,nrow = L*num_agents_Ck1 ,  ncol = num_agents_Ck2)
          energy_pdiff[[k_1]][[k_2]] <- matrix(0,nrow = L*num_agents_Ck1*d, ncol = num_agents_Ck2)
        }
        
        x_pdist_Ck1_Ck2 <- x_pdist[agents_Ck1, agents_Ck2]
        energy_pdist[[k_1]][[k_2]][row_ind_in_pdist,] <- x_pdist_Ck1_Ck2
        x_pdiff_Ck1_Ck2 <- x_pdiff[row_ind_Phi, agents_Ck2]
        energy_pdiff[[k_1]][[k_2]][row_ind_in_pdiff,] <- x_pdiff_Ck1_Ck2
        
      }
      
    }

  }
  
  
  
  result <- list()
  result[["energy_pdiff"]] <- energy_pdiff
  result[["energy_pdist"]] <- energy_pdist
  result[["Psi_row_ind"]] <- Psi_row_ind

  
  print("finish part on x and v")
  
  return(result)
}


# test

#rrr <- partition_x_and_v(x=obs$traj[[1]],v=c(),sys_info, learn_info)
# dim(rrr$energy_pdiff[[2]][[2]])


