## assemble_the_learning_matrix_on_x_and_v

assemble_the_learning_matrix_on_x_and_v <- function(energy_pdist, energy_pdiff, energy_basis,time_vec,
                                                    agents_info, Phi_row_ind, sys_info, learn_info){
  print("start assm the learning matrix on x and v")
  
  L <- length(time_vec)
  
  # mimic the cellfun in matlab
  energy_each_class <- vector("list",sys_info$K)
  for(i in 1:length(energy_each_class)){
    energy_each_class[[i]] <- vector("list",sys_info$K)
  }
  
  for(i in 1:sys_info$K){
    for(j in 1:sys_info$K){
      energy_each_class[[i]][[j]] <- length(energy_basis[[i]][[j]]$f)
    }
  }
  
  #print("num energy each class")
  #print(energy_each_class)

  total_energy_basis <- sum(unlist(energy_each_class))
  
  energy_Phi <- matrix(0, nrow = L*sys_info$N*sys_info$d, ncol = total_energy_basis)
  
  num_prev_energy <- 0
  
  for(k_1 in 1:sys_info$K){
    for(k_2 in 1:sys_info$K){
      ind_1 <- num_prev_energy + 1
      ind_2 <- num_prev_energy + energy_each_class[[k_1]][[k_2]]
      
      # print("ind_1")
      # print(ind_1)
      # print("ind_2")
      # print(ind_2)
      
      
      num_prev_energy <- num_prev_energy + energy_each_class[[k_1]][[k_2]]
      if((k_1==k_2) && (agents_info$num_agents[[k_2]] == 1)){
        energy_Phi[Phi_row_ind[[k_1]],ind_1:ind_2] <- 0
      }else{
        class_influence <- find_class_influence(energy_basis[[k_1]][[k_2]], energy_pdist[[k_1]][[k_2]],
                                                energy_pdiff[[k_1]][[k_2]],sys_info$d, agents_info$num_agents[k_2],
                                                sys_info$kappa[k_2],learn_info$Ebasis_info)
        
        
        # print("dim class_influence")
        # print(dim(class_influence))
        # print("dim energy_Phi[Phi_row_ind[[k_1]], ind_1:(ind_1+size(class_influence,2)-1)]")
        # print(dim(energy_Phi[Phi_row_ind[[k_1]], ind_1:(ind_1+size(class_influence,2)-1)]))
        
        
        #print(dim(energy_Phi))
        #print(sum((Phi_row_ind[[k_1]])!=c(1:1274)))
        #print(dim(energy_Phi[Phi_row_ind[[k_1]], ind_1:(ind_1+dim(class_influence)[2]-1)]))
        energy_Phi[Phi_row_ind[[k_1]], ind_1:(ind_1+dim(class_influence)[2]-1)] <- class_influence
      }

    }
  }
 
  result <- list()
  result[["energy_Phi"]] <- energy_Phi
  result["align_Phi"] <- list(NULL)
  
  print("finish assm the learning matrix on x and v")
  
  return(result)
}