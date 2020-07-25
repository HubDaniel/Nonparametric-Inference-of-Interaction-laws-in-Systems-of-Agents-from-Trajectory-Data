## serial_assemble_x_and_v


serial_assemble_x_and_v <- function(x, v, xi, dot_xv, time_vec, sys_info, learn_info,Riemann){
  
  print("serial ass x and v")
  
  Estimator <- list()
  
  # step 1: find maximum range in all trajectories (all data) Rs for basis construction
  M <- dim(x)[3]
  agents_info <- getAgentInfo(sys_info)
  N_per_class <- histc(sys_info$type_info,seq(from = 0.5, to = sys_info$K+1, by = 1))$cnt
  N_per_class[length(N_per_class)-1] <- N_per_class[length(N_per_class)-1] + N_per_class[length(N_per_class)]  # add the last one to the second last one because of the mechanism of histc
  N_per_class <- N_per_class[-length(N_per_class)]  # then remove the last one
  
  Mtrajs <- x
  max_rs <- array(0,dim = c(sys_info$K,sys_info$K,M))
  
  for (m in 1:M){
    traj <- Mtrajs[,,m]
    output <- find_maximums(traj, sys_info)   
    max_rs[,,m] <- output
  }
  Rs <- apply(max_rs,c(1,2),max)                      
  
  
  # step 2: construct basis
  Estimator[["Ebasis"]] <-  uniform_basis_by_class( Rs, sys_info$K, learn_info$Ebasis_info  )
  
  
  # step 3: construct estimation vectors and matrices (phi matrix // d^m // ...)
  num_Estimator.Phi_cols <- 0   # total number of basis functions
  
  for(i in 1:sys_info$K){
    num_Estimator.Phi_cols <- num_Estimator.Phi_cols + sum(sapply(Estimator$Ebasis[[i]],lengths)["f",])
  }
  
  print("num_Estimator.Phi_cols")
  print(num_Estimator.Phi_cols)
  
  Estimator.Phi <- matrix(0,ncol = num_Estimator.Phi_cols, nrow = num_Estimator.Phi_cols)
  
  Estimator.rhs_vec <- rep(0, num_Estimator.Phi_cols)
  
  rhs_in_l2_norm_sq <- 0
  
  for (m in 1:M){
    print(paste0("serial assemble -- constr estimate matrix and vec...  m = ",m,"/",M))
    
    one_x <- x[,,m]
    one_v <- c()
    one_xi <- c()
    one_dot_xv <- c()
    
    result <- one_step_assemble_on_x_and_v(one_x, one_v, one_xi,one_dot_xv, time_vec, agents_info, Estimator$Ebasis, 
                                          sys_info, learn_info, Riemann)

    
    
    
    one_energy_Estimator.Psi <- result$energy_Psi
    one_align_Estimator.Psi <- result$align_Psi
    #one_F_vec <- result$F_vec
    one_d_vec <- result$d_vec
    
    one_rhs_vec <- one_d_vec
    one_Psi <- one_energy_Estimator.Psi

    
    rhs_in_l2_norm_sq <- rhs_in_l2_norm_sq + norm(one_rhs_vec,type="2")^2
    Estimator.Phi <- Estimator.Phi + t(one_Psi) %*% one_Psi
    Estimator.rhs_vec <- Estimator.rhs_vec + t(one_Psi) %*% one_rhs_vec
  }
  
  

  
  extra.rhs_in_l2_norm_sq     = rhs_in_l2_norm_sq                                                
  extra.rhoLTemp              = c()
  
  
  
  #print(Estimator)
  

  Estimator[["Phi"]] <- Estimator.Phi
  Estimator[["rhs_vec"]] <- Estimator.rhs_vec
  
  extra <- list()
  extra[["rhs_in_l2_norm_sq"]] <- extra.rhs_in_l2_norm_sq
  extra[["rhoLTemp"]] <- extra.rhoLTemp
  
  final_result <- list()
  final_result[["Estimator"]] <- Estimator
  final_result[["extra"]] <- extra
  
  return(final_result)
}













