## one_step_assemble_on_x_and_v

one_step_assemble_on_x_and_v <- function(x, v, xi, dot_xv, time_vec,agents_info, energy_basis, sys_info,
                                         learn_info,Riemann){
  print("start one step ass on x and v")
  #print(dim(one_x))
  # step 1:  find d^m
  d_vec <- approximate_derivative_of_x_or_v( x,v,time_vec,sys_info)

  
  # step 2:
  # compute pairwise distance, pairwise differences, regulator, maximum interaction radii, 
  # row indices in the learning matrix to put the assembled data
  
  result <- partition_x_and_v( x,v,sys_info,learn_info )                 
  energy_pdist <- result$energy_pdist
  energy_pdiff <- result$energy_pdiff
  Psi_row_ind <- result$Psi_row_ind
  
  #print(Psi_row_ind)
  
  # step 3: find F vector
  L <- dim(x)[2]
  dn <- dim(x)[1]
  F_vec <- rep(0,dn*L)
  
  # step 4: find Phi (i.e. Psi) matrix
  result2 <- assemble_the_learning_matrix_on_x_and_v(energy_pdist, energy_pdiff, energy_basis,time_vec,
                                                     agents_info, Psi_row_ind, sys_info, learn_info)
  energy_Psi <- result2[["energy_Phi"]]
  align_Psi <- result2[["align_Phi"]]
  
  # print("energy_Psi")
  # print(energy_Psi)
  
  
  # step 5: scale the quantities
  result3 <- scale_the_quantities(F_vec, d_vec, energy_Psi, sys_info$N, sys_info$K,sys_info$type_info,time_vec,Riemann)
  #F_vec <- result3$F_vec
  d_vec <- result3$d_vec
  energy_Psi <- result3$energy_Psi
  
  
  final_result <- list()
  #final_result[["F_vec"]] <- F_vec
  final_result[["d_vec"]] <- d_vec
  final_result[["energy_Psi"]] <- energy_Psi
  final_result["align_Psi"] <- list(NULL)
  
  print("finish one step ass on x and v")
 # print(energy_Psi)

  return(final_result)
  
}