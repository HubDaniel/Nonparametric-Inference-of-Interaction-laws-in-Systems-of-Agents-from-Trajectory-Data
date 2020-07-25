## uniform_learn_interactions_on_x_and_v    step1: assemble Phi matrix, d^m etc.  step2: solve alpha
# input: those through learn_interactions_on_x_and_v
# output: list of results. 


uniform_learn_interactions_on_x_and_v <- function(x, v, xi, dot_xv, time_vec, sys_info, learn_info,Riemann){
  print("uniform learn interact on x and v")
  
  M <- dim(x)[3]
  
  result <- serial_assemble_x_and_v(x, v, xi, dot_xv, time_vec, sys_info, learn_info,Riemann)
  Estimator <- result$Estimator
  extra <- result$extra


  T_L <- time_vec[length(time_vec)]

  result2 <- solve_for_interactions_on_x_and_v_by_others( Estimator[["Phi"]],Estimator[["rhs_vec"]],
                                                          T_L,M,learn_info$solver_type)
  Estimator[["alpha"]] <- result2$alpha_vec
  
  
 
  
  
  opt_val <- result2$opt_val
  
  Estimator[["Info.phiSingVals"]] <- svd(Estimator[["Phi"]])$d
  Estimator[["Info.phiCond"]] <- Estimator[["Info.phiSingVals"]][1] / 
    Estimator[["Info.phiSingVals"]][length(Estimator[["Info.phiSingVals"]])]
  
  if (T_L == 0){
    Estimator[["emp_err"]] <- opt_val + extra[["rhs_in_l2_norm_sq"]]/M
  }else{
    Estimator[["emp_err"]] <- opt_val + extra[["rhs_in_l2_norm_sq"]]/(T_L * M)
  }
    
  final_result <- list()
  final_result[["Estimator"]] <- Estimator
  final_result[["extra"]] <- extra
  
  return(final_result)
}