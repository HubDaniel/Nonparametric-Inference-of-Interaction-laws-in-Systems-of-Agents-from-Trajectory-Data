## learn interactions and construct phiEhat
# input: the same as those through learn_from_dynamics
# output: list of results

learn_interactions_on_x_and_v <- function(x, v, xi, dot_xv, time_vec, sys_info, learn_info,Riemann){
  print("learn interact on x and v")
  
  
  result <- uniform_learn_interactions_on_x_and_v(x,  v, xi, dot_xv, time_vec, sys_info, learn_info,Riemann)
  Estimator <- result$Estimator
  extra <- result$extra
  
  
  #print(Estimator[["Ebasis"]])
  
  result2  <- LinearCombinationBasis( sys_info$K, Estimator[["Ebasis"]], Estimator[["alpha"]] )
  
  Estimator[["phiEhat"]] <- result2[["f_ptr"]]   # phiEhat is a list with "y" and "yp"
  
  final_result <- list()
  final_result[["Estimator"]] <- Estimator
  final_result[["extra"]] <- extra
  
  return(final_result)
}