## 
# input: obs_data: 3 dimension array of dimension dN*L*M. first dimension dN pooled data. second dimension time t1,...,tL.
#        third dimension trajectory 1,...,M
# output: list of learning result: list[[1]] -- Estimator(Ebasis,Abasis,phi,rhs_vec,alpha,emp_error,phiEhat,..)
#                                  list[[2]] -- extra_xv(rhoLTemp,rhs_in_l2_norm_sq)

library(MASS)
library(rdist)
library(pracma)

source("getAgentInfo.r")
source("learn_interactions_on_x_and_v.r")
source("uniform_learn_interactions_on_x_and_v.r")
source("LinearCombinationBasis.r")
source("serial_assemble_x_and_v.r")
source("solve_for_interactions_on_x_and_v_by_others.r")
source("standard_basis.r")
source("find_maximums.r")
source("one_step_assemble_on_x_and_v.r")
source("eval_basis_functions.r")
source("uniform_basis.r")
source("uniform_basis_by_class.r")
source("find_pair_diff.r")
source("assemble_the_learning_matrix_on_x_and_v.r")
source("find_class_influence.r")
source("partition_x_and_v.r")
source("approximate_derivative_of_x_or_v.r")
source("scale_the_quantities.r")
source("get_Riemann_based_time_modifier.r")



learn_from_dynamics <- function(sys_info, obs_info, learn_info, obs_data,Riemann){
  
  one_block <- sys_info$d * sys_info$N
  x <- obs_data$traj  # this is a list of length M. We need to transform this into a 3-D array
  x <- array(unlist(x),dim=c(sys_info$d*sys_info$N,obs_info$L,obs_info$M))  # 3rd dim of list is each traj
  #x <- obs_data
  
  v <- c()
  xi <- c()
  extra_xi <- c()
  dot_xv <- c()
  Estimator_xi <- c()
  result <- learn_interactions_on_x_and_v(x, v, xi, dot_xv, obs_info$obs_time_vec, sys_info, learn_info,Riemann)
  
  Estimator <- result$Estimator
  extra_xv <- result$extra
  Estimator["phiXihat"]     = list(NULL)
  Estimator["Xibasis"]      = list(NULL)
  Estimator["emp_err_xi"]   = list(NULL)
  

  
  learn_out <- list()
  learn_out[["extra_xv"]] <- extra_xv
  learn_out[["rhoLTemp"]] <- extra_xv$rhoLTemp
  learn_out[["Estimator"]] <- Estimator
  learn_out["extra_xi"]     = list(NULL)
  
  return(learn_out)
}