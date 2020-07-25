# self_organized_dynamics

self_organized_dynamics <- function(y_init, sys_info, obs_info,true_kernel=T){
  
  #print("self_organized_dynamics")

  constrFun <- function(x,sys_info){         # N,d can also be passed through function odefun in argument parms
    
    d <- sys_info$d
    N <- sys_info$N
    
    # reshape x into a matrix d*N for later functions to use
    x <- matrix(x, ncol=N, nrow=d)
    
    # step 1: construct part 1 use find_pair_diff.R
    #print("step 1 : start find pair diff")
    pair_diff <- find_pair_diff(x)
    
    # step 2: construct part 2 use sqdist_mod.R
    #print("step 2 : start sqdist_mod")
    pair_dist <- sqdist_mod(x)
    
    # step 3: construct part 3 use find_phis_of_pdist.R
    #print("step 3 : find_phis_of_pdist")
    
    phi_r <- find_phis_of_pdist(pair_dist,sys_info,true_kernel)
    #print(phi_r)

    
    
    # step 4: sum up part 1 and part 3 use find_collective_change.R
    #print("step 4")
    
    
    phi_r_rep <- rep_matrix(phi_r,d)             # extend phi_r to be dN*N matrix, matches the dim of pair_diff
    rhs <- rowSums(phi_r_rep * pair_diff)
    
    # step 5: return as a list as required in odefun
    #print("step 5")
    result <- list(rhs)
    return(result)
  }
  
  
  odefun <- function(t,x,parms){   # x is a d*N vector since we are solving d*N coupled ODE's
    return(constrFun(x,sys_info))
  }
  
  print("start ode")
  try(out <- ode(times = obs_info$time_vec, y = y_init, func = odefun, parms = NULL,method = "ode45"))
  # out is a matrix with first column the time_vec. column "1" and "2" are coordinates for agent 1, etc.
  print("finish ode")
  #print(dim(out))
  #print(head(out))

  
  # only extract those time points that we used to learn the dynamic (similar to observed_dynamics.m in matlab)
  indx <- out[,1] %in% obs_info$obs_time_vec
  out <- out[indx,]
  
  #print(dim(out))
  
  # reshape the data for later use
  # reshape to a matrix with 20*200, each col is a time stamp
  out1 <- t(out[,-1])
  
  # any noise?
  if (sys_info$noise_on_obs){
    if(sys_info$noise_type=="add"){
      noise_matrix <- sys_info$noise_sigma * matrix(runif(length(out1),-0.1,0.1),nrow = nrow(out1),ncol = ncol(out1))
      out1 <- out1 + noise_matrix
    }else{
      noise_matrix <- 1+(sys_info$noise_sigma * matrix(runif(length(out1),-0.1,0.1),nrow = nrow(out1),ncol = ncol(out1)))
      out1 <- out1 * noise_matrix
    }
  }
  
  #print(dim(out1))
  #print(head(out1))
  
  return(out1)
}


# test
#ICs <- generateICs(obs_info$M,sys_info)
#r <- self_organized_dynamics(ICs[,1], sys_info, obs_info)



