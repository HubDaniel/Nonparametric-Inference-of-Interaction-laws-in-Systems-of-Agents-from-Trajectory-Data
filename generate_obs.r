## generate observations
# input: N: number of agents
#        d: dimension of each agents
#        xini: initial conditions, a vector of length d*N, first d element corresponding to the first agent
#        timespan: time intervals we want to draw obs.
# output: matrix with dimension nrow*ncol where ncol=1+dN and nrow is specified in timespan
#         first d columns are d-dimensions for the first agents, etc.


library(deSolve)
source("find_pair_diff.r")
source("sqdist_mod.r")
source("find_phis_of_pdist.r")
source("rep_matrix.r")
source("generateICs.r")
source("PS_init_config.r")
source("LJ_init_config.r")
source("uniform_dist.r")
source("self_organized_dynamics.r")


generateObs <- function(sys_info, obs_info){
  
  obs_data <- list()
  
  sys_var_len <- sys_info$d * sys_info$N
  x_obs <- array(0,dim=c(sys_var_len,obs_info$L,obs_info$M))
  ICs <- generateICs(obs_info$M,sys_info)  # matrix of dim = dN*M
  
  traj <- vector("list",obs_info$M)
  
  for (m in 1:obs_info$M){
    traj[[m]] <- self_organized_dynamics(ICs[,m], sys_info, obs_info)   # we use obs_info$time_vec
  }
  
  obs_data$ICs <- ICs
  obs_data$traj <- traj
  

  return(obs_data)
}

# # test
#  rr <- generateObs(sys_info, obs_info)
# 
#  #rr$traj is a list of length 50. each is a trajectory of matrix 20 * 200. each col is a time stamp
# 
#  length(rr$traj) # 50 traj, stored in a list
#  dim(rr$traj[[1]]) # dim = 20 * 200
# 
#  rplot <- as.data.frame(t(rr$traj[[3]]))
#  plot(rplot$`1`,rplot$`2`,xlim=c(min(rplot),max(rplot)),ylim=c(min(rplot),max(rplot)),type="l")
#  lines(rplot$`3`,rplot$`4`)
#  lines(rplot$`5`,rplot$`6`)
#  lines(rplot$`7`,rplot$`8`)
#  lines(rplot$`9`,rplot$`10`)
#  lines(rplot$`11`,rplot$`12`)
#  lines(rplot$`13`,rplot$`14`)
#  lines(rplot$`15`,rplot$`16`)
#  lines(rplot$`17`,rplot$`18`)
#  lines(rplot$`19`,rplot$`20`,col="red")   # predator


 