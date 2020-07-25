# compute trajectory error

# computeTrajErr

source("self_organized_dynamics.r")
source("find_phis_of_pdist.r")
source("traj_norm.r")
source("state_norm.r")


computeTrajErr <- function(sys_info, obs_info, y_init){
  # compute trajectory error for both the learning period and future period
  
  # idea: use the phihat estimated in the learning process
  # 1. generate a traj with true kernel and initial condition 
  # (this initial condition is the same as in above learning process)
  # 2. generate a traj with estimated kernel (phihat) and the same initial condition
  # 3. error is the maximum of norm of their diff.
  
  
  # step 1: generate true traj
  
  print("step 1")
  M <- dim(y_init)[2]
  sys_var_len                 <- sys_info$d * sys_info$N
  trajt                       <- array(0,dim = c(sys_var_len, obs_info$L, M))
  trajtfut                    <- array(0,dim = c(sys_var_len, obs_info$L, M))
  
  for (m in 1:M){
    print(paste0("m=",m,"/",M))
    
    out.T_L <- self_organized_dynamics(y_init[,m], sys_info, obs_info)   # obs generated from T_1 to T_L
    # matrix of dN*M
    
    temp <- obs_info$time_vec
    temp2 <- obs_info$obs_time_vec
    obs_info$time_vec <- c(0,seq(from = obs_info$T_L, to = sys_info$T_f, length.out = obs_info$L)) #add time 0 to satisfy initial condition
    obs_info$obs_time_vec <- seq(from = obs_info$T_L, to = sys_info$T_f, length.out = obs_info$L)
    out.T_F <- self_organized_dynamics(y_init[,m], sys_info, obs_info)  # obs generated from T_L to T_F
    # matrix of dN*M
    obs_info$time_vec <- temp
    obs_info$obs_time_vec <- temp2

    trajt[,,m] <- out.T_L    # obs from 0 to T_L from true kernel
    trajtfut[,,m] <- out.T_F # obs from T_L to T_F from true kernel
    
  }
  
  
  
  # step 2: generate estimated traj
  print("step 2: generate estimated traj")
  
  trajhat                        <- array(0,dim = c(sys_var_len, obs_info$L, M))
  trajfuthat                     <- array(0,dim = c(sys_var_len, obs_info$L, M))
  
  
  tempkernel <- sys_info$phiE
  sys_info$phiE <- learn$Estimator$phiEhat
  
  #print("kernel")
  #print(sys_info$phiE)
  
  flaghat <- c()
  flaghatw <- c()
  flagfuthat <- c()
  flagfuthatw <- c()
  
  for (m in 1:M){
    print(paste0("m=",m,"/",M))
    
    print("start self org dyn for hat")
    
    outhat.T_L <- self_organized_dynamics(y_init[,m], sys_info, obs_info,true_kernel=F)   # obs generated from T_1 to T_L
    # matrix of dN*L
    print("finish self org dyn for hat")
    
    
    temp <- obs_info$time_vec
    temp2 <- obs_info$obs_time_vec
    obs_info$time_vec <- c(0,seq(from = obs_info$T_L, to = sys_info$T_f, length.out = obs_info$L))
    obs_info$obs_time_vec <- seq(from = obs_info$T_L, to = sys_info$T_f, length.out = obs_info$L)
    
    print("start self org dyn for future hat")
    outhat.T_F <- self_organized_dynamics(y_init[,m], sys_info, obs_info,true_kernel=F)  # obs generated from T_L to T_F
    print("finish self org dyn for future hat")
    
    
    # matrix of dN*L
    obs_info$time_vec <- temp
    obs_info$obs_time_vec <- temp2
    
    #trajhat[,,m] <- outhat.T_L
    #trajfuthat[,,m] <- outhat.T_F
    #try(trajhat[,,m] <- outhat.T_L)    # obs from 0 to T_L from estimated kernel
    #try(trajfuthat[,,m] <- outhat.T_F) # obs from T_L to T_F from estimated kernel
    
    
    tryCatch(trajhat[,,m] <- outhat.T_L,
             error = function(e){
               flaghat <<- c(flaghat,m)
             },
             warning = function(w){
               flaghatw <<- c(flaghatw,m)
             })


    tryCatch(trajfuthat[,,m] <- outhat.T_F,
             error = function(e){
               flagfuthat <<- c(flagfuthat,m)
             },
             warning = function(w){
               flagfuthatw <<- c(flagfuthatw,m)
             })
    
    print(paste0("flaghat = ",flaghat))
    print(paste0("flaghatw = ",flaghatw))
    print(paste0("flagfuthat = ",flagfuthat))
    print(paste0("flagfuthatw = ",flagfuthatw))
    
    
    
  }
  
  sys_info$phiE <- tempkernel
  
  
  # step 3: Remove flagged traj
  print("step 3")
  allFlags <- unique(c(flaghat,flaghatw,flagfuthat,flagfuthatw))
  
  print("allFlags")
  print(allFlags)
  
  if(!is.null(allFlags)){
    trajt <- trajt[,,-allFlags]
    trajtfut <- trajtfut[,,-allFlags]
    trajhat <- trajhat[,,-allFlags]
    trajfuthat <- trajfuthat[,,-allFlags]
  }


  # step 4: Compute trajectory error
  print("step 4")
  sup_err                     <- rep(0,dim(trajt)[3])
  sup_err_fut                 <- rep(0,dim(trajt)[3])
  for (m in 1:dim(trajt)[3]){
    print(paste0("m=",m,"/",dim(trajt)[3]))
    
    sup_err[m] <- traj_norm(trajt[,,m], trajhat[,,m],sys_info)
    sup_err_fut[m] <- traj_norm(trajtfut[,,m], trajfuthat[,,m],sys_info)
  }
  
  result <- list()
  result[["sup_err"]] <- sup_err
  result[["sup_err_fut"]] <- sup_err_fut
  result$trajhat <- trajhat
  result$trajfuthat <- trajfuthat
  result$trajt <- trajt
  result$trajtfut <- trajtfut
  result$flag <- allFlags
  
  return(result)
  
  
  
}