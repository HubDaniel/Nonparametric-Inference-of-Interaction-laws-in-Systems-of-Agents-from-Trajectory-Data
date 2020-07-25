# relativeErrorInfluenceFunction

source("L2_rho_T_energy.r")

relativeErrorInfluenceFunction <- function(phihat, phi, sys_info, basis, rhoLT){
  
  print("start rEIF")
  
  print("defining phidiff")
  
  phidiff <- vector("list",sys_info$K)
  
  for(i in 1:length(phidiff)){
    phidiff[[i]] <- vector("list",sys_info$K)
  }
  
  
  # function for force evaluation to prevent lazy evaluation
  
  h <- function(ind1,ind2){
    force(ind1)
    force(ind2)
    function(r) {phi[[ind1]][[ind2]](r) - phihat[[ind1]][[ind2]](r)$y} 
  }
  
  
  
  for(k_1 in seq(from=sys_info$K,to=1,by=-1)){
    for(k_2 in seq(from=sys_info$K,to=1,by=-1)){
      phidiff[[k_1]][[k_2]] <- h(k_1,k_2)
    }
  }

  

  
 
  
  print("defining L2_rho_T")
  
  L2_rho_T <- function(a,b,c,d) {L2_rho_T_energy(a,b,c,d)}
  
  print("start phi_L2norms")
  
  phi_L2norms              = L2_rho_T( phi, sys_info, rhoLT, basis )     # matrix of k*K
  
  print("start  phihat_L2rhoTdiff")
  phihat_L2rhoTdiff        = L2_rho_T( phidiff, sys_info, rhoLT, basis ) # matrix of k*K
  
  
  
  
  
  Err.Abs <- Err.Rel <- matrix(0,nrow=sys_info$K,ncol=sys_info$K)
  
  for(k_1 in 1:sys_info$K){
    for(k_2 in 1:sys_info$K){
      Err.Abs[k_1,k_2] <- phihat_L2rhoTdiff[k_1,k_2]
      if(phi_L2norms[k_1,k_2]!=0){
        Err.Rel[k_1,k_2] <- phihat_L2rhoTdiff[k_1,k_2]/ phi_L2norms[k_1,k_2]
      }else{
        if(phihat_L2rhoTdiff[k_1,k_2] == 0){
          Err.Rel[k_1,k_2] <- 0
        }else{
          Err.Rel[k_1,k_2] <- Inf
        }
      }
    }
  }
  
  
  
  Err <- list()
  Err[["Err.Abs"]] <- Err.Abs
  Err[["Err.Rel"]] <- Err.Rel
  
  return(Err)
  
}