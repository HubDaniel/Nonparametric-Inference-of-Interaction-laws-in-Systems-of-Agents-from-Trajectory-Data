

# function with sigma^2 solved analytically
likelihoodFn <- function(param,ns=30){
 
  # par <- exp(param)
  # Gamma <- par[1]
  # Eta <- par[2]

  
  Gamma <- param[1]    # we used gamma^2 below, no need to constraint it to be > 0, it's already > 0
  Eta <- exp(param[2]) # it is sigma_0^2 / sigma^2, need to control it to be > 0

  # construct psi(d) matrix (N*n) for training data
  psimat_exp <- vector("list",obs_info$M)   # list of length LMN, each is a matrix of psi(d_ijk)
  for(i in 1:length(psimat_exp)){
    psimat_exp[[i]] <- vector("list",ncol(train[[i]]))
  }
  for(i in 1:length(psimat_exp)){
    for(j in 1:length(psimat_exp[[i]])){
      psimat_exp[[i]][[j]] <- vector("list",sys_info$N)
    }
  }
  
  for(i in 1:length(psimat_exp)){
    for(j in 1:length(psimat_exp[[i]])){
      for(k in 1:length(psimat_exp[[i]][[j]])){
        psimat_exp[[i]][[j]][[k]] <- matrix(NA,nrow = sys_info$N,ncol = ns)
        pdisttemp <- Dmat[[i]][[j]][k,]  # distance between agent k with other agents on traj i at time j: d_ijk in notes page: star 1
        for(l in 1:ns){  # do not evaluates the formula since Rstudio evaluates factorial(large_number)=Inf, we evaluates the ratio term by term
          if(l==1){
            psimat_exp[[i]][[j]][[k]][,l] <- exp(-(Gamma*pdisttemp)^2) * sqrt((2*Gamma^2)^(l-1)/factorial(l-1)) * pdisttemp^(l-1)
          }else{
            temp <- 1
            for(z in 2:l){
              temp <- sqrt( (2*Gamma^2) ) / sqrt( z-1 ) * pdisttemp * temp
            }
            psimat_exp[[i]][[j]][[k]][,l] <- exp(-(Gamma*pdisttemp)^2) * temp   # provide NaN if temp is Inf (i.e, gamma too large)
          }
        }
      }
    }
  }
  
  
  # construct \tilde_psi(d_ijk) for training data
  
  psitildemat_exp <- vector("list",obs_info$M)   # list of length MLN, each is a matrix of psi(d_ijk)
  for(i in 1:length(psitildemat_exp)){
    psitildemat_exp[[i]] <- vector("list",ncol(train[[i]]))
  }
  for(i in 1:length(psitildemat_exp)){
    for(j in 1:length(psitildemat_exp[[i]])){
      psitildemat_exp[[i]][[j]] <- vector("list",sys_info$N)
    }
  }
  
  
  for(i in 1:length(psitildemat_exp)){
    for(j in 1:length(psitildemat_exp[[i]])){
      for(k in 1:length(psitildemat_exp[[i]][[j]])){
        indx2 <- k*sys_info$d
        indx1 <- indx2 - sys_info$d + 1   # [indx1,indx2] = [1,d],[d+1,2d],[2d+1,3d],...
        psitildemat_exp[[i]][[j]][[k]] <- 1/sys_info$N *Pmat[[i]][[j]][indx1:indx2,] %*% psimat_exp[[i]][[j]][[k]]
      }
    }
  }
  
  
  # construct \tilde_psi matrix from \tilde_psi(d_ijk) matrices
  tilde_psi_exp <- c()
  
  for(i in 1:length(psitildemat_exp)){
    for(j in 1:length(psitildemat_exp[[i]])){
      for(k in 1:length(psitildemat_exp[[i]][[j]])){
        tilde_psi_exp <- rbind(tilde_psi_exp,psitildemat_exp[[i]][[j]][[k]])
      }
    }
  }
  
  
  # here, we use the square of reciprocal of eta. If etaMLE has too many digits, it may indicate that the MLE result
  # is not accurate.
  
  LF1 <- length(Vtrain) * log(Eta) + determinant(diag(ncol(tilde_psi_exp)) + 1/Eta * t(tilde_psi_exp) %*% tilde_psi_exp)$mod[1]
  LF2 <- length(Vtrain) * 
    log(1/Eta * t(Vtrain) %*% Vtrain - 
                                (1/Eta)^2 * (t(Vtrain) %*% tilde_psi_exp) %*% solve(diag(ncol(tilde_psi_exp))+1/Eta * t(tilde_psi_exp) %*% tilde_psi_exp) %*% (t(tilde_psi_exp) %*% Vtrain))
  LF <- LF1 + LF2
    
  
  
  
  return(LF)
  
}














