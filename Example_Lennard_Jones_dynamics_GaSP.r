# Lennard_Jones_dynamics with GaSP

setwd("C:/Users/x_r_m/Desktop/Learn/RA/V2")

########################################################################################
########################################################################################
########################################################################################

sys_info                <- list()
sys_info$d              <- 2
sys_info$N              <- 7
sys_info$K              <- 1
sys_info$type_info      <- rep(1,sys_info$N)  
sys_info$kappa          <- rep(1,sys_info$K)
sys_info$T_f            <- 0.05
sys_info$mu0            <- function() LJ_init_config(sys_info$d, sys_info$N, 2)
sys_info$phiE           <- list()
sys_info$phiE[[1]]      <- list()
# sys_info$phiE[[1]][[1]] <- function(r,epsilon = 10,sigma = 1){
#   result <- rep(0,length(r))     # result holder
#   ind <- r!=0
#   rinv <- sigma/r[ind]
#   result[ind] <- 24 * epsilon/(sigma^2) * (-2*(rinv^14) + (rinv^8))
# 
#   ind <- r==0
#   result[ind] <- -Inf
# 
#   return(result)
# }


sys_info$phiE[[1]][[1]] <- function(r) {
                              sin(r)
}


obs_info                <- list()
obs_info$L              <- 50
obs_info$M              <- 3
obs_info$M_rhoT         <- 2000
obs_info$T_L            <- 0.01
obs_info$time_vec       <- c(0,seq(from = 1e-3, to = obs_info$T_L, length.out = obs_info$L))
obs_info$obs_time_vec   <- seq(from = 1e-3, to = obs_info$T_L, length.out = obs_info$L)   # time vec for obs.

obs_info$hist_num_bins  <- 10000

basis_info <- list()
basis_info$n            <- matrix(2*ceiling(150*(obs_info$M/log(obs_info$M))^(1/5)),nrow = 1,ncol = 1)
basis_info$type         <- 'standard'                                                             
basis_info$degree       <- matrix(1,ncol=1,nrow=1)

learn_info              <- list()
learn_info$Ebasis_info  <- basis_info
learn_info$solver_type  <- "svd"

# following part is for mimic the result of PNAS paper. out model need to add noise directly on V not on X (here is on X)
sys_info$noise_on_obs <- F                   
sys_info$noise_dist <- "uniform"
sys_info$noise_type <- "add"  # or mult
sys_info$noise_sigma <- 1     # also used as the standard deviation of our noise sigma_0

########################################################################################
########################################################################################
########################################################################################

# load packages used below
source("generate_obs.r")
source("learn_from_dynamics.r")
library(beepr)

#########################################################################################
########              Section 0 - Prepare the Data              #########
#########################################################################################

# generate data
obs <- generateObs(sys_info, obs_info)  # list of length 2, ICs: init conditions; traj: list of 1000 traj

# separate into train and test
p <- 0.7  # proportion of train
train <- obs$traj
for(i in 1:length(train)){
  train[[i]] <- obs$traj[[i]][,1:round(ncol(obs$traj[[i]])*p)]
}

test <- obs$traj
for(i in 1:length(test)){
  test[[i]] <- obs$traj[[i]][,(ncol(train[[i]])+1):ncol(obs$traj[[i]])]
}    # train and test are both a list of length M, each is a matrix of dimension dN*(time_length)

# Pairwise distance for training data
Dmat <- vector("list",obs_info$M)
for(i in 1:length(Dmat)){
  Dmat[[i]] <- vector("list",ncol(train[[i]])) 
}

for(i in 1:length(train)){
  for(j in 1:ncol(train[[i]])){
    D <- train[[i]][,j]
    Dmattemp <- matrix(D,nrow=sys_info$d)
    pdist <- sqdist_mod(Dmattemp) # matrix of N*N containing pairwise distances
    Dmat[[i]][[j]] <- pdist  
  }
}

# Pairwise distance for testing data
Dmat2 <- vector("list",obs_info$M)
for(i in 1:length(Dmat2)){
  Dmat2[[i]] <- vector("list",ncol(test[[i]])) 
}

for(i in 1:length(test)){
  for(j in 1:ncol(test[[i]])){
    D <- test[[i]][,j]
    Dmattemp <- matrix(D,nrow=sys_info$d)
    pdist <- sqdist_mod(Dmattemp) # matrix of N*N containing pairwise distances
    Dmat2[[i]][[j]] <- pdist  
  }
}

# Pairwise difference for training data
Pmat <- vector("list",obs_info$M)
for(i in 1:length(Pmat)){
  Pmat[[i]] <- vector("list",ncol(train[[i]])) 
}

for(i in 1:length(train)){
  for(j in 1:ncol(train[[i]])){
    P <- train[[i]][,j]
    Pmattemp <- matrix(P,nrow=sys_info$d)
    pdiff <- find_pair_diff(Pmattemp) # matrix of dN*N containing pairwise diff
    Pmat[[i]][[j]] <- pdiff  
  }
}

# Pairwise difference for testing data
Pmat2 <- vector("list",obs_info$M)
for(i in 1:length(Pmat2)){
  Pmat2[[i]] <- vector("list",ncol(test[[i]])) 
}

for(i in 1:length(test)){
  for(j in 1:ncol(test[[i]])){
    P <- test[[i]][,j]
    Pmattemp <- matrix(P,nrow=sys_info$d)
    pdiff <- find_pair_diff(Pmattemp) # matrix of dN*N containing pairwise diff
    Pmat2[[i]][[j]] <- pdiff  
  }
}

# True V's (with noise)
V <- vector("list",length(obs$traj))
for(i in 1:length(V)){
  # one traj at a time
  D <- obs$traj[[i]]  # matrix of dN*L
  
  # the following is our observed time stamps
  obs_info$obs_time_vec   <- seq(from = 1e-3, to = obs_info$T_L, length.out = obs_info$L)   # time vec for obs.
  
  # V is calculated using location X since we do not have V.
  # start from 1e-3, denote this as time t1, we use v1=(X2-X1)/(t2-t1)
  
  D1 <- D
  D1 <- cbind(D1,D1[,obs_info$L])
  Vtemp <- (D1[,2:(obs_info$L+1)]-D1[,1:obs_info$L])/(obs_info$obs_time_vec[2]-obs_info$obs_time_vec[1])
  # last col is 0 since we have no info on the later location
  # reshape V to be a long vector of length dNL
  #V1 <- as.vector(Vtemp) # in the form: (V_1,V_2,...,V_N) at t1, (V_1,V_2,...,V_N) at t2, ... , (V_1,V_2,...,V_N) at tL
  V[[i]] <- Vtemp
}

# add noise to V
for(i in 1:length(V)){
  V[[i]] <- V[[i]] + matrix(rnorm(dim(V[[i]])[1]*dim(V[[i]])[2],0,sys_info$noise_sigma),nrow = dim(V[[i]])[1],ncol = dim(V[[i]])[2])
}

Vtrain <- V
for(i in 1:length(Vtrain)){
  Vtrain[[i]] <- Vtrain[[i]][,1:round(ncol(Vtrain[[i]])*p)]
}

Vtest <- V
for(i in 1:length(Vtest)){
  Vtest[[i]] <- Vtest[[i]][,(ncol(Vtrain[[i]])+1):(ncol(Vtest[[i]])-1)]
}
# since we use difference to calculate the V, the last column of V is zero, not accurate, remove it

Vtrain <- unlist(Vtrain) # transform to a long vector to be used later



########################################################################################
########################################################################################
########################################################################################


#########################################################################################
########              Section 1 - Hyperparameter Estimation              #########
#########################################################################################


# use the formula in section 7 in notes
# find function Q w.r.t.   eta   and   tilde(psi)

source("Example_Lennard_Jones_dynamics_GaSP-LF.r")   # this contains the likelihood function we're going to min

ns=5
est_result <- optim(par = c(1,0), fn = likelihoodFn,ns=ns )  # ns: number of basis functions passed into likelihoodFn

gamma2MLE <- (est_result$par[1])^2
etaMLE <- exp(est_result$par[2])
# 1st one is gamma, 2nd one is eta

# sigma^2
# recalculate tilde_psi_exp to calculate sigmaMLE
psimat_exp_MLE <- vector("list",obs_info$M)   # list of length LMN, each is a matrix of psi(d_ijk)
for(i in 1:length(psimat_exp_MLE)){
  psimat_exp_MLE[[i]] <- vector("list",ncol(train[[i]]))
}
for(i in 1:length(psimat_exp_MLE)){
  for(j in 1:length(psimat_exp_MLE[[i]])){
    psimat_exp_MLE[[i]][[j]] <- vector("list",sys_info$N)
  }
}

for(i in 1:length(psimat_exp_MLE)){
  for(j in 1:length(psimat_exp_MLE[[i]])){
    for(k in 1:length(psimat_exp_MLE[[i]][[j]])){
      psimat_exp_MLE[[i]][[j]][[k]] <- matrix(NA,nrow = sys_info$N,ncol = ns)
      pdisttemp <- Dmat[[i]][[j]][k,]  # distance between agent k with other agents on traj i at time j: d_ijk in notes page: star 1
      for(l in 1:ns){  # do not evaluates the formula since Rstudio evaluates factorial(large_number)=Inf, we evaluates the ratio term by term
        if(l==1){
          psimat_exp_MLE[[i]][[j]][[k]][,l] <- exp(-(gamma2MLE*pdisttemp)^2) * sqrt((2*gamma2MLE^2)^(l-1)/factorial(l-1)) * pdisttemp^(l-1)
        }else{
          temp <- 1
          for(z in 2:l){
            temp <- sqrt( (2*gamma2MLE^2) ) / sqrt( z-1 ) * pdisttemp * temp
          }
          psimat_exp_MLE[[i]][[j]][[k]][,l] <- exp(-(gamma2MLE*pdisttemp)^2) * temp
        }
      }
    }
  }
}


# construct \tilde_psi(d_ijk) for training data

psitildemat_exp_MLE <- vector("list",obs_info$M)   # list of length MLN, each is a matrix of psi(d_ijk)
for(i in 1:length(psitildemat_exp_MLE)){
  psitildemat_exp_MLE[[i]] <- vector("list",ncol(train[[i]]))
}
for(i in 1:length(psitildemat_exp_MLE)){
  for(j in 1:length(psitildemat_exp_MLE[[i]])){
    psitildemat_exp_MLE[[i]][[j]] <- vector("list",sys_info$N)
  }
}


for(i in 1:length(psitildemat_exp_MLE)){
  for(j in 1:length(psitildemat_exp_MLE[[i]])){
    for(k in 1:length(psitildemat_exp_MLE[[i]][[j]])){
      indx2 <- k*sys_info$d
      indx1 <- indx2 - sys_info$d + 1   # [indx1,indx2] = [1,d],[d+1,2d],[2d+1,3d],...
      psitildemat_exp_MLE[[i]][[j]][[k]] <- 1/sys_info$N *Pmat[[i]][[j]][indx1:indx2,] %*% psimat_exp_MLE[[i]][[j]][[k]]
    }
  }
}


# construct \tilde_psi matrix from \tilde_psi(d_ijk) matrices
tilde_psi_exp_MLE <- c()

for(i in 1:length(psitildemat_exp_MLE)){
  for(j in 1:length(psitildemat_exp_MLE[[i]])){
    for(k in 1:length(psitildemat_exp_MLE[[i]][[j]])){
      tilde_psi_exp_MLE <- rbind(tilde_psi_exp_MLE,psitildemat_exp_MLE[[i]][[j]][[k]])
    }
  }
}


longterm <- (1/etaMLE) * diag(nrow(tilde_psi_exp_MLE)) - 
  (1/etaMLE)^2 * tilde_psi_exp_MLE %*% solve(diag(ncol(tilde_psi_exp_MLE)) + 
                                        (1/etaMLE) * t(tilde_psi_exp_MLE) %*% tilde_psi_exp_MLE) %*% t(tilde_psi_exp_MLE)

sigma2 <- (1/(sys_info$d*sys_info$N*obs_info$L*obs_info$M)) * t(Vtrain) %*% longterm %*% Vtrain   
sigma_02 <- etaMLE*sigma2  # should be close to 1 (set it to be 1)  

est_result
gamma2MLE
etaMLE
sigma_02
sigma2
est_result$value




########################################################################################
########################################################################################
########################################################################################


#########################################################################################
########              Section 2 - Piecewise Linear Basis Functions              #########
#########################################################################################

# assuming sigma_0 and sigma are both 1; basis functions are used as below.


# interval we are working on [0,Rmax].
Rmax <- max(unlist(Dmat))
interval <- seq(from=0,to=Rmax,length.out = basis_info$n[[1]][[1]]+1)  # use piecewise constant so num of inter is n

# construct \psi matrix of dim N*n
psimat <- vector("list",obs_info$M)   # list of length LMN, each is a matrix of psi(d_ijk)
for(i in 1:length(psimat)){
  psimat[[i]] <- vector("list",ncol(train[[i]]))
}
for(i in 1:length(psimat)){
  for(j in 1:length(psimat[[i]])){
    psimat[[i]][[j]] <- vector("list",sys_info$N)
  }
}

for(i in 1:length(psimat)){
  for(j in 1:length(psimat[[i]])){
    for(k in 1:length(psimat[[i]][[j]])){
      psimat[[i]][[j]][[k]] <- matrix(0,nrow = sys_info$N,ncol = basis_info$n[[1]][[1]])
      pdisttemp <- Dmat[[i]][[j]][k,]  # distance between agent k with other agents on traj i at time j: d_ijk in notes page: star 1
      for(l in 1:length(pdisttemp)){
        bin <- histc(pdisttemp[l],interval)$bin
        if(bin>basis_info$n){
          bin <- bin - 1 
        }
        psimat[[i]][[j]][[k]][l,bin] <- 1
      }
    }
  }
}

# as in notes, check \psi * \psi^T has identity diagonal (diag = 1)
bigpsi <- c()
for(i in 1:length(psimat)){
  for(j in 1:length(psimat[[i]])){
    for(k in 1:length(psimat[[i]][[j]])){
      bigpsi <- rbind(bigpsi,psimat[[i]][[j]][[k]])
    }
  }
}

diagonal_bigpsi <- (diag(bigpsi %*% t(bigpsi)))
table(diagonal_bigpsi)



# construct psi tilde

psitildemat <- vector("list",obs_info$M)   # list of length LMN, each is a matrix of psi(d_ijk)
for(i in 1:length(psitildemat)){
  psitildemat[[i]] <- vector("list",ncol(train[[i]]))
}
for(i in 1:length(psitildemat)){
  for(j in 1:length(psitildemat[[i]])){
    psitildemat[[i]][[j]] <- vector("list",sys_info$N)
  }
}


for(i in 1:length(psitildemat)){
  for(j in 1:length(psitildemat[[i]])){
    for(k in 1:length(psitildemat[[i]][[j]])){
      indx2 <- k*sys_info$d
      indx1 <- indx2 - sys_info$d + 1
      psitildemat[[i]][[j]][[k]] <- 1/sys_info$N *Pmat[[i]][[j]][indx1:indx2,] %*% psimat[[i]][[j]][[k]]
    }
  }
}

# construct \tilde_psi matrix from \tilde_psi(d_ijk) matrices
tilde_psi <- c()

for(i in 1:length(psitildemat)){
  for(j in 1:length(psitildemat[[i]])){
    for(k in 1:length(psitildemat[[i]][[j]])){
      tilde_psi <- rbind(tilde_psi,psitildemat[[i]][[j]][[k]])
    }
  }
}


# construct \tilde_psi_* matrix
# our prediction model is to predict on N agents on the same trajectory at the same time.


# 1. find each \psi(d_ijk*), exactly the same approach as above but with test data


# interval we are working on [0,Rmax2].
Rmax2 <- max(unlist(Dmat2))
interval2 <- seq(from=0,to=Rmax2,length.out = basis_info$n[[1]][[1]]+1)  # use piecewise constant so num of inter is n

# construct \psi matrix of dim N*n
psimat2 <- vector("list",obs_info$M)   # list of length LMN, each is a matrix of psi(d_ijk)
for(i in 1:length(psimat2)){
  psimat2[[i]] <- vector("list",ncol(Vtest[[i]]))
}
for(i in 1:length(psimat2)){
  for(j in 1:length(psimat2[[i]])){
    psimat2[[i]][[j]] <- vector("list",sys_info$N)
  }
}

for(i in 1:length(psimat2)){
  for(j in 1:length(psimat2[[i]])){
    for(k in 1:length(psimat2[[i]][[j]])){
      psimat2[[i]][[j]][[k]] <- matrix(0,nrow = sys_info$N,ncol = basis_info$n[[1]][[1]])
      pdisttemp <- Dmat2[[i]][[j]][k,]  # distance between agent k with other agents on traj i at time j: d_ijk in notes page: star 1
      for(l in 1:length(pdisttemp)){
        bin <- histc(pdisttemp[l],interval2)$bin
        if(bin>basis_info$n){
          bin <- bin - 1 
        }
        psimat2[[i]][[j]][[k]][l,bin] <- 1
      }
    }
  }
}


# construct psi tilde
psitildemat2 <- vector("list",obs_info$M)   # list of length LMN, each is a matrix of psi(d_ijk)
for(i in 1:length(psitildemat2)){
  psitildemat2[[i]] <- vector("list",ncol(Vtest[[i]]))
}
for(i in 1:length(psitildemat2)){
  for(j in 1:length(psitildemat2[[i]])){
    psitildemat2[[i]][[j]] <- vector("list",sys_info$N)
  }
}


for(i in 1:length(psitildemat2)){
  for(j in 1:length(psitildemat2[[i]])){
    for(k in 1:length(psitildemat2[[i]][[j]])){
      indx2 <- k*sys_info$d
      indx1 <- indx2 - sys_info$d + 1
      psitildemat2[[i]][[j]][[k]] <- 1/sys_info$N *Pmat2[[i]][[j]][indx1:indx2,] %*% psimat2[[i]][[j]][[k]]
    }
  }
}

# psitildemat2[[i]][[j]][[k]] contains \tilde\psi(d_ijk*)

tilde_psi_star <- vector("list",obs_info$M)
for(i in 1:length(tilde_psi_star)){
  tilde_psi_star[[i]] <- vector("list",ncol(Vtest[[i]]))
}

for(i in 1:length(psitildemat2)){
  for(j in 1:length(psitildemat2[[i]])){
    for(k in 1:length(psitildemat2[[i]][[j]])){
      tilde_psi_star[[i]][[j]] <- rbind(tilde_psi_star[[i]][[j]],psitildemat2[[i]][[j]][[k]])
    }
  }
}



## now we can make predictions use the results above

###
# now let's predict all test cases
# all cases means we need to use all matrices in the list tilde_psi_star[[i]][[j]] 
###
ninside <- 0
for(i in 1:length(tilde_psi_star)){
  for(j in 1:length(tilde_psi_star[[i]])){
    
    print(paste0("i = ",i,"; j = ",j,"  (total: i = ",length(tilde_psi_star)," j = ",length(tilde_psi_star[[i]]),")"))
    
    
    tilde_psi_star_temp <- tilde_psi_star[[i]][[j]]
    # real value
    real_V <- Vtest[[i]][,j]
    
    # predictive mean
    midterm <- diag(nrow(tilde_psi))-tilde_psi %*% solve(diag(ncol(tilde_psi))+t(tilde_psi)%*%tilde_psi) %*% t(tilde_psi)
    exp <- tilde_psi_star_temp %*% t(tilde_psi) %*% midterm %*% Vtrain
    pmean <- as.vector(exp)
    
    # predictive covariance
    pcov <- tilde_psi_star_temp %*% t(tilde_psi_star_temp) + diag(nrow(tilde_psi_star_temp)) -
      tilde_psi_star_temp %*% t(tilde_psi) %*% tilde_psi %*% t(tilde_psi_star_temp) + 
      tilde_psi_star_temp %*% t(tilde_psi) %*% tilde_psi %*% solve(diag(ncol(tilde_psi))+t(tilde_psi)%*%tilde_psi) %*% t(tilde_psi_star_temp %*% t(tilde_psi) %*% tilde_psi)
    pvar <- diag(pcov)
    
    # CI
    lowers <- pmean-(pvar)^(0.5)*qnorm(p = 0.975)
    uppers <- pmean+(pvar)^(0.5)*qnorm(p = 0.975)
    CI <- cbind(lowers,pmean,uppers,real_V)
    
    # calculate how many predictions landed inside the CI
    ninside <- ninside + sum( (CI[,"real_V"] > CI[,"lowers"]) & (CI[,"real_V"] < CI[,"uppers"]) )
    
  }
}

# percentage of inside CI
ntotal <- length(Vtest) * length(Vtest[[1]])
ntotal
ninside
ninside/ntotal     # 773/1218 = 63.46%





########################################################################################
########################################################################################
########################################################################################

#########################################################################################
########              Section 3 - SE Covariance Functions              #########
#########################################################################################

# now try to use the kernel on page 11 on https://www.csie.ntu.edu.tw/~cjlin/talks/kuleuven_svm.pdf
# use approximation of this kernel

# The only difference is how to construct \psi matrix, rest are the same.
coverage <- c()
avelength <- c()
RMSE <- c()
RMSE.sdPmean <- c()

nsv <- c(15,30,50)
betav <- c(0.8,1,3,5,6,7) # it's the gamma

for(ii in nsv){
  for(jj in betav){
    
    print(paste0("ns = ",ii,"; beta = ",jj))
    
    ns <- ii
    beta <- jj
    
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
              psimat_exp[[i]][[j]][[k]][,l] <- exp(-(beta*pdisttemp)^2) * sqrt((2*beta^2)^(l-1)/factorial(l-1)) * pdisttemp^(l-1)
            }else{
              temp <- 1
              for(z in 2:l){
                temp <- sqrt( (2*beta^2) ) / sqrt( z-1 ) * pdisttemp * temp
              }
              psimat_exp[[i]][[j]][[k]][,l] <- exp(-(beta*pdisttemp)^2) * temp
            }
          }
        }
      }
    }
    
    
    
    # # as in notes, check \psi * \psi^T has identity diagonal (diag = 1). want to see ratio_diag_is_1 = 1 as ns-->Inf
    # bigpsi_exp <- c()
    # for(i in 1:length(psimat_exp)){
    #   for(j in 1:length(psimat_exp[[i]])){
    #     for(k in 1:length(psimat_exp[[i]][[j]])){
    #       bigpsi_exp <- rbind(bigpsi_exp,psimat_exp[[i]][[j]][[k]])
    #     }
    #   }
    # }
    # 
    # diagonal_bigpsi_exp <- (diag(bigpsi_exp %*% t(bigpsi_exp)))
    # ratio_diag_is_1 <- mean(diagonal_bigpsi_exp == 1)
    
    
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
    
    
    ###############################################################################################################
    # construct \tilde_psi_* matrix
    # our prediction model is to predict on N agents on the same trajectory at the same time.
    
    
    # 1. find each \psi(d_ijk*), exactly the same approach as above but with test data
    
    
    psimat_exp2 <- vector("list",obs_info$M)   # list of length MLN, each is a matrix of psi(d_ijk*)
    for(i in 1:length(psimat_exp2)){
      psimat_exp2[[i]] <- vector("list",ncol(Vtest[[i]]))
    }
    for(i in 1:length(psimat_exp2)){
      for(j in 1:length(psimat_exp2[[i]])){
        psimat_exp2[[i]][[j]] <- vector("list",sys_info$N)
      }
    }
    
    for(i in 1:length(psimat_exp2)){
      for(j in 1:length(psimat_exp2[[i]])){
        for(k in 1:length(psimat_exp2[[i]][[j]])){
          psimat_exp2[[i]][[j]][[k]] <- matrix(NA,nrow = sys_info$N,ncol = ns)
          pdisttemp <- Dmat2[[i]][[j]][k,]  # distance between agent k with other agents on traj i at time j: d_ijk in notes page: star 1
          for(l in 1:ns){
            if(l==1){
              psimat_exp2[[i]][[j]][[k]][,l] <- exp(-(beta*pdisttemp)^2) * sqrt((2*beta^2)^(l-1)/factorial(l-1)) * pdisttemp^(l-1)
            }else{
              temp <- 1
              for(z in 2:l){
                temp <- sqrt( (2*beta^2) ) / sqrt( z-1 ) * pdisttemp * temp
              }
              psimat_exp2[[i]][[j]][[k]][,l] <- exp(-(beta*pdisttemp)^2) * temp
            }
          }
        }
      }
    }
    
   
    # construct \tilde\psi(d_ijk*)
    
    psitildemat_exp2 <- vector("list",obs_info$M)   # list of length LMN, each is a matrix of psi(d_ijk)
    for(i in 1:length(psitildemat_exp2)){
      psitildemat_exp2[[i]] <- vector("list",ncol(Vtest[[i]]))
    }
    for(i in 1:length(psitildemat_exp2)){
      for(j in 1:length(psitildemat_exp2[[i]])){
        psitildemat_exp2[[i]][[j]] <- vector("list",sys_info$N)
      }
    }
    
    
    for(i in 1:length(psitildemat_exp2)){
      for(j in 1:length(psitildemat_exp2[[i]])){
        for(k in 1:length(psitildemat_exp2[[i]][[j]])){
          indx2 <- k*sys_info$d
          indx1 <- indx2 - sys_info$d + 1
          psitildemat_exp2[[i]][[j]][[k]] <- 1/sys_info$N *Pmat2[[i]][[j]][indx1:indx2,] %*% psimat_exp2[[i]][[j]][[k]]
        }
      }
    }
    
    # psitildemat_exp2[[i]][[j]][[k]] contains \tilde\psi(d_ijk*)
    
    tilde_psi_star_exp <- vector("list",obs_info$M)
    for(i in 1:length(tilde_psi_star_exp)){
      tilde_psi_star_exp[[i]] <- vector("list",ncol(Vtest[[i]]))
    }
    
    for(i in 1:length(psitildemat_exp2)){
      for(j in 1:length(psitildemat_exp2[[i]])){
        for(k in 1:length(psitildemat_exp2[[i]][[j]])){
          tilde_psi_star_exp[[i]][[j]] <- rbind(tilde_psi_star_exp[[i]][[j]],psitildemat_exp2[[i]][[j]][[k]])
        }
      }
    }
    
    ###
    # now let's predict all test cases
    # all cases means we need to use all matrices in the list tilde_psi_star[[i]][[j]] 
    ###
    ninside2 <- 0
    sum_interval <- 0
    sumerror <- 0
    all_pmean <- c()
    all_length <- c()
    
    for(i in 1:length(tilde_psi_star_exp)){
      for(j in 1:length(tilde_psi_star_exp[[i]])){
        
        print(paste0("i = ",i,"; j = ",j,"  (total: i = ",length(tilde_psi_star_exp),"; j = ",length(tilde_psi_star_exp[[i]]),")"))
        
        
        tilde_psi_star_temp <- tilde_psi_star_exp[[i]][[j]]
        # real value
        real_V <- Vtest[[i]][,j]
        
        # # predictive mean
        # midterm <- diag(nrow(tilde_psi_exp))-tilde_psi_exp %*% solve(diag(ncol(tilde_psi_exp))+t(tilde_psi_exp)%*%tilde_psi_exp) %*% t(tilde_psi_exp)
        # exp <- tilde_psi_star_temp %*% t(tilde_psi_exp) %*% midterm %*% Vtrain
        # pmean <- as.vector(exp)
        # 
        # # predictive covariance
        # pcov <- tilde_psi_star_temp %*% t(tilde_psi_star_temp) + diag(nrow(tilde_psi_star_temp)) -
        #   tilde_psi_star_temp %*% t(tilde_psi_exp) %*% tilde_psi_exp %*% t(tilde_psi_star_temp) + 
        #   tilde_psi_star_temp %*% t(tilde_psi_exp) %*% tilde_psi_exp %*% solve(diag(ncol(tilde_psi_exp))+
        #                                                                          t(tilde_psi_exp)%*%tilde_psi_exp) %*% t(tilde_psi_star_temp %*% t(tilde_psi_exp) %*% tilde_psi_exp)
        # pvar <- diag(pcov)
        
        
        
        # predictive mean
        midterm <- tilde_psi_star_temp %*% t(tilde_psi_exp) %*% tilde_psi_exp %*% (solve( diag(ncol(tilde_psi_exp)) + t(tilde_psi_exp) %*% tilde_psi_exp ) %*% (t(tilde_psi_exp) %*% Vtrain))
        exp <- tilde_psi_star_temp %*% (t(tilde_psi_exp) %*% Vtrain) - midterm
        pmean <- as.vector(exp)
        
        # predictive covariance
        pcov <- (tilde_psi_star_temp %*% t(tilde_psi_star_temp) + (diag(nrow(tilde_psi_star_temp)))) -
          tilde_psi_star_temp %*% t(tilde_psi_exp) %*% tilde_psi_exp %*% t(tilde_psi_star_temp) + 
          tilde_psi_star_temp %*% t(tilde_psi_exp) %*% tilde_psi_exp %*% solve(diag(ncol(tilde_psi_exp))+
                                                                                                        t(tilde_psi_exp)%*%tilde_psi_exp) %*% (t(tilde_psi_exp) %*% (tilde_psi_exp %*% t(tilde_psi_star_temp)))
        pvar <- diag(pcov)
        
        
        
        
        
        
        
        
        
        # 95% CI
        lowers <- pmean-(pvar)^(0.5)*qnorm(p = 0.975)
        uppers <- pmean+(pvar)^(0.5)*qnorm(p = 0.975)
        CI <- cbind(lowers,pmean,uppers,real_V)
        CI <- cbind(CI,CI[,"uppers"]-CI[,"lowers"])
        # calculate how many predictions landed inside the CI
        ninside2 <- ninside2 + sum( (CI[,"real_V"] > CI[,"lowers"]) & (CI[,"real_V"] < CI[,"uppers"]) )
        # sum up interval length
        sum_interval <- sum_interval + sum(CI[,ncol(CI)])
        # RMSE
        sumerror <- sum((CI[,"pmean"]-CI[,"real_V"])^2)
        # store all post mean and all CI length
        all_pmean <- c(all_pmean,pmean)
        all_length <- c(all_length,CI[,ncol(CI)])
      }
    }
    ntotal <- length(Vtest) * length(Vtest[[1]])
    coverage <- c(coverage,ninside2/ntotal)
    avelength <- c(avelength,sum_interval/ntotal)
    RMSE <- c(RMSE,sqrt(sumerror/ntotal))
    RMSE.sdPmean <- c(RMSE.sdPmean,RMSE[length(RMSE)]/sd(all_pmean))
  }
}

result <- data.frame(ns=rep(nsv,each=length(betav)),gamma=rep(betav,times=length(nsv)),
                     coverage=coverage*100,avelength=avelength,RMSE=RMSE,"RMSE/sd(pmean)"=RMSE.sdPmean)
result


# sine function result for V
#     ns gamma coverage avelength       RMSE RMSE.sd.pmean.
# 1   50   0.8 94.90969  3.923100 0.09763126      0.2144010
# 2   50   1.0 94.90969  3.923596 0.09769796      0.2145293
# 3   50   3.0 94.25287  3.922965 0.09966892      0.2775932
# 4   50   5.0 93.76026  3.921539 0.11072756      0.8250313
# 5   50   6.0 93.59606  3.921035 0.10889798      1.0285410
# 6   50   7.0 93.43186  3.920593 0.10806027      9.8937349
# 7  150   0.8 94.90969  3.923100 0.09763126      0.2144010
# 8  150   1.0 94.90969  3.923597 0.09769796      0.2145293
# 9  150   3.0 95.07389  3.926065 0.09427750      0.2089319
# 10 150   5.0 94.41708  3.924472 0.09873475      0.2622624
# 11 150   6.0 93.18555  3.923570 0.10965750      0.4291253
# 12 150   7.0 93.43186  3.922882 0.11194785      0.5668784
# 13 300   0.8 94.90969  3.923100 0.09763126      0.2144010
# 14 300   1.0 94.90969  3.923597 0.09769796      0.2145293
# 15 300   3.0 95.15599  3.928059 0.09595322      0.2103229
# 16 300   5.0 94.99179  3.927129 0.09215638      0.2082940
# 17 300   6.0 94.66338  3.926500 0.09512647      0.2217082
# 18 300   7.0 94.41708  3.925668 0.09874614      0.2622856


plot(result$coverage[1:6]~result$gamma[1:6],type="l",xlab="gamma",ylab="coverage",main="sine function result for V",ylim=c(min(result$coverage),max(result$coverage)))
lines(result$coverage[7:12]~result$gamma[7:12],col="blue")
lines(result$coverage[13:18]~result$gamma[13:18],col="red")
legend("bottomleft", title="Number of basis functions",col = c("black","blue","red"),legend = c(50,150,300),lty=1)

plot(result$avelength[1:6]~result$gamma[1:6],type="l",xlab="gamma",ylab="avelength",main="sine function result for V",ylim=c(min(result$avelength),max(result$avelength)))
lines(result$avelength[7:12]~result$gamma[7:12],col="blue")
lines(result$avelength[13:18]~result$gamma[13:18],col="red")
legend("bottomleft", title="Number of basis functions",col = c("black","blue","red"),legend = c(50,150,300),lty=1)

plot(result$RMSE[1:6]~result$gamma[1:6],type="l",xlab="gamma",ylab="RMSE",main="sine function result for V",ylim=c(min(result$RMSE),max(result$RMSE)))
lines(result$RMSE[7:12]~result$gamma[7:12],col="blue")
lines(result$RMSE[13:18]~result$gamma[13:18],col="red")
legend("topleft", title="Number of basis functions",col = c("black","blue","red"),legend = c(50,150,300),lty=1)

plot(result$RMSE.sd.pmean.[1:6]~result$gamma[1:6],type="l",xlab="gamma",ylab="RMSE/sd(pmean)",main="sine function result for V",ylim=c(min(result$RMSE.sd.pmean.),1.2))
lines(result$RMSE.sd.pmean.[7:12]~result$gamma[7:12],col="blue")
lines(result$RMSE.sd.pmean.[13:18]~result$gamma[13:18],col="red")
legend("topleft", title="Number of basis functions",col = c("black","blue","red"),legend = c(50,150,300),lty=1)


# original function result for V

#     ns gamma coverage avelength      RMSE RMSE.sd.pmean.
# 19  50   0.8 44.17077  3.923652 0.5315356   1.316891e-01
# 20  50   1.0 45.32020  3.924253 0.5343588   1.269433e-01
# 21  50   3.0 53.28407  3.921958 0.3881982   9.824698e-02
# 22  50   5.0 55.33662  3.920233 0.3628032   1.035412e-01
# 23  50   6.0 61.57635  3.919928 0.2535274   3.148502e+01
# 24  50   7.0 61.49425  3.919928 0.2543213   4.103756e+07
# 25 150   0.8 44.17077  3.923652 0.5315364   1.316890e-01
# 26 150   1.0 45.56650  3.924318 0.5359400   1.270542e-01
# 27 150   3.0 52.29885  3.925365 0.3829124   9.530590e-02
# 28 150   5.0 71.10016  3.923082 0.2281168   7.241239e-02
# 29 150   6.0 83.41544  3.922387 0.1517812   5.465495e-02
# 30 150   7.0 79.39245  3.921672 0.1246621   4.869928e-02
# 31 300   0.8 44.17077  3.923652 0.5315364   1.316890e-01
# 32 300   1.0 45.56650  3.924318 0.5359400   1.270542e-01
# 33 300   3.0 52.13465  3.928576 0.3814643   9.407275e-02
# 34 300   5.0 71.83908  3.926291 0.2203067   7.180023e-02
# 35 300   6.0 79.88506  3.924981 0.1864488   6.403101e-02
# 36 300   7.0 86.12479  3.923967 0.1418386   5.355756e-02






#########################################################################################
########              Section 4 - Prediction on Phi with SE Covariance             #########
#########################################################################################

coverage <- c()
avelength <- c()
RMSE <- c()
RMSE.sdPmean <- c()

nsv <- 50 #c(50,150,300)
betav <- 1  #c(0.8,1,3,5,6,7)

for(ii in nsv){
  for(jj in betav){
    
    print(paste0("ns = ",ii,"; beta = ",jj))
    
    ns <- ii
    beta <- jj
    
    psimat_exp <- vector("list",obs_info$M)   # list of length MLN, each is a matrix of psi(d_ijk)
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
          for(l in 1:ns){
            if(l==1){
              psimat_exp[[i]][[j]][[k]][,l] <- exp(-(beta*pdisttemp)^2) * sqrt((2*beta^2)^(l-1)/factorial(l-1)) * pdisttemp^(l-1)
            }else{
              temp <- 1
              for(z in 2:l){
                temp <- sqrt( (2*beta^2) ) / sqrt( z-1 ) * pdisttemp * temp
              }
              psimat_exp[[i]][[j]][[k]][,l] <- exp(-(beta*pdisttemp)^2) * temp
            }
          }
        }
      }
    }
    
    # construct psi tilde
    psitildemat_exp <- vector("list",obs_info$M)   # list of length LMN, each is a matrix of psi(d_ijk)
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
          indx1 <- indx2 - sys_info$d + 1
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
    
    # same as above, on test set
    psimat_exp2 <- vector("list",obs_info$M)   # list of length MLN, each is a matrix of psi(d_ijk)
    for(i in 1:length(psimat_exp2)){
      psimat_exp2[[i]] <- vector("list",ncol(Vtest[[i]]))
    }
    for(i in 1:length(psimat_exp2)){
      for(j in 1:length(psimat_exp2[[i]])){
        psimat_exp2[[i]][[j]] <- vector("list",sys_info$N)
      }
    }
    
    for(i in 1:length(psimat_exp2)){
      for(j in 1:length(psimat_exp2[[i]])){
        for(k in 1:length(psimat_exp2[[i]][[j]])){
          psimat_exp2[[i]][[j]][[k]] <- matrix(NA,nrow = sys_info$N,ncol = ns)
          pdisttemp <- Dmat2[[i]][[j]][k,]  # distance between agent k with other agents on traj i at time j: d_ijk in notes page: star 1
          for(l in 1:ns){
            if(l==1){
              psimat_exp2[[i]][[j]][[k]][,l] <- exp(-(beta*pdisttemp)^2) * sqrt((2*beta^2)^(l-1)/factorial(l-1)) * pdisttemp^(l-1)
            }else{
              temp <- 1
              for(z in 2:l){
                temp <- sqrt( (2*beta^2) ) / sqrt( z-1 ) * pdisttemp * temp
              }
              psimat_exp2[[i]][[j]][[k]][,l] <- exp(-(beta*pdisttemp)^2) * temp
            }
          }
        }
      }
    }
    
    #########################################################################################################################
    # \tilde\psi(d_ijk*)
    # do not need \tilde\psi(d_ijk*) for predicting \phi, we only need \psi_*. dont need this part
    psitildemat_exp2 <- vector("list",obs_info$M)   # list of length LMN, each is a matrix of psi(d_ijk)
    for(i in 1:length(psitildemat_exp2)){
      psitildemat_exp2[[i]] <- vector("list",ncol(Vtest[[i]]))
    }
    for(i in 1:length(psitildemat_exp2)){
      for(j in 1:length(psitildemat_exp2[[i]])){
        psitildemat_exp2[[i]][[j]] <- vector("list",sys_info$N)
      }
    }
    
    
    for(i in 1:length(psitildemat_exp2)){
      for(j in 1:length(psitildemat_exp2[[i]])){
        for(k in 1:length(psitildemat_exp2[[i]][[j]])){
          indx2 <- k*sys_info$d
          indx1 <- indx2 - sys_info$d + 1
          psitildemat_exp2[[i]][[j]][[k]] <- 1/sys_info$N *Pmat2[[i]][[j]][indx1:indx2,] %*% psimat_exp2[[i]][[j]][[k]]
        }
      }
    }
    
    tilde_psi_star_exp <- vector("list",obs_info$M)
    for(i in 1:length(tilde_psi_star_exp)){
      tilde_psi_star_exp[[i]] <- vector("list",ncol(Vtest[[i]]))
    }
    
    for(i in 1:length(psitildemat_exp2)){
      for(j in 1:length(psitildemat_exp2[[i]])){
        for(k in 1:length(psitildemat_exp2[[i]][[j]])){
          tilde_psi_star_exp[[i]][[j]] <- rbind(tilde_psi_star_exp[[i]][[j]],psitildemat_exp2[[i]][[j]][[k]])
        }
      }
    }
    #########################################################################################################################
    # According to the fomulae in section 6 in the report, we need to calculate \psi_* and use \tilde\psi calculated above
    # still assume \sigma^2 = 1
    
    # construct \psi_*
    psi_star_exp <- vector("list",length(psimat_exp2))
    for(i in 1:length(psi_star_exp)){
      psi_star_exp[[i]] <- vector("list",length(psimat_exp2[[i]]))
    }
    
    for(i in 1:length(psi_star_exp)){
      for(j in 1:length(psi_star_exp[[i]])){
        temp <- c()
        for(k in 1:length(psimat_exp2[[i]][[j]])){
          temp <- rbind(temp,psimat_exp2[[i]][[j]][[k]])
        }
        psi_star_exp[[i]][[j]] <- temp
      }
    }
    
    # calculate real phi's
    real_phi <- vector("list",length(Dmat2))
    for(i in 1:length(real_phi)){
      real_phi[[i]] <- vector("list",length(Dmat2[[i]]))
    }
    
    for(i in 1:length(real_phi)){
      for(j in 1:length(real_phi[[i]])){
        pairdist <- as.vector(Dmat2[[i]][[j]])
        real_phi[[i]][[j]] <- sys_info$phiE[[1]][[1]](pairdist)
      }
    }
    
    ###
    # now let's predict all test cases
    # all cases means we need to use all matrices in the list psi_star_exp[[i]][[j]] 
    ###
    ninsidephi <- 0
    sum_interval <- 0
    sumerror <- 0
    ntotal <- 0
    all_pmean <- c()
    all_length <- c()
    
    for(i in 1:length(psi_star_exp)){
      for(j in 1:length(psi_star_exp[[i]])){
        
        print(paste0("i = ",i,"; j = ",j,"  (total: i = ",length(psi_star_exp),"; j = ",length(psi_star_exp[[i]]),")"))
        
        
        psi_star_temp <- psi_star_exp[[i]][[j]]
        # real value
        real_phi_temp <- real_phi[[i]][[j]]
        
        # predictive mean
        midterm <- diag(nrow(tilde_psi_exp))-tilde_psi_exp %*% solve(diag(ncol(tilde_psi_exp))+t(tilde_psi_exp)%*%tilde_psi_exp) %*% t(tilde_psi_exp)
        exp <- psi_star_temp %*% t(tilde_psi_exp) %*% midterm %*% Vtrain
        pmean <- as.vector(exp)
        
        # predictive covariance
        pcov <- psi_star_temp %*% t(psi_star_temp) - psi_star_temp %*% t(tilde_psi_exp) %*% midterm %*% t(psi_star_temp %*% t(tilde_psi_exp))
        pvar <- diag(pcov)
        
        # 95% CI
        lowers <- pmean-(pvar)^(0.5)*qnorm(p = 0.975)
        uppers <- pmean+(pvar)^(0.5)*qnorm(p = 0.975)
        CI <- cbind(lowers,pmean,uppers,real_phi_temp)
        CI <- cbind(CI,CI[,"uppers"]-CI[,"lowers"])
        # only consider those real_phi which are not -Inf
        indx <- which(CI[,"real_phi_temp"] != -Inf)
        CI <- CI[indx,]
        # calculate how many predictions landed inside the CI
        ninsidephi <- ninsidephi + sum( (CI[,"real_phi_temp"] > CI[,"lowers"]) & (CI[,"real_phi_temp"] < CI[,"uppers"]) )
        # sum up interval length
        sum_interval <- sum_interval + sum(CI[,ncol(CI)])
        # RMSE
        sumerror <- sum((CI[,"pmean"]-CI[,"real_phi_temp"])^2)  # do not include those Inf cases
        # total number of predictions (different from section 3, here need to calculate n manually each iteration since we may remove some results -- -Inf)
        ntotal <- ntotal + nrow(CI)
        # store all post mean and all CI length
        all_pmean <- c(all_pmean,pmean)
        all_length <- c(all_length,CI[,ncol(CI)])
      }
    }
    
    coverage <- c(coverage,ninsidephi/ntotal)
    avelength <- c(avelength,sum_interval/ntotal)
    RMSE <- c(RMSE,sqrt(sumerror/ntotal))
    RMSE.sdPmean <- c(RMSE.sdPmean,RMSE[length(RMSE)]/sd(all_pmean))
  }
}

result2 <- data.frame(ns=rep(nsv,each=length(betav)),gamma=rep(betav,times=length(nsv)),
                      coverage=coverage*100,avelength=avelength,RMSE=RMSE,"RMSE/sd(pmean)"=RMSE.sdPmean)
result2




# sine function result for phi
#     ns gamma  coverage avelength       RMSE RMSE.sd.pmean.
# 1   50   0.8 100.00000 0.7413978 0.01876677     0.05810367
# 2   50   1.0 100.00000 0.8383011 0.01300148     0.03731592
# 3   50   3.0  51.02041 0.9548302 0.04712745     0.08486943
# 4   50   5.0  22.44898 0.8256710 0.09636630     0.16582739
# 5   50   6.0  18.36735 0.7582080 0.09817906     0.20619524
# 6   50   7.0  18.36735 0.7094760 0.07847651     1.08622342
# 7  150   0.8 100.00000 0.7413978 0.01876677     0.05810367
# 8  150   1.0 100.00000 0.8383050 0.01300152     0.03731608
# 9  150   3.0  85.03401 1.1451220 0.02179074     0.05192536
# 10 150   5.0  36.89890 1.0527645 0.07145442     0.10283728
# 11 150   6.0  29.25170 0.9903441 0.09736073     0.13123382
# 12 150   7.0  27.89116 0.9265913 0.09439378     0.14980170
# 13 300   0.8 100.00000 0.7413978 0.01876677     0.05810367
# 14 300   1.0 100.00000 0.8383050 0.01300152     0.03731608
# 15 300   3.0  93.19728 1.2573010 0.01645439     0.03861771
# 16 300   5.0  62.28008 1.2053507 0.03393718     0.06128494
# 17 300   6.0  53.64767 1.1715852 0.04647234     0.07796000
# 18 300   7.0  37.55571 1.1089106 0.07953828     0.10666707




# original function result for phi
#     ns gamma   coverage    avelength      RMSE RMSE.sd.pmean.
# 1   50   0.8  0.5473454 4.536525e-01 1.5766376   7.044322e-02
# 2   50   1.0  2.2441160 4.848383e-01 1.5141937   7.580448e-02
# 3   50   3.0  0.6020799 2.878054e-01 1.1847660   8.431954e-02
# 4   50   5.0  0.2189381 1.024040e-01 1.2432440   9.373161e-02
# 5   50   6.0  0.0000000 2.740096e-03 1.5172394   4.336432e+01
# 6   50   7.0  0.0000000 9.672399e-07 1.5221393   5.445524e+07
# 7  150   0.8  0.5473454 4.536577e-01 1.5766288   7.044314e-02
# 8  150   1.0  3.6672140 4.872649e-01 1.5039428   7.556051e-02
# 9  150   3.0 10.6732348 5.425693e-01 1.0327377   8.143735e-02
# 10 150   5.0  1.8609743 3.301962e-01 0.5512032   4.356472e-02
# 11 150   6.0  1.3683634 2.799375e-01 0.4084547   3.295754e-02
# 12 150   7.0  1.2041598 2.281227e-01 0.4567296   3.815058e-02
# 13 300   0.8  0.5473454 4.536577e-01 1.5766288   7.044314e-02
# 14 300   1.0  3.6672140 4.872649e-01 1.5039428   7.556051e-02
# 15 300   3.0  3.5577449 6.758333e-01 1.0361406   8.299048e-02
# 16 300   5.0 12.0963328 5.751752e-01 0.4101111   3.554851e-02
# 17 300   6.0  6.6776136 4.535209e-01 0.2899306   2.628570e-02
# 18 300   7.0  6.0755337 3.627634e-01 0.3174778   2.682621e-02

plot(result2$coverage[1:6]~result2$gamma[1:6],type="l",xlab="gamma",ylab="coverage",main="sine function result for phi",ylim=c(min(result2$coverage),max(result2$coverage)))
lines(result2$coverage[7:12]~result2$gamma[7:12],col="blue")
lines(result2$coverage[13:18]~result2$gamma[13:18],col="red")
legend("bottomleft", title="Number of basis functions",col = c("black","blue","red"),legend = c(50,150,300),lty=1)

plot(result2$avelength[1:6]~result2$gamma[1:6],type="l",xlab="gamma",ylab="avelength",main="sine function result for phi",ylim=c(min(result2$avelength),max(result2$avelength)))
lines(result2$avelength[7:12]~result2$gamma[7:12],col="blue")
lines(result2$avelength[13:18]~result2$gamma[13:18],col="red")
legend("topleft", title="Number of basis functions",col = c("black","blue","red"),legend = c(50,150,300),lty=1)

plot(result2$RMSE[1:6]~result2$gamma[1:6],type="l",xlab="gamma",ylab="RMSE",main="sine function result for phi",ylim=c(min(result2$RMSE),max(result2$RMSE)))
lines(result2$RMSE[7:12]~result2$gamma[7:12],col="blue")
lines(result2$RMSE[13:18]~result2$gamma[13:18],col="red")
legend("topleft", title="Number of basis functions",col = c("black","blue","red"),legend = c(50,150,300),lty=1)

plot(result2$RMSE.sd.pmean.[1:6]~result2$gamma[1:6],type="l",xlab="gamma",ylab="RMSE/sd(pmean)",main="sine function result for phi",ylim=c(min(result2$RMSE.sd.pmean.),max(result2$RMSE.sd.pmean.)))
lines(result2$RMSE.sd.pmean.[7:12]~result2$gamma[7:12],col="blue")
lines(result2$RMSE.sd.pmean.[13:18]~result2$gamma[13:18],col="red")
legend("topleft", title="Number of basis functions",col = c("black","blue","red"),legend = c(50,150,300),lty=1)


