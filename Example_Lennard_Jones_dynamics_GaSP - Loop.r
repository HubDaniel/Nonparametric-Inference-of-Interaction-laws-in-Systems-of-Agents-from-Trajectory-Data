# Lennard_Jones_dynamics with GaSP

setwd("C:/Users/x_r_m/Desktop/Learn/RA/V2")
# load packages used below
source("generate_obs.r")
source("learn_from_dynamics.r")
source("Example_Lennard_Jones_dynamics_GaSP-LF.r") # this contains the likelihood function we're going to min

########################################################################################
########################################################################################
########################################################################################


# need sigma0 = 0.1 N = 7 M = 10 L = 200 ns = 50

sigma_0v <- c(sqrt(0.1),1,sqrt(10)) # standard deviation not variance
# if sigma0 is too small, may result in etaMLE=sigma02/sigma2 very small and later result in unstable computation
Nv <- c(3,7,10)
Mv <- c(3,10)
Lv <- c(30,150)
nsv <- c(15,30,50,100)   # when ns is very small, we may have negative posterior variance, which gives NA in coverage


gammaMLEv <- c()
etaMLEv <- c()
sigma2MLEv <- c()
sigma02MLEv <- c()
LFvalue <- c()
coverage <- c()
avelength <- c()
RMSE <- c()
RMSE.sdPmean <- c()
coverage_phiv <- c()
avelength_phiv <- c()
RMSE_phiv <- c()
RMSE.sdPmean_phiv <- c()



for(kk in sigma_0v){
  for(ll in Nv){
    for(ii in Mv){
      for(jj in Lv){
        
        
        sys_info                <- list()
        sys_info$d              <- 2
        sys_info$N              <- ll
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
        obs_info$M              <- ii
        obs_info$L              <- jj
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
        sys_info$noise_sigma <- kk 
        
        
        
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
        
        Vtrain <- unlist(Vtrain)
        
        
        
        
        #########################################################################################
        ########              Section 1 - Hyperparameter Estimation              #########
        #########################################################################################
        
        
        for(ns in nsv){
          print(paste0("estimating hyperparameter... sigma0 = ",kk," N = ",ll," M = ",ii," L = ",jj," ns = ",ns))
          est_result <- optim(par = c(1,0), fn = likelihoodFn,ns=ns )  # ns: number of basis functions passed into likelihoodFn
          
          gammaMLE <- est_result$par[1]
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
                    psimat_exp_MLE[[i]][[j]][[k]][,l] <- exp(-(gammaMLE*pdisttemp)^2) * sqrt((2*gammaMLE^2)^(l-1)/factorial(l-1)) * pdisttemp^(l-1)
                  }else{
                    temp <- 1
                    for(z in 2:l){
                      temp <- sqrt( (2*gammaMLE^2) ) / sqrt( z-1 ) * pdisttemp * temp
                    }
                    psimat_exp_MLE[[i]][[j]][[k]][,l] <- exp(-(gammaMLE*pdisttemp)^2) * temp
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
          
          if(nchar(etaMLE-round(etaMLE)) > 7){    # if etaMLE has too many decimals, the reciprocal of etaMLE results in unstable calculation
            sigma2MLE <- (1/length(Vtrain)) * t(Vtrain) %*% solve(tilde_psi_exp_MLE %*% t(tilde_psi_exp_MLE) + etaMLE * diag(nrow(tilde_psi_exp_MLE))) %*% Vtrain
          }else{                # if etaMLE is not very small, the following calculation is faster
            sigma2MLE <- (1/length(Vtrain)) *( (1/etaMLE) * t(Vtrain) %*% Vtrain - 
                                                 (1/etaMLE)^2 * (t(Vtrain) %*% tilde_psi_exp_MLE) %*% solve(diag(ncol(tilde_psi_exp_MLE))+1/etaMLE * t(tilde_psi_exp_MLE) %*% tilde_psi_exp_MLE) %*% (t(tilde_psi_exp_MLE) %*% Vtrain)  ) 
          }

          sigma2MLE <- as.vector(sigma2MLE)
          
          sigma02MLE <- etaMLE*sigma2MLE  # should be close to 1 (set it to be 1)
          
          sigma02MLE <- as.vector(sigma02MLE)
          
          gammaMLEv <- c(gammaMLEv,gammaMLE)
          etaMLEv <- c(etaMLEv,etaMLE)
          sigma2MLEv <- c(sigma2MLEv,sigma2MLE)
          sigma02MLEv <- c(sigma02MLEv,sigma02MLE)
          LFvalue <- c(LFvalue,est_result$value)
          
          
          beta <- gammaMLE # it's the gammaMLE
          
          
          tilde_psi_exp <- tilde_psi_exp_MLE
          
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
              
              print(paste0("V"," sigma0 = ",kk, " N = ", ll ," M = ",ii," L = ",jj,
                           " ns = ",ns," i = ",i,"; j = ",j,"  (total: i = ",length(tilde_psi_star_exp),"; j = ",length(tilde_psi_star_exp[[i]]),")"))
              
              
              tilde_psi_star_temp <- tilde_psi_star_exp[[i]][[j]]
              # real value
              real_V <- Vtest[[i]][,j]
              
              
              # if etaMLE has too many digits, the following calculation is unstable since it involves the 
              # square of the reciprocal of etaMLE. Also, since we used the square of the reciprocal of eta
              # in the likelihood function, an etaMLE with too many digits also indicates an inaccurate estimation
              
              
              # for further usage, we can force the result to be NA if:
              # nchar(etaMLE-round(etaMLE)) > 7
              
              
              
              
              # if etaMLE does not has too many digits, the following calculation is faster
              # predictive mean
              midterm <- tilde_psi_star_temp %*% t(tilde_psi_exp) %*% tilde_psi_exp %*% (solve( diag(ncol(tilde_psi_exp)) + (1/etaMLE) * t(tilde_psi_exp) %*% tilde_psi_exp ) %*% (t(tilde_psi_exp) %*% Vtrain))
              exp <- 1/etaMLE * tilde_psi_star_temp %*% (t(tilde_psi_exp) %*% Vtrain) - ((1/etaMLE)^2) * midterm
              pmean <- as.vector(exp)
              
              # predictive covariance
              pcov <- sigma2MLE * (tilde_psi_star_temp %*% t(tilde_psi_star_temp) + (etaMLE * diag(nrow(tilde_psi_star_temp)))) -
                sigma2MLE/etaMLE * tilde_psi_star_temp %*% t(tilde_psi_exp) %*% tilde_psi_exp %*% t(tilde_psi_star_temp) + 
                sigma2MLE/(etaMLE^2) * tilde_psi_star_temp %*% t(tilde_psi_exp) %*% tilde_psi_exp %*% solve(diag(ncol(tilde_psi_exp))+
                                                                                                              1/etaMLE * t(tilde_psi_exp)%*%tilde_psi_exp) %*% (t(tilde_psi_exp) %*% (tilde_psi_exp %*% t(tilde_psi_star_temp)))
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
          
          
          
          
          
          #########################################################################################
          ########              Prediction on Phi with SE Covariance             #########
          #########################################################################################
          
          
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
              
              print(paste0("Phi"," sigma0 = ",kk, " N = ", ll , " M = ",ii," L = ",jj,
                           " ns = ",ns," i = ",i,"; j = ",j,"  (total: i = ",length(psi_star_exp),"; j = ",length(psi_star_exp[[i]]),")"))
              
              
              psi_star_temp <- psi_star_exp[[i]][[j]]
              # real value
              real_phi_temp <- real_phi[[i]][[j]]
              
              
              
              # if etaMLE has too many digits, the following calculation is unstable since it involves the 
              # square of the reciprocal of etaMLE. Also, since we used the square of the reciprocal of eta
              # in the likelihood function, an etaMLE with too many digits also indicates an inaccurate estimation
              
              
              # for further usage, we can force the result to be NA if:
              # nchar(etaMLE-round(etaMLE)) > 7
              
              
              
              
              # if etaMLE does not has too many digits, the following calculation is faster
              # predictive mean
              midterm <- psi_star_temp %*% t(tilde_psi_exp) %*% tilde_psi_exp %*% (solve( diag(ncol(tilde_psi_exp)) + 1/etaMLE * t(tilde_psi_exp) %*% tilde_psi_exp ) %*% (t(tilde_psi_exp) %*% Vtrain))
              exp <- 1/etaMLE * psi_star_temp %*% (t(tilde_psi_exp) %*% Vtrain) - (1/etaMLE)^2 *midterm
              pmean <- as.vector(exp)
              
              # predictive covariance
              pcov1 <- sigma2MLE * psi_star_temp %*% t(psi_star_temp)
              pcov2 <- (sigma2MLE/etaMLE) * psi_star_temp %*% t(tilde_psi_exp) %*% tilde_psi_exp %*% t(psi_star_temp)
              pcov3inv <- solve( diag(ncol(tilde_psi_exp)) + (1/etaMLE) * t(tilde_psi_exp) %*% tilde_psi_exp )
              pcov3 <- (sigma2MLE/(etaMLE^2)) * psi_star_temp %*% t(tilde_psi_exp) %*% tilde_psi_exp %*% pcov3inv %*% (t(tilde_psi_exp) %*% (tilde_psi_exp %*% t(psi_star_temp)))
              pcov <- pcov1-pcov2+pcov3
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
          
          coverage_phiv <- c(coverage_phiv,ninsidephi/ntotal)
          avelength_phiv <- c(avelength_phiv,sum_interval/ntotal)
          RMSE_phiv <- c(RMSE_phiv,sqrt(sumerror/ntotal))
          RMSE.sdPmean_phiv <- c(RMSE.sdPmean_phiv,RMSE_phiv[length(RMSE_phiv)]/sd(all_pmean))
          
        }
        
        
        
      }
    }
  }
  
}




result <- data.frame(real_sigma_02=rep(sigma_0v^2,each=(length(Nv)*length(Mv)*length(Lv)*length(nsv))),
                     N=rep(Nv,each=(length(Mv)*length(Lv)*length(nsv))),
                     M=rep(Mv,each=(length(Lv)*length(nsv))),
                     L=rep(Lv,each=length(nsv)),
                     ns=nsv,gamma=gammaMLEv,gamma2=gammaMLEv^2,eta=etaMLEv,sigma2=sigma2MLEv,sigma02=sigma02MLEv,LFvalue=LFvalue,
                     coverage=coverage*100,avelength=avelength,RMSE=RMSE,"RMSE/sd(pmean)"=RMSE.sdPmean,coverage_phi=coverage_phiv*100,
                     avelength_phi=avelength_phiv,RMSE_phi=RMSE_phiv,RMSE.sdPmean_phi=RMSE.sdPmean_phiv)
result  
# NA resulted from negative posterior variance if etaMLE has too many digits, especially when its too small,
# since we used the square of the reciprocal of etaMLE, the calculation is unstable. I checked the code, it's okay.

# note that a very small etaMLE may also indicates an inaccurate MLE result since we used the square of 
# the reciprocal of eta in the likelihood funciton.

# simply ignore the rows with NA results.




















