## predator_swarm example

setwd("C:/Users/xrma/Desktop/ucsb/A ѧϰ/RA/V2")

########################################################################################
########################################################################################
########################################################################################
# set ups for Predator Prey
N_preys                 <- 9
N_predators             <- 1
prey_attract_prey       <- 1
predator_repulse_prey   <- 2
prey_attract_predator   <- 3.5
predator_sense_prey     <- 3

sys_info                <- list()
sys_info$d              <- 2
sys_info$N              <- N_preys + N_predators
sys_info$K              <- 2
sys_info$type_info      <- c(rep(1,N_preys),rep(2,N_predators))  
sys_info$kappa          <- rep(1,sys_info$K)
sys_info$T_f            <- 10
sys_info$mu0            <- function() PS_init_config(sys_info$N, sys_info$type_info, 1)
sys_info$phiE           <- list()
sys_info$phiE[[1]]      <- list()
sys_info$phiE[[2]]      <- list()
sys_info$phiE[[1]][[1]] <- function(r,prey_attract_prey=1) {prey_attract_prey - r^(-2)}
sys_info$phiE[[1]][[2]] <- function(r,predator_repulse_prey=2) {-predator_repulse_prey * r^(-2)}
sys_info$phiE[[2]][[1]] <- function(r, prey_attract_predator=3.5, predator_sense_prey=3) {prey_attract_predator * r^(-predator_sense_prey)}
sys_info$phiE[[2]][[2]] <- function(r) {rep(0,length(r))}

obs_info                <- list()
obs_info$L              <- 200
obs_info$M              <- 20
obs_info$M_rhoT         <- 2000
obs_info$T_L            <- sys_info$T_f/2
obs_info$time_vec       <- seq(from = 0, to = obs_info$T_L, length.out = obs_info$L)   # time vec for generation
obs_info$obs_time_vec   <- seq(from = 0, to = obs_info$T_L, length.out = obs_info$L)   # time vec for obs.
obs_info$hist_num_bins  <- 10000

basis_info <- list()
basis_info$n            <- pmax(64 * matrix(1,ncol = sys_info$K,nrow = sys_info$K), 
                              matrix(c(ceiling(obs_info$L * obs_info$M * N_preys * sys_info$d/500),
                               ceiling(obs_info$L * obs_info$M * sqrt(N_preys * N_predators) * sys_info$d/500),
                               ceiling(obs_info$L * obs_info$M * sqrt(N_preys * N_predators) * sys_info$d/500),
                               ceiling(obs_info$L * obs_info$M * N_predators * sys_info$d/500)),ncol=2))
basis_info$type         <- 'standard'                                                             
basis_info$degree       <- matrix(c(1, 1, 1, 0),ncol=2) 

learn_info              <- list()
learn_info$Ebasis_info  <- basis_info
learn_info$solver_type  <- "svd"


sys_info$noise_on_obs <- F
sys_info$noise_type <- "add"  # or mult
sys_info$noise_sigma <- 0.001

########################################################################################
########################################################################################
########################################################################################
# generate data

source("generate_obs.r")
obs <- generateObs(sys_info, obs_info)

# list of length 2, ICs: init conditions; traj: list of 50 traj

########################################################################################
########################################################################################
########################################################################################

# learn_from_dynamics use obs
source("learn_from_dynamics.r")
learn <- learn_from_dynamics(sys_info, obs_info, learn_info, obs_data=obs,Riemann=0)


########################################################################################
########################################################################################
########################################################################################

# estimateRhoLT use obs
source("estimateRhoLT.r")
rhoLT <- estimateRhoLT(obs_data=obs, sys_info, obs_info)


########################################################################################
########################################################################################
########################################################################################

# estimated true RhoLT use many traj
temp <- obs_info$M
obs_info$M <- 100   # in paper they used 2000
obs2 <- generateObs(sys_info, obs_info)
TrhoLT <- estimateRhoLT(obs_data=obs2, sys_info, obs_info)
obs_info$M <- temp


########################################################################################
########################################################################################
########################################################################################

# plot all 3 together
# 1 on 1 (prey on prey)
x_values <- seq(from=0.1,to=learn$Estimator$Ebasis[[1]][[1]]$Rmax,length.out = 10000)
rr <- learn$Estimator$phiEhat[[1]][[1]](x_values)
plot(x_values,sys_info$phiE[[1]][[1]](x_values),type="l",ylab = "kernel",xlab="r (pairwise distance)")
lines(x_values,rr$y,type = "l",col="blue",lty=2)

par(new = TRUE)

plot(stepfun(rhoLT$histedgesR[[1]][[1]][2:(length(rhoLT$histedgesR[[1]][[1]])-1)],rhoLT$hist[[1]][[1]]),
     col="red",main="",xlim=c(x_values[1],x_values[length(x_values)]),axes=F,xlab="",ylab="")
axis(side=4,at = pretty(range(rhoLT$hist[[1]][[1]])))
mtext("density",side=4)
lines(stepfun(TrhoLT$histedgesR[[1]][[1]][2:(length(TrhoLT$histedgesR[[1]][[1]])-1)],TrhoLT$hist[[1]][[1]]),
      col="orange")


# 1 on 2 (prey on predator)
x_values <- seq(from=0.4,to=learn$Estimator$Ebasis[[1]][[2]]$Rmax,length.out = 10000)
rr <- learn$Estimator$phiEhat[[1]][[2]](x_values)
plot(x_values,sys_info$phiE[[1]][[2]](x_values),type="l",ylab = "kernel",xlab="r (pairwise distance)")
lines(x_values,rr$y,type = "l",col="blue",lty=2)

par(new = TRUE)

plot(stepfun(rhoLT$histedgesR[[1]][[2]][2:(length(rhoLT$histedgesR[[1]][[2]])-1)],rhoLT$hist[[1]][[2]]),
     col="red",main = "estimated rhoLT",xlim=c(x_values[1],x_values[length(x_values)]),axes=F,
     xlab="",ylab="")
axis(side=4,at = pretty(range(rhoLT$hist[[1]][[2]])))
mtext("density",side=4)
lines(stepfun(TrhoLT$histedgesR[[1]][[2]][2:(length(TrhoLT$histedgesR[[1]][[2]])-1)],TrhoLT$hist[[1]][[2]]),
      col="orange")


# 2 on 1 (predator on prey)
x_values <- seq(from=0.4,to=learn$Estimator$Ebasis[[2]][[1]]$Rmax,length.out = 10000)
rr <- learn$Estimator$phiEhat[[2]][[1]](x_values)
plot(x_values,sys_info$phiE[[2]][[1]](x_values),type="l",ylab = "kernel",xlab="r (pairwise distance)")
lines(x_values,rr$y,type = "l",col="blue",lty=2)

par(new = TRUE)

plot(stepfun(rhoLT$histedgesR[[2]][[1]][2:(length(rhoLT$histedgesR[[2]][[1]])-1)],rhoLT$hist[[2]][[1]]),
     col="red",main = "estimated rhoLT",xlim=c(x_values[1],x_values[length(x_values)]),axes=F,
     xlab="",ylab="")
axis(side=4,at = pretty(range(rhoLT$hist[[2]][[1]])))
mtext("density",side=4)
lines(stepfun(TrhoLT$histedgesR[[2]][[1]][2:(length(TrhoLT$histedgesR[[2]][[1]])-1)],TrhoLT$hist[[2]][[1]]),
      col="orange")

# 2 on 2 (predator on predator)
x_values <- seq(from=0,to=learn$Estimator$Ebasis[[2]][[2]]$Rmax,length.out = 10000)
rr <- learn$Estimator$phiEhat[[2]][[2]](x_values)
plot(x_values,sys_info$phiE[[2]][[2]](x_values),type="l",ylab = "kernel",xlab="r (pairwise distance)")
lines(x_values,rr$y,type = "l",col="blue",lty=2)
mtext("density",side=4)


########################################################################################
########################################################################################
########################################################################################

# compute relative error

source("relativeErrorInfluenceFunction.r")    # mother function is computeL2rhoTErr.m (not used here)

# use estimated rhoLTM to calculate. 
# in their code, they use estimated true rhoLT. The only difference is to re-generate data, but this time 
# with many more traj M. Then use this data to estimate rhoLT and put the result in below
Err <- relativeErrorInfluenceFunction(learn$Estimator$phiEhat,sys_info$phiE,sys_info,
                                      learn$Estimator$Ebasis,rhoLT)
Err

# below is re-generating obs with a large M and estimated true rhoLT

Err2 <- relativeErrorInfluenceFunction(learn$Estimator$phiEhat,sys_info$phiE,sys_info,
                                      learn$Estimator$Ebasis,TrhoLT)
Err2



########################################################################################
########################################################################################
########################################################################################

# compute trajectory error for both the learning period and future period

# idea: use the phihat estimated above, 
# 1. generate a traj with true kernel and initial condition 
# (this initial condition is the same as in above learning process)
# 2. generate a traj with estimated kernel (phihat) and the same initial condition
# 3. error is the maximum of norm of their diff.

source("computeTrajErr.r")  
TrajErr <- computeTrajErr(sys_info, obs_info, obs$ICs)  
# errors may occur but won't stop the program. will remove the error traj before calculating error


which.max(TrajErr$sup_err)
sd(TrajErr$sup_err)
which.max(TrajErr$sup_err_fut)
sd(TrajErr$sup_err_fut)
# the error bound in the paper (I think) is after running ten learnings, each time with diff initial condition



########################################################################################
########################################################################################
########################################################################################
# plot estimated trajectory and true trajectory
# TrajErr will return the estimated traj


traj_to_plot <- 3

rplot <- as.data.frame(t(TrajErr$trajhat[,,traj_to_plot]))
rplot2 <- as.data.frame(t(TrajErr$trajt[,,traj_to_plot]))
rplot3 <- as.data.frame(t(TrajErr$trajfuthat[,,traj_to_plot]))
rplot4 <- as.data.frame(t(TrajErr$trajtfut[,,traj_to_plot]))

rplotodd <- rplot[,seq(from=1,to=ncol(rplot),by=2)]     # columns of coordinate x
rploteve <- rplot[,seq(from=2,to=ncol(rplot),by=2)]     # columns of coordinate y
rplot2odd <- rplot2[,seq(from=1,to=ncol(rplot),by=2)]
rplot2eve <- rplot2[,seq(from=2,to=ncol(rplot),by=2)]
rplot3odd <- rplot3[,seq(from=1,to=ncol(rplot),by=2)]
rplot3eve <- rplot3[,seq(from=2,to=ncol(rplot),by=2)]
rplot4odd <- rplot4[,seq(from=1,to=ncol(rplot),by=2)]
rplot4eve <- rplot4[,seq(from=2,to=ncol(rplot),by=2)]

xmin <- min(min(rplotodd),min(rplot2odd),min(rplot3odd),min(rplot4odd))
xmax <- max(max(rplotodd),max(rplot2odd),max(rplot3odd),max(rplot4odd))

ymin <- min(min(rploteve),min(rplot2eve),min(rplot3eve),min(rplot4eve))
ymax <- max(max(rploteve),max(rplot2eve),max(rplot3eve),max(rplot4eve))


plot(rplot$V1,rplot$V2,xlim=c(xmin,xmax),ylim=c(ymin,ymax),type="l",xlab="x1",ylab="x2",lty=2) # trajhat from t0~tL
for(i in seq(from=3,to=ncol(rplot),by=2)){
  if(i==19){
    lines(rplot[,i],rplot[,i+1],col="red",lty=2)    # predator
  }else{
    lines(rplot[,i],rplot[,i+1],lty=2)
  }
}

for(i in seq(from=1,to=ncol(rplot2),by=2)){                           # trajt from t0~tL
  if(i==19){
    lines(rplot2[,i],rplot2[,i+1],col="red")    # predator
  }else{
    lines(rplot2[,i],rplot2[,i+1])
  }
}

for(i in seq(from=1,to=ncol(rplot3),by=2)){                           # trajfuthat from tL~tF
  if(i==19){
    lines(rplot3[,i],rplot3[,i+1],col="orange",lty=2)    # predator
  }else{
    lines(rplot3[,i],rplot3[,i+1],lty=2,col="brown")
  }
}

for(i in seq(from=1,to=ncol(rplot4),by=2)){                           # trajtfut from tL~tF
  if(i==19){
    lines(rplot4[,i],rplot4[,i+1],col="orange")    # predator
  }else{
    lines(rplot4[,i],rplot4[,i+1],col="brown")
  }
}


