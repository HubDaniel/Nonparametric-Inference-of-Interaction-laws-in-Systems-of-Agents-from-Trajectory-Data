# L2_rho_T_energy

# input: 
# f: interaction function / estimated interaction function
# list: rhoLTE: rhoLTE[["supp"]], rhoLTE[["histedgesR"]], rhoLTE[["hist"]]



L2_rho_T_energy <- function(f, sys_info, rhoLTE, basis){
  
  print("start L2_rho_T_energy")
  
  L2rhoTnorm  <- matrix(0,nrow = sys_info$K,ncol = sys_info$K)
  
  for(k1 in 1:sys_info$K){
    N_k1 <- sum(sys_info$type_info==k1)
    for(k2 in 1:sys_info$K){
      if(k1==k2 && N_k1==1){
        L2rhoTnorm[k1, k2] <- 0
      }else{
        range       <- c(basis[[k1]][[k2]]$knots[1],basis[[k1]][[k2]]$knots[length(basis[[k1]][[k2]]$knots)])
        range[1]    <- max(range[1],rhoLTE$supp[[k1]][[k2]][1])
        range[2]    <- min(range[2],rhoLTE$supp[[k1]][[k2]][2])
        edges       <- rhoLTE$histedges[[k1]][[k2]]
        e_idxs      <- which(range[1]<=edges & edges<range[2])
        centers     <- ( edges[e_idxs[1 : (length(e_idxs) - 1)]] + edges[e_idxs[2 : length(e_idxs)]] ) / 2
        # midpoint of each interval in the range
        weights     <- centers
        histdata    <- rhoLTE$hist[[k1]][[k2]][e_idxs[1 : (length(e_idxs) - 1)]]  # hist count in the range
        f_vec       <- f[[k1]][[k2]](centers)  # function evaluated at midpoints
        f_integrand <- ((f_vec * weights)^2) * histdata
        L2rhoTnorm[k1,k2]  <- sqrt(sum(f_integrand))
      }
    }
  }

  print("finish L2_rho_T_energy")
  
  return(L2rhoTnorm)
}



