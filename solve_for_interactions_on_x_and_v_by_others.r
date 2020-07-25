##  solve_for_interactions_on_x_and_v_by_others

solve_for_interactions_on_x_and_v_by_others <- function(Phi, rhs_vec, T_L, M,solver_type){
  
  print("solve interact on x and v by others")
  #print(Phi)
  #print(rhs_vec)
  
  if(solver_type == "svd"){
    
    #print(dim(Phi))

    
    U <- svd(Phi)$u
    lambdas <- svd(Phi)$d
    V <- svd(Phi)$v
    
    #print(dim(U))
    #print(dim(V))
    
    
    the_zero <- 10^(-12)
    ind <- which(lambdas <= the_zero)[1]
    print("first ind less than the_zero = 10^(-12)")
    print(ind)
    if(!is.na(ind)){
      if(ind > 2){
        lambdas <- lambdas[1:(ind-1)]
        lam_len <- length(lambdas)
        ratios <- log(lambdas[1:(lam_len-1)]) - log(lambdas[2:lam_len])   # sharp drop off
        indd <- which.max(ratios)
        lambdas_cut <- lambdas[1:indd]
        U_cut <- U[,1:indd]
        V_cut <- V[,1:indd]
        
        if(length(lambdas_cut)==1){
          D_cut <- matrix(1/lambdas_cut,ncol=1,nrow=1)
        }else{
          D_cut <- diag(1/lambdas_cut)
        }
        


        alpha_vec  <- V_cut %*% D_cut %*% t(U_cut) %*% rhs_vec
      }else if(ind == 2){
        lambdas_cut <- lambdas[1]
        U_cut <- U[,1]
        V_cut <- V[,1]
        D_cut <- diag(1/lambdas_cut)
        alpha_vec  <- V_cut %*% D_cut %*% t(U_cut) %*% rhs_vec
      }else{
        alpha_vec <- rep(0,length(rhs_vec))
      }
    }else{
      lam_len          = length(lambdas)
      ratios           = log(lambdas[1 : (lam_len - 1)]) - log(lambdas[2 : lam_len])
      indd             = which.max(ratios)
      lambdas_cut      = lambdas[1 : indd]
      U_cut            = U[, 1 : indd]
      V_cut            = V[, 1 : indd]
      D_cut            = diag(1/lambdas_cut)
      alpha_vec        = V_cut %*% D_cut %*% t(U_cut) %*% rhs_vec
    }
  }else if(solver_type == "pseudo_inverse"){
    alpha_vec <- ginv(Phi) %*% rhs_vec
  }else if(solver_type == "inverse"){
    alpha_vec <- solve(Phi,rhs_vec)
  }
  
  if (T_L ==0){
    opt_val <- (t(Phi %*% alpha_vec) %*% alpha_vec - 2 * t(alpha_vec) %*% rhs_vec)/M
  } else {
    opt_val <- (t(Phi %*% alpha_vec) %*% alpha_vec - 2 * t(alpha_vec) %*% rhs_vec)/(T_L * M)
  }
  fric_coef <- c()
  

  
  result <- list()
  result[["alpha_vec"]] <- alpha_vec
  result[["opt_val"]] <- opt_val
  result["fric_coef"] <- list(NULL)
  
  return(result)
}