## scale_the_quantities

scale_the_quantities <- function(F_vec, d_vec, energy_Psi, N, K,class_info, time_vec,Riemann){
  
  print("start scale the quantities")
  
  L <- length(time_vec)
  d <- length(d_vec)/(L*N)
  # class_mod <- get_class_modifier(...)
  class_mod <- rep(0,N)
  
  #print(class_mod)
  
  for(k in 1:K){
    ind  <-  class_info == k
    class_mod[as.logical(ind)] <- sqrt(sum(ind))
  }
  
  class_mod <- rep(class_mod,each=d)
  class_mod <- rep(class_mod,L)
  
  #print(head(class_mod))

  
  #print("start get_Riemann_based_time_modifier")
  time_mod <- get_Riemann_based_time_modifier(N, d, time_vec, Riemann)
  #print("finish get_Riemann_based_time_modifier")
  #print(head(time_mod))
  
  D_mod <- diag(time_mod / class_mod)
  
  #print((D_mod[1:6,1:6]))
  
  # print(length(time_vec))
  # print(length(d_vec))
  # print(dim(energy_Psi))
  # print(dim(D_mod))

  #F_vec <- D_mod %*% F_vec
  d_vec <- D_mod %*% d_vec
  energy_Psi <- D_mod %*% energy_Psi
  
  result <- list()
  #result[["F_vec"]] <- F_vec
  result[["d_vec"]] <- d_vec
  result[["energy_Psi"]] <- energy_Psi
  
  print("finish scale the quantities")
  
  return(result)
}