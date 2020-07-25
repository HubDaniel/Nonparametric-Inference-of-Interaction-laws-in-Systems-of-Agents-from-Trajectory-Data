## uniform_basis_by_class

uniform_basis_by_class <- function(Rs, K, basis_info){
  
  # create a K*K list
  basis <- vector("list",K)
  for(i in 1:length(basis)){
    basis[[i]] <- vector("list",K)
  }
  
  
  for(k_1 in 1:K){
    for(k_2 in 1:K){
      b <- Rs[k_1,k_2]
      p <- basis_info$degree[k_1,k_2]
      num_basis_funs <- basis_info$n[k_1,k_2]
      basis[[k_1]][[k_2]] <- uniform_basis(b, p, num_basis_funs)
      basis[[k_1]][[k_2]]$degree <- p
      basis[[k_1]][[k_2]]$Rmax <- Rs[k_1,k_2]
    }
  }
  
  return(basis)
}

# # test
#  ll <- uniform_basis_by_class(3,2)
#  ll
# # num_basis_fun
#  length(ll$basis$f)
 