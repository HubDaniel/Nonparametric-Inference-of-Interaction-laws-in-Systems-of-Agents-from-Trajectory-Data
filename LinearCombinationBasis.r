## LinearCombinationBasis

LinearCombinationBasis <- function(K,basis,alpha){
  print("start linear com basis")
  
  num_classes <- K
  
  
  f_ptr <- vector("list",K)
  for(i in 1:length(f_ptr)){
    f_ptr[[i]] <- vector("list",K)
  }
  
  sum_prev_num_basis <- 0
  
  # function to force evaluation of a function. prevent lazy evaluation
  h <- function(a,b){
    force(a)
    force(b)
    function(x) eval_basis_functions(x,a,b)
  }
  
  
  
  for(k_1 in 1:num_classes){
    for(k_2 in 1:num_classes){
      one_basis <- basis[[k_1]][[k_2]]
      num_basis <- length(one_basis$f)
      ind_1 <- sum_prev_num_basis + 1
      ind_2 <- sum_prev_num_basis + num_basis
      one_alphas <- alpha[ind_1 : ind_2]
 
      f_ptr[[k_1]][[k_2]] <- h(a=one_alphas,b=one_basis)


      sum_prev_num_basis <- sum_prev_num_basis + num_basis
    }
  }
  
  lastIdx = sum_prev_num_basis
  
  final_result <- list()
  
  final_result[["f_ptr"]] <- f_ptr
  final_result[["lastIdx"]] <- lastIdx
  
  print("finish linear com basis")
  
  return(final_result)
}
