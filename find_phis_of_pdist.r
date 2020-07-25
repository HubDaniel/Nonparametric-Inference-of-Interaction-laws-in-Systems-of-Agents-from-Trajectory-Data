## Find phi(r)
# input: pair_dist matrix of dimension N*N with each element r_ii'
#        kernel_function: true kernel we use to generate data
# output: result: same N*N dimension but with each element as phi(r_ii')/N

find_phis_of_pdist <- function(pair_dist,sys_info,true_kernel=T){
  
  #print("start find_phis_of_pdist")
  
  phis_of_pdist <- matrix(0,nrow = sys_info$N,ncol = sys_info$N)
  
  phis <- sys_info$phiE
  kappa <- sys_info$kappa
  
  for(k1 in 1:sys_info$K){
    row_ind <- which(sys_info$type_info == k1)
    N_k1 <- length(row_ind)
    for(k2 in 1:sys_info$K){
      col_ind <- which(sys_info$type_info == k2)
      N_k2 <- length(col_ind)
      pair_dist_Ck1_Ck2_mat <- pair_dist[row_ind,col_ind]
      pdist_Ck1_Ck2_vec <- as.vector(pair_dist_Ck1_Ck2_mat)
      
      if(true_kernel){
        phi_of_pdist_Ck1_Ck2 <- matrix(phis[[k1]][[k2]](pdist_Ck1_Ck2_vec), nrow = N_k1, ncol = N_k2)
      }else{
        phi_of_pdist_Ck1_Ck2 <- matrix(phis[[k1]][[k2]](pdist_Ck1_Ck2_vec)$y, nrow = N_k1, ncol = N_k2)
      }
      
      # 
      # print("phi_of_pdist_Ck1_Ck2")
      # print(phi_of_pdist_Ck1_Ck2)
      
      
      
      phis_of_pdist[row_ind,col_ind] <- phi_of_pdist_Ck1_Ck2 * (kappa[k2]/N_k2)
    }
  }
  
  diag(phis_of_pdist) <- 0
  
  
  return(phis_of_pdist)
}