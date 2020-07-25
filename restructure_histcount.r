# restructure_histcount

restructure_histcount <- function(HCR, M, K, type_info, num_bins){
  
  new_HCR <- array(0,dim=c(K,K,num_bins,M))
  
  for(m in 1:M){
    HCR_m <- HCR[[m]]    
    # a list of K*K, [[k_1]][[k_2]] is the histcount of the pariwise distance between type k1 and k2 
    
    
    for(k1 in 1:K){
      num_Ck1 <- sum(type_info==k1)
      
      for(k2 in 1:K){
        if(k1==k2){
          if(num_Ck1==1){
            to_update <- FALSE
          }else{
            to_update <- TRUE
          }
        }else{
          to_update <- TRUE
        }
        
        if(to_update){
          new_HCR[k1,k2,,m] <- HCR_m[[k1]][[k2]] # histcount info is stored through dimension 3
        }
      }
    }
  }
  
  return(new_HCR)   
  # array with 4 dimensions. last dim indicates traj, 1st & 2nd indicate type between k1,k2
  # 3rd dim is used to store histcounts between type k1 and k2
  
}