# prepare_hist_items

prepare_hist_items <- function(K, hist_num_bins, M, max_rs){
  
  histedgesR <- vector("list",K)
  
  for(i in 1:length(histedgesR)){
    histedgesR[[i]] <- vector("list",K)
  }
  
  histbinwidthR <- matrix(0,ncol=K,nrow=K)
  
  for(k1 in 1:K){
    for(k2 in 1:K){
      histedgesR[[k1]][[k2]] <- seq(from=0,to=max_rs[k1,k2],length.out = hist_num_bins+1)
      histbinwidthR[k1,k2] <- max_rs[k1,k2]/hist_num_bins
    }
  }
  
  histcountR <- vector("list",M)
  
  result <- list()
  result$histedgesR <- histedgesR
  result$histbinwidthR <- histbinwidthR
  result$histcountR <- histcountR
  
  return(result)
  
}