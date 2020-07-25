## package_rhoLT

package_rhoLT <- function(histedgesR, histcountR, histbinwidthR, M, sys_info, obs_info,  max_rs, min_rs){
  
  
  
  
  # step 1:
  # restructure_histcount(histcountR, M, K, hist_num_bins) in the same function
  print("package rhoLT -- step 1")
  output <- restructure_histcount(histcountR, M, sys_info$K, sys_info$type_info, obs_info$hist_num_bins)
  
  histcountR <- output
  # array with 4 dimensions. last dim indicates traj, 1st & 2nd indicate type between k1,k2
  # 3rd dim is used to store histcounts between type k1 and k2
  
  
  # step 2:
  print("package rhoLT -- step 2")
  max_rs <- apply(max_rs,c(1,2),max) # max pairwise distance between type k1,k2 for all traj across all time
  min_rs <- apply(min_rs,c(1,2),min) # min pairwise distance between type k1,k2 for all traj across all time
  histcountR <- apply(histcountR,c(1,2,3),sum) # count in each bin between type k1,k2 across all traj

  histcount_R <- vector("list",sys_info$K)
  for(i in 1:length(histcount_R)){
    histcount_R[[i]] <- vector("list",sys_info$K)
  }
  
  hist_R <- vector("list",sys_info$K)
  for(i in 1:length(hist_R)){
    hist_R[[i]] <- vector("list",sys_info$K)
  }
  
  supp_R <- vector("list",sys_info$K)
  for(i in 1:length(supp_R)){
    supp_R[[i]] <- vector("list",sys_info$K)
  }
  
  
  for(k1 in 1:sys_info$K){
    for(k2 in 1:sys_info$K){
      supp_R[[k1]][[k2]] <- c(min_rs[k1,k2],max_rs[k1,k2]) # support between type k1,k2
      histcount_R[[k1]][[k2]] <- histcountR[k1, k2,]       # hist counts between type k1,k2
      hist_R[[k1]][[k2]] <- histcount_R[[k1]][[k2]]/(sum(histcount_R[[k1]][[k2]]) * histbinwidthR[k1, k2])
      # height of histograms to be plot between type k1,k2
    }
  }
  
  
  

  
  result <- list()
  result[["supp"]] <- supp_R            # list of K*K, support between type k1,k2
  result[["histcount"]] <- histcount_R  # list of K*K, hist counts between type k1,k2
  result[["hist"]] <- hist_R            # list of K*K, height of histograms to be plot between type k1,k2
  result[["histedgesR"]] <- histedgesR  # list of K*K, hist edges between type k1,k2
  result[["histbinwidth"]] <- histbinwidthR
  
  return(result)
}