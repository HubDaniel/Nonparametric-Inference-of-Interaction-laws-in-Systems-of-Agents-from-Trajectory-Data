## estimate RhoLT

source("find_maximums.r")
source("partition_traj.r")
source("sqdist_mod.r")
source("package_rhoLT.r")
source("squareform_tovector.r")
source("restructure_histcount.r")
source("prepare_hist_items.r")


estimateRhoLT <- function(obs_data, sys_info, obs_info){
  
  print("estimateRhoLT start...")
  
  obs_data <- obs_data$traj  # this is a list of length M. We need to transform this into a 3-D array
  obs_data <- array(unlist(obs_data),dim=c(sys_info$d*sys_info$N,obs_info$L,obs_info$M))  
  # 3rd dim of list is each traj
  
  max_rs <- array(0,dim=c(sys_info$K,sys_info$K,dim(obs_data)[3]))   
  # to store the max and min value of obs_data to 0
  min_rs <- array(0,dim=c(sys_info$K,sys_info$K,dim(obs_data)[3]))
  
  # step 1: find the maximum of all data (maximum dist to 0) in M trajectory during L times
  print("step 1")
  for( m in 1:dim(obs_data)[3] ){
    traj <- obs_data[,,m]
    output <- find_maximums(traj,sys_info)
    max_rs[,,m] <- output
  }
  max_rs1 <- apply(max_rs,c(1,2),max)
  
  
  # print("max_rs1 on all traj")
  # print(max_rs1)
  
  # step 2: prepare the bins for hist count using above maximum as right end point
  # prepare prepare_hist_items.m in the same function
  
  
  print("step 2")
  result <- prepare_hist_items(sys_info$K, obs_info$hist_num_bins, dim(obs_data)[3], max_rs1)
  
  histedgesR <- result$histedgesR       # list of dim K*K: edges of each type
  histbinwidthR <- result$histbinwidthR # matrix of dim K*K: bin width of each type
  histcountR <- result$histcountR       # empty list of length M
  
  
  # step 3: find the max and min pairwise distance
  print("step 3")
  max_rs2 <- array(0,dim=c(sys_info$K,sys_info$K,dim(obs_data)[3]))   
  # to store the max and min value of pairwise difference
  min_rs2 <- array(0,dim=c(sys_info$K,sys_info$K,dim(obs_data)[3]))
  
  for (m in 1:obs_info$M){
    
    print(paste0("m=",m,"/",obs_info$M))
    
    traj <- obs_data[,,m]
    
    # print("traj")
    # print(dim(traj))
    
    pdist_out <- partition_traj(traj,sys_info)
    # pdist_out is a list of 3. max_r, min_r and pdist_x
    # max_r and min_r are matrix of K*K, with each element as the value of max or min between type ck1 & ck2 
    # across all times
    # pdist_x is a list of K*K, with each element containing pairwise diff. between type ck1 & ck2
    
    
    # print("pdist_out$pdist_x")
    # print(dim(pdist_out$pdist_x))
    
    max_rs2[,,m] <- pdist_out$max_r
    min_rs2[,,m] <- pdist_out$min_r
    
    # print("max_rs2 vec")
    # print(max_rs2)
    # print("min_rs2 vec")
    # print(min_rs2)

    pdist_x <- pdist_out$pdist_x
    
    # print("pdist_x")
    # print(dim(pdist_x))
    
    # print("pdist_x_vec")
    # print(length(as.vector(pdist_x)))
    
    histcountR_m <- vector("list",sys_info$K)
    
    for(k1 in 1:sys_info$K){
      for(k2 in 1:sys_info$K){
        pdist_x_Ck1_Ck2 <- pdist_x[[k1]][[k2]]
        if(!is.null(pdist_x_Ck1_Ck2)){
          histcountR_m[[k1]][[k2]] <- 
            histc(as.vector(pdist_x_Ck1_Ck2), histedgesR[[k1]][[k2]])$cnt
          
          histcountR_m[[k1]][[k2]][length(histcountR_m[[k1]][[k2]])-1] <- 
            histcountR_m[[k1]][[k2]][length(histcountR_m[[k1]][[k2]])-1] + 
            histcountR_m[[k1]][[k2]][length(histcountR_m[[k1]][[k2]])]
          
          histcountR_m[[k1]][[k2]] <- histcountR_m[[k1]][[k2]][-length(histcountR_m[[k1]][[k2]])]
          
        }
      }
    }
   
    
    
    # print("histcountR_m")
    # print(length(histcountR_m))
    
    histcountR[[m]] <- histcountR_m  
    # this is a list of length M, each is a list of K*K, 
    # [[k_1]][[k_2]] is the histcount of the pariwise distance between type k1 and k2
    
    # print("histcountR[[m]]")
    # print(length(histcountR[[m]]))
  }
  
  # step 4: package the data
  print("step 4")
  rhoLT <- package_rhoLT(histedgesR, histcountR, histbinwidthR, obs_info$M, sys_info, obs_info, max_rs2, min_rs2)
  
  return(rhoLT)
  
  
}


