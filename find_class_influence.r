## find_class_influence

find_class_influence <- function(energy_basis, energy_pdist, energy_pdiff, d, num_agents, kappa,basis_info){
  
  
  #print("find class influence")
  #print(energy_basis)
  
  #print(energy_basis$f)
  
  n <- length(energy_basis$f)
  
  # print("total num of basis")
  # print(n)
  
  
  
  #print("num basis")
  #print(n)
  
  num_rows <- dim(energy_pdist)[1]
  num_cols <- dim(energy_pdist)[2]
  
  #print("num_rows")
  #print(num_rows)
  #print("num_cols")
  #print(num_cols)

  # print("dim energy_pdist")
  # print(head(energy_pdist))
  
  #print("dim energy_pdiff")
  #print(dim(energy_pdiff))
  
  energy_pdist <- as.vector(energy_pdist)
  pdist_data_sorted <- sort(energy_pdist)
  pdist_data_sorted_idxs <- order(energy_pdist)
  
  #print("pdist_data_sorted_idxs")
  #print(pdist_data_sorted_idxs)
  
  
  pdist_counts <- histc(pdist_data_sorted,energy_basis$knots)$cnt
  pdist_counts[length(pdist_counts)-1] <- pdist_counts[length(pdist_counts)-1] + 
    pdist_counts[length(pdist_counts)]
  pdist_counts <- pdist_counts[-length(pdist_counts)]
  
  
  #print("pdist_counts")
  #print(pdist_counts)
  
  populatedBins <- which(pdist_counts > 0)
  
  #print("populatedBins")
  #print(populatedBins)
  
  
  pdist_counts_cumsum <- cumsum(pdist_counts)
  
  #print("pdist_counts_cumsum")
  #print(pdist_counts_cumsum)
  
  
  class_influence <- matrix(0,nrow = dim(energy_pdiff)[1], ncol = n)
  
  
  #print("dim class influence")
  #print(dim(class_influence))
  
  #print("head(energy_pdiff)")
  #print(head(energy_pdiff))
  
  
  pdiff_data_d <- list()
  for (p in seq(from=d,to=1,by=-1)){
    pdiff_data_d_tmp <- energy_pdiff[seq(from=p,to=nrow(energy_pdiff),by=d),]
    
    #print("head(pdiff_data_d_tmp)")
    #print(head(pdiff_data_d_tmp))
    
    if(is.null(ncol(pdiff_data_d_tmp))){
      pdiff_data_d[[p]] <- pdiff_data_d_tmp[pdist_data_sorted_idxs]
    }else{
      pdiff_data_d[[p]] <- matrix(pdiff_data_d_tmp[pdist_data_sorted_idxs],ncol = ncol(pdiff_data_d_tmp))
    }
    
    
  }

  
  for (k in 1:length(populatedBins)){
    
    if (k==1){
      idxs <- 1: (pdist_counts_cumsum[populatedBins[k]])

      
    }else{
      idxs <- (pdist_counts_cumsum[populatedBins[k-1]]) : (pdist_counts_cumsum[populatedBins[k]])

    }
    f_idxs <- which(energy_basis$knotIdxs == populatedBins[k])

    
    for (kp in 1:length(f_idxs)){
      
      psi_of_r_idxs <- energy_basis$f[[f_idxs[kp]]](pdist_data_sorted[idxs])$psi * kappa/num_agents
      

      # mimic the ind2sub() function in matlab
      # [idxsI,idxsJ] = ind2sub([num_rows, num_cols],pdist_data_sorted_idxs(idxs))
      idxsI <- ((pdist_data_sorted_idxs[idxs]-1) %% num_rows) + 1
      idxsJ <- floor((pdist_data_sorted_idxs[idxs]-1) / num_rows) + 1

      for (p in 1:d){
        infl_tmp <- as.vector(psi_of_r_idxs) * as.vector(pdiff_data_d[[p]][idxs])
        
        infl_tmp_mat <- matrix(0,ncol = num_cols,nrow = num_rows)
        
        # mimic the sparse matrix part but without using sparse matrix
        for (i in 1:length(idxsI)){
          infl_tmp_mat[idxsI[i],idxsJ[i]] <- infl_tmp[i]
         
        }
        
        class_influence[seq(from=p,to=nrow(class_influence),by=d),f_idxs[kp]] <- rowSums(infl_tmp_mat)
        
      }
    }
    
  }

  return(class_influence)
}




