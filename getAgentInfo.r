# getAgentInfo

getAgentInfo <- function( sys_info ){
  agents_info <- list()
  agents_info$idxs <- vector("list",sys_info$K)
  agents_info$num_agents <- rep(0,sys_info$K)
  for(k in 1:sys_info$K){
    agents_info$idxs[[k]] <- which(sys_info$type_info == k)
    agents_info$num_agents[k] <- length(agents_info$idxs[[k]])
  }
  
  return(agents_info)
}