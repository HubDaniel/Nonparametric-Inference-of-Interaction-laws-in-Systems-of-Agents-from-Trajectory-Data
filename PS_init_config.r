# PS_init_config

PS_init_config <- function(N, type_info,kind){
  if (kind == 1){
    d                    = 2   # only for d=2
    
    y_init               = matrix(0,nrow=d,ncol=N)
    
    preys_ind            = type_info == 1
    num_preys            = sum(preys_ind)
    dist_type            = 'annulus'
    domain               = c(0.05, 0.15) * num_preys    
    y_init[, preys_ind]  = uniform_dist(d, num_preys, dist_type, domain)
    
    dist_type            = 'disk'
    domain               = 0.1
    preds_ind            = type_info == 2
    num_preds            = sum(preds_ind)
    y_init[, preds_ind]  = uniform_dist(d, num_preds, dist_type, domain)

    y_init               = as.vector(y_init)
  }
  return(y_init)
}