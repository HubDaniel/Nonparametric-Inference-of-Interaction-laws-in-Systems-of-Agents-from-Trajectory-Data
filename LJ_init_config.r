# LJ_init_config

LJ_init_config <- function(d, N, kind){

  if(kind == 2){
    y_init = matrix(rnorm(d*N),nrow=d)
    y_init = as.vector(y_init)
  }
  return(y_init)


}