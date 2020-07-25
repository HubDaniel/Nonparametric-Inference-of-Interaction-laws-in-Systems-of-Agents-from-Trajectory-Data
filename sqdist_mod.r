## Find Euclidean distance
# input: x: the data matrix, matrix of dimension d * N, where d is the size of each state vector and N is the total 
#           number of agents (observations). This is only one column of the trajectory matrix where it has L
#           columns of this d * N vectors.
# output: result: matrix of dimension N*N. For example, the [1,1] element of this matrix contains the distance
#                 between agent x1 and x1. The [1,2] element contains the distance between x1 and x2...
sqdist_mod <- function(x){
  d <- dim(x)[1]
  qn <- pn <- dim(x)[2]
  qmag <- pmag <- colSums(x * x)
  reshape_pmag <- matrix(pmag,ncol = pn,nrow = pn)
  reshape_qmag <- t(reshape_pmag)
  result <- abs(reshape_qmag+reshape_pmag-2*t(x)%*%x)
  diag(result) <- 0
  return(sqrt(result))
}