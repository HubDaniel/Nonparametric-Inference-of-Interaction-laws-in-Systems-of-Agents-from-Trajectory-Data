## find pair difference
# input: x: the data matrix, matrix of dimension d * N, where d is the size of each state vector and N is the total 
#           number of agents (observations). This is only one column of the trjectory matrix where it has L
#           columns of this d * N vectors.
# Output:
#   pair_diff : pair wise different of the following form:
#               |x_1 - x_1, x_2 - x_1, \ldots, x_N - x_1|
#               |x_1 - x_2, x_2 - x_2, \ldots, x_N - x_2|
#               |   \vdots,    \vdots, \ddots,    \vdots|
#               |x_1 - x_N, x_2 - x_N, \ldots, x_N - x_N|

find_pair_diff <- function(x){
  d <- dim(x)[1]
  N <- dim(x)[2]
  x_col_vec <- as.vector(x)      # convert into a long vector of length d*N
  x_col_mat <- matrix(x_col_vec,ncol=N,nrow=d*N)
  x_row_mat <- matrix(rep(t(x),N),ncol=ncol(x),byrow=T)
  pair_diff <- x_row_mat - x_col_mat
  return(pair_diff)
}


# test
#  aaa <- matrix(1:6,ncol=3)
#  find_pair_diff(aaa)
# 
#  
# 
