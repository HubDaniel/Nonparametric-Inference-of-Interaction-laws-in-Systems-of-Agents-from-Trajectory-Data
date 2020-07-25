## repeat a matrix n times
# input: a matrix  [a,b;c,d]
# output: a matrix repeated n times: 
#         |a  b|
#         |c  d|
# becomes (n=2)
#         |a  b|
#         |a  b|
#         |c  d|
#         |c  d|


rep_matrix <- function(mat,n){
  v <- rep(n,nrow(mat))
  result <- mat[rep(1:nrow(mat), times = v), ]
  return(result)
}

# test
# aaa <- matrix(1:4,ncol=2)
# rep_matrix(aaa,2)
