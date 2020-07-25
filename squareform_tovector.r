# squareform_tovector

squareform_tovector <- function(Y){
  n <- dim(Y)[2]
  Z <- Y[ lower.tri(matrix(TRUE,nrow = n, ncol = n)) ]
  return(Z)
}

# test
# Y = matrix(1:9,ncol=3)
# Y
# squareform_tovector(Y)
# Y[lower.tri(Y)]
# 
# Y = matrix(1:6,ncol=3)
# Y
# squareform_tovector(Y)
