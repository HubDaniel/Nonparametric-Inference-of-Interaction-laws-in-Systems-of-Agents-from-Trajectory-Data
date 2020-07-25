## construct standard basis (works for degree >=0)


standard_basis <- function(nn, x, xspan){
  # find the interval where the function is supported: [x_k, x_kp1]
  
  #print("xspan")
  #print(xspan)
  
  x_k         = xspan[1]
  x_kp1       = xspan[2]

  # only evaluate the points within this interval
  ind         = which((x_k <= x) & (x < x_kp1))     
  
  # initiliaze the quantities
  psi         = rep(0,length(x))
  dpsi        = rep(0,length(x))
  
  # the function value is x^n
  psi[ind]    = x[ind]^nn
  
  # the derivative is n * x^(n - 1) when n > 0
  if (nn > 0){
    dpsi[ind] = nn * x[ind]^(nn - 1)
  }
  
  result <- list()
  result[["psi"]] <- psi
  result[["dpsi"]] <- dpsi
  return(result)
}

# test
# standard_basis(1,c(0.5,1.2,2.2),c(0,2))














