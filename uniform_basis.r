## uniform_basis

uniform_basis <- function(Rs, degree, num_basis_funs){
  
  #print("uniform basis")
  
  
  # step 1: set_up_piecewise_polynomial_info in the same function
  
  # the left end point of the support of B-spline
  polynomial_info.left_end_pt   = 0
  # the right end point of the support of B-spline
  polynomial_info.right_end_pt  = Rs
  # the degree of B-splines (starting from 0)
  polynomial_info.degree        = degree
  # the number of basis functions
  polynomial_info.num_basis     = num_basis_funs

  
  
  # step 2: construct_piecewise_polynomial_basis in the same function
  
  a <- polynomial_info.left_end_pt
  b <- polynomial_info.right_end_pt
  n <- polynomial_info.num_basis
  num_sub_inter           <- n/(degree + 1)
  num_knots               <- num_sub_inter + 1
  # 
  # print(Rs)
  # print(num_sub_inter)
  # print(num_knots)
  
  
  basis <- list()
  basis[["knots"]]        <- seq(a, b, length.out = num_knots)
  basis[["knotIdxs"]]     <- rep(0,n)
  
  #print("basis[[knotIdxs]]")
  #print(basis[["knotIdxs"]])
  
  for (ind in 1:n){
    ell                   <- ind - 1
    #print("ell")
    #print(ell)
    nn                    <- ell %% (degree+ 1)
    #print("nn")
    #print(nn)
    k                     <- (ell - nn)/(degree + 1) + 1
    #print("k")
    #print(k)
    basis[["knotIdxs"]][ind]   <- k
  }
  
  
  
  my <- function(ind,degree.=degree,basis.=basis){
    function(x){
      ell                   <- ind - 1
      nn                    <- ell %% (degree.+ 1)
      k                     <- (ell - nn)/(degree. + 1) + 1
      xspan                 <- c(basis.[["knots"]][k], basis.[["knots"]][k + 1])
      
      standard_basis(nn=nn,x,xspan=xspan)
    }
  }
 
  
  
  basis[["f"]] <- lapply(1:n,my)
  
  return(basis)
  
}

# test
# l <- uniform_basis(3,1,6)
# l
# l$f[[1]](0.1)
# l$f[[2]](c(1,2.1,2.5,4))
# l$f[[3]](c(1,2.1,2.5,4))
# l$f[[4]](c(1,2.1,2.5,4))
# l$f[[5]](c(1,2.1,2.5,4))
# l$f[[6]](c(1,2.1,2.5,4))



