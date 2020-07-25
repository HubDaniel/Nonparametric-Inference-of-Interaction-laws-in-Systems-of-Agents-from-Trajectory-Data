## eval_basis_functions

eval_basis_functions <- function(x,alphas,basis){
  
  #print("eval basis function")
  # print("alphas")
  # print(alphas)
  
  D <- length(alphas)
  
  # print("D")
  # print(D)
  
  y <- rep(0,length(x))
  yp <- rep(0,length(x))
  
  # print("length f")
  # print(length(basis$f))
  # 
  # print("before for loop")
  
  for (i in 1:D){
    if(alphas[i]!=0){
      
      v <- basis$f[[i]](x)$psi
      vp <- basis$f[[i]](x)$dpsi
      if( sum(v!=0)!=0 ){
        y <- y+alphas[i]*v
      }
      if( sum(vp!=0)!=0 ){
        yp <- yp + alphas[i]*vp
      }
    }
  }

  # print("after for loop")
  
  
  result <- list()
  result[["y"]] <- y
  result[["yp"]] <- yp
  
  return(result)
}


# test
#eval_basis_functions(x=c(1,2),alphas = 1:10,basis)
# a <- 1:8
# b <- list()
# b$f <- list()
# b$f[[1]] <- function(x) x^2
# te <- function(x) eval_basis_functions(x,a,b)
# te(1:4)




