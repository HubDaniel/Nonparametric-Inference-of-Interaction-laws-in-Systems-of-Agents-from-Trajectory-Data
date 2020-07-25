## approximate_derivative_of_x_or_v

approximate_derivative_of_x_or_v <- function(x,v,time_vec,sys_info){
  L <- length(time_vec)
  print("approx derivative of x or v")
  # print(dim(x))
  data <- x
  
  N <- sys_info$N
  d <- sys_info$d
  K <- sys_info$K
  
  
  # print(L)
  # print(dim(data))
  delta_ts             = rep(0,length(time_vec))

  delta_ts[2 : L]      = time_vec[2 : L] - time_vec[1 : (L - 1)]
  # for the initial time, t_1 - t_2
  delta_ts[1]          = time_vec[1] - time_vec[2]
  delta_ts <- rep(delta_ts,each = N*d)
  # the difference between data
  the_diff             = matrix(0,ncol = L, nrow = d*N)
  # data_k - data_{k - 1}
  the_diff[, 2 : L]   = data[, 2 : L] - data[, 1 : (L - 1)]
  # for the data at initial time, data_1 - data_2
  the_diff[, 1]      = data[, 1] - data[, 2]
  # make it into a vector
  the_diff             = as.vector(the_diff)
  # the derivative: (data_k - data_{k - 1})/(t_k - t_{k -1}) except at the
  # initial time: (data_1 - data_2)/(t_1 - t_2)
  d_vec                = the_diff / delta_ts
  
  print("finish approx derivative of x or v")
  
  return(d_vec)
}



# test
# x <- matrix(1:8,nrow=2)
# time_vec <- c(0,0.5,1.5,3.5)
# approximate_derivative_of_x_or_v(x=x,v=c(),time_vec = time_vec,d=2,K=1,N=1)
# x <- obs[,,1]
