# uniform_dist

uniform_dist <- function(d, N, type, domain){
  if(d==2){
    if(type == 'annulus'){
      R_1              = domain[1] 
      R_2              = domain[2] 
      theta            = 2 * pi * runif(N, min=0, max=1) 
      r                = sqrt((R_2^2 - R_1^2) * runif(N, min=0, max=1)  + R_1^2) 
      point_dist       = matrix(0,nrow=2,ncol=N)
      point_dist[1,]   = r * cos(theta) 
      point_dist[2,]   = r * sin(theta) 
    }else if(type == 'disk'){
      R                = domain[1]
      theta            = 2 * pi * runif(N, min=0, max=1)
      r                = R * sqrt(runif(N, min=0, max=1))
      point_dist       = matrix(0,nrow=2,ncol=N)
      point_dist[1,]   = r * cos(theta)
      point_dist[2,]   = r * sin(theta)
    }
  }
  return(point_dist)
}