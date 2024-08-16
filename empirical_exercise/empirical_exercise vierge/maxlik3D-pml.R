pseudo.log.lik3D.super      <- function(c) {
  c11 <- c[1]
  c21 <- c[2]
  c31 <- c[3]
  c12 <- c[4]
  c22 <- c[5]
  c32 <- c[6]
  c13 <- c[7]
  c23 <- c[8]
  c33 <- c[9]
  
  c.1y <- c11*tu_chol[,1] + c21*tu_chol[,2] + c31*tu_chol[,3]
  c.2y <- c12*tu_chol[,1] + c22*tu_chol[,2] + c32*tu_chol[,3]
  c.3y <- c13*tu_chol[,1] + c23*tu_chol[,2] + c33*tu_chol[,3]
  
  #P <- 1/(2 * p^(1/p - 1))*gamma(1/p)
  g <- rep(NA,1)
  g[1] <- -1*(n*log(1/sqrt(2*pi))-sum((c.1y)^2/2) + 
                n*log(gamma((5+1)/2)/(sqrt(5*pi)*gamma(5/2)))-((5+1)/2)*sum(log(1+(c.2y)^2/5)) + 
                n*log(gamma((12+1)/2)/(sqrt(12*pi)*gamma(12/2)))-((12+1)/2)*sum(log(1+(c.3y)^2/12)))
  g
}
pseudo.log.lik3D.sub        <- function(c) {
  c11 <- c[1]
  c21 <- c[2]
  c31 <- c[3]
  c12 <- c[4]
  c22 <- c[5]
  c32 <- c[6]
  c13 <- c[7]
  c23 <- c[8]
  c33 <- c[9]
  
  c.1y <- c11*tu_chol[,1] + c21*tu_chol[,2] + c31*tu_chol[,3]
  c.2y <- c12*tu_chol[,1] + c22*tu_chol[,2] + c32*tu_chol[,3]
  c.3y <- c13*tu_chol[,1] + c23*tu_chol[,2] + c33*tu_chol[,3]
  
  #P <- 1/(2 * p^(1/p - 1))*gamma(1/p)
  g <- rep(NA,1)
  g[1] <- -1*(n*log(1/sqrt(2*pi)) + sum(pi*(c.1y)^2+log(cosh((pi/2)*c.1y))) + 
                n*log(1/sqrt(2*pi)) + sum(pi*(c.2y)^2+log(cosh((pi/2)*c.2y))) + 
                n*log(1/sqrt(2*pi)) + sum(pi*(c.3y)^2+log(cosh((pi/2)*c.3y))))
  g
}
orth.constr3D               <- function(c){
  c11 <- c[1]
  c21 <- c[2]
  c31 <- c[3]
  c12 <- c[4]
  c22 <- c[5]
  c32 <- c[6]
  c13 <- c[7]
  c23 <- c[8]
  c33 <- c[9]
  
  g <- rep(NA,6)
  g[1] <- c11*c12 + c21*c22 + c31*c32 # c.1*c.2
  g[2] <- c11*c13 + c21*c23 + c31*c33 # c.1*c.3
  g[3] <- c12*c13 + c22*c23 + c32*c33 # c.2*c.3
  g[4] <- c11^2 + c21^2 + c31^2 - 1         # c.1*c.1
  g[5] <- c12^2 + c22^2 + c32^2 - 1         # c.2*c.2
  g[6] <- c13^2 + c23^2 + c33^2 - 1         # c.3*c.3
  g
}