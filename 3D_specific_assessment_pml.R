## General assessment PML
## Load packege here for managing workflow

if(length(which(installed.packages() %in% c("rstudioapi") == T)) == 0){
  install.packages("rstudioapi")
  library(rstudioapi)
}else{library(rstudioapi)}

if(length(which(installed.packages() %in% c("here") == T)) == 0){
  install.packages("here")
  library(here)
}else{library(here)}

## Open ica-svars-comparative-analysis.Rproj from window ####
if(!base::grepl(x = here(),pattern = "ica-svars-comparative-analysis")){
  file_path <- file.choose()
}

if(!base::grepl(x = here(),pattern = "ica-svars-comparative-analysis")){
  openProject(file_path)
}

## Check that working directory contains the ica-svars-comparative-analysis.Rproj file
library(here)
dr_here()
here()
## Independent Component Analysis with different methods (MC analysis)
source(paste0(here(),"/Rpackages.R"),local = T)

## File with variables definitions and parameter settings

n              <- 400 # sample size
choose_cluster <- 3
## Number of shape parameter values (in the specific assessment and for inference table we only focuse on 4 distributional scenario - See paper)
seq1          <- seq(from = 0.5, to = 3.5, length.out = 15) ## Sequence of p-Generalized parameter for general evaluation exercise (1st part)
seq2          <- seq(from = 4, to = 100, length.out = 5) ## Sequence of p-Generalized parameter for general evaluation exercise (2nd part)
SEQ           <- as.numeric(c(seq1, seq2))

source(here("confg_file.R"),local = T)

#### Codes to be updated and experiment setting ##
source(here("functions_to_load.R"),local = T)


## Number of initialization Hypercybe sampling
if (n_lhs == 1) {
  lhs <- matrix(data = c(0,0),nrow = 1)
}else{
  set.seed(55)
  lhs            <- 2*pi*improvedLHS(n = n_lhs, k = 3)
}

for(z in 1:length(SEQ)){
  print(paste0("experiment ", z, " - first MC ", Sys.time()))
  print(paste0("sample size = ",n," initialization =  ", n_lhs, " MC = ",mc))
  ## First MC
  i <- 1
  set.seed(11)
  e1     <- stdEpsrnd(type = 'general',dof = 5, n = n, p = SEQ[z])
  set.seed(22)
  e2     <- stdEpsrnd(type = 'general',dof = 5, n = n, p = SEQ[z])
  set.seed(33)
  e3     <- stdEpsrnd(type = 'general',dof = 5, n = n, p = SEQ[z])
  
  E      <- cbind(e1,e2,e3)
  Atrue            <- Btrue3D
  U                <- Atrue %*% t(E)
  sigg             <- cov(t(U))
  C[[i]]           <- t(chol(sigg))
  u_chol[[i]]      <- t(solve(C[[i]]) %*% (U))
  Orth.true[[i]]   <- solve(C[[i]]) %*% Atrue
  
  ## Run methods and record best initialization matrix
  for (kk in 1:nrow(lhs)){
    ## Initialization
    winitx             <- matrix(data = c(1, 0, 0,
                                          0, cos(lhs[kk,1]), -sin(lhs[kk,1]),
                                          0, sin(lhs[kk,1]),  cos(lhs[kk,1])), ncol = 3, byrow = T)
    winity             <- matrix(data = c(cos(lhs[kk,2]), 0, -sin(lhs[kk,2]), 
                                          0,              1,  0,
                                          sin(lhs[kk,2]), 0,  cos(lhs[kk,2])), ncol = 3, byrow = T)
    winitz             <- matrix(data = c(cos(lhs[kk,3]), -sin(lhs[kk,3]), 0,
                                          sin(lhs[kk,3]),  cos(lhs[kk,3]), 0,
                                          0,               0,              1), ncol = 3, byrow = T)
    ## Store initial conditions
    
    orth_matrix_init  <- winitx %*% winity %*% winitz
    p.start[[kk]]           <- as.vector(C[[i]] %*% orth_matrix_init)
    ## Constrained optimization with more parameters than equations
    if (kurtosis(e1)>3 & kurtosis(e2)>3 | kurtosis(e1)>3 & kurtosis(e3)>3 | kurtosis(e2)>3 & kurtosis(e3)>3){
      ## super-gaussian case
      please           <- auglag(par = p.start[[kk]], fn = pseudo.log.lik3D.super, heq = orth.constr3D,
                                 control.outer = list(trace=F))$par 
      
    }else{
      ##sub-gaussian case
      please           <- auglag(par = p.start[[kk]], fn = pseudo.log.lik3D.sub, heq = orth.constr3D,
                                 control.outer = list(trace=F))$par 
    }
    Orth.hat           <- matrix(data = please,nrow = 3,byrow = F)
    B.init[[kk]]       <- C[[i]] %*% Orth.hat
    ilm.index[kk]      <- MD(W.hat = solve(B.init[[kk]]),A = Atrue)
  }
  ## Choose the best initialization
  best_init        <- p.start[[which.min(ilm.index)]]
  FR_pml_nt[i]     <- ilm.index[[which.min(ilm.index)]]
  A.pml[[i]]       <- B.init[[which.min(ilm.index)]]
  
  A.ID.pml[[i]]        <- maxfinder(A = A.pml[[i]])$A.id
  Orth.max[[i]]        <- solve(C[[i]]) %*% A.ID.pml[[i]]
  
  
  ## Statistics based on MLE
  ## Matrix A1 of Gorieroux with PML procedure, for both gaussian and sub-gaussian
  if (kurtosis(e1)>3 & kurtosis(e2)>3 | kurtosis(e1)>3 & kurtosis(e3)>3 | kurtosis(e2)>3 & kurtosis(e3)>3){
    ## super-gaussian case
    ## Matrix A1
    a12 <- (mean((-1)*(-1)) + mean((e2)*(-1*((e2)*(1+5)/(5-2+e2^2))))) * t(Orth.max[[i]][,2])
    a21 <- (mean((1+5)*((5-2-e2^2)/(5-2+e2^2)^2)) + mean((-e1)*e1)) *t(Orth.max[[i]][,1])
    a13 <- (mean((-1)*(-1)) + mean((e3)*(-((e3)*(1+12)/(12-2+e3^2))))) *t(Orth.max[[i]][,3])
    a31 <- (mean((1+12)*((12-2-e3^2)/(12-2+e3^2)^2)) + mean((-e1)*(e1))) * t(Orth.max[[i]][,1])
    a23 <- (mean((1+5)* ((5-2-e2^2)/(5-2+e2^2)^2)) + mean((e3)*(-((e3)*(1+12)/(12-2+e3^2)))))*t(Orth.max[[i]][,3])
    a32 <- (mean((1+12)*((12-2-e3^2)/(12-2+e3^2)^2)) + mean((e2)*(-((e2)*(1+5)/(5-2+e2^2)))))*t(Orth.max[[i]][,2])
    
    
    ## Matrix Omega
    om.1212 <- mean((-e1)^2) + mean(((-e2)*(1+5)/(5-2+e2^2))^2) - (2*mean((-e1)*e1)*mean(e2*((-e2)*(1+5)/(5-2+e2^2))))
    om.1313 <- mean((-e1)^2) + mean(((-e3)*(1+12)/(12-2+e3^2))^2) - (2*mean((-e1)*e1)*mean(e3*((-e3)*(1+12)/(12-2+e3^2))))
    om.2323 <- mean(((-e2)*(1+5)/(5-2+e2^2))^2) + mean(((-e3)*(1+12)/(12-2+e3^2))^2)-
      (2*mean((e2)*((-e2)*(1+5)/(5-2+e2^2)))*mean((e3)*((-e3)*(1+12)/(12-2+e3^2))))
    om.1213 <- mean((-e2)*(1+5)/(5-2+e2^2))*mean((-e3)*(1+12)/(12-2+e3^2))
    om.1223 <- -mean(-e1)*mean((-e3)*(1+12)/(12-2+e3^2))
    om.1312 <- om.1213
    om.1323 <- mean(-e1)*mean((-e2)*(1+5)/(5-2+e2^2))
    om.2312 <- -mean(((-e3)*(1+12)/(12-2+e3^2)))*mean(-e1)
    om.2313 <- mean((-e2)*(1+5)/(5-2+e2^2))*mean(-e1)
    
    ## Condition for local concavity
    #cond1[i]   <- mean(-1-e1*(-1))
    #cond2[i]   <- mean((-(1+5)*((5-2-e2^2)/(5-2+e2^2)^2))-((e2)*(-(e2)*(1+5)/(5-2+e2^2))))
    #cond3[i]   <- mean((-(1+12)*((12-2-e2^2)/(12-2+e2^2)^2))-((e2)*(-(e2)*(1+12)/(12-2+e2^2))))
  }else{
    ## sub-gaussian case
    ## MAtrix A1
    dl.g.sub <- function(e){
      2*pi*e + (pi/2)*(tanh(e*(pi/2)))
    }
    d2l.g.sub <- function(e){
      2*pi + ((pi/2)^2)*(1-tanh(e*(pi/2))^2)
    }
    ## MATRIX A1 for SUB GAUSSIAN
    a12 <- (mean((-1)*d2l.g.sub(e1)) + mean((e2)*dl.g.sub(e2))) * t(Orth.max[[i]][,2])
    a21 <- (mean((-1)*d2l.g.sub(e2)) + mean((e1)*dl.g.sub(e1))) * t(Orth.max[[i]][,1])
    a13 <- (mean((-1)*d2l.g.sub(e1)) + mean((e3)*dl.g.sub(e3))) * t(Orth.max[[i]][,3])
    a31 <- (mean((-1)*d2l.g.sub(e3)) + mean((e1)*dl.g.sub(e1))) * t(Orth.max[[i]][,1])
    a23 <- (mean((-1)*d2l.g.sub(e2)) + mean((e3)*dl.g.sub(e3))) * t(Orth.max[[i]][,3])
    a32 <- (mean((-1)*d2l.g.sub(e3)) + mean((e2)*dl.g.sub(e2))) * t(Orth.max[[i]][,2])
    
    ## Matrix Omega
    om.1212 <- mean(dl.g.sub(e1)^2) + mean(dl.g.sub(e2)^2) - (2*mean((e1)*dl.g.sub(e1))*mean((e2)*dl.g.sub(e2)))
    om.1313 <- mean(dl.g.sub(e1)^2) + mean(dl.g.sub(e3)^2) - (2*mean((e1)*dl.g.sub(e1))*mean((e3)*dl.g.sub(e3)))
    om.2323 <- mean(dl.g.sub(e2)^2) + mean(dl.g.sub(e3)^2) - (2*mean((e2)*dl.g.sub(e2))*mean((e3)*dl.g.sub(e3)))
    om.1213 <- mean(dl.g.sub(e2))*mean(dl.g.sub(e3))
    om.1223 <- -mean(dl.g.sub(e1))*mean(dl.g.sub(e3))
    om.1312 <- om.1213
    om.1323 <- mean(dl.g.sub(e1))*mean(dl.g.sub(e2))
    om.2312 <- -mean(dl.g.sub(e3))*mean(dl.g.sub(e1))
    om.2313 <- mean(dl.g.sub(e2))*mean(dl.g.sub(e1))
  }
  ## Matrix A1
  A1 <- matrix(data = c(a12,-a21,0,0,0,
                        a13,0,0,0,-a31,
                        0,0,0,a23,-a32),ncol = 9,byrow = T)
  
  ## Matrix A2
  A2 <- matrix(data = c(Orth.max[[i]][,2],Orth.max[[i]][,1],0,0,0,
                        Orth.max[[i]][,3],0,0,0,Orth.max[[i]][,1],
                        0,0,0,Orth.max[[i]][,3],Orth.max[[i]][,2]), ncol = 9, byrow = T)
  
  ## MAtrix A3
  A3 <- matrix(data = c(Orth.max[[i]][,1],0,0,0,0,0,0,
                        0,0,0,Orth.max[[i]][,2],0,0,0,
                        0,0,0,0,0,0,Orth.max[[i]][,3]),ncol = 9,byrow = T)
  
  A <- matrix(data = rbind(A1,A2,A3),ncol = 9)
  
  ## Matrix OMEGA blocks
  om.hat <- matrix(data = c(om.1212, om.1213, om.1223,
                            om.1312, om.1313, om.1323,
                            om.2312, om.2313, om.2323),ncol = 3,byrow = T)
  om.hat.1 <- solve(om.hat)
  om2      <- matrix(data=0,ncol = 6,nrow = 3)
  om.first.row <- cbind(om.hat.1,om2)
  OMEGA1   <- rbind(om.first.row,matrix(data=0,nrow = 6,ncol = 9))
  OMEGA    <- matrix(data = c(om.1212, om.1213, om.1223,0,0,0,0,0,0, 
                              om.1312, om.1313, om.1323,0,0,0,0,0,0,
                              om.2312, om.2313, om.2323,0,0,0,0,0,0,
                              0,0,0,0,0,0,0,0,0,
                              0,0,0,0,0,0,0,0,0,
                              0,0,0,0,0,0,0,0,0,
                              0,0,0,0,0,0,0,0,0,
                              0,0,0,0,0,0,0,0,0,
                              0,0,0,0,0,0,0,0,0),ncol = 9,byrow = T)
  var.matrix[[i]] <- solve(A) %*% OMEGA %*% t(solve(A))
  
  ## t-statistic
  ## This should be the the chisq
  #t.stat[i] <-n*t((vec(Bmax[[i]])-vec(Atrue))) %*% t(A) %*% OMEGA1 %*% A %*% (vec(Bmax[[i]])-vec(Atrue))
  gamma.hat <- as.vector(Orth.max[[i]])
  gamma.0   <- as.vector(Orth.true[[1]])
  t.stat[i] <- n*((t(gamma.hat - gamma.0) %*% t(A) %*% OMEGA1 %*% A %*% (gamma.hat - gamma.0)))
  
  
  first_mc_res <-  list(statistics      = t.stat[i],
                        Choleski        = C[[i]],
                        U               = U,
                        pml_init_mixing = B.init,
                        Ort.est         = Orth.max[[i]],
                        Ort.true        = Orth.true[[i]],
                        B.est           = A.ID.pml[[i]],
                        variance        = var.matrix[[i]],
                        perfm_plm       = FR_pml_nt[i])
  
  rm(i)
  
  ## Run parallel MCs ####
  nr_cluster <- choose_cluster
  cl <- makeCluster(nr_cluster)
  registerDoParallel(cl)
  set.seed(123+z)
  ciao <- foreach(i = 2:mc,
                  .packages = c('pgnorm', 'normtest', 'steadyICA', "ProDenICA",
                                'fastICA', 'JADE', 'tidyverse', 'magrittr', 'gtools', 
                                'DEoptim', "matrixcalc", "copula","moments",
                                "alabama")) %dorng% {
                                  
                                  Sys.sleep(1)  
                                  
                                  e1     <- stdEpsrnd(type = 'general',dof = 5, n = n, p = SEQ[z])
                                  e2     <- stdEpsrnd(type = 'general',dof = 5, n = n, p = SEQ[z])
                                  e3     <- stdEpsrnd(type = 'general',dof = 5, n = n, p = SEQ[z])
                                  
                                  E      <- cbind(e1,e2,e3)
                                  Atrue            <- Btrue3D
                                  U                <- Atrue %*% t(E)
                                  sigg             <- cov(t(U))
                                  C[[i]]           <- t(chol(sigg))
                                  u_chol[[i]]      <- t(solve(C[[i]]) %*% (U))
                                  Orth.true[[i]]   <- solve(C[[i]]) %*% Atrue
                                  
                                  if (kurtosis(e1)>3 & kurtosis(e2)>3 | kurtosis(e1)>3 & kurtosis(e3)>3 | kurtosis(e2)>3 & kurtosis(e3)>3){
                                    ## super-gaussian case
                                    please           <- auglag(par = best_init, fn = pseudo.log.lik3D.super, heq = orth.constr3D,
                                                               control.outer = list(trace=F))$par
                                  } else{
                                    ## sub-gaussian case
                                    please           <- auglag(par = best_init, fn = pseudo.log.lik3D.sub, heq = orth.constr3D,
                                                               control.outer = list(trace=F))$par
                                  }
                                  
                                  
                                  Orth.hat           <- matrix(data = please,nrow = 3,byrow = F)
                                  A.pml[[i]]         <- C[[i]] %*% Orth.hat
                                  FR_pml_nt[i]       <- MD(W.hat = solve(A.pml[[i]]),A = Atrue)
                                  
                                  
                                  A.ID.pml[[i]]        <- maxfinder(A = A.pml[[i]])$A.id
                                  Orth.max[[i]]        <- solve(C[[i]]) %*% A.ID.pml[[i]]
                                  
                                  ## Statistics based on MLE
                                  ## Matrix A1 of Gorieroux with PML procedure, for both gaussian and sub-gaussian
                                  if (kurtosis(e1)>3 & kurtosis(e2)>3 | kurtosis(e1)>3 & kurtosis(e3)>3 | kurtosis(e2)>3 & kurtosis(e3)>3){
                                    ## super-gaussian case
                                    ## Matrix A1
                                    a12 <- (mean((-1)*(-1)) + mean((e2)*(-1*((e2)*(1+5)/(5-2+e2^2))))) * t(Orth.max[[i]][,2])
                                    a21 <- (mean((1+5)*((5-2-e2^2)/(5-2+e2^2)^2)) + mean((-e1)*e1)) *t(Orth.max[[i]][,1])
                                    a13 <- (mean((-1)*(-1)) + mean((e3)*(-((e3)*(1+12)/(12-2+e3^2))))) *t(Orth.max[[i]][,3])
                                    a31 <- (mean((1+12)*((12-2-e3^2)/(12-2+e3^2)^2)) + mean((-e1)*(e1))) * t(Orth.max[[i]][,1])
                                    a23 <- (mean((1+5)* ((5-2-e2^2)/(5-2+e2^2)^2)) + mean((e3)*(-((e3)*(1+12)/(12-2+e3^2)))))*t(Orth.max[[i]][,3])
                                    a32 <- (mean((1+12)*((12-2-e3^2)/(12-2+e3^2)^2)) + mean((e2)*(-((e2)*(1+5)/(5-2+e2^2)))))*t(Orth.max[[i]][,2])
                                    
                                    
                                    ## Matrix Omega
                                    om.1212 <- mean((-e1)^2) + mean(((-e2)*(1+5)/(5-2+e2^2))^2) - (2*mean((-e1)*e1)*mean(e2*((-e2)*(1+5)/(5-2+e2^2))))
                                    om.1313 <- mean((-e1)^2) + mean(((-e3)*(1+12)/(12-2+e3^2))^2) - (2*mean((-e1)*e1)*mean(e3*((-e3)*(1+12)/(12-2+e3^2))))
                                    om.2323 <- mean(((-e2)*(1+5)/(5-2+e2^2))^2) + mean(((-e3)*(1+12)/(12-2+e3^2))^2)-
                                      (2*mean((e2)*((-e2)*(1+5)/(5-2+e2^2)))*mean((e3)*((-e3)*(1+12)/(12-2+e3^2))))
                                    om.1213 <- mean((-e2)*(1+5)/(5-2+e2^2))*mean((-e3)*(1+12)/(12-2+e3^2))
                                    om.1223 <- -mean(-e1)*mean((-e3)*(1+12)/(12-2+e3^2))
                                    om.1312 <- om.1213
                                    om.1323 <- mean(-e1)*mean((-e2)*(1+5)/(5-2+e2^2))
                                    om.2312 <- -mean(((-e3)*(1+12)/(12-2+e3^2)))*mean(-e1)
                                    om.2313 <- mean((-e2)*(1+5)/(5-2+e2^2))*mean(-e1)
                                    
                                    ## Condition for local concavity
                                    #cond1[i]   <- mean(-1-e1*(-1))
                                    #cond2[i]   <- mean((-(1+5)*((5-2-e2^2)/(5-2+e2^2)^2))-((e2)*(-(e2)*(1+5)/(5-2+e2^2))))
                                    #cond3[i]   <- mean((-(1+12)*((12-2-e2^2)/(12-2+e2^2)^2))-((e2)*(-(e2)*(1+12)/(12-2+e2^2))))
                                  }else{
                                    ## sub-gaussian case
                                    ## MAtrix A1
                                    dl.g.sub <- function(e){
                                      2*pi*e + (pi/2)*(tanh(e*(pi/2)))
                                    }
                                    d2l.g.sub <- function(e){
                                      2*pi + ((pi/2)^2)*(1-tanh(e*(pi/2))^2)
                                    }
                                    ## MATRIX A1 for SUB GAUSSIAN
                                    a12 <- (mean((-1)*d2l.g.sub(e1)) + mean((e2)*dl.g.sub(e2))) * t(Orth.max[[i]][,2])
                                    a21 <- (mean((-1)*d2l.g.sub(e2)) + mean((e1)*dl.g.sub(e1))) * t(Orth.max[[i]][,1])
                                    a13 <- (mean((-1)*d2l.g.sub(e1)) + mean((e3)*dl.g.sub(e3))) * t(Orth.max[[i]][,3])
                                    a31 <- (mean((-1)*d2l.g.sub(e3)) + mean((e1)*dl.g.sub(e1))) * t(Orth.max[[i]][,1])
                                    a23 <- (mean((-1)*d2l.g.sub(e2)) + mean((e3)*dl.g.sub(e3))) * t(Orth.max[[i]][,3])
                                    a32 <- (mean((-1)*d2l.g.sub(e3)) + mean((e2)*dl.g.sub(e2))) * t(Orth.max[[i]][,2])
                                    
                                    ## Matrix Omega
                                    om.1212 <- mean(dl.g.sub(e1)^2) + mean(dl.g.sub(e2)^2) - (2*mean((e1)*dl.g.sub(e1))*mean((e2)*dl.g.sub(e2)))
                                    om.1313 <- mean(dl.g.sub(e1)^2) + mean(dl.g.sub(e3)^2) - (2*mean((e1)*dl.g.sub(e1))*mean((e3)*dl.g.sub(e3)))
                                    om.2323 <- mean(dl.g.sub(e2)^2) + mean(dl.g.sub(e3)^2) - (2*mean((e2)*dl.g.sub(e2))*mean((e3)*dl.g.sub(e3)))
                                    om.1213 <- mean(dl.g.sub(e2))*mean(dl.g.sub(e3))
                                    om.1223 <- -mean(dl.g.sub(e1))*mean(dl.g.sub(e3))
                                    om.1312 <- om.1213
                                    om.1323 <- mean(dl.g.sub(e1))*mean(dl.g.sub(e2))
                                    om.2312 <- -mean(dl.g.sub(e3))*mean(dl.g.sub(e1))
                                    om.2313 <- mean(dl.g.sub(e2))*mean(dl.g.sub(e1))
                                  }
                                  ## Matrix A1
                                  A1 <- matrix(data = c(a12,-a21,0,0,0,
                                                        a13,0,0,0,-a31,
                                                        0,0,0,a23,-a32),ncol = 9,byrow = T)
                                  
                                  ## Matrix A2
                                  A2 <- matrix(data = c(Orth.max[[i]][,2],Orth.max[[i]][,1],0,0,0,
                                                        Orth.max[[i]][,3],0,0,0,Orth.max[[i]][,1],
                                                        0,0,0,Orth.max[[i]][,3],Orth.max[[i]][,2]), ncol = 9, byrow = T)
                                  
                                  ## MAtrix A3
                                  A3 <- matrix(data = c(Orth.max[[i]][,1],0,0,0,0,0,0,
                                                        0,0,0,Orth.max[[i]][,2],0,0,0,
                                                        0,0,0,0,0,0,Orth.max[[i]][,3]),ncol = 9,byrow = T)
                                  
                                  A <- matrix(data = rbind(A1,A2,A3),ncol = 9)
                                  
                                  ## Matrix OMEGA blocks
                                  om.hat <- matrix(data = c(om.1212, om.1213, om.1223,
                                                            om.1312, om.1313, om.1323,
                                                            om.2312, om.2313, om.2323),ncol = 3,byrow = T)
                                  om.hat.1 <- solve(om.hat)
                                  om2      <- matrix(data=0,ncol = 6,nrow = 3)
                                  om.first.row <- cbind(om.hat.1,om2)
                                  OMEGA1   <- rbind(om.first.row,matrix(data=0,nrow = 6,ncol = 9))
                                  OMEGA    <- matrix(data = c(om.1212, om.1213, om.1223,0,0,0,0,0,0, 
                                                              om.1312, om.1313, om.1323,0,0,0,0,0,0,
                                                              om.2312, om.2313, om.2323,0,0,0,0,0,0,
                                                              0,0,0,0,0,0,0,0,0,
                                                              0,0,0,0,0,0,0,0,0,
                                                              0,0,0,0,0,0,0,0,0,
                                                              0,0,0,0,0,0,0,0,0,
                                                              0,0,0,0,0,0,0,0,0,
                                                              0,0,0,0,0,0,0,0,0),ncol = 9,byrow = T)
                                  var.matrix[[i]] <- solve(A) %*% OMEGA %*% t(solve(A))
                                  
                                  ## t-statistic
                                  ## This should be the the chisq
                                  #t.stat[i] <-n*t((vec(Bmax[[i]])-vec(Atrue))) %*% t(A) %*% OMEGA1 %*% A %*% (vec(Bmax[[i]])-vec(Atrue))
                                  gamma.hat <- as.vector(Orth.max[[i]])
                                  gamma.0   <- as.vector(Orth.true[[1]])
                                  t.stat[i] <- n*((t(gamma.hat - gamma.0) %*% t(A) %*% OMEGA1 %*% A %*% (gamma.hat - gamma.0)))
                                  
                                  
                                  mc_results   <- list(statistics      = t.stat[i],
                                                       Choleski        = C[[i]],
                                                       U               = U,
                                                       Ort.est         = Orth.max[[i]],
                                                       Ort.true        = Orth.true[[i]],
                                                       B.est           = A.ID.pml[[i]],
                                                       variance        = var.matrix[[i]],
                                                       perfm_plm       = FR_pml_nt[i])
                                }
  stopCluster(cl)
  
  ciao[[mc]] <- first_mc_res
  ## Store results for each experiment
  pippo[[z]]           <- list(mc_results = ciao)
}

names(inference_results) <- paste0("p_",round(SEQ,2))

if (!dir.exists(here("results"))) {
  dir.create(here("results"))
}

save(inference_results,file = here("results",paste0("3D_specific_pml_n",n,".RData")))
rm(pippo)
