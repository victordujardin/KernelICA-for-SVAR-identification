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
  lhs            <- 2*pi*improvedLHS(n = n_lhs, k = 2)
}


for(z in 1:length(SEQ)){
  print(paste0("experiment ", z, " - first MC ", Sys.time()))
  print(paste0("sample size = ",n," initialization =  ", n_lhs, " MC = ",mc))
  i <- 1
  set.seed(11)
  e1     <- stdEpsrnd(type = 'general',dof = 5, n = n, p = SEQ[z])
  set.seed(22)
  e2     <- stdEpsrnd(type = 'general',dof = 5, n = n, p = SEQ[z])
  E        <- cbind(e1,e2)
  
  U                <- Btrue2Dz %*% t(E)
  sigg             <- cov(t(U))
  C[[i]]           <- t(chol(sigg))
  u_chol[[i]]      <- t(solve(C[[i]]) %*% (U))
  Orth.true[[i]]   <- solve(C[[i]]) %*% Btrue2Dz
  
  ## Run methods and record best initialization matrix
  for (kk in 1:nrow(lhs)){
    orth_matrix_init  <- matrix(data = c(cos(lhs[kk]), -sin(lhs[kk]), 
                                         sin(lhs[kk]),  cos(lhs[kk])), ncol = 2,byrow = T)
    p.start[[kk]]           <- as.vector(C[[i]] %*% orth_matrix_init)
    #p.start[[kk]]           <- as.vector(Btrue2Dz)
    ## Constrained optimization with more parameters than equations
    
    if (kurtosis(e1)>3 & kurtosis(e2)>3){
      ## super-gaussian case
      please           <- auglag(par = p.start[[kk]], fn = pseudo.log.lik2D.super, heq = orth.constr2D,
                                 control.outer = list(trace=F))$par
    }else{
      ## sub-gaussian case
      please           <- auglag(par = p.start[[kk]], fn = pseudo.log.lik2D.sub, heq = orth.constr2D,
                                 control.outer = list(trace=F))$par
    }
    
    Orth.hat           <- matrix(data = please,nrow = 2,byrow = F)
    B.init[[kk]]       <- C[[i]] %*% Orth.hat
    ilm.index[kk]      <- MD(W.hat = solve(B.init[[kk]]),A = Btrue2Dz)
    
  }
  
  ## Choose the best initialization
  best_init        <- p.start[[which.min(ilm.index)]]
  FR_pml_nt[i]     <- ilm.index[[which.min(ilm.index)]]
  A.pml[[i]]       <- B.init[[which.min(ilm.index)]]
  
  A.ID.pml[[i]]        <- maxfinder(A = A.pml[[i]])$A.id
  Orth.max[[i]]        <- solve(C[[i]]) %*% A.ID.pml[[i]]
  
  
  ## Try to built statistic based on MLE
  if (kurtosis(e1)>3 & kurtosis(e2)>3){
    ## Matrix A1 SUPER GUASSIAN)
    a12 <- (mean((-1)*(-1)) + mean((e2)*(-1*((e2)*(1+5)/(5-2+e2^2))))) * t(Orth.max[[i]][,2])
    a21 <- (mean((1+5)*((5-2-e2^2)/(5-2+e2^2)^2)) + mean((-e1)*e1)) *t(Orth.max[[i]][,1])
    
    ## Matrix OMEGA SUPER GAUSSIAN
    om.1212 <- mean((-e1)^2) + mean(((-e2)*(1+5)/(5-2+e2^2))^2) - 
      (2*mean((-e1)*e1)*mean(e2*((-e2)*(1+5)/(5-2+e2^2))))
    
    
  }else{
    ## Matrix A1 for SUB GAUSSIAN
    dl.g.sub <- function(e){
      2*pi*e + (pi/2)*(tanh(e*(pi/2)))
    }
    d2l.g.sub <- function(e){
      2*pi + ((pi/2)^2)*(1-tanh(e*(pi/2))^2)
    }
    a12 <- (mean((-1)*d2l.g.sub(e1)) + mean((e2)*dl.g.sub(e2))) * t(Orth.max[[i]][,2])
    a21 <- (mean((-1)*d2l.g.sub(e2)) + mean((e1)*dl.g.sub(e1))) * t(Orth.max[[i]][,1])
    
    ## Matrix OMEGA SUB GAUSSIAN
    om.1212 <- mean(dl.g.sub(e1)^2) + mean(dl.g.sub(e2)^2) - (2*mean((e1)*dl.g.sub(e1))*mean((e2)*dl.g.sub(e2)))
    
  }
  A1 <- matrix(data = c(a12,-a21),ncol = 4,byrow = T)
  
  ## Matrix A2
  A2 <- matrix(data = c(Orth.max[[i]][,2],Orth.max[[i]][,1]), ncol = 4, byrow = T)
  
  ## MAtrix A3
  A3 <- matrix(data = c(Orth.max[[i]][,1],0,0,
                        0,0,Orth.max[[i]][,2]),ncol = 4,byrow = T)
  
  A <- matrix(data = rbind(A1,A2,A3),ncol = 4)
  
  om.hat <- om.1212
  om.hat.1 <- om.hat^(-1)
  om2      <- matrix(data=0,ncol = 3,nrow = 1)
  om.first.row <- cbind(om.hat.1,om2)
  OMEGA1 <- rbind(om.first.row,matrix(data=0,nrow = 3,ncol = 4))
  OMEGA <- matrix(data = c(om.1212,0,0,0, 
                           0,0,0,0,
                           0,0,0,0,
                           0,0,0,0) ,ncol = 4,byrow = T)
  var.matrix[[i]] <- solve(A) %*% OMEGA %*% t(solve(A))
  
  
  ## t-statistic
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
                                  E        <- cbind(e1,e2)
                                  
                                  U                <- Btrue2Dz %*% t(E)
                                  sigg             <- cov(t(U))
                                  C[[i]]           <- t(chol(sigg))
                                  u_chol[[i]]      <- t(solve(C[[i]]) %*% (U))
                                  Orth.true[[i]]   <- solve(C[[i]]) %*% Btrue2Dz
                                  
                                  
                                  if (kurtosis(e1)>3 & kurtosis(e2)>3){
                                    ## super-gaussian case
                                    please           <- auglag(par = best_init, fn = pseudo.log.lik2D.super, heq = orth.constr2D,
                                                               control.outer = list(trace=F))$par
                                  }else{
                                    ## sub-gaussian case
                                    please           <- auglag(par = best_init, fn = pseudo.log.lik2D.sub, heq = orth.constr2D,
                                                               control.outer = list(trace=F))$par
                                  }
                                  
                                  Orth.hat           <- matrix(data = please,nrow = 2,byrow = F)
                                  A.pml[[i]]         <- C[[i]] %*% Orth.hat
                                  FR_pml_nt[i]       <- MD(W.hat = solve(A.pml[[i]]),A = Btrue2Dz)
                                  
                                  
                                  A.ID.pml[[i]]        <- maxfinder(A = A.pml[[i]])$A.id
                                  Orth.max[[i]]        <- solve(C[[i]]) %*% A.ID.pml[[i]]
                                  
                                  ## Try to built statistic based on MLE
                                  if (kurtosis(e1)>3 & kurtosis(e2)>3){
                                    ## Matrix A1 SUPER GUASSIAN)
                                    a12 <- (mean((-1)*(-1)) + mean((e2)*(-1*((e2)*(1+5)/(5-2+e2^2))))) * t(Orth.max[[i]][,2])
                                    a21 <- (mean((1+5)*((5-2-e2^2)/(5-2+e2^2)^2)) + mean((-e1)*e1)) *t(Orth.max[[i]][,1])
                                    
                                    ## Matrix OMEGA SUPER GAUSSIAN
                                    om.1212 <- mean((-e1)^2) + mean(((-e2)*(1+5)/(5-2+e2^2))^2) - 
                                      (2*mean((-e1)*e1)*mean(e2*((-e2)*(1+5)/(5-2+e2^2))))
                                    
                                    
                                  }else{
                                    ## Matrix A1 for SUB GAUSSIAN
                                    dl.g.sub <- function(e){
                                      2*pi*e + (pi/2)*(tanh(e*(pi/2)))
                                    }
                                    d2l.g.sub <- function(e){
                                      2*pi + ((pi/2)^2)*(1-tanh(e*(pi/2))^2)
                                    }
                                    a12 <- (mean((-1)*d2l.g.sub(e1)) + mean((e2)*dl.g.sub(e2))) * t(Orth.max[[i]][,2])
                                    a21 <- (mean((-1)*d2l.g.sub(e2)) + mean((e1)*dl.g.sub(e1))) * t(Orth.max[[i]][,1])
                                    
                                    ## Matrix OMEGA SUB GAUSSIAN
                                    om.1212 <- mean(dl.g.sub(e1)^2) + mean(dl.g.sub(e2)^2) - (2*mean((e1)*dl.g.sub(e1))*mean((e2)*dl.g.sub(e2)))
                                    
                                  }
                                  A1 <- matrix(data = c(a12,-a21),ncol = 4,byrow = T)
                                  
                                  ## Matrix A2
                                  A2 <- matrix(data = c(Orth.max[[i]][,2],Orth.max[[i]][,1]), ncol = 4, byrow = T)
                                  
                                  ## MAtrix A3
                                  A3 <- matrix(data = c(Orth.max[[i]][,1],0,0,
                                                        0,0,Orth.max[[i]][,2]),ncol = 4,byrow = T)
                                  
                                  A <- matrix(data = rbind(A1,A2,A3),ncol = 4)
                                  
                                  om.hat <- om.1212
                                  om.hat.1 <- om.hat^(-1)
                                  om2      <- matrix(data=0,ncol = 3,nrow = 1)
                                  om.first.row <- cbind(om.hat.1,om2)
                                  OMEGA1 <- rbind(om.first.row,matrix(data=0,nrow = 3,ncol = 4))
                                  OMEGA <- matrix(data = c(om.1212,0,0,0, 
                                                           0,0,0,0,
                                                           0,0,0,0,
                                                           0,0,0,0) ,ncol = 4,byrow = T)
                                  var.matrix[[i]] <- solve(A) %*% OMEGA %*% t(solve(A))
                                  
                                  
                                  ## t-statistic
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

inference_results <- pippo

names(inference_results) <- paste0("p_",round(SEQ,2))

if (!dir.exists(here("results"))) {
  dir.create(here("results"))
}

save(inference_results,file = here("results",paste0("2D_specific_pml_n",n,".RData")))
