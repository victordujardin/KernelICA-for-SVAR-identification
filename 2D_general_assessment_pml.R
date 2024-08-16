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

## Open ica-svars-comparative-analysis.Rproj
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

## Number of shape parameter values (in the specific assessment and for inference table we only focuse on 4 distributional scenario - See paper)
seq1   <- seq(from = 0.5, to = 3.5, length.out = 15) ## Sequence of p-Generalized parameter for general evaluation exercise (1st part)
seq2   <- seq(from = 4, to = 100, length.out = 5) ## Sequence of p-Generalized parameter for general evaluation exercise (2nd part)
#SEQ     <- c(0.5,1.57,2.43,100) ## Sequence of p-Generalized parameter for statistical inference exercise
SEQ    <- as.numeric(c(seq1, seq2))

## File with variables definitions and parameter settings

n              <- 100 # sample size
choose_cluster <- 3

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
  print(paste0("sample size = ",n," initialization =  ", n_lhs, " MC = ",mc, ", size = ",ncol(lhs)))
  nr_cluster <- choose_cluster
  cl <- makeCluster(nr_cluster)
  registerDoParallel(cl)
  ## Results reproducible 
  set.seed(123+z)
  ciao <- foreach(i = 1:mc,
                  .packages = c('pgnorm', 'normtest', 'steadyICA', "ProDenICA",
                                'fastICA', 'JADE', 'tidyverse', 'magrittr', 'gtools', 
                                'DEoptim', "matrixcalc", "copula","moments",
                                "alabama")) %dorng% {
                                  
                                  Sys.sleep(1)  
                                  
                                  e1     <- stdEpsrnd(type = 'general',dof = 5, n = n, p = SEQ[z])
                                  e2     <- stdEpsrnd(type = 'general',dof = 5, n = n, p = SEQ[z])
                                  E        <- cbind(e1,e2)
                                  
                                  Atrue            <- mixmat(p = ncol(lhs))
                                  U                <- Atrue %*% t(E)
                                  sigg             <- cov(t(U))
                                  C[[i]]           <- t(chol(sigg))
                                  u_chol[[i]]      <- t(solve(C[[i]]) %*% (U))
                                  Orth.true[[i]]   <- solve(C[[i]]) %*% Atrue
                                  
                                  ## Run methods and record best initialization matrix
                                  for (kk in 1:nrow(lhs)){
                                    orth_matrix_init  <- matrix(data = c(cos(lhs[kk]), -sin(lhs[kk]), 
                                                                         sin(lhs[kk]),  cos(lhs[kk])), ncol = 2,byrow = T)
                                    
                                    p.start[[kk]]           <- as.vector(C[[i]] %*% orth_matrix_init)
                                    
                                    ## Constrained optimization with more parameters than equations
                                    
                                    if (kurtosis(e1)-3>0 & kurtosis(e2)-3>0){
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
                                    ilm.index[kk]      <- MD(W.hat = solve(B.init[[kk]]),A = Atrue)
                                    
                                  }
                                  
                                  ## Choose the best initialization
                                  best_init        <- p.start[[which.min(ilm.index)]]
                                  FR_pml_nt[i]     <- ilm.index[[which.min(ilm.index)]]
                                  A.pml[[i]]       <- B.init[[which.min(ilm.index)]]
                                  
                                  if (kurtosis(e1)-3>0 & kurtosis(e2)-3>0){
                                    ## super-gaussian case
                                    please           <- auglag(par = best_init, fn = pseudo.log.lik2D.super, heq = orth.constr2D,
                                                               control.outer = list(trace=F))$par
                                  }else{
                                    ## sub-gaussian case
                                    please           <- auglag(par = best_init, fn = pseudo.log.lik2D.sub, heq = orth.constr2D,
                                                               control.outer = list(trace=F))$par
                                  }
                                  
                                  Orth.true[[i]]     <- matrix(data = please,nrow = 2,byrow = F)
                                  A.pml[[i]]         <- C[[i]] %*% Orth.true[[i]]
                                  FR_pml_nt[i]       <- MD(W.hat = solve(A.pml[[i]]),A = Atrue)
                                  
                                  
                                  A.ID.pml[[i]]        <- maxfinder(A = A.pml[[i]])$A.id
                                  
                                  
                                  mc_results   <- list(Choleski        = C[[i]],
                                                       TrueMixing      = Atrue,
                                                       U               = U,
                                                       pml_init_mixing = B.init,
                                                       Ort.true        = Orth.true[[i]],
                                                       B.est           = A.ID.pml[[i]],
                                                       perfm_plm       = FR_pml_nt[i])
                                  
                                }
  stopCluster(cl)
  ## Store results for each experiment
  pippo[[z]]           <- list(mc_results = ciao)
}


inference_results <- pippo

names(inference_results) <- paste0("p_",round(SEQ,2))

if (!dir.exists(here("results"))) {
  dir.create(here("results"))
}
save(inference_results,file = here("results", paste0("2D_pml_general_n",n,".RData")))
