#source("fuctions-needed/testlik.R")
set.seed(55) ## To make example reproducible
lhs    <- 2*pi*improvedLHS(n = n_lhs, k = 3)

## Sample size
n <- nrow(resid(var.interactions))

Wscaled.dc <- vector("list", n_lhs)
Wscaled.ica <- vector("list", n_lhs)
w0          <- vector("list", n_lhs)
A.ica       <- vector("list", n_lhs)
A.dc        <- vector("list", n_lhs)
A.ID.pml    <- vector("list", n_lhs)
A.ID.dc     <- vector("list", n_lhs)
A.ID.ica    <- vector("list", n_lhs)
A.ID.pml    <- vector("list", n_lhs)

ts_dim <- (1 + sqrt(1 + 8 * ncol(lhs))) / 2
combns <- combn(ts_dim, 2, simplify = FALSE)

for (kk in 1:nrow(lhs)){
  if (kk %% 100 == 0) {
    print(paste0("Initialization # ",kk))
  }
  #### Initialization ####
  we <- diag(ts_dim)
  pv <- lhs[kk,]
  for (i in seq_along(pv)) {
    
    tmp <- diag(ts_dim)
    tmp[combns[[i]][1], combns[[i]][1]] <- cos(pv[i])
    tmp[combns[[i]][2], combns[[i]][2]] <- cos(pv[i])
    tmp[combns[[i]][1], combns[[i]][2]] <- - sin(pv[i])
    tmp[combns[[i]][2], combns[[i]][1]] <- sin(pv[i])
    we <- we %*% tmp
  }
  w0[[kk]] <- we
  
  ## PML
  p.start           <- as.vector(w0[[kk]])
  
  please           <- auglag(par = p.start, fn = pseudo.log.lik3D.super, heq = orth.constr3D,
                             control.outer = list(trace=F))$par 
  Apml.raw         <- matrix(data = please,nrow = 3,byrow = F)
  Bpml             <- (C) %*% (Apml.raw)
  A.ID.pml[[kk]]   <- maxfinder(A = Bpml)$A.id
  Apml             <- solve(C) %*% A.ID.pml[[kk]]
  
  ## FastICA
  icares             <- fastICA(tu_chol, n.comp = ncol(sigg),tol=1e-14, maxit=3000, verbose=FALSE, w.init = w0[[kk]])
  W                  <- t((icares$K) %*% (icares$W)) 
  Wscal.ica          <- rescaleVar(W_hat = W,ut = u_chol)$Ws 
  Wscaled.ica[[kk]]  <- Wscal.ica %*% solve(C) 
  A.ica[[kk]]        <- solve(Wscaled.ica[[kk]]) # A is the mixing matrix
  
  ## Distance Covariance ##
  DC                 <- steadyICA(X = tu_chol,n.comp = ncol(sigg),w.init = w0[[kk]])
  Wscaled.dc[[kk]]   <- solve(DC$W) %*% solve(C)
  A.dc[[kk]]         <- solve(Wscaled.dc[[kk]])
  
  #Identification scheme
  A.ID.dc[[kk]]       <- maxfinder(A = A.dc[[kk]])$A.id
  A.ID.ica[[kk]]      <- maxfinder(A = A.ica[[kk]])$A.id
}
dd <- NULL
##CVM seems not to require any initialization
## the file dd contains the copula under the null hypothesis of statistical independence
if (is.null(dd)) {
  dd <- indepTestSim(n, ncol(tu_chol), verbose = T)
}
#save(dd, file = 'indepNULL3D_n400.RData')
lower <- rep(0, ncol(tu_chol) * (ncol(tu_chol) - 1)/2)
upper <- rep(pi, ncol(tu_chol) * (ncol(tu_chol) - 1)/2)
de_control <- list(itermax = 500, steptol = 100, 
                   trace = FALSE)
## First step of optimization with DEoptim
de_res <- DEoptim(testlik, lower = lower, upper = upper, 
                  control = de_control, faklow = C, u = u.ica, dd = dd)


k <- ncol(tu_chol)
iter2  <- 150
## Second step of optimization. Creating randomized starting angles around the optimized angles
## here maybe you insert the LHS stuff
theta_rot <- matrix(rnorm(n = (k * (k - 1)/2) * iter2, mean = de_res$optim$bestmem, 
                          sd = 0.3), (k * (k - 1)/2), iter2)
theta_rot <- cbind(theta_rot, de_res$optim$bestmem)
# Start vectors for iterative optimization approach
startvec_list <- as.list(as.data.frame(theta_rot))
erg_list <- lapply(X = startvec_list, FUN = optim, fn = testlik, 
                   gr = NULL, faklow = C, u = u.ica, dd = dd, 
                   method = ifelse(k == 2, "Brent", "Nelder-Mead"), 
                   lower = ifelse(k == 2, -.Machine$double.xmax/2, -Inf), 
                   upper = ifelse(k == 2, .Machine$double.xmax/2, Inf), 
                   control = list(maxit = 1000), 
                   hessian = FALSE)
# Print log-likelihood values from local maxima
logliks <- sapply(erg_list, "[[", "value")
if (min(logliks) < de_res$optim$bestval) {
  params <- sapply(erg_list, "[[", "par", simplify = FALSE)
  par_o <- params[[which.min(logliks)]]
  logs <- min(logliks)
  inc <- 1
} else {
  par_o <- de_res$optim$bestmem
  logs <- de_res$optim$bestval
  inc <- 0
}
Acvm <- rotmat(par_o, C)   ## Estimated mixing matrix

c.ica <- vector("list", n_lhs)
c.dc  <- vector("list", n_lhs)
c.pml <- vector("list",n_lhs)