## All auxiliary functions you need

# works fine up to 7 (maybe 8) variables, for more the function "permutations" gets slow (factorial grows quickly!!)
# need package "gtools" for function "permuations"

best_permutation <- function(W,W_boot) {
  
  library(gtools)
  nr <- nrow(W)
  
  allperm <- gtools::permutations(n=nr,r=nr)
  
  # to switch signs in ICs: first all possibilities with 1 IC switched, than with two, ...
  ls <- list()
  s <- array(0,dim=c(nr,1))
  for (i in 1:nr) {
    ls[[i]] <- combinations(n=nr,r=i) # in how many ways I can switch i signs
    s[i] <- dim(ls[[i]])[1]
  }
  
  ss <- sum(s)
  
  # for each permutation: original signs + ss different combinations of switched signs
  
  # get diagonal matrices with 1 and -1 in the columns to switch sign by multiplication
  Sgn <- list()
  Sgn[[1]] <- diag(1,nr)
  cnt <- 1
  for (i in 1:nr) {
    for (j in 1:s[i]) {
      cnt <- cnt+1
      Sgn[[cnt]] <- diag(1,nr)
      # here put -1 in right places:
      Sgn[[cnt]][ls[[i]][j,],] <- Sgn[[cnt]][ls[[i]][j,],]*(-1)
    }
  }
  
  
  norm <- 1e5 #just some big number
  #norm <- array(0,dim=c(nrow(allperm),ss+1))
  
  for (i in 1:nrow(allperm)) {
    
    # get one permutation
    Perm <- array(0,dim=c(nr,nr))
    for (j in 1:nr) Perm[allperm[i,j],j] <- 1
    Temp <- Perm%*%W_boot
    
    # try all different sign-combinations for the permutation
    for (j in 1:(ss+1)) {
      if (elp_norm(W-(Sgn[[j]]%*%Temp),2) < norm ) {
        norm <- elp_norm(W-(Sgn[[j]]%*%Temp),2)
        Perm_save <- Perm
        Sgn_save <- Sgn[[j]]
        #print(i)
        #print(j)
      }
    }
    
    #norm[i,] <- Linfnorm(W-Perm%*%W_boot)
    #norm[i,] <- elp_norm(W-(Perm%*%W_boot),2)
    
  }
  
  Sgn_save%*%Perm_save%*%W_boot
  #Perm_save%*%W_boot - like that it prunes all away!!
}

Linfnorm <- function(M) {
  
  # infinity norm of a squared Matrix M (max. absolute row sum)
  norm <- array(0,dim=nrow(M))
  for (i in 1:nrow(M)) norm[i] <- sum(abs(M[i,]))
  max(norm)
}

elp_norm <- function(M,p) {
  
  # elementwise p-norm
  (sum(abs(M)^p))^(1/p)
}

entries.matrix <- function(A.ID, n_lhs, d){
  
  B.entry <- array(data = 0,dim = c(n_lhs, d,d))
  #B.factors <- vector("list",length = d*d)
  l <- 0
  B.factors = list()
  names     = list()
  for (i in 1:d){
    for (j in 1:d){
      l <- l+1
      B.entry[,i,j]    <- do.call(rbind, lapply(A.ID,'[',i,j))
      B.factors[[l]]   <- B.entry[,i,j]
      names[[l]]       <- paste0("b",i,j)
    }
  }
  big_data <- as.data.frame(do.call(cbind, B.factors))
  colnames(big_data) <- names
  ddf      <- as.data.frame(gather(data = big_data))
  return(ddf)
}

getcoefSVARS <- function(x, B_hat) 
{
  if (inherits(x, "var.boot")) {
    u <- x$residuals
    Tob <- nrow(u)
    k <- ncol(u)
    residY <- u
  } else {
    u <- residuals(x)
    Tob <- nrow(u)
    k <- ncol(u)
    residY <- u
  }
  if (inherits(x, "var.boot")) {
    p <- x$p
    y <- x$y
    yOut = x$y
    type = x$type
    coef_x = x$coef_x
  } else if (inherits(x, "varest")) {
    p <- x$p
    y <- t(x$y)
    yOut = x$y
    type = x$type
    coef_x = coef(x)
  } else if (inherits(x, "nlVar")) {
    p <- x$lag
    y <- t(x$model[, 1:k])
    yOut <- x$model[, 1:k]
    coef_x <- t(coef(x))
    if (inherits(x, "VECM")) {
      coef_x <- t(VARrep(x))
    }
    if (rownames(coef_x)[1] %in% c("Intercept", "constant")) {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
    } else if (rownames(coef_x)[1] == "Trend") {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
    }
    if (rownames(coef_x)[1] %in% c("Intercept", "constant", 
                                   "Trend")) {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
    }
    type <- x$include
    coef_x <- split(coef_x, rep(1:ncol(coef_x), each = nrow(coef_x)))
    coef_x <- lapply(coef_x, as.matrix)
  } else if (inherits(x, "list")) {
    p <- x$order
    y <- t(x$data)
    yOut <- x$data
    coef_x <- t(coef(x))
    coef_x <- x$coef
    if (x$cnst == TRUE) {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
      type = "const"
    }
    coef_x <- split(coef_x, rep(1:ncol(coef_x), each = nrow(coef_x)))
    coef_x <- lapply(coef_x, as.matrix)
  } else if (inherits(x, "vec2var")) {
    coef_x <- vector("list", length = k)
    names(coef_x) <- colnames(x$y)
    p <- x$p
    y <- t(x$y)
    yOut <- x$y
    for (i in seq_len(k)) {
      for (j in seq_len(p)) coef_x[[i]] <- c(coef_x[[i]], 
                                             x$A[[j]][i, ])
      coef_x[[i]] <- c(coef_x[[i]], x$deterministic[i, 
      ])
    }
    coef_x <- lapply(coef_x, matrix)
    type <- "const"
  } else {
    stop("Object class is not supported")
  }
  
  ## HERE THERE WAS THE METHOD
  
  if (inherits(x, "var.boot")) {
    A_hat <- coef_x
  } else {
    A <- matrix(0, nrow = k, ncol = k * p)
    for (i in 1:k) {
      A[i, ] <- coef_x[[i]][1:(k * p), 1]
    }
    A_hat <- A
    if (type == "const") {
      v <- rep(1, k)
      for (i in 1:k) {
        v[i] <- coef_x[[i]][(k * p + 1), 1]
      }
      A_hat <- cbind(v, A)
    } else if (type == "trend") {
      trend <- rep(1, k)
      for (i in 1:k) {
        trend[i] <- coef_x[[i]][(k * p + 1), 1]
      }
      A_hat <- cbind(trend, A)
    } else if (type == "both") {
      v <- rep(1, k)
      for (i in 1:k) {
        v[i] <- coef_x[[i]][(k * p + 1), 1]
      }
      trend <- rep(1, k)
      for (i in 1:k) {
        trend[i] <- coef_x[[i]][(k * p + 2), 1]
      }
      A_hat <- cbind(v, trend, A)
    }
  }
  result <- list(B = B_hat, A_hat = A_hat, 
                 n = Tob, type = type, y = yOut, p = unname(p), K = k)
  class(result) <- "svars"
  return(result)
}


id.cvm.gp <- function (x, dd = NULL, itermax = 500, steptol = 100, iter2 = 75,ordering) 
{
  if (inherits(x, "var.boot")) {
    u <- x$residuals
    Tob <- nrow(u)
    k <- ncol(u)
    residY <- u
  }
  else {
    u <- residuals(x)
    Tob <- nrow(u)
    k <- ncol(u)
    residY <- u
  }
  if (inherits(x, "var.boot")) {
    p <- x$p
    y <- t(x$y)
    yOut <- x$y
    type = x$type
    coef_x = x$coef_x
  }
  else if (inherits(x, "varest")) {
    p <- x$p
    y <- t(x$y)
    yOut <- x$y
    type = x$type
    coef_x = coef(x)
  }
  else if (inherits(x, "nlVar")) {
    p <- x$lag
    y <- t(x$model[, 1:k])
    yOut <- x$model[, 1:k]
    coef_x <- t(coef(x))
    if (inherits(x, "VECM")) {
      coef_x <- t(VARrep(x))
    }
    if (rownames(coef_x)[1] %in% c("Intercept", "constant")) {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
    }
    else if (rownames(coef_x)[1] == "Trend") {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
    }
    if (rownames(coef_x)[1] %in% c("Intercept", "constant", 
                                   "Trend")) {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
    }
    type <- x$include
    coef_x <- split(coef_x, rep(1:ncol(coef_x), each = nrow(coef_x)))
    coef_x <- lapply(coef_x, as.matrix)
  }
  else if (inherits(x, "list")) {
    p <- x$order
    y <- t(x$data)
    yOut <- x$data
    coef_x <- x$coef
    if (x$cnst == TRUE) {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
      type = "const"
    }
    coef_x <- split(coef_x, rep(1:ncol(coef_x), each = nrow(coef_x)))
    coef_x <- lapply(coef_x, as.matrix)
  }
  else if (inherits(x, "vec2var")) {
    coef_x <- vector("list", length = k)
    names(coef_x) <- colnames(x$y)
    p <- x$p
    y <- t(x$y)
    yOut <- x$y
    for (i in seq_len(k)) {
      for (j in seq_len(p)) coef_x[[i]] <- c(coef_x[[i]], 
                                             x$A[[j]][i, ])
      coef_x[[i]] <- c(coef_x[[i]], x$deterministic[i, 
      ])
    }
    coef_x <- lapply(coef_x, matrix)
    type <- "const"
  }
  else {
    stop("Object class is not supported")
  }
  sigg1 <- crossprod(u)/(Tob - 1 - k * p)
  faklow1 <- t(chol(sigg1))
  if (is.null(dd)) {
    dd <- indepTestSim(Tob, k, N = 100, verbose = FALSE)
  }
  lower <- rep(0, k * (k - 1)/2)
  upper <- rep(pi, k * (k - 1)/2)
  if (is.null(dd)) {
    dd <- indepTestSim(Tob, k, verbose = F)
  }
  de_control <- list(itermax = itermax, steptol = steptol, 
                     trace = FALSE)
  de_res <- DEoptim(testlik, lower = lower, upper = upper, 
                    control = de_control, faklow = faklow1, u = u, dd = dd)
  theta_rot <- matrix(rnorm(n = (k * (k - 1)/2) * iter2, mean = de_res$optim$bestmem, 
                            sd = 0.3), (k * (k - 1)/2), iter2)
  theta_rot <- cbind(theta_rot, de_res$optim$bestmem)
  startvec_list <- as.list(as.data.frame(theta_rot))
  erg_list <- lapply(X = startvec_list, FUN = optim, fn = testlik, 
                     gr = NULL, faklow = faklow1, u = u, dd = dd, method = ifelse(k == 
                                                                                    2, "Brent", "Nelder-Mead"), lower = ifelse(k == 2, 
                                                                                                                               -.Machine$double.xmax/2, -Inf), upper = ifelse(k == 
                                                                                                                                                                                2, .Machine$double.xmax/2, Inf), control = list(maxit = 1000), 
                     hessian = FALSE)
  logliks <- sapply(erg_list, "[[", "value")
  if (min(logliks) < de_res$optim$bestval) {
    params <- sapply(erg_list, "[[", "par", simplify = FALSE)
    par_o <- params[[which.min(logliks)]]
    logs <- min(logliks)
    inc <- 1
  }
  else {
    par_o <- de_res$optim$bestmem
    logs <- de_res$optim$bestval
    inc <- 0
  }
  B_hat  <- rotmat(par_o, faklow1)
  
  if (ordering == "frob") {
    B_hat      <- B_hat %*% myfrob(A.hat = B_hat,A = x$B_original)
  }else if(ordering == "maxfinder"){
    B_hat      <- maxfinder(A = B_hat)$A.id
  }else{
    stop("Specifiy method for column-permutation indeterminacy")
  }
  
  
  if (inherits(x, "var.boot")) {
    A_hat <- coef_x
  }
  else {
    A <- matrix(0, nrow = k, ncol = k * p)
    for (i in 1:k) {
      A[i, ] <- coef_x[[i]][1:(k * p), 1]
    }
    A_hat <- A
    if (type == "const") {
      v <- rep(1, k)
      for (i in 1:k) {
        v[i] <- coef_x[[i]][(k * p + 1), 1]
      }
      A_hat <- cbind(v, A)
    }
    else if (type == "trend") {
      trend <- rep(1, k)
      for (i in 1:k) {
        trend[i] <- coef_x[[i]][(k * p + 1), 1]
      }
      A_hat <- cbind(trend, A)
    }
    else if (type == "both") {
      v <- rep(1, k)
      for (i in 1:k) {
        v[i] <- coef_x[[i]][(k * p + 1), 1]
      }
      trend <- rep(1, k)
      for (i in 1:k) {
        trend[i] <- coef_x[[i]][(k * p + 2), 1]
      }
      A_hat <- cbind(v, trend, A)
    }
  }
  result <- list(B = B_hat, A_hat = A_hat, method = "Cramer-von Mises distance", 
                 n = Tob, type = type, y = yOut, p = p, K = k, rotation_angles = par_o, 
                 inc = inc, test.stats = logs, iter1 = de_res$optim$iter, 
                 test1 = de_res$optim$bestval, test2 = min(logliks))
  class(result) <- "svars"
  return(result)
}

id.dc.gp <- function (x, PIT= T, ww, ordering) 
{
  if (inherits(x, "var.boot")) {
    u <- x$residuals
    Tob <- nrow(u)
    k <- ncol(u)
    residY <- u
  }
  else {
    u <- residuals(x)
    Tob <- nrow(u)
    k <- ncol(u)
    residY <- u
  }
  if (inherits(x, "var.boot")) {
    p <- x$p
    y <- t(x$y)
    yOut = x$y
    type = x$type
    coef_x = x$coef_x
  }
  else if (inherits(x, "varest")) {
    p <- x$p
    y <- t(x$y)
    yOut = x$y
    type = x$type
    coef_x = coef(x)
  }
  else if (inherits(x, "nlVar")) {
    p <- x$lag
    y <- t(x$model[, 1:k])
    yOut <- x$model[, 1:k]
    coef_x <- t(coef(x))
    if (inherits(x, "VECM")) {
      coef_x <- t(VARrep(x))
    }
    if (rownames(coef_x)[1] %in% c("Intercept", "constant")) {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
    }
    else if (rownames(coef_x)[1] == "Trend") {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
    }
    if (rownames(coef_x)[1] %in% c("Intercept", "constant", 
                                   "Trend")) {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
    }
    type <- x$include
    coef_x <- split(coef_x, rep(1:ncol(coef_x), each = nrow(coef_x)))
    coef_x <- lapply(coef_x, as.matrix)
  }
  else if (inherits(x, "list")) {
    p <- x$order
    y <- t(x$data)
    yOut <- x$data
    coef_x <- t(coef(x))
    coef_x <- x$coef
    if (x$cnst == TRUE) {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
      type = "const"
    }
    coef_x <- split(coef_x, rep(1:ncol(coef_x), each = nrow(coef_x)))
    coef_x <- lapply(coef_x, as.matrix)
  }
  else if (inherits(x, "vec2var")) {
    coef_x <- vector("list", length = k)
    names(coef_x) <- colnames(x$y)
    p <- x$p
    y <- t(x$y)
    yOut <- x$y
    for (i in seq_len(k)) {
      for (j in seq_len(p)) coef_x[[i]] <- c(coef_x[[i]], 
                                             x$A[[j]][i, ])
      coef_x[[i]] <- c(coef_x[[i]], x$deterministic[i, 
      ])
    }
    coef_x <- lapply(coef_x, matrix)
    type <- "const"
  }
  else {
    stop("Object class is not supported")
  }
  sigg <- crossprod(u)/(Tob - 1 - k * p)
  P_chol <- t(chol(sigg))
  u_chol <- t(solve(P_chol) %*% t(u))
  ICA    <- steadyICA(u_chol, symmetric = F, PIT = PIT, w.init = ww)
  Praw   <- P_chol %*% ICA$W
  if (ordering == "frob") {
    P      <- Praw %*% myfrob(A.hat = Praw,A = x$B_original)
  }else if(ordering == "maxfinder"){
    P      <- maxfinder(A = Praw)$A.id
  }else{
    stop("Specifiy method for column-permutation indeterminacy")
  }
  
  
  if (inherits(x, "var.boot")) {
    A_hat <- coef_x
  }
  else {
    A <- matrix(0, nrow = k, ncol = k * p)
    for (i in 1:k) {
      A[i, ] <- coef_x[[i]][1:(k * p), 1]
    }
    A_hat <- A
    if (type == "const") {
      v <- rep(1, k)
      for (i in 1:k) {
        v[i] <- coef_x[[i]][(k * p + 1), 1]
      }
      A_hat <- cbind(v, A)
    }
    else if (type == "trend") {
      trend <- rep(1, k)
      for (i in 1:k) {
        trend[i] <- coef_x[[i]][(k * p + 1), 1]
      }
      A_hat <- cbind(trend, A)
    }
    else if (type == "both") {
      v <- rep(1, k)
      for (i in 1:k) {
        v[i] <- coef_x[[i]][(k * p + 1), 1]
      }
      trend <- rep(1, k)
      for (i in 1:k) {
        trend[i] <- coef_x[[i]][(k * p + 2), 1]
      }
      A_hat <- cbind(v, trend, A)
    }
  }
  result <- list(B = P, A_hat = A_hat, method = "Distance covariances", 
                 n = Tob, type = type, y = yOut, p = unname(p), K = k, 
                 PIT = PIT)
  class(result) <- "svars"
  return(result)
}

id.fastICA.gp <- function (x, ww, ordering) 
{
  if (inherits(x, "var.boot")) {
    u <- x$residuals
    Tob <- nrow(u)
    k <- ncol(u)
    residY <- u
  } else {
    u <- residuals(x)
    Tob <- nrow(u)
    k <- ncol(u)
    residY <- u
  }
  if (inherits(x, "var.boot")) {
    p <- x$p
    y <- t(x$y)
    yOut = x$y
    type = x$type
    coef_x = x$coef_x
  } else if (inherits(x, "varest")) {
    p <- x$p
    y <- t(x$y)
    yOut = x$y
    type = x$type
    coef_x = coef(x)
  } else if (inherits(x, "nlVar")) {
    p <- x$lag
    y <- t(x$model[, 1:k])
    yOut <- x$model[, 1:k]
    coef_x <- t(coef(x))
    if (inherits(x, "VECM")) {
      coef_x <- t(VARrep(x))
    }
    if (rownames(coef_x)[1] %in% c("Intercept", "constant")) {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
    } else if (rownames(coef_x)[1] == "Trend") {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
    } 
    if (rownames(coef_x)[1] %in% c("Intercept", "constant", 
                                   "Trend")) {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
    }
    type <- x$include
    coef_x <- split(coef_x, rep(1:ncol(coef_x), each = nrow(coef_x)))
    coef_x <- lapply(coef_x, as.matrix)
  } else if (inherits(x, "list")) {
    p <- x$order
    y <- t(x$data)
    yOut <- x$data
    coef_x <- t(coef(x))
    coef_x <- x$coef
    if (x$cnst == TRUE) {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
      type = "const"
    }
    coef_x <- split(coef_x, rep(1:ncol(coef_x), each = nrow(coef_x)))
    coef_x <- lapply(coef_x, as.matrix)
  } else if (inherits(x, "vec2var")) {
    coef_x <- vector("list", length = k)
    names(coef_x) <- colnames(x$y)
    p <- x$p
    y <- t(x$y)
    yOut <- x$y
    for (i in seq_len(k)) {
      for (j in seq_len(p)) coef_x[[i]] <- c(coef_x[[i]], 
                                             x$A[[j]][i, ])
      coef_x[[i]] <- c(coef_x[[i]], x$deterministic[i, 
      ])
    }
    coef_x <- lapply(coef_x, matrix)
    type <- "const"
  } else {
    stop("Object class is not supported")
  }
  sigg <- crossprod(u)/(Tob - 1 - k * p)
  P_chol <- t(chol(sigg))
  u_chol <- t(solve(P_chol) %*% t(u))
  fast.ica           <- fastICA(u_chol, n.comp = k, tol = 1e-14, w.init = ww, 
                                maxit = 1000,verbose = FALSE)
  W.ica              <- t((fast.ica$K) %*% (fast.ica$W))
  W.scal.ica         <- rescaleVar(W_hat = W.ica, ut = t(u_chol))$Ws
  Wstar              <- W.scal.ica %*% solve(P_chol)
  Praw               <- solve(Wstar)
  if (ordering == "frob") {
    P      <- Praw %*% myfrob(A.hat = Praw,A = x$B_original)
  }else if(ordering == "maxfinder"){
    P      <- maxfinder(A = Praw)$A.id
  }else{
    stop("Specifiy method for column-permutation indeterminacy")
  }
  
  if (inherits(x, "var.boot")) {
    A_hat <- coef_x
  }
  else {
    A <- matrix(0, nrow = k, ncol = k * p)
    for (i in 1:k) {
      A[i, ] <- coef_x[[i]][1:(k * p), 1]
    }
    A_hat <- A
    if (type == "const") {
      v <- rep(1, k)
      for (i in 1:k) {
        v[i] <- coef_x[[i]][(k * p + 1), 1]
      }
      A_hat <- cbind(v, A)
    }
    else if (type == "trend") {
      trend <- rep(1, k)
      for (i in 1:k) {
        trend[i] <- coef_x[[i]][(k * p + 1), 1]
      }
      A_hat <- cbind(trend, A)
    }
    else if (type == "both") {
      v <- rep(1, k)
      for (i in 1:k) {
        v[i] <- coef_x[[i]][(k * p + 1), 1]
      }
      trend <- rep(1, k)
      for (i in 1:k) {
        trend[i] <- coef_x[[i]][(k * p + 2), 1]
      }
      A_hat <- cbind(v, trend, A)
    }
  }
  result <- list(B = P, A_hat = A_hat, method = "FastICA", 
                 n = Tob, type = type, y = yOut, p = unname(p), K = k)
  class(result) <- "svars"
  return(result)
}


id.KernelICA.gp <- function (x, ww, ordering) 
{
  if (inherits(x, "var.boot")) {
    u <- x$residuals
    Tob <- nrow(u)
    k <- ncol(u)
    residY <- u
  } else {
    u <- residuals(x)
    Tob <- nrow(u)
    k <- ncol(u)
    residY <- u
  }
  if (inherits(x, "var.boot")) {
    p <- x$p
    y <- t(x$y)
    yOut = x$y
    type = x$type
    coef_x = x$coef_x
  } else if (inherits(x, "varest")) {
    p <- x$p
    y <- t(x$y)
    yOut = x$y
    type = x$type
    coef_x = coef(x)
  } else if (inherits(x, "nlVar")) {
    p <- x$lag
    y <- t(x$model[, 1:k])
    yOut <- x$model[, 1:k]
    coef_x <- t(coef(x))
    if (inherits(x, "VECM")) {
      coef_x <- t(VARrep(x))
    }
    if (rownames(coef_x)[1] %in% c("Intercept", "constant")) {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
    } else if (rownames(coef_x)[1] == "Trend") {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
    } 
    if (rownames(coef_x)[1] %in% c("Intercept", "constant", 
                                   "Trend")) {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
    }
    type <- x$include
    coef_x <- split(coef_x, rep(1:ncol(coef_x), each = nrow(coef_x)))
    coef_x <- lapply(coef_x, as.matrix)
  } else if (inherits(x, "list")) {
    p <- x$order
    y <- t(x$data)
    yOut <- x$data
    coef_x <- t(coef(x))
    coef_x <- x$coef
    if (x$cnst == TRUE) {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
      type = "const"
    }
    coef_x <- split(coef_x, rep(1:ncol(coef_x), each = nrow(coef_x)))
    coef_x <- lapply(coef_x, as.matrix)
  } else if (inherits(x, "vec2var")) {
    coef_x <- vector("list", length = k)
    names(coef_x) <- colnames(x$y)
    p <- x$p
    y <- t(x$y)
    yOut <- x$y
    for (i in seq_len(k)) {
      for (j in seq_len(p)) coef_x[[i]] <- c(coef_x[[i]], 
                                             x$A[[j]][i, ])
      coef_x[[i]] <- c(coef_x[[i]], x$deterministic[i, 
      ])
    }
    coef_x <- lapply(coef_x, matrix)
    type <- "const"
  } else {
    stop("Object class is not supported")
  }
  sigg <- crossprod(u)/(Tob - 1 - k * p)
  P_chol <- t(chol(sigg))
  u_chol <- t(solve(P_chol) %*% t(u))
  kica_init <- KernelICA::kernel_ica(u_chol, variant = "kgv", kernel = "gauss", eps = 1e-14, init = list(ww))
  Wscaled_kica <- rescaleVar(W_hat = kica_init$W, ut = t(u_chol))$Ws
  Wstar <- Wscaled_kica %*% solve(P_chol)
  Praw <- solve(Wstar)
  if (ordering == "frob") {
    P      <- Praw %*% myfrob(A.hat = Praw,A = x$B_original)
  }else if(ordering == "maxfinder"){
    P      <- maxfinder(A = Praw)$A.id
  }else{
    stop("Specifiy method for column-permutation indeterminacy")
  }
  
  if (inherits(x, "var.boot")) {
    A_hat <- coef_x
  }
  else {
    A <- matrix(0, nrow = k, ncol = k * p)
    for (i in 1:k) {
      A[i, ] <- coef_x[[i]][1:(k * p), 1]
    }
    A_hat <- A
    if (type == "const") {
      v <- rep(1, k)
      for (i in 1:k) {
        v[i] <- coef_x[[i]][(k * p + 1), 1]
      }
      A_hat <- cbind(v, A)
    }
    else if (type == "trend") {
      trend <- rep(1, k)
      for (i in 1:k) {
        trend[i] <- coef_x[[i]][(k * p + 1), 1]
      }
      A_hat <- cbind(trend, A)
    }
    else if (type == "both") {
      v <- rep(1, k)
      for (i in 1:k) {
        v[i] <- coef_x[[i]][(k * p + 1), 1]
      }
      trend <- rep(1, k)
      for (i in 1:k) {
        trend[i] <- coef_x[[i]][(k * p + 2), 1]
      }
      A_hat <- cbind(v, trend, A)
    }
  }
  result <- list(B = P, A_hat = A_hat, method = "KernelICA", 
                 n = Tob, type = type, y = yOut, p = unname(p), K = k)
  class(result) <- "svars"
  return(result)
}


id.spec.gp <- function (x, ordering) 
{
  if (inherits(x, "var.boot")) {
    u <- x$residuals
    Tob <- nrow(u)
    k <- ncol(u)
    residY <- u
  } else {
    u <- residuals(x)
    Tob <- nrow(u)
    k <- ncol(u)
    residY <- u
  }
  if (inherits(x, "var.boot")) {
    p <- x$p
    y <- t(x$y)
    yOut = x$y
    type = x$type
    coef_x = x$coef_x
  } else if (inherits(x, "varest")) {
    p <- x$p
    y <- t(x$y)
    yOut = x$y
    type = x$type
    coef_x = coef(x)
  } else if (inherits(x, "nlVar")) {
    p <- x$lag
    y <- t(x$model[, 1:k])
    yOut <- x$model[, 1:k]
    coef_x <- t(coef(x))
    if (inherits(x, "VECM")) {
      coef_x <- t(VARrep(x))
    }
    if (rownames(coef_x)[1] %in% c("Intercept", "constant")) {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
    } else if (rownames(coef_x)[1] == "Trend") {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
    } 
    if (rownames(coef_x)[1] %in% c("Intercept", "constant", 
                                   "Trend")) {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
    }
    type <- x$include
    coef_x <- split(coef_x, rep(1:ncol(coef_x), each = nrow(coef_x)))
    coef_x <- lapply(coef_x, as.matrix)
  } else if (inherits(x, "list")) {
    p <- x$order
    y <- t(x$data)
    yOut <- x$data
    coef_x <- t(coef(x))
    coef_x <- x$coef
    if (x$cnst == TRUE) {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
      type = "const"
    }
    coef_x <- split(coef_x, rep(1:ncol(coef_x), each = nrow(coef_x)))
    coef_x <- lapply(coef_x, as.matrix)
  } else if (inherits(x, "vec2var")) {
    coef_x <- vector("list", length = k)
    names(coef_x) <- colnames(x$y)
    p <- x$p
    y <- t(x$y)
    yOut <- x$y
    for (i in seq_len(k)) {
      for (j in seq_len(p)) coef_x[[i]] <- c(coef_x[[i]], 
                                             x$A[[j]][i, ])
      coef_x[[i]] <- c(coef_x[[i]], x$deterministic[i, 
      ])
    }
    coef_x <- lapply(coef_x, matrix)
    type <- "const"
  } else {
    stop("Object class is not supported")
  }
  
  sigg <- crossprod(u)/(Tob - 1 - k * p)
  
  ## Spectral Decomposition
  spec_res             <- eigen(sigg)
  P                    <- spec_res$vectors
  D                    <- diag(spec_res$values)
  W                    <- P %*% sqrt(D) %*% t(P)
  Wscal.spec           <- rescaleVar(W_hat = W, ut = u)$Ws
  Wscaled.spec         <- Wscal.spec %*% solve(P)
  Praw                 <- solve(Wscaled.spec)
  
  if (ordering == "frob") {
    P <- Praw %*% myfrob(A.hat = Praw, A = x$B_original)
  } else if (ordering == "maxfinder") {
    P <- maxfinder(A = Praw)$A.id
  } else {
    stop("Specify method for column-permutation indeterminacy")
  }
  
  if (inherits(x, "var.boot")) {
    A_hat <- coef_x
  } else {
    A <- matrix(0, nrow = k, ncol = k * p)
    for (i in 1:k) {
      A[i, ] <- coef_x[[i]][1:(k * p), 1]
    }
    A_hat <- A
    if (type == "const") {
      v <- rep(1, k)
      for (i in 1:k) {
        v[i] <- coef_x[[i]][(k * p + 1), 1]
      }
      A_hat <- cbind(v, A)
    } else if (type == "trend") {
      trend <- rep(1, k)
      for (i in 1:k) {
        trend[i] <- coef_x[[i]][(k * p + 1), 1]
      }
      A_hat <- cbind(trend, A)
    } else if (type == "both") {
      v <- rep(1, k)
      for (i in 1:k) {
        v[i] <- coef_x[[i]][(k * p + 1), 1]
      }
      trend <- rep(1, k)
      for (i in 1:k) {
        trend[i] <- coef_x[[i]][(k * p + 2), 1]
      }
      A_hat <- cbind(v, trend, A)
    }
  }
  
  result <- list(B = P, A_hat = A_hat, method = "Spectral Decomposition", 
                 n = Tob, type = type, y = yOut, p = unname(p), K = k)
  class(result) <- "svars"
  return(result)
}


id.pml.gp <- function (x, ww, ordering) 
{
  p.start <- as.vector(ww)
  if (inherits(x, "var.boot")) {
    u <- x$residuals
    Tob <- nrow(u)
    k <- ncol(u)
    residY <- u
  }
  else {
    u <- residuals(x)
    Tob <- nrow(u)
    k <- ncol(u)
    residY <- u
  }
  if (inherits(x, "var.boot")) {
    p <- x$p
    y <- t(x$y)
    yOut = x$y
    type = x$type
    coef_x = x$coef_x
  }
  else if (inherits(x, "varest")) {
    p <- x$p
    y <- t(x$y)
    yOut = x$y
    type = x$type
    coef_x = coef(x)
  }
  else if (inherits(x, "nlVar")) {
    p <- x$lag
    y <- t(x$model[, 1:k])
    yOut <- x$model[, 1:k]
    coef_x <- t(coef(x))
    if (inherits(x, "VECM")) {
      coef_x <- t(VARrep(x))
    }
    if (rownames(coef_x)[1] %in% c("Intercept", "constant")) {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
    }
    else if (rownames(coef_x)[1] == "Trend") {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
    }
    if (rownames(coef_x)[1] %in% c("Intercept", "constant", 
                                   "Trend")) {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
    }
    type <- x$include
    coef_x <- split(coef_x, rep(1:ncol(coef_x), each = nrow(coef_x)))
    coef_x <- lapply(coef_x, as.matrix)
  }
  else if (inherits(x, "list")) {
    p <- x$order
    y <- t(x$data)
    yOut <- x$data
    coef_x <- t(coef(x))
    coef_x <- x$coef
    if (x$cnst == TRUE) {
      coef_x <- coef_x[c(2:nrow(coef_x), 1), ]
      type = "const"
    }
    coef_x <- split(coef_x, rep(1:ncol(coef_x), each = nrow(coef_x)))
    coef_x <- lapply(coef_x, as.matrix)
  }
  else if (inherits(x, "vec2var")) {
    coef_x <- vector("list", length = k)
    names(coef_x) <- colnames(x$y)
    p <- x$p
    y <- t(x$y)
    yOut <- x$y
    for (i in seq_len(k)) {
      for (j in seq_len(p)) coef_x[[i]] <- c(coef_x[[i]], 
                                             x$A[[j]][i, ])
      coef_x[[i]] <- c(coef_x[[i]], x$deterministic[i, 
      ])
    }
    coef_x <- lapply(coef_x, matrix)
    type <- "const"
  }
  else {
    stop("Object class is not supported")
  }
  sigg <- crossprod(u)/(Tob - 1 - k * p)
  P_chol <- t(chol(sigg))
  
  u_cholb <- t(solve(P_chol) %*% t(u))
  tu_cholb<- t(u_cholb)
  #rm(u_chol,tu_chol)
  n <- length(x = u_cholb[,1])
  if(k==5){
    pseudo.log.lik5D.super        <- function(c){
      c11 <- c[1]
      c21 <- c[2]
      c31 <- c[3] 
      c41 <- c[4]
      c51 <- c[5]
      c12 <- c[6]
      c22 <- c[7]
      c32 <- c[8]
      c42 <- c[9]
      c52 <- c[10]
      c13 <- c[11]
      c23 <- c[12]
      c33 <- c[13]
      c43 <- c[14]
      c53 <- c[15]
      c14 <- c[16]
      c24 <- c[17]
      c34 <- c[18]
      c44 <- c[19]
      c54 <- c[20]
      c15 <- c[21]
      c25 <- c[22]
      c35 <- c[23]
      c45 <- c[24]
      c55 <- c[25]
      
      c.1y <- c11*u_chol[,1] + c21*u_chol[,2] + c31*u_chol[,3] + c41*u_chol[,4] + c51*u_chol[,5]
      c.2y <- c12*u_chol[,1] + c22*u_chol[,2] + c32*u_chol[,3] + c42*u_chol[,4] + c52*u_chol[,5]
      c.3y <- c13*u_chol[,1] + c23*u_chol[,2] + c33*u_chol[,3] + c43*u_chol[,4] + c53*u_chol[,5]
      c.4y <- c14*u_chol[,1] + c24*u_chol[,2] + c34*u_chol[,3] + c44*u_chol[,4] + c54*u_chol[,5]
      c.5y <- c15*u_chol[,1] + c25*u_chol[,2] + c35*u_chol[,3] + c45*u_chol[,4] + c55*u_chol[,5]
      #P <- 1/(2 * p^(1/p - 1))*gamma(1/p)
      g <- rep(NA,1)
      g[1] <- -1*(n*log(1/sqrt(2*pi))-sum((c.1y)^2/2) + 
                    n*log(gamma((5+1)/2)/(sqrt(5*pi)*gamma(5/2)))-((5+1)/2)*sum(log(1+(c.2y)^2/5)) +
                    n*log(gamma((6+1)/2)/(sqrt(6*pi)*gamma(6/2)))-((6+1)/2)*sum(log(1+(c.3y)^2/6)) +
                    n*log(gamma((7+1)/2)/(sqrt(7*pi)*gamma(7/2)))-((7+1)/2)*sum(log(1+(c.4y)^2/7)) + 
                    n*log(gamma((5+1)/2)/(sqrt(5*pi)*gamma(5/2)))-((5+1)/2)*sum(log(1+(c.5y)^2/5)))  
      g
    }
    pseudo.log.lik5D.sub          <- function(c) {
      c11 <- c[1]
      c21 <- c[2]
      c12 <- c[3]
      c22 <- c[4]
      
      c.1y <- c11*u_chol[[i]][,1] + c21*u_chol[[i]][,2]
      c.2y <- c12*u_chol[[i]][,1] + c22*u_chol[[i]][,2]
      
      #P <- 1/(2 * p^(1/p - 1))*gamma(1/p)
      g <- rep(NA,1)
      g[1] <- -1*(n*log(1/sqrt(2*pi)) + sum(pi*(c.1y)^2+log(cosh((pi/2)*c.1y))) + 
                    n*log(1/sqrt(2*pi)) + sum(pi*(c.2y)^2+log(cosh((pi/2)*c.2y))))
      g
    }
    orth.constr5D                 <- function(c){
      c11 <- c[1]
      c21 <- c[2]
      c31 <- c[3] 
      c41 <- c[4]
      c51 <- c[5]
      c12 <- c[6]
      c22 <- c[7]
      c32 <- c[8]
      c42 <- c[9]
      c52 <- c[10]
      c13 <- c[11]
      c23 <- c[12]
      c33 <- c[13]
      c43 <- c[14]
      c53 <- c[15]
      c14 <- c[16]
      c24 <- c[17]
      c34 <- c[18]
      c44 <- c[19]
      c54 <- c[20]
      c15 <- c[21]
      c25 <- c[22]
      c35 <- c[23]
      c45 <- c[24]
      c55 <- c[25]
      g <- rep(NA,15)
      g[1] <- c11*c12 + c21*c22 + c31*c32 + c41*c42 + c51*c52 # c.1*c.2
      g[2] <- c11*c13 + c21*c23 + c31*c33 + c41*c43 + c51*c53 # c.1*c.3
      g[3] <- c11*c14 + c21*c24 + c31*c34 + c41*c44 + c51*c54 # c.1*c.4
      g[4] <- c11*c15 + c21*c25 + c31*c35 + c41*c45 + c51*c55 # c.1*c.5
      g[5] <- c12*c13 + c22*c23 + c32*c33 + c42*c43 + c52*c53 # c.2*c.3
      g[6] <- c12*c14 + c22*c24 + c32*c34 + c42*c44 + c52*c54 # c.2*c.4
      g[7] <- c12*c15 + c22*c25 + c32*c35 + c42*c45 + c52*c55 # c.2*c.5
      g[8] <- c13*c14 + c23*c24 + c33*c34 + c43*c34 + c52*c54 # c.3*c.4
      g[9] <- c13*c15 + c23*c25 + c33*c35 + c43*c35 + c52*c55 # c.3*c.5
      g[10]<- c14*c15 + c24*c25 + c34*c35 + c44*c35 + c54*c55 # c.4*c.5
      
      g[11] <- c11^2 + c21^2 + c31^2 + c41^2 + c51^2- 1 # c.1*c.1
      g[12] <- c12^2 + c22^2 + c32^2 + c41^2 + c52^2- 1 # c.2*c.2
      g[13] <- c13^2 + c23^3 + c33^2 + c43^2 + c53^2- 1
      g[14] <- c14^2 + c24^2 + c34^2 + c44^2 + c54^2- 1
      g[15] <- c15^2 + c25^2 + c35^2 + c45^2 + c55^2- 1 
      g
    }
    
    please <- auglag(par = c(improvedLHS(n = 1,k = 25)), fn = pseudo.log.lik5D.super, heq = orth.constr5D,
                     control.outer = list(trace=F))$par
    A            <- matrix(data = please,nrow = 5,byrow = F)
  }else if(k==3){
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
      
      c.1y <- c11*tu_cholb[,1] + c21*tu_cholb[,2] + c31*tu_cholb[,3]
      c.2y <- c12*tu_cholb[,1] + c22*tu_cholb[,2] + c32*tu_cholb[,3]
      c.3y <- c13*tu_cholb[,1] + c23*tu_cholb[,2] + c33*tu_cholb[,3]
      
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
      
      c.1y <- c11*tu_cholb[,1] + c21*tu_cholb[,2] + c31*tu_cholb[,3]
      c.2y <- c12*tu_cholb[,1] + c22*tu_cholb[,2] + c32*tu_cholb[,3]
      c.3y <- c13*tu_cholb[,1] + c23*tu_cholb[,2] + c33*tu_cholb[,3]
      
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
    please                      <- auglag(par = p.start, fn = pseudo.log.lik3D.super, heq = orth.constr3D,
                                          control.outer = list(trace=F))$par 
    A_mixing                    <- matrix(data = please,nrow = 3,byrow = F)
  }
  
  Praw         <- (P_chol) %*% (A_mixing)
  if (ordering == "frob") {
    P      <- Praw %*% myfrob(A.hat = Praw,A = x$B_original)
  }else if(ordering == "maxfinder"){
    P      <- maxfinder(A = Praw)$A.id
  }else{
    stop("Specifiy method for column-permutation indeterminacy")
  }
  
  if (inherits(x, "var.boot")) {
    A_hat <- coef_x
  }
  else {
    A <- matrix(0, nrow = k, ncol = k * p)
    for (i in 1:k) {
      A[i, ] <- coef_x[[i]][1:(k * p), 1]
    }
    A_hat <- A
    if (type == "const") {
      v <- rep(1, k)
      for (i in 1:k) {
        v[i] <- coef_x[[i]][(k * p + 1), 1]
      }
      A_hat <- cbind(v, A)
    }
    else if (type == "trend") {
      trend <- rep(1, k)
      for (i in 1:k) {
        trend[i] <- coef_x[[i]][(k * p + 1), 1]
      }
      A_hat <- cbind(trend, A)
    }
    else if (type == "both") {
      v <- rep(1, k)
      for (i in 1:k) {
        v[i] <- coef_x[[i]][(k * p + 1), 1]
      }
      trend <- rep(1, k)
      for (i in 1:k) {
        trend[i] <- coef_x[[i]][(k * p + 2), 1]
      }
      A_hat <- cbind(v, trend, A)
    }
  }
  result <- list(B = P, A_hat = A_hat, method = "PML", 
                 n = Tob, type = type, y = yOut, p = unname(p), K = k)
  class(result) <- "svars"
  return(result)
}

imrf.gianluca <- function (y, B_hat, A_hat, type,horizon = 20) 
{## y must be n*k observations
  # A_hat is taken from getcoefSVARS
  # B_hat the mixing matrix
  "%^%" <- function(A, n) {
    if (n == 1) {
      A
    }
    else {
      A %*% (A %^% (n - 1))
    }
  }
  IrF <- function(A_hat, B_hat, horizon) {
    k <- nrow(A_hat)
    p <- ncol(A_hat)/k
    if (p == 1) {
      irfa <- array(0, c(k, k, horizon))
      irfa[, , 1] <- B_hat
      for (i in 1:horizon) {
        irfa[, , i] <- (A_hat %^% i) %*% B_hat
      }
      return(irfa)
    }
    else {
      irfa <- array(0, c(k, k, horizon))
      irfa[, , 1] <- B_hat
      Mm <- matrix(0, nrow = k * p, ncol = k * p)
      Mm[1:k, 1:(k * p)] <- A_hat
      Mm[(k + 1):(k * p), 1:((p - 1) * k)] <- diag(k * 
                                                     (p - 1))
      Mm1 <- diag(k * p)
      for (i in 1:(horizon - 1)) {
        Mm1 <- Mm1 %*% Mm
        irfa[, , (i + 1)] <- Mm1[1:k, 1:k] %*% B_hat
      }
      return(irfa)
    }
  }
  if (type == "const") {
    A_hat <- A_hat[, -1]
  }
  else if (type == "trend") {
    A_hat <- A_hat[, -1]
  }
  else if (type == "both") {
    A_hat <- A_hat[, -c(1, 2)]
  }
  else {
    A_hat <- A_hat
  }
  B_hat <- B_hat
  IR <- IrF(A_hat, B_hat, horizon)
  impulse <- matrix(0, ncol = dim(IR)[2]^2 + 1, nrow = dim(IR)[3])
  colnames(impulse) <- rep("V1", ncol(impulse))
  cc <- 1
  impulse[, 1] <- seq(1, dim(IR)[3])
  for (i in 1:dim(IR)[2]) {
    for (j in 1:dim(IR)[2]) {
      cc <- cc + 1
      impulse[, cc] <- IR[i, j, ]
      colnames(impulse)[cc] <- paste("epsilon[", colnames(y)[j], 
                                     "]", "%->%", colnames(y)[i])
    }
  }
  impulse <- list(irf = as.data.frame(impulse))
  B_hat   <- B_hat
  class(impulse) <- "irf"
  results <- list(impulse = impulse, B_hat = B_hat)
  return(results)
}

maxfinder <- function(A){
  ## A is the estimated mixing matrix
  allperms <- gtools::permutations(n=nrow(A),r=ncol(A))
  nperms   <- nrow(allperms)
  A.perm   <- list()
  idx      <- list()
  bestval  <- Inf;
  besti    <- 0;
  failures <- 0
  for (i in 1:nperms){
    Pr          <- diag(nrow(A))
    Pr          <- Pr[,allperms[i,]]
    A.perm[[i]] <- A %*% Pr
    idx[[i]]      <- A.perm[[i]] %>% apply(MARGIN = 2, function(x) which.max(abs(x)))
  }
  daje <- seq(from = 1,by = 1,length.out = nrow(A))
  ## Finding the permutation matrix tha has the maximum in the rows 1,2,3 or
  ## alternatively, in the closest order possible [(1,2,1) or (1,1,1)]
  index <- idx %>% detect_index(.f = ~ all(. == daje))
  if (index == 0) {
    failures <-1
    warning('Maxima are not on different rows \n')
    for (i in 1:nperms){
      Pr          <- diag(nrow(A))
      Pr          <- Pr[allperms[i,],]
      Aperm        <- A %*% Pr
      c           <- nzdiagscore(Aperm)
      if (c<bestval){
        A.id <- Aperm
        bestval <- c
        besti <- i
      }
    }
  }else{A.id <- A.perm[[index]]}
  for (i in 1:nrow(A.id)){
    if (diag(A.id)[i]<0) A.id[,i] <- -A.id[,i]
  }
  ls <- list(fail = failures, A.id = A.id)
  return(ls)
}


mb.boot.gigi <- function (x, b.length = 15, n.ahead = 20, horizon, nboot, w00 = diag(x$K), method,nc = 1, dd = NULL, 
                          signrest = NULL, itermax = 300, steptol = 200, iter2 = 50,set_seed = F, ordering) 
{
  if (method == "CvM" & is.null(dd)) {
    dd <- copula::indepTestSim(x$n, x$K, verbose = F)
  }
  
  
  y_lag_cr <- function(y, lag_length) {
    y_lag <- matrix(NA, dim(y)[1], dim(y)[2] * lag_length)
    for (i in 1:lag_length) {
      y_lag[(1 + i):dim(y)[1], ((i * NCOL(y) - NCOL(y)) + 
                                  1):(i * NCOL(y))] <- y[1:(dim(y)[1] - i), (1:NCOL(y))]
    }
    y_lag <- as.matrix(y_lag[-(1:lag_length), ])
    out <- list(lags = y_lag)
  }
  
  
  sqrt.f <- function(Pstar, Sigma_u_star) {
    yy <- suppressMessages(sqrtm(Sigma_u_hat_old)) %*% solve(suppressMessages(sqrtm(Sigma_u_star))) %*% 
      Pstar
    return(yy)
  }
  
  
  y <- x$y
  p <- x$p
  obs <- x$n
  k <- x$K
  B <- x$B
  if (length(signrest) > k) {
    stop("too many sign restrictions")
  }
  A <- x$A_hat
  Z <- t(y_lag_cr(y, p)$lags)
  if (x$type == "const") {
    Z <- rbind(rep(1, ncol(Z)), Z)
  } else if (x$type == "trend") {
    Z <- rbind(seq(1, ncol(Z)), Z)
  } else if (x$type == "both") {
    Z <- rbind(rep(1, ncol(Z)), seq(1, ncol(Z)), Z)
  } else {
    Z <- Z
  }
  
  
  u <- t(y[-c(1:p), ]) - A %*% Z
  Sigma_u_hat_old <- tcrossprod(u)/(obs - 1 - k * p)
  errors <- list()
  N <- obs/b.length
  blocks <- array(NA, c(b.length, k, obs - b.length + 1))
  u <- t(u)
  for (i in 0:(obs - b.length)) {
    blocks[, , (i + 1)] <- u[(i + 1):(i + b.length), ]
  }
  for (i in 1:nboot) {
    epsilon.star <- matrix(0, b.length * ceiling(N), ncol(u))
    epsilon.star <- list()
    
    for (kk in 1:ceiling(N)) {
      
      epsilon.star[[kk]] <- blocks[, , floor(runif(1, 1, 
                                                   obs - b.length + 2))]
      
    }
    epsilon.star <- do.call("rbind", epsilon.star)
    for (s in 1:b.length) {
      b.mean <- colSums(epsilon.star[1:(s + (obs - b.length)), 
      ])/(obs - b.length + 1)
      for (j in 0:floor(N)) {
        epsilon.star[j * b.length + s, ] <- epsilon.star[j * 
                                                           b.length + s, ] - b.mean
      }
    }
    epsilon.star <- epsilon.star[1:obs, ]
    errors[[i]] <- t(epsilon.star)
  }
  
  
  bootf <- function(Ustar1){
    Ystar <- matrix(0, nrow(y), k)
    Ystar[1:p, ] <- y[1:p, ]
    if (x$type == "const" | x$type == "trend") {
      for (i in (p + 1):nrow(y)) {
        for (j in 1:k) {
          Ystar[i, j] <- A[j, 1] + A[j, -1] %*% c(t(Ystar[(i - 
                                                             1):(i - p), ])) + Ustar1[j, (i - p)]
        }
      }
    } else if (x$type == "both") {
      
      for (i in (p + 1):nrow(y)) {
        for (j in 1:k) {
          ## This works only with a VAR(4) in three variables with constant and trend
          Ystar[i,j] <- A[j, c(1)] + A[j, c(2)]*(i-p) + A[j,c(3:5)] %*% Ystar[i-1,] +
            A[j,c(6:8)] %*% Ystar[i-2,] +
            A[j,c(9:11)] %*% Ystar[i-3,] +
            A[j,c(12:14)] %*% Ystar[i-4,] + Ustar1[j,(i-p)]
        }
      }
    } else if (x$type == "none") {
      for (i in (p + 1):nrow(y)) {
        for (j in 1:k) {
          Ystar[i, j] <- A[j, ] %*% c(t(Ystar[(i - p):(i - 
                                                         1), ])) + Ustar1[j, (i - p)]
        }
      }
    }
    colnames(Ystar) <- colnames(y)
    varb <- suppressWarnings(VAR(Ystar, p = x$p, type = x$type))
    varb$B_original <- x$B
    class(varb) <- "varest"
    Ustar <- residuals(varb)
    Sigma_u_star <- crossprod(Ustar)/(obs - 1 - k * p)
    if (method == "Non-Gaussian maximum likelihood") {
      temp <- id.ngml_boot(varb, stage3 = x$stage3, Z = Z)
    } else if (method == "Changes in Volatility") {
      temp <- tryCatch(id.cv_boot(varb, SB = x$SB, Z = Z, 
                                  restriction_matrix = x$restriction_matrix), error = function(e) NULL)
    } else if (method == "CvM") {
      temp <- id.cvm.gp(varb, itermax = itermax, steptol = steptol, 
                        iter2 = iter2, dd,ordering = ordering)
      
    } else if (method == "Distance covariances") {
      ## here there is identification via maxfinder and same init
      temp <- id.dc.gp(varb, PIT = F,ww = w00,ordering = ordering) 
    } else if (method == "fastICA"){
      temp <- id.fastICA.gp(x = varb,ww = w00,ordering = ordering)
    } else if (method == "KernelICA"){
      temp <- id.KernelICA.gp(x = varb,ww = w00,ordering = ordering)
    } else if (method == "PML"){
      temp <- id.pml.gp(x = varb,ww = w00,ordering = ordering)
    } else if (method == "spec"){
    temp <- id.spec.gp(x = varb,ww = w00,ordering = ordering)
  }else {
      temp <- tryCatch(id.st_boot(varb, c_fix = x$est_c, 
                                  transition_variable = x$transition_variable, 
                                  restriction_matrix = x$restriction_matrix, gamma_fix = x$est_g, 
                                  max.iter = x$iteration, crit = 0.01, Z = Z), 
                       error = function(e) NULL)
    }
    if (!is.null(temp)) {
      Pstar <- temp$B
      if (!is.null(x$restriction_matrix)) {
        #Pstar1 <- Pstar
        #frobP <- frobICA_mod(Pstar1, B, standardize = TRUE)
        #Pstar1 <- maxfinder(A = Pstar)$A.id
      }
      else {
        #Pstar1 <- sqrt.f(Pstar, Sigma_u_star)
        #diag_sigma_root <- diag(diag(suppressMessages(sqrtm(Sigma_u_hat_old))))
        #frobP <- frobICA_mod(t(solve(diag_sigma_root) %*% 
        #                         Pstar1), t(solve(diag_sigma_root) %*% B), standardize = TRUE)
      }
      #Pstar <- Pstar1 %*% frobP$perm
      #temp$B <- Pstar
      #Pstar <- Pstar1
      ip.tot <- imrf.gianluca(y = temp$y,A_hat = temp$A_hat,
                              B_hat = Pstar, horizon = horizon, type = temp$type)
      ip     <- ip.tot$impulse
      Pstar  <- ip.tot$B_hat
      return(list(ip, Pstar))
    }
    else {
      return(NA)
    }
  }
  
  
  bootstraps <- pblapply(errors, bootf, cl = nc)
  
  
  delnull <- function(x) {
    x[unlist(lapply(x, length) != 0)]
  }
  bootstraps <- lapply(bootstraps, function(x) x[any(!is.na(x))])
  bootstraps <- delnull(bootstraps)
  Bs <- array(0, c(k, k, length(bootstraps)))
  ipb <- list()
  for (i in 1:length(bootstraps)) {
    Bs[, , i] <- bootstraps[[i]][[2]]
    ipb[[i]] <- bootstraps[[i]][[1]]
  }
  v.b <- matrix(Bs, ncol = k^2, byrow = T)
  cov.bs <- cov(v.b)
  if (method == "CvM" | method == 
      "Distance covariances" | method == "fastICA"| method == "PML"| method == "KernelICA") {
    SE <- matrix(sqrt(diag(cov.bs)), k, k)
    rownames(SE) <- rownames(x$B)
  }
  else {
    SE <- NULL
  }
  boot.mean <- matrix(colMeans(v.b), k, k)
  boot.median <- matrix(apply(v.b, 2, median), nrow = k, ncol = k)
  rownames(boot.mean) <- rownames(x$B)
  rownames(boot.median) <- rownames(x$B)
  if (!is.null(x$restriction_matrix)) {
    if (!is.null(signrest)) {
      cat("Testing signs only possible for unrestricted model \n")
    }
    sign.part <- NULL
    sign.complete <- NULL
  }
  else {
    if (is.null(signrest)) {
      sign.mat <- matrix(FALSE, nrow = k, ncol = k)
      sign.complete <- 0
      sign.part <- rep(0, times = k)
      for (i in 1:length(bootstraps)) {
        pBs <- permutation(Bs[, , i])
        sign.mat <- lapply(pBs, function(z) {
          sapply(1:k, function(ii) {
            all(z[, ii]/abs(z[, ii]) == x$B[, ii]/abs(x$B[, 
                                                          ii])) | all(z[, ii]/abs(z[, ii]) == x$B[, 
                                                                                                  ii]/abs(x$B[, ii]) * (-1))
          })
        })
        if (any(unlist(lapply(sign.mat, function(sign.mat) all(sign.mat == 
                                                               TRUE))))) {
          sign.complete <- sign.complete + 1
        }
        for (j in 1:k) {
          check <- rep(FALSE, k)
          for (l in 1:k) {
            check[l] <- any(all(pBs[[1]][, l]/abs(pBs[[1]][, 
                                                           l]) == x$B[, j]/abs(x$B)[, j]) | all(pBs[[1]][, 
                                                                                                         l]/abs(pBs[[1]][, l]) == x$B[, j]/abs(x$B)[, 
                                                                                                                                                    j] * (-1)))
          }
          if (sum(check) == 1) {
            sign.part[[j]] <- sign.part[[j]] + 1
          }
        }
      }
    }
    else {
      nrest <- length(signrest)
      sign.part <- rep(list(0), nrest)
      sign.complete <- 0
      for (j in 1:length(bootstraps)) {
        check.full <- 0
        for (i in 1:nrest) {
          check <- rep(FALSE, length(signrest[[i]][!is.na(signrest[[i]])]))
          for (l in 1:k) {
            check[l] <- any(all(Bs[!is.na(signrest[[i]]), 
                                   l, j]/abs(Bs[!is.na(signrest[[i]]), l, 
                                                j]) == signrest[[i]][!is.na(signrest[[i]])]) | 
                              all(Bs[!is.na(signrest[[i]]), l, j]/abs(Bs[!is.na(signrest[[i]]), 
                                                                         l, j]) == signrest[[i]][!is.na(signrest[[i]])] * 
                                    (-1)))
          }
          if (sum(check) == 1) {
            sign.part[[i]] <- sign.part[[i]] + 1
            check.full <- check.full + 1
          }
        }
        if (check.full == nrest) {
          sign.complete <- sign.complete + 1
        }
      }
      names(sign.part) <- names(signrest)
    }
  }
  ip <- imrf.gianluca(y = x$y ,B_hat = boot.mean, 
                      horizon = horizon,type = x$type,A_hat = A)$impulse
  result <- list(true = ip, bootstrap = ipb, SE = SE, nboot = nboot, 
                 point_estimate = x$B, boot_mean = boot.mean, boot_median = boot.median,
                 signrest = signrest, sign_complete = sign.complete, sign_part = sign.part, 
                 cov_bs = cov.bs, method = "Moving Block")
  class(result) <- "sboot"
  return(result)
}


mov.block.boot <- function(X, lagVar, B, C = C, w.init = diag(nrow(W)), 
                           nrboot, method, display = T, quantiles, b.length = 15, ordering,set_seed=F){
  
  u   <- X
  obs <- nrow(u)
  k   <- ncol(u)
  N <- obs/b.length
  blocks <- array(NA, c(b.length, k, obs - b.length + 1))
  Sigma_u_hat_old <- tcrossprod(X)/(obs - 1 - k * lagVar)
  #C <- t(chol(Sigma_u_hat_old))
  errors <- list()
  save_W <- array(NA,dim=c(nrboot,k^2))
  save_B <- array(NA,dim=c(nrboot,k^2))
  for (i in 0:(obs - b.length)) {
    blocks[, , (i + 1)] <- u[(i + 1):(i + b.length), ]
  }
  for (i in 1:nrboot) {
    if ((i%%100 == 0) & display) {cat("bootstrap run", i, "out of", nrboot, "\n")}
    epsilon.star <- matrix(0, b.length * ceiling(N), ncol(u))
    epsilon.star <- list()
    jj <- 1
    SEED <- seq(1,nrboot*ceiling(N))
    for (kk in 1:ceiling(N)) {
      if (set_seed == T) {
        set.seed(SEED[jj])
      }
      epsilon.star[[kk]] <- blocks[, , floor(runif(1, 1, 
                                                   obs - b.length + 2))]
      jj <- jj + 1
    }
    epsilon.star <- do.call("rbind", epsilon.star)
    for (s in 1:b.length) {
      b.mean <- colSums(epsilon.star[1:(s + (obs - b.length)), 
      ])/(obs - b.length + 1)
      for (j in 0:floor(N)) {
        epsilon.star[j * b.length + s, ] <- epsilon.star[j * 
                                                           b.length + s, ] - b.mean
      }
    }
    epsilon.star <- epsilon.star[1:obs, ]
    #errors[[i]] <- t(epsilon.star)
    X_boot      <- epsilon.star
    
    if (method=="fastICA"){
      icares              <- fastICA(X_boot, ncol(B), tol=1e-6, w.init=w.init )
      W.ica               <- t((icares$K) %*% (icares$W))
      W.scal.ica          <- rescaleVar(W_hat = W.ica, ut = t(X_boot))$Ws
      icares_MD           <- W.scal.ica %*% solve(C)
      W_boot              <- icares_MD
      if(ordering == "max"){
        B_boot              <- maxfinder(A = solve(W_boot))$A.id
      } else if(ordering == "frob") B_boot         <- solve(W_boot) %*% myfrob(A.hat = solve(W_boot),A = B)
      
      B_boot[,1]*1/B[1,1] -> B_boot[,1]
      B_boot[,2]*1/B[2,2] -> B_boot[,2]
      B_boot[,3]*1/B[3,3] -> B_boot[,3]
    }else if (method=="KernelICA"){
      
      
      
      ##KernelICA
      
      # Use kernel_ica method
      kica_init[[kk]] <- KernelICA::kernel_ica(u_chol, variant = "kgv", kernel = "gauss", eps = 1e-14, init = list(winit[[kk]]))
      
      
      # Given the nature of DCOV estimation, the matrix C%*%W is the estimate of the mixing matrix
      # So, in order to obtain the estimate of the unmixing matrix we have to calculate the following:
      
      
      Wscaled_kica[[kk]] <- rescaleVar(W_hat = kica_init[[kk]]$W, ut = t(u_chol))$Ws %*% solve(C)
      
      # If you have any other post-processing tasks similar to the one done in the provided code, 
      # like the rescaleVar for fastICA, you can do them here.
      
      
      
      
      W_boot              <- kicares_MD
      if(ordering == "max"){
        B_boot              <- maxfinder(A = solve(W_boot))$A.id
      } else if(ordering == "frob") B_boot         <- solve(W_boot) %*% myfrob(A.hat = solve(W_boot),A = B)
      
      B_boot[,1]*1/B[1,1] -> B_boot[,1]
      B_boot[,2]*1/B[2,2] -> B_boot[,2]
      B_boot[,3]*1/B[3,3] -> B_boot[,3]
    }else if(method == "DCov"){
      DCov                <- steadyICA(X_boot, n.comp = ncol(B), w.init = w.init, 
                                       PIT = T, symmetric = F,maxit = 1000,verbose = F)
      #DCov                <- dcovICA(Z = X_boot,theta.0 = w.init)
      DC_MD               <- solve(DCov$W) %*% solve(C)
      W_boot              <- DC_MD
      if(ordering == "max"){
        B_boot              <- maxfinder(A = solve(W_boot))$A.id
      } else if(ordering == "frob") B_boot         <- solve(W_boot) %*% myfrob(A.hat = solve(W_boot),A = B)
      B_boot[,1]*1/B[1,1] -> B_boot[,1]
      B_boot[,2]*1/B[2,2] -> B_boot[,2]
      B_boot[,3]*1/B[3,3] -> B_boot[,3]
      
    }else if(method == "CvM"){
      if (is.null(dd)) {
        dd <- indepTestSim(n, ncol(X), verbose = F)
      }
      lower <- rep(0, ncol(X) * (ncol(X) - 1)/2)
      upper <- rep(pi, ncol(X) * (ncol(X) - 1)/2)
      de_control <- list(itermax = 500, steptol = 100, 
                         trace = FALSE)
      ## First step of optimization with DEoptim
      
      de_res <- DEoptim(testlik, lower = lower, upper = upper, 
                        control = de_control, faklow = C, u = X_boot, dd = dd) #U or u_chol
      k <- ncol(X)
      iter2  <- 75
      ## Second step of optimization. Creating randomized starting angles around the optimized angles
      ## here maybe you insert the LHS stuff
      theta_rot <- matrix(rnorm(n = (k * (k - 1)/2) * iter2, mean = de_res$optim$bestmem, 
                                sd = 0.3), (k * (k - 1)/2), iter2)
      theta_rot <- cbind(theta_rot, de_res$optim$bestmem)
      # Start vectors for iterative optimization approach
      startvec_list <- as.list(as.data.frame(theta_rot))
      erg_list <- lapply(X = startvec_list, FUN = optim, fn = testlik, 
                         gr = NULL, faklow = C, u = X_boot, dd = dd, 
                         method = ifelse(k == 2, "Brent", "Nelder-Mead"), 
                         lower = ifelse(k == 2, -.Machine$double.xmax/2, -Inf), 
                         upper = ifelse(k == 2, .Machine$double.xmax/2, Inf), 
                         control = list(maxit = 1000), 
                         hessian = FALSE)
      #save(erg_list, file = 'erglist.RData')
      #load(file = 'erglist.RData')
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
      Acvm                <- rotmat(par_o, C)
      if(ordering == "max"){
        B_boot              <- maxfinder(A = Acvm)$A.id
      } else if(ordering == "frob") B_boot         <- Acvm %*% myfrob(A.hat = Acvm,A = B)
      B_boot[,1]*1/B[1,1] -> B_boot[,1]
      B_boot[,2]*1/B[2,2] -> B_boot[,2]
      B_boot[,3]*1/B[3,3] -> B_boot[,3]
    }
    W_boot <- solve(B_boot)
    
    
    # need to row-permute matrix for fitting it best to W
    #W_perm <- best_permutation(W,W_boot)
    #print(W_perm)
    
    #save_W[i,] <- W_perm
    save_B[i,] <- as.vector(B_boot)
    
    
  }
  B.boot <- vector("list",length = nrboot)
  test.t <- array(data = NA,dim = c(ncol(C),ncol(C),length(quantiles)),
                  dimnames = list(paste0(seq(1,ncol(C),1)),paste0(seq(1,ncol(C),1)), paste0(quantiles)))
  
  for (j in 1:nrboot){
    B.boot[[j]] <- (matrix(data = save_B[j,],nrow = ncol(C),ncol = ncol(C),byrow = F))
  }
  cnt <- 0
  for (alpha in quantiles){
    
    cnt <- cnt +1
    
    for (r in 1:nrow(C)){
      for (c in 1:ncol(C)){
        ## create columns with bootstrap estimates of the pmt
        eval(parse(text=paste("b",r,c,".boot<-do.call(rbind,lapply(B.boot,'[',",r,",",c,"))", sep="")))
        ## create columns with bootstrap estimates of b_boot - b_true, b_true = 0
        eval(parse(text=paste("XX",r,c," <- b",r,c,".boot",sep="")))
        a <- quantile(x = eval(parse(text=paste("XX",r,c,sep = ""))),probs = c(alpha/2, (1-alpha/2)))
        if (0 < min(a) | 0 > max(a)){
          test.t[r,c,cnt] <- TRUE ##  reject = 0
        }else test.t[r,c,cnt] <- FALSE 
        
      }
    }
  }
  ls   <- list(w.pmts = save_W, B_boot = B.boot, test = test.t)
}



myboot <- function(X, W, C = C, Bt, w.init = diag(nrow(W)), nrboot, method, display = T,seed = NULL) {
  ## Bootstrap on mixing matrix pmts
  ## X are the orthogonalized reduced-form residuals
  ## W is the unmixing matrix estimated
  ## C is the Choleski factor of U
  dims <- dim(X)
  nbootstrap <- nrboot # number or bootstrap samples
  save_W <- array(NA,dim=c(nbootstrap,dims[2]^2))
  
  nbootobs <- dims[1] # number of observations in each bootstrap sample
  
  
  for (i in 1:nbootstrap) {
    if ((i%%10 == 0) & display) {cat("bootstrap run", i, "out of", nbootstrap, "\n")}
    set.seed(seed)
    pick <- sample(dims[1], nbootobs, replace = TRUE)
    #pick <- ceiling(runif(nbootobs)*dims[1])
    X_boot <- X[pick,]
    
    if (method=="fastICA"){
      icares             <- fastICA(X_boot, ncol(W), tol=1e-6, w.init=w.init )
      W.ica              <- t((icares$K) %*% (icares$W))
      W.scal.ica         <- rescaleVar(W_hat = W.ica, ut = t(u_chol))$Ws
      icares_MD          <- W.scal.ica %*% solve(C)
      W_boot             <- icares_MD
      B_boot             <- solve(W_boot) %*% myfrob(A.hat = solve(W_boot), A = Bt)
      
      
    }else if(method == "KernelICA"){
      require(KernelICA)
      kernelica_res      <- KernelICA::kernel_ica(X_boot, init = list(w.init), variant = "kgv", kernel = "gauss")
      W_kernelica        <- kernelica_res$W
      #W_scal_kernelica   <- rescaleVar(W_hat = W_kernelica, ut = t(u_chol))$Ws  # Assurez-vous que la fonction rescaleVar est définie ou remplacez-la par une fonction équivalente
      kernelica_MD       <- rescaleVar(W_hat = W_kernelica, ut = t(u_chol))$Ws %*% solve(C)
      W_boot             <- kernelica_MD
      B_boot             <- solve(W_boot) %*% myfrob(A.hat = solve(W_boot), A = Bt)  # Assurez-vous que la fonction myfrob est définie ou remplacez-la par une fonction équivalente
      

 

    }else if(method == "DCov"){
      DCov               <- steadyICA(X_boot, n.comp = ncol(W), w.init = w.init, 
                                      PIT = T, symmetric = F,maxit = 1000,verbose = F)
      DC_MD              <- solve(DCov$W) %*% solve(C)
      W_boot             <- DC_MD
      B_boot             <- solve(W_boot) %*% myfrob(A.hat = solve(W_boot), A = Bt)
      
      
    }else if(method == "CvM"){
      if (is.null(dd)) {
        dd <- indepTestSim(n, ncol(X), verbose = F)
      }
      lower <- rep(0, ncol(X) * (ncol(X) - 1)/2)
      upper <- rep(pi, ncol(X) * (ncol(X) - 1)/2)
      de_control <- list(itermax = 500, steptol = 100, 
                         trace = FALSE)
      ## First step of optimization with DEoptim
      
      de_res <- DEoptim(testlik, lower = lower, upper = upper, 
                        control = de_control, faklow = C, u = t(C %*% t(X)), dd = dd) 
      k <- ncol(X)
      iter2  <- 75
      ## Second step of optimization. Creating randomized starting angles around the optimized angles
      ## here maybe you insert the LHS stuff
      theta_rot <- matrix(rnorm(n = (k * (k - 1)/2) * iter2, mean = de_res$optim$bestmem, 
                                sd = 0.3), (k * (k - 1)/2), iter2)
      theta_rot <- cbind(theta_rot, de_res$optim$bestmem)
      # Start vectors for iterative optimization approach
      startvec_list <- as.list(as.data.frame(theta_rot))
      erg_list <- lapply(X = startvec_list, FUN = optim, fn = testlik, 
                         gr = NULL, faklow = C, u = t(C %*% t(X)), dd = dd, 
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
      Acvm               <- rotmat(par_o, C)
      B_boot             <- Acvm %*% myfrob(A.hat = Acvm,A = Bt)
      
    }
    W_boot <- solve(B_boot)
    
    
    # need to row-permute matrix for fitting it best to W
    W_perm <- best_permutation(W,W_boot)
    
    
    save_W[i,] <- W_perm
    
  }
  B.boot <- vector("list",length = nrboot)
  test.t <- array(data = NA,dim = c(ncol(C),ncol(C),5), dimnames = 
                    list(paste0("row",seq(1,nrow(C))),
                         paste0("col",seq(1,nrow(C))),
                         c("a","a","a","a","a")))
  for (j in 1:nrboot){
    B.boot[[j]] <- solve(matrix(data = save_W[j,],nrow = ncol(C),ncol = ncol(C),byrow = F))
  }
  cnt <- 0
  for (alpha in c(0.01, 0.05, 0.1, 0.25, 0.5)){
    cnt <- cnt +1
    dimnames(test.t)[[3]][cnt] <- paste(text='alpha =', alpha,sep = "")
    for (r in 1:nrow(C)){
      for (c in 1:ncol(C)){
        ## create columns with bootstrap estimates of the pmt
        eval(parse(text=paste("b",r,c,".boot<-do.call(rbind,lapply(B.boot,'[',",r,",",c,"))", sep="")))
        ## create columns with bootstrap estimates of b_boot - b_true
        eval(parse(text=paste("XX",r,c," <- b",r,c,".boot-Atrue[",r,",",c,"]",sep="")))
        a <- quantile(x = eval(parse(text=paste("XX",r,c,sep = ""))),probs = c(alpha/2, (1-alpha/2)))
        if (0 < min(a) | 0 > max(a)){
          test.t[r,c,cnt] <- TRUE
        }else test.t[r,c,cnt] <- FALSE
        
      }
    }
  }
  ls   <- list(w.pmts = save_W, B_boot = B.boot, test = test.t)
  return(ls)
}


myfrob <- function(A.hat, A){
  ##sign-column permutation that minimizes the frobenius norm of A.hat - A_0
  a   <- gtools::permutations(ncol(A),ncol(A))
  b   <- gtools::permutations(2,ncol(A),repeats.allowed = T)
  b[b==2] <- -1
  PR   <- vector("list", length = nrow(a)*nrow(b))
  frob <- double(length = nrow(a)*nrow(b))
  pr  <- diag(ncol(A))
  cnt <- 0
  for (j in 1:nrow(a)){
    for (i in 1:nrow(b)){
      cnt <- cnt+1
      PR[[cnt]]  <- pr[,a[j,]] * b[i,]
      A.tilde    <- A.hat %*% PR[[cnt]]
      frob[cnt]  <- frobenius.norm(A.tilde-A)
    }
  }
  return(PR[[which.min(frob)]])
}


nzdiagscore <- function( W ) {
  
  res <- sum(1/diag(abs(W)))
  res
  
}

# computes all permutations of the columns of a matrix
# returns a list of matrices
permutation <- function(mat) {
  permutations <- permute_vector(1:ncol(mat))
  
  res_list <- list()
  for (i in 1:ncol(permutations)) {
    res_list[[i]] <- mat[, permutations[, i, drop = TRUE], drop = FALSE]
  }
  return(res_list)
}

permute_vector <- function(x) {
  # calculate all permutations of a vector x
  # naiive recursive implementation.
  # @param x a vector
  # @return matrix; each column is one unqiue permutation of x
  n <- length(x)
  num_permutation <- factorial(n)
  
  if (n == 1) { # base case
    return(as.matrix(x))
  } else {
    num_sub_permutation <- factorial(n - 1)
    res <- matrix(0.0, nrow = n, ncol = num_permutation)
    for (i in 1:n) {
      # swap first entry in x with i-th one
      y <- x
      tmp <- y[1]
      y[1] <- y[i]
      y[i] <- tmp
      
      # calculate all permutations of y[2:n]
      col_idx <- (i-1) * num_sub_permutation + seq(1:num_sub_permutation)
      res[1, col_idx] <- y[1]
      res[2:n, col_idx] <- permute_vector(y[2:n])
    }
    return(res)
  }
}

## Pseudo log-likelihood functions
## 2D case ####
## Super gaussian
pseudo.log.lik2D.super        <- function(c) {
  c11 <- c[1]
  c21 <- c[2]
  c12 <- c[3]
  c22 <- c[4]
  
  c.1y <- c11*u_chol[[i]][,1] + c21*u_chol[[i]][,2]
  c.2y <- c12*u_chol[[i]][,1] + c22*u_chol[[i]][,2]
  
  g <- rep(NA,1)
  g[1] <- -1*(n*log(1/sqrt(2*pi))-sum((c.1y)^2/2) + 
                n*log(gamma((5+1)/2)/(sqrt(5*pi)*gamma(5/2)))-((5+1)/2)*sum(log(1+(c.2y)^2/5)))  
  g
}

rescaleVar <- function(W_hat, ut){
  ## W_hat is the trasposded version of the unmixing matrix estimated by the method, E = W_hat * U
  ## ut is (kxn) matrix of structural shocks 
  AA_hat    <- solve(W_hat) # AA_hat is the mixing matrix
  eres_hat  <- t(W_hat %*% ut) # estimated ICAs n*k
  # we rescale the mixing matrix so that the independent components (i.e. structural shocks) have variance equal to one
  Ascaled_hat    <- AA_hat
  eresscaled_hat <- eres_hat
  for (i in 1:nrow(AA_hat)) {
    eresscaled_hat[, i] <- eres_hat[, i] / sd(eres_hat[, i])
    Ascaled_hat[, i]    <- AA_hat[, i] * sd(eres_hat[, i])
  }
  Wscaled_hat <- solve(Ascaled_hat)
  ## Wscaled is the unmixing matrix eps = Wscaled * u_t
  ## We want the mixing matrix A_hat s.t u_t = A_scaled*eps_t
  return(list(Ws = Wscaled_hat, As = Ascaled_hat))
}

stdEpsrnd <- function(type, n, dof, p = 2){
  ## type is the kind of distribution
  ## n is the sample size
  ## dof degrees of freedom
  ## seed 
  #set.seed(seed)
  if (type == "tstud"){
    eps_sim <- rt(n,dof)
  } else if(type == "laplace"){
    eps_sim <- rlaplace(n,m = 0, s = 1) 
    ## m,s are specific to resemble a nicer distrib. See different
  } else if(type == "logistic"){
    eps_sim <- rlogis(n, location = 0, scale = (3.14 ^ (0.5)))
  } else if(type == "general"){
    eps_sim <- rpgnorm(n = n,mean = 0, p = p,  method="pgenpolar")
  } else {
    print("Not valid distribution")
  }
  ## Standardize the random variable
  eps <- (eps_sim - mean(eps_sim)) / (sd(eps_sim))
  return(eps)
}

# function which rotates the B matrix and calculates the indpendence of the structural errors
testlik <- function(theta, faklow, u, dd) {
  
  temp_l <- rotmat(theta, faklow)
  
  ser_low <- tcrossprod(u, solve(temp_l))
  ddtest <- copula::indepTest(ser_low, dd)
  ddtest$global.statistic # * 10000000
  
}


# function to compute the optimal B matrix
rotmat <- function(pv, faklow) {
  
  ts_dim <- (1 + sqrt(1 + 8 * length(pv))) / 2
  combns <- combn(ts_dim, 2, simplify = FALSE)
  rotmat <- diag(ts_dim)
  
  for (i in seq_along(pv)) {
    
    tmp <- diag(ts_dim)
    tmp[combns[[i]][1], combns[[i]][1]] <- cos(pv[i])
    tmp[combns[[i]][2], combns[[i]][2]] <- cos(pv[i])
    tmp[combns[[i]][1], combns[[i]][2]] <- - sin(pv[i])
    tmp[combns[[i]][2], combns[[i]][1]] <- sin(pv[i])
    rotmat <- rotmat %*% tmp
    
  }
  
  tcrossprod(faklow, rotmat)
  
}


wild.boot.gp <- function(X, B, C = C, w.init = diag(nrow(W)), nrboot, method, display = T, quantiles, ordering) {
  ## Bootstrap on mixing matrix pmts
  ## X are the orthogonalized reduced-form residuals
  ## W is the unmixing matrix estimated
  ## C is the Choleski factor of U
  dims <- dim(X)
  nbootstrap <- nrboot # number or bootstrap samples
  save_W <- array(NA,dim=c(nbootstrap,dims[2]^2))
  save_B <- array(NA,dim=c(nbootstrap,dims[2]^2))
  
  nbootobs <- dims[1] # number of observations in each bootstrap sample
  
  
  for (i in 1:nbootstrap) {
    if ((i%%10 == 0) & display) {cat("bootstrap run", i, "out of", nbootstrap, "\n")}
    pick <- sample(dims[1], nbootobs, replace = TRUE)
    
    X_boot  <- X*rnorm(n = nrow(X))
    
    if (method=="fastICA"){
      icares <- fastICA(X_boot, ncol(B), tol=1e-6, w.init=w.init )
      W.ica              <- t((icares$K) %*% (icares$W))
      W.scal.ica         <- rescaleVar(W_hat = W.ica, ut = t(X_boot))$Ws
      icares_MD          <- W.scal.ica %*% solve(C)
      W_boot             <- icares_MD
      if(ordering == "maxfinder"){
        B_boot              <- maxfinder(A = solve(W_boot))$A.id
      } else if(ordering == "frob") B_boot         <- solve(W_boot) %*% myfrob(A.hat = solve(W_boot),A = B)
      B_boot[,1]*1/B[1,1] ->B_boot[,1]
      B_boot[,2]*1/B[2,2] ->B_boot[,2]
      B_boot[,3]*1/B[3,3] ->B_boot[,3]
    }else if(method == "DCov"){
      DCov               <- steadyICA(X_boot, n.comp = ncol(B), w.init = w.init, 
                                      PIT = F, symmetric = F,maxit = 1000,verbose = F)
      #DCov               <- dcovICA(Z = X_boot,theta.0 = w.init)
      DC_MD              <- solve(DCov$W) %*% solve(C)
      W_boot             <- DC_MD
      if(ordering == "max"){
        B_boot              <- maxfinder(A = solve(W_boot))$A.id
      } else if(ordering == "frob") B_boot         <- solve(W_boot) %*% myfrob(A.hat = solve(W_boot),A = B)
      B_boot[,1]*1/B[1,1] ->B_boot[,1]
      B_boot[,2]*1/B[2,2] ->B_boot[,2]
      B_boot[,3]*1/B[3,3] ->B_boot[,3]
      
    }else if(method == "CvM"){
      if (is.null(dd)) {
        dd <- indepTestSim(n, ncol(X), verbose = F)
      }
      lower <- rep(0, ncol(X) * (ncol(X) - 1)/2)
      upper <- rep(pi, ncol(X) * (ncol(X) - 1)/2)
      de_control <- list(itermax = 500, steptol = 100, 
                         trace = FALSE)
      ## First step of optimization with DEoptim
      
      de_res <- DEoptim(testlik, lower = lower, upper = upper, 
                        control = de_control, faklow = C, u = X_boot, dd = dd) #U or u_chol
      k <- ncol(X)
      iter2  <- 75
      ## Second step of optimization. Creating randomized starting angles around the optimized angles
      ## here maybe you insert the LHS stuff
      theta_rot <- matrix(rnorm(n = (k * (k - 1)/2) * iter2, mean = de_res$optim$bestmem, 
                                sd = 0.3), (k * (k - 1)/2), iter2)
      theta_rot <- cbind(theta_rot, de_res$optim$bestmem)
      # Start vectors for iterative optimization approach
      startvec_list <- as.list(as.data.frame(theta_rot))
      erg_list <- lapply(X = startvec_list, FUN = optim, fn = testlik, 
                         gr = NULL, faklow = C, u = X_boot, dd = dd, 
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
      Acvm               <- rotmat(par_o, C)
      if(ordering == "max"){
        B_boot              <- maxfinder(A = Acvm)$A.id
      } else if(ordering == "frob") B_boot         <- Acvm %*% myfrob(A.hat = Acvm,A = B)
      B_boot[,1]*1/B[1,1] ->B_boot[,1]
      B_boot[,2]*1/B[2,2] ->B_boot[,2]
      B_boot[,3]*1/B[3,3] ->B_boot[,3]
    }
    W_boot <- solve(B_boot)
    
  
    
    
    save_B[i,] <- as.vector(B_boot)
    
    
  }
  B.boot <- vector("list",length = nrboot)
  test.t <- array(data = NA,dim = c(ncol(C),ncol(C),length(quantiles)),
                  dimnames = list(paste0(seq(1,ncol(C),1)),paste0(seq(1,ncol(C),1)), paste0(quantiles)))
  
  for (j in 1:nbootstrap){
    B.boot[[j]] <- (matrix(data = save_B[j,],nrow = ncol(C),ncol = ncol(C),byrow = F))
  }
  cnt <- 0
  for (alpha in quantiles){
    
    cnt <- cnt +1
    
    for (r in 1:nrow(C)){
      for (c in 1:ncol(C)){
        ## create columns with bootstrap estimates of the pmt
        eval(parse(text=paste("b",r,c,".boot<-do.call(rbind,lapply(B.boot,'[',",r,",",c,"))", sep="")))
        ## create columns with bootstrap estimates of b_boot - b_true, b_true = 0
        eval(parse(text=paste("XX",r,c," <- b",r,c,".boot",sep="")))
        a <- quantile(x = eval(parse(text=paste("XX",r,c,sep = ""))),probs = c(alpha/2, (1-alpha/2)))
        if (0 < min(a) | 0 > max(a)){
          test.t[r,c,cnt] <- TRUE ##  reject = 0
        }else test.t[r,c,cnt] <- FALSE 
      }
    }
  }
  ls   <- list(w.pmts = save_W, B_boot = B.boot, test = test.t)
  return(ls)
}


## Sub gaussian
pseudo.log.lik2D.sub          <- function(c) {
  c11 <- c[1]
  c21 <- c[2]
  c12 <- c[3]
  c22 <- c[4]
  
  c.1y <- c11*u_chol[[i]][,1] + c21*u_chol[[i]][,2]
  c.2y <- c12*u_chol[[i]][,1] + c22*u_chol[[i]][,2]
  
  g <- rep(NA,1)
  g[1] <- -1*(n*log(1/sqrt(2*pi)) + sum(pi*(c.1y)^2+log(cosh((pi/2)*c.1y))) + 
                n*log(1/sqrt(2*pi)) + sum(pi*(c.2y)^2+log(cosh((pi/2)*c.2y))))
  g
}

## Contraint
orth.constr2D                 <- function(c){
  c11 <- c[1]
  c21 <- c[2]
  c12 <- c[3]
  c22 <- c[4]
  g <- rep(NA,3)
  g[1] <- c11*c12 + c21*c22 # c.1*c.2
  g[2] <- c11^2 + c21^2 - 1 # c.1*c.1
  g[3] <- c12^2 + c22^2 - 1 # c.2*c.2
  g
}

## 3D ####
## Super gaussian
pseudo.log.lik3D.super        <- function(c) {
  c11 <- c[1]
  c21 <- c[2]
  c31 <- c[3]
  c12 <- c[4]
  c22 <- c[5]
  c32 <- c[6]
  c13 <- c[7]
  c23 <- c[8]
  c33 <- c[9]
  
  c.1y <- c11*u_chol[[i]][,1] + c21*u_chol[[i]][,2] + c31*u_chol[[i]][,3]
  c.2y <- c12*u_chol[[i]][,1] + c22*u_chol[[i]][,2] + c32*u_chol[[i]][,3]
  c.3y <- c13*u_chol[[i]][,1] + c23*u_chol[[i]][,2] + c33*u_chol[[i]][,3]
  
  
  g <- rep(NA,1)
  g[1] <- -1*(n*log(1/sqrt(2*pi))-sum((c.1y)^2/2) + 
                n*log(gamma((5+1)/2)/(sqrt(5*pi)*gamma(5/2)))-((5+1)/2)*sum(log(1+(c.2y)^2/5)) + 
                n*log(gamma((12+1)/2)/(sqrt(12*pi)*gamma(12/2)))-((12+1)/2)*sum(log(1+(c.3y)^2/12)))
  g
}

## Sub gaussian
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
  
  c.1y <- c11*u_chol[[i]][,1] + c21*u_chol[[i]][,2] + c31*u_chol[[i]][,3]
  c.2y <- c12*u_chol[[i]][,1] + c22*u_chol[[i]][,2] + c32*u_chol[[i]][,3]
  c.3y <- c13*u_chol[[i]][,1] + c23*u_chol[[i]][,2] + c33*u_chol[[i]][,3]
  
  
  g <- rep(NA,1)
  g[1] <- -1*(n*log(1/sqrt(2*pi)) + sum(pi*(c.1y)^2+log(cosh((pi/2)*c.1y))) + 
                n*log(1/sqrt(2*pi)) + sum(pi*(c.2y)^2+log(cosh((pi/2)*c.2y))) + 
                n*log(1/sqrt(2*pi)) + sum(pi*(c.3y)^2+log(cosh((pi/2)*c.3y))))
  g
}
## Constraint
orth.constr3D    <- function(c){
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



