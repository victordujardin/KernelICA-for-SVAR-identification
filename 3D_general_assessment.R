## General assessment 
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

n              <- 200 # sample size
choose_cluster <- 15

source(here("confg_file.R"),local = T)

#### Codes to be updated and experiment setting ##
source(here("functions_to_load.R"),local = T)

## Create the distribution of the empirical copula under the null of inpendence(for CvM)
dd <- NULL
if (is.null(dd)) {
  dd <- indepTestSim(n, ncol(Btrue3D), verbose = T)
}
## Number of initialization Hypercybe sampling
if (n_lhs == 1) {
  lhs <- matrix(data = c(0,0,0),nrow = 1)
}else{
  set.seed(55)
  lhs            <- 2*pi*improvedLHS(n = n_lhs, k = 3)
}

for(z in 1:length(SEQ)){
  print(paste0("experiment ", z, " - first MC ", Sys.time()))
  print(paste0("sample size = ",n," initialization =  ", n_lhs, " MC = ",mc))
  nr_cluster <- choose_cluster
  cl <- makeCluster(nr_cluster)
  registerDoParallel(cl)
  ## Results reproducible 
  set.seed(123+z)
  ciao <- foreach(j = 1:mc,
                  .packages = c('pgnorm', 'normtest', 'steadyICA', "ProDenICA",
                                'fastICA', 'JADE', 'tidyverse', 'magrittr', 'gtools', 
                                'DEoptim', "matrixcalc", "copula")) %dorng% {
                                  
                                  Sys.sleep(1)  
                                  
                                  
                                  ## Run methods and record best initialization matrix
                                  eps1             <- stdEpsrnd(type = 'general', n = n, p = SEQ[z])
                                  eps2             <- stdEpsrnd(type = 'general', n = n, p = SEQ[z])
                                  eps3             <- stdEpsrnd(type = 'general', n = n, p = SEQ[z])
                                  pvalJB1[j]       <- as.numeric(ajb.norm.test(x = eps1, nrepl = n)$p.value)
                                  pvalJB2[j]       <- as.numeric(ajb.norm.test(x = eps2, nrepl = n)$p.value)
                                  pvalJB3[j]       <- as.numeric(ajb.norm.test(x = eps3, nrepl = n)$p.value)
                                  E                <- cbind(eps1, eps2, eps3)             # n*k matrix of structural shocks (iid)
                                  Atrue            <- mixmat(p = ncol(lhs))
                                  U                <- Atrue %*% t(E)
                                  
                                  # Sphering data
                                  sigg                <- cov(t(U))
                                  C                   <- t(chol(sigg))
                                  u_chol              <- t(solve(C) %*% U)
                                  
                                  # Find best Initial condition for W (Dcov & fastICA)
                                  for (kk in 1:nrow(lhs)){
                                    if((kk %% 20)==0) print(paste0("initialization #",kk))
                                    ## initialization
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
                                    
                                    winit[[kk]]        <- winitx %*% winity %*% winitz
                                    
                                    ## Dist Cov 
                                    dc_init[[kk]]      <- steadyICA(X = u_chol, n.comp = ncol(lhs), w.init = winit[[kk]], 
                                                                    PIT = F, symmetric = F,maxit = 1000,verbose = F)
                                    
                                    
                                    ## Given the nature of DCOV estimation, the matrix C%*%W is the estimate of the mixing matrix
                                    ## So, in order to obtain the estimate of the umixing matrix we have to calculate the following:
                                    ## The one and only!
                                    Wscaled_dc[[kk]]   <- solve(dc_init[[kk]]$W) %*% solve(C)
                                    ## fastICA  
                                    ica_init[[kk]]     <- fastICA(u_chol, n.comp = ncol(lhs), tol = 1e-14, w.init = winit[[kk]], 
                                                                  maxit = 1000,verbose = FALSE)
                                    Wica[[kk]]         <- t((ica_init[[kk]]$K) %*% (ica_init[[kk]]$W))
                                    Wscal_ica[[kk]]    <- rescaleVar(W_hat = Wica[[kk]], ut = t(u_chol))$Ws
                                    Wscaled_ica[[kk]]  <- Wscal_ica[[kk]]%*%solve(C)
                                    
                                    
                                    ##KernelICA
                                    
                                    # Use kernel_ica method
                                    kica_init[[kk]] <- KernelICA::kernel_ica(u_chol, variant = "kgv", kernel = "gauss", eps = 1e-14, init = list(winit[[kk]]))
                                    
                                    
                                    # Given the nature of DCOV estimation, the matrix C%*%W is the estimate of the mixing matrix
                                    # So, in order to obtain the estimate of the unmixing matrix we have to calculate the following:
                                    
                                    
                                    Wscaled_kica[[kk]] <- rescaleVar(W_hat = kica_init[[kk]]$W, ut = t(u_chol))$Ws %*% solve(C)
                                    
                                    ## Performance measures
                                    
                                    ## The more coherent measure for the raw unmixing matrix given the underyining model
                                    ## and the specification of the package function
                                    ## The one and only measure. Taken from svars. Proved and confirmed NT NT
                                    FR_dc_nt[[kk]]     <- MD(W.hat = (Wscaled_dc[[kk]]), A = (Atrue))
                                    FR_ica_nt[[kk]]    <- MD(W.hat = (Wscaled_ica[[kk]]), A = (Atrue))
                                    FR_kica_nt[[kk]]    <- MD(W.hat = (Wscaled_kica[[kk]]), A = (Atrue))
                                  }
                                  ## Select Best initial conditions (in 2D these are angles for steadyICA)
                                  best.W0.MD.dc     <- winit[[which.min(x = FR_dc_nt)]]
                                  best.W0.MD.ica    <- winit[[which.min(x = FR_ica_nt)]]
                                  best.W0.MD.kica    <- winit[[which.min(x = FR_kica_nt)]]
                                  
                                  
                                  # Select the best first W matrix
                                  DC_MD[[j]]        <- Wscaled_dc[[which.min(x = FR_dc_nt)]]
                                  icares_MD[[j]]    <- Wscaled_ica[[which.min(x = FR_ica_nt)]]
                                  kicares_MD[[j]]    <- Wscaled_kica[[which.min(x = FR_kica_nt)]]
                                  
                                  # 
                                  # #### Cvm method: Adapted from svars R package
                                  # if (is.null(dd)) {
                                  #   dd <- indepTestSim(n, ncol(E), verbose = F)
                                  # }
                                  # lower <- rep(0, ncol(E) * (ncol(E) - 1)/2)
                                  # upper <- rep(pi, ncol(E) * (ncol(E) - 1)/2)
                                  # de_control <- list(itermax = 500, steptol = 100, 
                                  #                    trace = FALSE)
                                  # ## First step of optimization with DEoptim
                                  # 
                                  # de_res <- DEoptim(testlik, lower = lower, upper = upper, 
                                  #                   control = de_control, faklow = C, u = t(U), dd = dd) #U or u_chol
                                  # k <- ncol(E)
                                  # iter2  <- 75
                                  # ## Second step of optimization. Creating randomized starting angles around the optimized angles
                                  # ## here maybe you insert the LHS stuff
                                  # theta_rot <- matrix(rnorm(n = (k * (k - 1)/2) * iter2, mean = de_res$optim$bestmem, 
                                  #                           sd = 0.3), (k * (k - 1)/2), iter2)
                                  # theta_rot <- cbind(theta_rot, de_res$optim$bestmem)
                                  # # Start vectors for iterative optimization approach
                                  # startvec_list <- as.list(as.data.frame(theta_rot))
                                  # erg_list <- lapply(X = startvec_list, FUN = optim, fn = testlik, 
                                  #                    gr = NULL, faklow = C, u = t(U), dd = dd, 
                                  #                    method = ifelse(k == 2, "Brent", "Nelder-Mead"), 
                                  #                    lower = ifelse(k == 2, -.Machine$double.xmax/2, -Inf), 
                                  #                    upper = ifelse(k == 2, .Machine$double.xmax/2, Inf), 
                                  #                    control = list(maxit = 1000), 
                                  #                    hessian = FALSE)
                                  # 
                                  # # Print log-likelihood values from local maxima
                                  # logliks <- sapply(erg_list, "[[", "value")
                                  # if (min(logliks) < de_res$optim$bestval) {
                                  #   params <- sapply(erg_list, "[[", "par", simplify = FALSE)
                                  #   par_o <- params[[which.min(logliks)]]
                                  #   logs <- min(logliks)
                                  #   inc <- 1
                                  # } else {
                                  #   par_o <- de_res$optim$bestmem
                                  #   logs <- de_res$optim$bestval
                                  #   inc <- 0
                                  # }
                                  # A.cvm[[j]] <- rotmat(par_o, C)   ## Estimated Mixing matrix
                                  
                                  
                                  ## Performances on the ICA given best W0 and same dgp
                                  perfm_dc[j]        <- MD(W.hat = DC_MD[[j]],A = Atrue)
                                  perfm_ica[[j]]     <- MD(W.hat = icares_MD[[j]],A = Atrue)
                                  perfm_kica[[j]]     <- MD(W.hat = kicares_MD[[j]],A = Atrue)
                                  # perfm_cvm[[j]]     <- MD(W.hat = solve(A.cvm[[j]]),A = Atrue) #solve because we compute first mixing
                                  # 
                                  A.dc[[j]]          <- solve(DC_MD[[j]])
                                  A.ica[[j]]         <- solve(icares_MD[[j]])
                                  
                                  
                                  
                                  
                                  mc_results <- list(
                                    # A.dc     = A.dc[[j]], 
                                    # A.ica    = A.ica[[j]],
                                    # A.cvm    = A.cvm[[j]],
                                    # C        = C,
                                    # U        = U,
                                    # Atrue    = Atrue,
                                    perfm_dc  = perfm_dc[j],  
                                    perfm_ica = perfm_ica[[j]],
                                    perfm_kica = perfm_kica[[j]]
                                    # ,
                                    # perfm_cvm = perfm_cvm[[j]]
                                    # ,
                                    # FR_ica_nt  = FR_ica_nt,
                                    # pvalJB1   = pvalJB1[j],
                                    # pvalJB2   = pvalJB2[j],
                                    # pvalJB3   = pvalJB3[j],
                                    # Wscaled_ica = Wscaled_ica,
                                    # Wscaled_dc  = Wscaled_dc,
                                    # icares_MD = icares_MD[[j]],
                                    # DC_MD     = DC_MD[[j]]
                                    )
                                }
  
  stopCluster(cl)
  
  pippo[[z]]           <- list(mc_results = ciao)
  
  
}

inference_results <- pippo

names(inference_results) <- paste0("p_",round(SEQ,2))

if (!dir.exists(here("results"))) {
  dir.create(here("results"))
}
save(inference_results,file = here("results", paste0("2D_other_ica_general_n",n,".RData")))
rm(pippo)





inference_results






































data_list <- list()

for (name in names(inference_results)) {
  p_data <- inference_results[[name]]$mc_results
  perf_metrics <- c(
    perfm_dc = unlist(lapply(p_data, function(x) x$perfm_dc)),
    perfm_ica = unlist(lapply(p_data, function(x) x$perfm_ica)),
    perfm_kica = unlist(lapply(p_data, function(x) x$perfm_kica))
    # ,
    # perfm_cvm = unlist(lapply(p_data, function(x) x$perfm_cvm))
  )
  data_list[[name]] <- perf_metrics
}

data_list



data_agg <- list()
# you can change this value to any other number as needed

for (p in names(data_list)) {
  # Create a named vector to store aggregated values
  aggregated_values <- c(
    perfm_dc = mean(data_list[[p]][paste0('perfm_dc', 1:mc)]),
    perfm_ica = mean(data_list[[p]][paste0('perfm_ica', 1:mc)]),
    perfm_kica = mean(data_list[[p]][paste0('perfm_kica', 1:mc)])
    # ,
    # perfm_cvm = mean(data_list[[p]][paste0('perfm_cvm', 1:mc)])
  )
  # Assign the named vector to data_agg
  data_agg[[p]] <- aggregated_values
}

data_agg












# Load necessary library
library(tidyverse)

# Convert data_agg to a data frame
data_df <- bind_rows(lapply(names(data_agg), function(p) {
  data.frame(
    p = gsub("p_", "", p),  # Extract p value from the name
    method = names(data_agg[[p]]),  # Extract method names from the names of x
    perfm = unlist(data_agg[[p]])  # Extract performance values
  )
}))

# Convert method to a factor for proper coloring
data_df$method <- factor(data_df$method,
                         levels = c("perfm_dc", "perfm_ica", "perfm_kica"),
                         labels = c("Distance-covariance", "FastICA", "KernelICA"))
# , "perfm_cvm"))

# Convert the 'p' column to a factor with ordered levels
data_df$p_factor <- factor(data_df$p, levels = sort(unique(data_df$p)))



# Convertir le facteur en numérique
vecteur_numerique <- as.numeric(as.character(data_df$p_factor))

# Trier le vecteur numérique
vecteur_numerique_ordonne <- sort(vecteur_numerique)

# Convertir le vecteur numérique trié en facteur avec les niveaux correctement ordonnés
facteur_ordonne <- factor(vecteur_numerique_ordonne, levels = sort(unique(vecteur_numerique_ordonne)))


# Plot the graph with lines connecting the points
ggplot(data_df, aes(x = facteur_ordonne, y = perfm, color = method)) +
  geom_line(aes(group = method)) +  # Group by method to connect lines
  geom_point() +
  labs( #title = paste("Performance by Method and value of p", "init : ", n_lhs, "sample size :", n, " mc : ",mc),
    x = "p",
    y = "Performance",
    color = "Method") +
  theme_minimal()


# 
# # Assuming inference_results is already loaded and available
# 
# # ... [Your initial package loading and setup code remains unchanged]
# 
# # Create an empty data frame to store the results
# result_df <- data.frame(
#   estimator = character(),
#   minimum_distance = numeric(),
#   scenario = numeric(),
#   dimension = numeric(),
#   mc = numeric(),
#   p_shape = numeric(),
#   stringsAsFactors = FALSE
# )
# 
# # Mapping the performance metrics to their respective estimator names
# estimator_map <- list(
#   perfm_dc = "DCov",
#   perfm_ica = "fastICA",
#   perfm_kica = "KernelICA",
#   perfm_cvm = "CvM"
# )
# 
# # Iterate over each name in inference_results
# for (name in names(inference_results)) {
#   p_data <- inference_results[[name]]$mc_results
#   
#   # Extract performance metrics for each estimator
#   for (estimator_name in names(estimator_map)) {
#     distances <- unlist(lapply(p_data, function(x) x[[estimator_name]]))
#     
#     tmp_df <- data.frame(
#       estimator = estimator_map[[estimator_name]],
#       minimum_distance = distances,
#       scenario = length(distances),  # Assuming scenario is the sample size, you might need to adjust this
#       dimension = 2,  # Since k=2
#       mc = 1:length(distances),  # Assuming mc represents the iteration of the monte carlo
#       p_shape = as.numeric(gsub("p_", "", name))  # Extracting p value from the name
#     )
#     
#     # Append to the result data frame
#     result_df <- rbind(result_df, tmp_df)
#   }
# # }
# 
# # Write to CSV
# write.csv(result_df, "desired_output.csv", row.names = FALSE)
# 
# 
# 
# 
# 









