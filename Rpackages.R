## Packages to load
packages_needed <- c("lhs","gtools","copula","pgnorm","normtest",
                         "DEoptim","steadyICA","fastICA","JADE","tidyverse",
                         "matrixcalc","doParallel","doSNOW","doRNG","moments",
                         "alabama","ProDenICA","gtools","pbapply","kableExtra","gridExtra","Hmisc","vars")

installed_packages  <- which(sapply(X = packages_needed,FUN = function(x) require(x,character.only = T))==F)
packages_to_install <- sapply(X = installed_packages, FUN = function(x) install.packages(x,lib =  .libPaths(),
                                                                                         repos = getOption("repos"),
                                                                                         contriburl = contrib.url(repos, type)))
# lhs
# gtools
# copula # required for indepTest: Test Independence of Continuous Random Variables via Empirical Copula
# pgnorm # The p-Generalized Normal Distribution
# normtest # To perform the augmented Jarque-Bera test
# DEoptim # Global Optimization by Differential Evolution, adapted from R package svars
# steadyICA # Distance Covariance estimator from Matteson Tsay
# fastICA # FastICA Algorithms to Perform ICA and Projection Pursuit
# JADE #required for function MD that computes the Minimum Distance index MD to evaluate the performance of an ICA algorithm.
# tidyverse # data handling
# matrixcalc # A collection of functions to support matrix calculations for probability, econometric and numerical analysis
# doParallel #Foreach Parallel Adaptor for the 'parallel' Package
# doSNOW # Performing simulation in parallel
# doRNG #Provides functions to perform reproducible parallel foreach loop
# moments ## kurtosis
# alabama ## Constrained optimization
# ProDenICA ## Use mixmat function to extrapolate random mixing matrix
# gtools ## Permutations
# pbapply ## Add progress bar in some functions
# gridExtra

