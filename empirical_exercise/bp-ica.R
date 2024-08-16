## Blanchard - Perotti exercise with ICA identification

## General assessment 
## Load packege here for managing workflow


library(MASS) 
library(reshape2) 
library(reshape) 

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

## Load auxiliary functions
source(here("functions_to_load.R"),local = T)

## Blanchard Perotti Database
BP02 <- read.csv(here("final_databases","blanchardperottius2002.csv"), row.names=1)
bp02 <- BP02[57:208,]


### Variables preparation ####
## Order of variables
yy <- data.frame(GG =  bp02$GG,
                 TT =  bp02$TT,
                 XX =  bp02$XX)
## Create time series object
tsy  <- ts(data = yy,start = c(1960,1),frequency = 4,names = c("GG","TT","XX"))

## Create Exogenous Variables matrix
source(here("empirical_exercise","exog-vars.R"),local = T)

#### Estimation VAR ####
## VAR a la Blanchard Perotti
var.interactions    <- VAR(y = tsy, p = 4, type = "both",exogen = EXOG)

## VAR for IRFs
var.irf             <- VAR(y = tsy, p = 4, type = "both",exogen = EXOG2)

var             <- VAR(y = tsy, p = 4, type = "both")





JB.test <- normality.test(x = var.interactions,
                          multivariate.only = F) ## Only Spending Reduced-form residuals are normal, ICA can be applied
JB.test <- normality.test(x = var.irf,
                          multivariate.only = F) ## Only Spending Reduced-form residuals are normal, ICA can be applied

#### ICA IDENTIFICATION ####
u.ica  <- residuals(var.interactions)
#u.ica  <- residuals(var.irf)
kurtosis(u.ica[,1])
kurtosis(u.ica[,2])
kurtosis(u.ica[,3])

ajb.norm.test(u.ica)

ajb.norm.test(u.ica[,1])
ajb.norm.test(u.ica[,2])
ajb.norm.test(u.ica[,3])

# plot(density(u.ica[,1]))
# plot(density(u.ica[,2]))
# plot(density(u.ica[,3]))


### Orthogonalize Reduced-form residuals
sigg   <- cov((u.ica)) ## Variance-Covariance Matrix
C      <- t(chol(sigg)) ## Lower-triangular Choleski factor
u_chol <- solve(C) %*% t(u.ica) ### Orthogonalized reduced form residuals (Preliminary Rotation)
tu_chol<- t(solve(C) %*% t(u.ica))

## Call function to maximize
source(here("empirical_exercise","/maxlik3D-pml.R"),local = T)


## Select Best initial conditions
n_lhs  <- 100 ## number of initializations

## Run initialization for FastICA, DCov and PML + CVM estimation
source(here("empirical_exercise","find-init-and-cvm.R"),local = T)

## Identified impact matrix under CvM
Aidcvm <- maxfinder(A = Acvm)$A.id

## Find the maximum entry per each column
for(i in 1:n_lhs){
  c.ica[[i]] <- A.ID.ica[[i]] %>% apply(MARGIN = 2, function(x) which.max(abs(x)))
  c.kica[[i]] <- A.ID.kica[[i]] %>% apply(MARGIN = 2, function(x) which.max(abs(x)))
  c.dc[[i]]  <- A.ID.dc[[i]]  %>% apply(MARGIN = 2, function(x) which.max(abs(x)))
  c.pml[[i]] <- A.ID.pml[[i]] %>% apply(MARGIN = 2, function(x) which.max(abs(x)))
}

entries <- bind_rows(
unlist(lapply(c.ica, function(x) paste(as.character(x),collapse = "-"))) %>% as_tibble() %>% 
  group_by(value) %>% summarise(tot_type = n()) %>% 
  mutate(method = "fastica"),
unlist(lapply(c.kica, function(x) paste(as.character(x),collapse = "-"))) %>% as_tibble() %>% 
  group_by(value) %>% summarise(tot_type = n()) %>% 
  mutate(method = "KernelICA"),
unlist(lapply(c.dc, function(x) paste(as.character(x),collapse = "-"))) %>% as_tibble() %>% 
  group_by(value) %>% summarise(tot_type = n()) %>% 
  mutate(method = "DCov"),
unlist(lapply(c.pml, function(x) paste(as.character(x),collapse = "-"))) %>% as_tibble() %>% 
  group_by(value) %>% summarise(tot_type = n()) %>% 
  mutate(method = "pml"))



## Which Maxfinder scheme is more frequent given the n_lhs initializations?
entries %>% 
  ggplot()+
  geom_bar(aes(x = value,weight = tot_type, fill = method))+
  facet_grid(~method)+
  coord_flip()+
  theme_bw()+
  theme(legend.position = "bottom")

## Find abs max by choosing the most frequent Maxfinder scheme
max_entry <- entries %>% 
  group_by(method) %>%
  slice_max(order_by = tot_type) %>% 
  ungroup()

## In this table you have convergence among initial conditions of row positions of maximum column for each entry
## and for each estimator
max_entry 

daje.fastica <- max_entry %>% filter(method == "fastica") %>% dplyr::select(value) %>%
  str_split(string = .,pattern = "-") %>% unlist()
daje.kernelica <- max_entry %>% filter(method == "KernelICA") %>% dplyr::select(value) %>%
  str_split(string = .,pattern = "-") %>% unlist()
daje.dcov <- max_entry %>% filter(method == "DCov") %>% dplyr::select(value) %>%
  str_split(string = .,pattern = "-") %>% unlist()
daje.pml <- max_entry %>% filter(method == "pml") %>% dplyr::select(value) %>%
  str_split(string = .,pattern = "-") %>% unlist()




## Among estimated mixing matrices, select those that have the greatest frequency of maximum rowposition for each column 
index.ica <- which(lapply(c.ica, function(x) x) %>% 
  map(.f = ~ all(. == daje.fastica)) %>% do.call(rbind,.)==TRUE)
index.kica <- which(lapply(c.kica, function(x) x) %>% 
  map(.f = ~ all(. == daje.kernelica)) %>% do.call(rbind,.)==TRUE)
index.dc <- which(lapply(c.dc, function(x) x) %>% 
  map(.f = ~ all(. == daje.dcov)) %>% do.call(rbind,.)==TRUE)
index.pml <- which(lapply(c.pml, function(x) x) %>% 
  map(.f = ~ all(. == daje.pml)) %>% do.call(rbind,.)==TRUE)


## Among those mixing matrices, Find negative impact of Tax on Outuput and of Spending on Tax (labeling shock)
source(here("empirical_exercise","fin-neg.R"),local = T)

unlist(lapply(neg.ica, function(x) paste(x,collapse = "-"))) %>% as_tibble() %>% 
  group_by(value) %>% summarise(tot_type = n()) %>% 
  arrange(desc(tot_type))
unlist(lapply(neg.kica, function(x) paste(x,collapse = "-"))) %>% as_tibble() %>% 
  group_by(value) %>% summarise(tot_type = n()) %>% 
  arrange(desc(tot_type))
unlist(lapply(neg.pml, function(x) paste(x,collapse = "-"))) %>% as_tibble() %>% 
  group_by(value) %>% summarise(tot_type = n()) %>% 
  arrange(desc(tot_type))
unlist(lapply(neg.dc, function(x) paste(x,collapse = "-"))) %>% as_tibble() %>%
  group_by(value) %>% summarise(tot_type = n()) %>% 
  arrange(desc(tot_type))



## We select impact matrices that delivers negative effects on gdp and spending (we "label" this shock as increases in taxes)
boom              <- as.integer(c(4,6)) #Index of negative impact coef

## Initializations that delivers this pattern
index.neg.ica <- do.call(rbind,lapply(neg.ica, function(x) paste0(x,collapse = "-"))) %>% as.data.frame() %>% 
  mutate(neg_index = 1:n()) %>% 
  filter(str_detect(V1,"4-6")) %>% #slice(1) %>% 
  dplyr::select(neg_index) #%>% as.double()

index.neg.kica <- do.call(rbind,lapply(neg.kica, function(x) paste0(x,collapse = "-"))) %>% as.data.frame() %>% 
  mutate(neg_index = 1:n()) %>% 
  filter(str_detect(V1,"4-6")) %>% #slice(1) %>% 
  dplyr::select(neg_index) #%>% as.double()


index.neg.dc <- do.call(rbind,lapply(neg.dc, function(x) paste0(x,collapse = "-"))) %>% as.data.frame() %>% 
  mutate(neg_index = 1:n()) %>% 
  filter(str_detect(V1,"4-6")) %>% #slice(1) %>% 
  dplyr::select(neg_index) #%>% as.double()

index.neg.pml <- do.call(rbind,lapply(neg.pml, function(x) paste0(x,collapse = "-"))) %>% as.data.frame() %>% 
  mutate(neg_index = 1:n()) %>% 
  filter(str_detect(V1,"4-6")) %>% #slice(1) %>% 
  dplyr::select(neg_index) #%>% as.double()




## Chose the MATRIX!

final.index.pml <- index.pml[index.neg.pml$neg_index[1]]
final.index.ica <- index.ica[index.neg.ica$neg_index[1]]
final.index.dc <-  index.dc[index.neg.dc$neg_index[1]]
final.index.kica <-  index.kica[index.neg.kica$neg_index[1]]

Aidpml         <- A.ID.pml[[final.index.pml]]
Aidica         <- A.ID.ica[[final.index.ica]]
Aidkica         <- A.ID.kica[[final.index.kica]]
Aiddc          <- A.ID.dc[[final.index.dc]]



round(Aidica * 100, 2)
round(Aidkica * 100, 2) 
round(Aiddc * 100, 2)






### IMPACT COEFFICIENTS SIGNIFICANCE ####

### Impulse Responses #### 
## SPending Shock



coef.dc <- getcoefSVARS(x = var,B_hat = Aiddc)
irf.dc <- imrf.gianluca(y = coef.dc$y,B_hat = coef.dc$B,
                           A_hat = coef.dc$A_hat,type = coef.dc$type,
                           horizon = 20)$impulse

coef.ica <- getcoefSVARS(x = var,B_hat = Aidica)
irf.ica <- imrf.gianluca(y = coef.ica$y,B_hat = coef.ica$B,
                            A_hat = coef.ica$A_hat,type = coef.ica$type,
                            horizon = 20)$impulse

coef.kica <- getcoefSVARS(x = var,B_hat = Aidkica)
irf.kica <- imrf.gianluca(y = coef.kica$y,B_hat = coef.kica$B,
                         A_hat = coef.kica$A_hat,type = coef.kica$type,
                         horizon = 20)$impulse

coef.pml <- getcoefSVARS(x = var,B_hat = Aidpml)
irf.pml <- imrf.gianluca(y = coef.pml$y,B_hat = coef.pml$B,
                         A_hat = coef.pml$A_hat,type = coef.pml$type,
                         horizon = 20)$impulse

coef.cvm <- getcoefSVARS(x = var,B_hat = Aidcvm)
irf.cvm <- imrf.gianluca(y = coef.cvm$y,B_hat = coef.cvm$B,
                         A_hat = coef.cvm$A_hat,type = coef.cvm$type,
                         horizon = 20)$impulse


point_irf <- list(fastICA = irf.ica,
                  KernelICA = irf.kica,
                  DCov = irf.dc,
                  #PML = irf.pml,
                  CvM = irf.cvm)

save(point_irf,file = here("empirical_exercise","/point_irf_rep.RData"))



print(paste0("Bootstrap IRF: Frobenius ordering"))

boot_replications <- 500
bootstrap.ica <- mb.boot.gigi(x = coef.ica,horizon = 20,nboot = boot_replications,
                              w00 = w0[[final.index.ica]],method = "fastICA",set_seed = F,
                              ordering = "frob")

bootstrap.kica <- mb.boot.gigi(x = coef.kica,horizon = 20,nboot = boot_replications,
                               w00 = w0[[final.index.kica]],method = "KernelICA",set_seed = F,
                               ordering = "frob")


bootstrap.dc <- mb.boot.gigi(x = coef.dc,horizon = 20,nboot = boot_replications,
                             w00 = w0[[final.index.dc]],method = "Distance covariances",set_seed = F,
                             ordering = "frob")

bootstrap.pml<- mb.boot.gigi(x = coef.pml,horizon = 20,nboot = boot_replications,
                             w00 = w0[[final.index.pml]],method = "PML",set_seed = F,
                             ordering = "frob")

print(paste0("Bootstrap started for CvM"))


# bootstrap.cvm<- mb.boot.gigi(x = coef.cvm,horizon = 20,nboot = boot_replications,
#                              w00 = w0[[final.index.cvm]],method = "CvM",set_seed = F,
#                              ordering = "frob")

my_boot <- list(fastICA = bootstrap.ica,
                DCov = bootstrap.dc,
                PML = bootstrap.pml,
                kernelICA = bootstrap.kica)
#                CvM = bootstrap.cvm)

save(my_boot,file = here("empirical_exercise","bootstrap_irfs_ordering_frob_temp.RData"))







# internal function to extract and assign all information from a reduced form var object
# to be used in svar identification functions

get_var_objects <- function(x){
  
  u <- Tob <- p <- k <- residY <- coef_x <- yOut <- type <- y <-  NULL
  if(inherits(x, "var.boot")){
    assign("u", x$residuals,envir = parent.frame())
    assign("Tob", nrow(x$residuals), envir = parent.frame())
    assign("k", ncol(x$residuals), envir = parent.frame())
    assign("residY", x$residuals, envir = parent.frame())
  }else{
    assign("u", residuals(x), envir = parent.frame())
    assign("Tob", nrow(residuals(x)), envir = parent.frame())
    assign("k", ncol(residuals(x)), envir = parent.frame())
    k <- ncol(residuals(x))
    assign("residY", residuals(x), envir = parent.frame())
    
  }
  
  if(inherits(x, "var.boot")){
    assign("p", x$p, envir = parent.frame())
    assign("y", t(x$y), envir = parent.frame())
    assign("yOut", x$y, envir = parent.frame())
    assign("type", x$type, envir = parent.frame())
    assign("coef_x", x$coef_x, envir = parent.frame())
    assign("A_hat", x$coef_x, envir = parent.frame())
    
  }else if(inherits(x, "varest")){
    assign("p", x$p, envir = parent.frame())
    p <- x$p
    assign("y", t(x$y), envir = parent.frame())
    assign("yOut", x$y, envir = parent.frame())
    assign("type", x$type, envir = parent.frame())
    type <- x$type
    assign("coef_x", coef(x), envir = parent.frame())
    if(type == "none"){
      assign("A_hat", vars::Bcoef(x), envir = parent.frame())
    }else{
      assign("A_hat", vars::Bcoef(x)[, c((k * p+1):ncol(vars::Bcoef(x)),1:(k * p))], envir = parent.frame())
    }
    
  }else if(inherits(x, "nlVar")){
    assign("p", x$lag, envir = parent.frame())
    p <- x$lag
    assign("k", x$k, envir = parent.frame())
    k <- x$k
    assign("y", t(x$model[, 1:x$k]), envir = parent.frame())
    assign("yOut", x$model[, 1:x$k], envir = parent.frame())
    assign("coef_x", t(coef(x)), envir = parent.frame())
    assign("A_hat", t(coef(x)), envir = parent.frame())
    if(x$include == "const"){
      assign("A_hat", A_hat[-1,], envir = parent.frame())
    }
    
  }else if(inherits(x, "vec2var")){
    assign("k", ncol(x$resid), envir = parent.frame())
    k <- ncol(x$resid)
    coef_x <- vector("list", length = k)
    names(coef_x) <- colnames(x$y)
    assign("coef_x", coef_x, envir = parent.frame())
    assign("yOut", x$y, envir = parent.frame())
    p <- x$p
    assign("p", x$p, envir = parent.frame())
    assign("y", t(x$y), envir = parent.frame())
    
    for (i in seq_len(k)) {
      for (j in seq_len(p)) coef_x[[i]] <- c(coef_x[[i]], x$A[[j]][i,])
      coef_x[[i]] <- c(coef_x[[i]], x$deterministic[i,])
    }
    assign("coef_x", lapply(coef_x, matrix), envir = parent.frame())
    assign("type", "const", envir = parent.frame())
    
  }else{
    stop("Object class is not supported")
  }
}







id.spec = function (x, order_k = NULL) 
{
  # Extracting necessary variables from the VAR model
  u <- Tob <- p <- k <- residY <- coef_x <- yOut <- type <- y <- A_hat <- NULL
  get_var_objects(x)
  
  # Getting variable names and covariance matrix
  names_k <- colnames(yOut)
  sigg <- crossprod(u)/(Tob - 1 - k * p)
  
  # Handling the order_k parameter
  if (is.null(order_k)) {
    order_k <- 1:k
  }
  else if (is.character(order_k)) {
    order_k <- sapply(order_k, FUN = function(x) which(names_k == x))
    if (is.list(order_k)) {
      stop("Check variable names given as characters in 'order_k' or use integers instead!")
    }
  }
  names(order_k) <- names_k
  
  # Spectral decomposition
  eig <- eigen(sigg)
  P <- eig$vectors
  D <- diag(eig$values)
  B <- P %*% sqrt(D) %*% t(P)
  
  # Preparing the result
  rownames(B) <- names_k
  result <- list(B = B, A_hat = A_hat, method = "Spectral Decomposition", 
                 order_k = order_k, n = Tob, type = type, y = yOut, p = unname(p), 
                 K = k, VAR = x)
  class(result) <- "svars"
  return(result)
}


# Assuming 'irf' is your impulse response function
# 'n' is the number of variables in the VAR system
# 'horizon' is the number of periods for the forecast

# Assuming irf.kica is your IRF object from SVAR with the structure you provided

library(vars)
library(svars)




round(solve(id.chol(var.interactions)$B * 100), 2)
round(solve(id.spec(var.interactions)$B * 100), 2)



chol = id.chol(var)
spec = id.spec(var)
garch = id.garch(var)
ngml = id.ngml(var)
st = id.st(var)


point_irf$chol = irf(chol, n.ahead = 20)
point_irf$spectral = irf(spec, n.ahead = 20)
#point_irf$garch = irf(garch, n.ahead = 20)
#point_irf$ngml = irf(ngml, n.ahead = 20)
#point_irf$st = irf(st, n.ahead = 20)





library(ggplot2)
library(tidyr)
library(cowplot)

# Convert the impulse response data to a long-format data frame for ggplot
irf_data_GGTT <- bind_rows(
  data.frame(response = point_irf$fastICA$irf$`epsilon[ GG ] %->% TT`, time = 1:20, method = "FastICA"),
  data.frame(response = point_irf$KernelICA$irf$`epsilon[ GG ] %->% TT`, time = 1:20, method = "KernelICA"),
  # data.frame(response = point_irf$DCov$irf$`epsilon[ GG ] %->% TT`, time = 1:20, method = "DCoV"),
  # data.frame(response = point_irf$PML$irf$`epsilon[ GG ] %->% TT`, time = 1:20, method = "Pseudo-ML"),
  # data.frame(response = point_irf$CvM$irf$`epsilon[ GG ] %->% TT`, time = 1:20, method = "CvM"),
  data.frame(response = point_irf$chol$irf$`epsilon[ GG ] %->% TT`, time = 1:20, method = "Cholesky"),
  data.frame(response = point_irf$spectral$irf$`epsilon[ GG ] %->% TT`, time = 1:20, method = "Spectral"),
  #data.frame(response = point_irf$garch$irf$`epsilon[ GG ] %->% TT`, time = 1:20, method = "Garch"),
  #data.frame(response = point_irf$ngml$irf$`epsilon[ GG ] %->% TT`, time = 1:20, method = "Non-gaussian ML"),
  #data.frame(response = point_irf$st$irf$`epsilon[ GG ] %->% TT`, time = 1:20, method = "Smooth transition") 
)



# Plot
p1 = ggplot(data = irf_data_GGTT, aes(x = time, y = response, color = method)) +
  geom_line() +
  labs(title = "Tax",
       x = "Quarters", y = "Response") +
  theme_minimal() +
  theme(legend.position = "right", # Place legend on right
        legend.direction = "vertical") # Make legend vertical

p1 <- p1 + theme(legend.position = "bottom", legend.direction = "horizontal")
legend = cowplot::get_legend(p1)


# Convert the impulse response data to a long-format data frame for ggplot
irf_data_GGGG <- bind_rows(
  data.frame(response = point_irf$fastICA$irf$`epsilon[ GG ] %->% GG`, time = 1:20, method = "fastICA"),
  data.frame(response = point_irf$KernelICA$irf$`epsilon[ GG ] %->% GG`, time = 1:20, method = "KernelICA"),
  # data.frame(response = point_irf$DCov$irf$`epsilon[ GG ] %->% GG`, time = 1:20, method = "DCov"),
  # data.frame(response = point_irf$PML$irf$`epsilon[ GG ] %->% GG`, time = 1:20, method = "PML"),
  # data.frame(response = point_irf$CvM$irf$`epsilon[ GG ] %->% GG`, time = 1:20, method = "CvM"),
  data.frame(response = point_irf$chol$irf$`epsilon[ GG ] %->% GG`, time = 1:20, method = "chol"),
  data.frame(response = point_irf$spectral$irf$`epsilon[ GG ] %->% GG`, time = 1:20, method = "spectral"),
  #data.frame(response = point_irf$garch$irf$`epsilon[ GG ] %->% GG`, time = 1:20, method = "garch"),
  #data.frame(response = point_irf$ngml$irf$`epsilon[ GG ] %->% GG`, time = 1:20, method = "ngml"),
  #data.frame(response = point_irf$st$irf$`epsilon[ GG ] %->% GG`, time = 1:20, method = "st") 
)

# Plot
p2 = ggplot(data = irf_data_GGGG, aes(x = time, y = response, color = method)) +
  geom_line() +
  labs(title = "Spending",
       x = "Quarters", y = "Response") +
  theme_minimal() +
  theme(legend.position = "bottom")




# Convert the impulse response data to a long-format data frame for ggplot
irf_data_GGXX <- bind_rows(
  data.frame(response = point_irf$fastICA$irf$`epsilon[ GG ] %->% XX`, time = 1:20, method = "fastICA"),
  data.frame(response = point_irf$KernelICA$irf$`epsilon[ GG ] %->% XX`, time = 1:20, method = "KernelICA"),
  # data.frame(response = point_irf$DCov$irf$`epsilon[ GG ] %->% XX`, time = 1:20, method = "DCov"),
  # data.frame(response = point_irf$PML$irf$`epsilon[ GG ] %->% XX`, time = 1:20, method = "PML"),
  # data.frame(response = point_irf$CvM$irf$`epsilon[ GG ] %->% XX`, time = 1:20, method = "CvM"),
  data.frame(response = point_irf$chol$irf$`epsilon[ GG ] %->% XX`, time = 1:20, method = "chol"),
  data.frame(response = point_irf$spectral$irf$`epsilon[ GG ] %->% XX`, time = 1:20, method = "spectral"),
  #data.frame(response = point_irf$garch$irf$`epsilon[ GG ] %->% XX`, time = 1:20, method = "garch"),
  # data.frame(response = point_irf$ngml$irf$`epsilon[ GG ] %->% XX`, time = 1:20, method = "ngml"),
  #data.frame(response = point_irf$st$irf$`epsilon[ GG ] %->% XX`, time = 1:20, method = "st") 
)

# Plot
p3 = ggplot(data = irf_data_GGXX, aes(x = time, y = response, color = method)) +
  geom_line() +
  labs(title = "GDP",
       x = "Quarters", y = "Response") +
  theme_minimal() +
  theme(legend.position = "bottom")



# Convert the impulse response data to a long-format data frame for ggplot
irf_data_XXTT <- bind_rows(
  data.frame(response = point_irf$fastICA$irf$`epsilon[ XX ] %->% TT`, time = 1:20, method = "fastICA"),
  data.frame(response = point_irf$KernelICA$irf$`epsilon[ XX ] %->% TT`, time = 1:20, method = "KernelICA"),
  # data.frame(response = point_irf$DCov$irf$`epsilon[ XX ] %->% TT`, time = 1:20, method = "DCov"),
  # data.frame(response = point_irf$PML$irf$`epsilon[ XX ] %->% TT`, time = 1:20, method = "PML"),
  # data.frame(response = point_irf$CvM$irf$`epsilon[ XX ] %->% TT`, time = 1:20, method = "CvM"),
  data.frame(response = point_irf$chol$irf$`epsilon[ XX ] %->% TT`, time = 1:20, method = "chol"),
  data.frame(response = point_irf$spectral$irf$`epsilon[ XX ] %->% TT`, time = 1:20, method = "spectral"),
  #data.frame(response = point_irf$garch$irf$`epsilon[ XX ] %->% TT`, time = 1:20, method = "garch"),
  #data.frame(response = point_irf$ngml$irf$`epsilon[ XX ] %->% TT`, time = 1:20, method = "ngml"),
  #data.frame(response = point_irf$st$irf$`epsilon[ XX ] %->% TT`, time = 1:20, method = "st") 
)

# Plot
p4 = ggplot(data = irf_data_XXTT, aes(x = time, y = response, color = method)) +
  geom_line() +
  labs(title = "GDP on tax",
       x = "Quarters", y = "Response") +
  theme_minimal() +
  theme(legend.position = "bottom")





# Convert the impulse response data to a long-format data frame for ggplot
irf_data_XXXX <- bind_rows(
  data.frame(response = point_irf$fastICA$irf$`epsilon[ XX ] %->% XX`, time = 1:20, method = "fastICA"),
  data.frame(response = point_irf$KernelICA$irf$`epsilon[ XX ] %->% XX`, time = 1:20, method = "KernelICA"),
  # data.frame(response = point_irf$DCov$irf$`epsilon[ XX ] %->% XX`, time = 1:20, method = "DCov"),
  # data.frame(response = point_irf$PML$irf$`epsilon[ XX ] %->% XX`, time = 1:20, method = "PML"),
  # data.frame(response = point_irf$CvM$irf$`epsilon[ XX ] %->% XX`, time = 1:20, method = "CvM"),
  data.frame(response = point_irf$chol$irf$`epsilon[ XX ] %->% XX`, time = 1:20, method = "chol"),
  data.frame(response = point_irf$spectral$irf$`epsilon[ XX ] %->% XX`, time = 1:20, method = "spectral"),
  #data.frame(response = point_irf$garch$irf$`epsilon[ XX ] %->% XX`, time = 1:20, method = "garch"),
  #data.frame(response = point_irf$ngml$irf$`epsilon[ XX ] %->% XX`, time = 1:20, method = "ngml"),
  #data.frame(response = point_irf$st$irf$`epsilon[ XX ] %->% XX`, time = 1:20, method = "st") 
)

# Plot
p5 = ggplot(data = irf_data_XXXX, aes(x = time, y = response, color = method)) +
  geom_line() +
  labs(title = "GDP on GDP",
       x = "Quarters", y = "Response") +
  theme_minimal() +
  theme(legend.position = "bottom")



# Convert the impulse response data to a long-format data frame for ggplot
irf_data_XXGG <- bind_rows(
  data.frame(response = point_irf$fastICA$irf$`epsilon[ XX ] %->% GG`, time = 1:20, method = "fastICA"),
  data.frame(response = point_irf$KernelICA$irf$`epsilon[ XX ] %->% GG`, time = 1:20, method = "KernelICA"),
  # data.frame(response = point_irf$DCov$irf$`epsilon[ XX ] %->% GG`, time = 1:20, method = "DCov"),
  # data.frame(response = point_irf$PML$irf$`epsilon[ XX ] %->% GG`, time = 1:20, method = "PML"),
  # data.frame(response = point_irf$CvM$irf$`epsilon[ XX ] %->% GG`, time = 1:20, method = "CvM"),
  data.frame(response = point_irf$chol$irf$`epsilon[ XX ] %->% GG`, time = 1:20, method = "chol"),
  data.frame(response = point_irf$spectral$irf$`epsilon[ XX ] %->% GG`, time = 1:20, method = "spectral"),
  #data.frame(response = point_irf$garch$irf$`epsilon[ XX ] %->% GG`, time = 1:20, method = "garch"),
  #data.frame(response = point_irf$ngml$irf$`epsilon[ XX ] %->% GG`, time = 1:20, method = "ngml"),
  #data.frame(response = point_irf$st$irf$`epsilon[ XX ] %->% GG`, time = 1:20, method = "st") 
)

# Plot
p6 = ggplot(data = irf_data_XXGG, aes(x = time, y = response, color = method)) +
  geom_line() +
  labs(title = "GDP on spending ",
       x = "Quarters", y = "Response") +
  theme_minimal() +
  theme(legend.position = "bottom")




# Convert the impulse response data to a long-format data frame for ggplot
irf_data_TTTT <- bind_rows(
  data.frame(response = point_irf$fastICA$irf$`epsilon[ TT ] %->% TT`, time = 1:20, method = "fastICA"),
  data.frame(response = point_irf$KernelICA$irf$`epsilon[ TT ] %->% TT`, time = 1:20, method = "KernelICA"),
  # data.frame(response = point_irf$DCov$irf$`epsilon[ TT ] %->% TT`, time = 1:20, method = "DCov"),
  # data.frame(response = point_irf$PML$irf$`epsilon[ TT ] %->% TT`, time = 1:20, method = "PML"),
  # data.frame(response = point_irf$CvM$irf$`epsilon[ TT ] %->% TT`, time = 1:20, method = "CvM"),
  data.frame(response = point_irf$chol$irf$`epsilon[ TT ] %->% TT`, time = 1:20, method = "chol"),
  data.frame(response = point_irf$spectral$irf$`epsilon[ TT ] %->% TT`, time = 1:20, method = "spectral"),
  #data.frame(response = point_irf$garch$irf$`epsilon[ TT ] %->% TT`, time = 1:20, method = "garch"),
  #data.frame(response = point_irf$ngml$irf$`epsilon[ TT ] %->% TT`, time = 1:20, method = "ngml"),
  #data.frame(response = point_irf$st$irf$`epsilon[ TT ] %->% TT`, time = 1:20, method = "st") 
)

# Plot
p7 = ggplot(data = irf_data_TTTT, aes(x = time, y = response, color = method)) +
  geom_line() +
  labs(title = "Tax",
       x = "Quarters", y = "Response") +
  theme_minimal() +
  theme(legend.position = "bottom")



# Convert the impulse response data to a long-format data frame for ggplot
irf_data_TTXX <- bind_rows(
  data.frame(response = point_irf$fastICA$irf$`epsilon[ TT ] %->% XX`, time = 1:20, method = "fastICA"),
  data.frame(response = point_irf$KernelICA$irf$`epsilon[ TT ] %->% XX`, time = 1:20, method = "KernelICA"),
  # data.frame(response = point_irf$DCov$irf$`epsilon[ TT ] %->% XX`, time = 1:20, method = "DCov"),
  # data.frame(response = point_irf$PML$irf$`epsilon[ TT ] %->% XX`, time = 1:20, method = "PML"),
  # data.frame(response = point_irf$CvM$irf$`epsilon[ TT ] %->% XX`, time = 1:20, method = "CvM"),
  data.frame(response = point_irf$chol$irf$`epsilon[ TT ] %->% XX`, time = 1:20, method = "chol"),
  data.frame(response = point_irf$spectral$irf$`epsilon[ TT ] %->% XX`, time = 1:20, method = "spectral"),
  #data.frame(response = point_irf$garch$irf$`epsilon[ TT ] %->% XX`, time = 1:20, method = "garch"),
  #data.frame(response = point_irf$ngml$irf$`epsilon[ TT ] %->% XX`, time = 1:20, method = "ngml"),
  #data.frame(response = point_irf$st$irf$`epsilon[ TT ] %->% XX`, time = 1:20, method = "st") 
)

# Plot
p8 = ggplot(data = irf_data_TTXX, aes(x = time, y = response, color = method)) +
  geom_line() +
  labs(title = "GDP",
       x = "Quarters", y = "Response") +
  theme_minimal() +
  theme(legend.position = "bottom")



# Convert the impulse response data to a long-format data frame for ggplot
irf_data_TTGG <- bind_rows(
  data.frame(response = point_irf$fastICA$irf$`epsilon[ TT ] %->% GG`, time = 1:20, method = "fastICA"),
  data.frame(response = point_irf$KernelICA$irf$`epsilon[ TT ] %->% GG`, time = 1:20, method = "KernelICA"),
  # data.frame(response = point_irf$DCov$irf$`epsilon[ TT ] %->% GG`, time = 1:20, method = "DCov"),
  # data.frame(response = point_irf$PML$irf$`epsilon[ TT ] %->% GG`, time = 1:20, method = "PML"),
  # data.frame(response = point_irf$CvM$irf$`epsilon[ TT ] %->% GG`, time = 1:20, method = "CvM"),
  data.frame(response = point_irf$chol$irf$`epsilon[ TT ] %->% GG`, time = 1:20, method = "chol"),
  data.frame(response = point_irf$spectral$irf$`epsilon[ TT ] %->% GG`, time = 1:20, method = "spectral"),
  #data.frame(response = point_irf$garch$irf$`epsilon[ TT ] %->% GG`, time = 1:20, method = "garch"),
  #data.frame(response = point_irf$ngml$irf$`epsilon[ TT ] %->% GG`, time = 1:20, method = "ngml"),
  #data.frame(response = point_irf$st$irf$`epsilon[ TT ] %->% GG`, time = 1:20, method = "st") 
)

# Plot
p9 = ggplot(data = irf_data_TTGG, aes(x = time, y = response, color = method)) +
  geom_line() +
  labs(title = "Spending",
       x = "Quarters", y = "Response") +
  theme_minimal() +
  theme(legend.position = "bottom")



legend = cowplot::get_legend(p1)

p1 <- p1 + theme(legend.position="none", plot.title = element_text(hjust = 0.5))
p2 <- p2 + theme(legend.position="none", plot.title = element_text(hjust = 0.5))
p3 <- p3 + theme(legend.position="none", plot.title = element_text(hjust = 0.5))
p4 <- p4 + theme(legend.position="none", plot.title = element_text(hjust = 0.5))
p5 <- p5 + theme(legend.position="none", plot.title = element_text(hjust = 0.5))
p6 <- p6 + theme(legend.position="none", plot.title = element_text(hjust = 0.5))
p7 <- p7 + theme(legend.position="none", plot.title = element_text(hjust = 0.5))
p8 <- p8 + theme(legend.position="none", plot.title = element_text(hjust = 0.5))
p9 <- p9 + theme(legend.position="none", plot.title = element_text(hjust = 0.5))





combined_plots = grid.arrange(grobs = list(p1, p2, p3, p7, p8, p9), ncol = 3)





# Assume combined_plots is your combined ggplot objects without legends
# and legend is the extracted legend

# Define the layout matrix for arranging the combined plot and the legend
layout_matrix <- rbind(
  c(1, 1, 1),
  c(2, 2, 2)
)

# Arrange everything in a grid with the legend at the bottom
# Adjust the heights to control the space allocated to the plots vs. the legend
grid.arrange(combined_plots, legend, layout_matrix = layout_matrix, heights = c(5, 1))




########## FEVD ###########################


########## fved for kica


# Extracting the IRFs
irfs_kica <- irf.kica$irf[, -1]  # Remove the first column as it is not needed

# Number of variables and time periods
n_vars <- 3  # Number of variables in the SVAR model
n_periods <- nrow(irfs_kica)

# Reshape the IRFs for computation
irfs_matrix_kica <- array(data = NA, dim = c(n_periods, n_vars, n_vars))
for (i in 1:n_vars) {
  for (j in 1:n_vars) {
    irfs_matrix_kica[, i, j] <- irfs_kica[, (i - 1) * n_vars + j]
  }
}







# Compute Forecast Variance Error Decomposition (FVED)
fevd_kica <- array(data = NA, dim = c(n_periods, n_vars, n_vars))
for (t in 1:n_periods) {
  for (i in 1:n_vars) {
    for (j in 1:n_vars) {
      irf_squared_cumsum_kica <- sum(irfs_matrix_kica[1:t, i, j]^2)
      total_variance_kica <- sum(irfs_matrix_kica[1:t, i, ]^2)
      fevd_kica[t, i, j] <- irf_squared_cumsum_kica / total_variance_kica
    }
  }
}

fevd_kica_XX = cbind(fevd_kica[,,1][,3],fevd_kica[,,2][,3],fevd_kica[,,3][,3])|> as.data.frame()

colnames(fevd_kica_XX) = c("GG", "TT", "XX")


########## fved for ica


# Extracting the IRFs
irfs_ica <- irf.ica$irf[, -1]  # Remove the first column as it is not needed

# Number of variables and time periods
n_vars <- 3  # Number of variables in the SVAR model
n_periods <- nrow(irfs_ica)

# Reshape the IRFs for computation
irfs_matrix_ica <- array(data = NA, dim = c(n_periods, n_vars, n_vars))
for (i in 1:n_vars) {
  for (j in 1:n_vars) {
    irfs_matrix_ica[, i, j] <- irfs_ica[, (i - 1) * n_vars + j]
  }
}







# Compute Forecast Variance Error Decomposition (FVED)
fevd_ica <- array(data = NA, dim = c(n_periods, n_vars, n_vars))
for (t in 1:n_periods) {
  for (i in 1:n_vars) {
    for (j in 1:n_vars) {
      irf_squared_cumsum_ica <- sum(irfs_matrix_ica[1:t, i, j]^2)
      total_variance_ica <- sum(irfs_matrix_ica[1:t, i, ]^2)
      fevd_ica[t, i, j] <- irf_squared_cumsum_ica / total_variance_ica
    }
  }
}

fevd_ica_XX = cbind(fevd_ica[,,1][,3],fevd_ica[,,2][,3],fevd_ica[,,3][,3])|> as.data.frame()

colnames(fevd_ica_XX) = c("GG", "TT", "XX")



########## fved for cvm


# Extracting the IRFs
irfs_cvm <- irf.cvm$irf[, -1]  # Remove the first column as it is not needed

# Number of variables and time periods
n_vars <- 3  # Number of variables in the SVAR model
n_periods <- nrow(irfs_cvm)

# Reshape the IRFs for computation
irfs_matrix_cvm <- array(data = NA, dim = c(n_periods, n_vars, n_vars))
for (i in 1:n_vars) {
  for (j in 1:n_vars) {
    irfs_matrix_cvm[, i, j] <- irfs_cvm[, (i - 1) * n_vars + j]
  }
}







# Compute Forecast Variance Error Decomposition (FVED)
fevd_cvm <- array(data = NA, dim = c(n_periods, n_vars, n_vars))
for (t in 1:n_periods) {
  for (i in 1:n_vars) {
    for (j in 1:n_vars) {
      irf_squared_cumsum_cvm <- sum(irfs_matrix_cvm[1:t, i, j]^2)
      total_variance_cvm <- sum(irfs_matrix_cvm[1:t, i, ]^2)
      fevd_cvm[t, i, j] <- irf_squared_cumsum_cvm / total_variance_cvm
    }
  }
}

fevd_cvm_XX = cbind(fevd_cvm[,,1][,3],fevd_cvm[,,2][,3],fevd_cvm[,,3][,3])|> as.data.frame()

colnames(fevd_cvm_XX) = c("GG", "TT", "XX")




########## fved for dc


# Extracting the IRFs
irfs_dc <- irf.dc$irf[, -1]  # Remove the first column as it is not needed

# Number of variables and time periods
n_vars <- 3  # Number of variables in the SVAR model
n_periods <- nrow(irfs_dc)

# Reshape the IRFs for computation
irfs_matrix_dc <- array(data = NA, dim = c(n_periods, n_vars, n_vars))
for (i in 1:n_vars) {
  for (j in 1:n_vars) {
    irfs_matrix_dc[, i, j] <- irfs_dc[, (i - 1) * n_vars + j]
  }
}







# Compute Forecast Variance Error Decomposition (FVED)
fevd_dc <- array(data = NA, dim = c(n_periods, n_vars, n_vars))
for (t in 1:n_periods) {
  for (i in 1:n_vars) {
    for (j in 1:n_vars) {
      irf_squared_cumsum_dc <- sum(irfs_matrix_dc[1:t, i, j]^2)
      total_variance_dc <- sum(irfs_matrix_dc[1:t, i, ]^2)
      fevd_dc[t, i, j] <- irf_squared_cumsum_dc / total_variance_dc
    }
  }
}

fevd_dc_XX = cbind(fevd_dc[,,1][,3],fevd_dc[,,2][,3],fevd_dc[,,3][,3])|> as.data.frame()

colnames(fevd_dc_XX) = c("GG", "TT", "XX")












n.ahead = 20

fevd_chol = fevd(chol, n.ahead = n.ahead)
fevd_spectral = fevd(spec, n.ahead = n.ahead)
fevd_garch = fevd(garch, n.ahead = n.ahead)
fevd_ngml = fevd(ngml, n.ahead = n.ahead)
fevd_st = fevd(st, n.ahead = n.ahead)



data_chol <- as.data.frame(fevd_chol$XX)
data_spectral <- as.data.frame(fevd_spectral$XX)
data_garch <- as.data.frame(fevd_garch$XX)
data_ngml <- as.data.frame(fevd_ngml$XX)
data_st <- as.data.frame(fevd_st$XX)



combined_data <- data.frame(Time = rownames(data_chol), Chol = data_chol[,1], Spectral = data_spectral[,1], Garch = data_garch[,1], Ngml = data_ngml[,1], St = data_st[,1], kica = fevd_kica_XX[,1]*100, ica = fevd_ica_XX[,1]*100, cvm = fevd_cvm_XX[,1]*100, dc = fevd_dc_XX[,1]*100)
combined_data <- data.frame(Time = rownames(data_chol), Chol = data_chol[,1], Spectral = data_spectral[,1], kica = fevd_kica_XX[,1]*100, ica = fevd_ica_XX[,1]*100)

long_data <- melt(combined_data, id = "Time")

long_data$Time <- as.numeric(as.character(long_data$Time))


colnames(long_data)[which(colnames(long_data) == "variable")] <- "method"


# Create the plot
p1 = ggplot(long_data, aes(x = Time, y = value, color = method, group = method)) +
  geom_line() +
  labs(x = "Forecast Horizon", y = "Variance Proportion", title = "Forecast Error variance of GDP explained by spending") +
  theme_minimal()+
  scale_color_discrete(labels = c("Chol" = "Cholesky", 
                                  "Spectral" = "Spectral", 
                                  "kica" = "KernelICA", 
                                  "ica" = "FastICA"))


combined_data <- data.frame(Time = rownames(data_chol), Chol = data_chol[,2], Spectral = data_spectral[,2], Garch = data_garch[,2], Ngml = data_ngml[,2], St = data_st[,2], kica = fevd_kica_XX[,2]*100, ica = fevd_ica_XX[,2]*100, cvm = fevd_cvm_XX[,2]*100, dc = fevd_dc_XX[,2]*100)
combined_data <- data.frame(Time = rownames(data_chol), Chol = data_chol[,2], Spectral = data_spectral[,2], kica = fevd_kica_XX[,2]*100, ica = fevd_ica_XX[,2]*100)
long_data <- melt(combined_data, id = "Time")

long_data$Time <- as.numeric(as.character(long_data$Time))


colnames(long_data)[which(colnames(long_data) == "variable")] <- "method"


# Create the plot
p2 = ggplot(long_data, aes(x = Time, y = value, color = method, group = method)) +
  geom_line() +
  labs(x = "Forecast Horizon", y = "Variance Proportion", title = "Forecast Error variance of GDP explained by taxes") +
  theme_minimal() +
  scale_color_discrete(labels = c("Chol" = "Cholesky", 
                                  "Spectral" = "Spectral Method", 
                                  "kica" = "Kernel ICA", 
                                  "ica" = "Independent Component Analysis"))

long_data$Time <- as.numeric(as.character(long_data$Time))

# Create the plot
p3 = ggplot(long_data, aes(x = Time, y = value, color = variable, group = variable)) +
  geom_line() +
  labs(x = "Forecast Horizon", y = "Variance Proportion", title = "Forecast Error variance of GDP explained by GDP") +
  theme_minimal()


legend = cowplot::get_legend(p1)

p1 = p1 + theme(legend.position="none")
p2 = p2 + theme(legend.position="none")



combined_plot = grid.arrange(p1, p2)

# Arrange the combined plot and the legend
grid.arrange(combined_plot, legend, ncol=2, widths=c(4, 1))








######################



# 
# # Plotting
# ggplot(fved_long, aes(x = Time, y = FVED, colour = Shock)) +
#   geom_line() +
#   facet_wrap(~ Variable, ncol = 1, scales = "free_y") +
#   theme_minimal() +
#   labs(title = "Forecast Variance Error Decomposition (FVED)",
#        x = "Time",
#        y = "FVED",
#        colour = "Shock")
# 
# 

















print(paste0("Bootstrap IRF: Frobenius ordering"))

boot_replications <- 500
bootstrap.ica <- mb.boot.gigi(x = coef.ica,horizon = 20,nboot = boot_replications,
                              w00 = w0[[final.index.ica]],method = "fastICA",set_seed = F,
                              ordering = "frob")

bootstrap.kica <- mb.boot.gigi(x = coef.kica,horizon = 20,nboot = boot_replications,
                              w00 = w0[[final.index.kica]],method = "KernelICA",set_seed = F,
                              ordering = "frob")


bootstrap.dc <- mb.boot.gigi(x = coef.dc,horizon = 20,nboot = boot_replications,
                             w00 = w0[[final.index.dc]],method = "Distance covariances",set_seed = F,
                             ordering = "frob")

bootstrap.pml<- mb.boot.gigi(x = coef.pml,horizon = 20,nboot = boot_replications,
                             w00 = w0[[final.index.pml]],method = "PML",set_seed = F,
                             ordering = "frob")




#print(paste0("Bootstrap started for CvM"))


# bootstrap.cvm<- mb.boot.gigi(x = coef.cvm,horizon = 20,nboot = boot_replications,
#                              w00 = w0[[final.index.cvm]],method = "CvM",set_seed = F,
#                              ordering = "frob")

my_boot <- list(fastICA = bootstrap.ica,
                KernelICA = bootstrap.kica,
                DCov = bootstrap.dc,
                PML = bootstrap.pml
                #,CvM = bootstrap.cvm
                )

save(my_boot,file = here("empirical_exercise","bootstrap_irfs_ordering_frob_temp.RData"))











