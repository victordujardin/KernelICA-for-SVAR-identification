## Blanchard - Perotti exercise with ICA identification

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

plot(density(u.ica[,1]))
plot(density(u.ica[,2]))
plot(density(u.ica[,3]))


### Orthogonalize Reduced-form residuals
sigg   <- cov((u.ica)) ## Variance-Covariance Matrix
C      <- t(chol(sigg)) ## Lower-triangular Choleski factor
u_chol <- solve(C) %*% t(u.ica) ### Orthogonalized reduced form residuals (Preliminary Rotation)
tu_chol<- t(solve(C) %*% t(u.ica))

## Call function to maximize
source(here("empirical_exercise","/maxlik3D-pml.R"),local = T)


## Select Best initial conditions
n_lhs  <- 500 ## number of initializations

## Run initialization for FastICA, DCov and PML + CVM estimation
source(here("empirical_exercise","find-init-and-cvm.R"),local = T)

## Identified impact matrix under CvM
Aidcvm <- maxfinder(A = Acvm)$A.id

## Find the maximum entry per each column
for(i in 1:n_lhs){
  c.ica[[i]] <- A.ID.ica[[i]] %>% apply(MARGIN = 2, function(x) which.max(abs(x)))
  c.dc[[i]]  <- A.ID.dc[[i]]  %>% apply(MARGIN = 2, function(x) which.max(abs(x)))
  c.pml[[i]] <- A.ID.pml[[i]] %>% apply(MARGIN = 2, function(x) which.max(abs(x)))
}

entries <- bind_rows(
unlist(lapply(c.ica, function(x) paste(as.character(x),collapse = "-"))) %>% as_tibble() %>% 
  group_by(value) %>% summarise(tot_type = n()) %>% 
  mutate(method = "fastica"),
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
daje.dcov <- max_entry %>% filter(method == "DCov") %>% dplyr::select(value) %>%
  str_split(string = .,pattern = "-") %>% unlist()
daje.pml <- max_entry %>% filter(method == "pml") %>% dplyr::select(value) %>%
  str_split(string = .,pattern = "-") %>% unlist()


## Among estimated mixing matrices, select those that have the greatest frequency of maximum rowposition for each column 
index.ica <- which(lapply(c.ica, function(x) x) %>% 
  map(.f = ~ all(. == daje.fastica)) %>% do.call(rbind,.)==TRUE)
index.dc <- which(lapply(c.dc, function(x) x) %>% 
  map(.f = ~ all(. == daje.dcov)) %>% do.call(rbind,.)==TRUE)
index.pml <- which(lapply(c.pml, function(x) x) %>% 
                    map(.f = ~ all(. == daje.pml)) %>% do.call(rbind,.)==TRUE)

## Among those mixing matrices, Find negative impact of Tax on Outuput and of Spending on Tax (labeling shock)
source(here("empirical_exercise","fin-neg.R"),local = T)

unlist(lapply(neg.ica, function(x) paste(x,collapse = "-"))) %>% as_tibble() %>% 
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
Aidpml         <- A.ID.pml[[final.index.pml]]
Aidica         <- A.ID.ica[[final.index.ica]]
Aiddc          <- A.ID.dc[[final.index.dc]]



### IMPACT COEFFICIENTS SIGNIFICANCE ####

### Impulse Responses #### 
## SPending Shock

coef.dc <- getcoefSVARS(x = var.irf,B_hat = Aiddc)
irf.dc <- imrf.gianluca(y = coef.dc$y,B_hat = coef.dc$B,
                           A_hat = coef.dc$A_hat,type = coef.dc$type,
                           horizon = 20)$impulse

coef.ica <- getcoefSVARS(x = var.irf,B_hat = Aidica)
irf.ica <- imrf.gianluca(y = coef.ica$y,B_hat = coef.ica$B,
                            A_hat = coef.ica$A_hat,type = coef.ica$type,
                            horizon = 20)$impulse

coef.pml <- getcoefSVARS(x = var.irf,B_hat = Aidpml)
irf.pml <- imrf.gianluca(y = coef.pml$y,B_hat = coef.pml$B,
                         A_hat = coef.pml$A_hat,type = coef.pml$type,
                         horizon = 20)$impulse

coef.cvm <- getcoefSVARS(x = var.irf,B_hat = Aidcvm)
irf.cvm <- imrf.gianluca(y = coef.cvm$y,B_hat = coef.cvm$B,
                         A_hat = coef.cvm$A_hat,type = coef.cvm$type,
                         horizon = 20)$impulse


point_irf <- list(fastICA = irf.ica,
                  DCov = irf.dc,
                  PML = irf.pml,
                  CvM = irf.cvm)

save(point_irf,file = here("empirical_exercise","/point_irf_rep.RData"))


print(paste0("Bootstrap IRF: Frobenius ordering"))

boot_replications <- 500
bootstrap.ica <- mb.boot.gigi(x = coef.ica,horizon = 20,nboot = boot_replications,
                              w00 = w0[[final.index.ica]],method = "fastICA",set_seed = F,
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
                PML = bootstrap.pml)
#                CvM = bootstrap.cvm)

save(my_boot,file = here("empirical_exercise","bootstrap_irfs_ordering_frob_temp.RData"))

