## Code to replicate Figure 2 of the paper
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


## Colors
palette <- c('CvM'     = '#66c2a5',
             'fastICA' = '#fc8d62',
             'DCov'    = '#8da0cb',
             'PML'     = '#e78ac3')

## Parameters of the p-generalized distribution
## Number of shape parameter values (in the specific assessment and for inference table we only focuse on 4 distributional scenario - See paper)
seq1   <- seq(from = 0.5, to = 3.5, length.out = 15) ## Sequence of p-Generalized parameter for general evaluation exercise (1st part)
seq2   <- seq(from = 4, to = 100, length.out = 5) ## Sequence of p-Generalized parameter for general evaluation exercise (2nd part)
#SEQ     <- c(0.5,1.57,2.43,100) ## Sequence of p-Generalized parameter for statistical inference exercise
SEQ    <- as.numeric(c(seq1, seq2))
rm(seq1,seq2)

## An arbitrary 2x2 Mixing Matrix
Btrue2Dz <- matrix(c(1.14, -0.38,
                     0,    1.26), nrow = 2, ncol = 2)

## An arbitrary 3x3 Mixing Matrix
## This Btrue3D has an inverse (the unmixing matrix containing contemporaneous relationships) 
## with an entry equal to 0
Btr3D  <- matrix(c(0.9,   0.150, 0.65,
                   -0.75,   1.13, 0.22,
                   0.21,  -0.53, 1.5)  , nrow = 3, ncol = 3)

Wtrue3D   <- t(solve(Btr3D))
Wtrue3D[3,1]<-0
Btrue3D   <- solve(Wtrue3D)


## Load Databases with simulated data
average_performance <- read_csv(file = here("final_databases","average_performance.csv")) #%>% 
  #dplyr::select(-"X1")


Mixing_True <- list(Btrue2Dz = Btrue2Dz,
                    Btrue3D  = Btrue3D)

percentiles <- NULL
for (cc in 1:2) {
  MDI <- double(length =  length(unique(average_performance$mc)))
  for (i in 1:length(unique(average_performance$mc))) {
    ## Function taken from ProDenICA::mixmat to generate
    ## a "well-behaving" mixing matrix with condition number C
    p <- cc+1
    set.seed(1+i)
    a <- matrix(rnorm(p * p), p, p)
    sa <- svd(a)
    #set.seed(1000+i)
    d <- sort(runif(p) + 1)
    mat <- sa$u %*% (sa$v * d)
    cond_number <- d[p]/d[1]
    
    MDI[i] <- MD(W.hat = solve(mat),A = Mixing_True[[cc]])
  }
  percentiles <- bind_rows(as.data.frame(MDI) %>% 
                             mutate(dimension = paste0("k=",cc+1)),percentiles)
  
}



## Bind files for CVM

rm(p,a,sa,d,mat,cond_number,Btr3D)

distr <- percentiles %>% as_tibble() %>% 
  ggplot()+
  geom_density(aes(x = MDI, group = dimension))+
  facet_grid(~dimension)

## Generate quantile of MDI distribution
quantile_rnd_matrix <- percentiles %>% as_tibble() %>% 
  group_by(dimension) %>% 
  summarise(percentile_5 = quantile(MDI,probs = 0.1))


## Colors
palette <- c('CvM'     = '#66c2a5',
             'fastICA' = '#fc8d62',
             'DCov'    = '#8da0cb',
             'PML'     = '#e78ac3')

shapes <- c('CvM'     = 21,
            'fastICA' = 22,
            'DCov'    = 24,
            'PML'     = 23)

average_performance %>% 
  mutate(scenario = str_replace(scenario,"n","T")) %>% filter(scenario == "T=400") %>% 
  rename(MDI = "Minimum Distance") %>% 
  group_by(estimator,p_shape,scenario,dimension) %>% 
  summarise_at(vars(MDI), list(mean = ~ mean(.),
                               sd   = ~ sd(.))) %>% 
  ungroup() %>% 
  left_join(.,quantile_rnd_matrix) %>% 
  mutate(p_shape = as.factor(p_shape),
         estimator = factor(estimator,levels = c("PML","fastICA","DCov","CvM"))) %>% 
  
  ggplot(.,aes(x = p_shape, group = estimator))+
  geom_line(aes(y = mean,color = estimator), size = 1)+
  geom_point(aes(y = mean,color = estimator, shape = estimator, fill = estimator))+
  geom_hline(aes(yintercept = percentile_5),linetype = "dotted")+
  geom_ribbon(aes(ymin = mean - sd,ymax = mean + sd, fill = estimator), alpha = 0.3)+
  facet_grid(dimension ~ estimator)+
  scale_color_manual(values = palette,breaks = c("PML","fastICA","DCov","CvM"))+
  scale_fill_manual(values = palette , breaks = c("PML","fastICA","DCov","CvM"))+
  scale_shape_manual(values = shapes , breaks = c("PML","fastICA","DCov","CvM")) + 
  scale_x_discrete(breaks = c("0.5","1.57","2","2.43","4","52","100"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "white"))+
  labs(title = "General Assessment",  y = "MDI", 
       fill = " ",color = " ",shape = " ")


different_figure <- average_performance %>% 
  filter(scenario == "n=400") %>% 
  rename(MDI = "Minimum Distance") %>% 
  group_by(estimator,p_shape,scenario,dimension) %>% 
  summarise_at(vars(MDI), list(mean = ~ mean(.),
                               sd   = ~ sd(.))) %>% 
  ungroup() %>% 
  left_join(.,quantile_rnd_matrix) %>% 
  mutate(p_shape = as.factor(p_shape),
         estimator = factor(estimator,levels = c("PML","fastICA","DCov","CvM"))) %>% 
  
  ggplot(.,aes(x = p_shape, group = estimator))+
  geom_line(aes(y = mean,color = estimator), size = 1)+
  geom_point(aes(y = mean,color = estimator, shape = estimator, fill = estimator), size = 3)+
  geom_hline(aes(yintercept = percentile_5),linetype = "dotted")+
  facet_grid(dimension~scenario)+
  scale_color_manual(values = palette,breaks = c("PML","fastICA","DCov","CvM"))+
  scale_fill_manual(values = palette , breaks = c("PML","fastICA","DCov","CvM"))+
  scale_shape_manual(values = shapes , breaks = c("PML","fastICA","DCov","CvM")) + 
  scale_x_discrete(breaks = c("0.5","1.57","2","2.43","4","52","100"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "white"))+
  labs(title = "General Assessment", y = "MDI", 
       fill = " ",color = " ",shape = " ")
different_figure
