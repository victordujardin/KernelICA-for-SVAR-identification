## Code to replicate Figure 3 of the paper
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


## Parameters of the p-generalized distribution
## Number of shape parameter values (in the specific assessment and for inference table we only focuse on 4 distributional scenario - See paper)
seq1   <- seq(from = 0.5, to = 3.5, length.out = 15) ## Sequence of p-Generalized parameter for general evaluation exercise (1st part)
seq2   <- seq(from = 4, to = 100, length.out = 5) ## Sequence of p-Generalized parameter for general evaluation exercise (2nd part)
#SEQ     <- c(0.5,1.57,2.43,100) ## Sequence of p-Generalized parameter for statistical inference exercise
SEQ    <- as.numeric(c(seq1, seq2))

rm(seq1,seq2)


## Load Databases with simulated data
average_performance_spec <- read_csv(file = here("final_databases","average_performance_specific.csv")) #%>% 
  #dplyr::select(-"X1")


## Colors
palette <- c('CvM'     = '#66c2a5',
             'fastICA' = '#fc8d62',
             'DCov'    = '#8da0cb',
             'PML'     = '#e78ac3')

shapes <- c('CvM'     = 21,
            'fastICA' = 22,
            'DCov'    = 24,
            'PML'     = 23)


different_figure <- average_performance_spec %>% 
  filter(scenario != "n=200") %>% 
  rename(MDI = "Minimum Distance") %>% 
  group_by(estimator,p_shape,scenario,dimension) %>% 
  summarise_at(vars(MDI), list(mean = ~ mean(.),
                               sd   = ~ sd(.))) %>% 
  ungroup() %>% 
  mutate(p_shape = as.factor(p_shape),
         estimator = factor(estimator,levels = c("PML","fastICA","DCov","CvM"))) %>% 
  ggplot(.,aes(x = p_shape, group = estimator, fill = estimator))+
  geom_line(aes(y = mean,color = estimator), size = 1)+
  geom_point(aes(y = mean, shape = estimator, col = estimator), size = 3)+
  facet_grid(dimension~scenario)+
  scale_color_manual(values = palette, breaks = c("PML","fastICA","DCov","CvM"))+
  scale_fill_manual(values = palette, breaks = c("PML","fastICA","DCov","CvM"))+
  scale_shape_manual(values = shapes, breaks = c("PML","fastICA","DCov","CvM"))+
  scale_x_discrete(breaks = c("0.5","1.57","2","2.43","4","52","100"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "white"))+
  labs(title = "Specific Assessment", y = "MDI", fill = " ", color = " ", 
       shape = " ")


average_performance_spec %>% 
  mutate(scenario = str_replace(scenario,"n","T")) %>% filter(scenario == "T=400") %>% 
  rename(MDI = "Minimum Distance") %>% 
  group_by(estimator,p_shape,scenario,dimension) %>% 
  summarise_at(vars(MDI), list(mean = ~ mean(.),
                               sd   = ~ sd(.))) %>% 
  ungroup() %>% 
  mutate(p_shape = as.factor(p_shape),
         estimator = factor(estimator,levels = c("PML","fastICA","DCov","CvM"))) %>% 
  ggplot(.,aes(x = p_shape, group = estimator))+
  geom_line(aes(y = mean,color = estimator), size = 1)+
  geom_point(aes(y = mean,color = estimator, shape = estimator, fill = estimator), size = 3)+
  geom_ribbon(aes(ymin = mean - sd,ymax = mean + sd, fill = estimator), alpha = 0.3)+
  facet_grid(dimension ~ estimator)+
  scale_color_manual(values = palette, breaks = c("PML","fastICA","DCov","CvM"))+
  scale_fill_manual(values = palette, breaks = c("PML","fastICA","DCov","CvM"))+
  scale_shape_manual(values = shapes, breaks = c("PML","fastICA","DCov","CvM"))+
  scale_x_discrete(breaks = c("0.5","1.57","2","2.43","4","52","100"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "white"))+
  labs(title = "Specific Assessment", y = "MDI")
#ss

