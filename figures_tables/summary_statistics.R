## Code to replicate Table 1 of the paper
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

coefficient_distribution <- read_csv(file = here("final_databases","coefficient_distribution.csv"))

coefficient_distribution %>% filter(scenario == "n=400") %>% 
  filter(entry %in% c("b11","b22","b12","b21")) %>% 
  ungroup() %>% 
  dplyr::select(dimension,estimator,entry,mean,sd,p_shape) %>% distinct() %>% 
  mutate(p_shape = round(p_shape,2)) %>% 
  pivot_wider(names_from = c("p_shape","dimension"),values_from = c("mean","sd")) %>% 
  dplyr::select(estimator,entry,
                "mean_0.5_k=2","sd_0.5_k=2",
                "mean_1.57_k=2","sd_1.57_k=2",
                "mean_2.43_k=2","sd_2.43_k=2",
                "mean_100_k=2","sd_100_k=2",
                "mean_0.5_k=3","sd_0.5_k=3",
                "mean_1.57_k=3","sd_1.57_k=3",
                "mean_2.43_k=3","sd_2.43_k=3",
                "mean_100_k=3","sd_100_k=3") %>% 
  arrange(entry,factor(estimator, levels = c('PML','fastICA','DCov','CvM'))) %>% 
  mutate(entry = paste0("$\\hat{b_{",str_sub(entry,2,3),"}} - b_{",str_sub(entry,2,3),"}$"),
         estimator = case_when(
           estimator == "PML" ~ "$PML$",
           T ~ paste0("$",estimator,"$")
         )) %>% 
  mutate(entry = if_else(estimator == "$PML$",entry,"")) %>% 
  print() %>%
  kable(.,format = "latex",escape = F,col.names = NULL,
        booktabs = T, digits = 3) %>% 
  add_header_above(c(" " = 2, 
                     rep(c("$mean$" = 1, "$sd$" = 1),8)),escape = F) %>% 
  add_header_above(c(" " = 2 ,rep(c("$p=0.5$" = 2, "$p=1.57$" = 2,"$p=2.47$" = 2,"$p=100$" = 2),2)),escape = F) %>% 
  add_header_above(c(" " = 2, "$k=2$" = 8, "$k=3$" = 8), escape = F)
