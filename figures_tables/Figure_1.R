## Code to replicate Figure 1 of the paper
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

source(here("functions_to_load.R"),local = T)
n <- 500
eps05 <- stdEpsrnd(type = "general",n = n,p = 0.5) %>% 
  as.data.frame() %>% 
  mutate(p = 0.5)

eps1  <- stdEpsrnd(type = "general",n = n,p = 1) %>% 
  as.data.frame() %>% 
  mutate(p = 1)

eps2  <- stdEpsrnd(type = "general",n = n,p = 2) %>% 
  as.data.frame() %>% 
  mutate(p = 2)

eps4  <- stdEpsrnd(type = "general",n = n,p = 4) %>% 
  as.data.frame() %>% 
  mutate(p = 4)

eps100  <- stdEpsrnd(type = "general",n = n,p = 100) %>% 
  as.data.frame() %>% 
  mutate(p = 100)

eps <- bind_rows(eps1,eps05,eps100,eps2,eps4) %>%
  arrange(p) %>% 
  mutate(p = as.character(p)) %>% 
  mutate(p = fct_inorder(p)) %>% 
  rename_all(~c("value","shape parameter"))

eps %>%  filter(value < 5 & value > -5) %>% 
  ggplot(.,aes(x = value, group = `shape parameter`)) +
  geom_density(aes(fill = `shape parameter`,col = `shape parameter`),adjust = 1.5,alpha = 0.8)+
  facet_grid(.~`shape parameter`)+
  scale_fill_brewer(palette = "Set2")+
  scale_color_brewer(palette = "Set2")+
  theme_bw()+
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())+
  labs(x="")
















 