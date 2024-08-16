## Inference of elements of mixing matrix
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

table_inference <- read_csv(here("final_databases","table_inference.csv"))

## Dimension ####
df_graph <- table_inference %>% 
  pivot_longer(4:23) %>%
  separate(col = "name",into = c("p_shape","nominal_size"),sep = "_") %>% 
  mutate(nominal_size = as.numeric(str_remove(nominal_size,"\\%"))/100) %>%
  mutate(size_distortion = value / nominal_size,
         nominal_size = as.factor(nominal_size)) %>%
  mutate(dimension = paste0(str_sub(dimension,1,1), " = ",str_sub(dimension,3,3)),
         p_shape = as.factor(p_shape),
         p_shape = fct_relevel(p_shape,"p=0.5","p=1.57","p=2.43","p=100"))


p1 <- df_graph %>% 
  mutate(estimator = case_when(
    estimator == "Distance Covariance" ~ "DCov",
    T ~ estimator
  )) %>% 
  arrange(desc(estimator)) %>% 
  mutate(estimator = fct_inorder(estimator)) %>% 
  ggplot(aes(x = nominal_size))+
  geom_boxplot(aes(x = nominal_size,y = size_distortion, 
                   group = interaction(dimension,nominal_size),
                   col = dimension,fill = dimension), alpha = 0.2)+
  geom_hline(yintercept = 1, linetype = "dashed")+
  facet_grid(estimator~p_shape)+
  coord_flip()+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill = "white"))+
  labs(y = "Size distortion", x = "", fill = "Dimension", col = "Dimension")
p1


## P-shape parameter ####
df_graph <- table_inference %>% 
  pivot_longer(4:23) %>%
  separate(col = "name",into = c("p_shape","nominal_size"),sep = "_") %>% 
  mutate(nominal_size = as.numeric(str_remove(nominal_size,"\\%"))/100) %>% 
  mutate(size_distortion = value / nominal_size,
         nominal_size = as.factor(nominal_size)) %>% 
  mutate(dimension = paste0(str_sub(dimension,1,1), " = ",str_sub(dimension,3,3)))


p2 <- df_graph %>% 
  mutate(p_shape = as.numeric(str_remove(str_extract(p_shape,"\\=.*"),"=")),
         p_shape = as.factor(p_shape)) %>% 
  ggplot(aes(x = nominal_size))+
  geom_boxplot(aes(x = nominal_size,y = size_distortion, 
                   group = interaction(p_shape,nominal_size),
                   col = p_shape,fill = p_shape), alpha = 0.2,
               position = position_dodge(width = 1))+
  geom_hline(yintercept = 1, linetype = "dashed")+
  facet_grid(~estimator)+
  theme_bw()+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_rect(fill = "white"))+
  labs(y = "Size distortion", x = "Nominal size", fill = "p-shape", col = "p-shape")

p2


grid.arrange(p1,p2,nrow =2)
