## Point estimates
## Code to replicate Table 1 of the paper
## Load packege here for managing workflow


library(dplyr)
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


## Load database with IRFs
load(file = here("final_databases","bootstrap_irfs_ordering_frob.RData"))

## Boostrap IRFS

bootstrap.irfs <- list(fastICA = my_boot$fastICA,
                       DCov    = my_boot$DCov,
                       PML     = my_boot$PML,
                       CvM     = my_boot$CvM)
irf_db <- NULL
for (i in 1:length(bootstrap.irfs)) {
  boot <- bootstrap.irfs[[i]]$bootstrap
  name <- names(bootstrap.irfs)[i]
  
  ## Calculate the median response for each irf of interest
  median_response <- do.call(rbind,lapply(boot, function(x){
    x$irf
  })) %>% as_tibble() %>% group_by(V1) %>% 
    mutate(boot_rep = 1:n()) %>% 
    pivot_longer(2:10,names_to = "irf_name",values_to = "value") %>% 
    dplyr::rename(horizon = V1) %>% 
    group_by(horizon,irf_name) %>% #filter(irf_name == "epsilon[ GG ] %->% GG" & horizon == 1) %>% print()
    nest() %>% 
    summarise(median_response = map_dbl(data, ~ median(.$value)),
              quantiles = map_df(data, ~quantile(.$value,probs = c(0.025,0.05,0.16,0.84,0.95,0.975))),
              standard_dev = map_dbl(data, ~sd(.$value)),
              mean_response = map_dbl(data,~mean(.$value))) %>%  
    unnest(quantiles) %>% 
    group_by(horizon,irf_name) %>% 
    mutate(shock = str_remove(str_extract(irf_name,"\\[ .."),"\\[ "),
           variable = str_remove(str_extract(irf_name,"\\% .."),"\\% ")) %>% 
    ungroup() %>% 
    mutate(estimator = name)
  
  irf_db <- bind_rows(median_response, irf_db)
}

## Mixing matrix
preliminary_table <- irf_db %>% 
  filter(horizon == 1) %>% 
  dplyr::select(-irf_name) %>% 
  dplyr::select(estimator,shock,variable,everything()) %>% 
  group_by(shock,estimator) %>% 
  mutate(value_norm = median_response[horizon == 1 & shock == variable]) %>% 
  mutate_at(vars(median_response,mean_response,contains("%")), ~ ./value_norm) %>% 
  mutate(entry_table = case_when(
    sign(median_response) == sign(`2.5%`)  & sign(median_response) == sign(`97.5%`) & shock != variable
    ~ paste0(round(median_response,3),"***"),
    sign(median_response) == sign(`5%`)  & sign(median_response) == sign(`95%`) & shock != variable
    ~ paste0(round(median_response,3),"**"),
    sign(median_response) == sign(`16%`)  & sign(median_response) == sign(`84%`) & shock != variable
    ~ paste0(round(median_response,3),"*"),
    T ~ as.character(round(median_response,3))
  )) %>% 
  mutate_at(vars(variable), ~ fct_relabel(.,~c("GG" = "$G_t$", "TT" = "$Tax_t$", "XX" = "$GDP_t$"))) %>% 
  dplyr::select(estimator,shock,variable,entry_table) %>% 
  pivot_wider(values_from = c("entry_table"), names_from = c("shock")) %>% 
  ungroup()

preliminary_table %>% pivot_wider(names_from = "estimator", values_from = c("GG","TT","XX")) %>% 
  dplyr::select(variable,contains("PML"),contains("fastICA"),contains("DCov"),contains("CvM")) %>% 
  kable(format = "latex",booktabs = T,escape = F,col.names = NULL) %>% 
  add_header_above(c(" " = 1, rep(c("\\\\varepsilon_1" = 1 ,"\\\\varepsilon_2" = 1,"\\\\varepsilon_3" = 1),4)),
                     escape = F, line = F) %>% 
  add_header_above(c(" " = 1, "$PML$" = 3, "$fastICA$" = 3, "$DCov$" = 3,"$CvM$" = 3),escape = F)

bind_cols(preliminary_table %>% filter(estimator %in% c("PML","DCov")) %>% dplyr::select(-estimator),
          preliminary_table %>% filter(estimator %in% c("fastICA","CvM")) %>% 
            arrange(desc(estimator)) %>% mutate(estimator = fct_inorder(estimator)) %>%  dplyr::select(-estimator)) %>% 
  kable(format = "latex",booktabs = T,escape = F,col.names = NULL) %>% 
  add_header_above(c(" " = 1, "\\\\varepsilon_1" = 1 ,"\\\\varepsilon_2" = 1,"\\\\varepsilon_3" = 1,
                     " " = 1,"\\\\varepsilon_1" = 1 ,"\\\\varepsilon_2" = 1,"\\\\varepsilon_3" = 1),escape = F, line = F) %>% 
  add_header_above(c(" " = 1, "$PML$" = 3, " " = 1, "$fastICA$" = 3),escape = F)
  

## Normalized IRF Spending
estimators <- c("PML","fastICA","DCov","CvM")
plots_irf <- vector("list",4)
for (p in 1:length(estimators)) {
  plots_irf[[p]] <- irf_db %>% 
    dplyr::select(horizon,irf_name,estimator,shock,variable,median_response,standard_dev,mean_response,contains("%")) %>% 
    group_by(shock,estimator) %>% 
    mutate(value_norm = median_response[horizon == 1 & variable == "GG"]) %>% 
    mutate_at(vars(median_response,standard_dev,mean_response,contains("%")), ~ ./value_norm) %>% 
    pivot_longer(median_response:mean_response) %>% 
    filter(shock == "GG" & estimator == estimators[p]) %>% 
    mutate(variable = fct_relabel(variable, ~c("GG" = "Spending","TT" = "Tax", "XX" = "Output"))) %>% 
    ggplot(aes(x = horizon))+
    geom_line(data = . %>% filter(name %in% c("median_response")),
              aes(y = value)) +
    geom_line(aes(y = `16%`),linetype = 2)+
    geom_line(aes(y = `84%`),linetype = 2)+
    #geom_ribbon(aes(ymin = `5%`,ymax = `95%`),alpha = 0.6, fill = "grey")+
    geom_hline(yintercept = 0,linetype = "dotted")+
    facet_wrap(~variable, ncol = 3)+
    theme_minimal()+
    theme(panel.grid.major = element_blank(),
          #panel.grid.minor.x = element_blank(),
          legend.position = "none")+
    ylim(-0.7,3)+
    labs(x = "", y = "", title = paste0(estimators[p]))
  
}

grid.arrange(grobs = plots_irf,nrow = 2)



## Normalized IRF Tax
estimators <- c("PML","fastICA","DCov","CvM")
plots_irf <- vector("list",4)
for (p in 1:length(estimators)) {
  plots_irf[[p]] <- irf_db %>% 
    dplyr::select(horizon,irf_name,estimator,shock,variable,median_response,standard_dev,mean_response,contains("%")) %>% 
    group_by(shock,estimator) %>% 
    mutate(value_norm = median_response[horizon == 1 & variable == "TT"]) %>% 
    mutate_at(vars(median_response,standard_dev,mean_response,contains("%")), ~ ./value_norm) %>% 
    pivot_longer(median_response:mean_response) %>% 
    filter(shock == "TT" & estimator == estimators[p]) %>% 
    mutate(variable = fct_relabel(variable, ~c("GG" = "Spending","TT" = "Tax", "XX" = "Output"))) %>%
    mutate(variable = fct_relevel(variable, c("Tax","Spending","Output"))) %>% 
    ggplot(aes(x = horizon))+
    geom_line(data = . %>% filter(name %in% c("median_response")),
              aes(y = value)) +
    geom_line(aes(y = `16%`),linetype = 2)+
    geom_line(aes(y = `84%`),linetype = 2)+
    #geom_ribbon(aes(ymin = `5%`,ymax = `95%`),alpha = 0.6, fill = "grey")+
    geom_hline(yintercept = 0,linetype = "dotted")+
    facet_wrap(~variable, ncol = 3)+
    theme_minimal()+
    theme(panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          legend.position = "none")+
    ylim(-0.6,1.5)+
    labs(x = "", y = "", title = paste0(estimators[p]))
  
}

grid.arrange(grobs = plots_irf,nrow = 2)
