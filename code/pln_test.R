library(tidyverse)
library(corrplot)
library(lattice)
library(GGally)
library(PLNmodels)



main_theme = theme_bw()+
  theme(line = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks =  element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size=22),
        axis.text.y = element_text(colour = "black", size=22),
        legend.title = element_text(colour = "black", size=20),
        legend.title.align=0.5,
        legend.text = element_text(colour = "black", size=18),
        axis.title=element_text(size=28))


vignette("PLNmodels")
##### Load data

read.table("data/data_clean/Metadata_grid_CAM.txt") -> metadata_grid
read.table("data/data_clean/OTU_plant_CAM.txt") -> otu_plant
read.table("data/data_clean/OTU_virus_CAM.txt") -> otu_virus

metadata_grid%>%select(pHwater, lime_tot, MO, Phos, K, Mg, Ca, Na, N, C, 
                                        CN, clay, SiltF, SiltC, SandF, SandC, Cl, Res, Cond,
                                        dep_oxy_num, Fauche, Pature, natural_landscape,wetland, cultivated, artificial, non_emitting)%>%
  scale()%>%
  as.data.frame()-> covariate_grid
metadata_grid %>%
  dplyr::select(pHwater, lime_tot, MO, Phos, K, Mg, Ca, Na, N, C, 
                CN, clay, SiltF, SiltC, SandF, SandC, Cl, Res, Cond,
                dep_oxy_num)%>%
  scale()%>%
  data.frame()-> metadata_grid_scale

######### PCA analysis of env variables
pca_env <-PCA(metadata_grid_scale, graph =F)
CAH_env = HCPC(pca_env, metric = "euclidean", method = "ward", nb.clust =3)
covariate_grid = cbind(covariate_grid, clust = CAH_env$data.clust$clust )
str(covariate_grid)
#Abundance pant
read.table("data/data_clean/abund_plant_grid.txt", header = T) -> abund_plant_grid 

grid_plant_mat <- prepare_data(abund_plant_grid,covariate_grid)
grid_plant_mat
##### PNL model
## Null model 
PLN_plant_null <- PLN(Abundance ~ 1, grid_plant_mat)
PLN_plant_null$criteria
data.frame(
  fitted   = as.vector(fitted(PLN_plant_null)),
  observed = as.vector(grid_plant_mat$Abundance)
) %>% 
  ggplot(aes(x = observed, y = fitted)) + 
  geom_point(size = .5, alpha =.25 ) + 
  scale_x_log10() + 
  scale_y_log10() + 
  annotation_logticks() +
  theme_bw() 

PLN_plant_null %>% sigma() %>% cov2cor() %>% corrplot()
PLN_plant_null%>% coef() %>% corrplot()

## Full model
PLN_plant_full <- PLN(Abundance ~ Res + non_emitting + Fauche + Pature + Phos + 
                        K + Cl + Ca + wetland + dep_oxy_num , grid_plant_mat)

data.frame(
  fitted   = as.vector(fitted(PLN_plant_full)),
  observed = as.vector(grid_plant_mat$Abundance)
) %>% 
  ggplot(aes(x = observed, y = fitted)) + 
  geom_point(size = .5, alpha =.25 ) + 
  scale_x_log10() + 
  scale_y_log10() + 
  annotation_logticks() +
  theme_bw() 

PLN_plant_full%>% coef()%>% levelplot()


PLN_plant_null$criteria
PLN_plant_full$criteria

##### PNL LDA
PLNLDA_plant_null <- PLNLDA(Abundance ~ 1, group = clust,grid_plant_mat)
data("trichoptera")
hist(as.matrix(trichoptera$Abundance))
plot(PLNLDA_plant_null)
