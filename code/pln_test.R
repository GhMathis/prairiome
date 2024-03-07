library(tidyverse)
library(corrplot)
library(lattice)
library(GGally)
library(PLNmodels)
library(FactoMineR)


main_theme = theme_bw()+
  theme(line = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks =  element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size=22),
        axis.text.y = element_text(colour = "black", size=22),
        legend.title = element_text(colour = "black", size=20, hjust = 0.5),
        legend.text = element_text(colour = "black", size=18),
        axis.title=element_text(size=28))


vignette("PLNmodels")
##### Load data

read.table("data/data_clean/Metadata_grid_CAM.txt") -> metadata_grid
read.table("data/data_clean/OTU_plant_CAM.txt") -> otu_plant
read.table("data/data_clean/OTU_virus_CAM.txt") -> otu_virus

metadata_grid%>%select(Grid_code,pHwater, lime_tot, MO, Phos, K, Mg, Ca, Na, N, C, 
                                        CN, clay, SiltF, SiltC, SandF, SandC, Cl, Res, Cond,
                                        dep_oxy_num, Fauche, Pature, natural_landscape,wetland, cultivated, artificial, non_emitting)%>%
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

#Abundance plant
read.table("data/data_clean/abund_plant_grid.txt", header = T) -> abund_plant_grid 



##### PNL model quadra
otu_plant%>%
  mutate(across(where(is.numeric), ~.*20),
         Grid_code = str_extract(Host_code, ".._CAM_.."))->abund_plant_quadra

abund_plant_quadra%>%
  select(Grid_code)%>%
  full_join(covariate_grid, by = "Grid_code") -> covariate_quad
abund_plant_quadra$Grid_code
abund_plant_quadra%>%
  select(Grid_code)%>%
  left_join(cbind(as.data.frame(pca_env$ind$coord),
                          Grid_code = rownames(abund_plant_grid)), by = "Grid_code")%>%
  select(-Grid_code)-> axis_quad

abund_plant_quadra%>%
  select(-c(Grid_code, Host_code)) -> abund_plant_quadra

#quad_plant_mat <- prepare_data(abund_plant_quadra,covariate_quad)
quad_plant_mat <- prepare_data(abund_plant_quadra,axis_quad)
quad_plant_mat


PLN_plant_null <- PLN(Abundance ~ 1, quad_plant_mat)
data.frame(
  fitted   = as.vector(fitted(PLN_plant_null)+1),
  observed = as.vector(quad_plant_mat$Abundance+1)
) %>% 
  ggplot(aes(x = observed, y = fitted)) + 
  geom_point(size = .5, alpha =.25 ) + 
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks()+
  theme_bw() 

PLN_plant_null %>% sigma() %>% cov2cor() %>% corrplot()

PLN_plant_null %>% sigma() -> temp
rownames(temp) = NULL
colnames(temp) = NULL
temp%>% cov2cor()%>% corrplot(
                 is.corr = T)
PLN_plant_full %>% coef() %>%exp()%>% levelplot()

PLN_plant_null %>% coef() %>% levelplot()


PLN_plant_full <- PLN(Abundance ~ ., quad_plant_mat)
data.frame(
  fitted   = as.vector(fitted(PLN_plant_full)+1),
  observed = as.vector(quad_plant_mat$Abundance+1)
) %>% 
  ggplot(aes(x = observed, y = fitted)) + 
  geom_point(size = .5, alpha =.25 ) + 
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks()+
  theme_bw() 

PLN_plant_full %>% sigma() -> temp2
rownames(temp2) = NULL
colnames(temp2) = NULL
temp2%>% corrplot(hclust.method = "ward.D2",is.corr = F)

PLN_plant_full %>% coef() %>%exp()%>% levelplot()

PLN_plant_full %>% coef() %>% t() %>% as.data.frame()%>%
  rownames_to_column(var = "sp")%>%
  ggplot()+
  geom_point(aes(Dim.1, Dim.2))+
  geom_text_repel(aes(Dim.1, Dim.2, label = sp))
  # data.frame(
#   fitted(PLN_plant_null))%>%
#   setNames(colnames(quad_plant_mat$Abundance))%>%
#   pivot_longer(everything(), names_to = "sp", values_to =  "Ab_esti") -> pred_long
# cbind(covariate_quad, data.frame(quad_plant_mat$Abundance))%>%
#   pivot_longer(-names(covariate_quad), names_to = "sp", values_to =  "Ab_obs") -> obs_long
# cbind(obs_long, pred_long[,2])-> all_data

##### PNL model
grid_plant_mat <- prepare_data(abund_plant_grid,covariate_grid%>%column_to_rownames(Grid_code))
grid_plant_mat
## Null model 
PLN_plant_null <- PLN(Abundance ~ 1, grid_plant_mat)
PLN_plant_null$criteria
data.frame(
  fitted   = as.vector(fitted(PLN_plant_null)+1),
  observed = as.vector(grid_plant_mat$Abundance+1)
) %>% 
  ggplot(aes(x = observed, y = fitted)) + 
  geom_point(size = .5, alpha =.25 ) + 
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() + annotation_logticks()
length( as.vector(grid_plant_mat$Abundance))
plot(as.vector(cor(grid_plant_mat$Abundance)), as.vector(PLN_plant_null %>% sigma() %>% cov2cor()),
     xlab= "cor observed", ylab = "cor estimated")
str(as.vector(PLN_plant_null %>% sigma() %>% cov2cor()))


PLN_plant_null %>% sigma() %>% cov2cor() %>% corrplot()
hclust()
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
  annotation_logticks() +
  theme_bw() 

PLN_plant_full%>% coef()%>% levelplot()
PLN_plant_null%>% coef()%>% levelplot()

PLN_plant_null$criteria
PLN_plant_full$criteria

##### PNL LDA
PLNLDA_plant_null <- PLNLDA(Abundance ~ 1, group = clust,grid_plant_mat)
data("trichoptera")
hist(as.matrix(trichoptera$Abundance))
plot(PLNLDA_plant_null)
