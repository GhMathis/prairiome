library(tidyverse)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(vegan)
library(ade4)
library(car)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(sf)
# objet pour mettre en forme les graphiques ggplot
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


##### CCA #####

##### Set up data
read.table("data/data_clean/Metadata_grid_CAM.txt", header = T) -> metadata_grid


read.table("data/data_clean/OTU_plant_CAM.txt") -> otu_plant 

read.table("data/data_clean/abund_plant_grid.txt", header = T, row.names = "Grid_code")-> abund_plant_grid 

read.table("data/data_clean/TAX_plant_CAM.txt")-> tax_plant
colnames(metadata_grid) 

metadata_grid%>%
  dplyr::select(pHwater, lime_tot, MO, Phos, K, Mg, Ca, Na, N, C, 
                CN, clay, SiltF, SiltC, SandF, SandC, Cl, Res, Cond,
                dep_oxy_num, Pature, Fauche, natural_landscape, cultivated, artificial, non_emitting,wetland)%>%
  mutate(across(c(Pature,Fauche), as.factor),
         across(where(is.numeric), ~decostand(., method = "standardize")))-> metadata_grid_standard
# select species represented at only one site. 1 is the lower abundance possible for each species across all site
abund_plant_grid %>%
  select(where( function(x) sum(x) >1))%>%
  decostand(method = "hellinger")-> abund_plant_stantdard 
table(colSums(abund_plant_stantdard)<=1)

##### Full model


cca_full = cca(abund_plant_stantdard ~ .  ,metadata_grid_standard)
cca_full$call

p_var_full = summary(cca_full)
###With colinearity
RsquareAdj(cca_full)$adj.r.squared
RsquareAdj(cca_full)$r.squared

##### Variable colinearity
## Varaible are removed one by one(remove order in the select()) until all VIF are below 10.
round(vif.cca(cca_full),2)
metadata_grid_standard%>%
  select(-c(N, Cond, SiltC, SandC, SiltF, Mg, Na, clay, artificial, cultivated, MO)) -> metadata_grid_selected
cca_selected = cca(abund_plant_stantdard ~ .  ,metadata_grid_selected)
round(vif.cca(cca_selected),2)

##### Variable selection

summary(cca_selected)

cca_ordistep <- ordiR2step(cca(abund_plant_stantdard ~ 1, metadata_grid_selected),
                           scope = formula(cca_selected),
                           direction = "both",
                           R2scope = F,
                           pstep = 10000,
                           trace = FALSE)
p_var_full3 = summary(cca_ordistep)
p_var_full3  

##### Anova
cca_ordistep$tot.chi
anova.cca(cca_ordistep, permutations = 10000)

anova_cca_ordistep = anova.cca(cca_ordistep, permutations = 10000, by ="term")

## conversion en pourcentage pour faciliter l'interprétation
str(anova_cca_ordistep)
anova_cca_ordistep$ChiSquare = (anova_cca_ordistep$ChiSquare/ cca_ordistep$tot.chi)*100

RsquareAdj(cca_ordistep)$r.square

##### extract axis for graph
pl_CCA <- ordiplot(cca_ordistep,scaling = 1)# Type 1 scaling
perc <- round(100*(summary(cca_ordistep)$cont$importance[2, 1:3]), 2)

## Continuous variables
env_var_CCA = as.data.frame(pl_CCA$biplot*attr(pl_CCA$biplot,"arrow.mul"))
env_var_CCA$type = rownames(env_var_CCA )

## Discret variables
env_var_discret_CCA = as.data.frame(pl_CCA$centroids[c(2,4),])
env_var_discret_CCA$type = rownames(env_var_discret_CCA)

## Sites
site_CCA = as.data.frame(pl_CCA$site)
site_CCA$type = rownames(site_CCA)

## Species
species_CCA = as.data.frame(pl_CCA$species)


##### Clustering

## Species
cluster2_sp_CCA = hcut(species_CCA, hc_metric = "euclidean", hc_method = "ward.D2", k =2)
cluster3_sp_CCA = hcut(species_CCA, hc_metric = "euclidean", hc_method = "ward.D2", k =3)
cluster4_sp_CCA = hcut(species_CCA, hc_metric = "euclidean", hc_method = "ward.D2", k =4)
cluster5_sp_CCA = hcut(species_CCA, hc_metric = "euclidean", hc_method = "ward.D2", k =5)
fviz_silhouette(cluster2_sp_CCA)
fviz_silhouette(cluster3_sp_CCA)
fviz_silhouette(cluster4_sp_CCA)
fviz_silhouette(cluster5_sp_CCA)

## Sites
cluster2_site_CCA = hcut(site_CCA[,1:2], hc_metric = "euclidean", hc_method = "ward.D2", k =2)
cluster3_site_CCA = hcut(site_CCA[,1:2], hc_metric = "euclidean", hc_method = "ward.D2", k =3)
cluster4_site_CCA = hcut(site_CCA[,1:2], hc_metric = "euclidean", hc_method = "ward.D2", k =4)
cluster5_site_CCA = hcut(site_CCA[,1:2], hc_metric = "euclidean", hc_method = "ward.D2", k =5)
fviz_silhouette(cluster2_site_CCA)
fviz_silhouette(cluster3_site_CCA)
fviz_silhouette(cluster4_site_CCA)
fviz_silhouette(cluster5_site_CCA)

cluster_site_df_CCA = data.frame(cluster3_site_CCA$data, clust_CCA_plant = as.factor(cluster3_site_CCA$cluster))

cluster_sp_df_CCA = data.frame(cluster3_sp_CCA$data, clust_CCA_plant = as.factor(cluster3_sp_CCA$cluster))

cluster_sp_df_CCA$Plant_species = rownames(cluster_sp_df_CCA)

cluster_sp_df_CCA%>%
  full_join(tax_plant, by = join_by("Plant_species"))-> cluster_sp_df_CCA

env_var_CCA$type
env_var_CCA$type = c("Res","Non em.", "Cl","SanF","Wetlands", "Ca", "Lime", "K", "Nat. lands")
env_var_discret_CCA$type = c("Mow", "Pasture")
ggplot()+
  geom_segment(data = env_var_CCA ,aes(x=CCA1,y = CCA2, xend =0, yend = 0), arrow = arrow(length=unit(0.20,"cm"), ends="first", type = "open"), size = 1.5, col = "gray3")+
  geom_label_repel(data = env_var_CCA ,aes(x=CCA1,y = CCA2,label = type),  cex = 5, alpha = 0.75)+
  geom_point(data = cluster_sp_df_CCA, aes(x=CCA1,y = CCA2, col = clust_CCA_plant), size = 2, stroke=1.2, shape =4)+
  geom_label_repel(data = cluster_sp_df_CCA,aes(x=CCA1,y = CCA2, col = clust_CCA_plant, label = Plant_species), cex = 2,alpha = 0.80)+
  geom_point(data = cluster_site_df_CCA, aes(x=CCA1,y = CCA2, shape = clust_CCA_plant), cex = 2.4, col = "black")+
  geom_point(data = cluster_site_df_CCA, aes(x=CCA1,y = CCA2, shape = clust_CCA_plant), cex = 2, col = "darkgray")+
  geom_segment(data = env_var_discret_CCA,aes(x=CCA1,y = CCA2, xend =0, yend = 0), arrow = arrow(length=unit(0.20,"cm"), ends="first", type = "open"), size = 1, col = "gray3")+
  geom_label_repel(data = env_var_discret_CCA,aes(x=CCA1,y = CCA2,label = type),  cex = 5, alpha = 0.75)+
  labs(x = paste0("CCA1 (", perc[1], "%)"), y = paste0("CCA2 (", perc[2], "%)"), col = "Cluster species", shape = "Cluster sites")+
  scale_color_brewer(palette ="Set1")+
  main_theme+
  theme(line = element_line())


##### SBM Plants #####
library(alluvial)
library(sbm)
library(pals)
read.table("data/data_clean/Metadata_quadra_CAM.txt", header = T) -> metadata_quadra

read.table("data/data_clean/OTU_plant_CAM.txt", header = T, row.names = "Host_code")%>%
  mutate(across(where(is.numeric), ~ifelse(. != 0,1,0 ))) %>%
  select(where( function(x) sum(x) >1))%>%
  t() -> bin_plant_quadra

plotMyMatrix(as.matrix(bin_plant_quadra), dimLabels = c("Plants", "Quadrats"), plotOptions = list(rowNames = T,colNames = F))+
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 4, face="italic"),
    axis.text.y = element_text(size = 4))

if(file.exists("Results/LBM_bin.Rdata")){
  load(file="Results/LBM_bin.Rdata")
  
}else{
  LBM_bin <-
    as.matrix(bin_plant_quadra) %>%
    estimateBipartiteSBM(
      model = 'bernoulli', 
      dimLabels = c(row = "Plants", col = "Quadrats"))
  save(LBM_bin, file="Results/LBM_bin.Rdata")
}

LBM_bin$storedModels %>% arrange(ICL)
plot(LBM_bin, dimLabels = list(row = "Plants", col= "Quadrats"))

##### Aluvial plot #####
library(alluvial)

cluster_site_df_CCA%>%
  rownames_to_column(var = "Grid_code") ->cluster_site_df_CCA

data.frame(Host_code = colnames(bin_plant_quadra),
           clust_smb_quad_plant = LBM_bin$memberships$Quadrats)%>%
  mutate( Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
  left_join(cluster_site_df_CCA%>%select(-c(CCA1,CCA2)), by = "Grid_code" ) ->alluvial_site_df

alluvial_site_df$clust_smb_quad_plant
B <-
  as.data.frame(
    table(
      alluvial_site_df$Grid_code,
      alluvial_site_df$clust_smb_quad_plant,
      alluvial_site_df$clust_CCA_plant)
    )


colnames(B) = c( "Grid", "Quadra profiles",  "CCAclust", "Freq")

w   <- which(B$Freq != 0)
B <- B[w, ]
col_clust = brewer.pal(3, "Set1")
alluvial(B[, c(1, 2, 3)],
         freq = B$Freq,
         col = case_when(B$CCAclust == 1 ~ col_clust[1],
                         B$CCAclust == 2 ~ col_clust[2],
                         .default = col_clust[3]))

metadata_quadra%>%
  select(-c(Grid_code,Host_community))%>%
  left_join(alluvial_site_df, by = "Host_code") -> metadata_quadra2

write.table(metadata_quadra2, "data/data_clean/Metadata_quadra2_CAM.txt")
str(metadata_quadra2)
B <-
  as.data.frame(
    table(
      alluvial_site_df$Grid_code,
      alluvial_site_df$Ecosystem,
      alluvial_site_df$clust_smb_quad_plant,
      alluvial_site_df$clust_CCA_plant))


colnames(B) = c( "Grid","Ecosystem", "Quadra profiles",  "CCAclust", "Freq")

w   <- which(B$Freq != 0)
B <- B[w, ]
data.frame(Ecosystem = unique(B$Ecosystem),
           Col_Ecosys = watlington(10)) -> col_ecosys
B%>%
  left_join(col_ecosys, by = "Ecosystem") -> B

alluvial(B[, c(1, 2, 3,4)],
         freq = B$Freq,
         col = B$Col_Ecosys)
##### Species

data.frame(Plant_species = rownames(bin_plant_quadra),
           memb_plant_obs = LBM_bin$memberships$Plants)%>%
  left_join(cluster_sp_df_CCA, by = "Plant_species" ) -> alluvial_sp_df
B <-
  as.data.frame(
    table(
      alluvial_sp_df$Plant_species,
      alluvial_sp_df$Plant_order,
      alluvial_sp_df$memb_plant_obs,
      alluvial_sp_df$clust_CCA_plant))


colnames(B) = c( "Plant","Plant_order", "Plant profiles",  "CCAclust", "Freq")

w   <- which(B$Freq != 0)
B <- B[w, ]
col_clust = brewer.pal(3, "Set1")

alluvial(B[, c(1, 2, 3, 4)],
         freq = B$Freq,
         col = case_when(B$CCAclust == 1 ~ col_clust[1],
                         B$CCAclust == 2 ~ col_clust[2],
                         .default = col_clust[3]))

data.frame(Plant_order = unique(B$Plant_order),
           Col_order = watlington()) -> col_order
B%>%
  left_join(col_order, by = "Plant_order") -> B
alluvial(B[, c(1, 2, 3, 4)],
         freq = B$Freq,
         col = B$Col_order)





##### SBM Virus #####

read.table("data/data_clean/OTU_virus_CAM.txt", header = T, row.names = "Host_code")%>%
  select(where( function(x) sum(x) >1))%>%
  t() -> bin_virus_quadra
str(read.table("data/data_clean/OTU_virus_CAM.txt", header = T, row.names = "Host_code"))
plotMyMatrix(as.matrix(bin_virus_quadra), dimLabels = c("Virus", "Quadrats"), plotOptions = list(rowNames = T,colNames = F))+
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 4, face="italic"),
    axis.text.y = element_text(size = 4))

if(file.exists("Results/LBM_bin_virus.Rdata")){
  load(file="Results/LBM_bin_virus.Rdata")
  
}else{
  LBM_bin_virus <-
    as.matrix(bin_virus_quadra) %>%
    estimateBipartiteSBM(
      model = 'bernoulli', 
      dimLabels = c(row = "Virus", col = "Quadrats"))
  save(LBM_bin_virus, file="Results/LBM_bin_virus.Rdata")
}

LBM_bin_virus$storedModels %>% arrange(ICL)
plot(LBM_bin_virus, dimLabels = list(row = "Virus", col= "Quadrats"))

##### Alluvile plot #####

data.frame(Host_code = colnames(bin_virus_quadra),
           memb_qdt_virus_obs = LBM_bin_virus$memberships$Quadrats)%>%
  left_join(alluvial_site_virus_df, by = "Host_code" ) ->alluvial_site_virus_df




B <-
  as.data.frame(
    table(
      alluvial_site_virus_df$Grid_code,
      alluvial_site_virus_df$clust_smb_quad_plant,
      alluvial_site_virus_df$memb_qdt_virus_obs,
      alluvial_site_virus_df$clust_CCA_plant,
      alluvial_site_virus_df$Ecosystem))


colnames(B) = c( "Grid", "Quadra profiles Plt",  "Quadra profiles Vir",
                 "CCAclust","Ecosystems", "Freq")

w   <- which(B$Freq != 0)
B <- B[w, ]
col_clust = brewer.pal(3, "Set1")
alluvial(B[, c(5, 1, 3, 5)],
         freq = B$Freq,
         col = case_when(B[,3] == 1 ~ col_clust[1],
                         #B$memb_qdt_virus_obs == 2 ~ col_clust[2],
                         .default = col_clust[3]))
alluvial_site_virus_df%>%
  left_join(metadata_quadra%>%select(-c(Grid_code)), by = "Host_code") -> alluvial_site_virus_df
B <-
  as.data.frame(
    table(
      alluvial_site_virus_df$Grid_code,
      alluvial_site_virus_df$Ecosystem,
      alluvial_site_virus_df$clust_smb_quad_plant,
      alluvial_site_virus_df$clust_CCA_plant))


colnames(B) = c( "Grid","Ecosystem", "Quadra profiles",  "CCAclust", "Freq")

w   <- which(B$Freq != 0)
B <- B[w, ]
data.frame(Ecosystem = unique(B$Ecosystem),
           Col_Ecosys = watlington(10)) -> col_ecosys
B%>%
  left_join(col_ecosys, by = "Ecosystem") -> B

alluvial(B[, c(1, 2, 3,4)],
         freq = B$Freq,
         col = B$Col_Ecosys)

alluvial_site_virus_df$clust_smb_quad_plant

##### Species

data.frame(Virus_species = rownames(bin_virus_quadra),
           memb_vir_obs = LBM_bin_virus$memberships$Virus) -> alluvial_sp_virus_df
  #left_join(cluster_sp_df_CCA, by = "Virus_species" )
B <-
  as.data.frame(
    table(
      alluvial_sp_virus_df$Virus_species,
      #alluvial_sp_virus_df$Plant_order,
      alluvial_sp_virus_df$memb_vir_obs))
      #alluvial_sp_virus_df$clust_CCA_plant))


colnames(B) = c( "Virus", "Virus profiles", "Freq")

w   <- which(B$Freq != 0)
B <- B[w, ]
col_clust = brewer.pal(3, "Set1")

alluvial(B[, c(1, 2)],
         freq = B$Freq
         # col = case_when(B$CCAclust == 1 ~ col_clust[1],
         #                 B$CCAclust == 2 ~ col_clust[2],
         #                 .default = col_clust[3])
         )

data.frame(Plant_order = unique(B$Plant_order),
           Col_order = watlington()) -> col_order
B%>%
  left_join(col_order, by = "Plant_order") -> B
alluvial(B[, c(1, 2, 3, 4)],
         freq = B$Freq,
         col = B$Col_order)

alluvial_site_virus_df%>%
  filter(memb_qdt_virus_obs == 1)%>%
  pull(Host_code) -> temp
alluvial_sp_virus_df%>%
  filter(memb_vir_obs == 1)%>%
  pull(Virus_species ) -> temp2
sort(rowSums(bin_virus_quadra[temp2, temp]))

##### PLN Virus ~ plant_community_cluster #####

library(Hmsc)
vignette("Hmsc")
str(as.data.frame(t(bin_virus_quadra)))
str(alluvial_site_df%>%
      select(clust_smb_quad_plant, Host_code, Grid_code))
length(alluvial_site_df$clust_smb_quad_plant)
quad_mat <- prepare_data(as.data.frame(t(bin_virus_quadra)),
                         alluvial_site_df%>%
                           select(clust_smb_quad_plant, Host_code, Grid_code)%>%
                           mutate(clust_smb_quad_plant = as.factor(clust_smb_quad_plant))%>%
                           column_to_rownames("Host_code")
                         )
quad_mat
str(alluvial_site_df)
alluvial_site_df%>%
  select(clust_smb_quad_plant)%>%
  as.matrix() -> covar
hmsc_plant_null <- Hmsc(Y = as.data.frame(t(bin_virus_quadra)),
                        X = covar,
                        XFormula = ~clust_smb_quad_plant,
                        distr = "probit")
nChains = 2
thin = 1 
samples = 100 
transient = 50

verbose = 0 
m = sampleMcmc(hmsc_plant_null, thin = thin, samples = samples,
                           transient = transient, nChains = nChains, nParallel = nChains, verbose = verbose)
postBeta = getPostEstimate(m, parName = "Beta") 
plotBeta(m, post = postBeta, param = "Support",
         supportLevel = 0.95, split=.4, spNamesNumbers = c(F,F))


data.frame(
  fitted   = as.vector(fitted(PLN_plant_null)+1),
  observed = as.vector(quad_mat$Abundance+1)
) %>% 
  ggplot(aes(x = observed, y = fitted)) + 
  geom_point(size = .5, alpha =.25 ) + 
  theme_bw() 

PLN_plant_null %>% sigma() %>% cov2cor() %>% corrplot()

PLN_plant_null %>% sigma() -> temp
rownames(temp) = NULL
colnames(temp) = NULL
temp%>% cov2cor()%>% corrplot(
  is.corr = T)
PLN_plant_null %>% coef() %>%exp()%>% levelplot()

PLN_plant_null %>% coef() %>% levelplot()


PLN_plant_full <- PLN(Abundance ~ ., quad_mat)
data.frame(
  fitted   = as.vector(fitted(PLN_plant_full)+1),
  observed = as.vector(quad_mat$Abundance+1)
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