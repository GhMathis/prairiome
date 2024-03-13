#SBM
library(alluvial)
library(sbm)

#general
library(tidyverse)
library(vegan)
library(GGally)
library(pals)
library(RColorBrewer)

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
##### SBM Plants #####
read.table("data/data_clean/Metadata_quadra_CAM.txt", header = T) -> metadata_quadra

read.table("data/data_clean/OTU_plant_CAM.txt", header = T)%>%
  filter(!str_detect("22_CAM_15....", Host_code))%>%
  column_to_rownames(var = "Host_code")%>%
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


##### SBM Virus #####

read.table("data/data_clean/OTU_virus_CAM.txt", header = T)%>%
  filter(!str_detect("22_CAM_15....", Host_code))%>%
  column_to_rownames(var = "Host_code")%>%
  select(where( function(x) sum(x) >1))%>%
  t() -> bin_virus_quadra

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


##### SBM environement #####
library(car)
library(FactoMineR)
library(factoextra)
read.table("data/data_clean/Metadata_grid_CAM.txt", header = T)%>%
  filter(!str_detect("22_CAM_15", Grid_code)) -> metadata_grid

metadata_grid%>%
  dplyr::select(pHwater, lime_tot, MO, Phos, K, Mg, Ca, Na, N, C, 
                CN, clay, SiltF, SiltC, SandF, SandC, Cl, Res, Cond,
                dep_oxy_num)%>%
         mutate(across(where(is.numeric), ~decostand(., method = "standardize")))-> metadata_grid_standard

pca_env <-PCA(metadata_grid_standard,graph = F)
CAH_env = HCPC(pca_env, metric = "euclidean", method = "ward", nb.clust =3)
x11()
fviz_screeplot(pca_env)

x11()
fviz_pca_biplot(pca_env,col.var = "orange2", repel = TRUE, habillage = CAH_env$data.clust$clust,
                addEllipses=TRUE, ellipse.level=0.95)+
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1")+
  main_theme



ggpairs(metadata_grid_standard%>%as.matrix()%>%as.data.frame()) 

metadata_grid_standard%>%
  select(-c(Na, Cond, MO,SandC, C)) -> metadata_grid_standard


metadata_grid_standard%>%
  t() -> metadata_grid_standard

plotMyMatrix(as.matrix(metadata_grid_standard), dimLabels = c("metadata", "Grid"), plotOptions = list(rowNames = T,colNames = F))+
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 4, face="italic"),
    axis.text.y = element_text(size = 4))

LBM_bin_env <-
  as.matrix(metadata_grid_standard)%>%
  estimateBipartiteSBM(
    model = 'gaussian', 
    dimLabels = c(row = "Env", col = "Grid"))


df=as.data.frame(metadata_grid_standard)
df$env=rownames(df)
df_long<-pivot_longer(df, cols=names(df)[-ncol(df)], names_to =
                        "grid", values_to = "values")

env_list<-rownames(metadata_grid_standard)
grid_list<-colnames(metadata_grid_standard)

df_tax_cluster=data.frame(env=env_list,env_cluster=LBM_bin_env$memberships$Env)
df_long<- df_long%>%left_join(df_tax_cluster, by="env")

df_long %>% arrange(env_cluster) -> df_long
correct_env_order <- unique(df_long$env)

df_grid_cluster=data.frame(grid = grid_list,spl_cluster=LBM_bin_env$memberships$Grid)
df_long<-df_long %>% left_join(df_grid_cluster, by="grid")

df_long %>% arrange(spl_cluster) -> df_long
correct_grid_order <- unique(df_long$grid)

ggplot(df_long,aes(x=grid, y=env, fill=values))+
  geom_tile()+
  theme_minimal()+
  theme(legend.position="none",
        axis.text.y = element_text(vjust = 0, hjust = 1, size = 6),
        axis.text.x = element_text(angle = 90, size = 6))+
  scale_fill_gradient(low="white", high="Black")+
  scale_x_discrete(limits = correct_grid_order)+
  scale_y_discrete(limits = correct_env_order)+
  geom_hline(yintercept=0.5+cumsum(table(LBM_bin_env$memberships$Env))[-length(table(LBM_bin_env$memberships$Env))],
             color = "red")+
  geom_vline(xintercept=0.5+cumsum(table(LBM_bin_env$memberships$Grid))[-length(table(LBM_bin_env$memberships$Grid))],
             color = "red")+
  ylab("Environnement")+
  xlab("Grid")

data.frame(Host_code = colnames(bin_plant_quadra),
           clust_smb_quad_plant = LBM_bin$memberships$Quadrats,
           clust_smb_quad_virus = LBM_bin_virus$memberships$Quadrats)%>%
  mutate( Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
  left_join(data.frame(Grid_code = metadata_grid$Grid_code,
                       clust_smb_grid_env = LBM_bin_env$memberships$Grid,
                       clus_pca_grid_env = CAH_env$data.clust$clust), by = "Grid_code" ) -> alluvial_site_df


alluvial_site_df$clust_smb_quad_plant
B <-
  as.data.frame(
    table(
      alluvial_site_df$Grid_code,
      alluvial_site_df$clust_smb_quad_plant,
      alluvial_site_df$clust_smb_quad_virus,
      alluvial_site_df$clust_smb_grid_env,
      alluvial_site_df$clus_pca_grid_env)
  )


colnames(B) = c( "Grid", "Quadra profiles plant",
                 "Quadra profiles Virus",
                 "Quadra profiles Env",
                 "Quadra profiles Env PCA", "Freq")

w   <- which(B$Freq != 0)
B <- B[w, ]
col_clust = brewer.pal(7, "Set1")
alluvial(B[, c(5, 2, 3)],
         freq = B$Freq,
         col = case_when(B$`Quadra profiles Env PCA` == 1 ~ col_clust[1],
                         B$`Quadra profiles Env PCA` == 2 ~ col_clust[2],
                         .default = col_clust[7]))
# 
# lapply(metadata_grid_standard%>%t(), function(x) dist(cbind(x,x)))
# dist(metadata_grid_standard%>%t())
# MultiplexFitIndep <- estimateMultiplexSBM(list(netA, netB), dependent = FALSE,
#                                           estimOptions = list(verbosity = 0))
# 


