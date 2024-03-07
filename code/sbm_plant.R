library(sbm)
library(tidyverse)

#SBM----
read.table("data/data_clean/OTU_plant_CAM.txt") -> otu_plant
otu_plant%>%
  mutate(across(where(is.numeric), ~.*2),
         Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
  select(-c(Host_code, Grid_code)) %>%
  select(where( function(x) sum(x) >10))->abund_plant_quadra
table(colSums(abund_plant_quadra) ==1)
plotMyMatrix(as.matrix(abund_plant_quadra), dimLabels = c("Quadrats", "Plants"), plotOptions = list(rowNames = T,colNames = F))+
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 4, face="italic"),
    axis.text.y = element_text(size = 4))

LBM <-
  as.matrix(log(abund_plant_quadra+1)) %>%
  estimateBipartiteSBM(
    model = 'gaussian', 
    dimLabels = c(row = "Quadrats", col = "Plants"))
str(LBM_bernoulli)
# save(LBM_bernoulli, ps_rare, file = "results/LBM_bernouilli_24000.Rdata")
# ggsave("results/LBM_bernouilli_24000.svg")


abund_plant_grid %>%
  select(where( function(x) sum(x) >1))-> abund_plant_stantdard

N_SEUIL = 0
1*(abund_plant_quadra>N_SEUIL)%>%
  as.data.frame()%>%
  select(where( function(x) sum(x) >1))%>%
  t()->data_bin


table(rowSums(data_bin) == 1)

plotMyMatrix(as.matrix(data_bin), dimLabels = c("Plants", "Quadrats"), plotOptions = list(rowNames = T,colNames = F))+
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 4, face="italic"),
    axis.text.y = element_text(size = 4))

LBM_bin <-
  as.matrix(data_bin) %>%
  estimateBipartiteSBM(
    model = 'bernoulli', 
    dimLabels = c(row = "Plants", col = "Quadrats"))

str(LBM_bin)


plot(LBM_bin, dimLabels = list(row = "Plants", col= "Quadrats"))

save(LBM_bin,abund_plant_stantdard, file="Results/LBM_bin.Rdata")

LBM_bin$storedModels %>% arrange(ICL)

memb_qdt_obs <- LBM_bin$memberships$Quadrats


memb_plt_obs <- LBM_bin$memberships$Plants



ggsave("results/one_binary_matrix.svg")

library(alluvial)
cluster_sp_df_CCA
tit <- as.data.frame(Titanic)

# 2d
tit2d <- aggregate( Freq ~ Class + Survived, data=tit, sum)
alluvial( tit2d[,1:2], freq=tit2d$Freq, xw=0.0, alpha=0.8,
          gap.width=0.1, col= "steelblue", border="white",
          layer = tit2d$Survived != "Yes" )

#Alluvial plot----
library(alluvial)
cluster_site_df_CAA%>%
  rownames_to_column(var = "Grid_code") ->cluster_site_df_CAA
otu_plant%>%
  mutate( Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
            pull(Grid_code)->grid_names_vect
otu_plant%>%
  mutate( Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
  left_join(cluster_site_df_CAA, by = "Grid_code" )%>%
  pull(clust)->grid_clustCCA_vect
B <-
  as.data.frame(
    table(
      grid_names_vect,
      memb_qdt_obs,
      grid_clustCCA_vect))
colnames(B) = c( "Grid", "Plant profiles",  "CCA", "Freq")

w   <- which(B$Freq != 0)
B <- B[w, ]
alluvial(B[, c(1, 2, 3)], freq = B$Freq)

alluvial(B[, c(1, 2, 3)],
         freq = B$Freq,
         col = case_when(B$CCA == 1 ~ "red",
                         B$CCA == 2 ~ "blue",
                         .default = "darkgreen"))

B <-
  as.data.frame(
    table(
      sample_data(ps_bin)$fruit,
      sample_data(ps_bin)$locality,
      memb_spl_obs))
colnames(B) = c( "Fruit", "Locality", "Gut bacteria", "Freq")

w   <- which(B$Freq != 0)
B <- B[w, ]
alluvial(B[, c(1, 2, 3)], freq = B$Freq)
alluvial(B[, c(1, 3, 2)], freq = B$Freq)

cbind(sp = rownames(data_bin), memb_plt_obs) %>% 
  as.data.frame() %>% 
  mutate(memb_plt_obs=as.numeric(memb_plt_obs)) %>% 
  arrange(memb_plt_obs)

