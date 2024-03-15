#SBM
library(alluvial)
library(sbm)

#multi analy
library(car)
library(FactoMineR)
library(factoextra)
library(PerformanceAnalytics)
library(ade4)


#general
library(tidyverse)
library(vegan)
library(GGally)
library(pals)
library(RColorBrewer)
library(doParallel)

#graph
library(igraph)
library(incidentally)

source("code/functions_modified.R")

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
read.table("data/data_clean/Metadata_quadra_CAM.txt", header = T)%>%
  filter(!str_detect(Host_code, "22_CAM_15....")) -> metadata_quadra

read.table("data/data_clean/OTU_plant_CAM.txt", header = T)%>%
  filter(!str_detect(Host_code, "22_CAM_15...."))%>%
  column_to_rownames(var = "Host_code")%>%
  mutate(across(where(is.numeric), ~ifelse(. != 0,1,0 )))%>%
  dplyr::select(where( function(x) sum(x) >1))%>%
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
  filter(!str_detect(Host_code, "22_CAM_15...."))%>%
  column_to_rownames(var = "Host_code")%>%
  dplyr::select(where( function(x) sum(x) >1))%>%
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

read.table("data/data_clean/Metadata_grid_CAM.txt", header = T)%>%
  filter(!str_detect(Grid_code,"22_CAM_15")) -> metadata_grid
str(metadata_grid)
metadata_grid%>%
  dplyr::select(pHwater, lime_tot, MO, Phos, K, Mg, Ca, Na, N, C, 
                CN, clay, SiltF, SiltC, SandF, SandC, Cl, Res, Cond,
                wetland, non_emitting, cultivated,artificial,natural_landscape, depth_oxy,Fauche, Pature)%>%
         mutate(across(c(depth_oxy,Fauche, Pature), as.factor),
           across(where(is.numeric), ~decostand(., method = "standardize")))-> metadata_grid_standard

chart.Correlation(metadata_grid_standard%>%dplyr::select(where(is.numeric),-c(N,C, Na, Cond, SandC, Mg, SiltF)), histogram=TRUE, pch=19)
FAMD_env <- FAMD(metadata_grid_standard, graph = F)

CAH_env_FAMD = HCPC(FAMD_env, metric = "euclidean", method = "ward", nb.clust =3)
x11()
plot(CAH_env_FAMD, choice = "map", title = "FAMD")

pca_env <- PCA(metadata_grid_standard%>%dplyr::select(where(is.numeric)), graph = F)

CAH_env_pca = HCPC(pca_env, metric = "euclidean", method = "ward", nb.clust =3)
x11()
plot(CAH_env_pca, choice = "map", title = "PCA")

metadata_grid_standard%>%
  dplyr::select(-c(N,C, Na, Cond, SandC, Mg, SiltF)) -> metadata_grid_standard
# 
# metadata_grid_standard%>%
#   t() -> metadata_grid_standard_t
# 
# plotMyMatrix(as.matrix(metadata_grid_standard_t), dimLabels = c("metadata", "Grid"), plotOptions = list(rowNames = T,colNames = F))+
#   theme(
#     legend.position="none",
#     axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 4, face="italic"),
#     axis.text.y = element_text(size = 4))
# 
# LBM_bin_env <-
#   as.matrix(metadata_grid_standard_t)%>%
#   estimateBipartiteSBM(
#     model = 'gaussian', 
#     dimLabels = c(row = "Env", col = "Grid"))
# 
# 
# df=as.data.frame(metadata_grid_standard_t)
# df$env=rownames(df)
# df_long<-pivot_longer(df, cols=names(df)[-ncol(df)], names_to =
#                         "grid", values_to = "values")
# 
# env_list<-rownames(metadata_grid_standard_t)
# grid_list<-colnames(metadata_grid_standard_t)
# 
# df_tax_cluster=data.frame(env=env_list,env_cluster=LBM_bin_env$memberships$Env)
# df_long<- df_long%>%left_join(df_tax_cluster, by="env")
# 
# df_long %>% arrange(env_cluster) -> df_long
# correct_env_order <- unique(df_long$env)
# 
# df_grid_cluster=data.frame(grid = grid_list,spl_cluster=LBM_bin_env$memberships$Grid)
# df_long<-df_long %>% left_join(df_grid_cluster, by="grid")
# 
# df_long %>% arrange(spl_cluster) -> df_long
# correct_grid_order <- unique(df_long$grid)
# 
# ggplot(df_long,aes(x=grid, y=env, fill=values))+
#   geom_tile()+
#   theme_minimal()+
#   theme(legend.position="none",
#         axis.text.y = element_text(vjust = 0, hjust = 1, size = 6),
#         axis.text.x = element_text(angle = 90, size = 6))+
#   scale_fill_gradient(low="white", high="Black")+
#   scale_x_discrete(limits = correct_grid_order)+
#   scale_y_discrete(limits = correct_env_order)+
#   geom_hline(yintercept=0.5+cumsum(table(LBM_bin_env$memberships$Env))[-length(table(LBM_bin_env$memberships$Env))],
#              color = "red")+
#   geom_vline(xintercept=0.5+cumsum(table(LBM_bin_env$memberships$Grid))[-length(table(LBM_bin_env$memberships$Grid))],
#              color = "red")+
#   ylab("Environnement")+
#   xlab("Grid")
# 
data.frame(Host_code = colnames(bin_plant_quadra),
           clust_smb_quad_plant = LBM_bin$memberships$Quadrats,
           clust_smb_quad_virus = LBM_bin_virus$memberships$Quadrats)%>%
  mutate( Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
  left_join(data.frame(Grid_code = metadata_grid$Grid_code,
                       clus_famd_grid_env = CAH_env_FAMD$data.clust$clust), by = "Grid_code" ) -> clust_quad_df


clust_quad_df$clust_smb_quad_plant
B <-
  as.data.frame(
    table(
      clust_quad_df$Grid_code,
      clust_quad_df$clust_smb_quad_plant,
      clust_quad_df$clust_smb_quad_virus,
      clust_quad_df$clus_famd_grid_env)
  )


colnames(B) = c( "Grid", "Quadra profiles plant",
                 "Quadra profiles Virus",
                 "Quadra profiles Env FAMD", "Freq")

w   <- which(B$Freq != 0)
B <- B[w, ]
col_clust = brewer.pal(7, "Set1")
alluvial(B[, c(4, 2, 3)],
         freq = B$Freq,
         col = case_when(B$`Quadra profiles Env FAMD` == 1 ~ col_clust[1],
                         B$`Quadra profiles Env FAMD` == 2 ~ col_clust[2],
                         .default = col_clust[7]))
# 
# lapply(metadata_grid_standard%>%t(), function(x) dist(cbind(x,x)))
# dist(metadata_grid_standard%>%t())
# MultiplexFitIndep <- estimateMultiplexSBM(list(netA, netB), dependent = FALSE,
#                                           estimOptions = list(verbosity = 0))
# 

##### CCA on plant community #####

##### Full model

FAMD_env$ind$coord[,1:3]%>%as.data.frame()%>%
  mutate(grid_id = rownames(.)) %>% slice(rep(1:n(), each = 9)) -> coord_FAMD_long

modules<-tab.disjonctif(as.factor(clust_quad_df$clust_smb_quad_plant))

sbm_clustering = only.F.alt(memberships = clust_quad_df$clust_smb_quad_plant, var1 = coord_FAMD_long[,1], 
                            var2 = coord_FAMD_long[,2],var3 = coord_FAMD_long[,3])

test = analysis.function.alt(incid = bin_plant_quadra%>%t(),memb_obs = clust_quad_df$clust_smb_quad_plant, var1 = coord_FAMD_long[,1], 
                      var2 = coord_FAMD_long[,2],var3 = coord_FAMD_long[,3],NPERM = 1000)
bin_plant_quadra%>%
  as.matrix()%>%
  plotMyMatrix( dimLabels = c("Plants", "Quadrats"), plotOptions = list(rowNames = T,colNames = F))+
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 4, face="italic"),
    axis.text.y = element_text(size = 4))

##### Randomisation exemple #####
bin_plant_quadra%>%
  curveball()%>%
  as.data.frame()%>%
  set_names(colnames(bin_plant_quadra))%>%
  mutate(sp = rownames(bin_plant_quadra))%>%
  column_to_rownames("sp")%>%
  as.matrix() -> bin_plant_quadra_rand

bin_plant_quadra_rand%>%
  plotMyMatrix( dimLabels = c("Plants", "Quadrats"), plotOptions = list(rowNames = T,colNames = F))+
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 4, face="italic"),
    axis.text.y = element_text(size = 4))

bin_plant_quadra_rand %>%
  estimateBipartiteSBM(
    model = 'bernoulli', 
    dimLabels = c(row = "Plants", col = "Quadrats")) -> LBM_bin_plant_rand
LBM_bin_plant_rand$memberships$Quadrats

##### Parallelisation #####

cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)
trials = 100
test1 <- foreach(icount(trials), .combine=cbind,
                 .packages = c("tidyverse", "sbm","incidentally")) %dopar% {
  bin_plant_quadra%>%
    curveball()%>%
    as.data.frame()%>%
    set_names(colnames(bin_plant_quadra))%>%
    mutate(sp = rownames(bin_plant_quadra))%>%
    column_to_rownames("sp")%>%
    as.matrix()%>%
    estimateBipartiteSBM(
      model = 'bernoulli', 
      dimLabels = c(row = "Plants", col = "Quadrats"),
      estimOptions  = list(plot = F)) -> LBM_bin_plant_rand
  LBM_bin_plant_rand$memberships$Quadrats
    }
stopCluster(cl)

test_env_ana = analysis.function.alt(incid = bin_plant_quadra%>%t(),memb_obs = clust_quad_df$clust_smb_quad_plant, var1 = coord_FAMD_long[,1], 
                      var2 = coord_FAMD_long[,2],var3 = coord_FAMD_long[,3],NPERM = 1000, configs = test1)

##### CCA on virus community #####


cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)
trials = 1000
metadata_quadra
permutation_virus <- foreach(icount(trials), .combine=cbind,
                 .packages = c("tidyverse", "sbm","incidentally")) %dopar% {
                   bin_virus_quadra%>%
                     curveball()%>%
                     as.data.frame()%>%
                     set_names(colnames(bin_virus_quadra))%>%
                     mutate(sp = rownames(bin_virus_quadra))%>%
                     column_to_rownames("sp")%>%
                     as.matrix()%>%
                     estimateBipartiteSBM(
                       model = 'bernoulli', 
                       dimLabels = c(row = "Virus", col = "Quadrats"),
                       estimOptions  = list(plot = F)) -> LBM_bin_rand
                   LBM_bin_rand$memberships$Quadrats
                 }
stopCluster(cl)

permutation_virus[,1:10]
test_env_ana = analysis.function.alt(incid = bin_virus_quadra%>%t(),memb_obs = as.factor(clust_quad_df$clust_smb_quad_virus),
                                     var1 = as.factor(clust_quad_df$clust_smb_quad_plant), 
                                     var2 = as.factor(clust_quad_df$clus_famd_grid_env),
                                     var3 = as.factor(metadata_quadra$Collection_date), 
                                     NPERM = 1000, configs = permutation_virus%>%as.data.frame())

test_env_ana

alluvial_df <-
  as.data.frame(
    table(
      clust_quad_df$Grid_code,
      clust_quad_df$clust_smb_quad_plant,
      clust_quad_df$clust_smb_quad_virus,
      clust_quad_df$clus_famd_grid_env,
      metadata_quadra$Collection_date)
  )


colnames(alluvial_df) = c( "Grid", "Quadra profiles plant",
                 "Quadra profiles Virus",
                 "Quadra profiles Env FAMD",
                 "Collection date", "Freq")

w   <- which(alluvial_df$Freq != 0)
alluvial_df <- alluvial_df[w, ]
col_clust = brewer.pal(7, "Set1")
alluvial(alluvial_df[, c(3,1, 5, 2, 4)],
         freq = alluvial_df$Freq,
         col = case_when(alluvial_df$`Quadra profiles Virus` == 1 ~ col_clust[1],
                         .default = col_clust[2]))
##### RDA on virus community without 0 #####
library(ROCR)
rgpd.function =function(incid_mat, side = "L"){
  svd_mat = svd(bin_virus_quadra%>%t())
  
  svd_L = svd_mat$u%*%diag(sqrt(svd_mat$d))
  svd_R = diag(sqrt(svd_mat$d)) %*% svd_mat$v
  if(side == "L"){
    return(svd_L)
  }else if (side =="R"){
    return(svd_R)
  }
}
rgpd.function(bin_virus_quadra%>%t())
svd_virus = svd(bin_virus_quadra%>%t())

svd_L_virus = svd_virus$u%*%diag(sqrt(svd_virus$d))
svd_R_virus = diag(sqrt(svd_virus$d)) %*% svd_virus$v

##### Why need to be binary ?
approx.from.svd(incid = bin_virus_quadra%>%t() ,svd.L = svd_L_virus,svd.R = svd_R_virus,nvect = 20)%>%
  t()%>%
  plotMyMatrix( dimLabels = c("Virus", "Quadrats"), plotOptions = list(rowNames = T,colNames = F))+
  theme(
    legend.position="none",
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 4, face="italic"),
    axis.text.y = element_text(size = 4))

##### RGPD approx

approx_virus_quadra = svd_L_virus[,1:50] %*% t(svd_R_virus[,1:50])

analysis.function.rda<-function(incid,var1,var2,var3,configs=NULL){#industrial processing of rda
  n<-dim(incid)[1]
  modules<-incid
  one.1<-rda(modules ~ var1, na.action=na.omit)
  one.2<-rda(modules ~ var2,  na.action=na.omit)
  one.3<-rda(modules ~ var3,   na.action=na.omit)
  two.1<-rda(modules ~ var1 + var2, na.action=na.omit)
  two.2<-rda(modules ~ var1 + var3,  na.action=na.omit)
  two.3<-rda(modules ~ var2 + var3,  na.action=na.omit)
  three<-rda(modules ~ var1 + var2 + var3, na.action=na.omit)
  one.1.alone<-rda(modules ~ var1 + Condition(var2) + Condition(var3), na.action=na.omit)
  one.2.alone<-rda(modules ~ Condition(var1) + var2 + Condition(var3), na.action=na.omit)
  one.3.alone<-rda(modules ~ Condition(var1) + Condition(var2) + var3, na.action=na.omit)
  one.1.vs.3<-rda(modules ~ var1 +  Condition(var3), na.action=na.omit)
  one.1.vs.2<-rda(modules ~ var1 +  Condition(var2), na.action=na.omit)
  one.2.vs.3<-rda(modules ~ var2 +  Condition(var3), na.action=na.omit)
  one.2.vs.1<-rda(modules ~ var2 +  Condition(var1), na.action=na.omit)
  one.3.vs.1<-rda(modules ~ var3 +  Condition(var1), na.action=na.omit)
  one.3.vs.2<-rda(modules ~ var3 +  Condition(var2), na.action=na.omit)
  xxx<-list(one.1,one.2,one.3,two.1,two.2,two.3,three,one.1.alone,one.2.alone,one.3.alone,one.1.vs.3,one.1.vs.2,one.2.vs.3,one.2.vs.1,one.3.vs.1,one.3.vs.2)
  res<-t(sapply(1:16,function(x) DSFP(anova(xxx[[x]],permutations = how(nperm=NPERM-1)))))
  F_stat<-res[,3]
  Forms<-sapply(1:16,function(x) extractTerms(xxx[[x]]))
  
  if(is.null(configs)) {
    results<-data.frame("Formulas"=Forms,"df"=res[,1],"SSq"=res[,2],"F"=res[,3],"P"=res[,4])
  }
  else {
    depth<-length(configs)
    cl <- makeCluster(detectCores()-2)
    registerDoParallel(cl)
    all.F <- foreach(x=icount(depth),.combine = cbind,
                                      .packages = c("tidyverse", "vegan","permute")) %dopar% {
                                        NPERM = 1000
                                        DSFP<-function(anovaobject){
                                          c(anovaobject$D[1],anovaobject$Variance[1],anovaobject$F[1],anovaobject$P[1])
                                        }
                                        
                                        only.F.rda<-function(memberships,var1,var2,var3){
                                          modules<-memberships
                                          one.1<- rda(modules ~ var1, na.action=na.omit)
                                          one.2<- rda(modules ~ var2, na.action=na.omit)
                                          one.3<- rda(modules ~ var3, na.action=na.omit)
                                          two.1<- rda(modules ~ var1 + var2, na.action=na.omit)
                                          two.2<- rda(modules ~ var1 + var3, na.action=na.omit)
                                          two.3<- rda(modules ~ var2 + var3, na.action=na.omit)
                                          three<- rda(modules ~ var1 + var2 + var3, na.action=na.omit)
                                          one.1.alone<- rda(modules ~ var1 + Condition(var2) + Condition(var3), na.action=na.omit)
                                          one.2.alone<- rda(modules ~ Condition(var1) + var2 + Condition(var3), na.action=na.omit)
                                          one.3.alone<- rda(modules ~ Condition(var1) + Condition(var2) + var3, na.action=na.omit)
                                          one.1.vs.3<- rda(modules ~ var1 +  Condition(var3), na.action=na.omit)
                                          one.1.vs.2<- rda(modules ~ var1 +  Condition(var2), na.action=na.omit)
                                          one.2.vs.3<- rda(modules ~ var2 +  Condition(var3), na.action=na.omit)
                                          one.2.vs.1<- rda(modules ~ var2 +  Condition(var1), na.action=na.omit)
                                          one.3.vs.1<- rda(modules ~ var3 +  Condition(var1), na.action=na.omit)
                                          one.3.vs.2<- rda(modules ~ var3 +  Condition(var2), na.action=na.omit)
                                          xxx<-list(one.1,one.2,one.3,two.1,two.2,two.3,three,one.1.alone,one.2.alone,one.3.alone,one.1.vs.3,one.1.vs.2,one.2.vs.3,one.2.vs.1,one.3.vs.1,one.3.vs.2)
                                          res<-t(sapply(1:16,function(x) DSFP(anova(xxx[[x]],permutations = how(nperm=NPERM-1)))))
                                          res[,3]
                                        }
                                        
                                        
                                        
                                        only.F.rda(configs[[x]],var1,var2,var3)
                                      }
    stopCluster(cl)
    #all.F<-sapply(1:depth,function(x) only.F.rda(configs[[x]],var1,var2,var3))
    all.ecdf<-apply(all.F,1,ecdf)
    P2<-sapply(1:16,function(x) 1-((all.ecdf[[x]])(F_stat[x])))
    results<-data.frame("Formulas"=Forms,"df"=res[,1],"SSq"=res[,2],"F"=res[,3],"P"=res[,4],"P2"=P2)
  }
  results
  
}


cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)
trials = 100
metadata_quadra
permutation_virus_rgpd <- foreach(icount(trials),
                             .packages = c("tidyverse", "sbm","incidentally")) %dopar% {
                               bin_virus_quadra%>%
                                 curveball()%>%
                                 as.data.frame()%>%
                                 set_names(colnames(bin_virus_quadra))%>%
                                 mutate(sp = rownames(bin_virus_quadra))%>%
                                 column_to_rownames("sp")%>%
                                 as.matrix()%>%
                                 t()%>%
                                 rgpd.function(side = "L") -> rgpd_virus_rand
                               rgpd_virus_rand
                             }
stopCluster(cl)
str(permutation_virus_rgpd)

clust_quad_df%>%
  mutate(Collection_date = metadata_quadra$Collection_date ) ->clust_quad_df

rda_variance = analysis.function.rda(incid = approx_virus_quadra,
                      var1 = as.factor(clust_quad_df$clust_smb_quad_plant), 
                      var2 = as.factor(clust_quad_df$clus_famd_grid_env),
                      var3 = as.factor(metadata_quadra$Collection_date), 
                      configs = permutation_virus_rgpd)

mod1 = varpart(Y = as.data.frame(approx_virus_quadra),clust_quad_df%>%select(clust_smb_quad_plant),
        clust_quad_df%>%select(clus_famd_grid_env),
        clust_quad_df%>%select(Collection_date))
plot(mod1)

mod2 = varpart(Y = clust_quad_df%>%select(clust_smb_quad_virus),
               clust_quad_df%>%select(clust_smb_quad_plant),
              clust_quad_df%>%select(clus_famd_grid_env),
              clust_quad_df%>%select(Collection_date))
str(clust_quad_df)
clust_quad_df%>%
  ggplot()+
  geom_bin_2d(aes(Collection_date,clust_smb_quad_plant))+
  geom_jitter(aes(Collection_date,clust_smb_quad_plant, col = as.factor(clust_smb_quad_virus)),width = 0.1, height = 0.1)+
  labs(x ="date quadra cluster", y = "plant quadra cluster")
hist(clust_quad_df$Collection_date)

clust_quad_df%>%
  ggplot()+
  #facet_wrap(~as.factor(clust_smb_quad_virus))+
  geom_bin_2d(aes(clus_famd_grid_env,clust_smb_quad_plant,
                  fill = after_stat(density), group),
              alpha =0.8)+
  geom_jitter(aes(clus_famd_grid_env,clust_smb_quad_plant, col = as.factor(clust_smb_quad_virus)),width = 0.1, height = 0.1,
              cex = 2, alpha =0.6)+
  labs(x ="Env quadra cluster", y = "plant quadra cluster", col = "clust  vir")+
  scale_fill_distiller(palette = "OrRd", direction = 1)+
  main_theme

##### NMI #####
library(aricode)
clust_quad_df$clust_smb_quad_virus

obs_nmi_plant  = NMI(clust_quad_df$clust_smb_quad_virus, clust_quad_df$clust_smb_quad_plant)
prem_nmi_egde = apply(permutation_virus,MARGIN = 2, function(x) NMI(x, clust_quad_df$clust_smb_quad_plant))
hist(prem_nmi_egde, breaks = 100)
abline(v = obs_nmi_plant ,col = "red")

prem_nmi_row = sapply(1:1000, function(x) NMI(sample(clust_quad_df$clust_smb_quad_virus), clust_quad_df$clust_smb_quad_plant))
hist(prem_nmi_row, breaks = 100)
abline(v = obs_nmi_plant ,col = "red")




obs_nmi_date  = NMI(clust_quad_df$clust_smb_quad_virus, clust_quad_df$Collection_date)
prem_nmi_egde2 = apply(permutation_virus,MARGIN = 2, function(x) NMI(x, clust_quad_df$Collection_date))
hist(prem_nmi_egde2, breaks = 100)
abline(v = obs_nmi_date ,col = "red")

prem_nmi_row2 = sapply(1:1000, function(x) NMI(sample(clust_quad_df$clust_smb_quad_virus), clust_quad_df$Collection_date))
hist(prem_nmi_row2, breaks = 100)
abline(v = obs_nmi_date ,col = "red")

