library(PLNmodels)
library(tidyverse)


library(corrplot)
library(vegan)
library(proxy)
library(ade4)
library(FactoMineR)
library(factoextra)
library(phyloseq)


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

read.table("data/data_clean/Metadata_grid_CAM.txt") -> metadata_grid
read.table("data/data_clean/OTU_plant_CAM.txt") -> otu_plant
read.table("data/data_clean/OTU_virus_CAM.txt") -> otu_virus

covariate_grid = metadata_grid%>%select(pHwater, lime_tot, MO, Phos, K, Mg, Ca, Na, N, C, 
                                        CN, clay, SiltF, SiltC, SandF, SandC, Cl, Res, Cond,
                                        dep_oxy_num, Pature, Fauche, natural_landscape,wetland, cultivated, artificial, non_emitting)

metadata_grid %>%
  dplyr::select(pHwater, lime_tot, MO, Phos, K, Mg, Ca, Na, N, C, 
                CN, clay, SiltF, SiltC, SandF, SandC, Cl, Res, Cond,
                dep_oxy_num)%>%
  scale()%>%
  data.frame()-> metadata_grid_scale

######### PCA analysis of env variables
pca_env <-PCA(metadata_grid_scale, graph =F)
CAH_env = HCPC(pca_env, metric = "euclidean", method = "ward", nb.clust =3)
x11()
fviz_screeplot(pca_env)

x11()
fviz_pca_biplot(pca_env,col.var = "orange2", repel = TRUE, habillage = CAH_env$data.clust$clust,
                addEllipses=TRUE, ellipse.level=0.95)+
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1")+
  main_theme



colSums(otu_plant[-1] !=0)>1-> plant_filter
table(plant_filter) 
  
colSums(otu_virus[-1])>1-> virus_filter
table(virus_filter) 

# remove rare species and stantadization
otu_plant%>%
  select(-Host_code) -> otu_plant_standa
otu_plant_standa[,plant_filter] -> otu_plant_standa
hist(as.matrix(otu_plant_standa))

otu_virus%>%
  select(-Host_code)-> otu_virus_standa
otu_virus_standa[virus_filter] -> otu_virus_standa

##### CCA virus
str(otu_virus_standa)
virus_cca = dudi.coa(df = otu_virus_standa, scannf = FALSE, nf = 4)
s.label(virus_cca$co)
##### PCA plant
plant_pca = dudi.pca(df = otu_plant_standa,row.w = virus_cca$lw, scannf = FALSE, nf = 5)
s.corcircle(plant_pca$co, label = NULL)

##### Coinertia
coin_virus_plant <- coinertia(virus_cca,plant_pca,scannf=F,nf=4) #perfom the coinertia
coin_virus_plant
coin_virus_plant$RV
plot(coin_virus_plant)

##### PCoA
virus_pcoa =dudi.pca(d =otu_virus_standa, scannf = FALSE, nf = 3) # proxy::dist(otu_virus_standa, method = "Jaccard"),
virus_pcoa$lw
scatter(virus_pcoa)
s.label(virus_pcoa$li)

##### PCA plant
plant_pcoa = dudi.pca(d = otu_plant_standa, scannf = FALSE, nf = 3)
plant_pcoa$lw
s.corcircle(plant_pca$co, label = NULL)
scatter(plant_pcoa)
plant_pcoa$lw
##### Coinertia
coin_virus_plant <- coinertia(virus_pcoa,plant_pcoa,scannf=F,nf=3) #perfom the coinerti
coin_virus_plant$RV
randtest(x=coin_virus_plant,nrepet=1000)
coin_virus_plant

plot(coin_virus_plant)

##### grid coinertia
# Plant
otu_plant%>%
  mutate( Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
  group_by(Grid_code) %>%
  mutate(across(where(is.numeric), ~ifelse(. != 0,1,0 )))%>%
  summarise_if(is.numeric, sum)%>%
  select(-Grid_code)%>%
  decostand(method = "hellinger")-> otu_plant_grid_standa
otu_plant_grid_standa = otu_plant_grid_standa[plant_filter]

#Virus
otu_virus%>%
  mutate( Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
  group_by(Grid_code) %>%
  mutate(across(where(is.numeric), ~ifelse(. != 0,1,0 )))%>%
  summarise_if(is.numeric, sum)%>%
  select(-Grid_code)%>%
  decostand(method = "hellinger")-> otu_virus_grid_standa
otu_virus_grid_standa = otu_virus_grid_standa[virus_filter]

##### PCA virus
virus_pca_grid =dudi.pca(d =otu_virus_grid_standa, scannf = FALSE, nf = 3) # proxy::dist(otu_virus_standa, method = "Jaccard"),
virus_pca_grid$lw
s.corcircle(virus_pca_grid$co, label = NULL)
scatter(virus_pca_grid)
s.label(virus_pca_grid$li)

##### PCA plant
plant_pca_grid = dudi.pca(d = otu_plant_grid_standa, scannf = FALSE, nf = 3)
plant_pca_grid$lw
s.corcircle(plant_pca_grid$co, label = NULL)
scatter(plant_pca_grid)
s.label(plant_pca_grid$li)

##### Coinertia
coin_virus_plant_grid <- coinertia(virus_pca_grid, plant_pca_grid,scannf=F,nf=3) #perfom the coinerti
coin_virus_plant_grid$RV
coinrand_test = randtest(x=coin_virus_plant_grid,nrepet=10000)
plot(coinrand_test, nclass = 10, coeff = 1)
coin_virus_plant_grid

plot(coin_virus_plant_grid)

vir_coord = coin_virus_plant_grid$mX[,1:2]
colnames(vir_coord) = c("virX", "virY")
plant_coord = coin_virus_plant_grid$mY[,1:2]
colnames(plant_coord) = c("pltX", "pltY")

df_coine = cbind(cbind(plant_coord, vir_coord), clust = CAH_env$data.clust$clust)
ggplot(df_coine)+
  geom_point(aes(virX, virY))+
  geom_segment(data = df_coine,aes(x=pltX,y = pltY, xend =virX, yend = virY, col = clust), arrow = arrow(length=unit(0.30,"cm"), ends="first",
                                                                                                         type = "open"), linewidth = 1)+
  main_theme
##### PCoA virus

ggplot(df_coine)+
  geom_point(aes(virX, virY))
ggplot(df_coine)+
  geom_point(aes(pltX, pltY))
otu_virus%>%
  mutate( Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
  group_by(Grid_code) %>%
  mutate(across(where(is.numeric), ~ifelse(. != 0,1,0 )))%>%
  summarise_if(is.numeric, sum)%>%
  select(-Grid_code)%>%
  proxy::dist( method = "Jaccard") -> jaccard_grid_virus
virus_pcoa_grid =dudi.pco(d =jaccard_grid_virus, scannf = FALSE, nf = 3) # proxy::dist(otu_virus_standa, method = "Jaccard"),
virus_pcoa_grid$lw
scatter(virus_pcoa_grid)
s.label(virus_pcoa_grid$li)

##### PCoA plant
otu_plant%>%
  mutate( Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
  group_by(Grid_code) %>%
  mutate(across(where(is.numeric), ~ifelse(. != 0,1,0 )))%>%
  summarise_if(is.numeric, sum)%>%
  select(-Grid_code)%>%
  proxy::dist( method = "Jaccard") -> jaccard_grid_plant

plant_pcoa_grid = dudi.pco(d = jaccard_grid_plant, scannf = FALSE, nf = 3)
plant_pcoa_grid$lw
scatter(plant_pcoa_grid)
s.label(plant_pcoa_grid$li)

##### Coinertia
coin_virus_plant_grid <- coinertia(virus_pcoa_grid, plant_pcoa_grid,scannf=F,nf=3) #perfom the coinerti
plot(coin_virus_plant_grid)
coin_virus_plant_grid$RV
coinrand_test = randtest(x=coin_virus_plant_grid,nrepet=1000)
plot(coinrand_test, nclass = 10, coeff = 1)
coin_virus_plant_grid

###### Bray curtis
##### PCoA virus
otu_virus%>%
  mutate( Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
  group_by(Grid_code) %>%
  mutate(across(where(is.numeric), ~ifelse(. != 0,1,0 )))%>%
  summarise_if(is.numeric, sum)%>%
  select(-Grid_code)%>%
  proxy::dist( method = "bray") -> bray_grid_virus
virus_pcoa_grid =dudi.pco(d =bray_grid_virus, scannf = FALSE, nf = 3) # proxy::dist(otu_virus_standa, method = "Jaccard"),
virus_pcoa_grid$eig
scatter(virus_pcoa_grid)
s.label(virus_pcoa_grid$li)


##### PCoA plant
otu_plant%>%
  mutate( Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
  group_by(Grid_code) %>%
  mutate(across(where(is.numeric), ~ifelse(. != 0,1,0 )))%>%
  summarise_if(is.numeric, sum)%>%
  select(-Grid_code)%>%
  proxy::dist( method = "bray")-> bray_grid_plant

plant_pcoa_grid = dudi.pco(d = bray_grid_plant, scannf = FALSE, nf = 3)
plant_pcoa_grid$eig
scatter(plant_pcoa_grid)
s.label(plant_pcoa_grid$li)

##### Coinertia
coin_virus_plant_grid <- coinertia(virus_pcoa_grid, plant_pcoa_grid,scannf=F,nf=3) #perfom the coinerti
plot(coin_virus_plant_grid)
coin_virus_plant_grid$RV
coinrand_test = randtest(x=coin_virus_plant_grid,nrepet=1000)
plot(coinrand_test, nclass = 10, coeff = 1)
coin_virus_plant_grid

vir_coord = coin_virus_plant_grid$mX[,1:2]
colnames(vir_coord) = c("virX", "virY")
plant_coord = coin_virus_plant_grid$mY[,1:2]
colnames(plant_coord) = c("pltX", "pltY")

df_coine = cbind(cbind(plant_coord, vir_coord), clust = CAH_env$data.clust$clust)
ggplot(df_coine)+
  geom_point(aes(virX, virY))+
  geom_segment(data = df_coine,aes(x=pltX,y = pltY, xend =virX, yend = virY, col = clust), arrow = arrow(length=unit(0.30,"cm"), ends="first",
                                                                                            type = "open"), linewidth = 1)+
  main_theme

##### Load data
read.table("data/data_clean/Metadata_Grid_CAM.txt", header = T) -> metadata_grid
vignette("PNL")
metadata_grid %>%
  dplyr::select(pHwater, lime_tot, MO, Phos, K, Mg, Ca, Na, N, C, 
                CN, clay, SiltF, SiltC, SandF, SandC, Cl, Res, Cond,
                dep_oxy_num)%>%
  scale()%>%
  data.frame()-> metadata_grid_scale
metadata_grid
cor(metadata_grid_scale)%>%
  corrplot( method = "number")
PerformanceAnalytics::chart.Correlation(metadata_grid_scale, histogram=TRUE, pch=19,
                                        labels = names(metadata_grid))



X = data.frame(x1 = c(2,3,3), x2 = c(2, 0, 1) ) 
Y = data.frame(y1 = c(1,0,1), y2 = c(1,1,0) ) 
test1 = dudi.pca(X, scannf = FALSE, nf = 2)
test2 =  dudi.pca(Y, scannf = FALSE, nf = 2)
s.corcircle(test1$co)
s.label(test1$li)
s.label(test2$li)
biplot(test1)
biplot(test2)
co_test = coinertia(test1, test2,scannf=F,nf=2)
plot(co_test)
biplot(co_test)
test3 =  dudi.pca(cbind(X,Y),scannf = FALSE, nf = 2)
biplot(test3)
s.corcircle(test3$co)
s.label(test3$li)
dev.off()
plot(1:100, log(dpois(1:100, 11)), "l")
rand_pt = rpois(100, 11)
plot(seq(1,10,length.out = 100), log(rand_pt))
points(seq(1,10,length.out = 100), rand_pt, col="red")

