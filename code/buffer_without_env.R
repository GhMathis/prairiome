##### Pacakges
# General
library(RColorBrewer)
library(tidyverse)
library(ggrepel)
library(readxl)
#library(open)
# Multivar
library(FactoMineR)
library(factoextra)
library(tidyverse)
library(corrplot)
library(ade4)
library(vegan)

library(PerformanceAnalytics)

# Geographical
library(sp)
library(sf)
library(siland)
library(terra)

main_theme = theme_bw()+
  theme(line = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks =  element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size=15),
        axis.text.y = element_text(colour = "black", size=15),
        axis.title.x = element_text(colour = "black", size=18),
        axis.title.y = element_text(colour = "black", size=18),
        legend.title = element_text(colour = "black", size=20),
        legend.title.align=0.5,
        legend.text = element_text(colour = "black", size=18),
        axis.title=element_text(size=10),
        strip.background = element_rect(fill="cornsilk"))

##### Load data
read.table("data/data_clean/Metadata_Grid_CAM.txt", header = T) -> metadata_grid

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
######### PCA analysis of env variables
pca_env <-PCA(metadata_grid_scale, graph =F)
str(pca_env)
rowSums(pca_env$var$contrib/colSums(pca_env$var$contrib)*100)
as.data.frame(round(pca_env$var$contrib,1))[,1:3]%>%
  arrange(desc(Dim.1))
rowSums(pca_env$var$cos2)
fviz_pca_var(pca_env, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
fviz_pca_ind(pca_env, col.ind = "cos2",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
round(pca_env$eig,1)
# clustering
CAH_env = HCPC(pca_env, metric = "euclidean", method = "ward", nb.clust =3)
x11()
fviz_screeplot(pca_env)

x11()
fviz_pca_biplot(pca_env,col.var = "orange2", repel = TRUE, habillage = CAH_env$data.clust$clust,
                addEllipses=TRUE, ellipse.level=0.95)+
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1")+
  main_theme

x11()
fviz_pca_biplot(pca_env,col.var = "orange2", repel = TRUE,axes = c(1,3), habillage = CAH_env$data.clust$clust,
                addEllipses=TRUE, ellipse.level=0.95)+
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1")+
  main_theme
dev.set(0)

## test 

pca_env$ind$coord
lamb = diag(pca_env$svd$vs, nrow = nrow(pca_env$svd$U))[,1:5]


pca_env$svd$U %*% diag(pca_env$svd$vs[1:5])
pca_env$var$coord


##### Compute buffer size

st_read("data/shapefiles/sql_statement_d551600.shp")  %>% 
  #rename(lib = lib4_16)%>%
  # mutate(lib1_16 = str_replace_all(lib1_16, "[,/ ']","_"),
  #        lib1_16 = str_replace_all(lib1_16, "[.]","_"),
  #        lib1_16 = str_replace_all(lib1_16, "[êéè]","e"),
  #        lib1_16 = str_remove(lib1_16, "[*]$"),
  #        lib1_16 = str_remove(lib1_16, "_$"))%>%
  mutate(lib = case_when(
    lib4_16 %in% c("Bâti diffus", "Bâti individuel dense", "Tissu urbain continu", "Bâti individuel lâche",            
                    "Bâti isolé", "Décharge",  "Bâti collectif",  "Terrain vague en zone urbaine",     
                    "Équipement sportif et de loisirs", "Bâti léger ou informel", "Zone d'activité et équipement","Extraction de matériaux",           
                    "Chantier", "Place", "Jardins familiaux", "Espace vert urbain",                
                    "Zone portuaire", "Bâti individuel dans parc paysager", "Cimetière" ) ~ "artificial",
    
    lib4_16 %in% c("Marais ouvert", "Feuillu", "Formation arbustive et arborée semi-ouverte", "Ripisylve",
                   "Conifère",  "Formation arbustive et arborée fermée","Plage",
                   "Etang et/ou lagune", "Canal", "Forêt mélangée") ~ "non_emitting",
    
    lib4_16 %in% c("Riz", "Luzerne", "Prairie temporaire", "Blé", "Tournesol", "Verger, oliveraie", "Friche récente",    
                   "Culture maraichère", "Sorgho, soja", "Vignoble", "Colza", "Maïs" ) ~ "cultivated",
    
    lib4_16 %in% c("Prairie naturelle","Coussoul", "Dune embryonnaire",                         
                   "Dune végétalisée", "Dune à végétation arbustive") ~ "natural_landscape",
    
    lib4_16 %in% c("Roselière", "Sansouire basse", "Marais ouvert", "Sansouire haute",                  
                  "Jonchaie", "Autre marais à végétation émergée", "Marais à marisque", "Sol nu",                           
                  "Lagune de pré-concentration", "Table saunante", "Friche salicole récente" ) ~ "wetland",
    
    
    .default = lib4_16
  ))%>%
  dplyr::select(lib)%>%
  st_transform( 2154)%>%
  
  mutate(area =  st_area(.), #area of each polygone
         n =1,
         lib2 = as.factor(lib))%>%
  pivot_wider(names_from = lib2, values_from = n, values_fill = 0)-> soil_occu


##### Change positions of grids to lambert projection (for siland package)
metadata_grid%>%
  st_as_sf(coords = c("X", "Y"), crs = 4326)%>%
  st_transform(grid_pos_metadata, crs = st_crs(soil_occu))%>%
  mutate(X = as.numeric(st_coordinates(geometry)[, 1]),
         Y =  as.numeric(st_coordinates(geometry)[, 2]))%>%
  st_join(.,
          soil_occu%>%
            dplyr::select(lib,area)
  )%>%
  distinct(Grid_code, .keep_all = TRUE)%>% # if duplication is append in the st_joit above... just in case
  as.data.frame(.)%>%
  bind_cols(as.data.frame(pca_env$ind$coord))-> metadata_grid

metadata_grid%>%
  #dplyr::select(X,Y)%>%
  st_as_sf(coords = c("X", "Y"), crs = st_crs(soil_occu))%>%
  st_buffer(4000)%>%
  st_union()%>%
  st_cast(to = "POLYGON") ->limit_map
  
soil_occu_crop = st_crop(x = soil_occu, y =limit_map)
soil_occu_crop%>%
  mutate(lib = factor(lib,levels = c("artificial", "cultivated", "wetland", "non_emitting", "natural_landscape"),
                      labels =  c("artificial", "cultivated", "wetland", "non emitting", "natural landscape")))%>%
ggplot()+
  geom_sf(col = "gray",aes(fill = lib) )+
  geom_point(data = metadata_grid,aes(X,Y, col = CAH_env$data.clust$clust),  cex = 2.5)+
  geom_text_repel(data = metadata_grid,aes(X,Y, col = CAH_env$data.clust$clust, label = 1:42),  cex = 3.5)+
  scale_fill_brewer(palette = "Pastel2")+
  scale_color_brewer(palette = "Set1")+
  labs(colour = "Grid cluster", fill ="Cover type")+
  main_theme

##### Area per class cover on the overall landscape
as.data.frame(soil_occu_crop)%>%
  dplyr::select(lib, area)%>%
  group_by(lib)%>%
  summarise(sum_area = sum(area, na.rm = T))%>%
  mutate(perc_area = as.numeric(sum_area/sum(sum_area)*100))%>%
  arrange(desc(perc_area)) -> area_per_class

area_per_class%>%
  #filter(perc_area >2.5)%>%
  pull(lib) -> cover_names 

soil_occu_crop_forma = soil_occu_crop%>%
  filter(lib %in% cover_names)

##### Compute buffers sizes estimations
#explication de la richesse par des combinaison linéaire de variables environnementales
# et par les variables paysagéres.
(fmla <- as.formula(paste("Richness_grid ~ Dim.1 + Dim.2 + Dim.3 +", paste(cover_names, collapse= "+"))))
(fmla <- as.formula(paste("Richness_grid ~depth_oxy+  ", paste(cover_names, collapse= "+"))))
mod1pois = siland(fmla, land = soil_occu_crop, init = c(100, 100, 100, 100, 100), data = metadata_grid, wd = 10, family = "poisson")
summary(mod1pois)
str(mod1pois)


(358.03-123.97)/358.03

covar = names(metadata_grid)[c(7:25,28:30)]
(fmla_full <- as.formula(paste("Richness_grid ~  ", paste(covar, collapse= "+"))))
mod2 = glm(fmla_full ,family=poisson(link="log"), data =metadata_grid  )
outmod = step(mod2, list(lower = ~1, upper = mod2), direction  = "both", steps = 5000)
summary(mod2)
mod_depth = glm(Richness_grid ~ depth_oxy ,family=poisson(link="log"), data =metadata_grid  )
summary(mod_depth)
hist(mod2$residuals, breaks = 20)
(358.03-327.54)/358.03
summary(glm(Richness_grid~1 ,family=poisson(link="log"), data =metadata_grid  ))

#### Shanon
plot(metadata_grid$Dim.1, metadata_grid$Simpson_grid)
summary(glm(Shannon_grid~1 ,family=gaussian(link="identity"), data =metadata_grid  ))
mod2_shanon = glm(Shannon_grid~  Dim.1 + Dim.2 + Dim.3,family=gaussian(link="identity"), data =metadata_grid  )
summary(mod2_shanon)
(16.065    -14.968  )/16.065    
(fmla <- as.formula(paste("Richness_grid ~ Dim.1 + Dim.2 + Dim.3 +", paste(cover_names, collapse= "+"))))
mod3_shanon = siland(fmla, land = soil_occu_crop, init = c(100, 100, 100, 100, 100), data = metadata_grid, wd = 10, family = "poisson")

#### Deviance >> df => surdispestion => binomial negative
library(MASS)
glm.negbin = glm.nb(Richness_grid~ Dim.1 + Dim.2 + Dim.3 + Dim.4 +Dim.5 , data =metadata_grid )
summary(glm.negbin)
(47.581  -42.963)/47.581
summary(glm.nb(Richness_grid~ 1 , data =metadata_grid ))
glm.negbin = glm.nb(Richness_grid~ depth_oxy , data =metadata_grid )
summary(glm.negbin)

glm.negbin = glm.nb(Richness_grid~ clay + MO + Cond + pHwater , data =metadata_grid)


round(exp(mod1pois$result$coefficients),1)
ggplot(metadata_grid)+
  geom_text(aes(Dim.2, Richness_grid , label = 1:42))+
  geom_smooth(aes(Dim.2, Richness_grid ), method = "lm")
  geom_abline(intercept = -4.05, slope = -0.06)
summary(mod)

likres = siland.lik(mod1pois,land = soil_occu_crop,data = metadata_grid, varnames = cover_names)
likres+
  main_theme
siland.quantile(mod1pois)

plotsiland.sif(mod1pois)+
  main_theme

#####
library(raster)
library(scalescape)
library(doParallel)
library(nlme) 


vignette("scalescape")
ext <- floor(extent(soil_occu_crop))

raster_list = list()
rr <- raster(ext, res=50)
names_lands = c("wetland", "cultivated", "natural_landscape", "non_emitting", "artificial")

temp = soil_occu_crop%>%dplyr::select(contains(names_lands[1]))
raster::rasterize(temp,rr)
cl <- makeCluster(5)
registerDoParallel(cl)    
list_raster = foreach(i=1:5, .packages = c("tidyverse", "raster","sf")) %dopar% {
  temp = soil_occu_crop%>%dplyr::select(contains(names_lands[i]))
  
  raster::rasterize(temp,rr)
  
}
stopCluster()

metadata_grid%>%
  dplyr::select(X,Y)%>%
  st_as_sf(coords = c("X","Y")) -> site
str(rr)
wetland_matrix <- landscape_matrix(list_raster[[1]], site, max.radius = 2000) 
cultivated_matrix <- landscape_matrix(list_raster[[2]], site, max.radius = 2000) 
natural_matrix <- landscape_matrix(list_raster[[3]], site, max.radius = 2000) 
non_emitting_matrix <- landscape_matrix(list_raster[[4]], site, max.radius = 2000) 
artificial_matrix <- landscape_matrix(list_raster[[5]], site, max.radius = 2000)

mod0.gls <- gls(Richness_grid ~ Dim.1 + Dim.2 + Dim.3, metadata_grid, method = "ML", correlation=corGaus(form = ~ X + Y, nugget=T))
summary(mod0.gls)
# Fit null model without landscape variables 
mod0 <- glm(Richness_grid~Dim.1 + Dim.2 + Dim.3, family=poisson(link="log"), data = metadata_grid)
#run dist weight 
mod <- dist_weight(mod0 = mod0, landscape.vars = list(wetland = wetland_matrix, 
                                                      cultivated = cultivated_matrix,
                                                      natural = natural_matrix, 
                                                      non_emitting = non_emitting_matrix,
                                                      artificial = artificial_matrix),
    landscape.formula = '~ . + wetland + cultivated + natural + non_emitting + artificial', data = metadata_grid) 

##### Bootstrap of landscape
# library(doParallel)
# library(foreach)
# land_bootstrap = function(soil_occu_forma, metadata_grid){
#     # soil_occu_forma %>%
#     #   select(lib)%>%
#     #   mutate(lib = sample(lib))%>%
#     #   mutate(area =  st_area(.), #area of each polygone
#     #          n =1,
#     #          lib2 = as.factor(lib))%>%
#     #   pivot_wider(names_from = lib2, values_from = n, values_fill = 0) -> soil_occu_forma
#     # 
#     soil_occu_forma%>%
#       pull(lib)%>%
#       unique() -> cover_names
#   ID = sample(1:42)
#   temp= metadata_grid
#   temp[,names(temp)%in%c("Richness_grid","Dim.1", "Dim.2", "Dim.3")] = temp[ID,names(temp)%in%c("Richness_grid","Dim.1", "Dim.2", "Dim.3")]
#     # metadata_grid %>%
#     #   mutate(Richness_grid = sample((Richness_grid))) -> temp
#   (fmla <- as.formula(paste("Richness_grid ~ Dim.1 + Dim.2 + Dim.3 +", paste(cover_names, collapse= "+"))))
#   
#   mod = siland(fmla, land = soil_occu_forma, init = c(100, 100, 100, 100, 100), data = temp, wd = 10, family = "poisson")
#   return(mod)
# }#mod_test = siland(fmla, land = soil_occu_forma, init = c(100, 100, 100, 100, 100), data = metadata_grid, wd = 20)
# 
# #vignette("foreach")
# 
# #Create cluster with desired number of cores (careful need lot of RAM)  
# #core processes
# cl <- makeCluster(4)
# registerDoParallel(cl)
# 
# test <- foreach(n_iter = 1:80, .packages = c("tidyverse", "siland", "sf")) %dopar% {
#   tryCatch({
#     land_bootstrap(soil_occu_forma, metadata_grid)
#   }, error = function(e) {
#     cat("Error in iteration", n_iter, ":", conditionMessage(e), "\n")
#     NULL  # Returning NULL to avoid issues with gathering results
#   })
# }
# stopCluster(cl)
# 
# save(test, file = "outputs/richenss_bootstrap.RData")
# 
# vignette("siland")
# summary(test[[4]])
# test[[6]]$pval0
# siland.lik(test[[4]],land = soil_occu_forma,data = metadata_grid, varnames = cover_names)
# test[[6]]$loglik
# test[[6]]$loglik0
# test[[6]]$landcontri
# str(test[[2]]$result$R)
# test[[4]]$result$converged
# 
# AIC_mod = c()
# AIC_null = c()
# nolandscape_effect = c()
# for(i in 1:length(test)){
#   AIC_mod = c(AIC_mod, test[[i]]$AIC)
#   AIC_null = c(AIC_null,test[[i]]$AIC0)
#   nolandscape_effect = c(nolandscape_effect, test[[i]]$pval0)
#   
# }
# par(mfrow = c(1,1))
# hist(AIC_null-AIC_mod, breaks = 60)
# abline(v =( mod1pois$AIC0-mod1pois$AIC), col= "red")
# hist(nolandscape_effect, breaks = 60)
# abline(v =mod1$pval0, col= "red")


buffer.around.point = function(pt, geo_data, layer_name, buffer_size){
  pt = as.matrix(pt)
  if(nrow(pt) == 2){
    pt = t(pt)
  }
  pt_vect = terra::vect(pt, type="points", atts=NULL, crs=terra::crs(geo_data))
  buffer_vec <- terra::buffer(pt_vect, buffer_size)
  crop_suface = terra::crop(geo_data, buffer_vec)
  
  names(crop_suface)[names(crop_suface) == layer_name] = "focal_layer"
  
  sufaces_class = data.frame( class= crop_suface$focal_layer, surface = expanse(crop_suface, unit="m", transform=TRUE))
  sufaces_per_class = sufaces_class%>%
    group_by(class)%>%
    summarise(sum_area = sum(surface, na.rm = T))%>%
    mutate(perc_area = sum_area/sum(sum_area)*100,
           X = pt[1],
           Y = pt[2])%>%
    arrange(desc(perc_area))
  
  return(list(sufaces_per_class = sufaces_per_class, crop_suface = st_as_sf(crop_suface) ))
}

##### Compute area percentage within a buffer
buffer_size = 500 # buffer radius
croped_surfaces = apply(metadata_grid[,2:3], 1, function(x)
  buffer.around.point(pt = x, geo_data = vect(soil_occup_forma), layer_name = "lib", buffer_size = buffer_size))

# croped_surfaces contain percentage area dataframes and also a sp objects of 
# each buffer (if needed for graphical representation of the landscape for exemple)

##### Extract and group percentage area dataframes into one dataframe

area_per_buffer = croped_surfaces[[1]][[1]]
area_per_buffer$Grid_code =metadata_grid$Grid_code[1]
for (i in 2:length(croped_surfaces)){
  temp = croped_surfaces[[i]][[1]]
  temp$Grid_code = metadata_grid$Grid_code[i]
  area_per_buffer = rbind(area_per_buffer ,temp)
}

ggplot(area_per_buffer)+
  geom_point(aes(x=class, perc_area),cex = 2)+
  facet_wrap(~Grid_code)

area_per_buffer$buffer_size = buffer_size
area_per_buffer%>%
  select(class, perc_area, Grid_code,buffer_size)%>%
  pivot_wider(names_from = class, values_from = perc_area, values_fill = 0) -> area_per_buffer_wide 


area_per_buffer_wide%>%
  write.table(file = "data/clean_data/Buffer_cover_CAM.txt")
