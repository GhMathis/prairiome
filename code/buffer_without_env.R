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
library(car)
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
(fmla <- as.formula(paste("Plant_richness_grid ~ Dim.1 + Dim.2 + Dim.3 +", paste(cover_names, collapse= "+"))))
mod1pois = siland(fmla, land = soil_occu_crop, init = c(100, 100, 100, 100, 100), data = metadata_grid, wd = 10, family = "poisson")

mod1pois
summary(mod1pois)


likres = siland.lik(mod1pois,land = soil_occu_crop,data = metadata_grid, varnames = cover_names)
likres+
  main_theme
siland.quantile(mod1pois)

plotsiland.sif(mod1pois)+
  main_theme
######

metadata_grid
metadata_grid%>%
  mutate(Na_class= cut(Na,
                       breaks = c(0, 500,2000, Inf),
                       labels = c("0-500", "500-2000", ">2000"),
                       include.lowest = TRUE),
         clay_class= cut(clay,
                         breaks = c(0, 25, Inf),
                         labels = c("0-25", "25-50"),
                         include.lowest = TRUE)) -> metadata_grid


mod_no_land <- glm(Plant_richness_grid~Pature*Fauche + Na_class +clay_class , family=poisson(link="log"), data = metadata_grid)
summary(mod_no_land)
hist(mod_no_land$residuals, breaks = 20)
Anova(mod_no_land,3)

(mod_no_land$null.deviance - mod_no_land$deviance)/mod_no_land$null.deviance

metadata_grid$richness_resid = mod_no_land$residual

(fmla_res <- as.formula(paste("richness_resid~1 +", paste(cover_names, collapse= "+"))))
mod_res = siland(fmla_res, land = soil_occu_crop, data = metadata_grid,init = c(100,500,300,100,500), wd = 5, family = "gaussian")
mod_res
summary(mod_res)
likres2 = siland.lik(mod_res,land = soil_occu_crop,data = metadata_grid, varnames = cover_names)
likres2+
  main_theme
siland.quantile(mod_res)
#####
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
buffer_size = 1500 # buffer radius
croped_surfaces = apply(data.frame(metadata_grid$X, metadata_grid$Y), 1, function(x)
  buffer.around.point(pt = x, geo_data = vect(soil_occu_crop), layer_name = "lib", buffer_size = buffer_size))

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
  write.table(file = "data/data_clean/Buffer_5ldscp_CAM.txt")
