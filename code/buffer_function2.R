library(readxl)
library(sp)
library(sf)
library(tidyverse)
library(siland)
library(terra)
library(ade4)
library(FactoMineR)

vignette("siland")

##### graphics setup #####
main_theme = theme_bw()+
  theme(line = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks =  element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size=10),
        axis.text.y = element_text(colour = "black", size=10),
        legend.title = element_text(colour = "black", size=20),
        legend.title.align=0.5,
        legend.text = element_text(colour = "black", size=15),
        axis.title=element_text(size=10),
        strip.background = element_rect(fill="cornsilk"))
##### Load positions of grids and convert to lambert projection (for siland package)
read_excel("data/Metadata_Grid.xlsx") %>%
  filter(Locality == "Arles") %>%
  mutate(
    X = as.numeric(str_split_fixed(GPS_localisation, pattern = ",", n = 2)[, 2]),
    Y = as.numeric(str_split_fixed(GPS_localisation, pattern = ",", n = 2)[, 1])
  ) %>%
  select(X, Y, Grid_code) %>%
  st_as_sf(coords = c("X", "Y"), crs = 4326)%>%
  st_transform(grid_pos, crs = 2154)%>%
  mutate(X = st_coordinates(geometry)[, 1],
         Y =  st_coordinates(geometry)[, 2])%>%
  left_join(read.table("data/grid_abudance.txt",header = T),by = "Grid_code") -> grid_pos


# 
# soil_occu <- st_read("data/shapefiles/bd_ocsol_evol_2006_2014_v2/soil_occupation.shp")%>%
#   head(100)
# data_lambert <- st_transform(grid_pos, crs = 2154) 
# data_lambert$X <- st_coordinates(data_lambert)[, 1]
# soil_occu_binary = soil_occu %>% 
#   select(geometry)
# 
# soil_occu = soil_occu %>%
#   #select(geometry, LIB3_14)%>%
#   mutate(LIB3_14 = str_replace_all(LIB3_14, "[,/ ')()]","_"),
#          LIB3_14 = str_replace_all(LIB3_14, "[.]","_"),
#          LIB3_14 = str_remove(LIB3_14, "[*]$"),
#          LIB3_14 = str_remove(LIB3_14, "_$"))%>%
#   group_by(row = row_number(), LIB3_14) %>%
#   tally() %>%
#   pivot_wider(names_from = LIB3_14, values_from = n, values_fill = 0)%>%
#   select(-row)


##### Load shapefile of soil occupation
st_read("data/shapefiles/sql_statement_d551600.shp") %>% 
  select(geometry,lib = lib3_16)%>%
  mutate(lib = str_replace_all(lib, "[,/ ']","_"),
         lib = str_replace_all(lib, "[.]","_"),
         lib = str_remove(lib, "[*]$"),
         lib = str_remove(lib, "_$")) -> soil_occu

soil_occu <- st_transform(soil_occu, 2154) # convert to Lambert

soil_occu %>%
  mutate(n =1) %>%
  pivot_wider(names_from = lib, values_from = n, values_fill = 0) -> soil_occu # make a complete disjunctive table


##### Limit tthe soil occupation shape file given a buffer around point to limit computation time with siland
limit_map = st_buffer(grid_pos$geometry, 5000)%>%
  st_union()%>%
  st_cast(to = "POLYGON")
soil_occu_crop = st_crop(x = soil_occu, y =limit_map)

map = ggplot()+
  geom_sf(data = soil_occu_crop,col = "gray", fill = NA)+
  geom_sf(data = limit_map, col ="red", fill = NA)+
  geom_sf(data = grid_pos)+
  main_theme

map
##### Format df of grid position to respect for siland
soil_occu_crop = soil_occu_crop %>%
  mutate(across(where(is.integer), as.numeric))

grid_pos_forma = as.data.frame(grid_pos)%>%
  rename( x1= "year", obs = "plant_abundance_grid")%>%
  mutate(Id = 1:n())%>%
  select(X,Y,x1,Id,obs)%>%
  mutate(across(c(Id,obs,x1), as.numeric))


cover_names = names(soil_occu_crop)[names(soil_occu_crop) != "geometry"]




mod0 = siland("obs ~ x1", land = soil_occu_crop, data = grid_pos_forma,wd = 10, maxD=1000)
for(i in 1: ){
  (fmla <- as.formula(paste("obs ~ x1 + ", paste(cover_names, collapse= "+"))))
  mod=siland(fmla, land = soil_occu_crop, data = grid_pos_forma,wd = 10, maxD=1000)
  AIC_prev = mod$AIC
}

##### Create buffer around grid
buffer.aroud.point = function(pt, geo_data, layer_name, buffer_size){
  pt = as.matrix(pt)
  if(nrow(pt) == 2){
    pt = t(pt)
  }
  pt_vect = terra::vect(pt, type="points", atts=NULL, crs=crs(geo_data))
  buffer_vec <- terra::buffer(pt_vect, buffer_size)
  crop_suface = terra::crop(geo_data, buffer_vec)
  
  names(crop_suface)[names(crop_suface) == layer_name] = "focal_layer"
  sufaces_class = crop_suface$focal_layer
  sufaces_class = cbind(sufaces_class, surface = expanse(crop_suface, unit="m", transform=TRUE))

  sufaces_per_class = sufaces_class%>%
    group_by(focal_layer)%>%
    summarise(sum_area = sum(surface, na.rm = T))%>%
    mutate(perc_area = sum_area/sum(sum_area)*100,
           X = pt[1],
           Y = pt[2])%>%
    arrange(desc(perc_area))
  
  return(list(sufaces_per_class = sufaces_per_class, crop_suface = crop_suface))
}


croped_surfaces = apply(grid_pos[1:2,1:2], 1, function(x) buffer.aroud.point(pt = x, geo_data = soil_occu2, layer_name = "lib3_16", buffer_size = 300))
buffer.aroud.point(pt = grid_pos[1,1:2], geo_data = soil_occu2, layer_name = "lib3_16", buffer_size = 300)
for(i in 1:lenght(croped_surfaces)){
  croped_surfaces[[i]]$surface_per_class
}
