library(readxl)
library(sp)
library(sf)
library(tidyverse)
library(siland)
library(terra)

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


##### Load shapefile of soil occupation
st_read("data/shapefiles/sql_statement_d551600.shp")  %>% 
  select(lib = lib1_16)%>%
  mutate(lib = str_replace_all(lib, "[,/ ']","_"),
         lib = str_replace_all(lib, "[.]","_"),
         lib = str_replace_all(lib, "[êéè]","e"),
         lib = str_remove(lib, "[*]$"),
         lib = str_remove(lib, "_$"))%>%
  st_transform( 2154)%>%
  
  mutate(area =  st_area(.), #area of each polygone
         n =1,
         lib2 = as.factor(lib))%>%
  pivot_wider(names_from = lib2, values_from = n, values_fill = 0)-> soil_occu
(soil_occu)
unique(soil_occu$lib)
##### Load positions of grids and convert to lambert projection (for siland package)
read_excel("data/Metadata_sample.xlsx")%>%
  
  filter(Locality == "Arles") %>%
  select(Biomass, Grid_code)%>%
  group_by(Grid_code)%>%
  summarise(biomass_mean = mean(Biomass),
            biomass_sd = sd(Biomass)) -> biomass

read_excel("data/Metadata_Grid.xlsx")%>%
  filter(Locality == "Arles") %>%
  mutate(
    X = as.numeric(str_split_fixed(GPS_localisation, pattern = ",", n = 2)[, 2]),
    Y = as.numeric(str_split_fixed(GPS_localisation, pattern = ",", n = 2)[, 1])
  ) %>%
  select(X, Y, Grid_code) %>%
  st_as_sf(coords = c("X", "Y"), crs = 4326)%>%
  st_transform(grid_pos_metadata, crs = st_crs(soil_occu))%>%
  mutate(X = as.numeric(st_coordinates(geometry)[, 1]),
         Y =  as.numeric(st_coordinates(geometry)[, 2]))%>%
  st_join(.,
    soil_occu%>%
            select(lib,area)
          )%>%
  distinct(Grid_code, .keep_all = TRUE)%>% # duplication is apprening in the st_joit above for some reason... Idk
  as.data.frame(.)%>%
  select(!geometry)%>%
  left_join(read.table("data/grid_abudance.txt",header = T),by = "Grid_code")%>%
  left_join(biomass,by = "Grid_code")-> grid_pos_metadata


##### Some visualisations
ggplot(grid_pos_metadata)+
  geom_point(aes(plant_abundance_grid, biomass_mean), cex = 3)+
  main_theme
ggplot(grid_pos_metadata)+
  geom_point(aes(plant_abundance_grid, biomass_sd), cex = 3)+
  main_theme
map = ggplot()+
  #geom_sf(data = soil_occu,col = "gray",aes(fill = lib))+
  #geom_sf(data = limit_map, col ="red", fill = NA)+
  geom_point(data = grid_pos_metadata, aes(X,Y), cex = 2.2, fill ="black")+
  geom_point(data = grid_pos_metadata, aes(X,Y, col =as.factor(year)), cex = 2)+
  main_theme

map

##### Format df of grid position to respect for siland
grid_pos_metadata%>%
  rename( obs = "plant_abundance_grid")%>%
  mutate(Id = 1:n())%>%
  select(X,Y, lib,area, year, biomass_mean, biomass_sd,Id,obs)%>%
  mutate(across(c(Id,obs), as.numeric),
         year = as.character(year)) -> grid_pos_metadata_forma


##### Compute area per cover type over all the zone observed
as.data.frame(soil_occu)%>%
  select(lib, area)%>%
  group_by(lib)%>%
  summarise(sum_area = sum(area, na.rm = T))%>%
  mutate(perc_area = as.numeric(sum_area/sum(sum_area)*100))%>%
  arrange(desc(perc_area)) -> area_per_class

str(area_per_class)
area_per_class%>%
  filter(perc_area >2.5)%>%
  pull(lib) -> cover_names 

soil_occu_forma = soil_occu%>%
  filter(lib %in% cover_names)

(fmla <- as.formula(paste("obs ~ biomass_mean +", paste(cover_names, collapse= "+"))))
mod = siland(fmla, land = soil_occu_forma, init = c(100, 100, 100, 100), data = grid_pos_metadata_forma, wd = 5)

summary(mod)
siland.quantile(mod)

plotsiland.sif(mod)
likres = siland.lik(mod,land= soil_occu_forma,data=grid_pos_metadata_forma,varnames=cover_names)
likres


###### Compute percetage of cover type within a buffer arround grid 


# Buffer size are estimated with the model in 01_compute_buffer_sizes. 


##### Function to compute area percentage 
buffer.around.point = function(pt, geo_data, layer_name, buffer_size){
  pt = as.matrix(pt)
  if(nrow(pt) == 2){
    pt = t(pt)
  }
  pt_vect = terra::vect(pt, type="points", atts=NULL, crs=crs(geo_data))
  buffer_vec <- terra::buffer(pt_vect, buffer_size)
  crop_suface = terra::crop(geo_data, buffer_vec)
  
  names(crop_suface)[names(crop_suface) == layer_name] = "focal_layer"

  sufaces_class = data.frame( class= crop_suface$focal_layer, surface = expanse(crop_suface, unit="m", transform=TRUE))
  print(sufaces_class)
  print(pt)
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
buffer_size = 300 # buffer radius
croped_surfaces = apply(grid_pos_metadata[,2:3], 1, function(x)
  buffer.around.point(pt = x, geo_data = vect(soil_occu_forma), layer_name = "lib", buffer_size = buffer_size))

# croped_surfaces contain percentage area dataframes and also a sp objects of 
# each buffer (if needed for graphical representation of the landscape for exemple)

##### Extract and group percentage area dataframes into one dataframe

area_per_buffer = croped_surfaces[[1]][[1]]
area_per_buffer$Grid_code =grid_pos_metadata$Grid_code[1]
for (i in 2:length(croped_surfaces)){
  temp = croped_surfaces[[i]][[1]]
  temp$Grid_code = grid_pos_metadata$Grid_code[i]
  area_per_buffer = rbind(area_per_buffer ,temp)
}

ggplot(area_per_buffer)+
  geom_point(aes(x=class, perc_area),cex = 2)+
  facet_wrap(~Grid_code)

area_per_buffer$buffer_size = buffer_size
area_per_buffer%>%
  select(class, perc_area, Grid_code,buffer_size)%>%
  pivot_wider(names_from = class, values_from = perc_area, values_fill = 0) -> area_per_buffer_wide 

grid_pos_metadata%>%
  left_join(area_per_buffer_wide,by = "Grid_code")%>%
  write.table(file = "data/Metadata_Grid_CAM.txt")

