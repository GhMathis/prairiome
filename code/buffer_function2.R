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
  select(geometry,lib = lib1_16)%>%
  mutate(lib = str_replace_all(lib, "[,/ ']","_"),
         lib = str_replace_all(lib, "[.]","_"),
         lib = str_remove(lib, "[*]$"),
         lib = str_remove(lib, "_$"))%>%
  mutate(lib = case_when(
    lib %in% c("Espace_bâti_diffus_et_autre_bâti","Tissu_urbain_discontinu",
               "Tissu_urbain_continu","Décharge","Espace_ouvert_urbain",
               "Équipement_sportif_et_de_loisirs","Zone_d_activité_et_équipement",
               "Extraction_de_matériaux", "Chantier","Zone_portuaire") ~ "zone_urbaine",
    lib %in% c("Plage", "Dune") ~ "littoral",
    lib %in% c("Verger__oliveraie", "Culture_maraichère", "Vignoble") ~ "zone_de_culture",
    #lib %in% c("Zones_humides", "Forêts_et_milieux_naturels_terrestres", "Zones_en_eau") ~ "milieux_naturel", #for lib 1 if nedded
    .default = lib
  ))%>%
  st_transform( 2154)-> soil_occu
unique(soil_occu$lib)

#area of each polygone
soil_occu$area = soil_occu%>%
  st_area()

soil_occu %>%
  mutate(n =1,
         lib2 = lib) %>% # duplicate lib (need it below to compute area per cover class)
  pivot_wider(names_from = lib2, values_from = n, values_fill = 0) -> soil_occu # make a complete disjunctive table



##### Load positions of grids and convert to lambert projection (for siland package)
read_excel("data/Metadata_Grid.xlsx") %>%
  filter(Locality == "Arles") %>%
  mutate(
    X = as.numeric(str_split_fixed(GPS_localisation, pattern = ",", n = 2)[, 2]),
    Y = as.numeric(str_split_fixed(GPS_localisation, pattern = ",", n = 2)[, 1])
  ) %>%
  select(X, Y, Grid_code, Ecosystem) %>%
  st_as_sf(coords = c("X", "Y"), crs = 4326)%>%
  st_transform(grid_pos, crs = 2154)%>%
  mutate(X = st_coordinates(geometry)[, 1],
         Y =  st_coordinates(geometry)[, 2])%>%
  st_join(soil_occu%>%
            select(lib,area))%>%
 
  left_join(read.table("data/grid_abudance.txt",header = T),by = "Grid_code") -> grid_pos



##### Limit the soil occupation shapefile given a buffer around point to limit computation time with siland
limit_map = st_buffer(grid_pos$geometry, 4000)%>%
  st_union()%>%
  st_cast(to = "POLYGON")
soil_occu_crop = st_crop(x = soil_occu, y =limit_map)

map = ggplot()+
  geom_sf(data = soil_occu_crop,col = "gray",aes(fill = lib))+
  #geom_sf(data = limit_map, col ="red", fill = NA)+
  geom_sf(data = grid_pos, cex = 2.2, fill ="black")+
  geom_sf(data = grid_pos, cex = 2)+
  main_theme
over
map

##### Format df of grid position to respect for siland

grid_pos_forma = as.data.frame(grid_pos)%>%
  rename( obs = "plant_abundance_grid")%>%
  mutate(Id = 1:n())%>%
  select(X,Y,Ecosystem, lib, area, year,Id,obs)%>%
  mutate(across(c(Id,obs), as.numeric),
         year = as.character(year))
str(grid_pos_forma)

##### Compute area per cover type over all the zone observed
area_per_class = as.data.frame(soil_occu_crop)%>%
  select(lib, area)%>%
  group_by(lib)%>%
  summarise(sum_area = sum(area, na.rm = T))%>%
  mutate(perc_area = as.numeric(sum_area/sum(sum_area)*100))%>%
  arrange(desc(perc_area))

str(area_per_class)
cover_names =  area_per_class%>%
  filter(perc_area >2.5)%>%
  pull(lib)


hist(grid_pos_forma$obs, breaks =40)
unique(grid_pos_forma$Ecosystem)
str(grid_pos_forma)
t1 = Sys.time()
TEST = glm(obs ~ lib, dat = grid_pos_forma)
par(mfrow = c(2,2))
hist(TEST$residuals, breaks = 20)
plot(TEST)
ggplot(grid_pos_forma)+
  geom_boxplot(aes(lib,obs))+
  main_theme

grid_pos_forma$year == 20
summary(TEST)
AIC(TEST)
grid_pos_forma$lib
(fmla <- as.formula(paste("obs ~ lib +", paste(cover_names, collapse= "+"))))
mod = siland(fmla, land = soil_occu_crop%>%select(cover_names), data = grid_pos_forma, init = c(100,100,100,100), wd = 10, maxD = 2000)
t2 = Sys.time()
mod
t2-t1
mod$coefficients
summary(mod)
plotsiland.sif(mod)
siland.quantile(mod)
# d = 1:500
# func.exp = function(d,delta){(2/(pi*delta**2))*exp(-2*d/delta)}
# plot(d, func.exp(d,250), "l", col = "red")
# for(delta in seq(50,200, by= 50)){
#   points(d, func.exp(d,delta), "l")
#   print(mean(func.exp(d,delta)))
# }
d = seq(0,2000,length.out = 1000)
func.gaus = function(d,delta,beta){beta*(1/(2*delta*sqrt(pi))*exp(-(pi*d)/(2*delta))**2)}
delta = mod$paramSIF
beta = mod$coefficients[4:7]
col = c("red", "blue", "green", "orange", "purple")
par(mfrow = c(1,1))
plot(d, func.gaus(d,delta[1], beta[1]), "l", col = col[1])
#abline(h=mean(func.gaus(d,delta[1], beta[1])), lty = 3)
#abline(h=median(func.gaus(d,delta[1], beta[1])), lty = 2)
abline(v=100,lty = 3)
for(i in 2:length(delta)){
  points(d, func.gaus(d,delta[i], beta[i]), "l", col = col[i])
  print(mean(func.gaus(d,delta[i], beta[i])))
}
par(mfrow = c(2,2))
plot(mod)
plot(grid_pos_forma$x1,grid_pos_forma$obs)
mod_test = lm(obs~x1, data = grid_pos_forma)
plot(mod_test)
plotsiland.sif(mod)
likres=siland.lik(mod,land= soil_occu_crop%>%select(cover_names),data=grid_pos_forma,varnames=cover_names)
likres

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
