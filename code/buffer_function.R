library(readxl)
library(stringr)
library(alm)
library(terra)
library(sp)
library(sf)

Metadata_grid = read_excel(
  "data/Metadata_Grid.xlsx"
)%>%
  filter(Locality == "Arles")%>%
  mutate(lat = as.numeric(str_split_fixed(GPS_localisation, pattern = ",",n=2)[,1]),
           lon = as.numeric(str_split_fixed(GPS_localisation, pattern = ",",n=2)[,2]))

data_habitat <- read.csv("data_habitat.csv", header = T, sep=";")%>%
  select(surf_ha,perimetre,contains("lib") ,id)


str(Metadata_grid)
str(data_habitat)

##### Create buffer around grid

# buffer.aroud.point = function(pt, geo_data, buffer_size){
#   
#   p = vect(t(as.matrix(pt)), type="points", atts=NULL, crs=crs(geo_data))
#   b <- terra::buffer(p, buffer_size)
# 
#   return(crop(geo_data, b))
# }

map_hab <- st_read("sql_statement_d551600.shp")
plot(map_hab)

# df lon/lat grid pos
grid_pos = as.data.frame(cbind(X = Metadata_grid$lon, Y = Metadata_grid$lat))
write.table(grid_pos, "data/grid_location.txt")
grid_pos$Buffer = paste(1:nrow(grid_pos),"_id",sep="")
# bff = buffer.aroud.point(pt = grid_pos[7,], geo_data = map_hab, buffer_size = 1000) #test

# Buffer creation, center on grid pos
list_buffer <- alm_create_buffer(grid_pos[1,], # the previous table
                                 set.endCapStyle = "ROUND", # type of Buffer: either ROUND or SQUARE
                                 set.dist = 200, # a distance in meters : radius (ROUND), side (SQUARE)
                                 set.crs = st_crs(map_hab)
)

par(mfrow = c(1,1))
str(map_hab)
plot(map_hab["lib1_16"], key.width = lcm(1))
for (i in 1:length(list_buffer)){
  points(list_buffer[[i]])
}



list_shape <- alm_set_shp_sources(c("shp1"), c("data/shp1.shp"),
                                  c("lib1_16"), 1,
                                  set.crs = st_crs(map_hab),
                                 #set.save = TRUE, # want to save the output ?
                                 trunc.with.buffers = TRUE, # only keep the useful part ?
                                 list.buffer = list_buffer, # created in the previous step
                                 check.validity = TRUE, # want to check shape validity ?
                                 check.fillable = TRUE) # want to know how much we can fill of all 
list_shape$shp1$lib1_16
str(list_shape)
plot(list_shape$shp1, key.width = lcm(2))
plot(map_hab["lib1_16"], key.width = lcm(2))
points(list_shape$shp1)
# buffers with our source shapefiles ?

# v = vect(grid_pos, type="points", atts=NULL, crs=crs(p))
# 
# b100 <- buffer(v, 100)
# b1000 <- buffer(v, 1000)
# 
# par(mfrow = c(1,2))
# plot(p,xlim = c(4.5,4.8), ylim = c(43.39,43.65))
# points(b100)
# points(v, col ="red")
# 
# plot(p,xlim = c(4.5,4.8), ylim = c(43.39,43.65))
# points(b1000)
# points(v, col ="red")
# 
# v1 = vect(grid_pos[,1], type="points", atts=NULL, crs=crs(p))
# mask()
# b1 <- buffer(v, 1000)
# cm = mask(p, b1000)
# str(cm)
# plot(cm)
#####

dim(list_shape[[1]])
str(list_buffer)
str(list_shape)

empty = list_buffer[[1]]
st_geometry(empty) <- "geometry"
str(d_geom)
st_geometry_type(empty$geometry)
d_geom <- empty %>% dplyr::filter(grepl("GEOMETRYCOLLECTION", sf::st_geometry_type(geometry)))

df_sans_d <- empty %>% dplyr::filter(grepl("POLYGON", sf::st_geometry_type(geometry)))

d_cor <- suppressWarnings(sf::st_cast(d_geom)[which(sf::st_is(sf::st_cast(d_geom), c("POLYGON", "MULTIPOLYGON"))),])

data_final_cor <- dplyr::bind_rows(df_sans_d, d_cor)
#empty = data_final_cor
x= list_shape[[1]]
sf::st_agr(x) <- "constant"
sf::st_agr(empty) <- "constant"
list_intersection = tryCatch(sf::st_intersection(x,empty),
         error = function(e) {
           message("\n [B1] Need cut for current shape : ", "sph1")
           empty <- sf::st_buffer(empty, dist = -0.0001)
           return(sf::st_intersection(x, empty))
         })

sf::st_agr(list_intersection) <- "constant"

# correction linéaire
#list_intersection <- alm_get_data_final(list_intersection)

d_geom <- list_intersection %>% dplyr::filter(grepl("GEOMETRYCOLLECTION", sf::st_geometry_type(geometry)))

df_sans_d <- list_intersection %>% dplyr::filter(grepl("POLYGON", sf::st_geometry_type(geometry)))
plot((d_geom))
st_make_valid(d_geom)$Buffer
d_cor <- suppressWarnings(sf::st_cast(d_geom)[which(sf::st_is(sf::st_cast(d_geom), c("POLYGON", "MULTIPOLYGON"))),])

data_final_cor <- dplyr::bind_rows(df_sans_d, d_cor)
# Second intersection to avoid bugs
st_geometry(list_intersection)
st_geometry(empty)
st_geometry(x)
st_is_valid(st_make_valid(list_intersection))

#st_make_valid
list_intersection <- tryCatch(sf::st_intersection(st_make_valid(list_intersection), empty),
                              error = function(e) {
                                message("\n [B1] Need cut for current shape : ", "shp1")
                                empty <- sf::st_buffer(empty, dist = -0.0001)
                                return(sf::st_intersection(list_intersection, empty))
                              })


data_final <- alm_run(to.fill = "all", # here we want to fill all previously specified buffers
                     list.buffer = list_buffer, # created in the first step
                     list.shape = list_shape, # create in the second step
                     log.save = FALSE, # do we want to save console logs ?
                     log.display = TRUE, # do we want to display them too ?
                     save = TRUE, # do we want to save last results ?
                     save.state = TRUE, # do we want to save state after each buffer has been processed ?
                     use.save.state = FALSE, # do we want to use some previous save states ?
                     set.dist = -0.0001) # a distance (in meter) parameter for some part of the function
e <- simpleError("test error")
test = tryCatch(stop(e), error = function(e) e, finally = print("Hello"))
test
str(test)
