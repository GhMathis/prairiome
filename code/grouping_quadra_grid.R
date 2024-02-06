metadata_sample = read_xlsx("data/Metadata_Sample.xlsx")
OTU_plant = read.table("data/otu_plant_grid_CAM.txt")
str(OTU_plant)
str(metadata_sample)
head(metadata_sample)
metadata_sample%>%
  filter(str_detect(Locality, "Arles")) -> metadata_sample_cam

otu_plant = read_xlsx("data/OTU_plant.xlsx")%>%
  filter(str_detect( Host_code, "CAM"),
         !str_detect(Host_code, "21_CAM_1[45]")) %>%
  replace(is.na(.),0)

library(GGally)
##########Host_code


otu_plant_long = otu_plant%>%
  
  mutate(year = as.factor(substring(Host_code,1,2)),
         grid = as.factor(substring(Host_code,8,9)),
         quadra = as.factor(substring(Host_code,10,11)),
         total_cover = rowSums(across(where(is.numeric))),
         plant_richness = rowSums(across(where(is.numeric)) != 0)
         
  )%>%
  select(c(year, grid, quadra, total_cover, plant_richness, Host_code))-> data_quadra
 


ggplot(metadata_sample_cam)+
  geom_point(aes(asBiomass, Soil))

data_quadra$Biommass = as.numeric(metadata_sample_cam$Biomass)
data_quadra$Soil = as.numeric(metadata_sample_cam$Soil)
data_quadra$Litter = as.numeric(metadata_sample_cam$Litter)
data_quadra$Rock = as.numeric(metadata_sample_cam$Rocks)

ggplot(data_quadra)+
  geom_boxplot(aes(as.factor(plant_richness), Biommass))

ggplot(data_quadra)+
  geom_boxplot(aes(as.factor(Soil), plant_richness))

ggpairs(data_quadra%>%
         select(Biommass, Soil, plant_richness))
