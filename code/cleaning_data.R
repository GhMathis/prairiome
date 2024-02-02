library(tidyverse)
library(readxl)

##### OTU Camargue
read_xlsx("data/OTU_plant.xlsx")%>%
  filter(str_detect( Host_code, "CAM"),
         !str_detect(Host_code, "21_CAM_1[45]")) %>%
  replace(is.na(.),0)%>%
  select(-Richness)%>%
  pivot_longer(-Host_code, names_to = "Plant", values_to = "cover")  %>%
  filter(cover != 0)%>%
  
  pivot_wider(values_from = cover, names_from = Plant, values_fill = 0) -> otu_plant_cam
write.table(otu_plant_cam, "data/data_clean/OTU_plant_CAM.txt")

##### Metadata Quadra 
read_xlsx("data/Metadata_Sample.xlsx")%>%
  filter(str_detect(Locality, "Arles"))%>%
  select(-c(Rocks, Additional_information))-> metadata_quad_cam

read_xlsx("data/OTU_plant.xlsx")%>%
  
  filter(str_detect( Host_code, "CAM"),
         !str_detect(Host_code, "21_CAM_1[45]")) %>%
  replace(is.na(.),0)%>%
  mutate(year = as.factor(substring(Host_code,1,2)),
         grid = as.factor(substring(Host_code,8,9)),
         quadra = as.factor(substring(Host_code,10,11)),
         total_cover = rowSums(across(where(is.numeric)))
         
  )%>%
  pivot_longer(-c(Host_code,Richness, year, grid, quadra, total_cover), names_to = "Plant", values_to = "cover") %>%
  filter(cover != 0)%>%
  pivot_wider(values_from = cover, names_from = Plant)%>%
  select(Richness,year,grid,quadra, total_cover,Host_code)%>%
  left_join(metadata_quad_cam, by ="Host_code")%>%
  mutate(across(c(Vegetation, Litter, Soil),~as.numeric(.x)))%>%
  mutate(Soil = case_when(str_detect(Host_code, "20_CAM_01")~ 0,
         .default = Soil))%>%
  mutate(Free_space = Litter +Soil) -> metadata_quad_cam
  

write.table(metadata_quad_cam, "data/data_clean/Metadata_quadra_CAM.txt")

##### Metadata grid


read_xlsx("data/donnee_sol.xlsx",col_names = F, skip = 3)%>%
  rename_with(~c("Grid_code", "Num", "pHwater", "lime_tot", "MO", "Phos", "K", "Mg",
                 "Ca", "Na", "N","C","CN","clay",  "SiltF", "SiltC","SandF", "SandC",
                 "Cl", "Res", "Cond"))%>%
  filter(str_detect(Grid_code,"CAM" ))%>%
  mutate(Grid_code = str_extract(Grid_code,".._CAM_.." ))%>%
  left_join(y = read_xlsx("data/distance_fer_sol_EDGG_CAM.xlsx"), by = join_by(Grid_code == EDGG))%>%
  mutate(Profondeur = as.numeric(str_remove(Profondeur, "[>]")),
         depth_oxy = cut(Profondeur,
                         breaks = c(0, 9, 19, 29, 39, Inf),
                         labels = c("0-9", "10-19", "20-29", "30-39", ">40"),
                         include.lowest = TRUE)) -> soil_data

read_xlsx("data/Metadata_grid.xlsx")%>%
  filter(Locality == "Arles") %>%
  mutate(
    X = as.numeric(str_split_fixed(GPS_localisation, pattern = ",", n = 2)[, 2]),
    Y = as.numeric(str_split_fixed(GPS_localisation, pattern = ",", n = 2)[, 1])
  )%>%
  select(-GPS_localisation)%>%
  left_join(soil_data, by= "Grid_code")-> Metadata_grid_cam

write.table(Metadata_grid_cam, "data/data_clean/Metadata_grid_CAM.txt")
