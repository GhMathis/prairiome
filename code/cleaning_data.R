library(tidyverse)
library(readxl)
library(iNEXT)

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

##### OTU Camargue

### PLANT
read_xlsx("data/OTU_plant.xlsx")%>%
  filter(str_detect( Host_code, "CAM"),
         !str_detect(Host_code, "21_CAM_1[45]")) %>%
  replace(is.na(.),0)%>%
  select(-Richness)%>%
  select(Host_code, where(~is.numeric(.x) && sum(.x) != 0 )) -> otu_plant_cam

write.table(otu_plant_cam, "data/data_clean/OTU_plant_CAM.txt",row.names = F)

otu_plant_cam%>%
  mutate( Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
  group_by(Grid_code) %>%
  mutate(across(where(is.numeric), ~.*2))%>%
  summarise_if(is.numeric, sum)-> abund_plant_grid 


write.table(abund_plant_grid, "data/data_clean/abund_plant_grid.txt",row.names = F)

### VIRUS
read.delim2("data/S1_Viral_OTU.txt",header = T)%>%
  rename(Host_code = "X")%>%
  filter(str_detect( Host_code, "CAM"),
         !str_detect(Host_code, "21_CAM_1[45]")) %>%
  select(Host_code, where(~is.numeric(.x) && sum(.x) != 0 ))-> otu_virus_cam

write.table(otu_virus_cam, "data/data_clean/OTU_virus_CAM.txt",row.names = F)


##### Metadata Quadra 
read_xlsx("data/Metadata_Sample.xlsx")%>%
  filter(str_detect(Locality, "Arles"))%>%
  select(-c(Rocks, Additional_information))-> metadata_quad_cam
str(metadata_quad_cam)

read.table("data/data_clean/OTU_plant_CAM.txt",header = T)%>%
  reframe(Host_code = Host_code,
         year = as.factor(substring(Host_code,1,2)),
         grid = as.factor(substring(Host_code,8,9)),
         quadra = as.factor(substring(Host_code,10,11)),
         total_cover = rowSums(across(where(is.numeric))),
         Plant_richness = rowSums(across(where(is.numeric), ~.x != 0))
         
  )%>%
  select(Plant_richness,year,grid,quadra, total_cover,Host_code)%>%
  left_join(metadata_quad_cam, by ="Host_code")%>%
  mutate(across(c(Vegetation, Litter, Soil),~as.numeric(.x)),
         Soil = case_when(str_detect(Host_code, "20_CAM_01")~ 0,
         .default = Soil),
         Free_space = Litter +Soil) -> metadata_quad_cam

read.table("data/data_clean/OTU_virus_CAM.txt",header = T)%>%
  reframe(Host_code = Host_code,
          Viral_richness = rowSums(across(where(is.numeric), ~.x != 0))
          )%>%
  left_join(metadata_quad_cam, by ="Host_code")-> metadata_quad_cam

write.table(metadata_quad_cam, "data/data_clean/Metadata_quadra_CAM.txt")

##### Metadata grid
otu_plant_cam%>%
  pivot_longer(-Host_code, names_to = "Plant", values_to = "cover")%>%
  mutate(Grid_code = substring(Host_code, 1,9))%>%
  filter(cover != 0) -> metadatat_grid_plant

read_xlsx("data/Habitats_EDGG_patures.xlsx")%>%
  select(Grid_code = Sampling_n, Pature,Fauche)-> paturage


read_xlsx("data/donnee_sol.xlsx",col_names = F, skip = 3)%>%
  rename_with(~c("Grid_code", "Num", "pHwater", "lime_tot", "MO", "Phos", "K", "Mg",
                 "Ca", "Na", "N","C","CN","clay",  "SiltF", "SiltC","SandF", "SandC",
                 "Cl", "Res", "Cond"))%>%
  filter(str_detect(Grid_code,"CAM" ))%>%
  mutate(Grid_code = str_extract(Grid_code,".._CAM_.." ))%>%
  left_join(y = read_xlsx("data/distance_fer_sol_EDGG_CAM.xlsx"), by = join_by(Grid_code == EDGG))%>%
  mutate(dep_oxy_num = as.numeric(str_remove(Profondeur, "[>]")),
         depth_oxy = cut(dep_oxy_num,
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
  left_join(soil_data, by= "Grid_code")%>%
  left_join(y = read.table("data/data_clean/Buffer_5ldscp_CAM.txt"), by = "Grid_code")%>%
  left_join(paturage, by = "Grid_code") -> Metadata_grid_cam

#####Abundance indices at grid levels

rbind(t(data.frame(sample_size = rep(9,42))),.)

read.table("data/data_clean/OTU_plant_CAM.txt")%>%
  mutate( Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
  group_by(Grid_code) %>%
  mutate(across(where(is.numeric), ~ifelse(. != 0,1,0 ))) %>%
  summarise_if(is.numeric, sum)%>%
  pivot_longer(-c(Grid_code), names_to = "Plant", values_to = "occur")%>%
  pivot_wider(values_from = occur, names_from = Grid_code )%>%
  column_to_rownames("Plant")-> OUT_plant_grid_binary_df

rbind(n = rep(9,42),OUT_plant_grid_binary_df) -> OUT_plant_grid_binary_df

read.table("data/data_clean/OTU_virus_CAM.txt")%>%
  mutate( Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
  group_by(Grid_code) %>%
  mutate(across(where(is.numeric), ~ifelse(. != 0,1,0 ))) %>%
  summarise_if(is.numeric, sum)%>%
  pivot_longer(-Grid_code, names_to = "Virus", values_to = "occur")%>%
  group_by(Grid_code)%>%
  filter(sum(occur) >=2 )%>%
  pivot_wider(values_from = occur, names_from = Grid_code )%>%
  column_to_rownames("Virus") -> OUT_virus_grid_binary_df

rbind(n = rep(9,ncol(OUT_virus_grid_binary_df)),OUT_virus_grid_binary_df) -> OUT_virus_grid_binary_df
OUT_virus_grid_binary_df$`21_CAM_16`
OUT_virus_grid_binary_df$`21_CAM_17`
colSums(OUT_virus_grid_binary_df)
# 
# otu_plant_cam%>%
#   pivot_longer(-Host_code, names_to = "Plant", values_to = "cover")  %>%
#   pivot_wider(values_from = cover, names_from = Host_code   , values_fill = 0) -> otu_plant_cam_trpose
# otu_plant_cam_trpose = as.data.frame(otu_plant_cam_trpose)
# rownames(otu_plant_cam_trpose) = otu_plant_cam_trpose$Plant
# as.data.frame(otu_plant_cam_trpose)%>%
#   dplyr::select(-Plant) -> otu_plant_cam_trpose
# otu_plant_cam_binary = otu_plant_cam_trpose
# otu_plant_cam_binary[otu_plant_cam_trpose != 0] = 1 
OUT_plant_grid_binary_df
if(file.exists("data/data_clean/Inext_plant_grid.RData")){
  load("data/data_clean/Inext_plant_grid.RData")
}else{
  hills_numbers_plant_df = iNEXT(OUT_plant_grid_binary_df, q = 0,nboot = 140,
                                 datatype="incidence_freq")
  save(hills_numbers_plant_df, file = "data/data_clean/Inext_plant_grid.RData")
}


hills_numbers_virus_df = iNEXT(OUT_virus_grid_binary_df,q = 1, nboot = 120,
                               datatype="incidence_freq")
hills_numbers_plant_df%>%
  ggiNEXT(., type = 1, color.var = "Order.q")+
  facet_wrap(~Assemblage)+
  main_theme

colSums(OUT_virus_grid_binary_df)
hills_numbers_virus_df%>%
  ggiNEXT(., type = 1, color.var = "Order.q")+
  facet_wrap(~Assemblage)+
  main_theme
hills_numbers_virus_df%>%
  ggiNEXT(., type = 3, color.var = "Order.q")+
  facet_wrap(~Assemblage)+
  main_theme

iNEXT(rowSums(OUT_virus_grid_binary_df),q = 0, nboot = 100,
      datatype="incidence_freq")%>%
  ggiNEXT(., type = 1, color.var = "Order.q")+
  ylab("Global Viral richness")+
  xlab("Sample size")+
  main_theme

temp_plant = iNEXT(rowSums(OUT_plant_grid_binary_df),q = 0, nboot = 100,
      datatype="incidence_freq")
temp_plant%>%  
ggiNEXT(., type = 3, color.var = "Order.q")+
  ylab("Global Plant richness")+
  xlab("Sample size")+
  main_theme

temp = iNEXT(rowSums(OUT_virus_grid_binary_df),q = 0, nboot = 100,
      datatype="incidence_freq")
temp%>%  
ggiNEXT(., type = 3, color.var = "Order.q")+
  ylab("Global Plant richness")+
  xlab("Sample size")+
  main_theme
# iNEXT(OUT_virus_grid_binary_list[[38]],q = 0, se=TRUE, size = c(1, 10, 20, 30, 40, 50, 80, 100), nboot = 100) %>% 
#   ggiNEXT(., type = 1, color.var = "Order.q")
# ChaoRichness(temp, datatype="incidence_raw")
# 
# # Loop through sample names and extract corresponding columns into list
# otu_virus_cam
# for(sample_name in unique(Metadata_grid_cam$Grid_code)[1]) {
#   ### Plant
#   temp <- otu_plant_cam_binary[, str_detect(names(otu_plant_cam_binary), sample_name)]
#   Plant_richness_grid = c(Plant_richness_grid , ChaoRichness(temp, datatype="incidence_raw")[[1]])
#   Plant_shannon_grid = c(Plant_shannon_grid, ChaoShannon(temp, datatype="incidence_raw")[[1]])
#   Plant_simpson_grid = c(Plant_simpson_grid, ChaoSimpson(temp, datatype="incidence_raw")[[1]])
#   
#   ### Virus
#   temp <- otu_virus_cam_trpose[, str_detect(names(otu_virus_cam_trpose), sample_name)]
#   if(sum(temp) != 0){
#    
#     Virus_richness_grid = c(Virus_richness_grid , ChaoRichness(temp, datatype="incidence_raw")[[1]])
#     Virus_shannon_grid = c(Virus_shannon_grid, ChaoShannon(temp, datatype="incidence_raw")[[1]])
#     Virus_simpson_grid = c(Virus_simpson_grid, ChaoSimpson(temp, datatype="incidence_raw")[[1]])
#   }else{
#     Virus_richness_grid = c(Virus_richness_grid , 0)
#     Virus_shannon_grid = c(Virus_shannon_grid, 0)
#     Virus_simpson_grid = c(Virus_simpson_grid, 0)
#   }
#   
# }
# 
Metadata_grid_cam%>%
  mutate(Plant_richness_grid = OUT_plant_grid_binary_df,
         Plant_shannon_grid = Plant_shannon_grid,
         Plant_simpson_grid = Plant_simpson_grid,
         Virus_richness_grid = Virus_richness_grid,
         Virus_shannon_grid = Virus_shannon_grid,
         Virus_simpson_grid = Virus_simpson_grid) -> Metadata_grid_cam

write.table(Metadata_grid_cam, "data/data_clean/Metadata_grid_CAM.txt")
