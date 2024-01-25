
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
                     include.lowest = TRUE)) -> Meta_data_ground 

read.table("data/Metadata_Grid_CAM.txt", header = T)%>%
  left_join(Meta_data_ground, by = "Grid_code")%>%
  write.table(file = "data/Metadata_Grid_CAM.txt")

