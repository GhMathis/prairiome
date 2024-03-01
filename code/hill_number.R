library(vegan)
read.table("data/data_clean/Metadata_quadra_CAM.txt", header = T) -> metadata_quadra 

read.table("data/data_clean/OTU_plant_CAM.txt", header = T)%>%
  # mutate( Grid_code = str_extract(Host_code, ".._CAM_.."))%>%
  # group_by(Grid_code)%>%
  # summarise_if(is.numeric, sum)%>%

  pivot_longer(-c(Host_code), names_to = "Plant", values_to = "occur")%>%
  pivot_wider(values_from = occur, names_from = Host_code )%>%
  column_to_rownames("Plant") -> OTU_plant_trspose

  mutate_all( ~.x/sum(.x) )
apply(abundance_plant_trspose, function(x) x * colSums(abundance_plant_trspose/100))
abundance_plant_trspose = as.data.frame(OTU_plant_trspose*2)
colSums(abundance_plant_trspose)
rbind(n = rep(200,ncol(abundance_plant_trspose)),abundance_plant_trspose) -> abundance_plant_trspose        
abundance_plant_trspose$`20_CAM_010200`
abundance_plant_trspose$`20_CAM_010100`
hills_numbers_plant_df = iNEXT(abundance_plant_trspose,q = 0,nboot = 50,
                               datatype="incidence_freq")
save(hills_numbers_plant_df,file = "data/data_clean/Inext_quadra.RData")
hills_numbers_plant_df$iNextEst$coverage_based%>%
  mutate(qD_rar = qD,
         qD_ext = qD, 
         qD.LCL_bot = qD.LCL,
         qD.UCL_bot = qD.UCL) -> tmp

hills_numbers_plant_df$iNextEst$coverage_based%>%
  dplyr::filter(Assemblage %in% c("20_CAM_010100", "20_CAM_010200", "20_CAM_010300", "20_CAM_010400"))%>%
  ggplot(., aes(x = SC, y = qD, group = Assemblage)) +
  geom_point(aes(col = Method), size = 4)+
  geom_ribbon(aes(ymin=qD.LCL, ymax=qD.UCL), alpha = 0.05, linetype = 0) +
  facet_wrap(~Assemblage)+
  ylab("Plant richness")+
  xlab("Coverage")+

  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15, face = "bold"),
    legend.title = element_text(size=12, face="bold"),
    legend.text = element_text(size=12, face="italic")
  )


data(spider)
out1 <- iNEXT(spider, q=c(0,1,2), datatype="abundance")
data(bird)
out2 <- iNEXT(bird, q=0, datatype="abundance")
colSums(bird)
