library(readxl)
library(tidyverse)
library(otuSummary)
library(ade4)       # used for plotting spatial data
library(spatialreg)
library(spdep)
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


read_xlsx("data/distance_quadrats.xlsx")%>%
  mutate(Quadrat = str_remove(Quadrat, "0"))-> dist_quad
xy <- as.matrix(dist_quad[, 2:3])
  
read.table("data/data_clean/Metadata_quadra_CAM.txt", header = T)%>%
  mutate(quadra = as.character(quadra))%>%
  left_join(dist_quad, by = join_by(quadra == Quadrat)) -> metadata_quadra

#read.table("data/data_clean/Metadata_grid_CAM.txt", header = T) -> metadata_grid
distmatrix = dist(xy, diag = F)
dist_df = matrixConvert(distmatrix, colname = c("Q", "QQ", "dist"))
unique(dist_df$dist)

metadata_quadra


div.knear4 <- knearneigh(x = xy,  # Coordinates (matrix form required)
                         k = 4) # Number of neighbors
distances <- nbdists(nb = knn2nb(div.knear4), coords = xy) 
invd1 <- lapply(distances, function(x) (1 / x))
pond4.stan <- nb2listw(neighbours = knn2nb(div.knear4),  # Enable compatibility by 
                       style = "B", # Binary values (0 / 1),
                       glist = invd1,
                       zero.policy = TRUE  # Set 0 when the observation has no neighbors
)

pond4.stan

str(metadata_quadra)
moran.test(x=metadata_quadra$Plant_richness, listw = pond4.stan,randomisation = T)
metadata_quadra%>%
  group_by(Grid_code)%>%
  summarise(Moran_index_Richness = moran.test(x=Plant_richness, listw = pond4.stan,randomisation = T)$estimate[1],
            Moran_expect = moran.test(x=Plant_richness, listw = pond4.stan,randomisation = T)$estimate[2],
            Moran_varariance = moran.test(x=Plant_richness, listw = pond4.stan,randomisation = T)$estimate[3],
            Moran_pval = moran.mc(x=Plant_richness, listw = pond4.stan, nsim = 1000)$p.value,
            significativity = ifelse(Moran_pval <0.05, "signif", "no signif")) -> moran_index_Richness

ggplot(moran_index_Richness)+
  geom_point(aes(Grid_code, Moran_index_Richness, col= significativity), cex =2)+
  geom_errorbar( aes(x= Grid_code,
                      ymin=Moran_index_Richness-sqrt(Moran_varariance), ymax=Moran_index_Richness+sqrt(Moran_varariance), col= significativity))+
  geom_point(aes(Grid_code, Moran_expect), col = "red", cex =2)+
  main_theme+
  theme(axis.text.x = element_text(angle = 90))

metadata_quadra%>%
  group_by(Grid_code)%>%
  filter(sum(Viral_richness)>=3)%>%
  summarise(Moran_index_Richness = moran.test(x=Viral_richness, listw = pond4.stan,randomisation = T)$estimate[1],
            Moran_expect = moran.test(x=Viral_richness, listw = pond4.stan,randomisation = T)$estimate[2],
            Moran_varariance = moran.test(x=Viral_richness, listw = pond4.stan,randomisation = T)$estimate[3],
            Moran_pval = moran.mc(x=Viral_richness, listw = pond4.stan, nsim = 1000)$p.value,
            significativity = ifelse(Moran_pval <0.05, "signif", "no signif")) -> moran_index_Richness_viral

ggplot(moran_index_Richness_viral)+
  geom_point(aes(Grid_code, Moran_index_Richness, col= significativity), cex =2)+
  geom_errorbar( aes(x= Grid_code,
                     ymin=Moran_index_Richness-sqrt(Moran_varariance), ymax=Moran_index_Richness+sqrt(Moran_varariance), col= significativity))+
  geom_point(aes(Grid_code, Moran_expect), col = "red", cex =2)+
  main_theme+
  theme(axis.text.x = element_text(angle = 90))

moran.test(x=metadata_quadra%>%filter(Grid_code == "20_CAM_02")%>%pull(Plant_richness), listw = pond4.stan,randomisation = T)
moran.mc(x=metadata_quadra%>%filter(Grid_code == "20_CAM_02")%>%pull(Plant_richness), listw = pond4.stan, nsim = 100)
cor4pp <-
  sp.correlogram(
    knn2nb(div.knear4),
    metadata_quadra%>%filter(Grid_code == "20_CAM_01")%>%pull(Plant_richness),
    order = 3,  # maximum of the scale (lag)
    method = "I",
    zero.policy = TRUE,
    style = "B"
  )
plot(cor4pp, main = "Moran autocorrelogram for bird richness (4 neighbors)")
