library(readxl)
library(tidyverse)
library(otuSummary)

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
distmatrix = dist(xy, diag = F)
dist_df = matrixConvert(distmatrix, colname = c("quadra_A", "quadra_B", "dist"))
unique(dist_df$dist)


read_xlsx("data/Metadata_Sample.xlsx") %>%
  dplyr::select(Host_code, Grid_code)%>%
  mutate(quadra_A = as.numeric(substr(Host_code, 10,11)))%>%
  group_by(Grid_code)%>%
  full_join(dist_df, by = join_by(quadra_A == quadra_A), relationship = "many-to-many" )-> dist_for_each_quadra
dist_for_each_quadra = na.omit(dist_for_each_quadra)

write.table(dist_for_each_quadra, "data/data_clean/dist_for_each_quadra.txt")

test = read.table("data/data_clean/dist_for_each_quadra.txt", header = T)
unique(test$dist)
