
##### Load package #####

library(tidyverse)
library(readxl)

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
##### Load data #####

otu_plant = read.delim2("data/otu_plant.txt", header = T)%>%
  filter(grepl("CAM", Host_code))
str(otu_plant)
otu_plant%>%
  select_if(is.character)%>%
  filter(grepl("-",Iris.reichenbachiana))

otu_plant$Iris.reichenbachiana[grepl("-",otu_plant$Iris.reichenbachiana)]  = "20"
otu_plant$Schedonorus.arundinaceus = str_replace_all(otu_plant$Schedonorus.arundinaceus, ",", ".")
otu_plant$Schedonorus.arundinaceus[grepl(">",otu_plant$Schedonorus.arundinaceus)]  = "62.5"

otu_plant[,2:ncol(otu_plant)]= data.frame(lapply(otu_plant[,2:ncol(otu_plant)], as.numeric))
head(otu_plant)
otu_plant[,50:70]
otu_plant = otu_plant[,c(T,colSums(otu_plant[,colnames(otu_plant) != "Host_code"]) != 0) ]
otu_plant[grepl("21_CAM_15", otu_plant$Host_code)| grepl("21_CAM_14", otu_plant$Host_code),]
otu_plant$Host_code
##########Host_code


otu_plant_long = otu_plant%>%
  mutate(year = as.factor(substring(Host_code,1,2)),
         grid = as.factor(substring(Host_code,8,9)),
         quadra = as.factor(substring(Host_code,10,11)),
         total_cover = rowSums(across(where(is.numeric))),
         plant_richness = rowSums(across(where(is.numeric)) != 0)
         
  )%>%
  pivot_longer(cols = -c(year, grid, quadra, total_cover, plant_richness, Host_code), names_to = "Plant", values_to = "cover")  %>%
  group_by(year, grid, quadra)%>%
  mutate(
    relative_abundance = cover/total_cover,
    )%>%
  filter(cover != 0)

table(otu_plant_long$Plant)

ggplot(otu_plant_long)+
  geom_point(aes(total_cover, plant_richness))+
  main_theme
str(otu_plant_long)
otu_plant_long%>%
  group_by(year, grid)%>%
  summarise(grid_richness = n())

otu_plant_long%>%
  
  group_by(Plant) %>%
  summarise(n = n())%>%
  mutate(Freq = n()/sum(n))%>%
ggplot(aes(n))+
    geom_histogram(bins = 60)+
  labs(x = "nb. obs.", y = "count of species" )+
  main_theme

grid_abudance = otu_plant_long%>%
  group_by(year, grid)%>%
  summarise(plant_abundance_grid = n_distinct(Plant))%>%
  mutate(Grid_code = paste(year,"_CAM_",as.character(grid),sep = ""))
  
write.table(grid_abudance, "data/grid_abudance.txt")

str(grid_abudance)
ggplot()+
  geom_point(aes(cover_mean, n))
ChaoRichness(otu_plant_long)
ChaoRichness(spider$Girdled, datatype="abundance")
##### Indice d'abondance
library(iNEXT)

library(sjmisc)
otu_plant__pivot = otu_plant%>%
  column_to_rownames(var = "Host_code")%>%
  rotate_df()
output = iNEXT(x = otu_plant__pivot, q = c(0))