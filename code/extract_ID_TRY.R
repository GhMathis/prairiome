library(readxl)
library(tidyverse)
try_data_base = read.delim("data/TryAccSpecies.txt", header = T)
tax_plant =read_xlsx("data/TAX_Plant.xlsx")
try_data_base
str(try_data_base)
tax_plant$Plant_species
tax_plant %>%
  mutate(Sp_names = str_extract(Plant_species, "^\\w+\\s\\w+")) -> tax_plant

tibble(try_data_base)%>%
  filter(AccSpeciesName %in% tax_plant$Plant_species)
tibble(try_data_base)%>%
  filter(AccSpeciesName %in% tax_plant$Sp_names) -> identified_sp

write.table(identified_sp$AccSpeciesID, row.names = F, col.names = F, file = "data/IDsp_TRY_EDGG.txt",  eol =", " )

