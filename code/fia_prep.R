

library(data.table)
library(dplyr)
library(tidyr)
library(raster)
library(ggplot2)
select <- dplyr::select

#### load, filter, and join FIA data

tree <- fread("E:/fia/raw_data/ENTIRE/TREE.csv", stringsAsFactors=F)
tree <- select(tree, PLT_CN, PLOT, SUBP, TREE, SPCD, DIA, HT)
gc()

plot <- fread("E:/fia/raw_data/ENTIRE/PLOT.csv", stringsAsFactors=F)
subplot <- fread("E:/fia/raw_data/ENTIRE/SUBPLOT.csv", stringsAsFactors=F)
species <- fread("E:/fia/raw_data/REF_SPECIES.csv", stringsAsFactors=F)

plot <- select(plot, CN, PLOT, LAT, LON, ELEV) %>%
      mutate(PLT_CN=CN) %>% select(-CN)
subplot <- select(subplot, PLT_CN, SUBP, SLOPE, ASPECT)
species <- select(species, SPCD, COMMON_NAME, GENUS, SPECIES)
gc()

d <- left_join(tree, subplot) %>%
      left_join(plot) %>%
      left_join(species)
colnames(d) <- tolower(colnames(d))

rm(tree)
rm(plot)
rm(subplot)
rm(species)
gc()

# correct outdated nomenclature
d$genus[d$genus=="Lithocarpus"] <- "Notholithocarpus"
d$gs <- paste(d$genus, d$species)

# filter to focal species list
spp <- c("Acer macrophyllum", "Aesculus californica", "Arbutus menziesii", 
         "Pseudotsuga menziesii", "Notholithocarpus densiflorus", "Quercus agrifolia", 
         "Quercus douglasii", "Quercus garryana","Quercus kelloggii", 
         "Quercus lobata", "Sequoia sempervirens", "Umbellularia californica")
d <- filter(d, gs %in% spp)

# export
write.csv(d, "data/fia.csv", row.names=F)



