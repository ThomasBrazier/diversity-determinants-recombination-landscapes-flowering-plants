############################################################################
#                    ECOBIO - PhD
#
#       Genomic ressources for Genetic Maps - Database reading and management
#
############################################################################
# by Thomas Brazier
# brazier.thomas@gmail.com

# Supervisor: Sylvain GLEMIN
#             ECOBIO Lab


#==========================================================
# LOADING ENVIRONMENT
#==========================================================

# clear global environment: remove all variables
rm(list=ls(all=TRUE))

#----------------------
# Loading packages
library(rstudioapi)
library(ggplot2)
library(MareyMap)
library(gdata)
library(rentrez)

#======================================
# Loading variables & objects
#======================================

# Get the directory of the file & set working directory
filedir=dirname(rstudioapi::getSourceEditorContext()$path)
wd=filedir
setwd(wd)

source("sources/MareyMap.R")

#======================================
# Import & compile bibliographic database (DB)
#======================================
# Import the handmade bibliographic database in 'data'
# metadata = read.xls(xls = paste(wd, "/data/Genetic_maps/Genetic_maps_ressources.xlsx", sep = ""), sheet = 1, header = TRUE)
metadata = read.csv(file = paste(wd, "/data/Genetic_maps/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")

#-------------------------------------
# Compute the number of markers per genetic map
#-------------------------------------
# Return a data frame with the map name and the number of markers associated
metadata = read.csv(file = paste(wd, "/data/Genetic_maps/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")
nbMkr = get.nbMkr(path = "data/Genetic_maps/")
# Completing the database
for (i in 1:nrow(nbMkr)) {
  metadata$nb_markers[which(metadata$id == nbMkr$set[i])] = nbMkr$nbMkr[i]
}
# Saving the updated database
write.table(metadata, file = paste(wd, "/data/Genetic_maps/Genetic_maps_ressources.csv", sep = ""), col.names = TRUE, row.names = FALSE, sep = ";")

#-------------------------------------
# Linkage map length
#-------------------------------------
# Return a data frame with the map name and the linkage map length associated
metadata = read.csv(file = paste(wd, "/data/Genetic_maps/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")
linkageMap.length = get.maplength(path = "data/Genetic_maps/")
# Completing the database
for (i in 1:nrow(linkageMap.length)) {
  if (is.na(metadata$linkage_map_length[i])) {
    metadata$linkage_map_length[which(metadata$id == linkageMap.length$set[i])] = linkageMap.length$linkage_map_length[i]
  }
}
# Saving the updated database
write.table(metadata, file = paste(wd, "/data/Genetic_maps/Genetic_maps_ressources.csv", sep = ""), col.names = TRUE, row.names = FALSE, sep = ";")



#-------------------------------------
# Phylogenetic informations
#-------------------------------------
# Retrieve taxid, taxon (division) and taxon_group
metadata = read.table("data/Genetic_maps/Genetic_maps_ressources.csv", header = TRUE, sep = ";")
sp = metadata$species
taxid = character(rep(NA, length(sp)))
group = character(rep(NA, length(sp)))
for (i in 1:length(sp)) {
  # Taxid in ncbi
  r_search = entrez_search(db="taxonomy", term=sp[i], retmax = 1)
  r_search$ids
  if (length(r_search$ids) > 0) {
    taxize_summ = entrez_summary(db="taxonomy", id=r_search$ids)
    taxid[i] = taxize_summ$taxid
    group[i] = taxize_summ$division
  }
}
taxid
group
metadata$taxid = taxid
metadata$taxon_group = group
taxon = metadata$taxon
taxon[group %in% c("eudicots", "monocots", "flowering plants")] = "Angiosperm"
metadata$taxon = taxon
# Update database with results
write.table(metadata, "data/Genetic_maps/Genetic_maps_ressources.csv", col.names = TRUE,
            row.names = FALSE, quote = FALSE, sep = ";")

#======================================
# Genetic maps filtering
#======================================

# Select maps with a sufficient marker resolution to be intergrated in the final dataset


# Compile the DB into a clean DB .txt file in 'data-cleaned'
# Only maps with a complete Marey map

#======================================
# Database summary
#======================================
# Min number of markers
min(metadata$nb_markers, na.rm = TRUE)
# Max number of markers
max(metadata$nb_markers, na.rm = TRUE)
# Mean number of markers
mean(metadata$nb_markers, na.rm = TRUE)


# Number of species with a genetic map
length(unique(metadata$species[metadata$id != "Data Not Published"]))

#======================================
# END
#======================================