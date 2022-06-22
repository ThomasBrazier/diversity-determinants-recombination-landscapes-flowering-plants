############################################################################
#                    ECOBIO - PhD
#
#       Generate and assemble physical maps for each reference genome
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

#======================================
# Loading variables & objects
#======================================

# Get the directory of the file & set working directory
wd=dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
# Custom functions
source("sources/PhysMap.R") # Get the physical positions of a marker given a list of marker names
source("sources/MareyMap.R") # Consolidate a MareyMap input file
source("sources/core_pipeline.R") # Consolidate a MareyMap input file

#########-------------------------------##########
# FROM THIS POINT, READ & SAVE DATA IN 'data'
#########-------------------------------##########

#======================================
# Loading variables & objects
#======================================
# Display the list of maps
list = system(paste("ls ", wd, "/data/Genetic_maps/", sep = ""), intern = TRUE)
list = list[grep("*.txt", list)]
# How many maps?
# length(list)
list

# Select the data you which to work on
set = list[21]
set
species = regmatches(set, regexpr("^[A-Z][a-z]+_[a-z]+", as.character(set)))

#======================================
# Blast for marker position in a marker dataset
#======================================
# install.packages("remotes")
# remotes::install_github("mhahsler/rBLAST")
# https://rdrr.io/github/mhahsler/rBLAST/man/BLAST.html

# "Malus_domestica", "Rosa_chinensis" à refaire en améliorant détection des chromosomes
# species = "Sorghum_bicolor"
# species = "Hordeum_vulgare"
# species = "Elaeis_guineensis" # Not working well, no alignment
# Load the marker database to blast
mkrDB = read.csv(paste("data/Physical_maps/Markers/", species, ".csv", sep = ""), sep = ";")
accession = levels(mkrDB$accession)[1]
# accession = "GCA_014319735.1"

#--------------------------------------
# If a marker file contains too much useless markers
# e.g. with Prunus avium and 143223 markers with sequences to blast
# Hence reduce the marker dataset to the markers already present in the genetic map to consolidate
# mkrDB = read.csv(paste("data/Physical_maps/Markers/", species, ".csv", sep = ""), sep = ";")
# genmap = read.table(paste("data/Genetic_maps/", set, sep = ""), header = TRUE)
# mkrDB$mkr = as.character(mkrDB$mkr)
# genmap$Marker = as.character(genmap$Marker)
# # Number of common markers
# sum(mkrDB$mkr %in% genmap$Marker)
# which(mkrDB$mkr %in% genmap$Marker)
# # mkrDB$mkr[545]
# # genmap$Marker[genmap$Marker == mkrDB$mkr[545]]
# # Trim mkrDB to common markers only and save the new file
# mkrDB = mkrDB[which(mkrDB$mkr %in% genmap$Marker),]
# write.table(mkrDB, file = paste("data/Physical_maps/Markers/",species ,".csv", sep = ""), col.names = TRUE, row.names = FALSE, sep = ";")
# 

#--------------------------------------
# Blast markers
# Given a marker database, the function retrieve the position of the markers
# An accession to a reference genome must be provided to narrow the search
# Return a marker database updated with positions
# tol: a tolerance parameter, the minimum percentage of identity with the reference genome accepted to accept a marker position
source("sources/PhysMap.R")
new_positions = blast.markers(species = species, accession = accession, db = mkrDB, tol = 90,
                              update.all = FALSE)


# Diagnostics
sum(is.na(new_positions$phys)) # Number of NA (markers not found with Blast)
# sum(is.na(new_positions$map)) # Number of NA (markers not found with Blast)
table(new_positions$map) # Number of markers per chromosome
# View(new_positions)
# If you retrieved a sufficient number of markers with a high confidence, then
# save your marker dataset
write.table(new_positions, file = paste("data/Physical_maps/Markers/",species ,".csv", sep = ""), col.names = TRUE, row.names = FALSE, sep = ";")



#--------------------------------------
# Blast scaffolds
source("sources/PhysMap.R")
new_positions = blast.scaffolds(species, update = TRUE, tol = 95, submitter = "tbrazier")
# Diagnostics
sum(is.na(new_positions$phys)) # Number of NA (markers not found with Blast)
table(as.numeric(new_positions$map)) # Number of markers per chromosome
# If you retrieved a sufficient number of markers with a high confidence, then
# save your marker dataset
write.table(new_positions, file = paste("data/Physical_maps/Markers/",species ,".csv", sep = ""), col.names = TRUE, row.names = FALSE, sep = ";")



#======================================
#   END
#======================================
