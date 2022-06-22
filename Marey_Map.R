########################################################################### #
#                    ECOBIO - PhD
#
#       Creating recombination landscapes with Marey Map
#
########################################################################### #
# by Thomas Brazier
# brazier.thomas@gmail.com

# Supervisor: Sylvain GLEMIN
#             ECOBIO Lab


#============================================================================#
# LOADING ENVIRONMENT ----
#============================================================================#

# clear global environment: remove all variables
rm(list=ls(all=TRUE))

# Loading packages
library(rstudioapi)
library(ggplot2)
library(MareyMap)
library(stringr)
library(myTAI)

#============================================================================#
# Loading variables & objects ----
#============================================================================#

# Get the directory of the file & set working directory
wd=dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

# Custom functions
source("sources/PhysMap.R") # Get the physical positions of a marker given a list of marker names
source("sources/MareyMap.R") # Consolidate a MareyMap input file
source("sources/core_pipeline.R") # Basic functions but essential for the pipeline (e.g. saving logs)
source("sources/stats.R") # Statistical functions implemented


#########-------------------------------######## ##
# FROM THIS POINT, READ & SAVE DATA IN 'data' ----
#########-------------------------------######## ##

# At each step, the variable 'set' define which genetic/physical/Marey map is processed
# hence it is important to be careful to your definition of 'set' all along the pipeline

# Display the list of maps
list = system(paste("ls ", wd, "/data/Genetic_maps/", sep = ""), intern = TRUE)
list = list[grep("*.txt", list)]
# How many maps?
length(list)
list

# Select the data you which to work on
set = list[36]
set
# Beginning with existing reference genome, on which we tried to apply genetic maps...

#============================================================================#
# Build input files for Marey Map - Marey map assembly ----
#============================================================================#

# Need a genetic and a physical maps sharing the same markers
# https://github.com/aursiber/MareyMap
# So, when we have access to the genetic map and physical maps separately (e.g. markers need to be blasted to a reference genome)
# we then need to consolidate both maps to a single Marey map dataset

# Map a physical map onto a genetic map
# For use only when data does not provide both genetic and physical maps
# Access to the markers database for the species and retrieve physical positions and chromosome number of the genetic map markers
# Otherwise, fail to recover markers position and give a warning
consolidateMareyMap(set = set, markerDB = "data/Physical_maps/Markers/",
                    input = "data/Genetic_maps/", output = "data/Marey_maps/")


#============================================================================#
# Map quality control & Data cleaning ----
#============================================================================#

# Diagnostic plots
# Genetic distance (cM) as a function of physical position (bp)
# Make a diagnostic plot of all maps in 'data/Marey_maps/'
# List of all maps
list = system(paste("ls ", wd, "/data/Marey_maps/", sep = ""), intern = TRUE)
list = list[grep("*.txt", list)]
# length(list)
list
set = list[107]
set
# Save diagnostic plots as one figure per dataset in 'output/marey_maps/diagnostic_plots/'
# Valid markers are in blue and unvalid markers are in red
# list = list[14]
# the whole list
# marey.diagplot(set = list, dir = "data/Marey_maps/", output = "output/marey_maps/diagnostic_plots/", display = TRUE)
# Or a single dataset
marey.diagplot(set = set, dir = "data/Marey_maps/", output = "output/marey_maps/diagnostic_plots/", display = TRUE)
# rm(list)

#-------------------------------------#
# Manual curating - verify if dataset open in MareyMap
# Open the GUI with the command:
startMareyMapGUI() # Long & tiedious; no automation of opening procedures

# Implement a graphic interface for fast point & click to unvalidate outliers & aberrant markers
# One click to unvalidate/validate
# + Zone selection

#-------------------------------------#
# Check & correct bad orientation of linkage groups
# Markers oriented in the direction of a reference map if necessary
# A reference genome with markers mapped onto it seems inevitable

# FULL MAP CORRECTION
# Select the data you wish to work on; ONE map ONLY
set = list[67]
set
# Plot the data
marey.diagplot(set = set, dir = "data/Marey_maps/", save.plot = TRUE, display = TRUE)
# Apply correction (automatic check not implemented yet)
# Set = name of a map
# id_chr = index of chromosomes to flip
new_map = check_map_orientation(set, id_chr = c(10,11,13,14,2,4,6,8))
# Plot the data to validate the procedure
marey.diagplot(set = new_map)
# Save the corrected map once validated
write.table(new_map, file = paste("data/Marey_maps/", set, sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
save.log(msg = paste("Check and correct bad orientation of linkage groups: Map saved.", sep = " "))


# CORRECTION OF A SEGMENT WITHIN A MAP
source("sources/MareyMap.R")
# Select the data you wish to work on; ONE map ONLY
set = list[13]
set
# Plot the data
marey.diagplot(set = set, dir = "data/Marey_maps/", save.plot = FALSE, display = TRUE)
# Apply correction (automatic check not implemented yet)
# Set = name of a map
# Work on a given chromosome for a given dataset
chr = "8"
# id_chr = index of chromosomes to flip
new_map = flip.segment(set, chr)
# Plot the data to validate the procedure
marey.diagplot(set = new_map, save.plot = FALSE)
# Save the corrected map once validated
write.table(new_map, file = paste("data/Marey_maps/", set, sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
save.log(msg = paste("Check and correct bad orientation of a segment in a chromosome: Map saved.", sep = " "))


#-------------------------------------#
# Filter out outliers automatically
# Selection of valid markers

######### STEP 1 ######### 

# A rough filtering of clear outliers, i.e. points clearly aberrant
# Before noise reduction

# TODO How to implement an automatic procedure for filtering outliers?
# 1/ Fit a curve and filter markers away from the curve (threshold): moving average (https://en.wikipedia.org/wiki/Moving_average) acting as a low-pass filter
# bootstrap the dist. of the mean departure from the curve, exclude points far away than 95% (2.5% each side)
# Advantage, selection is curve specific
# Yet, big problem: interpolation is also a curve fitting, circular approach !!
# 2/ Sample the subset of markers minimizing the RMSE (resampling approach)

# Plot the data
marey.diagplot(set = set, dir = "data/Marey_maps/",save.plot = FALSE, display = TRUE)
# Apply automatic check and correction
# TODO Yet to do, may be not necessary at this time
# Plot the data to validate the procedure
marey.diagplot(set = set, dir = "data/Marey_maps/",save.plot = FALSE, display = TRUE)


#-------------------------------------#
# Noise estimation and reduction
# See signal processing and information theory
# Keywords: noise reduction

# In this part, we can consider that
# (a) the noise of the measure of recombination distance (cM) explains the non-alignement of markers
# (b) quality of genome assembly (e.g. gaps) can influence the quality of the Marey map
# So we decided to treat this as a problem of signal processing/information theory
# Hence the noise reduction step for Marey curves that are not monotoneously perfectly aligned markers...


#-------------------------------------#
# Filter out outliers manually
# Visually check for outliers in the Marey map plot
source("sources/MareyMap.R")
# Select the data you wish to work on; ONE map ONLY
set = list[19]
set
# Plot the data
marey.diagplot(set = set, dir = "data/Marey_maps/", save.plot = FALSE, display = TRUE)
# Apply correction (automatic check not implemented yet)
# Set = name of a map
# id_chr = index of chromosomes to flip
chr = "15"
new_map = outliers.selection(set , chr)
# Plot the data to validate the procedure
marey.diagplot(set = new_map, save.plot = FALSE, display = TRUE)
# Save the corrected map once validated
write.table(new_map, file = paste("data/Marey_maps/", set, sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
save.log(msg = paste("Manually removed some outliers in a chromosome: Map saved.", sep = " "))

# Discard a complete chromosome
chr = "6"
new_map = read.table(file = paste("data/Marey_maps/", set, sep = ""), header = TRUE)
new_map$vld[new_map$map == chr] = FALSE
write.table(new_map, file = paste("data/Marey_maps/", set, sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
save.log(msg = paste("Chromosome", chr,"in", set,"discarded: Map saved.", sep = " "))


######### STEP 2 ######### 

# After noise reduction, apply a second step of a more subtle filtering, trying to keep only markers with a high confidence degree

#-------------------------------------#
# Filter out chromosomes with insufficient number of markers,
# based on an interval threshold or number of markers specified by the user
# mean interval distance
# Plot the data
marey.diagplot(set = set, dir = "data/Marey_maps/",save.plot = FALSE, display = TRUE)
# Apply automatic check and correction
# TODO, manual curating for the time
# Plot the data to validate the procedure
marey.diagplot(set = set, dir = "data/Marey_maps/",save.plot = FALSE, display = TRUE)



#============================================================================#
# Consensus map, if necessary ----
#============================================================================#

# Consensus maps is a way to rpoduce highly saturated linkage maps from multiple maps with gaps
# Yet, it has already be done in model species
# And we lack data for non-model species
# Mays be this aprt involves too much work for the small benefice expected...

# Merge multiple maps to obtain a high-density consensus map
## Step by step procedure inspired by Mapfuser (Van Muijen et al. 2017)
# The LPmerge package assumes sufficient marker overlap and identical orientation of linkage groups
#-------------------------------------#
# 1/ Correct map orientation



#-------------------------------------#
# 2/ Check anchor markers



#-------------------------------------#
# 3/ Evaluate overlap between maps in network visualisation
# igraph R package (Csardi and Nepusz, 2006)
# GeneticMapComparator (Holtz et al., 2017) for map visualization


#-------------------------------------#
# Merging maps with LPmerge (Endelman and Plomion, 2014) or LCSLCIS (De Matteo et al., 2018)


#-------------------------------------#
# Visual assessment of the Marey map, in order to manually exclude outliers and validate the map at the end of the quality control procedure
startMareyMapGUI()
# Reformat output from MareyMapGUI
system(paste("sed -i '' 's/\"//g' ", wd, "/data/Marey_maps/*.txt", sep = "")) # Remove "
system(paste("sed -i '' 's/ /\t/g' ", wd, "/data/Marey_maps/*.txt", sep = "")) # Change spaces to tab

# Ultimately save cleaned maps in '/data-cleaned/Marey_maps/'



#============================================================================#
# MAP METADATA ----
#============================================================================#
# Metadata are associated to each map dataset, identified by 'Genus_species_AuthorYear'
# Some metadata are (1) gathered in associated litterature (e.g. crossing, number of progenies, genome accession)
# While other are (2) specific to species and come from a 'Species_metadata.txt' file
# Finally, some metadata are associated to the dataset (e.g. number of markers)

# Metadata (1) are written directly into 'Genetic_maps_ressources.csv'

# Metadata (2) are written in a separate file 'Species_metadata.txt' then propagated to the 'Genetic_maps_ressources.csv'
# with following columns:
# species
# taxid
# taxon
# taxon_group
# family
# genome_size
# genome_Cvalue
# ploidy
# hcn
metadata = read.csv(file = paste(wd, "/data/Genetic_maps/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")
df = data.frame(species = sort(as.character(unique(metadata$species))),
  taxid = NA,
  #taxon = NA,
  taxon_group = NA,
  angiosperm_group = NA,
  family = NA,
  genome_size = NA,
  genome_Cvalue = NA,
  ploidy = NA,
  hcn = NA)
# Retrieve taxid and phylogenetic infos with R and myTAI
# taxonomy() retrieves taxonomic information from NCBI taxonomy
# taxo = taxonomy(organism = "Zea_mays", 
#                 db = "ncbi",
#                 output = "classification")
for (i in 1:nrow(df)) {
  cat(as.character(df$species[i]), "\n")
  # Since ncbi is limiting the number of requests to 3 per seconds
  Sys.sleep(0.7)
  taxo = taxonomy(organism = gsub("_", " ", df$species[i]), 
                  db = "ncbi",
                  output = "classification")
  if (length(taxo) > 1) {
    df$taxid[i] = taxo$id[which(taxo$rank == "species")]
    df$taxon_group[i] = taxo$name[12] # Mono- or Eudicotyledons
    df$family[i] = taxo$name[which(taxo$rank == "family")]
  }
}
# Cvalue database used to query genome informations
cvalue_db = read.csv(file = paste(wd, "/data-cleaned/CValueDataBase.csv", sep = ""), header = TRUE, sep = ";")
cvalue_db$sp = paste(cvalue_db$Genus, cvalue_db$Species, sep = "_")
for (i in 1:nrow(df)) {
  if (df$species[i] %in% cvalue_db$sp) {
    df$angiosperm_group[i] = as.character(cvalue_db$Angiosperm.Group[which(cvalue_db$sp == df$species[i])])
    # df$genome_size[i] = NA
    df$genome_Cvalue[i] = cvalue_db$Cval[which(cvalue_db$sp == df$species[i])]
    df$ploidy[i] = cvalue_db$PloidyLevel[which(cvalue_db$sp == df$species[i])]
    # Haploid chromosome number is the diploid chromosome numbers divided by ploidy level
    if (!is.na(cvalue_db$ChromosomeNumber[which(cvalue_db$sp == df$species[i])])) {
      df$hcn[i] = as.integer(min(cvalue_db$ChromosomeNumber[which(cvalue_db$sp == df$species[i])], na.rm = TRUE)/
        df$ploidy[i])
    } else {
      df$hcn[i] = NA
    }
  }
}

write.table(df, file = "data-cleaned/species_metadata.csv",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep=";")
df = read.table("data-cleaned/species_metadata.csv", header = TRUE, sep = ";")

# Update 'Genetic_maps_ressources.csv' with 'Species_metadata.txt'
df$species = as.character(df$species)
metadata$species = as.character(metadata$species)
metadata$taxon_group = as.character(metadata$taxon_group)
for (i in 1:nrow(metadata)) {
  if (metadata$species[i] %in% df$species) {
    metadata$taxid[i] = df$taxid[which(df$species == metadata$species[i])]
    metadata$taxon_group[i] = df$taxon_group[which(df$species == metadata$species[i])]
    metadata$family[i] = df$family[which(df$species == metadata$species[i])]
    metadata$genome_Cvalue[i] = df$genome_Cvalue[which(df$species == metadata$species[i])]
    metadata$ploidy[i] = df$ploidy[which(df$species == metadata$species[i])]
    metadata$hcn[i] = df$hcn[which(df$species == metadata$species[i])]
  }
}
write.table(metadata, file = "data/Genetic_maps/Genetic_maps_ressources.csv",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep=";")


# Estimate (3), based on data
# Display the list of maps
list = system(paste("ls ", wd, "/data/Marey_maps/", sep = ""), intern = TRUE)
list = list[grep("*.txt", list)]
# How many maps?
list = gsub(".txt", "", list)
length(list)
list

df = data.frame(set = as.character(list),
                nb_markers = NA,
                linkage_map_length = NA,
                physical_length = NA)
df$set = as.character(df$set)
for (i in 1:nrow(df)) {
  cat(as.character(df$set[i]), "\n")
  data = read.table(paste("data/Marey_maps/", df$set[i], ".txt", sep = ""), header = TRUE)
  # Number or markers per dataset
  # is the number of genetic positions not NA
  df$nb_markers[i] = sum(!is.na(data$gen))
  # Total linkage map length
  # the sum of all chromosome linkage length
  data$gen = as.numeric(data$gen)
  subset = split(data$gen, data$map)
  subset = subset[!(unlist(lapply(subset, function(x) max(x, na.rm = TRUE)))) == -Inf]
  df$linkage_map_length[i] = sum((unlist(lapply(subset, function(x) max(x, na.rm = TRUE)))), na.rm = TRUE)
  # Total genome size
  data$phys = as.numeric(data$phys)
  subset = split(data$phys, data$map)
  subset = subset[!(unlist(lapply(subset, function(x) max(x, na.rm = TRUE)))) == -Inf]
  df$physical_length[i] = sum((unlist(lapply(subset, function(x) max(x, na.rm = TRUE)))))
}


# Save the updated 'genetic_maps_ressources'
metadata = read.csv(file = paste(wd, "/data/Genetic_maps/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")
for (i in 1:nrow(df)) {
  cat(df$set[i], "\n")
  idx = which(metadata$id == df$set[i])
  metadata$nb_markers[idx] = df$nb_markers[i]
  metadata$linkage_map_length[idx] = df$linkage_map_length[i]
  metadata$physical_length[idx] = df$physical_length[i]
}
write.table(metadata, file = "data/Genetic_maps/Genetic_maps_ressources.csv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep=";")




#============================================================================#
# Map filtering - Save in data-cleaned ----
#============================================================================#
# Auto-update data-cleaned Marey maps and associated metadata
source("sources/MareyMap.R")
update_datacleaned()

# Before continuing, we need to validate our dataset and to discard some chromosomes that may be of poor quality
# A choice based on visual assessment of Marey maps produced and
# summary statistics (e.g. number of markers) to help the decision
# Load the joint database with all maps
source("sources/MareyMap.R")
data.final = read.table(file = "data-cleaned/marey_maps/AllMaps.txt", header = TRUE, sep = "\t")

map_stats = map.statistics(data.final)

# Initiate a vectrr of indexes of maps to discard
discard_idx = NA
# Be very conservative at this step, since interpolation quality relies on data quality

#-------------------------------------#
# Filtering by number of markers in the map
# Span calibration and interpolation fitting are very sensitive to low number of markers
# Hence use a minimal number of markers
summary(map_stats$nb.markers)
# How many under 40 markers?
sum(map_stats$nb.markers < 40)
paste(map_stats$set[which(map_stats$nb.markers < 40)], map_stats$chromosome[which(map_stats$nb.markers < 40)])
# How many under 30 markers?
sum(map_stats$nb.markers < 30)
paste(map_stats$set[which(map_stats$nb.markers < 30)], map_stats$chromosome[which(map_stats$nb.markers < 30)])
# Differences between 30 and 40 markers are explained almost exclusively by chromosomes in Camelina_sativa_Singh2015
# However, despite a low number of markers, maps are really good looking
# Consequently, choose to discard chromosome maps below 30 markers
# Discard chromosomes with less than 30 markers
discard_idx = which(map_stats$nb.markers < 30)

#-------------------------------------#
# Filtering by marker density (nb of markers/total genetic length)
# i.e. the averaged number of markers per cM, more is better
summary(map_stats$density.markers.cM)
# At least one marker every 4 centimorgans
# How many chromosomes have less?
tmp_idx = which(map_stats$density.markers.cM < 1/4)
length(tmp_idx)
paste(map_stats$set[tmp_idx], map_stats$chromosome[tmp_idx])

discard_idx = unique(c(discard_idx, tmp_idx))

# And marker density in bp
# At least 1 markers per 10Mb
summary(map_stats$density.markers.bp)*10000000
tmp_idx = which(map_stats$density.markers.bp*10000000 < 1)
length(tmp_idx)
paste(map_stats$set[tmp_idx], map_stats$chromosome[tmp_idx])

discard_idx = unique(c(discard_idx, tmp_idx))

#-------------------------------------#
# Filtering by gaps
# Gaps in genetic distances can indicate errors in genome assembly and will be falsely detected as recombination peaks
# e.g. Draba nivalis chromosome 6
summary(map_stats$largest.gap.cM)
tmp_idx = which(map_stats$largest.gap.cM > 20)
length(tmp_idx)
paste(map_stats$set[tmp_idx], map_stats$chromosome[tmp_idx])
# However, gaps are also dependent on map length (e.g. Triticum aestivum)
# Hence take into account the relative gap size instead of absolute value?
# Gaps larger than 20% of map length
summary(map_stats$largest.gap.cM)
tmp_idx = which(map_stats$largest.gap.cM > 0.2*map_stats$linkage.map.length)
length(tmp_idx)
paste(map_stats$set[tmp_idx], map_stats$chromosome[tmp_idx])

# This criterion discard too much good looking maps
# hence, do not use this criterion to reject maps
# But it can be a real help to identify errors in maps
# discard_idx = unique(discard_idx, tmp_idx)

#-------------------------------------#
# Map coverage (Chakravarti et al. (1991), Hall and Willis (2005))
#         "map coverage c. proportion c of the genome that is within distance d cM of a marker, assuming random distribution of markers,
#         was estimated using c = 1 - exp(-2dn/L), where n is the number of markers and L is the estimated genome length"
#   We estimated the map coverage for a distance d of 1 cM
summary(map_stats$map.coverage)
hist(map_stats$map.coverage)
tmp_idx = which(map_stats$map.coverage < 0.5)
length(tmp_idx)
paste(map_stats$set[tmp_idx], map_stats$chromosome[tmp_idx])

discard_idx = unique(c(discard_idx, tmp_idx))

#-------------------------------------#
# genome coverage
# Maps incomplete
# i.e. starting well after 0
# remove maps with a min physical position > 10% of the total physical length
for (i in 1:nrow(map_stats)) {
  map = read.table(file = paste("data/Marey_maps/",  map_stats$set[i], ".txt", sep = ""), header = TRUE)
  map$phys = as.numeric(map$phys)
  map = subset(map, map == map_stats$chromosome[i])
  map = subset(map, vld == TRUE)
  totallength = max(map$phys, na.rm = TRUE)
  if (nrow(map) > 0) {
    #if (min(map$phys, na.rm = TRUE) > 0.1*totallength | max(map$phys, na.rm = TRUE) < 0.9*totallength) {
    #For now, assuming that we have manually assessed coverage on the max side, because no evidence of the 'true' total length
    if (min(map$phys, na.rm = TRUE) > 0.1*totallength) {
      cat(map_stats$set[i], "chromosome", map_stats$chromosome[i], "\n")
      new_map$vld[new_map$map == map_stats$chromosome[i]] = FALSE
      write.table(new_map, file = paste("data/Marey_maps/", map_stats$set[i], ".txt", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
      set =  map_stats$chromosome[i]
      save.log(msg = paste("Chromosome", map_stats$chromosome[i],"in", map_stats$set[i],"discarded: Map saved.", sep = " "))
    }
  }
}


# Finally remove all chromosome maps in the 'discard_idx'
paste(map_stats$set[discard_idx], map_stats$chromosome[discard_idx])

for (i in discard_idx) {
  # Discard a complete chromosome
  new_map = read.table(file = paste("data/Marey_maps/", map_stats$set[i], ".txt", sep = ""), header = TRUE)
  new_map$vld[new_map$map == map_stats$chromosome[i]] = FALSE
  write.table(new_map, file = paste("data/Marey_maps/", map_stats$set[i], ".txt", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
  set =  map_stats$chromosome[i]
  save.log(msg = paste("Chromosome", map_stats$chromosome[i],"in", map_stats$set[i],"discarded: Map saved.", sep = " "))
}

# Re-iterate auto-update of data-cleaned Marey maps and associated metadata
# AFTER FINAL TRIMMING OF THE DATASET
source("sources/MareyMap.R")
update_datacleaned()



#----------------------------------------#
# Final cleaning step ----
# Final cleaning step - Manually remove maps that are not satysfying
source("sources/MareyMap.R")
data.final = read.table(file = "data-cleaned/marey_maps/AllMaps.txt", header = TRUE, sep = "\t")
map_stats = read.table(file = "tables/chromosome.stats.csv", header = TRUE, sep = ";")

# Initiate a vector of indexes of maps to discard
map_stats$set = as.character(map_stats$set)
map_stats$chromosome = as.character(map_stats$chromosome)
discard_idx = c(which(map_stats$set == "Brassica_napus_Yang2017" & map_stats$chromosome == "A08"),
                which(map_stats$set == "Brassica_napus_Yang2017" & map_stats$chromosome == "C07"),
                which(map_stats$set == "Camellia_sinensis_Xu2018" & map_stats$chromosome == "2"),
                which(map_stats$set == "Camellia_sinensis_Xu2018" & map_stats$chromosome == "4"),
                which(map_stats$set == "Camellia_sinensis_Xu2018" & map_stats$chromosome == "5"),
                which(map_stats$set == "Camellia_sinensis_Xu2018" & map_stats$chromosome == "7"),
                which(map_stats$set == "Citrus_sinensis_Huang2018" & map_stats$chromosome == "2"),
                which(map_stats$set == "Cucurbita_pepo_MonteroPau2017" & map_stats$chromosome == "18"),
                which(map_stats$set == "Cucurbita_pepo_MonteroPau2017" & map_stats$chromosome == "20"),
                which(map_stats$set == "Elaeis_guineensis_Yaakub2020" & map_stats$chromosome == "15"),
                which(map_stats$set == "Triticum_aestivum_GutierrezGonzalez2019" & map_stats$chromosome == "5A"))
discard_idx
paste(map_stats$set[discard_idx], map_stats$chromosome[discard_idx])

for (i in discard_idx) {
  # Discard a complete chromosome
  new_map = read.table(file = paste("data/Marey_maps/", map_stats$set[i], ".txt", sep = ""), header = TRUE)
  new_map$vld[new_map$map == map_stats$chromosome[i]] = FALSE
  write.table(new_map, file = paste("data/Marey_maps/", map_stats$set[i], ".txt", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
  set =  map_stats$chromosome[i]
  save.log(msg = paste("Chromosome", map_stats$chromosome[i],"in", map_stats$set[i],"discarded: Map saved.", sep = " "))
}

# Re-iterate auto-update of data-cleaned Marey maps and associated metadata
# AFTER FINAL TRIMMING OF THE DATASET
source("sources/MareyMap.R")
update_datacleaned()

############################################################################# #
#       POST PROCESSING ----
############################################################################# #

#########-------------------------------########## #
# FROM THIS POINT, READ & SAVE DATA IN 'data-cleaned' ----
#########-------------------------------########## #

#============================================================================#
#============================================================================#
#           DATA IS NOW CLEANED ----
#     From now, we worked on a single data frame containing
#       the concatenated database with all maps cleaned
#============================================================================#
#============================================================================#

# Load the joint database with all maps
data.final = read.table(file = "data-cleaned/marey_maps/AllMaps.txt", header = TRUE, sep = "\t")

#============================================================================#
# View maps ----
#============================================================================#

############# deprecated #####################
# First, use a Unix formatting of text files
#system(paste("dos2unix ", wd, "/data/Marey_maps/*.txt", sep=""))
# Open the GUI with the command:
#startMareyMapGUI()
############# deprecated #####################



# TODO
# Verify congruence chromosome size/LG size -> LG can be named differently from chromosomes
# Manual curating for now


#============================================================================#
# Estimate the local recombination rate ----
#============================================================================#

# Estimate the recombination rate at many positions
# (for instance all the genes of a genome)
# This can be done by up-loading a text file including all the positions.
# a text file (txt extension) containing at least a “map” column and a “phys” column
# indicating respectively the map and the physical position of each gene.

# Individual dataset
# set = "Arabidopsis_thaliana_Serin2017"
# set = "Oryza_sativa_Jiang2017"
# set = "Vitis_vinifera_Brault2020"
# set = "Coffea_canephora_Crouzillat2020"
# set = "Cenchrus_americanus_Pucher2017"
# Otherwise, automatic treatment of all maps
# get the list of maps
list.maps = unique(data.final$set)
source("sources/MareyMap.R")
# Clean directories from previous results
list_dirs = c("output/recombination_maps/loess/1Mbwind", "output/recombination_maps/loess/100kbwind", "output/recombination_maps/loess/bins",
              "output/recombination_maps/loess/pointwise",
              "output/recombination_maps_fig/loess/1Mbwind", "output/recombination_maps_fig/loess/100kbwind", "output/recombination_maps_fig/loess/bins", "output/recombination_maps_fig/loess/pointwise",
              "output/recombination_maps/smooth.spline/1Mbwind", "output/recombination_maps/smooth.spline/100kbwind", "output/recombination_maps/smooth.spline/bins",
              "output/recombination_maps/smooth.spline/pointwise",
              "output/recombination_maps_fig/smooth.spline/1Mbwind", "output/recombination_maps_fig/smooth.spline/100kbwind", "output/recombination_maps_fig/smooth.spline/bins", "output/recombination_maps_fig/smooth.spline/pointwise")
unlink(list_dirs, recursive = TRUE)
sapply(list_dirs, function(x) dir.create(x, recursive = TRUE))


#----------------------------------------------------------------------#
# LOESS method ----
#----------------------------------------------------------------------#
# ?loess # As in MareyMap

# Par ailleurs, dans la mesure où elle repose sur des régressions par les moindres carrés,
# la régression locale bénéficie aussi de la plupart des outils liés à ces méthodes de régression, notamment la théorie de calcul des incertitudes de prédiction et de calibrage.
# Beaucoup d'autres tests et procédures utilisés pour valider les modèles par les moindres carrés peuvent également être étendus aux modèles de régression locale. 
# https://fr.wikipedia.org/wiki/R%C3%A9gression_locale
# https://en.wikipedia.org/wiki/Local_regression
source("sources/MareyMap.R")
(list.maps = unique(data.final$set))

# For a chosen dataset
# set = "Panicum_hallii_Lovell2018"
# set = list.maps[3]
# recombination.map(set = set, chr = "all", method = "loess", K = 3, boot = 1000)

# Or all datasets at once
# get the list of maps
# list.maps = unique(data.final$set)
# list.maps
# list.maps = list.maps[-c(3)]
# for (i in 1:length(list.maps)) {
for (i in c(1:40, 42:52, 54:length(list.maps))) {
  cat("Processing ", as.character(list.maps[i]), "...\n")
  recombination.map(set = as.character(list.maps[i]), chr = "all", method = "loess", K = 3, boot = 1000)
}
# Avoid 41 - Panicum_hallii_Lovell2018 - Too much SNPs - Run on a machine with lots of memory
# Avoid 53 - Triticum_aestivum_GutierrezGonzalez2019 - Too much SNPs - Run on a machine with lots of memory
for (i in c(41, 53)) {
  cat("Processing ", as.character(list.maps[i]), "...\n")
  recombination.map(set = as.character(list.maps[i]), chr = "all", method = "loess", K = 3, boot = 100)
}


# Special cases
# In a few cases, the automatic calibration of loess regression
# produced sub-optimal recombination landscapes
# Hence, loess regression was manually adjusted (span) for this special cases
# Arachis duranensis (4) chr A10, span = 0.4
# Camelina sativa chr () 17
# Cenchrus americanus () chr 1
# Prunus persica (44) chr 5, span = 0.35
# Triticum aestivum (53) chr 1D, span = 0.3
list.maps
source("sources/MareyMap.R")
recombination.map(set = as.character(list.maps[44]), chr = "5",
                  method = "loess", K = 3, boot = 1000, span = 0.35)


#----------------------------------------------------------------------#
# Sliding window ----
#----------------------------------------------------------------------#
# "sliding window, which is the simplest and most widely used method. The idea is just sliding a window along a chromosome and getting the local 
# estimate with the slope of the best line fit to the data in the local window. Parameters are window size and shift." Rezvoy et al. 2007

# May be the worst interpolation method, so no implementation yet

#----------------------------------------------------------------------#
# Cubic spline ----
#----------------------------------------------------------------------#
# "cubic splines, which is probably the best method to estimate recombination rate with Marey map approach (Berloff et al., 2002; Yu et al., 2001)." Rezvoy et al. 2007
# ?smooth.spline # As in MareyMap
# Another implementation is the function qsreg for robust smoothing spline from the package fields (Oh, et al., 2004),
# For a chosen dataset
# recombination.map(set = set, chr = "all", method = "smooth.spline", K = 5, boot = 1000)
source("sources/MareyMap.R")

# Or all datasets at once
for (i in c(1:40, 42:52, 54:length(list.maps))) {
  cat("Processing ", as.character(list.maps[i]), "...\n")
  recombination.map(set = as.character(list.maps[i]), chr = "all", method = "smooth.spline", K = 3, boot = 1000)
}
# Avoid 41 - Panicum_hallii_Lovell2018 - Too much SNPs - Run on a machine with lots of memory
# Avoid 53 - Triticum_aestivum_GutierrezGonzalez2019 - Too much SNPs - Run on a machine with lots of memory
for (i in c(41, 53)) {
  cat("Processing ", as.character(list.maps[i]), "...\n")
  recombination.map(set = as.character(list.maps[i]), chr = "all", method = "smooth.spline", K = 3, boot = 100)
}
# Cubic smooth splines seems to be very sensitive to missing data and produces some artefacts
# Lower confidence in results
# Estimates are smoother
# But may be used as a confirmation of the global pattern for LOESS

# Below are methods of interest to improve esimates, but they have to be tested

#----------------------------------------------------------------------#
# Natural splines method of Berloff (2002) with penalized max likelihood estimation ----
#----------------------------------------------------------------------#
# Custom implementation of the penalized max likelihood estimation with natural splines developped by Berloff et al. (2002)

# TODO Important improvement of smoothing splines

#----------------------------------------------------------------------#
# Thin-plate regression splines ----
#----------------------------------------------------------------------#
# R package mgcv to fit thin-plate regression splines (Wood, 2003).
# TODO


#----------------------------------------------------------------------#
# GAM models with mgcv ----
#----------------------------------------------------------------------#
# TODO
# mod = gam(gen~s(phys), data = chr.data)
# summary(mod)
# gam.check(mod)
# 
# plot(mod,shade=TRUE,seWithMean=TRUE,scale=0)
# 
# predict(mod, newdata=chr.data$phys)

#----------------------------------------------------------------------#
# Bayesian estimates (@petitVariationRecombinationRate2017a) Petit et al. 2017's method ----
#----------------------------------------------------------------------#
# Method based directly on CO positions, not on genetic distances

# TODO


#============================================================================#
# Validate the recombination map ----
#============================================================================#
# A new round of validation for recombination maps
# (1) Maps must not present regions of inconsistent recombination rate (i.e. over 30 cM.Mb-1)
# (2) CI must be under a threshold of quality
# (3) loess and smooth.spline must show convergence
# (4) Maps must not show a high discrepancy between the chromosome-wide recombiantion rate and the estimated mean recombination rate

# Otherwise, maps are discarded in the original data and outputs/figures are removed
#discard.map() # Remove data in the original dataset

# Then analyses are runned again to update the results


#----------------------------------------------------------------------#
# RMSE method
#----------------------------------------------------------------------#



#============================================================================#
# VALIDATION OF ESTIMATED RECOMBINATION MAPS ----
#============================================================================#
# Recombination map selection step
# A function returning a data frame with summary statistics of selected maps
# AND saving a file with summary statistics if save.file = TRUE
source("sources/MareyMap.R")
# Load the joint database with all maps
data.final = read.table(file = "data-cleaned/marey_maps/AllMaps.txt", header = TRUE, sep = "\t")
chromosome.stats = recombination.map.statistics(data.final, save.file = FALSE)

#----------------------------------------------------------------------#
# Qualitative assessment - Visual validation

#----------------------------------------------------------------------#
# Quantitative assessment - Correlation between the chromosome wide recombination rate & the mean estimated recombination rate
# Show the linear relationship between chromosome wide rate and the mean estimated recombination rate as a cross-validation procedure
plot(chromosome.stats$chrwide.rate, chromosome.stats$mean.recrate)
abline(a = 0, b = 1, col = "Red")
# Same mean: no bias in interpolation, but may be noise produced (variance)
# Rejection of maps too far from the QQline
# i.e. points that are out of the confidence interval around the QQline (arbitrary interval of +- 0.2)

# Plotting the interval
# quantile_interval = data.frame(x = round(seq(0, 10, 0.0001), digits = 4), y.upper = seq(0, 10, 0.0001) + interval, y.lower = seq(0, 10, 0.0001) - interval)
# lines(quantile_interval$x, quantile_interval$y.upper, col = "Red")
# lines(quantile_interval$x, quantile_interval$y.lower, col = "Red")

# Interval calibration by exploring distributions
plot(density(abs(chromosome.stats$mean.recrate-chromosome.stats$chrwide.rate), na.rm = TRUE), col = "Red")
abline(v = 1, col = "Red", lty = 2)
# We see a peak around 0.1-0.2

# test if a map is outside the interval
# See interval calibration above
interval = 1
df_test = data.frame(x = round(chromosome.stats$chrwide.rate , digits = 4), y = chromosome.stats$mean.recrate) # A data frame of map coordinates to test
idx_reject = logical(nrow(df_test)) # the index of maps rejected
for (i in 1:nrow(df_test)) {
  idx_reject[i] = (df_test$y[i] > quantile_interval$y.upper[which(quantile_interval$x == df_test$x[i])]) | (df_test$y[i] < quantile_interval$y.lower[quantile_interval$x == df_test$x[i]])
}
sum(idx_reject)
# plot rejected maps in red
plot(chromosome.stats$chrwide.rate, chromosome.stats$mean.recrate)
abline(a = 0, b = 1, col = "Red")
points(chromosome.stats$chrwide.rate[idx_reject], chromosome.stats$mean.recrate[idx_reject], col = "Red", bg = "Red", pch = 16)

# Which maps are outliers
paste(chromosome.stats$set[which(idx_reject)], "chr", chromosome.stats$chromosome[which(idx_reject)], sep = " ")

# Then discard low quality maps
for (i in which(idx_reject)) {
  # Discard a complete chromosome
  new_map = read.table(file = paste("data/Marey_maps/", map_stats$set[i], ".txt", sep = ""), header = TRUE)
  new_map$vld[new_map$map == map_stats$chromosome[i]] = FALSE
  write.table(new_map, file = paste("data/Marey_maps/", map_stats$set[i], ".txt", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
  set =  map_stats$chromosome[i]
  save.log(msg = paste("Chromosome", map_stats$chromosome[i],"in", map_stats$set[i],"discarded: Map saved.", sep = " "))
}

#----------------------------------------------------------------------#
# Quantitative assessment - Confidence Intervals



#----------------------------------------------------------------------#
# Now re-iterate update_cleaned data and map recombination rate estimation...
#----------------------------------------------------------------------#






#============================================================================#
# Significance testing of the recombination map ----
#============================================================================#
# Which part of the recombination landscape have a recombination rate significantly higher than the average recombination rate?





#============================================================================#
# Save the recombination map ----
#============================================================================#

#----------------------------------------------------------------------#
# Saving in the directory '/output/MareyMaps'
# Maps can be saved to R data files (rda, Rda, rdata or Rdata) or to text files (txt).


#----------------------------------------------------------------------#
# Save in a txt file the estimated local recombination rate for each marker position



#----------------------------------------------------------------------#
# Indicate the centromere position
# BLAST the centromeric motif, since centromeres are not necessarily at the center of the chromosome (e.g. metacentric, acrocentric or telocentric)

#============================================================================#
# END ----
#============================================================================#