############################################################################
#                    ECOBIO - PhD
#
#           Exploring patterns
#       Correlation tests and Models
#
############################################################################
# by Thomas Brazier
# brazier.thomas@gmail.com

# Supervisor: Sylvain GLEMIN
#             ECOBIO Lab


#============================================================================
# LOADING ENVIRONMENT
#============================================================================

# clear global environment: remove all variables
rm(list=ls(all=TRUE))

#----------------------
# Loading packages
library(rstudioapi)
library(ggplot2)
library(MareyMap)
library(stringr) 

#============================================================================
# Loading variables & objects
#============================================================================

# Get the directory of the file & set working directory
filedir=dirname(rstudioapi::getSourceEditorContext()$path)
wd=filedir
setwd(wd)

# Custom functions
source("sources/PhysMap.R") # Get the physical positions of a marker given a list of marker names
source("sources/MareyMap.R") # Consolidate a MareyMap input file
source("sources/core_pipeline.R") # Basic functions but essential for the pipeline (e.g. saving logs)
source("sources/stats.R") # Statistical functions implemented

metadata = read.csv(file = paste(wd, "/data/Genetic_maps/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")
# Load the joint database with all maps
data.final = read.table(file = "data-cleaned/marey_maps/AllMaps.txt", header = TRUE, sep = "\t")
chromosome.stats = read.table("tables/chromosome.stats.csv", header = TRUE)

#============================================================================
# Figure. Mean recombination rate as a function of physical length for each chromosome
#============================================================================

# Log ~ log shows a linear relationship
plot(log(chromosome.stats$phys.map.length), log(chromosome.stats$mean.recrate))
mod = lm(log(mean.recrate) ~ log(phys.map.length), data = chromosome.stats)
abline(mod)
summary(mod)
cor.test(log(chromosome.stats$phys.map.length), log(chromosome.stats$mean.recrate))



#============================================================================
# END
#============================================================================