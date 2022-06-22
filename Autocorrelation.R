############################################################################
#                    ECOBIO - PhD
#
#       Broken stick model
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

#============================================================================
# Plot for autocorrelation in the recombination landscapes
#============================================================================
# Test a single map first
set = "Arabidopsis_thaliana_Serin2017"
set = "Zea_mays_IBM_MaizeSNP50"
chromosome = "1"
map = read.table(file = paste("output/recombination_maps/loess/100kbwind/", set, "_chromosome", chromosome,".txt", sep = ""), header = TRUE, sep = "\t")
# map = read.table(file = paste("output/recombination_maps/loess/1Mbwind/", set, "_chromosome", chromosome,".txt", sep = ""), header = TRUE, sep = "\t")

# use of the acf() function of R
# Lag is the number of windows in 100kb
# e.g. lag = 5 means 500kb (5 windows)
# Lag max adjusted to the length of the sequence (i.e. number of windows)
autocorr = acf(map$rec.rate[!is.na(map$rec.rate)], lag.max = length(map$rec.rate))

# Compare two species and two scales
par(mfrow = c(2,2))
set = "Arabidopsis_thaliana_Serin2017"
map = read.table(file = paste("output/recombination_maps/loess/100kbwind/", set, "_chromosome", chromosome,".txt", sep = ""), header = TRUE, sep = "\t")
autocorr = acf(map$rec.rate[!is.na(map$rec.rate)], lag.max = length(map$rec.rate),
               main = "A. thaliana chromosome 1, windows of 100kb")
set = "Zea_mays_IBM_MaizeSNP50"
map = read.table(file = paste("output/recombination_maps/loess/100kbwind/", set, "_chromosome", chromosome,".txt", sep = ""), header = TRUE, sep = "\t")
autocorr = acf(map$rec.rate[!is.na(map$rec.rate)], lag.max = length(map$rec.rate),
               main = "Z. mays chromosome 1, windows of 100kb")
set = "Arabidopsis_thaliana_Serin2017"
map = read.table(file = paste("output/recombination_maps/loess/1Mbwind/", set, "_chromosome", chromosome,".txt", sep = ""), header = TRUE, sep = "\t")
autocorr = acf(map$rec.rate[!is.na(map$rec.rate)], lag.max = length(map$rec.rate),
               main = "A. thaliana chromosome 1, windows of 1Mb")
set = "Zea_mays_IBM_MaizeSNP50"
map = read.table(file = paste("output/recombination_maps/loess/1Mbwind/", set, "_chromosome", chromosome,".txt", sep = ""), header = TRUE, sep = "\t")
autocorr = acf(map$rec.rate[!is.na(map$rec.rate)], lag.max = length(map$rec.rate),
               main = "Z. mays chromosome 1, windows of 1Mb")
par(mfrow = c(1,1))


# Is the autocorrelation lag similar between species
# Check autocorrelation in bins
# Compare two species and two scales
par(mfrow = c(1,2))
set = "Arabidopsis_thaliana_Serin2017"
map = read.table(file = paste("output/recombination_maps/loess/bins/", set, "_chromosome", chromosome,".txt", sep = ""), header = TRUE, sep = "\t")
autocorr = acf(map$rec.rate[!is.na(map$rec.rate)], lag.max = length(map$rec.rate),
               main = "A. thaliana chromosome 1, 1000 bins")
set = "Zea_mays_IBM_MaizeSNP50"
map = read.table(file = paste("output/recombination_maps/loess/bins/", set, "_chromosome", chromosome,".txt", sep = ""), header = TRUE, sep = "\t")
autocorr = acf(map$rec.rate[!is.na(map$rec.rate)], lag.max = length(map$rec.rate),
               main = "Z. mays chromosome 1, 1000 bins")
par(mfrow = c(1,1))

#============================================================================
# Test for autocorrelation in the recombination landscapes
#============================================================================

set = "Arabidopsis_thaliana_Serin2017"
set = "Zea_mays_IBM_MaizeSNP50"
chromosome = "1"
map = read.table(file = paste("output/recombination_maps/loess/100kbwind/", set, "_chromosome", chromosome,".txt", sep = ""), header = TRUE, sep = "\t")

map = map[!is.na(map$phys),]
map = map[!is.na(map$rec.rate),]



#============================================================================
# END 
#============================================================================