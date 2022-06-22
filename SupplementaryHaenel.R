########################################################################### #
#                    ECOBIO - PhD
#
#       Meta-analysis Haenel et al. 2019
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

#----------------------#
# Loading packages
library(rstudioapi)
library(ggplot2)
library(MareyMap)
library(stringr)
library(lme4)
library(caper)
library(car)
library(lmerTest) # add p-values to lme4 models

#============================================================================#
# Loading variables & objects ----
#============================================================================#

# Get the directory of the file & set working directory
wd=dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

# Custom functions
source("sources/PhysMap.R") # Get the physical positions of a marker given a list of marker names
source("sources/MareyMap.R")
source("sources/core_pipeline.R") # Basic functions but essential for the pipeline (e.g. saving logs)
source("sources/stats.R") # Statistical functions implemented

#============================================================================#
# Summary statistics ----
#============================================================================#
haenel = read.table("data/Supplementary_Haenel2019.csv", header = TRUE, sep = ",")

summary(haenel)

summary(haenel$Physical.map.length.Mb[which(haenel$Kingdom == "Plants")])
summary(haenel$Physical.map.length.Mb[which(haenel$Kingdom == "Fungi")])
summary(haenel$Physical.map.length.Mb[which(haenel$Kingdom == "Animals")])

boxplot(haenel$Physical.map.length.Mb ~ haenel$Kingdom)

summary(haenel$Genetic.map.length.cM[which(haenel$Kingdom == "Plants")])
# Ratio of change
max(haenel$Genetic.map.length.cM[which(haenel$Kingdom == "Plants")])/min(haenel$Genetic.map.length.cM[which(haenel$Kingdom == "Plants")])
# Range
max(haenel$Genetic.map.length.cM[which(haenel$Kingdom == "Plants")]) - min(haenel$Genetic.map.length.cM[which(haenel$Kingdom == "Plants")])

summary(haenel$Genetic.map.length.cM[which(haenel$Kingdom == "Fungi")])
# Ratio of change
max(haenel$Genetic.map.length.cM[which(haenel$Kingdom == "Fungi")])/min(haenel$Genetic.map.length.cM[which(haenel$Kingdom == "Fungi")])
max(haenel$Genetic.map.length.cM[which(haenel$Kingdom == "Fungi")]) - min(haenel$Genetic.map.length.cM[which(haenel$Kingdom == "Fungi")])


summary(haenel$Genetic.map.length.cM[which(haenel$Kingdom == "Animals")])
# Subset vertebrates
vertebrates = unique(haenel$Species)[c(2,6,8:11,14:15,17,27,29,32,35,41,48)]
summary(haenel$Genetic.map.length.cM[which(haenel$Species %in% 
                                             vertebrates)])
summary(haenel$Physical.map.length.Mb[which(haenel$Species %in% 
                                             vertebrates)])
# Subset mammals
mammals = unique(haenel$Species)[c(9:11,29,41)]
summary(haenel$Genetic.map.length.cM[which(haenel$Species %in% 
                                             mammals)])
summary(haenel$Physical.map.length.Mb[which(haenel$Species %in% 
                                              mammals)])

# Ratio of change
max(haenel$Genetic.map.length.cM[which(haenel$Kingdom == "Animals")])/min(haenel$Genetic.map.length.cM[which(haenel$Kingdom == "Animals")])
# Range
max(haenel$Genetic.map.length.cM[which(haenel$Kingdom == "Animals")]) - min(haenel$Genetic.map.length.cM[which(haenel$Kingdom == "Animals")])

boxplot(haenel$Genetic.map.length.cM ~ haenel$Kingdom)
ggplot(data = haenel, aes(x = Kingdom, y = Genetic.map.length.cM)) +
  geom_violin()

# ratio of change  and range in Stapley et al. 2017
# Plants
8184/309
8184 - 309
# Animals
5961/90
5961 - 90
# Fungi
5860/86


summary(haenel$Crossover.rate[which(haenel$Kingdom == "Plants")])
# Ratio of change
max(haenel$Crossover.rate[which(haenel$Kingdom == "Plants")])/min(haenel$Crossover.rate[which(haenel$Kingdom == "Plants")])
summary(haenel$Crossover.rate[which(haenel$Kingdom == "Fungi")])
# Ratio of change
max(haenel$Crossover.rate[which(haenel$Kingdom == "Fungi")])/min(haenel$Crossover.rate[which(haenel$Kingdom == "Fungi")])
summary(haenel$Crossover.rate[which(haenel$Kingdom == "Animals")])
# Ratio of change
max(haenel$Crossover.rate[which(haenel$Kingdom == "Animals")])/min(haenel$Crossover.rate[which(haenel$Kingdom == "Animals")])

boxplot(haenel$Crossover.rate ~ haenel$Kingdom)
ggplot(data = haenel, aes(x = Kingdom, y = Crossover.rate)) +
  geom_violin()


# Linkage map length of the largest plant genome
plant_subset = haenel[which(haenel$Kingdom == "Plants"),]
plant_subset$Genetic.map.length.cM[which(plant_subset$Species == plant_subset$Species[which.max(plant_subset$Physical.map.length.Mb)])]
summary(plant_subset$Genetic.map.length.cM[which(plant_subset$Species == plant_subset$Species[which.max(plant_subset$Physical.map.length.Mb)])])
sum(plant_subset$Genetic.map.length.cM[which(plant_subset$Species == plant_subset$Species[which.max(plant_subset$Physical.map.length.Mb)])])

# Linkage map length of the largest animal genome
animal_subset = haenel[which(haenel$Kingdom == "Animals"),]
animal_subset$Genetic.map.length.cM[which(animal_subset$Species == animal_subset$Species[which.max(animal_subset$Physical.map.length.Mb)])]
summary(animal_subset$Genetic.map.length.cM[which(animal_subset$Species == animal_subset$Species[which.max(animal_subset$Physical.map.length.Mb)])])
sum(animal_subset$Genetic.map.length.cM[which(animal_subset$Species == animal_subset$Species[which.max(animal_subset$Physical.map.length.Mb)])])
