############################################################################
#                    ECOBIO - PhD
#
#       Simulate Recombination maps (genetic vs physical maps)
#     To test sensitivyt of different methods of estimating recombination rates
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
# R/Qtl: a package for simulating a genetic map
# library(qtl)

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
# SIMULATIONS
#============================================================================

# Marey map can be modelled by the logit function, the reciprocal function of the sigmoid, or a beta CDF
# with different degrees of skewness toward the telomere/ flatness of the center/centromeres
# And then, noise and missing data can be applied on the theoretical ideal Marey map (i.e. a monotone increasing function)


#============================================================================
# SENSITIVITY
#============================================================================




#============================================================================
# Loading variables & objects
#============================================================================
