########################################################################## #
#                    ECOBIO - PhD
#
#                 Genome GC
#         Exploration of genome characteristics and recombination
#             as a function of GC
#
########################################################################## #
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
# library(MareyMap)
# library(stringr) 
require(gdata)
require(taxize)
require(taxizedb)
require(dplyr)
require(rentrez)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# install.packages("devtools")
# devtools::install_github("mhahsler/rBLAST")
library(rBLAST) #
# BiocManager::install("Biostrings")
library(Biostrings) #
#https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html
require(parallel)
require(MASS)
require(pbmcapply)


#============================================================================#
# Loading variables & objects ----
#============================================================================#
# Get the directory of the file & set working directory
wd=dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

# Custom functions
# source("sources/PhysMap.R") # Get the physical positions of a marker given a list of marker names
# source("sources/MareyMap.R") # Consolidate a MareyMap input file
source("sources/core_pipeline.R") # Basic functions but essential for the pipeline (e.g. saving logs)
# source("sources/Genome.R") # Essential functions to parse GFF files and compute genomic statistics

chromosome.stats = read.table(file = "tables/chromosome.stats.csv", header = TRUE, sep = ";")
metadata.clean = read.table("data-cleaned/Genetic_maps_ressources.csv", header = TRUE, sep = ";")

#============================================================================#
# List of genomes in the dataset ----
#============================================================================#
######## Work on genome in the accession list: i.e. reference genomes for validated recombination maps
# Load genome database
genomes = read.table("data/Genome/Genome_ressources.csv", header = TRUE, sep = ";")
# Import list of species
list_species = unique(genomes$species)
# Import list of accessions to process
list_accessions = unique(genomes$accession)

# The list of actual species in the final dataset
unique(chromosome.stats$species)
# The list of actual dataset
unique(chromosome.stats$set)



#============================================================================#
# ESTIMATE GC CONTENT AT A GENE LEVEL ----
#============================================================================#
# GFF-based files
# with supplementary columns
# col. 10 = exon rank (if pertinent)
# col. 11 = gene id
# col. 12 = recombination rate
# col. 13-16= GC/GC1/GC2/GC3

# Produce one file for each dataset, saved in the "gc_genes" directory
unique(chromosome.stats$set)
# TEST
source("sources/GC.R") # Essential functions to parse GFF files and compute genomic statistics
# Produce GC maps at a genomic feature level ----
results = GC_content(map, parallel = TRUE)
summary(results)

 
# Choose a map - GC maps estimated...
# One by one procedure - Data heterogeneity
  map = "Arabidopsis_thaliana_Serin2017"
map = "Arachis_duranensis_Bertioli2016"
map = "Arachis_hypogaea_Zhuang2019"
# map = "Boechera_stricta_Lee2020" # No genome
  map = "Brachypodium_distachyon_Huo2011"
map = "Camelina_sativa_King2019"
map = "Camellia_sinensis_Xu2018" # NO EXONS; Errors with IDs
map = "Brassica_rapa_CodyMarkelz2017" # Problem with exon rank; no info
map = "Capsella_rubella_Slotte2013"
map = "Cucumis_sativus_Zhu2016"
map = "Cucurbita_pepo_MonteroPau2017"
map = "Dioscorea_alata_Cormier2019"
map = "Elaeis_guineensis_Yaakub2020"
map = "Eucalyptus_grandis_Bertholome2015"
map = "Glycine_max_Watanabe2017"
map = "Brassica_napus_Yang2017"
map = "Lupinus_albus_Ksiazkiewicz2017" # No exons
map = "Lupinus_angustifolius_Zhou2017"
map = "Malus_domestica_DiPierro2016"
map = "Manihot_esculenta_ICGMC2015"
map = "Oryza_sativa_DeLeon2016"
map = "Panicum_hallii_Lovell2018"
  map = "Phaseolus_vulgaris_Song2015" # Recombination rates did not worked
map = "Prunus_mume_Zhang2015"
map = "Prunus_persica_Verde2017"
map = "Sesamum_indicum_Wang2016"
map = "Setaria_italica_Bennetzen2012" # NO EXONS
map = "Solanum_lycopersicum_Gonda2018"
  map = "Sorghum_bicolor_Zou2012"
map = "Theobroma_cacao_Royaert2016" # No exons
map = "Vigna_unguiculata_Lonardi2019"

map = "Zea_mays_MaizeGDBConsensus_v4"





 
map = "Hordeum_vulgare_MunozAmatriain2014" # Memory overflow; use parallel = FALSE




map = "Capsicum_annuum_Han2016" # Not working !!! No exon positions; nothing to do
  map = "Cenchrus_americanus_Pucher2017" # Not working : no GFF -> convert gbff to gff Error importing gff: Not a regular file
map = "Citrullus_lanatus_Ren2015" # Not working !!! No exon positions; nothing to do
  map = "Citrus_sinensis_Huang2018" # Not working !!! -> only scaffolds in gff -> search better annotation file : use GCA_000317415.1 Error importing gff: Not a regular file
map = "Coffea_canephora_Crouzillat2020" # Not working : no GFF; no solution
map = "Cucumis_melo_Pereira2018" # Not working !!!
map = "Cucurbita_maxima_Wang2020"
  map = "Draba_nivalis_Birkeland2020" # Error importing gff: Not a regular file

  map = "Gossypium_hirsutum_Zhang2019" # Erreur : $ operator is invalid for atomic vectors
  map = "Gossypium_raimondii_Wang2013" # Erreur : $ operator is invalid for atomic vectors
  map = "Helianthus_annuus_Talukder2014" # Erreur : $ operator is invalid for atomic vectors
  map = "Oryza_nivara_Ma2016"
  map = "Solanum_tuberosum_Endelman2016"
  map = "Triticum_aestivum_GutierrezGonzalez2019"
  map = "Vitis_vinifera_Brault2020"
  
  
  map = "Mangifera_indica_Luo2016" # GBFF to convert...
  
map = "Juglans_regia_Luo2015" # No genome or annotation found !
map = "Nelumbo_nucifera_Gui2018" # No GFF
map = "Quercus_sp_Bodenes2016" # le tableau de remplacement a 1 lignes, le tableau remplacé en a 0
map = "Raphanus_sativus_Mun2015" # Pas de gff
map = "Triticum_dicoccoides_Jorgensen2017" # Annotation does not contain genes, CDS or exons; problem with GFF conversion?
map = "Triticum_urartu_Ling2018" # Annotation does not contain genes, CDS or exons; problem with GFF conversion?


