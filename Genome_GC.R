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








#============================================================================#
# ESTIMATE GC CONTENT AT A GENOME LEVEL ----
#============================================================================#
# Recombination map like file
# with supplementary columns for 100kb windows
# GC/GC1/GC2/GC3






#============================================================================#
# DEPRECATED UNDER THIS LINE
#============================================================================#


# #============================================================================#
# # Recombination map level
# #============================================================================#
# (species = as.character(metadata.clean$species[which(metadata.clean$id == map)]))
# (accession = as.character(metadata.clean$accession[which(metadata.clean$id == map)]))
# source("sources/Genome.R") # Essential functions to parse GFF files and compute genomic statistics
# exonpos = exon.position(species, accession)
# write.table(exonpos, gzfile(paste("data-cleaned/genome/exon_positions/", species,"_exon_positions.txt.gz", sep = "")), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)
# 
# # Process a dataset
# # source("sources/Genome.R") # Essential functions to parse GFF files and compute genomic statistics
# gc.map(map)
# 
# # RECOMBINATION MAPS
# gcmap = read.table(file = gzfile(paste("data-cleaned/genome/gc_maps/", map, "_chromosome1.txt.gz", sep = "")), header = TRUE, sep = "\t")
# # Filter recombination rates
# # values at 0 are considered as errors in interpolation
# # Hence these windows are masked for analyses, by converting them to NAs
# gcmap$rec.rate[gcmap$rec.rate < 0] = NA
# 
# #============================================================================#
# # Gene level
# #============================================================================#
# 
# # GENES
# gcgenes = read.table(file = gzfile(paste("data-cleaned/genome/gc_genes/", map, ".txt.gz", sep = "")), header = TRUE, sep = "\t")
# 
# #----------------------------------------------------------------------------#
# # Gene size in bp
# #----------------------------------------------------------------------------#
# species = "Malus_domestica"
# genepos = read.table(file = gzfile(paste("data-cleaned/genome/gene_positions/", species, "_gene_positions.txt.gz", sep = "")), header = TRUE, sep = "\t")
# 
# gene_size = genepos$end - genepos$start
# 
# summary(gene_size)
# hist(gene_size, freq = FALSE, breaks = 100)
# lines(density(gene_size), col = "red")
# 
# #----------------------------------------------------------------------------#
# # Gene size - Number of exons
# #----------------------------------------------------------------------------#
# # Initial design for
# map = "Malus_domestica_DiPierro2016"
# 
# # Choose a map
# map = "Arabidopsis_thaliana_Serin2017" #DONE
# map = "Oryza_sativa_DeLeon2016" #DONE
# map = "Malus_domestica_DiPierro2016" #DONE
# map = "Brachypodium_distachyon_Huo2011" # DONE
# # map = "Oryza_sativa_Jiang2017" 
# # map = "Capsicum_annuum_Han2016" 
# map = "Cucumis_sativus_Zhu2016" # DONE
# map = "Gossypium_raimondii_Wang2013" # DONE
# map = "Phaseolus_vulgaris_Song2015" # DONE
# map = "Sorghum_bicolor_Zou2012" # DONE
# # map = "Theobroma_cacao_Royaert2016"
# map = "Zea_mays_MaizeGDBConsensus_v4" # DONE
# map = "Zea_mays_IBM_MaizeSNP50" # DONE
# map = "Setaria_italica_Ni2017" # DONE
# map = "Solanum_lycopersicum_Gonda2018" # DONE
# map = "Triticum_aestivum_GutierrezGonzalez2019"
# 
# # Produce data from genome annotations
# 
# # A function that take map name and return & save a file containing the list of genes and the number of exons associated
# # Plus additional informations such as recombination and GC3 (gene) + GC3 exon 1
# source("sources/Genome.R")
# df = genesize_nbexons(map)
# # Display results
# df = read.table(file = gzfile(paste("data-cleaned/genome/number_exons/", map, "_exonnumber.txt.gz", sep = "")), header = TRUE, sep = "\t")
# summary(df$number_exons)
# hist(df$number_exons, freq = FALSE, breaks = max(df$number_exons))
# lines(density(df$number_exons), col = "red")
# # Rec rate as a function of number of exons
# nbexons = unique(df$number_exons)
# recrate = numeric(length(nbexons))
# for (i in 1:length(nbexons)) {
#   recrate[i] = mean(df$rec.rate[which(df$number_exons == nbexons[i])], na.rm = TRUE)
# }
# plot(nbexons, recrate, xlim = c(0, 20))
# # GC3 as a function of number of exons
# nbexons = unique(df$number_exons)
# gc3_for_exons = numeric(length(nbexons))
# for (i in 1:length(nbexons)) {
#   gc3_for_exons[i] = mean(df$GC3[which(df$number_exons == nbexons[i])], na.rm = TRUE)
# }
# plot(nbexons, gc3_for_exons, xlim = c(0, 20))
# # GC3ex1 as a function of number of exons
# nbexons = unique(df$number_exons)
# gc3ex1_for_exons = numeric(length(nbexons))
# for (i in 1:length(nbexons)) {
#   gc3ex1_for_exons[i] = mean(df$GC3ex1[which(df$number_exons == nbexons[i])], na.rm = TRUE)
# }
# plot(nbexons, gc3ex1_for_exons, xlim = c(0, 20))
# 
# # GC3 gradients
# # Load data
# gcexons = read.table(file = gzfile(paste("data-cleaned/genome/gc_exons/", map, ".txt.gz", sep = "")), header = TRUE, sep = "\t")
# # Init a data frame with rank and the mean GC3 for the rank
# gradient = data.frame(rank = sort(unique(gcexons$rank)), GC3 = NA)
# for (i in 1:nrow(gradient)) {
#   gradient$GC3[i] = mean(gcexons$gc3[gcexons$rank == gradient$rank[i]], na.rm = TRUE)
# }
# # Figure
# plot(gradient$rank, gradient$GC3, xlim = c(0, 20))
# 
# #----------------------------------------------------------------------------#
# # FIGURES
# #----------------------------------------------------------------------------#
# 
# # Gather data for figures
# # Combine all species
# list_maps = c("Arabidopsis_thaliana_Serin2017",  "Oryza_sativa_DeLeon2016",
#  "Malus_domestica_DiPierro2016" ,  "Brachypodium_distachyon_Huo2011", "Cucumis_sativus_Zhu2016", "Gossypium_raimondii_Wang2013",
#  "Phaseolus_vulgaris_Song2015", "Sorghum_bicolor_Zou2012", "Zea_mays_MaizeGDBConsensus_v4", "Zea_mays_IBM_MaizeSNP50",
#  "Setaria_italica_Ni2017", "Solanum_lycopersicum_Gonda2018", "Triticum_aestivum_GutierrezGonzalez2019")
# 
# df_all = read.table(file = gzfile(paste("data-cleaned/genome/number_exons/", list_maps[1], "_exonnumber.txt.gz", sep = "")), header = TRUE, sep = "\t")
# df_all$species = gsub("_", " ",gsub("[A-Za-z0-9]*$", "",list_maps[1]))
# for (i in 2:length(list_maps)) {
#   df = read.table(file = gzfile(paste("data-cleaned/genome/number_exons/", list_maps[i], "_exonnumber.txt.gz", sep = "")), header = TRUE, sep = "\t")
#   df$species = gsub("_", " ",gsub("[A-Za-z0-9]*$", "",list_maps[i]))
#   df_all = rbind(df_all, df)
#   rm(df)
# }
# 
# #----------------------------------------------------------------------------#
# # Gene size distribution (number of exons)
# # Plot the density of number of exons per species
# Dist_exonsnumber = ggplot(data = df_all, aes(x = number_exons, color = species)) +
#   geom_density(aes(), bw = 1) +
#   labs(x = "Number of exons", y="Density", fill="Species") +
#   xlim(0, 20) +
#   theme(axis.line = element_line(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
#         plot.subtitle = element_text(color="black",size=28,hjust = 0.5),
#         axis.title.x = element_text(color="black", size=28),
#         axis.title.y = element_text(color="black", size=28),
#         axis.text=element_text(size=28, colour="black"),
#         strip.text=element_text(size=28, colour="black"),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.text=element_text(size=28),
#         legend.title=element_text(size=28)
#   )
# Dist_exonsnumber
# ggsave("figures/Dist_exonsnumber.png",
#        device="png",dpi=320,units="cm",width=70,height=40)
# 
# #----------------------------------------------------------------------------#
# # Shorter genes recombine more ?
# # Rec rate as a function of number of exons
# # Make a data frame with mean of recombination per number of exons and species
# # Only for 1 to 20 exons
# data_plot = data.frame(species = rep(unique(df_all$species), each = 20), number_exons = rep(1:20, times = length(unique(df_all$species))), rec.rate = NA, rec.rate.upper = NA, rec.rate.lower = NA)
# # Compute mean recombination and CI by bootstrap
# for (i in 1:nrow(data_plot)) {
#   cat("Row", i, "of", nrow(data_plot), "\n")
#   nboot = 1000
#   boot = numeric(nboot)
#   for (j in 1:nboot) {
#     boot[j] = mean(sample(df_all$rec.rate[df_all$species == data_plot$species[i] & df_all$number_exons == data_plot$number_exons[i]], replace = TRUE), na.rm = TRUE)
#   }
#   data_plot$rec.rate[i] = mean(boot)
#   data_plot$rec.rate.upper[i] = quantile(boot, 0.975)
#   data_plot$rec.rate.lower[i] = quantile(boot, 0.025)
# }
# # Save in a file, to avoid recomputing (bootstrapping is long)
# write.table(data_plot, file = paste("output/figures/Recombination~exonnumber.txt", sep = ""),
#             quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
# # Import file
# data_plot = read.table(file = paste("output/figures/Recombination~exonnumber.txt", sep = ""), header = TRUE, sep = "\t")
# # Add a small horizontal jitter
# data_plot$number_exons = data_plot$number_exons + rep(seq(0, 0.5, length.out = length(unique(data_plot$species))), each = 20)
# 
# # Plot
# Rec_exonsnumber = ggplot(data = data_plot, aes(x = number_exons, y = rec.rate, colour = species)) +
#   # One point per species and number of exon
#   geom_point() +
#   # Add confidence intervals
#   geom_errorbar(aes(ymin = rec.rate.lower, ymax = rec.rate.upper), width = 0) +
#   labs(x = "Number of exons", y="Recombination rate (cM/Mb)", fill="Species") +
#   theme(axis.line = element_line(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
#         plot.subtitle = element_text(color="black",size=28,hjust = 0.5),
#         axis.title.x = element_text(color="black", size=28),
#         axis.title.y = element_text(color="black", size=28),
#         axis.text=element_text(size=28, colour="black"),
#         strip.text=element_text(size=28, colour="black"),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.text=element_text(size=28),
#         legend.title=element_text(size=28)
#   )
# Rec_exonsnumber
# ggsave("figures/Recombination~exonnumber.png",
#        device="png",dpi=320,units="cm",width=70,height=40)
# 
# #----------------------------------------------------------------------------#
# # Shorter genes are supposed to be more GC-rich
# # GC3 as a function of number of exons
# # Make a data frame with mean of GC3 per number of exons and species
# # Only for 1 to 20 exons
# data_plot = data.frame(species = rep(unique(df_all$species), each = 20), number_exons = rep(1:20, times = length(unique(df_all$species))), GC3 = NA, GC3.upper = NA, GC3.lower = NA)
# # Compute mean recombination and CI by bootstrap
# for (i in 1:nrow(data_plot)) {
#   cat("Row", i, "of", nrow(data_plot), "\n")
#   nboot = 1000
#   boot = numeric(nboot)
#   for (j in 1:nboot) {
#     boot[j] = mean(sample(df_all$GC3[df_all$species == data_plot$species[i] & df_all$number_exons == data_plot$number_exons[i]], replace = TRUE), na.rm = TRUE)
#   }
#   data_plot$GC3[i] = mean(boot)
#   data_plot$GC3.upper[i] = quantile(boot, 0.975)
#   data_plot$GC3.lower[i] = quantile(boot, 0.025)
# }
# # Save in a file, to avoid recomputing (bootstrapping is long)
# write.table(data_plot, file = paste("output/figures/GC3~exonnumber.txt", sep = ""),
#             quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
# # Import file
# data_plot = read.table(file = paste("output/figures/GC3~exonnumber.txt", sep = ""), header = TRUE, sep = "\t")
# # Add a small horizontal jitter
# data_plot$number_exons = data_plot$number_exons + rep(seq(0, 0.5, length.out = length(unique(data_plot$species))), each = 20)
# 
# # Plot
# GC3_exonsnumber = ggplot(data = data_plot, aes(x = number_exons, y = GC3, colour = species)) +
#   # One point per species and number of exon
#   geom_point() +
#   # Add confidence intervals
#   geom_errorbar(aes(ymin = GC3.lower, ymax = GC3.upper), width = 0) +
#   labs(x = "Number of exons", y="GC3 proportion", fill="Species") +
#   theme(axis.line = element_line(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
#         plot.subtitle = element_text(color="black",size=28,hjust = 0.5),
#         axis.title.x = element_text(color="black", size=28),
#         axis.title.y = element_text(color="black", size=28),
#         axis.text=element_text(size=28, colour="black"),
#         strip.text=element_text(size=28, colour="black"),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.text=element_text(size=28),
#         legend.title=element_text(size=28)
#   )
# GC3_exonsnumber
# ggsave("figures/GC3~exonnumber.png",
#        device="png",dpi=320,units="cm",width=70,height=40)
# 
# #----------------------------------------------------------------------------#
# # If shorter genes recombine more, we can expect a GC3ex1 higher in shorter genes
# # GC3ex1 as a function of number of exons
# # Make a data frame with mean of GC3 exon 1 per number of exons and species
# # Only for 1 to 20 exons
# data_plot = data.frame(species = rep(unique(df_all$species), each = 20), number_exons = rep(1:20, times = length(unique(df_all$species))), GC3ex1 = NA, GC3ex1.upper = NA, GC3ex1.lower = NA)
# # Compute mean recombination and CI by bootstrap
# for (i in 1:nrow(data_plot)) {
#   cat("Row", i, "of", nrow(data_plot), "\n")
#   nboot = 1000
#   boot = numeric(nboot)
#   for (j in 1:nboot) {
#     boot[j] = mean(sample(df_all$GC3ex1[df_all$species == data_plot$species[i] & df_all$number_exons == data_plot$number_exons[i]], replace = TRUE), na.rm = TRUE)
#   }
#   data_plot$GC3ex1[i] = mean(boot)
#   data_plot$GC3ex1.upper[i] = quantile(boot, 0.975)
#   data_plot$GC3ex1.lower[i] = quantile(boot, 0.025)
# }
# # Save in a file, to avoid recomputing (bootstrapping is long)
# write.table(data_plot, file = paste("output/figures/GC3~exonnumber.txt", sep = ""),
#             quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
# # Import file
# data_plot = read.table(file = paste("output/figures/GC3~exonnumber.txt", sep = ""), header = TRUE, sep = "\t")
# # Add a small horizontal jitter
# data_plot$number_exons = data_plot$number_exons + rep(seq(0, 0.5, length.out = length(unique(data_plot$species))), each = 20)
# 
# # Plot
# GC3ex1_exonsnumber = ggplot(data = data_plot, aes(x = number_exons, y = GC3ex1, colour = species)) +
#   # One point per species and number of exon
#   geom_point() +
#   # Add confidence intervals
#   geom_errorbar(aes(ymin = GC3ex1.lower, ymax = GC3ex1.upper), width = 0) +
#   labs(x = "Number of exons", y="GC3 exon 1 proportion", fill="Species") +
#   theme(axis.line = element_line(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
#         plot.subtitle = element_text(color="black",size=28,hjust = 0.5),
#         axis.title.x = element_text(color="black", size=28),
#         axis.title.y = element_text(color="black", size=28),
#         axis.text=element_text(size=28, colour="black"),
#         strip.text=element_text(size=28, colour="black"),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.text=element_text(size=28),
#         legend.title=element_text(size=28)
#   )
# GC3ex1_exonsnumber
# ggsave("figures/GC3ex1~exonnumber.png",
#        device="png",dpi=320,units="cm",width=70,height=40)
# 
# 
# 
# 
# 
# #============================================================================#
# # Exon level
# #============================================================================#
# 
# # EXONS
# gcexons = read.table(file = gzfile(paste("data-cleaned/genome/gc_exons/", map, ".txt.gz", sep = "")), header = TRUE, sep = "\t")
# 
# 
# #----------------------------------------------------------------------------#
# # GC3 gradients for all species
# # Gather data for figures
# # Combine all species
# list_maps = c("Arabidopsis_thaliana_Serin2017",  "Oryza_sativa_DeLeon2016",
#               "Malus_domestica_DiPierro2016" ,  "Brachypodium_distachyon_Huo2011", "Cucumis_sativus_Zhu2016", "Gossypium_raimondii_Wang2013",
#               "Phaseolus_vulgaris_Song2015", "Sorghum_bicolor_Zou2012", "Zea_mays_MaizeGDBConsensus_v4", "Zea_mays_IBM_MaizeSNP50",
#               "Setaria_italica_Ni2017", "Solanum_lycopersicum_Gonda2018", "Triticum_aestivum_GutierrezGonzalez2019")
# 
# # Compute the mean GC3 and 95% CI for the gradient of GC3 along genes (up to the nth exon)
# # Take the name of the map in argument
# # Results are in a data frame with map, species, exon rank, mean GC3, upper bound of mean GC3, lower bound of mean GC3
# # map = list_maps[1]
# GC3_gradient = function(map = "", n_exons = 20) {
#   res = data.frame(map = rep(map, n_exons), species =  gsub("_", " ",gsub("[A-Za-z0-9]*$", "", map)), exon_rank = 1:n_exons, GC3 = NA, upper_GC3 = NA, lower_GC3 = NA)
#   # Import data for the map
#   cat("Loading", map, "\n")
#   gcexons = read.table(file = gzfile(paste("data-cleaned/genome/gc_exons/", map, ".txt.gz", sep = "")), header = TRUE, sep = "\t")
#   # Run a bootstrap of nboot iterations to compute the mean and 95% CI
#   cat("Bootstrapping...\n")
#   for (i in 1:nrow(res)) {
#     cat("Rank", i, "of", nrow(res), "\n")
#     nboot = 1000
#     boot = numeric(nboot)
#     for (j in 1:nboot) {
#       boot[j] = mean(sample(gcexons$gc3[gcexons$rank == res$exon_rank[i]], replace = TRUE), na.rm = TRUE)
#     }
#     # boot
#     res$GC3[i] = mean(boot)
#     res$upper_GC3[i] = quantile(boot, 0.975)
#     res$lower_GC3[i] = quantile(boot, 0.025)
#   }
#   # Return the results as a data frame
#   return(res)
# }
# # Run the gradient function for all species combined in a single dataset
# data_plot = GC3_gradient(list_maps[1])
# for (i in 2:length(list_maps)) {
#   df = GC3_gradient(list_maps[i])
#   data_plot = rbind(data_plot, df)
#   rm(df)
# }
# # Save in a file, to avoid recomputing (bootstrapping is long)
# write.table(data_plot, file = paste("output/figures/GC3_gradient.txt", sep = ""),
#             quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
# # Import file
# data_plot = read.table(file = paste("output/figures/GC3_gradient.txt", sep = ""), header = TRUE, sep = "\t")
# 
# # Plotting the GC3 gradient
# GC3_gradient = ggplot(data = data_plot, aes(x = exon_rank, y = GC3, colour = species)) +
#   # One point per species and number of exon
#   geom_point() +
#   # Add confidence intervals
#   geom_errorbar(aes(ymin = lower_GC3, ymax = upper_GC3), width = 0) +
#   geom_line() +
#   labs(x = "Exon rank", y="GC3 proportion", fill="Species") +
#   theme(axis.line = element_line(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
#         plot.subtitle = element_text(color="black",size=28,hjust = 0.5),
#         axis.title.x = element_text(color="black", size=28),
#         axis.title.y = element_text(color="black", size=28),
#         axis.text=element_text(size=28, colour="black"),
#         strip.text=element_text(size=28, colour="black"),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.text=element_text(size=28),
#         legend.title=element_text(size=28)
#   )
# GC3_gradient
# ggsave("figures/GC3_gradient.png",
#        device="png",dpi=320,units="cm",width=70,height=40)
# 
# 
# 
# #============================================================================#
# # Sampling effect (Sylvain Glémin)
# #============================================================================#
# # Salut Léo. Désolé je prends la conversation tard. Pour tes graphes du Nb d’exon en fonction du taux de GC tu as des choses non monotones étranges au premier abord car par prédit de façon évidente quand on sait que le GC décroît avec le rang de l’exon (cf stage de Maya et certains de tes graphes je crois). On s’attendrait donc à avoir une relation négative aussi entre taux de GC du gène et nb d’exons et pas ces courbes en cloche. En fait il s’agit d’un artefact. D’où l’importance de bien ploter ce qu’on veux tester et de faire attention à comment on regroupe les données. Ici tu veux tester que le taux de GC diminue avec la longueur des gènes (mesurée en nb d’exons). Donc ce que tu dois faire c’est le GC en fonction du nombre d’exons. Si tu fais ça tu verras que ça diminue sans doute chez toutes les espèces en étant plus plat chez les espèces pauvres en GC. Mais en tout cas je parie que tu ne trouveras pas une relation en cloche.
# # 22 h 09
# # D’où ça vient alors...
# # 22 h 11
# # Du fait que tu as beaucoup plus de gènes courts que de gènes longs. La distribution doit être plus ou moins exponentielle. Du coup comme il y a beaucoup plus de gène, la gamme de GC couverte est simplement plus grande. De ce fait quand tu regroupes sur le taux de GC, dans les valeurs extremes tu a beaucoup plus de chance de tomber sur des gènes courts que sur des gènes longs alors que pour les taux de GC intermédiaires tu peux tomber sur toutes les tailles.
# # 22 h 13
# # Pour t’en convaincre voici deux figures à partir de 3 cas simulés: aucun effet du nombre d’exons sur le GC  (bleu), un effet faiblement décroissant (rouge) et un effet fortement décroissant (violet)
# 
# set.seed(100)
# # Jeu de données simulé:
# # Distribution leptokurtique du nombre d'exon (comme observé la plus part du temps)
# NbExon <- c(rep(1,10000),rep(2,5000),rep(3,2500),rep(4,1250),rep(5,600),rep(6,300),rep(7,150),rep(8,75),rep(9,40),rep(10,20))
# # Aucun effet moyen du GC sur le nombre d'exon: on suppose une taille d'exon de 100 nucléotide et le taux de GC est tiré dans une binomial de paramètre n = 100*Nb exon et p = 0.5
# GC <- sapply(NbExon,FUN = function(x) rbinom(1,size = 100*x,prob = 0.5)/(100*x) )
# # Calcul du taux de GC moyen par nb d'exon
# GC_mean <- aggregate(x = GC, by = list(NbExon), FUN = mean)
# names(GC_mean) <- c("NbExon","gc")
# # Calcul du nb d'exon moyens par quantile de GC
# GCquant <- cut(GC,quantile(GC,probs = c(0:10)/10))
# NbExon_mean <- aggregate(x = list(GC,NbExon), by = list(GCquant), FUN = mean)
# names(NbExon_mean) <- c("Quantile","gc","NbExon")
# # Plot du GC en fonction du nb d'exons
# # Chaque gène (en noir) et la moyenne par nb d'exon (en rouge)
# plot(GC~NbExon)
# points(GC_mean$gc ~ GC_mean$NbExon, col = "red", pch = 16)
# # Plot du Nb d'exons moyen en fonction du GC
# plot(NbExon~GC)
# points(NbExon_mean$NbExon ~ NbExon_mean$gc, col = "red", pch = 16)
# # La même chose mais maintenant le taux de GC diminue faiblement avec le nombre d'exons
# GCbis <- sapply(NbExon,FUN = function(x) rbinom(1,size = 100*x,prob = 0.55*(1 - 0.01*x))/(100*x) )
# # Calcul du taux de GC moyen par nb d'exon
# GCbis_mean <- aggregate(x = GCbis, by = list(NbExon), FUN = mean)
# names(GCbis_mean) <- c("NbExon","gc")
# # Calcul du nb d'exon moyens par quantile de GC
# GCbisquant <- cut(GCbis,quantile(GCbis,probs = c(0:10)/10))
# NbExonbis_mean <- aggregate(x = list(GCbis,NbExon), by = list(GCbisquant), FUN = mean)
# names(NbExonbis_mean) <- c("Quantile","gc","NbExon")
# # Plot du GC en fonction du nb d'exons
# # Chaque gène (en noir) et la moyenne par nb d'exon (en rouge)
# plot(GCbis~NbExon)
# points(GCbis_mean$gc ~ GCbis_mean$NbExon, col = "red", pch = 16)
# # Plot du Nb d'exons moyen en fonction du GC
# plot(NbExon~GCbis)
# points(NbExonbis_mean$NbExon ~ NbExonbis_mean$gc, col = "red", pch = 16)
# # La même chose mais maintenant le taux de GC diminue fortement avec le nombre d'exons
# GCter <- sapply(NbExon,FUN = function(x) rbinom(1,size = 100*x,prob = 0.65*(1 - 0.05*x))/(100*x) )
# # Calcul du taux de GC moyen par nb d'exon
# GCter_mean <- aggregate(x = GCter, by = list(NbExon), FUN = mean)
# names(GCter_mean) <- c("NbExon","gc")
# # Calcul du nb d'exon moyens par quantile de GC
# GCterquant <- cut(GCter,quantile(GCter,probs = c(0:10)/10))
# NbExonter_mean <- aggregate(x = list(GCter,NbExon), by = list(GCterquant), FUN = mean)
# names(NbExonter_mean) <- c("Quantile","gc","NbExon")
# # Plot du GC en fonction du nb d'exons
# # Chaque gène (en noir) et la moyenne par nb d'exon (en rouge)
# plot(GCter~NbExon)
# points(GCter_mean$gc ~ GCter_mean$NbExon, col = "red", pch = 16)
# # Plot du Nb d'exons moyen en fonction du GC
# plot(NbExon~GCter)
# points(NbExonter_mean$NbExon ~ NbExonter_mean$gc, col = "red", pch = 16)
# # Comparaison des trois cas:
# # Nb exons en fonction du taux de GC (ce que tu as fait)
# plot(NbExon_mean$NbExon ~ NbExon_mean$gc,col="blue",pch=16,xlim=c(0.3,0.7),ylim=c(1,5),xlab="GC",ylab="NbExons")
# lines(NbExon_mean$NbExon ~ NbExon_mean$gc,col="blue")
# points(NbExonbis_mean$NbExon ~ NbExonbis_mean$gc,col="red",pch=16)
# lines(NbExonbis_mean$NbExon ~ NbExonbis_mean$gc,col="red")
# points(NbExonter_mean$NbExon ~ NbExonter_mean$gc,col="purple",pch=16)
# lines(NbExonter_mean$NbExon ~ NbExonter_mean$gc,col="purple")
# # Comparaison des trois cas:
# # GC en fonction du nb d'exons 
# plot(GC_mean$gc ~ GC_mean$NbExon,col="blue",pch=16,xlim=c(1,10),ylim=c(0.3,0.65),xlab="NbExons",ylab="GC")
# lines(GC_mean$gc ~ GC_mean$NbExon,col="blue")
# points(GCbis_mean$gc ~ GCbis_mean$NbExon,col="red",pch=16)
# lines(GCbis_mean$gc ~ GCbis_mean$NbExon,col="red")
# points(GCter_mean$gc ~ GCter_mean$NbExon,col="purple",pch=16)
# lines(GCter_mean$gc ~ GCter_mean$NbExon,col="purple")


