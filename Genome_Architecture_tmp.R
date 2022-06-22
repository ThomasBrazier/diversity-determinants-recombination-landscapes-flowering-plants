#============================================================================#
#                    ECOBIO - PhD
#
#                 Genome Architecture
#
#============================================================================#
# by Thomas Brazier
# brazier.thomas@gmail.com

# Supervisor: Sylvain GLEMIN
#             ECOBIO Lab


# Within this script you have functions to
# - retrieve gene/exon/CDS positions for a given map and reference genome
# - Compute genome statistics in windows (e.g. GC content, gene length)

#============================================================================#
# Loading environment ----
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
library(rBLAST)
library(Biostrings)
library(pals)

#https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html

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
source("sources/Genome.R") # Essential functions to parse GFF files and compute genomic statistics

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

# Or a custopm list of species
# list_species = sort(c("Oryza_sativa", "Oryza_indica", "Sesamum_indicum", "Triticum_aestivum", "Arabidopsis_thaliana",
#                       "Solanum_tuberosum", "Solanum_lycopersicum", "Brachypodium_distachyon", "Citrullus_lanatus", "Gossypium_hirsutum", "Gossypium_raimondii",
#                       "Malus_domestica", "Phaseolus_vulgaris", "Prunus_mume", "Theobroma_cacao", "Zea_mays"))
# list_species = c("Arabidopsis_thaliana")
# list_accession = c("GCA_902460315.1")
# chromosome = 2

#============================================================================#
# Chromosome characteristics - Chromosome names and size ----
#============================================================================#
# Make the list of all chromosome names in their NCBI annotation and litteral forms
# And get the size of the sequence at the same time

# Generate the list of NCBI annotation forms
metadata.clean = read.csv(file = paste(wd, "/data-cleaned/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")
chrnames_translationtable = data.frame(set = character(),
                                       accession = character(),
                                       annotname = character(),
                                       litteralname = character(),
                                       chrsize.bp = numeric(),
                                       annot = character())

getChrName = function(i) {
  cat(i, "on", nrow(metadata.clean), ":", metadata.clean$species[i], "\n")
  species = metadata.clean$species[i]
  acc = metadata.clean$accession[i]
  fastaname = system(paste("ls ", wd, "/data/Genome/", tolower(gsub(" ", "_", species)),
                           "/", acc, "/fasta/ | grep '.f[n]*a[sta]*.gz'", sep = ""), intern = TRUE)

  # Get chromosome names
  
  chrnames = system(paste("zcat ", wd, "/data/Genome/", tolower(gsub(" ", "_", species)),
                          "/", acc, "/fasta/", fastaname, " | grep '^>'", sep = ""), intern = TRUE)
  # Index which lines to keep (trim scaffolds)
  idx = which(unlist(lapply(chrnames, function(x) grepl(">CM", x))) |
                              unlist(lapply(chrnames, function(x) grepl(">NC", x))) |
                              unlist(lapply(chrnames, function(x) grepl(">BDM", x))) |
                              unlist(lapply(chrnames, function(x) grepl("^[chromse]*[0-9]+[A-Za-z]*$", x))))
  # Get chromosome sizes
  system(paste("gunzip ", wd, "/data/Genome/", tolower(gsub(" ", "_", species)),
               "/", acc, "/fasta/", fastaname, sep = ""), intern = TRUE)
  system(paste("bgzip ", wd, "/data/Genome/", tolower(gsub(" ", "_", species)),
               "/", acc, "/fasta/", gsub(".gz", "", fastaname), sep = ""), intern = TRUE)
  chrsizes = system(paste("faidx ", wd, "/data/Genome/", tolower(gsub(" ", "_", species)),
               "/", acc, "/fasta/", fastaname, " -i chromsizes | awk '{ print $2 }'", sep = ""), intern = TRUE)

  # Trim chromosome names and sequence size by index
  chrnames = chrnames[idx]
  chrsizes = chrsizes[idx]
  # Separate chromosome name and full sequence header/metadata
  annot = chrnames
  chrnames = unlist(lapply(chrnames, function(x) unlist(strsplit(x, " "))[1]))
  
  # Make a data frame of results
  annot = data.frame(set = metadata.clean$id[i],
                     accession = acc,
                     annotname = chrnames,
                     litteralname = NA,
                     chrsize.bp = chrsizes,
                     annot = annot)
  return(annot)
}
# i = 19
# chrnames = getChrName(i)

# BACTH
for (i in 1:nrow(metadata.clean)) {
  chrnames_translationtable = rbind(chrnames_translationtable,
                                    getChrName(i))
}
# Save file
write.table(chrnames_translationtable, "data-cleaned/genome/chrnames_translationtable.csv",
            colnames = TRUE, rownames = FALSE, sep = "\t", quote = FALSE)
chrnames_translationtable = read.table("data-cleaned/genome/chrnames_translationtable.csv",
                                       header = TRUE)




#============================================================================#
# C-values per chromosome (chromosome sizes in pg) ----
#============================================================================#
# Kew Garden database
Cvalues = read.csv("data-cleaned/CValueDataBase.csv", header = TRUE, sep = ";")
metadata = read.table("data/Genetic_maps/Genetic_maps_ressources.csv", header = TRUE, sep = ";")

# Retrieve a vector for species and c values
sp = metadata$species
cval = metadata$genome_Cvalue
HCN = metadata$hcn
ploidy = metadata$ploidy
# Get the Cvalue for each species
idx_species = paste(Cvalues$Genus, Cvalues$Species, sep = "_")
# Arbitrary choose the first value returned if multiple results
for (i in 1:length(sp)) {
  print(sp[i])
  if (length(Cvalues$Cval[which(idx_species == as.character(sp[i]))] > 0)) {
    cval[i] = (Cvalues$Cval[which(idx_species == as.character(sp[i]))])[1]
    print(cval[i])
    if (is.na(HCN[i])) {
      HCN[i] = (Cvalues$ChromosomeNumber[which(idx_species == as.character(sp[i]))])[1] / (Cvalues$PloidyLevel[which(idx_species == as.character(sp[i]))])[1]
    }
    if (is.na(ploidy[i])) {
      (ploidy[i] = Cvalues$PloidyLevel[which(idx_species == as.character(sp[i]))])[1]
    }
  }
}
# Update the metadata with HCN, ploidy level and C-values
metadata$genome_Cvalue = cval
metadata$hcn = HCN
metadata$ploidy = ploidy

# Inspect results
# tmp = metadata[,-c(2:5, 7:18, 23:26)]

write.table(metadata, "data/Genetic_maps/Genetic_maps_ressources.csv", col.names = TRUE,
            row.names = FALSE, quote = FALSE, sep = ";")

#============================================================================#

# Genome features positions -----

#============================================================================#

# One file per features, for all species
# Files saved in 'data-cleaned/genomes'
# Columns are
# (1) species
# (2) accession
# (3) chromosome
# (4) start
# (5) end
# (6) sequence name
# (7) score of the annotated sequence

#============================================================================#
# Gene positions ----
#============================================================================#
# Data saved in 'gene.positions.txt'
# All species and accessions in the same file

# Retrieve gene positions in GFF files and save in a tidy format
source("sources/Genome.R") # Essential functions to parse GFF files and compute genomic statistics

# A list of accession to process
metadata.clean = read.table("data-cleaned/Genetic_maps_ressources.csv", header = TRUE, sep = ";")
metadata.clean$accession
metadata.clean$id[is.na(metadata.clean$accession)]
list_accessions = metadata.clean$accession
list_accessions = as.character(list_accessions[!is.na(list_accessions)])
# list_accessions[c(3, 4, 12, 14, 22, 24, 25, 26, 27)] # GFF to repair
# list_accessions = list_accessions[-c(3, 12, 14, 25, 27)] # GFF not working/available for some accessions
list_accessions
# Metadata.clean contains all reference genomes to process, associated to recombination maps
metadata.clean$species
i = 3
(list_species = as.character(metadata.clean$species[i]))
(list_accessions = as.character(metadata.clean$accession[i]))
species = as.character(metadata.clean$species[i])

# BATCH PROCESS
# Process each genome dataset and save in separate files per species
source("sources/Genome.R")
batch.gene.position(list_accessions)

# Explore the dataset
genpos = read.table(gzfile(paste("data-cleaned/genome/gene_positions/", species,"_gene_positions.txt.gz", sep = "")), header = TRUE, sep = "\t")
View(genpos)

# Number of species with an annotation file and gene positions
list = system(paste("ls ", wd, "/data-cleaned/genome/gene_positions/", sep = ""), intern = TRUE)
list = gsub("_gene_positions.txt.gz", "", list)
length(list)
list

# c(3,33,17,4,11,16,27,29,44,48,38,32,43,56)
# How many alternative genome annotation to re-map?
# 33 - Nelumbo_nucifera GCA_003033695.1 None annotation for this version; annotation in GCF_000365185.1
# 3 - Capsicum, no gene annotation in the actual gff -> go to GCF_000710875.1 (latest)
# 17 - Raphanus_sativus GCA_002197605.1 None annotation -> scaffolds annotated in GCF_000801105.1


# Indexes with problems --> Problem avec biopython.convert when converting gbff to gff
# NO FEATURE FOUND IN DATA
# 4 - Citrullus_lanatus GCA_000238415.2, none annotation
# 11 - Mangifera_indica GCA_011075055.1, none annotation
# 16 - Quercus_sp GCA_900291515.1, none annotation
# 27 - Triticum_dicoccoides GCA_002575655.1 None annotation

# Otherwise data not downloaded or not a regular file
# 29 - Cenchrus_americanus GCA_002174835.2 None annotation
# 44 - Momordica_charantia GCA_013281855.1 None annotation
# 48 - Draba_nivalis PRJNA657155
# 38 - Triticum_urartu	GCA_003073215.1 None annotation

# No data
# 32 - Coffea canephora
# 43 - Boechera_stricta NA
# 56 - Juglans regia NA

# RESOLVED with success
# 30 - Brassica napus GCF_000686985.2
# 15 - Prunus_mume GCA_000346735.1
# 7 - Cucurbita_pepo GCA_002806865.2
# 18 - Sesamum_indicum GCA_000512975.1
# 19 - Solanum_tuberosum PGSC_DM_v4.03
# 41 - Prunus_persica GCA_000346465.2
# 45 - Elaeis_guineensis GCA_000442705.1
# 46 - Vigna_unguiculata GCA_004118075.1 
# 49 - Gossypium_hirsutum HAU_G.hirsutum_AD1genome_v1.1
# 57 - Cucurbita_maxima Cmaxima_v1.1
# 52 - Citrus_sinensis GCA_000695605.1 Assembly level to scaffold, not chromosome, but LGs correspond to nine first scaffolds
# 51 - Capsella_rubella GCA_000375325.1
# 42 - Arachis_hypogaea GCA_003713155.1 
# 39 - Glycine_max GCA_000004515.4
# 53 - Camelina_sativa GCA_000633955.1
# 37 - Arachis_duranensis GCF_000817695.2


#============================================================================#
# Exon positions ----
#============================================================================#
#Process a single dataset
source("sources/Genome.R")
# metadata.clean$species
# i = 28
# (species = as.character(metadata.clean$species[i]))
# (accession = as.character(metadata.clean$accession[i]))
(species = as.character(metadata.clean$species[which(metadata.clean$id == map)]))
(accession = as.character(metadata.clean$accession[which(metadata.clean$id == map)]))

exonpos = exon.position(species, accession)
write.table(exonpos, gzfile(paste("data-cleaned/genome/exon_positions/", species,"_exon_positions.txt.gz", sep = "")), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)


# Process each genome dataset and save in separate files per species
for (acc in 1:length(list_accessions)) { # LOOP ON ACCESSIONS
  # Getting infos on dataset to treat
  accession = list_accessions[acc]
  species = list_species[acc]
  cat(species, accession, "\n")
  # Retrieving feature position, i.e. genes
  exonpos = exon.position(species, accession)
  # Some data is not available or no features have been reported in annotation files (bad annotation)
  if (nrow(exonpos) > 0) {
    # if (((exonpos) != "Accession is not available in your dataset.")[1]) {
    # Add species and accession to the results
    exonpos = cbind(data.frame(species = rep(species, nrow(exonpos)), accession = rep(accession, nrow(exonpos))),
                    exonpos)
    #Load the dataset of the species if it exists, otherwise create a new one
    if (file.exists(paste("data-cleaned/genome/exon_positions/", species,"_exon_positions.txt.gz", sep = ""))) {
      # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
      dataset = read.table(gzfile(paste("data-cleaned/genome/exon_positions/", species,"_exon_positions.txt.gz", sep = "")),
                           header = TRUE, sep = "\t")
    } else {
      # Template of a new empty dataset
      dataset = data.frame(species = character(0), accession = character(0), id = character(0),
                           type = character(0), chromosome = character(0),
                           strand = character(0), rank = numeric(0), start = numeric(0), end = numeric(0))
      write.table(dataset, gzfile(paste("data-cleaned/genome/exon_positions/", species,"_exon_positions.txt.gz", sep = "")), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)      # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
      # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
      dataset = read.table(gzfile(paste("data-cleaned/genome/exon_positions/", species,"_exon_positions.txt.gz", sep = "")),
                           header = TRUE, sep = "\t")
    }
    
    # Remove older rows in the dataset
    dataset = dataset[!(dataset$accession == accession),]
    # Bind new rows
    dataset = rbind(dataset, exonpos)
    # Init an empty dataset
    # dataset = data.frame(species = character(0), accession = character(0), id = character(0),
    #                      type = character(0), chromosome = character(0), strand = character(0),
    #                      start = character(0), end = character(0)) 
    # Write updated dataset in 'data-cleaned'
    write.table(dataset, gzfile(paste("data-cleaned/genome/exon_positions/", species,"_exon_positions.txt.gz", sep = "")), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)
    rm(dataset)
    # } else {
    #   warning("Accession is not available in your dataset.")
    # }
  } else {
    warning(paste(species, accession,"No feature found in data.", sep = " "))
  }
  # END OF ACCESSION LOOP
}

# Explore a species dataset
species = ""
dataset = read.table(gzfile(paste("data-cleaned/genome/exon_positions/", species,"_exon_positions.txt.gz", sep = "")), header = TRUE, sep = "\t")

#============================================================================#
# CDS positions ----
#============================================================================#
# Retrieve gene positions in GFF files and save in a tidy format
source("sources/Genome.R") # Essential functions to parse GFF files and compute genomic statistics

# A list of accession to process
metadata.clean = read.table("data-cleaned/Genetic_maps_ressources.csv", header = TRUE, sep = ";")
metadata.clean$accession
metadata.clean$id[is.na(metadata.clean$accession)]
list_accessions = metadata.clean$accession
list_accessions = as.character(list_accessions[!is.na(list_accessions)])
list_accessions
# Metadata.clean contains all reference genomes to process, associated to recombination maps
metadata.clean$species
i = 38
(list_species = as.character(metadata.clean$species[i]))
(list_accessions = as.character(metadata.clean$accession[i]))
species = as.character(metadata.clean$species[i])

# BATCH PROCESS
# Process each genome dataset and save in separate files per species
source("sources/Genome.R")
batch.cds.position(list_accessions)

# Explore the dataset
cdspos = read.table(gzfile(paste("data-cleaned/genome/cds_positions/", species,"_cds_positions.txt.gz", sep = "")), header = TRUE, sep = "\t")
View(cdspos)

# Number of species with an annotation file and gene positions
list = system(paste("ls ", wd, "/data-cleaned/genome/cds_positions/", sep = ""), intern = TRUE)
list = gsub("_cds_positions.txt.gz", "", list)
length(list)
list

# # Process each genome dataset and save in separate files per species
# for (acc in 1:length(list_accessions)) { # LOOP ON ACCESSIONS
#   # Getting infos on dataset to treat
#   accession = list_accessions[acc]
#   species = list_species[acc]
#   cat(species, accession, "\n")
#   # Retrieving feature position, i.e. genes
#   cdspos = cds.position(species, accession)
#   # Some data is not available or no features have been reported in annotation files (bad annotation)
#   if (nrow(cdspos) > 0) {
#     # if (((cdspos) != "Accession is not available in your dataset.")[1]) {
#       # Add species and accession to the results
#       cdspos = cbind(data.frame(species = rep(species, nrow(cdspos)), accession = rep(accession, nrow(cdspos))),
#                      cdspos)
#       #Load the dataset of the species if it exists, otherwise create a new one
#       if (file.exists(paste("data-cleaned/genome/CDS_positions/", species,"_CDS_positions.txt.gz", sep = ""))) {
#         # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
#         dataset = read.table(gzfile(paste("data-cleaned/genome/CDS_positions/", species,"_CDS_positions.txt.gz", sep = "")),
#                              header = TRUE, sep = "\t")
#       } else {
#         # Template of a new empty dataset
#         dataset = data.frame(species = character(0), accession = character(0), id = character(0),
#                              type = character(0), chromosome = character(0),
#                              strand = character(0), phase = numeric(0), start = numeric(0), end = numeric(0))
#         write.table(dataset, gzfile(paste("data-cleaned/genome/CDS_positions/", species,"_CDS_positions.txt.gz", sep = "")), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)      # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
#         # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
#         dataset = read.table(gzfile(paste("data-cleaned/genome/CDS_positions/", species,"_CDS_positions.txt.gz", sep = "")),
#                              header = TRUE, sep = "\t")
#       }
#       
#       # Remove older rows in the dataset
#       dataset = dataset[!(dataset$accession == accession),]
#       # Bind new rows
#       dataset = rbind(dataset, cdspos)
#       # Init an empty dataset
#       # dataset = data.frame(species = character(0), accession = character(0), id = character(0),
#       #                      type = character(0), chromosome = character(0), strand = character(0),
#       #                      start = character(0), end = character(0)) 
#       # Write updated dataset in 'data-cleaned'
#       write.table(dataset, gzfile(paste("data-cleaned/genome/CDS_positions/", species,"_CDS_positions.txt.gz", sep = "")), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)
#       rm(dataset)
#     # } else {
#     #   warning("Accession is not available in your dataset.")
#     # }
#   } else {
#     warning(paste(species, accession,"No feature found in data.", sep = " "))
#   }
#   # END OF ACCESSION LOOP
# }

# Explore a species dataset
# species = ""
# dataset = read.table(gzfile(paste("data-cleaned/genome/CDS_positions/", species,"_CDS_positions.txt.gz", sep = "")), header = TRUE, sep = "\t")

#============================================================================#
#   Gene length ----
#============================================================================#
# Compute the length of each gene (stop - start),
# its coding length (sum of all CDS lengths),
# and the estimated recombination rate (the windows rec. rate)
# For each species with gene positions

estimate_gene_length = function() {
  # For each species with gene positions
  list_species = system(paste("ls ", wd, "/data-cleaned/genome/gene_positions/", sep = ""), intern = TRUE)
  list_species = gsub("_gene_positions.txt.gz", "", list_species)
  length(list_species)
  list_species
  # Import gene positions
  readall_genepos = function(x) {
    data = read.table(file = gzfile(paste("data-cleaned/genome/gene_positions/", x, "_gene_positions.txt.gz", sep = "")),
                      colClasses = c(rep("character", 7), rep("numeric", 2)),
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    return(data)
  }
  df_genepos = bind_rows(lapply(list_species, function(x) readall_genepos(x)))
  
  # Import also CDS positions
  readall_cdspos = function(x) {
    data = read.table(file = gzfile(paste("data-cleaned/genome/cds_positions/", x, "_cds_positions.txt.gz", sep = "")),
                      colClasses = c(rep("character", 7), rep("numeric", 2)),
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    return(data)
  }
  df_cdspos = bind_rows(lapply(list_species, function(x) readall_cdspos(x)))
  
  # Combine genes and CDS for further same procedure
  df = rbind(df_genepos, df_cdspos)
  rm(df_genepos)
  rm(df_cdspos)
  
  # Load recombination maps
  # List of all maps with a gene position annotation
  list_maps = system(paste("ls ", wd, "/output/recombination_maps/loess/100kbwind/", sep = ""), intern = TRUE)
  list_maps = list_maps[unlist(lapply(list_species, function (x) grep(x, list_maps)))]
  readall_recmaps = function(x) {
    data = read.table(file = gzfile(paste("output/recombination_maps/loess/100kbwind/", x, sep = "")),
                      colClasses = rep("numeric", 4),
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    data$chromosome = gsub(".txt", "", gsub(".*_chromosome", "", x))
    data$species = gsub("_MaizeGDBConsensus", "", gsub("_[A-Za-z0-9]*_[A-Za-z0-9]*.txt", "", x))
    return(data)
  }
  df_recmaps = bind_rows(lapply(list_maps, function(x) readall_recmaps(x)))
  
  # Now that data is built
  # Compute gene/CDS length (stop - start)
  df$length = df$end - df$start
  
  # Retrieve the recombination rate of the windows which contains the gene/CDS 
  # i.e. the window containing the start position
  df$recrate = as.numeric(NA)
  est_recrate = function(x) {
    x = df[x,]
    recrate = df_recmaps$rec.rate[which((x$start/1000000) >= (df_recmaps$phys-0.05) &
                                          (x$start/1000000) <= (df_recmaps$phys+0.05) &
                                          df_recmaps$species == x$species & df_recmaps$chromosome == x$chromosome)]
    recrate = recrate[1] # Take the first if multiple results
    if (length(recrate) == 0) {recrate = NA}
    return(recrate)
  }
  
  # for (i in 1:nrow(df)) {
  #   print(i/nrow(df)*100)
  #   df$recrate[i] = est_recrate(i)
  # }
  # x = df[2985955,]
  # x = df[1,]
  # est_recrate(x)
  library(pbmcapply)
  df$recrate = unlist(pbmclapply(1:nrow(df), function(x) est_recrate(x)))
  
  # Save the file
  write.table(df, file = gzfile("output/gene_length.txt.gz"), sep = "\t",
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  return(df)
}

estimate_gene_length()
df = read.table(file = gzfile("output/gene_length.txt.gz"), sep = "\t", header = TRUE)





#============================================================================#
#   GC content ----
#============================================================================#
# Global composition in GC, in windows of 100kb and 1Mb of same coordinates as recombination maps
# Composition in 1rst, 2nd and 3rd codon position in the CDS contained by the 100kb window (i.e. CDS with start position in the window)
# In order to be explanatory variables for recombination landscapes, windows need to be the same as used for recombination rates
# Hence, importing each recombination maps to compute GC proportions within defined windows

# Enhance data with a new species
list_species = as.character(metadata.clean$species[38])
list_accessions = as.character(metadata.clean$accession[38])
# Process each genome dataset and save in separate files per species
# Gene position
for (acc in 1:length(list_accessions)) { # LOOP ON ACCESSIONS
  # Getting infos on dataset to treat
  accession = list_accessions[acc]
  species = list_species[acc]
  cat(species, accession, "\n")
  # Retrieving feature position, i.e. genes
  genpos = gene.position(species, accession)
  # Some data is not available or no features have been reported in annotation files (bad annotation)
  if (nrow(genpos) > 0) {
    # if (((genpos) != "Accession is not available in your dataset.")[1]) {
    # Add species and accession to the results
    genpos = cbind(data.frame(species = rep(species, nrow(genpos)), accession = rep(accession, nrow(genpos))),
                   genpos)
    #Load the dataset of the species if it exists, otherwise create a new one
    if (file.exists(paste("data-cleaned/genome/gene_positions/", species,"_gene_positions.txt.gz", sep = ""))) {
      # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
      dataset = read.table(gzfile(paste("data-cleaned/genome/gene_positions/", species,"_gene_positions.txt.gz", sep = "")),
                           header = TRUE, sep = "\t")
    } else {
      # Template of a new empty dataset
      dataset = data.frame(species = character(0), accession = character(0), id = character(0),
                           type = character(0), chromosome = character(0),
                           strand = character(0), phase = numeric(0), start = numeric(0), end = numeric(0))
      write.table(dataset, gzfile(paste("data-cleaned/genome/gene_positions/", species,"_gene_positions.txt.gz", sep = "")), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)      # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
      # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
      dataset = read.table(gzfile(paste("data-cleaned/genome/gene_positions/", species,"_gene_positions.txt.gz", sep = "")),
                           header = TRUE, sep = "\t")
    }
    
    # Remove older rows in the dataset
    dataset = dataset[!(dataset$accession == accession),]
    # Bind new rows
    dataset = rbind(dataset, genpos)
    # Init an empty dataset
    # dataset = data.frame(species = character(0), accession = character(0), id = character(0),
    #                      type = character(0), chromosome = character(0), strand = character(0),
    #                      start = character(0), end = character(0)) 
    # Write updated dataset in 'data-cleaned'
    write.table(dataset, gzfile(paste("data-cleaned/genome/gene_positions/", species,"_gene_positions.txt.gz", sep = "")), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)
    rm(dataset)
    # } else {
    #   warning("Accession is not available in your dataset.")
    # }
  } else {
    warning(paste(species, accession,"No feature found in data.", sep = " "))
  }
  # END OF ACCESSION LOOP
}
# Exon position
for (acc in 1:length(list_accessions)) { # LOOP ON ACCESSIONS
  # Getting infos on dataset to treat
  accession = list_accessions[acc]
  species = list_species[acc]
  cat(species, accession, "\n")
  # Retrieving feature position, i.e. genes
  exonpos = exon.position(species, accession)
  # Some data is not available or no features have been reported in annotation files (bad annotation)
  if (nrow(exonpos) > 0) {
    # if (((exonpos) != "Accession is not available in your dataset.")[1]) {
    # Add species and accession to the results
    exonpos = cbind(data.frame(species = rep(species, nrow(exonpos)), accession = rep(accession, nrow(exonpos))),
                    exonpos)
    #Load the dataset of the species if it exists, otherwise create a new one
    if (file.exists(paste("data-cleaned/genome/exon_positions/", species,"_exon_positions.txt.gz", sep = ""))) {
      # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
      dataset = read.table(gzfile(paste("data-cleaned/genome/exon_positions/", species,"_exon_positions.txt.gz", sep = "")),
                           header = TRUE, sep = "\t")
    } else {
      # Template of a new empty dataset
      dataset = data.frame(species = character(0), accession = character(0), id = character(0),
                           type = character(0), chromosome = character(0),
                           strand = character(0), rank = numeric(0), start = numeric(0), end = numeric(0))
      write.table(dataset, gzfile(paste("data-cleaned/genome/exon_positions/", species,"_exon_positions.txt.gz", sep = "")), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)      # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
      # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
      dataset = read.table(gzfile(paste("data-cleaned/genome/exon_positions/", species,"_exon_positions.txt.gz", sep = "")),
                           header = TRUE, sep = "\t")
    }
    
    # Remove older rows in the dataset
    dataset = dataset[!(dataset$accession == accession),]
    # Bind new rows
    dataset = rbind(dataset, exonpos)
    # Init an empty dataset
    # dataset = data.frame(species = character(0), accession = character(0), id = character(0),
    #                      type = character(0), chromosome = character(0), strand = character(0),
    #                      start = character(0), end = character(0)) 
    # Write updated dataset in 'data-cleaned'
    write.table(dataset, gzfile(paste("data-cleaned/genome/exon_positions/", species,"_exon_positions.txt.gz", sep = "")), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)
    rm(dataset)
    # } else {
    #   warning("Accession is not available in your dataset.")
    # }
  } else {
    warning(paste(species, accession,"No feature found in data.", sep = " "))
  }
  # END OF ACCESSION LOOP
}
# CDS position
for (acc in 1:length(list_accessions)) { # LOOP ON ACCESSIONS
  # Getting infos on dataset to treat
  accession = list_accessions[acc]
  species = list_species[acc]
  cat(species, accession, "\n")
  # Retrieving feature position, i.e. genes
  cdspos = cds.position(species, accession)
  # Some data is not available or no features have been reported in annotation files (bad annotation)
  if (nrow(cdspos) > 0) {
    # if (((cdspos) != "Accession is not available in your dataset.")[1]) {
    # Add species and accession to the results
    cdspos = cbind(data.frame(species = rep(species, nrow(cdspos)), accession = rep(accession, nrow(cdspos))),
                   cdspos)
    #Load the dataset of the species if it exists, otherwise create a new one
    if (file.exists(paste("data-cleaned/genome/CDS_positions/", species,"_CDS_positions.txt.gz", sep = ""))) {
      # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
      dataset = read.table(gzfile(paste("data-cleaned/genome/CDS_positions/", species,"_CDS_positions.txt.gz", sep = "")),
                           header = TRUE, sep = "\t")
    } else {
      # Template of a new empty dataset
      dataset = data.frame(species = character(0), accession = character(0), id = character(0),
                           type = character(0), chromosome = character(0),
                           strand = character(0), phase = numeric(0), start = numeric(0), end = numeric(0))
      write.table(dataset, gzfile(paste("data-cleaned/genome/CDS_positions/", species,"_CDS_positions.txt.gz", sep = "")), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)      # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
      # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
      dataset = read.table(gzfile(paste("data-cleaned/genome/CDS_positions/", species,"_CDS_positions.txt.gz", sep = "")),
                           header = TRUE, sep = "\t")
    }
    
    # Remove older rows in the dataset
    dataset = dataset[!(dataset$accession == accession),]
    # Bind new rows
    dataset = rbind(dataset, cdspos)
    # Init an empty dataset
    # dataset = data.frame(species = character(0), accession = character(0), id = character(0),
    #                      type = character(0), chromosome = character(0), strand = character(0),
    #                      start = character(0), end = character(0)) 
    # Write updated dataset in 'data-cleaned'
    write.table(dataset, gzfile(paste("data-cleaned/genome/CDS_positions/", species,"_CDS_positions.txt.gz", sep = "")), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)
    rm(dataset)
    # } else {
    #   warning("Accession is not available in your dataset.")
    # }
  } else {
    warning(paste(species, accession,"No feature found in data.", sep = " "))
  }
  # END OF ACCESSION LOOP
}

# Then compute a new map
map = as.character(metadata.clean$id[38])
gc.map(map, ex = TRUE, method = 'gene')
gcgenes = read.table(file = gzfile(paste("data-cleaned/genome/gc_genes/", map, ".txt.gz", sep = "")), header = TRUE, sep = "\t")




# Then compute GC content
metadata.clean = read.table("data-cleaned/Genetic_maps_ressources.csv", header = TRUE, sep = ";")
list_maps = unique(metadata.clean$id)
# Process each recombination map with a genome dataset
for (map in list_maps[2]) {
  gc.map(map)
}

# EXAMPLE
# Problems with features: chromosomes names (e.g. CM014049.1 for Malus domestica), type -> character(0)

# Works with Arabidopsis, Oryza, Glycine max (long computational time)
map = "Arabidopsis_thaliana_Serin2017" # WORKING; DONE
map = "Oryza_sativa_DeLeon2016" # WORKING; DONE
map = "Malus_domestica_DiPierro2016" # WORKING; DONE
map = "Brachypodium_distachyon_Huo2011" # WORKING; DONE
map = "Oryza_sativa_Jiang2017" # WORKING; DONE
# map = "Capsicum_annuum_Han2016" # WORKING; DONE
map = "Cucumis_sativus_Zhu2016" # WORKING; DONE
map = "Gossypium_raimondii_Wang2013" # WORKING; DONE
map = "Phaseolus_vulgaris_Song2015" # TODO WORKING; DONE, beware... "Erreur : subscript contains out-of-bounds indices " -> verify GFF/Fasta congruence on accession
map = "Sorghum_bicolor_Zou2012" # TODO "Erreur : subscript contains out-of-bounds indices " -> verify GFF/Fasta congruence on accession
# map = "Theobroma_cacao_Royaert2016" # TODO "Erreur : subscript contains out-of-bounds indices " for method 'gene'  WORKING; DONE
map = "Zea_mays_MaizeGDBConsensus_v4" # WORKING; DONE
map = "Zea_mays_IBM_MaizeSNP50" # WORKING; DONE  Error on numerals
map = "Setaria_italica_Ni2017" # WORKING; DONE
map = "Solanum_lycopersicum_Gonda2018" # WORKING; DONE
map = "Triticum_aestivum_GutierrezGonzalez2019" #  WORKING; DONE, Very long computation, crashed sometimes
map = "Solanum_tuberosum_Endelman2016"
# map = list_maps[18]
# map = list_maps[16]
source("sources/Genome.R") # Essential functions to parse GFF files and compute genomic statistics

# GENES
gc.map(map, ex = TRUE, method = 'gene')
gcgenes = read.table(file = gzfile(paste("data-cleaned/genome/gc_genes/", map, ".txt.gz", sep = "")), header = TRUE, sep = "\t")
# View(gcgenes)
gcgenes$rec.rate[gcgenes$rec.rate < 0] = 0
plot(gcgenes$rec.rate, gcgenes$gc)
plot(gcgenes$rec.rate, gcgenes$gc3)
cor.test(gcgenes$rec.rate, gcgenes$gc3, method = "spearman")
# Downsampling to 1% of data to reduce autocorrelation
idx = sample(1:27445, 275)
cor.test(gcgenes$rec.rate[idx], gcgenes$gc3[idx], method = "spearman")

# EXONS
gc.map(map, ex = TRUE, method = 'exon')
gcexons = read.table(file = gzfile(paste("data-cleaned/genome/gc_exons/", map, ".txt.gz", sep = "")), header = TRUE, sep = "\t")
plot(gcexons$rec.rate, gcexons$gc)
plot(gcexons$rec.rate, gcexons$gc3)
cor.test(gcexons$rec.rate, gcexons$gc3, method = "spearman")

# RECOMBINATION MAPS
gc.map(map, ex = TRUE)
gcmap = read.table(file = gzfile(paste("data-cleaned/genome/gc_maps/", map, "_chromosome1.txt.gz", sep = "")), header = TRUE, sep = "\t")
# gcmap = mapdata
gcmap$rec.rate[gcmap$rec.rate < 0] = 0

# Figures of nucleotide landscapes
plot(gcmap$phys, gcmap$gc, ylim = c(0, 1), xlab = "Position along the chromosome (Mb)", ylab = "GC proportion (whole sequence)")
par(mfrow = c(1, 3))
plot(gcmap$phys, gcmap$gc1, ylim = c(0, 1), xlab = "Position along the chromosome (Mb)", ylab = "GC1 proportion (CDS)")
plot(gcmap$phys, gcmap$gc2, ylim = c(0, 1), xlab = "Position along the chromosome (Mb)", ylab = "GC2 proportion (CDS)")
plot(gcmap$phys, gcmap$gc3, ylim = c(0, 1), xlab = "Position along the chromosome (Mb)", ylab = "GC3 proportion (CDS)")
par(mfrow = c(1, 1))

plot(gcmap$gc1ex1, ylim = c(0, 1))
plot(gcmap$gc2ex1, ylim = c(0, 1))
plot(gcmap$gc3ex1, ylim = c(0, 1))
apply(cbind(gcmap$gc, gcmap$gc1, gcmap$gc2, gcmap$gc3, gcmap$gc1ex1, gcmap$gc2ex1, gcmap$gc3ex1), 2, function(x) mean(x, na.rm = TRUE))
# Figures of nucleotide landscapes
plot(gcmap$gc)
plot(gcmap$gc1)
plot(gcmap$gc2)
plot(gcmap$gc3)
plot(gcmap$gc1ex1)
plot(gcmap$gc2ex1)
plot(gcmap$gc3ex1)
# Relationship between GC content and recombination rate
plot(gcmap$rec.rate, gcmap$gc)
plot(gcmap$rec.rate, gcmap$gc1)
plot(gcmap$rec.rate, gcmap$gc2)
plot(gcmap$rec.rate, gcmap$gc3)
plot(gcmap$rec.rate, gcmap$gc1ex1)
plot(gcmap$rec.rate, gcmap$gc2ex1)
plot(gcmap$rec.rate, gcmap$gc3ex1)
# Correlation tests
cor.test(gcmap$rec.rate, gcmap$gc, method = "spearman")
cor.test(gcmap$rec.rate, gcmap$gc1, method = "spearman")
cor.test(gcmap$rec.rate, gcmap$gc2, method = "spearman")
cor.test(gcmap$rec.rate, gcmap$gc3, method = "spearman")
# GC proportion of first exon only
cor.test(gcmap$rec.rate, gcmap$gc1ex1, method = "spearman")
cor.test(gcmap$rec.rate, gcmap$gc2ex1, method = "spearman")
cor.test(gcmap$rec.rate, gcmap$gc3ex1, method = "spearman")



png(paste("output/comm/GC_", map,"_chr1.png",sep = ""), width = 800, height = 600)
plot(gcmap$phys, gcmap$gc, ylim = c(0, 1),
     xlab = "Position along the chromosome (Mb)", ylab = "GC proportion (whole sequence)", cex.lab = 1.5, cex.axis = 2)
dev.off()

png(paste("output/comm/GC_CDS_", map,"_chr1.png",sep = ""), width = 800, height = 600)
par(mfrow = c(1, 3))
plot(gcmap$phys, gcmap$gc1, ylim = c(0, 1), xlab = "Position along the chromosome (Mb)", ylab = "GC1 proportion (CDS)", cex.lab = 1.5, cex.axis = 2)
plot(gcmap$phys, gcmap$gc2, ylim = c(0, 1), xlab = "Position along the chromosome (Mb)", ylab = "GC2 proportion (CDS)", cex.lab = 1.5, cex.axis = 2)
plot(gcmap$phys, gcmap$gc3, ylim = c(0, 1), xlab = "Position along the chromosome (Mb)", ylab = "GC3 proportion (CDS)", cex.lab = 1.5, cex.axis = 2)
par(mfrow = c(1, 1))
dev.off()

# TEST
# mysequence <- s2c("agtctggggggccccttttaagtagatagatagctagtcgta")
# GC(mysequence)  # 0.4761905
# GC1(mysequence) # 0.6428571
# GC2(mysequence) # 0.3571429
# GC3(mysequence) # 0.4285714



# Batch mapping
# List of maps
list_maps = c("Arabidopsis_thaliana_Serin2017", # WORKING; DONE
              "Oryza_sativa_DeLeon2016", # WORKING; DONE
              "Malus_domestica_DiPierro2016", # WORKING; DONE
              "Brachypodium_distachyon_Huo2011", # WORKING; DONE
              "Oryza_sativa_Jiang2017", # WORKING; DONE
              # "Capsicum_annuum_Han2016", # TODO "Erreur : subscript contains out-of-bounds indices " # WORKING; DONE
              "Cucumis_sativus_Zhu2016", # WORKING; DONE
              "Gossypium_raimondii_Wang2013", # WORKING; DONE
              "Phaseolus_vulgaris_Song2015", # TODO WORKING; DONE, beware... "Erreur : subscript contains out-of-bounds indices " -> verify GFF/Fasta congruence on accession
              "Sorghum_bicolor_Zou2012", # TODO "Erreur : subscript contains out-of-bounds indices " -> verify GFF/Fasta congruence on accession
              # "Theobroma_cacao_Royaert2016", # TODO "Erreur : subscript contains out-of-bounds indices " for method 'gene'  WORKING; DONE
              "Zea_mays_MaizeGDBConsensus_v4", # WORKING; DONE
              "Zea_mays_IBM_MaizeSNP50", # WORKING; DONE  Error on numerals
              "Setaria_italica_Ni2017", # WORKING; DONE
              "Solanum_lycopersicum_Gonda2018", # WORKING; DONE
              "Triticum_aestivum_GutierrezGonzalez2019")

batch_gc = function(x) {
  source("sources/Genome.R") # Essential functions to parse GFF files and compute genomic statistics
  # GENES
  gc.map(x, ex = TRUE, method = 'gene')
  # EXONS
  gc.map(x, ex = TRUE, method = 'exon')
  # RECOMBINATION MAPS
  gc.map(x, ex = TRUE, method = 'map')
}
sapply(list_maps, batch_gc)


#============================================================================#
#   Nucleotidic diversity ----
#============================================================================#





#============================================================================#
# Promoters positions ----
#============================================================================#











#============================================================================#

#   Repeatitive elements & Transposons ----

#============================================================================#

#============================================================================#
# Repeat density ----
#============================================================================#
# A simple proxy of distribution of repeatitive elements in the genome
# The proportion of masked nucleotide sin windows of 100kb
# Masked nucleotides are written in lower case in fasta for NCBI/Ensemble genomes







#============================================================================#

#   Codon usage bias -----

#============================================================================#

# Codon distribution
# Degenerated codon may be biased in favor of some codon more easily translated - coevolution with ARNt


#============================================================================#
# End ----
#============================================================================#
