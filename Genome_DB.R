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

# Within this script you have functions to build the genome database and gather some statistics and metadata from the NCBI API

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

#==========================================================
# Loading variables & objects
#==========================================================

# Get the directory of the file & set working directory
filedir=dirname(rstudioapi::getSourceEditorContext()$path)
wd=filedir
setwd(wd)
source("sources/GenomeDB.R")

#==========================================================
# Import & compile bibliographic database (DB)
#==========================================================
# Compile the DB into a clean DB .txt file in 'data'
genomeDB.compile(data = "data/Genome/Genome_ressources.xlsx", output = "data/Genome/GenomeDB.txt")

#==========================================================
# Summary statistics
#==========================================================
# read the database
genDB = read.table("data/Genome/GenomeDB.txt", header = TRUE, sep = "\t")
genomeDB.summary(data = "data/Genome/GenomeDB.txt")



#==========================================================
# Download genomes
#==========================================================
# make a list of all genomes you want to download
list_species = unique(genDB$species)
# First, select the latest genome assembly in NCBI/Ensembl for each species
index = data.frame(species = rep(NA, length(list_species)),
                   accession = rep(NA, length(list_species)),
                   database = rep(NA, length(list_species)))  # A vector of accessions to keep in the database
genDB$accession = as.character(genDB$accession)
genDB$source = as.character(genDB$source)
# Retrieve all accessions for each species
for (i in 1:length(list_species)) {
  print(i)
  index$species[i] = as.character(list_species[i])
  print(list_species[i])
  # Retrieve accessions
  acc = genDB$accession[which(as.character(genDB$species) == list_species[i])]
  print(acc)
  # If more than one accession,
  # Which accession is the latest?
  t = as.Date(genDB$release_date[which(as.character(genDB$species) == list_species[i])])
  print(t)
  if(length(t) > 1) {
    # Which accession is the latest?
    if (length(acc[which(t == max(t, na.rm = TRUE))]) > 1) { # If different accession numbers
      # Keep the one with the motif 'GCA' or 'GCF'
      index$accession[i] = max(grep("GC", acc[which(t == max(t, na.rm = TRUE))], value = TRUE))
      index$database[i] = genDB$source[which(genDB$accession == max(grep("GC", acc[which(t == max(t, na.rm = TRUE))], value = TRUE)))]
    } else {
      index$accession[i] = acc[which(t == max(t, na.rm = TRUE))]
      index$database[i] = genDB$source[which(genDB$accession == acc[which(t == max(t, na.rm = TRUE))])]
    }
    
  } else { # If no accession for this species
    if(length(t) == 0) {
      index$accession[i] = NA
      index$database[i] = NA
    } else {  # Otherwise, get accession for the species (t == 1)
      index$accession[i] = acc
      index$database[i] = genDB$source[which(as.Date(genDB$release_date) == t)]
    }
  }
}
# remove NA in accessions or database
index = index[!is.na(index$accession),]
index = index[!is.na(index$database),]
# Save the list of accessions (i.e. genomes) to download...
write.table(index, file = "data/Genome/list_accessions.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Now you can download genomes with the script 'genome_dl.sh' to construct your own database

#==========================================================
# END
#==========================================================