############################################################################
#                    ECOBIO - PhD
#
#       Linkage Map analysis of A. lyrata wild populations
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
library(qtl)
library(ASMap)

# install.packages('argparser')
# install.packages('yaml')
# install.packages('abind')
# install.packages('viridis')
# install.packages('xlsx')
# BiocManager::install("BiocLite")
# BiocManager::install(c("VariantAnnotation", "rhdf5"))
# BiocManager::install(c('GenomeInfoDb', 'GenomicRanges', 'rtracklayer', 'VariantAnnotation'))
# install.packages("remotes")
# remotes::install_github("gact/shmootl")
library(shmootl)
library(utils)
# devtools::install_github("gact/utl")
library(utl)
library(readODS)
library(Rsamtools)
library(onemap)
library(vcfR)

#============================================================================
# Loading variables & objects
#============================================================================

# Get the directory of the file & set working directory
wd=dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

# Custom functions
source("sources/PhysMap.R") # Get the physical positions of a marker given a list of marker names
source("sources/MareyMap.R") # Consolidate a MareyMap input file
source("sources/core_pipeline.R") # Basic functions but essential for the pipeline (e.g. saving logs)
source("sources/stats.R") # Statistical functions implemented
source("sources/convert_vcf_to_genfile.R")
# source("sources/readvcf2geno.R")

#============================================================================
# Import data
#============================================================================
# Example data
# data("mapBCu")

# Read crossings infos with pedigree
crossings = read.table("data/Genetic_maps/rawdata/AlyrataCrossingsParents/takou_transcriptomes_familyInfo.csv", header = TRUE, sep = "\t")
# Population sizes
nrow(crossings)
# 133 descendants
length(unique(crossings$SP.parent))
length(unique(crossings$PL.parent))
# 6 SP parents, 5 PL parents
length(unique(paste(crossings$Dam, "x", crossings$Sire, sep = " ")))
# 20 different type of reciprocal exchanges (Dam x Sire)

# We should have 133 + 11 = 144 samples


# Samples IDs of progeny
sampleID = getSamplesVCF(paste(wd, "/data/Genetic_maps/rawdata/AlyrataCrossingsParents/Alyrata_F1transc.minQ60.minDP10.maxalleles2.noindels.recode.vcf.gz", sep = ""))
parent1ID = NA
parent2ID = NA
for (i in 1:length(sampleID)) {
  parent1ID[i] = as.character(crossings$Dam[which(crossings$PlantID == sampleID[i])])
  parent2ID[i] = as.character(crossings$Sire[which(crossings$PlantID == sampleID[i])])
}
parent1ID
parent2ID

# # Convert vcf to raw onemap format
# vcf2raw(input = gzfile("data/Genetic_maps/rawdata/AlyrataCrossingsParents/Alyrata_F1transc.minQ60.minDP10.maxalleles2.noindels.recode.vcf.gz"),
#         output = "data/Genetic_maps/rawdata/AlyrataCrossingsParents/Onemap_input.raw",
#         cross = c("outcross"),
#         parent1 = parent1ID,
#         parent2 = parent2ID)

vcfR_progenies = read.vcfR("data/Genetic_maps/rawdata/AlyrataCrossingsParents/Alyrata_F1transc.minQ60.minDP10.maxalleles2.noindels.recode.vcf.gz")
vcfR_parents = read.vcfR("data/Genetic_maps/rawdata/AlyrataCrossingsParents/Alyrata.SPPL.NoRepeats.NoHet.minQ30.minDP10.miss80.minGQ30.maf05.recode.vcf.gz")
# Subset SNPs that are both in parents and progenies
# Markers with the same chromosome and position
mkr_progenies = paste(vcfR_progenies@fix[,1], "_", vcfR_progenies@fix[,2], sep = "")
mkr_parents = paste(vcfR_parents@fix[,1], "_", vcfR_parents@fix[,2], sep = "")
# Number of shared markers
sum(mkr_progenies %in% mkr_parents)
sum(mkr_parents %in% mkr_progenies)
# Not the same number of markers in each selection
# Some markers are counted twice in mkr_progenies?

# Making lists with the same number of markers
idx_progenies = which(mkr_progenies %in% mkr_parents)
idx_parents = rep(NA, length(idx_progenies))
# for (i in 1:length(idx_progenies)) {
#   cat(i/length(idx_progenies)*100, "\n")
#   idx_parents[i] = which(mkr_parents == mkr_progenies[idx_progenies[i]])[1]
# }
# For loop is really too slow
# x = idx_progenies[1]
# which(mkr_parents == mkr_progenies[x])[1]
# idx_progenies_tmp = idx_progenies[1:100]
library(pbapply)
# As there is multiple answers for some markers, select only the first correspondance,
# although it is arbitrary
idx_parents = pblapply(idx_progenies, function(x) which(mkr_parents == mkr_progenies[x])[1])
idx_parents = unlist(idx_parents)

# Subset only shared markers
vcfR_progenies_subset = vcfR_progenies[idx_progenies,]
vcfR_parents_subset = vcfR_parents[idx_parents,]

# The subset represents % of the original dataset

# Save subsets in .Rda
save(vcfR_progenies_subset, file = "data/Genetic_maps/rawdata/AlyrataCrossingsParents/vcfR_progenies_subset.Rda")
load("data/Genetic_maps/rawdata/AlyrataCrossingsParents/vcfR_progenies_subset.Rda")
save(vcfR_parents_subset, file = "data/Genetic_maps/rawdata/AlyrataCrossingsParents/vcfR_parents_subset.Rda")
load("data/Genetic_maps/rawdata/AlyrataCrossingsParents/vcfR_parents_subset.Rda")

# Combine parents and progenies in a single vcfR object
vcfR_file = vcfR_progenies_subset
# vcfR_file@fix = rbind(vcfR_progenies_subset@fix, vcfR_parents_subset@fix)
vcfR_file@gt = cbind(vcfR_progenies_subset@gt, vcfR_parents_subset@gt)
vcfR_file
# 173 samples
# 8 CHROMs
# 457,006 variants
save(vcfR_file, file = "data/Genetic_maps/rawdata/AlyrataCrossingsParents/onemap_input.Rda")
load("data/Genetic_maps/rawdata/AlyrataCrossingsParents/onemap_input.Rda")

# Convert to onemap input
onemap_input = onemap_read_vcfR(vcfR.object = vcfR_file,
                    cross = c("outcross"),
                    parent1 = parent1ID,
                    parent2 = parent2ID)










#-----------------------------------------
#   DEPRECATED
# Convert vcf to ABH genotypye format
# Tassel GenosToABHPlugin
# This plugin was designed to convert genotypes of a biparental population that came out of the TASSEL GBS pipeline in a nucleotide-based format to a parent-based format. It also removes sites in which the parent genotype is missing, ambiguous or heterozygous. The output of this plugin is a .csv file with data arrangement specific to the R package ‘qtl’ [1].
# The parameters to this plugin are:
# -o <Output file>: output filename, will be written as .csv
# -parentA <Parent A>: A plain text file containing the name of all samples from parent A as they are found in the input filename. One name per line.
# -parentB <Parent B>: A plain text file containing the name of all samples from parent B as they are found in the input filename. One name per line.
# -outputFormat <output format>: if "c", output will be A,H,B for parent A, het, and parent B: if "i", output will be 0,1,2 for parent A, het, parent B, if "r", output will be 0,0.5,1 for parent A, het, parent B. This field it not required. If absent, the default is "c" - output will be in the form of A,H,B.

# Generate parent A

# Generate parent B



# Read data files
# Importing QTL mapping data into R is accomplished with the read.cross function.






# ?readGenoVCF()
genotypes = readGenoVCF(gzfile("data/Genetic_maps/rawdata/AlyrataCrossingsParents/Alyrata_F1transc.minQ60.minDP10.maxalleles2.noindels.recode.vcf.gz"),
                        samples = unique(crossings$PlantID),
                        founder = unique(c(crossings$Dam, crossings$Sire)))
genotypes = readGenoVCF(paste(wd, "/data/Genetic_maps/rawdata/AlyrataCrossingsParents/Alyrata_F1transc.minQ60.minDP10.maxalleles2.noindels.recode.vcf", sep = ""),
                        samples = unique(crossings$PlantID),
                        founder = unique(c(crossings$Dam, crossings$Sire)))
getSamplesVCF(paste(wd, "/data/Genetic_maps/rawdata/AlyrataCrossingsParents/Alyrata_F1transc.minQ60.minDP10.maxalleles2.noindels.recode.vcf", sep = ""))

#' infiles <- c('samples.vcf', 'founders.vcf')
#' genfile <- 'geno.csv'
# convert_vcf_to_genfile(c("data/Genetic_maps/rawdata/AlyrataCrossingsParents/Alyrata_F1transc.minQ60.minDP10.maxalleles2.noindels.recode.vcf",
#                          "data/Genetic_maps/rawdata/AlyrataCrossingsParents/Alyrata.SPPL.NoRepeats.NoHet.minQ30.minDP10.miss80.minGQ30.maf05.recode.vcf"),
#                        genfile = "data/Genetic_maps/rawdata/AlyrataCrossingsParents/geno.csv",
#                        samples = read_samples_from_vcf('data/Genetic_maps/rawdata/AlyrataCrossingsParents/Alyrata_F1transc.minQ60.minDP10.maxalleles2.noindels.recode.vcf'),
#                        founders = read_samples_from_vcf('data/Genetic_maps/rawdata/AlyrataCrossingsParents/Alyrata.SPPL.NoRepeats.NoHet.minQ30.minDP10.miss80.minGQ30.maf05.recode.vcf'))
# 
# convert_vcf_to_genfile(c("data/Genetic_maps/rawdata/AlyrataCrossingsParents/Alyrata_F1transc.minQ60.minDP10.maxalleles2.noindels.recode.vcf",
#                          "data/Genetic_maps/rawdata/AlyrataCrossingsParents/Alyrata.SPPL.NoRepeats.NoHet.minQ30.minDP10.miss80.minGQ30.maf05.recode.vcf"),
#                        genfile = "data/Genetic_maps/rawdata/AlyrataCrossingsParents/geno.csv",
#                        samples = c("1964", "1960", "1760"),
#                        founders = c("70535", "80936"))
genotypes = read.table("data/Genetic_maps/rawdata/AlyrataCrossingsParents/geno.csv", sep = ",")
genotypes = read.table("data/Genetic_maps/rawdata/AlyrataCrossingsParents/listeria.csv", header = TRUE, sep = ",")
genotypes = read.table("data/Genetic_maps/rawdata/AlyrataCrossingsParents/test_hapmap.csv", sep = "\t")


genotypes = read.cross(format = "csvs", genfile="data/Genetic_maps/rawdata/AlyrataCrossingsParents/geno.csv",
                       phefile="data/Genetic_maps/rawdata/AlyrataCrossingsParents/pheno.csv",
                       alleles = c("A", "B"),
                       estimate.map = FALSE, genotypes = c("A","H","B","D","C"), na.strings = "-")

genotypes = read.cross(format = "csv", file = "data/Genetic_maps/rawdata/AlyrataCrossingsParents/test_hapmap.csv")


cross = readCrossHDF5("data/Genetic_maps/rawdata/AlyrataCrossingsParents/test")




setwd(paste(wd, "/data/Genetic_maps/rawdata/AlyrataCrossingsParents", sep = ""))
genotypes = readGenoVCF(c("Alyrata_F1transc.minQ60.minDP10.maxalleles2.noindels.recode.vcf",
                          "Alyrata.SPPL.NoRepeats.NoHet.minQ30.minDP10.miss80.minGQ30.maf05.recode.vcf"),
                        samples = getSamplesVCF("Alyrata_F1transc.minQ60.minDP10.maxalleles2.noindels.recode.vcf"))
genotypes = readGenoVCF(c("Alyrata_F1transc.minQ60.minDP10.maxalleles2.noindels.recode.vcf",
                          "Alyrata.SPPL.NoRepeats.NoHet.minQ30.minDP10.miss80.minGQ30.maf05.recode.vcf"),
                        samples = getSamplesVCF("Alyrata_F1transc.minQ60.minDP10.maxalleles2.noindels.recode.vcf"))
# genotypes = readGenoVCF(c("data/Genetic_maps/rawdata/AlyrataCrossingsParents/Alyrata_F1transc.minQ60.minDP10.maxalleles2.noindels.recode.vcf",
#                           "data/Genetic_maps/rawdata/AlyrataCrossingsParents/Alyrata.SPPL.NoRepeats.NoHet.minQ30.minDP10.miss80.minGQ30.maf05.recode.vcf"),
#                         samples = getSamplesVCF("data/Genetic_maps/rawdata/AlyrataCrossingsParents/Alyrata_F1transc.minQ60.minDP10.maxalleles2.noindels.recode.vcf"),
#                         founder = getSamplesVCF("data/Genetic_maps/rawdata/AlyrataCrossingsParents/Alyrata.SPPL.NoRepeats.NoHet.minQ30.minDP10.miss80.minGQ30.maf05.recode.vcf"))
setwd(wd)


# readGenoVCF("./Merged_RIAILs_SNPs/RIAILs_SNPs_filtered_position.vcf.recode.vcf", samples = rils, founders = mapping(c(APS7 = 'A', APS14 = 'B')))

# Convert to ASMap dataset


#============================================================================
#
#============================================================================















#============================================================================
# Linkage map construction @davikDdRADBasedLinkage2015
# Data produced by Stacks were coded as an F1 segregating population using the genotypes of
# the parental lines to assign segregation. Data were filtered for all markers containing more than
# 50% missing values and a chi-squared analysis was performed to determine segregation distortion
# at the 5% level of significance (Chi-squared = 3.841 (1 d.f.), 5.991 (2 d.f.) or 7.815 (3 d.f.)).
# Initially, only robust markers for which no significant segregation distortion was observed at
# the 5% level were used for linkage mapping using JoinMap 4.1 (Kyazma, NL). An initial linkage
# map was constructed using the Maximum Likelihood mapping function and assessed for spurious
# linkages or inflated genetic distances, with individual genotypes being converted to missing
# values where necessary or loci being removed completely where they caused conflicts in the
# data. Following scrutiny of the Maximum Likelihood data, linkage mapping was performed
# with the initial marker set using regression mapping and v1.0 linkage groups were produced.
# Additional data were then added to the dataset for markers exhibiting significant segregation
# distortion. Marker placement was determined using regression mapping with a minimum logarithm
# of odds (LOD) score threshold of 3.0, a recombination fraction threshold of 0.35, ripple
# value of 1.0, jump threshold of 3.0 and a triplet threshold of 5.0, and mapping distances were


# @konarHighqualityGeneticMapping2017
# JoinMap for generating the map
# MapChart for visualization
