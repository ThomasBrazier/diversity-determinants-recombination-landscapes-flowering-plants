########################################################################### #
#                    ECOBIO - PhD
#
#           Genetic Shuffling
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
library(ggplot2)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
library(ggplotify)
library(pals)
library(cowplot)
library(lme4)
library(lmerTest)
library(car)

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
source("sources/geneticshuffling.R")

############################################################################ #
#============================================================================#
# ESTIMATE THE INTRA-CHROMOSOMAL GENETIC SHUFFLING ----
#============================================================================#
############################################################################ #

#============================================================================#
# Interpolate 1,000 evenly spaced markers with genetic distances ----
#============================================================================#
# For an evenly sampling along the chromosome
# Use the recombination-Map_Veller() function to interpolate a Marey map of 1,000 evenly distributed markers (according to physical positions)
# Interpolated maps of 1,000 loci are stored in 'data-cleaned/veller/'
# No bootstrap
source("sources/MareyMap.R") 
source("sources/geneticshuffling.R") 

# For a chosen dataset
set = "Arabidopsis_thaliana_Serin2017"
set = "Panicum_hallii_Lovell2018"
recombinationMap_Veller(set = set, chr = "all", K = 3)

# Or all datasets at once
# List of all maps
list = system(paste("ls ", wd, "/data-cleaned/marey_maps/", sep = ""), intern = TRUE)
list = list[grep("*.txt", list)]
list = gsub(".txt", "", list)
list = list[!(list == "AllMaps")]
list
for (set in list) {
  cat("===============================================\n")
  cat("Processing ", set, "...\n")
  recombinationMap_Veller(set = set, chr = "all", K = 3)
}

#============================================================================#
# Genetic shuffling rate of chromosomes ----
#============================================================================#
# Estimate the intra-chromosomal part of the genetic shuffling rate, according to Veller et al. 2019
# Open the dataset
df.veller = read.table("output/veller/GeneticDistances.txt", header = TRUE)

# Reversed mapping function: from the genetic distances to the recombination fraction
# Reversed Kosambi from the formula given in Veller et al. 2019
# r_ij = r(d_ij) = 1/2(tanh(2d_ij))
source("sources/geneticshuffling.R") 

# For a single chromosome
set = "Gossypium_raimondii_Wang2013"
chromosome = "13"
R_chr(df.veller, set = set, chromosome = chromosome)
R_chr_fromMatlab(df.veller, set = set, chromosome = chromosome)


# Compute the genetic shuffling rate per chromosome
list_chr = unique(paste(df.veller$set, df.veller$map, sep = ":::"))
df_geneticshuffling_physical_chr = data.frame(set = gsub(":::[A-Za-z0-9]*$", "", list_chr),
                                              chromosome = gsub("[_A-Za-z0-9]*:::", "", list_chr),
                                              geneticshuffling = NA)
for (i in 1:nrow(df_geneticshuffling_physical_chr)) {
  cat("==========================================\n")
  df_geneticshuffling_physical_chr$geneticshuffling[i] = R_chr(df.veller, set = as.character(df_geneticshuffling_physical_chr$set[i]),
                                                               chromosome = as.character(df_geneticshuffling_physical_chr$chromosome[i]))
}

# Save in tables
write.table(df_geneticshuffling_physical_chr, file = "tables/df_geneticshuffling_physical_chr.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
df_geneticshuffling_physical_chr = read.table("tables/df_geneticshuffling_physical_chr.txt", header = TRUE)


# # Compute the mean genetic shuffling rate per species
# list_set = unique(df.veller$set)
# df_geneticshuffling_physical = data.frame(set = list_set, geneticshuffling = NA)
# for (i in 1:nrow(df_geneticshuffling_physical)) {
#   print(i)
#   df_geneticshuffling_physical$geneticshuffling[i] = R_chr_mean(df.veller, df_geneticshuffling_physical$set[i])
# }
# 
# # Save in tables
# write.table(df_geneticshuffling_physical, file = "tables/df_geneticshuffling_physical.csv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
# df_geneticshuffling_physical = read.table("tables/df_geneticshuffling_physical.csv", header = TRUE)

# Alternative method to estimate r_bar
# Compute the genetic shuffling rate per chromosome
source("sources/geneticshuffling.R") 
list_chr = unique(paste(df.veller$set, df.veller$map, sep = ":::"))
df_geneticshuffling_physical_chr_matlab = data.frame(set = gsub(":::[A-Za-z0-9]*$", "", list_chr),
                                              chromosome = gsub("[_A-Za-z0-9]*:::", "", list_chr),
                                              geneticshuffling = NA)
for (i in 1:nrow(df_geneticshuffling_physical_chr_matlab)) {
  cat("==========================================\n")
  df_geneticshuffling_physical_chr_matlab$geneticshuffling[i] = R_chr_fromMatlab(df.veller, set = as.character(df_geneticshuffling_physical_chr_matlab$set[i]),
                                                               chromosome = as.character(df_geneticshuffling_physical_chr_matlab$chromosome[i]))
}
# Faster implementation
# Same results

# Assess convergence of the methods
plot(df_geneticshuffling_physical_chr$geneticshuffling, df_geneticshuffling_physical_chr_matlab$geneticshuffling,
     xlab = "Method 1 (equation)", ylab = "Method 2 (Matlab)")

#============================================================================#
# Exploratory analyses ----
#============================================================================#
df_geneticshuffling_physical_chr = read.table("tables/df_geneticshuffling_physical_chr.txt", header = TRUE, stringsAsFactors = FALSE)
metadata.clean = read.table("data-cleaned/Genetic_maps_ressources.csv", header = TRUE, sep = ";", stringsAsFactors = FALSE)
chromosome.stats = read.table("tables/chromosome.stats.csv", header = TRUE, sep = ";", stringsAsFactors = FALSE)
chromosome.stats = merge(chromosome.stats, df_geneticshuffling_physical_chr, by = c("set", "chromosome"))

#============================================================================#
# At a chromosome level ----
#============================================================================#
# "A CO toward the middle of a chromosome causes more shuffling than a CO toward the tips." veller et al. 2019
# (1) Hypothesis: Genetic shuffling is larger when COs ar more evenly spaced along the chromosome.
#       We expect a significant negative relationship with the periphery-bias ratio
#       and a positive relationship with chromosome size/recombination as an indirect proxy of the u-shaped landscapes

# (2) Hypothesis: When thinking in terms of gene distance instead of basepair distance, genetic shuffling should be more equally distributed among chromosomes

hist(chromosome.stats$geneticshuffling, breaks = 40, xlab = "Genetic shuffling")

#--------------------------------------------------------------------------#
# Is there a relationship between chromosome size and the genetic shuffling rate?
# Larger chromosomes have recombination concentrated in the tips
# Hence chromosome size should be an indirect proxy of the evenness of the distribution of the recombination
# We expect larger chromosomes to have lower genetic shuffling rate
plot(chromosome.stats$phys.map.length/1000000, chromosome.stats$geneticshuffling,
     log = "x", xlab = "Chromosome size", ylab = "Genetic shuffling")
# None clear relationship can be seen
# Dig furthermore with models...
mod = lm(geneticshuffling ~ log10(phys.map.length/1000000), data = chromosome.stats)
summary(mod)
AIC(mod)

# Linear Mixed Model with species effect
lmer.model = lmer(geneticshuffling ~ log10(phys.map.length/1000000) + (1|species), data = chromosome.stats)
lmer.model = lmer(geneticshuffling ~ log10(phys.map.length/1000000) + (log10(phys.map.length/1000000)|species), data = chromosome.stats)
summary(lmer.model)
plot(lmer.model)
AIC(lmer.model)

# Diagnostic plots
par(mfrow = c(1, 3))
# Residuals distribution
qqPlot(residuals(lmer.model), ylab = "Residuals", id = F)
# Residuals as a function of fitted
plot(fitted(lmer.model), residuals(lmer.model), xlab = "fitted", ylab = "residuals")
abline(h = 0)
# Observed~predicted qqplot
plot(log10(chromosome.stats$phys.map.length/1000000), predict(lmer.model), xlab = "observed", ylab = "predicted")
abline(a = 0, b = 1)
par(mfrow = c(1,1))

plot(chromosome.stats$phys.map.length/1000000, chromosome.stats$geneticshuffling,
     xlab = "Chromosome size",
     ylab = "Genetic shuffling",
     main = "LM (red) & LMM (blue)",
     log = "x")
abline(mod,
       untf = FALSE,
       col = "red")
abline(fixef(lmer.model)[1],
       fixef(lmer.model)[2],
       untf = FALSE,
       col = "blue")

# We cannot conclude on these results


#--------------------------------------------------------------------------#
# Now is there a relationship between recombination rates and the genetic shuffling rate?
plot(chromosome.stats$mean.recrate, chromosome.stats$geneticshuffling,
     log = "x", xlab = "Recombination rate", ylab = "Genetic shuffling")
# Interestingly, we can see a relationship, though there is a lot of variance
# Dig furthermore with models...
mod = lm(geneticshuffling ~ log10(mean.recrate), data = chromosome.stats)
summary(mod)
AIC(mod)

# Linear Mixed Model with species effect
lmer.model = lmer(geneticshuffling ~ log10(mean.recrate) + (1|species), data = chromosome.stats)
lmer.model = lmer(geneticshuffling ~ log10(mean.recrate) + (log10(mean.recrate)|species), data = chromosome.stats)
summary(lmer.model)
plot(lmer.model)
AIC(lmer.model)

# Diagnostic plots
par(mfrow = c(1, 3))
# Residuals distribution
qqPlot(residuals(lmer.model), ylab = "Residuals", id = F)
# Residuals as a function of fitted
plot(fitted(lmer.model), residuals(lmer.model), xlab = "fitted", ylab = "residuals")
abline(h = 0)
# Observed~predicted qqplot
plot(chromosome.stats$geneticshuffling[which(!is.na(chromosome.stats$mean.recrate))], predict(lmer.model), xlab = "observed", ylab = "predicted")
abline(a = 0, b = 1)
par(mfrow = c(1,1))

plot(chromosome.stats$mean.recrate, chromosome.stats$geneticshuffling,
     xlab = "Recombination rate",
     ylab = "Genetic shuffling",
     main = "LM (red) & LMM (blue)",
     log = "x")
abline(mod,
       untf = FALSE,
       col = "red")
abline(fixef(lmer.model)[1],
       fixef(lmer.model)[2],
       untf = FALSE,
       col = "blue")

#--------------------------------------------------------------------------#
# Is there a relationship between the periphery-bias ratio and the genetic shuffling rate?
# According to Veller et al. 2019, chromosomes that tend to concentrate recombination in their distal regions
# have a lower chromosomal genetic shuffling rate
# The higher is the bias toward the periphery, the lower is the genetic shuffling, according to (1)
hist(chromosome.stats$peripherybias_ratio, breaks = 40, xlab = "Periphery-bias ratio")
hist(sqrt(chromosome.stats$peripherybias_ratio), breaks = 40, xlab = "Periphery-bias ratio (sqrt)")

plot(sqrt(chromosome.stats$peripherybias_ratio), chromosome.stats$geneticshuffling,
     xlab = "Periphery-bias ratio", ylab = "Genetic shuffling")
# The result is not the expected. Chromosomes with a lower bias toward periphery

# Dig furthermore with models...
mod = lm(geneticshuffling ~ sqrt(peripherybias_ratio), data = chromosome.stats)
summary(mod)
AIC(mod)

# Linear Mixed Model with species effect
lmer.model = lmer(geneticshuffling ~ sqrt(peripherybias_ratio) + (1|species), data = chromosome.stats)
lmer.model = lmer(geneticshuffling ~ sqrt(peripherybias_ratio) + (sqrt(peripherybias_ratio)|species), data = chromosome.stats)
summary(lmer.model)
plot(lmer.model)
AIC(lmer.model)

# Diagnostic plots
par(mfrow = c(1, 3))
# Residuals distribution
qqPlot(residuals(lmer.model), ylab = "Residuals", id = F)
# Residuals as a function of fitted
plot(fitted(lmer.model), residuals(lmer.model), xlab = "fitted", ylab = "residuals")
abline(h = 0)
# Observed~predicted qqplot
plot(chromosome.stats$geneticshuffling[which(!is.na(chromosome.stats$peripherybias_ratio))], predict(lmer.model), xlab = "observed", ylab = "predicted")
abline(a = 0, b = 1)
par(mfrow = c(1,1))

plot(sqrt(chromosome.stats$peripherybias_ratio), chromosome.stats$geneticshuffling,
     xlab = "Periphery-bias ratio", ylab = "Genetic shuffling")
abline(mod,
       col = "red")
abline(fixef(lmer.model)[1],
       fixef(lmer.model)[2],
       col = "blue")


# Finally, we have to consider the impact of the genetic map length on the genetic shuffling (i.e. the number of COs)
# Larger genetic map should have a higher genetic shuffling, though this increases sublinearly according to Veller et al. 2019
#--------------------------------------------------------------------------#
# Is there a relationship between the genetic map length and the genetic shuffling rate?
hist(chromosome.stats$linkage.map.length.correctedHW, breaks = 40, xlab = "Genetic map length")

plot(chromosome.stats$linkage.map.length.correctedHW, chromosome.stats$geneticshuffling,
     xlab = "Genetic map length", ylab = "Genetic shuffling")
# Indeed, larger maps, with more COs, have a higher genetic shuffling rate.
# Dig furthermore with models...
mod = lm(geneticshuffling ~ linkage.map.length.correctedHW, data = chromosome.stats)
summary(mod)
AIC(mod)

# Linear Mixed Model with species effect
lmer.model = lmer(geneticshuffling ~ linkage.map.length.correctedHW + (1|species), data = chromosome.stats)
lmer.model = lmer(geneticshuffling ~ linkage.map.length.correctedHW + (linkage.map.length.correctedHW|species), data = chromosome.stats)
summary(lmer.model)
plot(lmer.model)
AIC(lmer.model)

# Diagnostic plots
par(mfrow = c(1, 3))
# Residuals distribution
qqPlot(residuals(lmer.model), ylab = "Residuals", id = F)
# Residuals as a function of fitted
plot(fitted(lmer.model), residuals(lmer.model), xlab = "fitted", ylab = "residuals")
abline(h = 0)
# Observed~predicted qqplot
plot(chromosome.stats$geneticshuffling, predict(lmer.model), xlab = "observed", ylab = "predicted")
abline(a = 0, b = 1)
par(mfrow = c(1,1))

plot(chromosome.stats$linkage.map.length.correctedHW, chromosome.stats$geneticshuffling,
     xlab = "Genetic map length", ylab = "Genetic shuffling")
abline(mod,
       col = "red")
abline(fixef(lmer.model)[1],
       fixef(lmer.model)[2],
       col = "blue")

# It seems trivial but the more COs they have, the more chromosomes shuffles their alleles at random
# Hence the obvious correlation between recombination rate and genetic shuffling

# Yet the periphery-bias ratio is unsensitive to chromosome size, number of COs (linkage map length) or recombination rate inter-specific variations
# And we still find a slight negative relationship between the bias toward the periphery and the genetic shuffling.

# Biologically it is interesting to note that larger chromosomes have lower recombination rates and a clear bias toward the tips
# hence they should be less efficient than smaller chromosomes in terms of genetic shuffling rate

# BUT we have seen that measuring distances in basepairs may not be the most appropriate measure:
# DNA sequences are linear, but the genome is not physically strictly linear
# DNA can be condensed or not (heterochromatin/euchromatin), and not all regions are qually accessible for recombination
# Besides, recombination seems to be preferentially initated in the begining of genes...

# A new measure of genomic distances could be the gene distances (cumulative number of genes along the genome)
# A functional measure, taking into account variation in gene density and other genomic features that make the genome heterogeneous




############################################################################ #
#============================================================================#
# GENE DISTANCES ----
#============================================================================#
############################################################################ #

#============================================================================#
# Genetic distance (cM) expressed as a function of gene distance (cumulative number of genes) ----
#============================================================================#
source("sources/Genome.R")

# Test the method on A. thaliana
# Retrieve gene distances (cumulative number of genes) in a given Marey map
# gene_distance.map(set = map, scale = 0.1)

#-----------------------------------------------#
# Estimate
#-----------------------------------------------#
# Reproduce to all species in dataset
# List of all maps (all maps with a gc_genes file giving the list of genes and their position in a reference genome)
list = system(paste("ls ", wd, "/data-cleaned/genome/gene_positions/", sep = ""), intern = TRUE)
list = list[grep("*.txt.gz", list)]
list = gsub("_gene_positions.txt.gz", "", list)
list

# BUGFIX 2021/11/19
sp = "Manihot_esculenta"
sp = "Theobroma_cacao"
sp = "Citrus_sinensis"
sp = "Arachis_duranensis"
source("sources/Genome.R")
gene_distance.map(species = sp, scale = 0.1)

# For all species at once
source("sources/Genome.R") 
for (sp in list) {
  gene_distance.map(species = sp, scale = 0.1)
}

#-----------------------------------------------#
# Visualize
#-----------------------------------------------#
# Load the map with gene distances and genetic distances (cM)
list
set = "Arabidopsis_thaliana_Serin2017"
chromosome = "1"
gene_distance_map = read.table(file = gzfile(paste("data-cleaned/genome/gene_distance/", set, ".txt.gz", sep = "")), header = TRUE, sep = "\t")
gene_distance_map = gene_distance_map[gene_distance_map$map == chromosome,]
# Figure
par(mfrow = c(1,2))
plot(gene_distance_map$phys/1000000, gene_distance_map$gen, xlab = "Physical distance (Mb)", ylab = "Genetic distance (cM)")
lines(x = c(min(gene_distance_map$phys/1000000), max(gene_distance_map$phys/1000000)), y = c(min(gene_distance_map$gen), max(gene_distance_map$gen)))
plot(gene_distance_map$gene_distance, gene_distance_map$gen, xlab = "Gene distance\n(cumulative number of genes)", ylab = "Genetic distance (cM)")
lines(x = c(min(gene_distance_map$gene_distance), max(gene_distance_map$gene_distance)), y = c(min(gene_distance_map$gen), max(gene_distance_map$gen)))
par(mfrow = c(1,1))
# Interestingly,
plot(gene_distance_map$phys/1000000, gene_distance_map$gene_distance, xlab = "Physical distance (Mb)", ylab = "Genetic distance (cM)")


set = "Zea_mays_MaizeGDBConsensus_v4"
chromosome = "1"
gene_distance_map = read.table(file = gzfile(paste("data-cleaned/genome/gene_distance/", set, ".txt.gz", sep = "")), header = TRUE, sep = "\t")
gene_distance_map = gene_distance_map[gene_distance_map$map == chromosome,]
# Figure
par(mfrow = c(1,2))
plot(gene_distance_map$phys/1000000, gene_distance_map$gen, xlab = "Physical distance (Mb)", ylab = "Genetic distance (cM)")
lines(x = c(min(gene_distance_map$phys/1000000), max(gene_distance_map$phys/1000000)), y = c(min(gene_distance_map$gen), max(gene_distance_map$gen)))
plot(gene_distance_map$gene_distance, gene_distance_map$gen, xlab = "Gene distance\n(cumulative number of genes)", ylab = "Genetic distance (cM)")
lines(x = c(min(gene_distance_map$gene_distance), max(gene_distance_map$gene_distance)), y = c(min(gene_distance_map$gen), max(gene_distance_map$gen)))
par(mfrow = c(1,1))
# Interestingly,
plot(gene_distance_map$phys/1000000, gene_distance_map$gene_distance, xlab = "Physical distance (Mb)", ylab = "Genetic distance (cM)")




############################################################################ #
#============================================================================#
# GENE DISTANCES - ESTIMATE THE INTRA-CHROMOSOMAL GENETIC SHUFFLING ----
#============================================================================#
############################################################################ #

#============================================================================#
# Interpolate 1,000 evenly spaced markers with genes distances ----
#============================================================================#
# Rethink the genome from the point of view of genes distribution

# For an evenly sampling along the chromosome
# Use the recombination.map() function with method = "veller" to interpolate a Marey map of 1,000 evenly distributed markers (according to physical positions)
# Interpolated maps of 1,000 loci are stored in 'data-cleaned/veller/'
# No bootstrap
source("sources/MareyMap.R") 

# For a chosen dataset
# set = "Lupinus_albus_Ksiazkiewicz2017"
# GeneDistancesShuffling(set = set, chr = "all", K = 5)

# Or all datasets at once
# List of all maps
list = system(paste("ls ", wd, "/data-cleaned/genome/gene_distance/", sep = ""), intern = TRUE)
list = list[grep("*.txt.gz", list)]
list = gsub(".txt.gz", "", list)
list = list[!(list == "AllMaps")]
list
# Citrus_sinensis_Huang2018 has too few markers to efficiently interpolate with gene distances
list = list[-which(list == "Citrus_sinensis_Huang2018")]
for (set in list) {
  GeneDistancesShuffling(set = set, chr = "all", K = 3)
}

# Open the dataset
df.veller = read.table("output/veller/GeneDistances.txt", header = TRUE, stringsAsFactors = FALSE)
unique(df.veller$set)

# Compute the genetic shuffling rate per chromosome
list_chr = unique(cbind(df.veller$set, df.veller$map))
df_geneticshuffling_genedist_chr = data.frame(set = list_chr[,1],
                                              chromosome = list_chr[,2],
                                              geneticshuffling = NA)
for (i in 1:nrow(df_geneticshuffling_genedist_chr)) {
  cat("===================================================\n")
  df_geneticshuffling_genedist_chr$geneticshuffling[i] = R_chr(df.veller, set = as.character(df_geneticshuffling_genedist_chr$set[i]),
                                                               chromosome = as.character(df_geneticshuffling_genedist_chr$chromosome[i]))
}

# Save in tables
write.table(df_geneticshuffling_genedist_chr, file = "tables/df_geneticshuffling_genedist_chr.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
df_geneticshuffling_genedist_chr = read.table("tables/df_geneticshuffling_genedist_chr.txt", header = TRUE, stringsAsFactors = FALSE)





#============================================================================#
# Exploratory analyses ----
#============================================================================#
df_geneticshuffling_physical_chr = read.table("tables/df_geneticshuffling_physical_chr.txt", header = TRUE, stringsAsFactors = FALSE)
df_geneticshuffling_genedist_chr = read.table("tables/df_geneticshuffling_genedist_chr.txt", header = TRUE, stringsAsFactors = FALSE)
metadata.clean = read.table("data-cleaned/Genetic_maps_ressources.csv", header = TRUE, sep = ";", stringsAsFactors = FALSE)
chromosome.stats = read.table("tables/chromosome.stats.csv", header = TRUE, sep = ";", stringsAsFactors = FALSE)
chromosome.stats = merge(chromosome.stats, df_geneticshuffling_physical_chr, by = c("set", "chromosome"))
df_geneticshuffling_genedist_chr = read.table("tables/df_geneticshuffling_genedist_chr.txt", header = TRUE, stringsAsFactors = FALSE)
colnames(df_geneticshuffling_genedist_chr)[3] = "geneticshuffling_genedist"
chromosome.stats = merge(chromosome.stats, df_geneticshuffling_genedist_chr, by = c("set", "chromosome"))

#============================================================================#
# At a chromosome level ----
#============================================================================#
# "A CO toward the middle of a chromosome causes more shuffling than a CO toward the tips." veller et al. 2019
# (2) Hypothesis: When thinking in terms of gene distance instead of basepair distance, genetic shuffling should be more equally distributed among chromosomes

par(mfrow = c(1,2))
hist(chromosome.stats$geneticshuffling, breaks = 40, xlab = "Genetic shuffling", xlim = c(0,1))
hist(chromosome.stats$geneticshuffling_genedist, breaks = 40, xlab = "Genetic shuffling in gene distances", xlim = c(0,1))
par(mfrow = c(1,1))

plot(chromosome.stats$geneticshuffling, chromosome.stats$geneticshuffling_genedist,
     xlab = "Genetic shuffling", ylab = "Genetic shuffling in gene distances")
abline(0,1)
cor.test(chromosome.stats$geneticshuffling, chromosome.stats$geneticshuffling_genedist, method = "spearman")

#--------------------------------------------------------------------------#
# Is there a relationship between chromosome size and the genetic shuffling rate?
# Larger chromosomes have recombination concentrated in the tips
# Hence chromosome size should be an indirect proxy of the evenness of the distribution of the recombination
# We expect larger chromosomes to have lower genetic shuffling rate
plot(chromosome.stats$phys.map.length/1000000, chromosome.stats$geneticshuffling_genedist,
     log = "x", xlab = "Chromosome size", ylab = "Genetic shuffling in gene distances")
# None clear relationship can be seen
# Dig furthermore with models...
mod = lm(geneticshuffling_genedist ~ log10(phys.map.length/1000000), data = chromosome.stats)
summary(mod)
AIC(mod)

# Linear Mixed Model with species effect
lmer.model = lmer(geneticshuffling_genedist ~ log10(phys.map.length/1000000) + (1|species), data = chromosome.stats)
lmer.model = lmer(geneticshuffling_genedist ~ log10(phys.map.length/1000000) + (log10(phys.map.length/1000000)|species), data = chromosome.stats)
summary(lmer.model)
plot(lmer.model)
AIC(lmer.model)

# Diagnostic plots
par(mfrow = c(1, 3))
# Residuals distribution
qqPlot(residuals(lmer.model), ylab = "Residuals", id = F)
# Residuals as a function of fitted
plot(fitted(lmer.model), residuals(lmer.model), xlab = "fitted", ylab = "residuals")
abline(h = 0)
# Observed~predicted qqplot
plot(chromosome.stats$geneticshuffling_genedist, predict(lmer.model), xlab = "observed", ylab = "predicted")
abline(a = 0, b = 1)
par(mfrow = c(1,1))

plot(chromosome.stats$phys.map.length/1000000, chromosome.stats$geneticshuffling_genedist,
     xlab = "Chromosome size",
     ylab = "Genetic shuffling in gene distances",
     main = "LM (red) & LMM (blue)",
     log = "x")
abline(mod,
       untf = FALSE,
       col = "red")
abline(fixef(lmer.model)[1],
       fixef(lmer.model)[2],
       untf = FALSE,
       col = "blue")

# We cannot conclude on these results

#--------------------------------------------------------------------------#
# Is there a relationship between recombination rates and the genetic shuffling rate?
plot(chromosome.stats$mean.recrate, chromosome.stats$geneticshuffling_genedist,
     log = "x", xlab = "Recombination rate", ylab = "Genetic shuffling in gene distances")
# None clear relationship can be seen
# Dig furthermore with models...
mod = lm(geneticshuffling_genedist ~ log10(mean.recrate), data = chromosome.stats)
summary(mod)
AIC(mod)

# Linear Mixed Model with species effect
lmer.model = lmer(geneticshuffling_genedist ~ log10(mean.recrate) + (1|species), data = chromosome.stats)
lmer.model = lmer(geneticshuffling_genedist ~ log10(mean.recrate) + (log10(mean.recrate)|species), data = chromosome.stats)
summary(lmer.model)
plot(lmer.model)
AIC(lmer.model)

# Diagnostic plots
par(mfrow = c(1, 3))
# Residuals distribution
qqPlot(residuals(lmer.model), ylab = "Residuals", id = F)
# Residuals as a function of fitted
plot(fitted(lmer.model), residuals(lmer.model), xlab = "fitted", ylab = "residuals")
abline(h = 0)
# Observed~predicted qqplot
plot(chromosome.stats$geneticshuffling_genedist, predict(lmer.model), xlab = "observed", ylab = "predicted")
abline(a = 0, b = 1)
par(mfrow = c(1,1))

plot(chromosome.stats$mean.recrate, chromosome.stats$geneticshuffling_genedist,
     xlab = "Recombination rate",
     ylab = "Genetic shuffling in gene distances",
     main = "LM (red) & LMM (blue)",
     log = "x")
abline(mod,
       untf = FALSE,
       col = "red")
abline(fixef(lmer.model)[1],
       fixef(lmer.model)[2],
       untf = FALSE,
       col = "blue")


#--------------------------------------------------------------------------#
# Is there a relationship between the periphery-bias ratio and the genetic shuffling rate?
# According to Veller et al. 2019, chromosomes that tend to concentrate recombination in their distal regions
# have a lower chromosomal genetic shuffling rate
# The higher is the bias toward the periphery, the lower is the genetic shuffling, according to (1)
plot(chromosome.stats$peripherybias_ratio, chromosome.stats$geneticshuffling_genedist,
     xlab = "Periphery-bias ratio", ylab = "Genetic shuffling in gene distance")
# The result is not the expected. Chromosomes with a lower bias toward periphery

# Dig furthermore with models...
mod = lm(geneticshuffling_genedist ~ (peripherybias_ratio), data = chromosome.stats)
summary(mod)
AIC(mod)

# Linear Mixed Model with species effect
lmer.model = lmer(geneticshuffling_genedist ~ (peripherybias_ratio) + (1|species), data = chromosome.stats)
lmer.model = lmer(geneticshuffling_genedist ~ (peripherybias_ratio) + (sqrt(peripherybias_ratio)|species), data = chromosome.stats)
summary(lmer.model)
plot(lmer.model)
AIC(lmer.model)

# Diagnostic plots
par(mfrow = c(1, 3))
# Residuals distribution
qqPlot(residuals(lmer.model), ylab = "Residuals", id = F)
# Residuals as a function of fitted
plot(fitted(lmer.model), residuals(lmer.model), xlab = "fitted", ylab = "residuals")
abline(h = 0)
# Observed~predicted qqplot
plot(chromosome.stats$geneticshuffling_genedist[which(!is.na(chromosome.stats$peripherybias_ratio))], predict(lmer.model), xlab = "observed", ylab = "predicted")
abline(a = 0, b = 1)
par(mfrow = c(1,1))

plot(sqrt(chromosome.stats$peripherybias_ratio), chromosome.stats$geneticshuffling_genedist,
     xlab = "Periphery-bias ratio", ylab = "Genetic shuffling in gene distance")
abline(mod,
       col = "red")
abline(fixef(lmer.model)[1],
       fixef(lmer.model)[2],
       col = "blue")

#--------------------------------------------------------------------------#
# Is there a relationship between the genetic map length and the genetic shuffling rate?
plot(chromosome.stats$linkage.map.length.correctedHW, chromosome.stats$geneticshuffling_genedist,
     xlab = "Genetic map length", ylab = "Genetic shuffling in gene distance")
# Indeed, larger maps, with more COs, have a higher genetic shuffling rate.
# Dig furthermore with models...
mod = lm(geneticshuffling_genedist ~ linkage.map.length.correctedHW, data = chromosome.stats)
summary(mod)
AIC(mod)

# Linear Mixed Model with species effect
lmer.model = lmer(geneticshuffling_genedist ~ linkage.map.length.correctedHW + (1|species), data = chromosome.stats)
lmer.model = lmer(geneticshuffling_genedist ~ linkage.map.length.correctedHW + (linkage.map.length.correctedHW|species), data = chromosome.stats)
summary(lmer.model)
plot(lmer.model)
AIC(lmer.model)

# Diagnostic plots
par(mfrow = c(1, 3))
# Residuals distribution
qqPlot(residuals(lmer.model), ylab = "Residuals", id = F)
# Residuals as a function of fitted
plot(fitted(lmer.model), residuals(lmer.model), xlab = "fitted", ylab = "residuals")
abline(h = 0)
# Observed~predicted qqplot
plot(chromosome.stats$geneticshuffling_genedist, predict(lmer.model), xlab = "observed", ylab = "predicted")
abline(a = 0, b = 1)
par(mfrow = c(1,1))

plot(chromosome.stats$linkage.map.length.correctedHW, chromosome.stats$geneticshuffling_genedist,
     xlab = "Genetic map length", ylab = "Genetic shuffling in gene distance")
abline(mod,
       col = "red")
abline(fixef(lmer.model)[1],
       fixef(lmer.model)[2],
       col = "blue")




############################################################################ #
#============================================================================#
# ESTIMATE THE SPECIES INTRA-CHROMOSOMAL GENETIC SHUFFLING ----
#============================================================================#
############################################################################ #


#============================================================================#
# Genetic shuffling rate of species ----
#============================================================================#
# Estimate the intra-chromosomal part of the genetic shuffling rate, according to Veller et al. 2019
df_geneticshuffling_physical_chr = read.table("tables/df_geneticshuffling_physical_chr.txt", header = TRUE)

# Sampling problem: we don't have complete sampling of all chromosomes for each species !!

# Compute the mean genetic shuffling rate per species
list_set = unique(df.veller$set)
df_geneticshuffling_physical = data.frame(set = list_set, geneticshuffling = NA)
for (i in 1:nrow(df_geneticshuffling_physical)) {
  print(i)
  df_geneticshuffling_physical$geneticshuffling[i] = R_chr_mean(df.veller, df_geneticshuffling_physical$set[i])
}

# Save in tables
write.table(df_geneticshuffling_physical, file = "tables/df_geneticshuffling_physical.csv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
df_geneticshuffling_physical = read.table("tables/df_geneticshuffling_physical.csv", header = TRUE)




#============================================================================#
# END ----
#============================================================================#

