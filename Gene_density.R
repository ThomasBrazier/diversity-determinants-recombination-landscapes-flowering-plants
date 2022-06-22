############################################################################ #
#                    ECOBIO - PhD
#
#       Gene Density
#
############################################################################ #
# by Thomas Brazier
# brazier.thomas@gmail.com

# Supervisor: Sylvain GLEMIN
#             ECOBIO Lab


# Within this script you have functions to
# - estimate the gene count/density along the genome for a given recombination map
# - analyses of the correlation between gene count/density and recombination
# - template of figures

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
library(nlme)
library(pals)
library(egg)
library(data.table)
library(car)
library(ape)
library(dplyr)
library(partR2)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)

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
source("sources/Genome.R") 

metadata.clean = read.csv(file = paste(wd, "/data-cleaned/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")
chromosome.stats = read.table(file = "tables/chromosome.stats.csv", header = TRUE, sep = ";")

#============================================================================#
# Consider only protein-coding genes ----
#============================================================================#
# Which species have other than protein-coding genes in annotations?
# List of species
list_species = system(paste("ls ", wd, "/data-cleaned/genome/gene_positions/", sep = ""), intern = TRUE)
list_species = list_species[grep("*.txt.gz", list_species)]
list_species = gsub("_gene_positions.txt.gz", "", list_species)
list_species
# List of accessions
df = metadata.clean[which(metadata.clean$species %in% list_species),c("id", "species", "accession")]
df$coding = NA
# In data/Genome/species/accession/gff3/.gff.gz
# i=28
# sp = df$species[i]
# acc = df$accession[i]
# filename = system(paste("ls ", wd, paste("/data/Genome/", tolower(sp), "/", acc, "/gff3/", sep = ""), sep = ""), intern = TRUE)
# filename = filename[grep("*.gff*", filename)]
# gff = rtracklayer::import.gff(paste("data/Genome/", tolower(sp), "/", acc, "/gff3/", filename, sep = ""))
# #class(gff)
# #colnames(mcols(gff))
# #gff$gene_biotype
# tab = as.data.frame(table(gff$gene_biotype))
# df$noncoding[i] = tab$Freq[which(tab$Var1 == "protein_coding")]/sum(tab$Freq)

df = df[-which(df$species == "Aegilops_speltoides"),]
# df = df[-which(df$species == "Camellia_sinensis"),]
 
for (i in 1:nrow(df)) {
  sp = df$species[i]
  acc = df$accession[i]
  cat(sp, acc, "\n")
  filename = system(paste("ls ", wd, paste("/data/Genome/", tolower(sp), "/", acc, "/gff3/", sep = ""), sep = ""), intern = TRUE)
  filename = filename[grep(".gff", filename)]
  # filename = filename[-grep(".chr", filename)]
  gff = rtracklayer::import.gff(paste("data/Genome/", tolower(sp), "/", acc, "/gff3/", filename, sep = ""))
  #class(gff)
  #colnames(mcols(gff))
  #gff$gene_biotype
  # tab = as.data.frame(table(gff$gene_biotype))
  # df$coding[i] = tab$Freq[which(tab$Var1 == "protein_coding")]/sum(tab$Freq)
  
  if (length(table(gff$biotype)) != 0) {
    tab = as.data.frame(table(gff$biotype))
    df$coding[i] = tab$Freq[which(tab$Var1 == "protein_coding")]/sum(tab$Freq)
  } else {
    if (length(table(gff$gene_biotype)) != 0) {
      tab = as.data.frame(table(gff$gene_biotype))
      df$coding[i] = tab$Freq[which(tab$Var1 == "protein_coding")]/sum(tab$Freq)
    } else {
      df$coding[i] = NA
    }
  }
}
write.table(df, "data/Genome/coding_proportion.csv", col.names = T, row.names = F, quote = F, sep = "\t")

# Plot gene landscapes with only protein coding genes
gff_gene = gff[which(gff$type == "gene")]
gff_gene = gff_gene[which(gff_gene$biotype == "protein_coding")]
gff_gene = gff_gene[which(gff_gene@seqnames == 1)]

unique(gff_gene$biotype)
sum(gff_gene$biotype == "protein_coding")
dfgene = data.frame(start = gff_gene@ranges@start,
                    chr = 1)
# Bin in 100kb windows
dfplot = data.frame(win = seq(50000, max(dfgene$start), by = 100000),
          nbgenes = NA)
for (i in 1:nrow(dfplot)) {
  dfplot$nbgenes[i] = sum((dfgene$start > (dfplot$win[i] - 50000)) & (dfgene$start < (dfplot$win[i] + 50000)))
}

ggplot(data = dfplot, aes(x = win, y = nbgenes)) +
  # ggtitle(paste(chromosome.stats.shuffling$species[x], "chromosome", chromosome.stats.shuffling$chromosome[x], sep = " ")) +
  # xlim(0, max(c(genemap$phys, xmax), na.rm = TRUE)) +
  geom_point(size = 1.2, alpha = 0.4) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  xlab("Genomic position (Mb)") + ylab("Gene count (100kb windows)") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"))




#============================================================================#
# Estimating gene count in recombination maps ----
#============================================================================#
# Test the method on A. thaliana
set = "Arabidopsis_thaliana_Serin2017"
# set = "Zea_mays_MaizeGDBConsensus_v4"
# set = "Brachypodium_distachyon_Huo2011"
# set = "Oryza_sativa_Jiang2017"
# The method need the gene positions already retrieved in the reference genome
# in a file "data-cleaned/genome/gene_positions/", species,"_gene_positions.txt.gz"
# AND an existing loess 100kb recombination map

# Correction 2021/11/19
set = "Theobroma_cacao_Royaert2016"
source("sources/Genome.R")
gene_count(set = set, scale = 0.1)

# Repoduce to all species in dataset
# List of all maps
list = system(paste("ls ", wd, "/data-cleaned/genome/gene_positions/", sep = ""), intern = TRUE)
list = list[grep("*.txt.gz", list)]
list = gsub("_gene_positions.txt.gz", "", list)
list
# Dataset with a gene positions map (i.e. a gff annotation file)
list.set = as.character(metadata.clean$id[which(metadata.clean$species %in% list)])
source("sources/Genome.R") 
for (set in list.set) {
  gene_count(set, scale = 0.1)
}

#============================================================================#
# CORRELATION BETWEEN GENE COUNT AND RECOMBINATION ----
#============================================================================#

# Correlation between gene density and recombination rate
# set = "Zea_mays_MaizeGDBConsensus_v4"
# set = "Aegilops_speltoides_Zhang2020"
chromosome = "1"
genecount_map = read.table(file = gzfile(paste("data-cleaned/genome/gene_count/", set, "_chromosome", chromosome, ".txt.gz", sep = "")),
                             header = TRUE, sep = "\t")

summary(genecount_map$gene_count)
hist(genecount_map$gene_count)
# TODO Low number of genes in each window, may be increase the scale !!!
# Distribution of genes along the chromosome
plot(genecount_map$phys, genecount_map$gene_count)

#============================================================================#
# Estimating the mean recombination rate for each gene count, for each species ----
#============================================================================#

# A function that estimate the mean recombination rate for each gene count in a given dataset
estimmeanrecrate_genecount = function(species, scaled = FALSE) {
  library(dplyr)
  library(data.table)
  metadata.clean = read.csv(file = paste(wd, "/data-cleaned/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";", stringsAsFactors = FALSE)
  chromosome.stats = read.table(file = "tables/chromosome.stats.csv", header = TRUE, sep = ";", stringsAsFactors = FALSE)
  # How many chromosomes for the species
  list_chromosomes = chromosome.stats$chromosome[which(as.character(chromosome.stats$species) == species)]
  set = unique(chromosome.stats$set[which(as.character(chromosome.stats$species) == species)])
  # Import data for ALL chromosomes
  files = unlist(lapply(set, function (x) list.files(path = "data-cleaned/genome/gene_count/", pattern = x)))
  df = bind_rows(lapply(files, function(x) read.table(file = gzfile(paste("data-cleaned/genome/gene_count/", x, sep = "")),
                                                                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)))
  # Estimate the mean recombination ± C.I. for each gene count
  genecount_levels = sort(unique(df$gene_count))
  # Define the bootstrap function, taking a vector of recombination rates values as argument
  bootmean = function(recrates, nboot = 1000) {
    boot = numeric(nboot)
    for (i in 1:nboot) {
      boot[i] = mean(sample(recrates, replace = TRUE), na.rm = TRUE)
    }
    results = data.frame(mean_rec = mean(boot, na.rm = TRUE),
                         lower_rec = quantile(boot, 0.025, na.rm = TRUE),
                         upper_rec = quantile(boot, 0.975, na.rm = TRUE))
    return(results)
  }
  # Estimate the recombination
  res = data.frame(species = rep(gsub(" ", "_", species), length(genecount_levels)), gene_count = genecount_levels)
  res = cbind(res, bind_rows(lapply(res$gene_count, function(x) bootmean(df$rec.rate[which(df$gene_count == x)]))))
  # Estimate the mean recombination rate (all data) and scale the recombination rates per gene count
  if (scaled == TRUE) {
    # meanrecrate = bootmean(df$rec.rate)
    # res$mean_rec = res$mean_rec/meanrecrate$mean_rec
    # res$lower_rec = res$lower_rec/meanrecrate$mean_rec
    # res$upper_rec = res$upper_rec/meanrecrate$mean_rec
    res$mean_rec = scale(res$mean_rec)
    res$lower_rec = scale(res$lower_rec)
    res$upper_rec = scale(res$upper_rec)
  }
  
  # Save in gene count
  if (scaled == TRUE) {
    savedfile = read.table("output/gene_count/meanrecrate_genecount_scaled.txt", header = TRUE)
    # Remove older results
    savedfile = savedfile[which(savedfile$species != gsub(" ", "_", species)),]
    savedfile = rbind(savedfile, res)
    write.table(savedfile, "output/gene_count/meanrecrate_genecount_scaled.txt",
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  } else {
    savedfile = read.table("output/gene_count/meanrecrate_genecount.txt", header = TRUE)
    # Remove older results
    savedfile = savedfile[which(savedfile$species != gsub(" ", "_", species)),]
    savedfile = rbind(savedfile, res)
    write.table(savedfile, "output/gene_count/meanrecrate_genecount.txt",
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  }
  
    return(res)
}

source("sources/Genome.R")
# species = "Arabidopsis thaliana"
species = "Theobroma cacao"
estimmeanrecrate_genecount(species)

# Batch for all species
list_species = system(paste("ls ", wd, "/data-cleaned/genome/gene_positions/", sep = ""), intern = TRUE)
list_species = list_species[grep("*.txt.gz", list_species)]
list_species = gsub("_gene_positions.txt.gz", "", list_species)
(list_species = gsub("_", " ", list_species))
# Compute the recombination rate per gene count
for (sp in list_species) {
  cat(sp, "\n")
  estimmeanrecrate_genecount(sp)
}
# Compute the scaled recombination rate per gene count
for (sp in list_species) {
  cat(sp, "\n")
  estimmeanrecrate_genecount(sp, scaled = TRUE)
}

#============================================================================#
# Rec rate ~ gene count, for each species, 2X2 plot ----
#============================================================================#
meanrecrate_genecount = read.table("output/gene_count/meanrecrate_genecount.txt", header = TRUE)
meanrecrate_genecount = read.table("output/gene_count/meanrecrate_genecount_scaled.txt", header = TRUE)

meanrecrate_genecount$species = gsub("_", " ", meanrecrate_genecount$species)
unique(meanrecrate_genecount$species)

# FIGURES
# Divide the figure in 4 plots for visualization (44 species)
# Divide in four subsets
d1 = meanrecrate_genecount[which(meanrecrate_genecount$species %in% unique(meanrecrate_genecount$species)[1:11]),]
d2 = meanrecrate_genecount[which(meanrecrate_genecount$species %in% unique(meanrecrate_genecount$species)[12:22]),]
d3 = meanrecrate_genecount[which(meanrecrate_genecount$species %in% unique(meanrecrate_genecount$species)[23:33]),]
d4 = meanrecrate_genecount[which(meanrecrate_genecount$species %in% unique(meanrecrate_genecount$species)[34:44]),]

# Add a x jitter per species
spe = unique(d1$species)
i = 0
for (s in spe) {
  cat(s, "\n")
  d1$gene_count[which(d1$species == s)] = d1$gene_count[which(d1$species == s)] + i
  i = i + 0.04
}
spe = unique(d2$species)
i = 0
for (s in spe) {
  cat(s, "\n")
  d2$gene_count[which(d2$species == s)] = d2$gene_count[which(d2$species == s)] + i
  i = i + 0.04
}
spe = unique(d3$species)
i = 0
for (s in spe) {
  cat(s, "\n")
  d3$gene_count[which(d3$species == s)] = d3$gene_count[which(d3$species == s)] + i
  i = i + 0.04
}
spe = unique(d4$species)
i = 0
for (s in spe) {
  cat(s, "\n")
  d4$gene_count[which(d4$species == s)] = d4$gene_count[which(d4$species == s)] + i
  i = i + 0.04
}

p1 = ggplot(data = d1, aes(x = gene_count, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  scale_x_continuous(limits = c(0, 20)) +
  # scale_y_continuous(limits = c(0, 7)) +
  labs(x = "Gene count", y = "Recombination rate", colour = "Species") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16))
p1
p2 = ggplot(data = d2, aes(x = gene_count, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  scale_x_continuous(limits = c(0, 20)) +
  # scale_y_continuous(limits = c(0, 7)) +
  labs(x = "Gene count", y = "Recombination rate", colour = "Species") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16))
# p2
p3 = ggplot(data = d3, aes(x = gene_count, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  scale_x_continuous(limits = c(0, 20)) +
  # scale_y_continuous(limits = c(0, 7)) +
  labs(x = "Gene count", y = "Recombination rate", colour = "Species") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16))
# p3
p4 = ggplot(data = d4, aes(x = gene_count, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  scale_x_continuous(limits = c(0, 20)) +
  # scale_y_continuous(limits = c(0, 7)) +
  labs(x = "Gene count", y = "Recombination rate", colour = "Species") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16))
# p4
ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

#============================================================================#
# Rec rate ~ gene count, all species pooled (dark grey transparent dots) + quadratic phylogenetic regression with C.I. ----
#============================================================================#
meanrecrate_genecount = read.table("output/gene_count/meanrecrate_genecount.txt", header = TRUE)
meanrecrate_genecount$species = gsub("_", " ", meanrecrate_genecount$species)
# Subset only for a gene count <= 20
meanrecrate_genecount = meanrecrate_genecount[which(meanrecrate_genecount$gene_count < 21),]

summary(meanrecrate_genecount$mean_rec)
# Quadratic linear regression
mod = lm(mean_rec ~ gene_count + I(gene_count^2), data = meanrecrate_genecount)
summary(mod)
plot(mod, which = c(1,2,3,4))
AIC(mod)
# Poor model quality with some strong outliers
# Yet if it is only for visualisation

# Quadratic Phylogenetic Regression
# Load the phylogenetic tree estimated in 'Phylogeny.R'
load("output/phylogeny/tree.Rdata")
tree
tree$tip.label
# Drop unsampled tips
library(ape)
tree = drop.tip(tree, tip = tree$tip.label[!(tree$tip.label %in% meanrecrate_genecount$species)])
tree

# PGLMM (model selection)
# pglmm.model = pglmm(mean_rec ~ gene_count + (1|species), data = meanrecrate_genecount, family = 'gaussian',
#                     cov_ranef = list(species = tree))
library(phyr)
pglmm.model = pglmm(mean_rec ~ gene_count + I(gene_count^2) + (1|species), data = meanrecrate_genecount, family = 'gaussian',
                    cov_ranef = list(species = tree))
pglmm.model = pglmm(mean_rec ~ gene_count + I(gene_count^2) + (1|species__), data = meanrecrate_genecount, family = 'gaussian',
                    cov_ranef = list(species = tree))
summary(pglmm.model)
# No clear differences in AIC with/without phylogeny
# Test the effect of the phylogeny
# Bayesian estimates for confidence interval (Bayesian INLA)
pglmm.model = pglmm(mean_rec ~ gene_count + I(gene_count^2) + (1|species__), data = meanrecrate_genecount, family = 'gaussian',
                    cov_ranef = list(species = tree), bayes = TRUE)
summary(pglmm.model)
# The model show a significant effect of the phylogeny (yet the model without phylogeny seems better)
# 31% of variance explained by phylogeny

# Hence phylogenetic autocorrelation should not be considered furthermore

# Diagnostic plots
# Residuals distribution
qqPlot(residuals(pglmm.model))
# Residuals as a function of fitted
plot(fitted(pglmm.model), residuals(pglmm.model))
abline(h = 0)
# Observed~predicted qqplot
plot(unlist(pglmm_predicted_values(pglmm.model)), pglmm.model$Y)
abline(a = 0, b = 1)

# Bayesian estimates for confidence interval (Bayesian INLA)
pglmm.model = pglmm(mean_rec ~ gene_count + I(gene_count^2) + (1|species__), data = meanrecrate_genecount, family = 'gaussian',
                    cov_ranef = list(species = tree), bayes = TRUE)
summary(pglmm.model)
# Quadratic function and C.I.
quadrafun = function(x, c) fixef(pglmm.model)[[1]][1] + fixef(pglmm.model)[[1]][2]*x + fixef(pglmm.model)[[1]][3]*x^2
lowerquadrafun = function(x, c) fixef(pglmm.model)[[2]][1] + fixef(pglmm.model)[[2]][2]*x + fixef(pglmm.model)[[2]][3]*x^2
upperquadrafun = function(x, c) fixef(pglmm.model)[[3]][1] + fixef(pglmm.model)[[3]][2]*x + fixef(pglmm.model)[[3]][3]*x^2
# The C.I. values
df.new = data.frame(gene_count = 0:20)
df.new$lwr.pred = lowerquadrafun(df.new$gene_count)
df.new$upr.pred = upperquadrafun(df.new$gene_count)

# FIGURES
p = ggplot(data = meanrecrate_genecount, aes(x = gene_count, y = mean_rec)) +
  geom_point(colour = "Darkgrey", alpha = 0.5) +
  stat_smooth(aes(y = mean_rec), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Blue", fill = "Blue", alpha = 0.3) +
  scale_x_continuous(limits = c(0, 20)) +
  # scale_y_continuous(limits = c(0, 6)) +
  labs(x = "Gene count", y = "Recombination rate") +
  # geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
  stat_function(fun = quadrafun, colour = "Black", size = 1) +
  geom_ribbon(data = df.new, aes(x = gene_count, ymin = lwr.pred, ymax = upr.pred), alpha = 0.4, inherit.aes = F, fill = "Black") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16))
p

#----------------------------------------------------------------------------#
# Diagnostic plots ----
#----------------------------------------------------------------------------#



#============================================================================#
# Rec rate ~ gene count, all species pooled (dark grey transparent dots), + bootstrapped mean ± C.I. (pooled species) ----
#============================================================================#
meanrecrate_genecount = read.table("output/gene_count/meanrecrate_genecount.txt", header = TRUE)
meanrecrate_genecount$species = gsub("_", " ", meanrecrate_genecount$species)

# Compute mean and C.I. for each gene count
df = data.frame(gene_count = c(0:40), mean_rec = NA, lower_rec = NA, upper_rec = NA)

for (i in 1:nrow(df)) {
  nboot = 1000
  boot = numeric(nboot)
  for (j in 1:nboot) {
    boot[j] = mean(sample(meanrecrate_genecount$mean_rec[which(meanrecrate_genecount$gene_count == df$gene_count[i])], replace = TRUE), na.rm = TRUE)
  }
  df$mean_rec[i] = mean(boot, na.rm = TRUE)
  df$lower_rec[i] = quantile(boot, 0.025, na.rm = TRUE)
  df$upper_rec[i] = quantile(boot, 0.975, na.rm = TRUE)
}

df

# FIGURES
p = ggplot(data = meanrecrate_genecount, aes(x = gene_count, y = mean_rec)) +
  geom_point(colour = "Darkgrey", alpha = 0.5) +
  # stat_smooth(aes(y = mean_rec), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Black", fill = "Grey12") +
  scale_x_continuous(limits = c(0, 40)) +
  scale_y_continuous(limits = c(0, 10)) +
  labs(x = "Gene count", y = "Recombination rate") +
  # geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
  geom_line(data = df, aes(x = gene_count, y = mean_rec)) +
  geom_ribbon(data = df, aes(x = gene_count, ymin = lower_rec, ymax = upper_rec), alpha = 0.4, inherit.aes = F, fill = "Black") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16))
p





#============================================================================#
# CORRELATION BETWEEN GENE DENSITY (PROPORTION OF GENE SEQUENCE) AND RECOMBINATION ----
#============================================================================#

#============================================================================#
# Estimating gene density in recombination maps ----
#============================================================================#
# Test the method on A. thaliana
set = "Arabidopsis_thaliana_Serin2017"
set = "Aegilops_speltoides_Zhang2019"
# The method need the gene positions already retrieved in the reference genome
# in a file "data-cleaned/genome/gene_positions/", species,"_gene_positions.txt.gz"
# AND an existing loess 100kb recombination map
source("sources/Genome.R")
gene_density(set = set, scale = 0.1)

#----------------------------------------------------------------------------#
# Repoduce to all species in dataset
# List of all maps
list = system(paste("ls ", wd, "/data-cleaned/genome/gene_positions/", sep = ""), intern = TRUE)
list = list[grep("*.txt.gz", list)]
list = gsub("_gene_positions.txt.gz", "", list)
list
# Dataset with a gene positions map (i.e. a gff annotation file)
list.set = as.character(metadata.clean$id[which(metadata.clean$species %in% list)])
source("sources/Genome.R") 
for (set in list.set) {
  gene_density(set, scale = 0.1)
}

#============================================================================#
# Computing recombination rates in 20 quantiles of gene density ----
#============================================================================#
# A function that estimate the mean recombination rate for each gene density quantile (N=20 quantiles) in a given dataset
# When divbygenecount is true, the windowed gene density is divided by the gene count of the same window
estimmeanrecrate_genedensity = function(species, rec.scaled = TRUE, density.scaled = TRUE, nquantiles = 20, sampling = c("quantiles", "systematic"), divbygenecount = FALSE) {
  library(dplyr)
  library(data.table)
  metadata.clean = read.csv(file = paste(wd, "/data-cleaned/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";", stringsAsFactors = FALSE)
  chromosome.stats = read.table(file = "tables/chromosome.stats.csv", header = TRUE, sep = ";", stringsAsFactors = FALSE)
  # How many chromosomes for the species
  list_chromosomes = chromosome.stats$chromosome[which(as.character(chromosome.stats$species) == species)]
  set = unique(chromosome.stats$set[which(as.character(chromosome.stats$species) == species)])
  # Import data for ALL chromosomes
  files = unlist(lapply(set, function (x) list.files(path = "data-cleaned/genome/gene_density/", pattern = x)))
  df = bind_rows(lapply(files, function(x) read.table(file = gzfile(paste("data-cleaned/genome/gene_density/", x, sep = "")),
                                                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)))
  # summary(df$gene_density)
  
  if (divbygenecount == TRUE) {
    # Import gene count
    files = unlist(lapply(set, function (x) list.files(path = "data-cleaned/genome/gene_count/", pattern = x)))
    df_genecount = bind_rows(lapply(files, function(x) read.table(file = gzfile(paste("data-cleaned/genome/gene_count/", x, sep = "")),
                                                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)))
    # Divide
    df$gene_density = df$gene_density/df_genecount$gene_count
    df$gene_density[which(df$gene_density == "NaN")] = NA
    df$gene_density[which(df$gene_density == "Inf")] = NA
  }
  
  if (density.scaled == TRUE) {
    df$gene_density = scale(df$gene_density, center = TRUE)
    # df$gene_density = df$gene_density/max(df$gene_density, na.rm = TRUE)
  }
  # hist(df$gene_density)
  
  #--------------------------------------------------#
  # Compute quantiles boundaries for gene density
  # N=20 quantiles within the distribution of gene density OR equally distributed from 0 to 1
  if (sampling == "quantiles") {
    quantiles = c(0, quantile(df$gene_density[which(df$gene_density >0)], probs = seq(0, 1, length.out = nquantiles), na.rm = TRUE))
  } else {
    if (sampling == "systematic") {
      quantiles = seq(0, 1, length.out = nquantiles+1)
    }
  }
  # Estimate the mean recombination ± C.I. for each quantile
  # Define the bootstrap function, taking a vector of recombination rates values as argument
  bootmean = function(recrates, nboot = 1000) {
    boot = numeric(nboot)
    for (i in 1:nboot) {
      boot[i] = mean(sample(recrates, replace = TRUE), na.rm = TRUE)
    }
    results = data.frame(mean_rec = mean(boot, na.rm = TRUE),
                         lower_rec = quantile(boot, 0.025, na.rm = TRUE),
                         upper_rec = quantile(boot, 0.975, na.rm = TRUE))
    return(results)
  }
  # Estimate the recombination
  res = data.frame(species = rep(gsub(" ", "_", species), nquantiles), quantile = c(1:20),
                   quantile.lower = quantiles[1:20], quantile.upper = quantiles[2:21])
  res = cbind(res, bind_rows(lapply(res$quantile,
                                    function(x) bootmean(df$rec.rate[which(df$gene_density >= res$quantile.lower[x] & df$gene_density <= res$quantile.upper[x])]))))
  # Estimate the mean recombination rate (all data) and scale the recombination rates per gene count
  if (rec.scaled == TRUE) {
    # meanrecrate = bootmean(df$rec.rate)
    # res$mean_rec = res$mean_rec/meanrecrate$mean_rec
    # res$lower_rec = res$lower_rec/meanrecrate$mean_rec
    # res$upper_rec = res$upper_rec/meanrecrate$mean_rec
    res$mean_rec = scale(res$mean_rec, center = TRUE)
    res$lower_rec = scale(res$lower_rec, center = TRUE)
    res$upper_rec = scale(res$upper_rec, center = TRUE)
  }
  
  # Save in gene count
  if(divbygenecount == FALSE) {
    savedfile = read.table("output/gene_density/meanrecrate_genedensity.txt", header = TRUE, sep = "\t")
    # Remove older results
    savedfile = savedfile[which(savedfile$species != gsub(" ", "_", species)),]
    savedfile = rbind(savedfile, res)
    write.table(savedfile, "output/gene_density/meanrecrate_genedensity.txt",
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  } else {
    savedfile = read.table("output/gene_density/meanrecrate_genedensity_divbygenecount.txt", header = TRUE, sep = "\t")
    # Remove older results
    savedfile = savedfile[which(savedfile$species != gsub(" ", "_", species)),]
    savedfile = rbind(savedfile, res)
    write.table(savedfile, "output/gene_density/meanrecrate_genedensity_divbygenecount.txt",
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  }

  return(res)
}

source("sources/Genome.R")
species = "Aegilops speltoides"
estimmeanrecrate_genedensity(species)

# Batch for all species
list_species = system(paste("ls ", wd, "/data-cleaned/genome/gene_positions/", sep = ""), intern = TRUE)
list_species = list_species[grep("*.txt.gz", list_species)]
list_species = gsub("_gene_positions.txt.gz", "", list_species)
(list_species = gsub("_", " ", list_species))
# Compute the scaled recombination rate per gene count
for (sp in list_species) {
  cat(sp, "\n")
  # estimmeanrecrate_genedensity(sp, rec.scaled = TRUE, density.scaled = TRUE, sampling = "quantiles", divbygenecount = TRUE)
  estimmeanrecrate_genedensity(sp, rec.scaled = TRUE, density.scaled = FALSE, sampling = "systematic", divbygenecount = TRUE)
  estimmeanrecrate_genedensity(sp, rec.scaled = TRUE, density.scaled = FALSE, sampling = "systematic", divbygenecount = FALSE)
  # estimmeanrecrate_genedensity(sp, rec.scaled = FALSE, density.scaled = FALSE, sampling = "systematic", divbygenecount = TRUE)
  # estimmeanrecrate_genedensity(sp, rec.scaled = FALSE, density.scaled = FALSE, sampling = "systematic", divbygenecount = FALSE)
}

#============================================================================#
# Descriptive model ----
#============================================================================#
# Gene density
meanrecrate_genedensity = read.table("output/gene_density/meanrecrate_genedensity.txt", header = TRUE, sep = "\t")
# Or gene density divided by the number of genes
meanrecrate_genedensity = read.table("output/gene_density/meanrecrate_genedensity_divbygenecount.txt", header = TRUE, sep = "\t")

meanrecrate_genedensity$species = gsub("_", " ", meanrecrate_genedensity$species)
unique(meanrecrate_genedensity$species)
meanrecrate_genedensity$gene_density = rowMeans(cbind(meanrecrate_genedensity$quantile.lower, meanrecrate_genedensity$quantile.upper), na.rm = TRUE)

# Correlation
cor.test(meanrecrate_genedensity$mean_rec, meanrecrate_genedensity$gene_density, method = "spearman")

# PGLMM regression
# Quadratic Phylogenetic Regression
# Load the phylogenetic tree estimated in 'Phylogeny.R'
load("output/phylogeny/tree.Rdata")
# Drop unsampled tips
library(ape)
tree = drop.tip(tree, tip = tree$tip.label[!(tree$tip.label %in% meanrecrate_genedensity$species)])
library(phyr)
# PGLMM (model selection)
# pglmm.model = pglmm(mean_rec ~ gene_density + (1|species), data = meanrecrate_genedensity, family = 'gaussian',
#                     cov_ranef = list(species = tree))
pglmm.model = pglmm(mean_rec ~ gene_density + I(gene_density^2) + (1|species), data = meanrecrate_genedensity, family = 'gaussian',
                    cov_ranef = list(species = tree))
pglmm.model = pglmm(mean_rec ~ gene_density + I(gene_density^2) + (1|species__), data = meanrecrate_genedensity, family = 'gaussian',
                    cov_ranef = list(species = tree))
summary(pglmm.model)
# No clear differences in AIC with/without phylogeny
# Test the effect of the phylogeny
# Bayesian estimates for confidence interval (Bayesian INLA)
pglmm.model = pglmm(mean_rec ~ gene_density + I(gene_density^2) + (1|species__), data = meanrecrate_genedensity, family = 'gaussian',
                    cov_ranef = list(species = tree), bayes = TRUE)
summary(pglmm.model)
# The model show a significant effect of the phylogeny (yet the model without phylogeny seems better)

# Hence phylogenetic autocorrelation should not be considered furthermore

# Diagnostic plots
# Residuals distribution
library(car)
qqPlot(residuals(pglmm.model))
# Residuals as a function of fitted
plot(fitted(pglmm.model), residuals(pglmm.model))
abline(h = 0)
# Observed~predicted qqplot
plot(unlist(pglmm_predicted_values(pglmm.model)), pglmm.model$Y)
abline(a = 0, b = 1)

# Bayesian estimates for confidence interval (Bayesian INLA)
pglmm.model = pglmm(mean_rec ~ gene_density + I(gene_density^2) + (1|species__), data = meanrecrate_genedensity, family = 'gaussian',
                    cov_ranef = list(species = tree), bayes = TRUE)
summary(pglmm.model)
# Quadratic function and C.I.
quadrafun = function(x, c) fixef(pglmm.model)[[1]][1] + fixef(pglmm.model)[[1]][2]*x + fixef(pglmm.model)[[1]][3]*x^2
lowerquadrafun = function(x, c) fixef(pglmm.model)[[2]][1] + fixef(pglmm.model)[[2]][2]*x + fixef(pglmm.model)[[2]][3]*x^2
upperquadrafun = function(x, c) fixef(pglmm.model)[[3]][1] + fixef(pglmm.model)[[3]][2]*x + fixef(pglmm.model)[[3]][3]*x^2
# The C.I. values
df.new = data.frame(gene_density = seq(0, 1, length.out = nquantiles+1))
df.new$lwr.pred = lowerquadrafun(df.new$gene_density)
df.new$upr.pred = upperquadrafun(df.new$gene_density)

# FIGURES
p = ggplot(data = meanrecrate_genedensity, aes(x = gene_density, y = mean_rec)) +
  geom_point(colour = "Darkgrey", alpha = 0.5) +
  stat_smooth(aes(y = mean_rec), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Blue", fill = "Blue", alpha = 0.3) +
  stat_smooth(aes(y = mean_rec), method = "lm", formula = y ~ x, size = 1, colour = "Red", fill = "Red", alpha = 0.3) +
  scale_x_continuous(limits = c(0, 1)) +
  # scale_y_continuous(limits = c(-5, 5)) +
  labs(x = "Gene density", y = "Recombination rate") +
  # geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
  stat_function(fun = quadrafun, colour = "Black", size = 1) +
  geom_ribbon(data = df.new, aes(x = gene_density, ymin = lwr.pred, ymax = upr.pred), alpha = 0.4, inherit.aes = F, fill = "Black") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16))
p


#============================================================================#
# Rec rate ~ gene density quantiles, for each species, 2X2 plot ----
#============================================================================#
# Gene density
meanrecrate_genedensity = read.table("output/gene_density/meanrecrate_genedensity.txt", header = TRUE, sep = "\t")
# Or gene density divided by the number of genes
meanrecrate_genedensity = read.table("output/gene_density/meanrecrate_genedensity_divbygenecount.txt", header = TRUE, sep = "\t")

meanrecrate_genedensity$species = gsub("_", " ", meanrecrate_genedensity$species)
unique(meanrecrate_genedensity$species)
meanrecrate_genedensity$gene_density = rowMeans(cbind(meanrecrate_genedensity$quantile.lower, meanrecrate_genedensity$quantile.upper), na.rm = TRUE)

# FIGURES
# Divide the figure in 4 plots for visualization (44 species)
# Divide in four subsets
d1 = meanrecrate_genedensity[which(meanrecrate_genedensity$species %in% unique(meanrecrate_genedensity$species)[1:11]),]
d2 = meanrecrate_genedensity[which(meanrecrate_genedensity$species %in% unique(meanrecrate_genedensity$species)[12:22]),]
d3 = meanrecrate_genedensity[which(meanrecrate_genedensity$species %in% unique(meanrecrate_genedensity$species)[23:33]),]
d4 = meanrecrate_genedensity[which(meanrecrate_genedensity$species %in% unique(meanrecrate_genedensity$species)[34:44]),]

p1 = ggplot(data = d1, aes(x = gene_density, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  # scale_x_continuous(limits = c(0, 20)) +
  # scale_y_continuous(limits = c(0, 7)) +
  labs(x = "Gene density", y = "Recombination rate", colour = "Species") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16))
# p1
p2 = ggplot(data = d2, aes(x = gene_density, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  # scale_x_continuous(limits = c(0, 20)) +
  # scale_y_continuous(limits = c(0, 7)) +
  labs(x = "Gene density", y = "Recombination rate", colour = "Species") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16))
# p2
p3 = ggplot(data = d3, aes(x = gene_density, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  # scale_x_continuous(limits = c(0, 20)) +
  # scale_y_continuous(limits = c(0, 7)) +
  labs(x = "Gene density", y = "Recombination rate", colour = "Species") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16))
# p3
p4 = ggplot(data = d4, aes(x = gene_density, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  # scale_x_continuous(limits = c(0, 20)) +
  # scale_y_continuous(limits = c(0, 7)) +
  labs(x = "Gene density", y = "Recombination rate", colour = "Species") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16))
# p4
ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)











#============================================================================#
# CORRELATION BETWEEN CODING DENSITY (PROPORTION OF CODING SEQUENCE) AND RECOMBINATION ----
#============================================================================#

#============================================================================#
# Computing recombination rates in 20 quantiles of gene density ----
#============================================================================#
# A function that estimate the mean recombination rate for each gene density quantile (N=20 quantiles) in a given dataset
estimmeanrecrate_codingdensity = function(species, rec.scaled = TRUE, sampling = c("quantiles", "systematic"), nquantiles = 20) {
  library(dplyr)
  library(data.table)
  metadata.clean = read.csv(file = paste(wd, "/data-cleaned/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";", stringsAsFactors = FALSE)
  chromosome.stats = read.table(file = "tables/chromosome.stats.csv", header = TRUE, sep = ";", stringsAsFactors = FALSE)
  # How many chromosomes for the species
  list_chromosomes = chromosome.stats$chromosome[which(as.character(chromosome.stats$species) == species)]
  set = unique(chromosome.stats$set[which(as.character(chromosome.stats$species) == species)])
  # Import data for ALL chromosomes
  files = unlist(lapply(set, function (x) list.files(path = "data-cleaned/genome/gene_density/", pattern = x)))
  df = bind_rows(lapply(files, function(x) read.table(file = gzfile(paste("data-cleaned/genome/gene_density/", x, sep = "")),
                                                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)))
  # summary(df$coding_density)
  #--------------------------------------------------#
  # Compute quantiles boundaries for gene density
  # N=20 quantiles within the distribution of gene density OR equally distributed from 0 to 1
  if (sampling == "quantiles") {
    quantiles = c(0, quantile(df$gene_density[which(df$gene_density >0)], probs = seq(0, 1, length.out = nquantiles), na.rm = TRUE))
  } else {
    if (sampling == "systematic") {
      quantiles = seq(0, 1, length.out = nquantiles+1)
    }
  }
  # Estimate the mean recombination ± C.I. for each quantile
  # Define the bootstrap function, taking a vector of recombination rates values as argument
  bootmean = function(recrates, nboot = 1000) {
    boot = numeric(nboot)
    for (i in 1:nboot) {
      boot[i] = mean(sample(recrates, replace = TRUE), na.rm = TRUE)
    }
    results = data.frame(mean_rec = mean(boot, na.rm = TRUE),
                         lower_rec = quantile(boot, 0.025, na.rm = TRUE),
                         upper_rec = quantile(boot, 0.975, na.rm = TRUE))
    return(results)
  }
  # Estimate the recombination
  res = data.frame(species = rep(gsub(" ", "_", species), nquantiles), quantile = c(1:20),
                   quantile.lower = quantiles[1:20], quantile.upper = quantiles[2:21])
  res = cbind(res, bind_rows(lapply(res$quantile,
                                    function(x) bootmean(df$rec.rate[which(df$coding_density >= res$quantile.lower[x] & df$coding_density <= res$quantile.upper[x])]))))
  # Estimate the mean recombination rate (all data) and scale the recombination rates per gene count
  if (rec.scaled == TRUE) {
    res$mean_rec = scale(res$mean_rec, center = TRUE)
    res$lower_rec = scale(res$lower_rec, center = TRUE)
    res$upper_rec = scale(res$upper_rec, center = TRUE)
  }
  
  # Save in gene count
  savedfile = read.table("output/gene_density/meanrecrate_codingdensity.txt", header = TRUE, sep = "\t")
  # Remove older results
  savedfile = savedfile[which(savedfile$species != gsub(" ", "_", species)),]
  savedfile = rbind(savedfile, res)
  write.table(savedfile, "output/gene_density/meanrecrate_codingdensity.txt",
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  return(res)
}

# Batch for all species
list_species = system(paste("ls ", wd, "/data-cleaned/genome/gene_positions/", sep = ""), intern = TRUE)
list_species = list_species[grep("*.txt.gz", list_species)]
list_species = gsub("_gene_positions.txt.gz", "", list_species)
(list_species = gsub("_", " ", list_species))
# Compute the scaled recombination rate per gene count
for (sp in list_species) {
  cat(sp, "\n")
  estimmeanrecrate_codingdensity(sp, rec.scaled = TRUE, sampling = "systematic")
}

#============================================================================#
# Descriptive model ----
#============================================================================#
meanrecrate_codingdensity = read.table("output/gene_density/meanrecrate_codingdensity.txt", header = TRUE, sep = "\t")
meanrecrate_codingdensity$species = gsub("_", " ", meanrecrate_codingdensity$species)
unique(meanrecrate_codingdensity$species)
meanrecrate_codingdensity$coding_density = rowMeans(cbind(meanrecrate_codingdensity$quantile.lower, meanrecrate_codingdensity$quantile.upper), na.rm = TRUE)

# Correlation
cor.test(meanrecrate_codingdensity$mean_rec, meanrecrate_codingdensity$coding_density, method = "spearman")

# Subset to coding density <= 0.5 (few values above 0.5)
meanrecrate_codingdensity = subset(meanrecrate_codingdensity, meanrecrate_codingdensity$coding_density <= 0.5)

# PGLMM regression
# Quadratic Phylogenetic Regression
# Load the phylogenetic tree estimated in 'Phylogeny.R'
load("output/phylogeny/tree.Rdata")
# Drop unsampled tips
library(ape)
tree = drop.tip(tree, tip = tree$tip.label[!(tree$tip.label %in% meanrecrate_codingdensity$species)])
library(phyr)
# PGLMM (model selection)
# pglmm.model = pglmm(mean_rec ~ coding_density + (1|species), data = meanrecrate_codingdensity, family = 'gaussian',
#                     cov_ranef = list(species = tree))
pglmm.model = pglmm(mean_rec ~ coding_density + I(coding_density^2) + (1|species), data = meanrecrate_codingdensity, family = 'gaussian',
                    cov_ranef = list(species = tree))
pglmm.model = pglmm(mean_rec ~ coding_density + I(coding_density^2) + (1|species__), data = meanrecrate_codingdensity, family = 'gaussian',
                    cov_ranef = list(species = tree))
summary(pglmm.model)
# No clear differences in AIC with/without phylogeny
# Test the effect of the phylogeny
# Bayesian estimates for confidence interval (Bayesian INLA)
pglmm.model = pglmm(mean_rec ~ coding_density + I(coding_density^2) + (1|species__), data = meanrecrate_codingdensity, family = 'gaussian',
                    cov_ranef = list(species = tree), bayes = TRUE)
summary(pglmm.model)
# The model show a significant effect of the phylogeny (yet the model without phylogeny seems better)

# Hence phylogenetic autocorrelation should not be considered furthermore

# Diagnostic plots
# Residuals distribution
library(car)
qqPlot(residuals(pglmm.model))
# Residuals as a function of fitted
plot(fitted(pglmm.model), residuals(pglmm.model))
abline(h = 0)
# Observed~predicted qqplot
plot(unlist(pglmm_predicted_values(pglmm.model)), pglmm.model$Y)
abline(a = 0, b = 1)

# Bayesian estimates for confidence interval (Bayesian INLA)
pglmm.model = pglmm(mean_rec ~ coding_density + I(coding_density^2) + (1|species__), data = meanrecrate_codingdensity, family = 'gaussian',
                    cov_ranef = list(species = tree), bayes = TRUE)
summary(pglmm.model)
# Quadratic function and C.I.
quadrafun = function(x, c) fixef(pglmm.model)[[1]][1] + fixef(pglmm.model)[[1]][2]*x + fixef(pglmm.model)[[1]][3]*x^2
lowerquadrafun = function(x, c) fixef(pglmm.model)[[2]][1] + fixef(pglmm.model)[[2]][2]*x + fixef(pglmm.model)[[2]][3]*x^2
upperquadrafun = function(x, c) fixef(pglmm.model)[[3]][1] + fixef(pglmm.model)[[3]][2]*x + fixef(pglmm.model)[[3]][3]*x^2
# The C.I. values
df.new = data.frame(coding_density = seq(0, 1, length.out = nquantiles+1))
df.new$lwr.pred = lowerquadrafun(df.new$coding_density)
df.new$upr.pred = upperquadrafun(df.new$coding_density)

# FIGURES
p = ggplot(data = meanrecrate_codingdensity, aes(x = coding_density, y = mean_rec)) +
  geom_point(colour = "Darkgrey", alpha = 0.5) +
  stat_smooth(aes(y = mean_rec), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Blue", fill = "Blue", alpha = 0.3) +
  stat_smooth(aes(y = mean_rec), method = "lm", formula = y ~ x, size = 1, colour = "Red", fill = "Red", alpha = 0.3) +
  scale_x_continuous(limits = c(0, 0.5)) +
  scale_y_continuous(limits = c(-5, 5)) +
  labs(x = "Coding density", y = "Standardized recombination rate") +
  # geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
  stat_function(fun = quadrafun, colour = "Black", size = 1) +
  geom_ribbon(data = df.new, aes(x = coding_density, ymin = lwr.pred, ymax = upr.pred), alpha = 0.4, inherit.aes = F, fill = "Black") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16))
p

#============================================================================#
# Figures -Rec rate ~ coding density quantiles, for each species, 2X2 plot ----
#============================================================================#

# FIGURES
# Divide the figure in 4 plots for visualization (44 species)
# Divide in four subsets
d1 = meanrecrate_codingdensity[which(meanrecrate_codingdensity$species %in% unique(meanrecrate_codingdensity$species)[1:11]),]
d2 = meanrecrate_codingdensity[which(meanrecrate_codingdensity$species %in% unique(meanrecrate_codingdensity$species)[12:22]),]
d3 = meanrecrate_codingdensity[which(meanrecrate_codingdensity$species %in% unique(meanrecrate_codingdensity$species)[23:33]),]
d4 = meanrecrate_codingdensity[which(meanrecrate_codingdensity$species %in% unique(meanrecrate_codingdensity$species)[34:44]),]
p1 = ggplot(data = d1, aes(x = coding_density, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  # scale_y_continuous(limits = c(0, 7)) +
  labs(x = "Coding density", y = "Recombination rate", colour = "Species") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16))
p2 = ggplot(data = d2, aes(x = coding_density, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  # scale_y_continuous(limits = c(0, 7)) +
  labs(x = "Coding density", y = "Recombination rate", colour = "Species") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16))
p3 = ggplot(data = d3, aes(x = coding_density, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  # scale_y_continuous(limits = c(0, 7)) +
  labs(x = "Coding density", y = "Recombination rate", colour = "Species") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16))
p4 = ggplot(data = d4, aes(x = coding_density, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  # scale_y_continuous(limits = c(0, 7)) +
  labs(x = "Coding density", y = "Recombination rate", colour = "Species") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16))
ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)





#============================================================================#
# DECOMPOSITION OF VARIANCE ----
#============================================================================#

# Analysis of variance of different effects: gene count, gene density, coding density

# Compare the % of variance explained by each effect


#============================================================================#
# ANALYSIS OF MULTICOLLINEARITY - Linear Model ----
#============================================================================#
# Collinearity diagnostics
# https://cran.r-project.org/web/packages/olsrr/vignettes/regression_diagnostics.html
# https://rdrr.io/cran/perturb/man/colldiag.html
# https://towardsdatascience.com/multi-collinearity-in-regression-fe7a2c1467ea

# install.packages("olsrr")
# library(olsrr)
# 
# # Import data for ALL chromosomes & species
# # Repoduce to all species in dataset
# # List of all maps
# list = system(paste("ls ", wd, "/data-cleaned/genome/gene_count/", sep = ""), intern = TRUE)
# readall = function(x) {
#   data = read.table(file = gzfile(paste("data-cleaned/genome/gene_count/", x, sep = "")),
#                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#   gendens = read.table(file = gzfile(paste("data-cleaned/genome/gene_density/", x, sep = "")),
#                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#   data$gene_density = gendens$gene_density
#   data$coding_density = gendens$coding_density
#   data$species = gsub("_", " ", rep(gsub("_[A-Za-z0-9]+_chromosome[A-Za-z0-9]*.txt.gz", "", x), nrow(data)))
#   data$species = gsub(" MaizeGDBConsensus", "", data$species)
#   return(data)
# }
# df = bind_rows(lapply(list, function(x) readall(x)))
# 
# # CHECK YOUR DATA FIRST
# # handle missing data previously to analyses
# sum(is.na(df$rec.rate))
# df = df[!is.na(df$rec.rate),]
# df = df[!is.na(df$gene_count),]
# df = df[!is.na(df$gene_density),]
# df = df[!is.na(df$coding_density),]
# 
# # EXPLANATORY VARIABLES
# hist(df$gene_count[df$gene_count < 30])
# sum(df$gene_count > 30)/length(df$gene_count)*100
# # 0.01% of very high values of gene count
# df = df[which(df$gene_count < 30),]
# hist(df$gene_count)
# 
# # Linear model for exploration
# lm.model = lm(log(rec.rate + 0.0001) ~ gene_count + gene_density + coding_density, data = df)
# summary(lm.model)
# plot(lm.model, which = c(1,2,3,5), ask = FALSE)
# 
# # Variance inflation factor, tolerance, eigenvalues and condition indices.
# # https://rdrr.io/cran/olsrr/man/ols_coll_diag.html
# coll_diag = ols_coll_diag(lm.model)
# coll_diag
# # Plot to detect non-linearity, influential observations and outliers.
# # Consists of side-by-side quantile plots of the centered fit and the residuals.
# # It shows how much variation in the data is explained by the fit and how much remains in the residuals.
# # For inappropriate models,
# # the spread of the residuals in such a plot is often greater than the spread of the centered fit.
# ols_plot_resid_fit_spread(model)
# 
# # Part & Partial correlations
# # Relative importance of independent variables in determining Y. How much each variable uniquely
# # contributes to R2 over and above that which can be accounted for by the other predictors.
# ols_correlations(model)

#============================================================================#
# GLOBAL MODEL - Build a global Linear Mixed model before decomposing variance ----
#============================================================================#
# Multicolinearity very strong in data
# Hence a new approach: multicolinearity analysis & decomposition of variance

# A GLMM model was built with all variables of interest for which we wanted to test partial effects on variance
# All three variables were highly correlated
# Hence the assumptions of the GLMM was poorly achieved (residuals)
# And model was a poor fit (predicted vs observed)

# BUT the full model was not intended for fit/predictions/inferences.
# It was a template full model for variables to compare in backward selection.

# Import data for ALL chromosomes & species
# Reproduce to all species in dataset
# List of all maps
list = system(paste("ls ", wd, "/data-cleaned/genome/gene_count/", sep = ""), intern = TRUE)
readall = function(x) {
  data = read.table(file = gzfile(paste("data-cleaned/genome/gene_count/", x, sep = "")),
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  gendens = read.table(file = gzfile(paste("data-cleaned/genome/gene_density/", x, sep = "")),
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  data$gene_density = gendens$gene_density
  data$coding_density = gendens$coding_density
  data$species = gsub("_", " ", rep(gsub("_[A-Za-z0-9]+_chromosome[A-Za-z0-9]*.txt.gz", "", x), nrow(data)))
  data$species = gsub(" MaizeGDBConsensus", "", data$species)
  return(data)
}
df = bind_rows(lapply(list, function(x) readall(x)))

# CHECK YOUR DATA FIRST
# handle missing data previously to analyses
sum(is.na(df$rec.rate))
df = df[!is.na(df$rec.rate),]
df = df[!is.na(df$gene_count),]
df = df[!is.na(df$gene_density),]
df = df[!is.na(df$coding_density),]

# EXPLANATORY VARIABLES
hist(df$gene_count[df$gene_count < 30])
sum(df$gene_count > 30)/length(df$gene_count)*100
# 0.01% of very high values of gene count
df = df[which(df$gene_count < 30),]
hist(df$gene_count)

# STANDARDISE PREDCITOR VARIABLES (CENTERED-REDUCED)
# Scale variables
# We are interested in comparing variance explained by each variable
# Hence centering-reducing to compare variables with different dimensions and units in the same model
# to compare their relative contributions ot the variance
df$scaled_gene_count = as.numeric(scale(df$gene_count))
df$scaled_gene_density = as.numeric(scale(df$gene_density))
df$scaled_coding_density = as.numeric(scale(df$coding_density))
# df$gene_density_divbygenecount = scale(df$gene_density_divbygenecount)
summary(df)

#----------------------------------------------------------------------------#
# Linear Mixed Model (Phylogenetic GLMM has too many individuals/variables, size too large)
# Decoupling null values and positive values of recombination
# Because they can be influenced by different processes (e.g. centromere suppression of recombination)
# And they can have very different sources of error and error rates
# i.e. null recombination is deeply impacted by the fit quality in regions of null recombination and can be less accurate

# Since there is an excess of zero values, I estimated a two-part model with
# 1/ Poisson component of 0/positive value (zero outcomes)
# zero-inflation and over-dispersion, see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3238139/
sum(df$rec.rate < 0.001)
sum(df$rec.rate < 0.001)/length(df$rec.rate)
# Consider very low recombination rates (i.e. < 0.001) as zero values
df$zero.rec = as.logical(df$rec.rate < 0.001)

# Regions of null recombination are genomically significantly different than regions with recombination
# Hence suggesting discriminating null recombination was relevant   
boxplot(df$gene_count~df$zero.rec)
boxplot(df$gene_density~df$zero.rec)
boxplot(df$coding_density~df$zero.rec)

# Test a binomial model for gene count
glmer.model.binomcomponent = glmer(zero.rec ~ gene_count + (1|species),
                                   data = df, family = binomial(link = logit))
summary(glmer.model.binomcomponent)
# Probability coefficient
plogis(fixef(glmer.model.binomcomponent))
# plot data
plot(x = df$gene_count, 
     y = df$zero.rec,
     main = "Probit Model for Gene count",
     xlab = "Gene count",
     ylab = "Probability",
     pch = 20,
     ylim = c(-0.4, 1.4),
     cex.main = 0.85)
# add horizontal dashed lines and text
abline(h = 1, lty = 2, col = "darkred")
abline(h = 0, lty = 2, col = "darkred")
# add estimated regression line
x = seq(0, 30, 1)
y = predict(glmer.model.binomcomponent, newdata = df, type = "response")
points(df$gene_count, y, lwd = 1.5, col = "steelblue")

# NOW RUN A FULL MODEL
glmer.model.binomcomponent = glmer(zero.rec ~ gene_count + gene_density + coding_density + (1|species),
                                   data = df, family = binomial(link = logit))
summary(glmer.model.binomcomponent)
# Probability coefficient, inverse of the logit
plogis(fixef(glmer.model.binomcomponent))

#----------------------------------------------------------------------------#
# 2/ Regression on non-zero values
# Consider very low recombination rates (i.e. < 0.001) as zero values
hist(log(df$rec.rate))
abline(v = log(0.001))
# New dataset of positive values
df_positive = df[which(df$rec.rate > 0.001),]

hist(df_positive$gene_count)
hist(df_positive$gene_density)
hist(df_positive$coding_density)
hist(df_positive$scaled_gene_count)
hist(df_positive$scaled_gene_density)
hist(df_positive$scaled_coding_density)
# Skewed distribution
# Log-transform

# DEPENDENT VARIABLE
hist(log(df_positive$rec.rate))
# Slightly skewed distribution toward short recombination rates

# LMER Model
lmm.model.positivecomponent = lmer(log(rec.rate) ~ gene_count + gene_density + coding_density + (1|species),
                                   data = df_positive)
# lmm.model.positivecomponent = glmer(rec.rate ~ gene_count + gene_density + coding_density + (1|species), data = df_positive,
#                                     family = Gamma(link = log))
summary(lmm.model.positivecomponent)
# Diagnostic plots
# Residuals distribution
par(mfrow = c(2,2))
qqPlot(residuals(lmm.model.positivecomponent), ylab = "Residuals")
# Residuals as a function of fitted
plot(fitted(lmm.model.positivecomponent), residuals(lmm.model.positivecomponent), xlab = "Fitted", ylab = "Residuals")
abline(h = 0)
# Observed~predicted qqplot
plot(df$rec.rate, predict(lmm.model.positivecomponent, newdata = df), xlab = "Observed", ylab = "Predicted")
abline(a = 0, b = 1)
par(mfrow = c(1,1))

# The full-set model has moderate fit and poor predictions,
# with slight violation of assumptions (normality of residuals)
# Yet it is only a template model to test for decomposition of variance

#============================================================================#
# DECOMPOSITION OF THE VARIANCE WITH PARTR2 ----
#============================================================================#
# In order to investigate structure of the variance and collinearity among independent variables

# library(MuMIn)
# install.packages("remotes")
# remotes::install_github("mastoffel/partR2")
library(partR2)
# ?partR2

#----------------------------------------------------------------------------#
# Binomial component
# Very long, even for 100 bootstraps (~1h for nboot = 10), "vecteurs de mémoire épuisés" when using parallel = TRUE
# R2_binom = partR2(glmer.model.binomcomponent, partvars = c("gene_count", "gene_density",
#                                                            "coding_density"), data = df,
#                   nboot = 100, CI = 0.95)
# summary(R2_binom)
# Save R2 decomposition of variance
# save(R2_binom, file = "output/models/glmer.model.binomcomponent_genedensity_R2.Rda")

load(file = "output/models/glmer.model.binomcomponent_genedensity_R2.Rda")
summary(R2_binom)
# Part R2 is the variance explained by each individual predictor, by iterative removal of each predictor
forestplot(R2_binom, type = "R2", line_size = 0.7, text_size = 14, point_size = 3)
# Check also inclusive R2 (gives information about the structure of variance explained; square of the structure coefficient, i.e. the correlation between predictors)
# and beta weights (standardized regression slopes)
forestplot(R2_binom, type = "IR2", line_size = 0.7, text_size = 14, point_size = 3)
forestplot(R2_binom, type = "BW", line_size = 0.7, text_size = 14, point_size = 3)
# Structure coefficients show the correlation between individual predictor variables
forestplot(R2_binom, type = "SC", line_size = 0.7, text_size = 14, point_size = 3)


#----------------------------------------------------------------------------#
# Positive component
# R2_positive = partR2(lmm.model.positivecomponent,
#                      partvars = c("gene_count", "gene_density", "coding_density"), data = df_positive,
#                      nboot = 1000, CI = 0.95)
# summary(R2_positive)
# Save R2 decomposition of variance
# save(R2_positive, file = "output/models/lmm.model.positivecomponent_genedensity_R2.Rda")

load(file = "output/models/lmm.model.positivecomponent_genedensity_R2.Rda")
summary(R2_positive)# Part R2 is the variance explained by each individual predictor, by iterative removal of each predictor
forestplot(R2_positive, type = "R2", line_size = 0.7, text_size = 14, point_size = 3)
# Check also inclusive R2 (gives information about the structure of variance explained; square of the structure coefficient, i.e. the correlation between predictors)
# and beta weights (standardized regression slopes)
forestplot(R2_positive, type = "IR2", line_size = 0.7, text_size = 14, point_size = 3)
forestplot(R2_positive, type = "BW", line_size = 0.7, text_size = 14, point_size = 3)
# Structure coefficients show the correlation between individual predictor variables
forestplot(R2_positive, type = "SC", line_size = 0.7, text_size = 14, point_size = 3)


#============================================================================#
# IS THERE AN EFFECT OF GENE LENGTH? ----
#============================================================================#

# Gene length were estimated in the 'Genome_architecture.R' file
df = read.table(file = gzfile("output/gene_length.txt.gz"), sep = "\t", header = TRUE)

# Now we explore if on average longer genes have a higher recombination rate.
# We take length of genes by start-end coordinates in the gff file
# And the estimated recombination rate associated to the genomic window where the gene lies (start position)

# Estimate 20 quantiles of gene length per species, and compute the mean recombination rate.

list_species = unique(df$species)
df_genelength = data.frame(species = rep(list_species, each = 20), gene_length = NA, rec_rate = NA)
for (l in list_species) {
  cat(l, "\n")
  # Compute gene length quantiles
  quantiles = quantile(df$length[which(df$species == l)], probs = seq(0, 1, length.out = 21), na.rm = TRUE)
  mean_length = c()
  mean_recrate = c()
  for (i in 1:20) {
    mean_length[i] = mean(df$length[which(df$species == l & df$length > quantiles[i] & df$length < quantiles[i+1])], na.rm = TRUE)
    mean_recrate[i] = mean(df$recrate[which(df$species == l & df$length > quantiles[i] & df$length < quantiles[i+1])], na.rm = TRUE)
  }
  df_genelength$gene_length[which(df_genelength$species == l)] = mean_length
  df_genelength$rec_rate[which(df_genelength$species == l)] = mean_recrate
}
# Save results
write.table(df_genelength, "output/recrate~genelength.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)
df_genelength = read.table("output/recrate~genelength.txt", header = TRUE)
# Quick visual exploration
hist(df_genelength$gene_length)
hist(log(df_genelength$gene_length))
hist(df_genelength$rec_rate)
plot(log(df_genelength$gene_length),
     df_genelength$rec_rate,
     xlab = "Gene length (log)", ylab = "Recombination rate")

# One line per species
list_species = unique(df_genelength$species)
plot(log(df_genelength$gene_length[which(df_genelength$species == list_species[1])]),
     df_genelength$rec_rate[which(df_genelength$species == list_species[1])],
     xlab = "Gene length (log)", ylab = "Recombination rate", type = "l",
     xlim = c(4, 12), ylim = c(0,20))
for (l in 2:length(list_species)) {
  points(log(df_genelength$gene_length[which(df_genelength$species == list_species[l])]),
         df_genelength$rec_rate[which(df_genelength$species == list_species[l])], type = "l")
}

# It seems that there is no effect of gene length on recombination rates.

############################################################################################# #
# DEPRECATED BELOW THIS LINE ----
############################################################################################# #










































# #============================================================================#
# # Genetic distance (cM) expressed as a function of gene distance (cumulative number of genes) ----
# #============================================================================#
# # Enhance data with a new species
list_species = "Eucalyptus_grandis"
list_accessions = "Egrandis_297_v2.0"
# # Process each genome dataset and save in separate files per species
# # Gene position
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

    # Clean genpos for scaffolds and NA chromosomes
    genpos = subset(genpos, !grepl("scaffold", genpos$chromosome))
    
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

# # Exon position
# for (acc in 1:length(list_accessions)) { # LOOP ON ACCESSIONS
#   # Getting infos on dataset to treat
#   accession = list_accessions[acc]
#   species = list_species[acc]
#   cat(species, accession, "\n")
#   # Retrieving feature position, i.e. genes
#   exonpos = exon.position(species, accession)
#   # Some data is not available or no features have been reported in annotation files (bad annotation)
#   if (nrow(exonpos) > 0) {
#     # if (((exonpos) != "Accession is not available in your dataset.")[1]) {
#     # Add species and accession to the results
#     exonpos = cbind(data.frame(species = rep(species, nrow(exonpos)), accession = rep(accession, nrow(exonpos))),
#                     exonpos)
#     #Load the dataset of the species if it exists, otherwise create a new one
#     if (file.exists(paste("data-cleaned/genome/exon_positions/", species,"_exon_positions.txt.gz", sep = ""))) {
#       # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
#       dataset = read.table(gzfile(paste("data-cleaned/genome/exon_positions/", species,"_exon_positions.txt.gz", sep = "")),
#                            header = TRUE, sep = "\t")
#     } else {
#       # Template of a new empty dataset
#       dataset = data.frame(species = character(0), accession = character(0), id = character(0),
#                            type = character(0), chromosome = character(0),
#                            strand = character(0), rank = numeric(0), start = numeric(0), end = numeric(0))
#       write.table(dataset, gzfile(paste("data-cleaned/genome/exon_positions/", species,"_exon_positions.txt.gz", sep = "")), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)      # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
#       # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
#       dataset = read.table(gzfile(paste("data-cleaned/genome/exon_positions/", species,"_exon_positions.txt.gz", sep = "")),
#                            header = TRUE, sep = "\t")
#     }
#     
#     # Remove older rows in the dataset
#     dataset = dataset[!(dataset$accession == accession),]
#     # Bind new rows
#     dataset = rbind(dataset, exonpos)
#     # Init an empty dataset
#     # dataset = data.frame(species = character(0), accession = character(0), id = character(0),
#     #                      type = character(0), chromosome = character(0), strand = character(0),
#     #                      start = character(0), end = character(0)) 
#     # Write updated dataset in 'data-cleaned'
#     write.table(dataset, gzfile(paste("data-cleaned/genome/exon_positions/", species,"_exon_positions.txt.gz", sep = "")), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)
#     rm(dataset)
#     # } else {
#     #   warning("Accession is not available in your dataset.")
#     # }
#   } else {
#     warning(paste(species, accession,"No feature found in data.", sep = " "))
#   }
#   # END OF ACCESSION LOOP
# }
# # CDS position
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
#     # Add species and accession to the results
#     cdspos = cbind(data.frame(species = rep(species, nrow(cdspos)), accession = rep(accession, nrow(cdspos))),
#                    cdspos)
#     #Load the dataset of the species if it exists, otherwise create a new one
#     if (file.exists(paste("data-cleaned/genome/CDS_positions/", species,"_CDS_positions.txt.gz", sep = ""))) {
#       # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
#       dataset = read.table(gzfile(paste("data-cleaned/genome/CDS_positions/", species,"_CDS_positions.txt.gz", sep = "")),
#                            header = TRUE, sep = "\t")
#     } else {
#       # Template of a new empty dataset
#       dataset = data.frame(species = character(0), accession = character(0), id = character(0),
#                            type = character(0), chromosome = character(0),
#                            strand = character(0), phase = numeric(0), start = numeric(0), end = numeric(0))
#       write.table(dataset, gzfile(paste("data-cleaned/genome/CDS_positions/", species,"_CDS_positions.txt.gz", sep = "")), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)      # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
#       # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
#       dataset = read.table(gzfile(paste("data-cleaned/genome/CDS_positions/", species,"_CDS_positions.txt.gz", sep = "")),
#                            header = TRUE, sep = "\t")
#     }
#     
#     # Remove older rows in the dataset
#     dataset = dataset[!(dataset$accession == accession),]
#     # Bind new rows
#     dataset = rbind(dataset, cdspos)
#     # Init an empty dataset
#     # dataset = data.frame(species = character(0), accession = character(0), id = character(0),
#     #                      type = character(0), chromosome = character(0), strand = character(0),
#     #                      start = character(0), end = character(0)) 
#     # Write updated dataset in 'data-cleaned'
#     write.table(dataset, gzfile(paste("data-cleaned/genome/CDS_positions/", species,"_CDS_positions.txt.gz", sep = "")), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)
#     rm(dataset)
#     # } else {
#     #   warning("Accession is not available in your dataset.")
#     # }
#   } else {
#     warning(paste(species, accession,"No feature found in data.", sep = " "))
#   }
#   # END OF ACCESSION LOOP
# }
# 
# # Then compute a new map
# map = as.character(metadata.clean$id[38])
# gc.map(map, ex = TRUE, method = 'gene')
# # gcgenes = read.table(file = gzfile(paste("data-cleaned/genome/gc_genes/", map, ".txt.gz", sep = "")), header = TRUE, sep = "\t")
# gene_distance.map(set = map, scale = 0.1)
# 
# 
# 
# 
# 
# 
# 
# 
# # Test the method on A. thaliana
# # Retrieve gene distances (cumulative number of genes) in a given Marey map
# source("sources/Genome.R") 
# # set = "Arabidopsis_thaliana_Serin2017"
# # set = "Zea_mays_MaizeGDBConsensus_v4"
# # set ="Triticum_aestivum_GutierrezGonzalez2019"
# # set = "Solanum_lycopersicum_Gonda2018"
# # set = "Phaseolus_vulgaris_Song2015"
# # A function that 
# gene_distance.map(set = map, scale = 0.1)
# 
# 
# # Reproduce to all species in dataset
# # List of all maps (all maps with a gc_genes file giving the list of genes and their position in a reference genome)
# list = system(paste("ls ", wd, "/data-cleaned/genome/gc_genes/", sep = ""), intern = TRUE)
# list = list[grep("*.txt.gz", list)]
# list = gsub(".txt.gz", "", list)
# list
# source("sources/Genome.R")
# for (sp in list) {
#   gene_distance.map(set = sp, scale = 0.1)
# }

# 
# 
# # Load the map with gene distances and genetic distances (cM)
# list
# set = list[11]
# chromosome = "1"
# gene_distance_map = read.table(file = gzfile(paste("data-cleaned/genome/gene_distance/", set, ".txt.gz", sep = "")), header = TRUE, sep = "\t")
# 
# gene_distance_map = gene_distance_map[gene_distance_map$map == chromosome,]
# 
# # Figure
# plot(gene_distance_map$gene_distance, gene_distance_map$gen)
# 
# 
# # Figure
# par(mfrow = c(1,2))
# plot(gene_distance_map$phys/1000000, gene_distance_map$gen, xlab = "Physical distance (Mb)", ylab = "Genetic distance (cM)")
# lines(x = c(min(gene_distance_map$phys/1000000), max(gene_distance_map$phys/1000000)), y = c(min(gene_distance_map$gen), max(gene_distance_map$gen)))
# plot(gene_distance_map$gene_distance, gene_distance_map$gen, xlab = "Gene distance (cumulative number of genes)", ylab = "Genetic distance (cM)")
# lines(x = c(min(gene_distance_map$gene_distance), max(gene_distance_map$gene_distance)), y = c(min(gene_distance_map$gen), max(gene_distance_map$gen)))
# par(mfrow = c(1,1))
# 
# # Interestingly,
# plot(gene_distance_map$phys/1000000, gene_distance_map$gene_distance, xlab = "Physical distance (Mb)", ylab = "Genetic distance (cM)")



#============================================================================#
# Mixed model to explain recombination rates as a function of gene density
#============================================================================#

#----------------------------------------------------------------------------#
# Prepare dataset
#----------------------------------------------------------------------------#
# List of all maps
list = system(paste("ls ", wd, "/data-cleaned/genome/gene_density/", sep = ""), intern = TRUE)
list = list[grep("*.txt.gz", list)]
list = gsub(".txt.gz", "", list)
list
# Import all gene density maps in a single data frame
tmp = read.table(file = gzfile(paste("data-cleaned/genome/gene_density/", list[1], ".txt.gz", sep = "")),
                                            header = TRUE, sep = "\t")
df = cbind(data.frame(set = rep(list, nrow(tmp))), tmp)
for (i in 2:length(list)) {
  tmp = read.table(file = gzfile(paste("data-cleaned/genome/gene_density/", list[i], ".txt.gz", sep = "")),
                   header = TRUE, sep = "\t")
  tmp = cbind(data.frame(set = rep(list[i], nrow(tmp))), tmp)
  df = rbind(df, tmp)
}
rm(tmp)
# Get the chromosome and species names
df$species = gsub("_", " ", (gsub("_[0-9A-Za-z]+_chromosome[0-9A-Za-z]+$", "", df$set)))
df$chromosome = gsub("[0-9A-Za-z_]+_chromosome", "", df$set)
df$species = gsub(" MaizeGDBConsensus", "", df$species)
df$species = gsub(" IBM", "", df$species)
# Save and re-upload data
write.table(df, gzfile("data-cleaned/genome/gene_density_all.txt.gz"), row.names = FALSE, col.names = TRUE)


#----------------------------------------------------------------------------#
# Load dataset
#----------------------------------------------------------------------------#
df = read.table(gzfile("data-cleaned/genome/gene_density_all.txt.gz"), header = TRUE)

#----------------------------------------------------------------------------#
#Models
#----------------------------------------------------------------------------#
# Plot the species-specific mean recombination rate for each gene count
# Chromosomes pooled
# i.e. one point per species and gene count (up to 30 genes)
summary(df$gene_count)
maxgene = 30
nboot = 1000 # Number of bootstrap iterations
# Create a new data frame
list_species = unique(df$species)
list_species
plotdf = data.frame(species = rep(list_species, each = (maxgene+1)), gene_count = rep(0:maxgene, length(list_species)), mean_recrate = NA,
                    lowerCI = NA, upperCI = NA)
# Compute the mean recombination rate for each class
# TODO compute by bootstrap the mean and CI
for (i in 1:nrow(plotdf)) {
  cat("Progress: ", i/nrow(plotdf)*100, "%\n")
  boot = numeric(nboot)
  for (j in 1:nboot) {
    boot[j] = mean(sample(df$rec.rate[which((df$gene_count == plotdf$gene_count[i]) & (df$species == as.character(plotdf$species[i])))], replace = TRUE), na.rm = TRUE)
  }
  plotdf$mean_recrate[i] = mean(boot, na.rm = TRUE)
  plotdf$lowerCI[i] = quantile(boot, 0.025, na.rm = TRUE)
  plotdf$upperCI[i] = quantile(boot, 0.975, na.rm = TRUE)
}

# A color per species
color.species = data.frame(species = list_species, color = as.character(polychrome(length(list_species))))
# show_col(as.character(polychrome(length(unique(species)))))
plotdf = merge(plotdf, color.species, by = "species")

# Save results
write.table(plotdf, gzfile("output/figures/Rec~Gene_count.txt.gz"), row.names = FALSE, col.names = TRUE)
plotdf = read.table(gzfile("output/figures/Rec~Gene_count.txt.gz"), header = TRUE)

# Plot it
# plot(plotdf$gene_count, plotdf$mean_recrate)
plot(plotdf$gene_count, plotdf$mean_recrate, xlab = "Number of genes (100kb windows)", ylab = "Recombination rate",
     bg = as.character(plotdf$color), pch = 21)
segments(plotdf$gene_count, plotdf$lowerCI, plotdf$gene_count, plotdf$upperCI, col = as.character(plotdf$color))



# TODO Compute the recombination rate for each gene count (all species pooled)
# And see if there is a differences with species specific patterns

#----------------------------------------------------------------------------#
# Mixed model with Linear Regression ----
#----------------------------------------------------------------------------#
# Rec ~ Gene density + (1|Species)
# With lme4 implementation
# lmer.genedensity = lmer(rec.rate ~ gene_count + (1|species), data = df)
# summary(lmer.genedensity)

# With nlme implementation
df$species = as.character
lme.genedensity = lme(fixed = rec.rate ~ gene_count, random = ~ 1|species, data = df, na.action = na.omit)
summary(lme.genedensity)
plot(lme.genedensity) # Residuals of the model
# "Remember, for a well fitting regression, we want the plot of our residuals to meet the following criteria: (1) they’re pretty symmetrically distributed (2) they’re relatively small and (3) they don’t follow a clear pattern"
# Source: https://crumplab.github.io/psyc7709/book/docs/a-tutorial-for-using-the-lme-function-from-the-nlme-package-.html
# Our results does not follow these prerequisites; model not validated
# Proportion of variance explained by the random factor
(0.8723962/(0.8723962+1.433164))
# Log likelihood
lme.genedensity$logLik
# 37% of variance explained by the random factor
lme.genedensity$coefficients
sumEffects = summary(lme.genedensity)
sumEffects$tTable
# Grahic representation of the model
plot(plotdf$gene_count, plotdf$mean_recrate, bg = as.character(plotdf$color), pch = 21)
abline(sumEffects$tTable[1], sumEffects$tTable[2])
# Yet, the effect does not seem to be linear (at least for all gene counts)
# And model is not validated

# However, a linear relationship can be imagined for gene count <10
lme.genedensity = lme(fixed = rec.rate ~ gene_count, random = ~ 1|species, data = df[which(df$gene_count < 11),], na.action = na.omit)
summary(lme.genedensity)
plot(lme.genedensity) # Residuals of the model
plot(plotdf$gene_count, plotdf$mean_recrate, bg = as.character(plotdf$color), pch = 21, xlim = c(0, 10))
abline(sumEffects$tTable[1], sumEffects$tTable[2])

#----------------------------------------------------------------------------#
#  Non Linear Model with Random Effect ----
#----------------------------------------------------------------------------#
# Which model will we use?

# Logistic regression
# Seems a reasonable choice, given the graphical representation of data
# use nlme instead of lme
nlme.genedensity = nlme(model = rec.rate ~ exp(gene_count)/(exp(gene_count)+1),
                        fixed = rec.rate ~ gene_count, random = rec.rate ~ 1|species, data = df,
                        start = c(mean(df$gene_count)), na.action = na.omit)
# Use glmer in lme4 for a Mixed effects logistic regression
# Tutorial at https://stats.idre.ucla.edu/r/dae/mixed-effects-logistic-regression/
# https://stats.idre.ucla.edu/other/mult-pkg/introduction-to-generalized-linear-mixed-models/
glmer.genedensity = glmer(rec.rate ~ gene_count + (1 | species), data = df, family = poisson, na.action = na.omit)
summary(glmer.genedensity)
plot(glmer.genedensity)





#============================================================================#
# Differences in distribution of recombination between physical distances and gene distances, ----
# explained by a coefficient of departure from evenness
#============================================================================#

# Gini's coefficient

# Is evenness of the distribution different in physical distances and gene distances?

# I need Marey maps with physical distances and gene distances and an evenly distribution of markers
# i.e. the map computed for t_bar, with an interpolated grid of 1,000 markers evenly distributed along the genome
# Compute the same grid for gene distances
# Data is in "data-cleaned/genome/gene_distance/" for physical & gene distances ~ genetic distances

# TODO Even sampling in a grid of 1,000 markers by interpolating a curve

# Otherwise, I could use the recombination maps with physical distances and gene distances
# And measure the evenness of the distribution with the coefficient of variation of squared deviation to the linear relationship


set = "Zea_mays_MaizeGDBConsensus_v4"
# Load the map with gene distances and genetic distances (cM)
chromosome = "1"
gene_distance_map = read.table(file = gzfile(paste(wd, "data-cleaned/genome/gene_distance/", set, ".txt.gz", sep = "")), header = TRUE, sep = "\t")
gene_distance_map = gene_distance_map[gene_distance_map$map == chromosome,]

par(mfrow = c(1,2))
plot(gene_distance_map$phys/1000000, gene_distance_map$gen, xlab = "Physical distance (Mb)", ylab = "Genetic distance (cM)")
lines(x = c(min(gene_distance_map$phys/1000000), max(gene_distance_map$phys/1000000)), y = c(min(gene_distance_map$gen), max(gene_distance_map$gen)))
plot(gene_distance_map$gene_distance, gene_distance_map$gen, xlab = "Gene distance (cumulative number of genes)", ylab = "Genetic distance (cM)")
lines(x = c(min(gene_distance_map$gene_distance), max(gene_distance_map$gene_distance)), y = c(min(gene_distance_map$gen), max(gene_distance_map$gen)))
par(mfrow = c(1,1))

# Compute Gini's index on the observed distribution (beware, sampling not even, but same individuals for gene or physical)
# x is a vector of numerical values of interest for the chromosome (i.e. )
Gini = function(x) {
  # Numeric values centered-reduced to avoid scale problems
  x = scale(x)
  # Mean difference of Gini, mean of all absolute deviations for all pairs of individuals
  # (expected deviation between values of two individuals randomly sampled)
  # The matrix of all pairs
  m = matrix(NA, ncol = length(x), nrow = length(x))
  for (i in 1:nrow(m)) {
    for (j in 1:ncol(m)) {
      m[i,j] = as.numeric(x[1])*as.numeric(x[3])
    }
  }
  E = sum(m, na.rm = TRUE)/(length(x)^2)
  # The mean of all x values
  M = mean(x, na.rm = TRUE)
  # The Gini's index is
  E/2*M
}
Gini(gene_distance_map$phys)
Gini(gene_distance_map$gene)
# Batch Gini's index of all mpas and chromosomes with gene distance data
# List of all maps (all maps with a gc_genes file giving the list of genes and their position in a reference genome)
list = system(paste("ls ", wd, "/data-cleaned/genome/gc_genes/", sep = ""), intern = TRUE)
list = list[grep("*.txt.gz", list)]
list = gsub(".txt.gz", "", list)
list
df = data.frame(set = character(0), map = character(0), Gini.phys = numeric(0), Gini.gene = numeric(0))
for (set in list) {
  cat("Processing", set, "\n")
  gene_distance_map = read.table(file = gzfile(paste(wd, "data-cleaned/genome/gene_distance/", set, ".txt.gz", sep = "")), header = TRUE, sep = "\t")
  # List of chromosomes
  chromosomes = unique(gene_distance_map$map)
  # Process each chromosome
  for (chr in chromosomes) {
    gene_distance_chr = gene_distance_map[gene_distance_map$map == chr,]
    Gini.phys = Gini(gene_distance_chr$phys)
    Gini.gene = Gini(gene_distance_chr$gene)
    df = rbind(df, data.frame(set, chr, Gini.phys, Gini.gene))
  }
}

# Save
write.table(df, file = paste("output/gene_distance/Gini.txt", sep = ""),
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
boxplot(df$Gini.phys, df$Gini.gene)

# Gini's index in this implementation did not showed differences between physical distances and gene distance
# Use another metric for the deviation from the null model (uniform distribution)
# Ideally, the metric must be scale-free, to compare physical and gene distances using very different scales and units
# Instead, use the Mean Absolute Error (MAE, https://en.wikipedia.org/wiki/Mean_absolute_error)

# First step, model the null hypothesis given a chromosome
# The line between the first and last markers
# Take the vector of X in argument
uniform_model = function() {
  
}

# Compute the MAE between the empirical distribution and the theoretical null model
# One value per marker
# Take the vectors of X and Y
# Y is constant , X is a scaled value being physical or gene distance




#============================================================================#
# Gene landscapes correlate to recombination landscapes? ----
#============================================================================#

#============================================================================#
#  Estimate mean number of genes ~ relative distance ----
#============================================================================#
chromosome.stats = read.table(file = "tables/chromosome.stats.csv", header = TRUE, sep = ";")
df_dist2telomere = data.frame(set = character(0), chromosome = character(0), gene_count = numeric(0), dist2telomere = numeric(0))

# After determining the center of the chromosome, the distance to the telomere is simply physical position for the first half
# and absolute(max - physical position) for the other half
# In other words, it is also the minimal value of distances taken from both tips of the chromosome
# YET this means we could do an averaging of two contrasted patterns, i.e. it can erase a true signal
# Hence we needed to consider both halves separately, but they were not independent
# Only one half for each chromosome was randomly sampled to avoid pseudo-replicates
list_maps = data.frame(set = as.character(chromosome.stats$set[which(!is.na(chromosome.stats$genecount))]), chromosome = as.character(chromosome.stats$chromosome[which(!is.na(chromosome.stats$genecount))]), stringsAsFactors = FALSE)

for(i in 1:nrow(list_maps)) {
  cat(list_maps[i,1], list_maps[i,2], "\n")
  gene_map = read.table(paste("data-cleaned/genome/gene_count/", list_maps[i,1], "_chromosome", list_maps[i,2], ".txt.gz", sep = ""), header = TRUE)
  gene_map = gene_map[order(gene_map$phys),]
  totalchrsize = max(gene_map$phys, na.rm = TRUE)
  # Randomly subset half of the chromosome
  side = sample(c("left", "right"), size = 1)
  if (side == "left") {
    gene_map = gene_map[which(gene_map$phys < totalchrsize/2),]
  } else {
    gene_map = gene_map[which(gene_map$phys > totalchrsize/2),]
  }
  # Estimate distance to telomere
  gene_map$dist2telomere = apply(data.frame(gene_map$phys, abs(gene_map$phys-totalchrsize)), 1, min)
  # Results
  df = data.frame(set = rep(list_maps[i,1], nrow(gene_map)), chromosome = rep(list_maps[i,2], nrow(gene_map)),
                  gene_count = gene_map$gene_count, dist2telomere = gene_map$dist2telomere)
  # remove older estimates
  df_dist2telomere = df_dist2telomere[!(df_dist2telomere$set == list_maps[i,1] & df_dist2telomere$chromosome == list_maps[i,2]),]
  # Append data
  df_dist2telomere = rbind(df_dist2telomere, df)
}
df_dist2telomere

write.table(df_dist2telomere, "output/gene_count/AllChromosomes_Dist2telomeres.txt", row.names = FALSE, col.names = TRUE)

#============================================================================#
#   Raw data - Gene count ~ distance to nearest chromosome end ----
#============================================================================#
df_dist2telomere = read.table("output/gene_count/AllChromosomes_Dist2telomeres.txt", header = TRUE)

# Remove NA and transform null values
dist2telomere = subset(df_dist2telomere, (!is.na(df_dist2telomere$gene_count)))
# dist2telomere$rec.rate = dist2telomere$rec.rate + 0.0001
# dist2telomere$dist2telomere = dist2telomere$dist2telomere + 0.0001
# Remove null values
dist2telomere = dist2telomere[(dist2telomere$gene_count != 0),]
dist2telomere = dist2telomere[(dist2telomere$dist2telomere != 0),]
# Species names
dist2telomere = merge(dist2telomere, unique(chromosome.stats[,c(1, 24)]), by = "set", all.x = TRUE)

# Distribution of data
hist(dist2telomere$dist2telomere, breaks = 20, xlab = "Distance to telomere (Mb)", main = "")
hist(log(dist2telomere$dist2telomere), breaks = 20, xlab = "Distance to telomere (log-transformed)", main = "")

hist(dist2telomere$gene_count, breaks = 20, xlab = "Gene count", main = "")
hist(log(dist2telomere$gene_count), breaks = 20, xlab = "Gene count (log-transformed)", main = "")


#============================================================================#
#   Rec. rate ~ relative distance to nearest chromosome end (quantiles of distances) ----
#============================================================================#

df_dist2telomere = read.table("output/gene_count/AllChromosomes_Dist2telomeres.txt", header = TRUE)
# Compute the quantiles of distances
# Prepare the dataset with species and chromosomes to gather
nquantiles = 20
df = data.frame(set = character(0), chromosome = character(0),
                gene_count = numeric(0), lower = numeric(0), upper = numeric(0),
                quantile.lower = numeric(0), quantile.upper = numeric(0), dist2telomere = numeric(0))
# Compute the quantiles of distances for each species and chromosome
for (i in 1:nrow(list_maps)) {
  cat(list_maps[i,1], list_maps[i,2], "\n")
  df_subset = df_dist2telomere[which(df_dist2telomere$set == list_maps[i,1] & df_dist2telomere$chromosome == list_maps[i,2]),]
  # Boundaries of 20 bins of equal distance
  # quantiles = seq(df_subset$dist2telomere, probs = seq(0, 1, length.out = nquantiles+1))
  bins = seq(from = 0, to = max(df_dist2telomere$dist2telomere[which(df_dist2telomere$set == list_maps[i,1] & df_dist2telomere$chromosome == list_maps[i,2])], na.rm = TRUE), length.out = nquantiles+1)
  quantile.lower = bins[1:nquantiles]
  quantile.upper = bins[2:(nquantiles+1)]
  newdata = data.frame(set = rep(as.character(unique(df_subset$set)), each = nquantiles), chromosome = rep(as.character(unique(df_subset$chromosome)), each = nquantiles),
                       gene_count = NA, lower = NA, upper = NA,
                       quantile.lower = quantile.lower, quantile.upper = quantile.upper, dist2telomere = (quantile.upper+quantile.lower)/2)
  df = rbind(df,
             newdata)
  # # Process each chromosome
  # list_chromosomes = unique(df_sp$chromosome)
  # for (chr in list_chromosomes) {
  #   # Boundaries of 20 quantiles
  #   quantiles = quantile(df_dist2telomere$dist2telomere[df_dist2telomere$set == set & df_dist2telomere$chromosome == chr], probs = seq(0, 1, length.out = 21))
  #   quantile.lower = quantiles[1:20]
  #   quantile.upper = quantiles[2:21]
  #   df = rbind(df, data.frame(set = rep(set, 20), chromosome = rep(chr, 20), rec.rate = NA, upper = NA, lower = NA, quantile.lower = quantile.lower, quantile.upper = quantile.upper))
  # }
  rm(df_subset)
}
df

#----------------------------------------------------------------------------#
# Compute the mean gene count in distance intervals
# Take a vector of indexes in the df dataframe
# and return a vector of gene count
genecount_in_distance_quantile = function(idx) {
  genecount.quantile = mean(df_dist2telomere$gene_count[df_dist2telomere$set == df$set[idx]
                                                     & df_dist2telomere$chromosome == df$chromosome[idx]
                                                     & df_dist2telomere$dist2telomere >= df$quantile.lower[idx]
                                                     & df_dist2telomere$dist2telomere <= df$quantile.upper[idx]], na.rm = TRUE)
  return(genecount.quantile)
}

idx = c(1:nrow(df))
require(pbmcapply)
gene_count = unlist(pbmclapply(idx, function(x) genecount_in_distance_quantile(x), mc.cores =  7))
df$gene_count = gene_count
hist(df$gene_count)


#----------------------------------------------------------------------------#
# Transform in relative distances
# relatives_distances = function(idx) {
#   reldist = data.frame(quantile.lower = df$quantile.lower[idx], quantile.upper = df$quantile.upper[idx], dist2telomere = df$dist2telomere[idx])
#   maxdist = max(df_dist2telomere$dist2telomere[df_dist2telomere$set == df$set[idx] & df_dist2telomere$chromosome == df$chromosome[idx]], na.rm = TRUE)
#   reldist$quantile.lower = reldist$quantile.lower/maxdist
#   reldist$quantile.upper = reldist$quantile.upper/maxdist
#   reldist$dist2telomere = reldist$dist2telomere/maxdist
#   return(reldist)
# }
# idx = c(1:nrow(df))
# require(pbmcapply)
# reldist = pbmclapply(idx, function(x) relatives_distances(x), mc.cores =  7)
# reldist = data.frame(matrix(unlist(reldist), nrow=length(reldist), byrow=T))
# reldist
# df$quantile.lower = reldist$X1
# df$quantile.upper = reldist$X2
# df$dist2telomere = reldist$X3
df$quantile.lower = rep(seq(0, 0.45, length.out = 20), nrow(df)/20)
df$quantile.upper = rep(seq(0.05, 0.5, length.out = 20), nrow(df)/20)
df$dist2telomere = rep(seq(0.025, 0.475, length.out = 20), nrow(df)/20)

#----------------------------------------------------------------------------#
write.table(df, "output/gene_count/DistancesRelative_bins.txt", row.names = FALSE, col.names = TRUE)

df = read.table("output/gene_count/DistancesRelative_bins.txt", header = TRUE)

#----------------------------------------------------------------------------#
# Scale recombination rates
for (i in 1:nrow(list_maps)) {
  cat(list_maps$set[i], list_maps$chromosome[i], "\n")
  df$gene_count[which(df$set == list_maps$set[i] & df$chromosome == list_maps$chromosome[i])] =
    scale(df$gene_count[which(df$set == list_maps$set[i] & df$chromosome == list_maps$chromosome[i])])
}
hist(df$gene_count)
write.table(df, "output/gene_count/DistancesRelativeScaled_bins.txt", row.names = FALSE, col.names = TRUE)

df = read.table("output/gene_count/DistancesRelativeScaled_bins.txt", header = TRUE)




#============================================================================#
# END ----
#============================================================================#