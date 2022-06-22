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
# - compute recombination rates as a function of the distance to telomere/centromere
# - infer the best model explaining distribution of the recombination (i.e. centromere or telomere driven)

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
library(rBLAST)
library(Biostrings)
library(pals)
library(egg)
library(lme4)
library(nlme)
library(lmerTest)
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
chromosome.stats = read.table(file = "tables/chromosome.stats.csv", header = TRUE, sep = ";")

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


#============================================================================#
#
#   TELOMERES: IS RECOMBINATION DRIVEN BY TELOMERES? ----
#
#============================================================================#

# The hypothesis is that recombination is a process beginning in telomeres and decreasing as physical distance increases


#============================================================================#
#  Estimate Rec. rate ~ relative distance ----
#============================================================================#
chromosome.stats = read.table(file = "tables/chromosome.stats.csv", header = TRUE, sep = ";")
df_dist2telomere = data.frame(set = character(0), chromosome = character(0), rec.rate = numeric(0), dist2telomere = numeric(0))

# After determining the center of the chromosome, the distance to the telomere is simply physical position for the first half
# and absolute(max - physical position) for the other half
# In other words, it is also the minimal value of distances taken from both tips of the chromosome
# YET this means we could do an averaging of two contrasted patterns, i.e. it can erase a true signal
# Hence we needed to consider both halves separately, but they were not independent
# Only one half for each chromosome was randomly sampled to avoid pseudo-replicates
list_maps = data.frame(set = as.character(chromosome.stats$set), chromosome = as.character(chromosome.stats$chromosome), stringsAsFactors = FALSE)

for(i in 1:nrow(list_maps)) {
  cat(list_maps[i,1], list_maps[i,2], "\n")
  recombination_map = read.table(paste("output/recombination_maps/loess/100kbwind/", list_maps[i,1], "_chromosome", list_maps[i,2], ".txt", sep = ""), header = TRUE)
  recombination_map = recombination_map[order(recombination_map$phys),]
  totalchrsize = max(recombination_map$phys, na.rm = TRUE)
  # Randomly subset half of the chromosome
  side = sample(c("left", "right"), size = 1)
  if (side == "left") {
    recombination_map = recombination_map[which(recombination_map$phys < totalchrsize/2),]
  } else {
    recombination_map = recombination_map[which(recombination_map$phys > totalchrsize/2),]
  }
  # Estimate distance to telomere
  recombination_map$dist2telomere = apply(data.frame(recombination_map$phys, abs(recombination_map$phys-totalchrsize)), 1, min)
  # plot(recombination_map$dist2telomere, recombination_map$rec.rate, main = map)
  # Results
  df = data.frame(set = rep(list_maps[i,1], nrow(recombination_map)), chromosome = rep(list_maps[i,2], nrow(recombination_map)),
                  rec.rate = recombination_map$rec.rate, dist2telomere = recombination_map$dist2telomere)
  # remove older estimates
  df_dist2telomere = df_dist2telomere[!(df_dist2telomere$set == list_maps[i,1] & df_dist2telomere$chromosome == list_maps[i,2]),]
  # Append data
  df_dist2telomere = rbind(df_dist2telomere, df)
}
df_dist2telomere

write.table(df_dist2telomere, "output/dist2telomere/AllChromosomes_Dist2telomeres.txt", row.names = FALSE, col.names = TRUE)

#============================================================================#
#   Raw data - Rec. rate ~ distance to nearest chromosome end ----
#============================================================================#
df_dist2telomere = read.table("output/dist2telomere/AllChromosomes_Dist2telomeres.txt", header = TRUE)

# Remove NA and transform null values
dist2telomere = subset(df_dist2telomere, (!is.na(df_dist2telomere$rec.rate)))
# dist2telomere$rec.rate = dist2telomere$rec.rate + 0.0001
# dist2telomere$dist2telomere = dist2telomere$dist2telomere + 0.0001
# Remove null values
dist2telomere = dist2telomere[(dist2telomere$rec.rate != 0),]
dist2telomere = dist2telomere[(dist2telomere$dist2telomere != 0),]
# Species names
dist2telomere = merge(dist2telomere, unique(chromosome.stats[,c(1, 24)]), by = "set", all.x = TRUE)

# Distribution of data
hist(dist2telomere$dist2telomere, breaks = 20, xlab = "Distance to telomere (Mb)", main = "")
hist(log(dist2telomere$dist2telomere), breaks = 20, xlab = "Distance to telomere (log-transformed)", main = "")

hist(dist2telomere$rec.rate, breaks = 20, xlab = "Recombination rate (cM/Mb)", main = "")
hist(log(dist2telomere$rec.rate), breaks = 20, xlab = "Recombination rate (log-transformed)", main = "")

# A figure with all data points (100kb windows)
# And LM (red) + LMM (blue, species effect) regression lines
plot(log(dist2telomere$dist2telomere), log(dist2telomere$rec.rate),
     xlab = "Distance to telomere (log-transformed)",
     ylab = "Recombination rate (log-transformed)")

# A simple regression model
mod_dist2telomere = lm(log(rec.rate) ~ log(dist2telomere), data = dist2telomere)
abline(mod_dist2telomere, col = "Red")
# Linear Mixed Model
lmm_dist2telomere = lmer(log(rec.rate) ~ log(dist2telomere) + (log(dist2telomere) | species), data = dist2telomere)
summary(lmm_dist2telomere)
# Significance of coefficients
anova(lmm_dist2telomere)

# Mixed regression line
abline(a = fixef(lmm_dist2telomere)[1], b = fixef(lmm_dist2telomere)[2], col = "Blue")

# Correlation
cor.test(log(dist2telomere$dist2telomere), log(dist2telomere$rec.rate), method = "spearman")




#============================================================================#
#   Rec. rate ~ relative distance to nearest chromosome end (quantiles of distances) ----
#============================================================================#

df_dist2telomere = read.table("output/dist2telomere/AllChromosomes_Dist2telomeres.txt", header = TRUE)
# Compute the quantiles of distances
# Prepare the dataset with species and chromosomes to gather
nquantiles = 20
df = data.frame(set = character(0), chromosome = character(0),
                rec.rate = numeric(0), lower = numeric(0), upper = numeric(0),
                quantile.lower = numeric(0), quantile.upper = numeric(0), dist2telomere = numeric(0))
# Compute the quantiles of distances for each species and chromosome
list_maps = data.frame(set = as.character(chromosome.stats$set), chromosome = as.character(chromosome.stats$chromosome), stringsAsFactors = FALSE)
for (i in 1:nrow(list_maps)) {
  cat(list_maps[i,1], list_maps[i,2], "\n")
  df_subset = df_dist2telomere[which(df_dist2telomere$set == list_maps[i,1] & df_dist2telomere$chromosome == list_maps[i,2]),]
  # Boundaries of 20 bins of equal distance
  # quantiles = seq(df_subset$dist2telomere, probs = seq(0, 1, length.out = nquantiles+1))
  bins = seq(from = 0, to = max(df_dist2telomere$dist2telomere[which(df_dist2telomere$set == list_maps[i,1] & df_dist2telomere$chromosome == list_maps[i,2])], na.rm = TRUE), length.out = nquantiles+1)
  quantile.lower = bins[1:nquantiles]
  quantile.upper = bins[2:(nquantiles+1)]
  newdata = data.frame(set = rep(as.character(unique(df_subset$set)), each = nquantiles), chromosome = rep(as.character(unique(df_subset$chromosome)), each = nquantiles),
                       rec.rate = NA, lower = NA, upper = NA,
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
# Compute the mean recombination rate in distance intervals
# Take a vector of indexes in the df dataframe
# and return a vector of rec.rate
# Standardized recombination rate among chromosomes (i.e. each chromosome standardized independently)
# for (i in 1:nrow(list_maps)) {
#   cat(list_maps$set[i], list_maps$chromosome[i], "\n")
#   df_dist2telomere$rec.rate[which(df_dist2telomere$set == list_maps$set[i] & df_dist2telomere$chromosome == list_maps$chromosome[i])] =
#     scale(df_dist2telomere$rec.rate[which(df_dist2telomere$set == list_maps$set[i] & df_dist2telomere$chromosome == list_maps$chromosome[i])])
# }
rec_in_distance_quantile = function(idx) {
  # subset = subset(df_dist2telomere, set == df$set[idx] & chromosome == df$chromosome[idx])
  # rec.rate.quantile = mean(subset$rec.rate, na.rm = TRUE)
  rec.rate.quantile = mean(df_dist2telomere$rec.rate[df_dist2telomere$set == df$set[idx]
                                                     & df_dist2telomere$chromosome == df$chromosome[idx]
                                                     & df_dist2telomere$dist2telomere >= df$quantile.lower[idx]
                                                     & df_dist2telomere$dist2telomere <= df$quantile.upper[idx]], na.rm = TRUE)
  return(rec.rate.quantile)
}

idx = c(1:nrow(df))
require(pbmcapply)
rec.rate = unlist(pbmclapply(idx, function(x) rec_in_distance_quantile(x), mc.cores =  7))
df$rec.rate = rec.rate
hist(df$rec.rate)


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
write.table(df, "output/dist2telomere/DistancesRelative_bins.txt", row.names = FALSE, col.names = TRUE)

df = read.table("output/dist2telomere/DistancesRelative_bins.txt", header = TRUE)

#----------------------------------------------------------------------------#
# Scale recombination rates
for (i in 1:nrow(list_maps)) {
  cat(list_maps$set[i], list_maps$chromosome[i], "\n")
  df$rec.rate[which(df$set == list_maps$set[i] & df$chromosome == list_maps$chromosome[i])] =
    scale(df$rec.rate[which(df$set == list_maps$set[i] & df$chromosome == list_maps$chromosome[i])])
}
hist(df$rec.rate)
write.table(df, "output/dist2telomere/DistancesRelativeScaled_bins.txt", row.names = FALSE, col.names = TRUE)

df = read.table("output/dist2telomere/DistancesRelativeScaled_bins.txt", header = TRUE)



#----------------------------------------------------------------------------#
# FIGURES ----
#----------------------------------------------------------------------------#
# Quantiles of rec. rate ~ distance to telomere for each chromosome
# one color per species (set)
DistancesRelative_bins = ggplot(data = df, aes(x = dist2telomere, y = rec.rate)) +
  geom_line(aes(group = interaction(set, chromosome), colour = set)) +
  geom_point(aes(colour = set)) +
  xlab("Relative distance to telomere") + ylab("Recombination rate (cM/Mb)") +
  stat_smooth(aes(y = rec.rate), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Black", fill = "DarkGrey", alpha = 0.3) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
DistancesRelative_bins

# Unreadable, too much lines

#----------------------------------------------------------------------------#
# Split the dataset to identify groups/patterns ----
#----------------------------------------------------------------------------#
# Figures per classes of chromosome sizes?
# Indeed there seem to have co-existence of different patterns
hist(chromosome.stats$phys.map.length, breaks = 60)
abline(v = 4e+08, col = "Red")
abline(v = 0.85e+08, col = "Red", lty = 2)
# Species with larger chromosomes
largerChr = chromosome.stats$set[which(chromosome.stats$phys.map.length > 4e+08)]
df_largerChr = subset(df, df$set %in% largerChr)
DistancesRelative_binsLargerChr = ggplot(data = df_largerChr, aes(x = dist2telomere, y = rec.rate)) +
  # geom_line(aes(group = interaction(set, chromosome), colour = set)) +
  geom_point(aes(colour = set, alpha = 0.5)) +
  stat_smooth(aes(y = rec.rate), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Black", fill = "DarkGrey", alpha = 0.3) +
  xlab("Relative distance to telomere") + ylab("Recombination rate (cM/Mb)") +
  xlim(0, 0.5) +
  # geom_abline(intercept = mod$coefficients[1], slope = mod$coefficients[2], linetype = "dashed") + # The linear regression line
  # geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
DistancesRelative_binsLargerChr

# Medium chromosomes
mediumChr = chromosome.stats$set[which(chromosome.stats$phys.map.length > 0.85e+08 & chromosome.stats$phys.map.length < 4e+08)]
df_mediumChr = subset(df, df$set %in% mediumChr)
DistancesRelative_binsMediumChr = ggplot(data = df_mediumChr, aes(x = dist2telomere, y = rec.rate)) +
  # geom_line(aes(group = interaction(set, chromosome), colour = set)) +
  geom_point(aes(colour = set, alpha = 0.5)) +
  stat_smooth(aes(y = rec.rate), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Black", fill = "DarkGrey", alpha = 0.3) +
  xlab("Relative distance to telomere") + ylab("Recombination rate (cM/Mb)") +
  xlim(0, 0.5) +
  # geom_abline(intercept = mod$coefficients[1], slope = mod$coefficients[2], linetype = "dashed") + # The linear regression line
  # geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
DistancesRelative_binsMediumChr
# And smaller chromosomes
smallerChr = chromosome.stats$set[which(chromosome.stats$phys.map.length < 0.85e+08)]
df_smallerChr = subset(df, (df$set %in% smallerChr))
DistancesRelative_binsSmallerChr = ggplot(data = df_smallerChr, aes(x = dist2telomere, y = rec.rate)) +
  # geom_line(aes(group = interaction(set, chromosome), colour = set)) +
  geom_point(aes(colour = set, alpha = 0.5)) +
  stat_smooth(aes(y = rec.rate), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Black", fill = "DarkGrey", alpha = 0.3) +
  xlab("Relative distance to telomere") + ylab("Recombination rate (cM/Mb)") +
  xlim(0, 0.5) +
  # geom_abline(intercept = mod$coefficients[1], slope = mod$coefficients[2], linetype = "dashed") + # The linear regression line
  # geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
DistancesRelative_binsSmallerChr

ggarrange(DistancesRelative_binsLargerChr, DistancesRelative_binsMediumChr, DistancesRelative_binsSmallerChr,
          ncol = 3)

#----------------------------------------------------------------------------#
# Chromosomes pooled per species ----
#----------------------------------------------------------------------------#
df_pooled = read.table("output/dist2telomere/DistancesRelativeScaled_bins.txt", header = TRUE)
# Quantiles of relative distances are pooled together
df_pooled$species = gsub("_[A-Za-z0-9]*$", "", df_pooled$set)
df_pooled$species = gsub("_MaizeGDBConsensus", "", df_pooled$species)
df_pooled$species = gsub("_", " ", df_pooled$species)

df_pooled = aggregate(df_pooled$rec.rate~df_pooled$dist2telomere+df_pooled$species, FUN=mean)
colnames(df_pooled) = c("dist2telomere", "species", "rec.rate")

# Quadratic Phylogenetic Regression
# Load the phylogenetic tree estimated in 'Phylogeny.R'
load("output/phylogeny/tree.Rdata")
tree
tree$tip.label
# Drop unsampled tips
library(ape)
tree = drop.tip(tree, tip = tree$tip.label[!(tree$tip.label %in% df_pooled$species)])
tree

# PGLMM (model selection, model minimizing AIC)
hist(df_pooled$rec.rate)
library(phyr)
# pglmm.model = pglmm(rec.rate ~ dist2telomere + (1|species), data = df_pooled, family = 'gaussian',
#                     cov_ranef = list(species = tree))
pglmm.model = pglmm(rec.rate ~ dist2telomere + I(dist2telomere^2) + (1|species), data = df_pooled, family = 'gaussian',
                    cov_ranef = list(species = tree))
summary(pglmm.model)

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
pglmm.model = pglmm(rec.rate ~ dist2telomere + I(dist2telomere^2) + (1|species), data = df_pooled, family = 'gaussian',
                    cov_ranef = list(species = tree), bayes = TRUE)
summary(pglmm.model)
# Quadratic function and C.I.
quadrafun = function(x, c) fixef(pglmm.model)[[1]][1] + fixef(pglmm.model)[[1]][2]*x + fixef(pglmm.model)[[1]][3]*x^2
lowerquadrafun = function(x, c) fixef(pglmm.model)[[2]][1] + fixef(pglmm.model)[[2]][2]*x + fixef(pglmm.model)[[2]][3]*x^2
upperquadrafun = function(x, c) fixef(pglmm.model)[[3]][1] + fixef(pglmm.model)[[3]][2]*x + fixef(pglmm.model)[[3]][3]*x^2
# The C.I. values
df.new = data.frame(dist2telomere = df_pooled$dist2telomere[1:20])
df.new$lwr.pred = lowerquadrafun(df.new$dist2telomere)
df.new$upr.pred = upperquadrafun(df.new$dist2telomere)

DistancesRelative_pooledChromosomes = ggplot(data = df_pooled, aes(x = dist2telomere, y = rec.rate)) +
  # geom_line(aes(group = species, colour = species)) +
  geom_point(aes(colour = species), alpha = 0.4) +
  xlab("Relative distance to telomere") + ylab("Recombination rate (cM/Mb)") +
  # stat_smooth(aes(y = rec.rate), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Black", fill = "DarkGrey", alpha = 0.3) +
  stat_function(fun = quadrafun, colour = "Black", size = 1) +
  geom_ribbon(data = df.new, aes(x = dist2telomere, ymin = lwr.pred, ymax = upr.pred), alpha = 0.4, inherit.aes = F, fill = "Black") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
DistancesRelative_pooledChromosomes

#----------------------------------------------------------------------------#
# Distribution of correlations between distance to telomere & recombination ----
#----------------------------------------------------------------------------#
# One correlation coefficient per chromosome
df_dist2telomere = read.table("output/dist2telomere/AllChromosomes_Dist2telomeres_relativedistances.txt", header = TRUE)
library(plyr)
# Trim for NA values
df_dist2telomere_trimmed = subset(df_dist2telomere,
                                  !(is.na(df_dist2telomere$rec.rate) | is.na(df_dist2telomere$dist2telomere)))
cor_dist2telomere = ddply(df_dist2telomere_trimmed, c("set", "chromosome"),
      function(x) cor(x$rec.rate, x$dist2telomere, method = "spearman"))
colnames(cor_dist2telomere)[3] = "correlation"

hist(cor_dist2telomere$correlation, freq = TRUE,
     breaks = 20, xlab = "Correlation", main = "")


#----------------------------------------------------------------------------#
# Test by random resampling ----
# Distribution of correlations between distance to telomere & recombination
#----------------------------------------------------------------------------#
# One correlation coefficient per chromosome
# Is recombination significantly correlated to distance to telomere
# H0: random reshuffling of recombination rates
# 1,000 iterations
df_dist2telomere_resampling = df_dist2telomere_trimmed
df_dist2telomere_resampling$rec.rate = sample(df_dist2telomere_resampling$rec.rate, replace = FALSE)
cor_dist2telomere_resampling = ddply(df_dist2telomere_resampling, c("set", "chromosome"),
                          function(x) cor(x$rec.rate, x$dist2telomere, method = "spearman"))
colnames(cor_dist2telomere_resampling)[3] = "correlation"
hist(cor_dist2telomere_resampling$correlation, breaks = 20, xlab = "Correlation", main = "")

hist(cor_dist2telomere$correlation, border = "red", freq = FALSE, breaks = 20,
     xlab = "Correlation", main = "", ylim = c(0, 5))
hist(cor_dist2telomere_resampling$correlation, border = "black", freq = FALSE,
     breaks = 20, add = TRUE)

# Is it really interesting to prove the obvious?






#----------------------------------------------------------------------------#
# Identify periphery-bias ----
#----------------------------------------------------------------------------#
# Take the problem in a different approach, less a priori on hypotheses (chromosome size)
# And more data driven
# Discriminate chromosomes with a bias toward periphery and others.
# Chromosomes with significantly more recombination in quantiles < 0.25 than > 0.25
list_maps$peripherybias = NA
# prop = 0.15
for (i in 1:nrow(list_maps)) {
  cat(list_maps$set[i], list_maps$chromosome[i], "\n")
  # mean_distal = mean(df$rec.rate[which(df$set == list_maps$set[i] & df$chromosome == list_maps$chromosome[i] & df$dist2telomere < prop)], na.rm = TRUE)
  # mean_proximal = mean(df$rec.rate[which(df$set == list_maps$set[i] & df$chromosome == list_maps$chromosome[i] & df$dist2telomere > prop)], na.rm = TRUE)
  # if (is.na(mean_proximal) | is.na(mean_distal)) {
  #   list_maps$peripherybias[i] = NA
  # } else {
  #   if (wilcox.test(mean_distal, mean_proximal)$p.value < 0.05 & (mean_distal > mean_proximal)) {
  #     list_maps$peripherybias[i] = TRUE
  #   } else {
  #     list_maps$peripherybias[i] = FALSE
  #   }
  # }
  
  # Mode of the distribution
  if (length(which.max(df$rec.rate[which(df$set == list_maps$set[i] & df$chromosome == list_maps$chromosome[i])])) > 0 ) {
    if ((which.max(df$rec.rate[which(df$set == list_maps$set[i] & df$chromosome == list_maps$chromosome[i])]) < 3)) {
      list_maps$peripherybias[i] = TRUE
    } else {
      list_maps$peripherybias[i] = FALSE
    }
  } else {
    list_maps$peripherybias[i] = FALSE
  }
  

   
}
table(list_maps$peripherybias)

# Species with a periphery-bias
df_peripherybias = subset(df, (df$set %in% list_maps$set[list_maps$peripherybias == TRUE] & df$chromosome %in% list_maps$chromosome[list_maps$peripherybias == TRUE]))
DistancesRelative_binsPeripheryBias = ggplot(data = df_peripherybias, aes(x = dist2telomere, y = rec.rate)) +
  # geom_line(aes(group = interaction(set, chromosome), colour = set)) +
  geom_point(aes(colour = set)) +
  xlab("Relative distance to telomere") + ylab("Recombination rate (cM/Mb)") +
  # geom_abline(intercept = mod$coefficients[1], slope = mod$coefficients[2], linetype = "dashed") + # The linear regression line
  # geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
DistancesRelative_binsPeripheryBias

# Species with larger chromosomes
df_noperipherybias = subset(df, !(df$set %in% list_maps$set[list_maps$peripherybias == TRUE] & df$chromosome %in% list_maps$chromosome[list_maps$peripherybias == TRUE]))
DistancesRelative_binsNoPeripheryBias = ggplot(data = df_noperipherybias, aes(x = dist2telomere, y = rec.rate)) +
  # geom_line(aes(group = interaction(set, chromosome), colour = set)) +
  geom_point(aes(colour = set)) +
  xlab("Relative distance to telomere") + ylab("Recombination rate (cM/Mb)") +
  # geom_abline(intercept = mod$coefficients[1], slope = mod$coefficients[2], linetype = "dashed") + # The linear regression line
  # geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
DistancesRelative_binsNoPeripheryBias

#----------------------------------------------------------------------------#
# Compute the mode of each line (chromosome) as a proxy of the spatial bias toward the periphery ----
#---------------------------------------------------------------------------- #
chromosome.stats$modeDist2telomere = NA
for (i in 1:nrow(chromosome.stats)) {
  subset = subset(df, df$set == chromosome.stats$set[i] & df$chromosome == chromosome.stats$chromosome[i])
  chromosome.stats$modeDist2telomere[i] = subset$dist2telomere[which.max(subset$rec.rate)]
}
hist(chromosome.stats$modeDist2telomere, breaks = 20)
summary(chromosome.stats$modeDist2telomere)
cor.test(chromosome.stats$modeDist2telomere, chromosome.stats$phys.map.length/1000000, method = "spearman")
plot(log(chromosome.stats$phys.map.length/1000000), chromosome.stats$modeDist2telomere)
# Is there a statistical relationships when taking into account species autocorrelation?
# lmer.model = lmer(modeDist2telomere~log(phys.map.length/1000000) + (1|species), data = chromosome.stats)
# summary(lmer.model)
# (lmer.model)
# plot(lmer.model)
load("output/phylogeny/tree_reduced.Rdata")
# Branch lengths all set to 1
tree$edge.length = rep(1, length(tree$edge))
pglmm.model = pglmm(modeDist2telomere~log(phys.map.length/1000000) + (1|species__), data = chromosome.stats, family = 'gaussian',
                    cov_ranef = list(species = tree))
summary(pglmm.model)
# Diagnostic plots
qqPlot(residuals(pglmm.model))
plot(fitted(pglmm.model), residuals(pglmm.model))
abline(h = 0)
# Plot the fitted regression line over data
plot(log(chromosome.stats$phys.map.length/1000000), chromosome.stats$modeDist2telomere)
# abline(a = fixef(lmer.model)[1], b = fixef(lmer.model)[2])
abline(a = fixef(pglmm.model)[[1]][1], b = fixef(pglmm.model)[[1]][2])

# Yet this gives only the information of where the recombination is maximal in the genome
# Apparently it is more often in the end of the chromosome
# BUT gives no information on the strenght of this bias
# How much greater is it in the end that in the rest of the chromosome?
# A ratio end of the chromosome/rest of the chromosome (i.e. ratio bin of the tip/mean of all other bins)
# OR a regression line on the bins, per chromosome, where the slope is the strength of the bias
# Compare both estimators between them and with observed data
# and make correlations with chromosome size

# For a chromosome only
i = 8
subs = subset(df_dist2telomere, set == list_maps[i,1] & chromosome == list_maps[i,2])
plot(subs$dist2telomere, subs$rec.rate)
# What can we dot with this?
# Problem is there is two sides with different landscapes
# Hence the use of avergaing in 20 bins
subs = subset(df, set == list_maps[i,1] & chromosome == list_maps[i,2])
plot(subs$dist2telomere, subs$rec.rate)
abline(v = chromosome.stats$modeDist2telomere[which(chromosome.stats$set == list_maps[i,1] & chromosome.stats$chromosome == list_maps[i,2])], col = "Red") # Display the mode
mod = lm(rec.rate~dist2telomere, data = subs)
abline(mod)
abline(h = mean(subs$rec.rate[-1], na.rm = TRUE), col = "Blue", lty = 2)
# The ratio of end/rest of the chromosome
subs$rec.rate[1]/mean(subs$rec.rate, na.rm = TRUE)
# The higher the ratio, the more the recombination is concentrated in telomeric regions

#----------------------------------------------------------------------------#
# Compute the periphery recombination bias ----
# as a ratio of recombination at the end of the chromosome divided by total mean recombination
#---------------------------------------------------------------------------- #
# chromosome.stats$peripherybias_regression = NA
chromosome.stats$peripherybias_ratio = NA
for (i in 1:nrow(chromosome.stats)) {
  subs = subset(df, set == chromosome.stats$set[i] & chromosome == chromosome.stats$chromosome[i])
  # plot(subs$dist2telomere, subs$rec.rate)
  # abline(v = chromosome.stats$modeDist2telomere[which(chromosome.stats$set == list_maps[i,1] & chromosome.stats$chromosome == list_maps[i,2])], col = "Red") # Display the mode
  # mod = lm(rec.rate~dist2telomere, data = subs)
  # chromosome.stats$peripherybias_regression[i] = as.numeric(coefficients(mod)[2])
  # abline(mod)
  # abline(h = mean(subs$rec.rate[-1], na.rm = TRUE), col = "Blue", lty = 2)
  # The ratio of end/rest of the chromosome
  N = 1
  chromosome.stats$peripherybias_ratio[i] = (mean(subs$rec.rate[c(1:N)], na.rm = TRUE)/mean(subs$rec.rate, na.rm = TRUE))
  # chromosome.stats$peripherybias_proportionrec[i] = (sum(subs$rec.rate[1:N], na.rm = TRUE)/sum(subs$rec.rate, na.rm = TRUE)) # Proportion of the total recombination in the N first bin,
  # that is exactly the same thing as the ratio; the proportion of recombination concentrated in the end of the chromosome
}

chromosome.stats$peripherybias_ratio[which(chromosome.stats$peripherybias_ratio == 0)] = NA
write.table(chromosome.stats, file = "tables/chromosome.stats.csv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ";")

chromosomes = read.table("tables/chromosome.stats.csv", header = TRUE, sep = ";")
# Remove NAs
chromosomes = chromosomes[which(!is.na(chromosomes$peripherybias_ratio)),]

# plot(log(chromosomes$peripherybias_ratio), chromosomes$peripherybias_regression)
# plot(log(chromosomes$phys.map.length/1000000), chromosomes$peripherybias_regression) # Regression DOES NOT WORK
# plot(log(chromosomes$phys.map.length/1000000), log(chromosomes$peripherybias_ratio))
# Regression DOES NOT WORK
# BUT the peripheric ratio seems to work well
cor.test(log(chromosomes$phys.map.length/1000000), log(chromosomes$peripherybias_ratio), method = "spearman")
mod = lm(log(chromosomes$peripherybias_ratio)~log(chromosomes$phys.map.length/1000000))
# lmer.model = lmer(log(chromosomes$peripherybias_ratio)~log(chromosomes$phys.map.length/1000000) + (1|chromosomes$species))
# summary(lmer.model)
# (lmer.model)
load("output/phylogeny/tree_reduced.Rdata")
# Branch lengths all set to 1
tree$edge.length = rep(1, length(tree$edge))
pglmm.model = pglmm(log(peripherybias_ratio)~log(phys.map.length/1000000) + (1|species__), data = chromosome.stats, family = 'gaussian',
                    cov_ranef = list(species = tree))
summary(pglmm.model)
# Diagnostic plots
# qqPlot(residuals(pglmm.model))
# plot(fitted(pglmm.model), residuals(pglmm.model))
# abline(h = 0)

# plot(lmer.model)
# Plot the fitted regression line over data
plot(log(chromosomes$phys.map.length/1000000), log(chromosomes$peripherybias_ratio))
abline(mod, lty = 2)
# abline(a = fixef(lmer.model)[1], b = fixef(lmer.model)[2])
abline(a = fixef(pglmm.model)[[1]][1], b = fixef(pglmm.model)[[1]][2])
# the ratio seems independant of the phylogeny
# Yet, we need to verify this method carefully
# If recombination higher in distal regions, then log(ratio) would be > 0
abline(h = 0, col = "Blue")
# We observed that log(ratio) are distributed both above and below 0
# No universal feature, even if there was correlation


# Problems were:
# Averaging two halves of a chromosome (consider chromosomes halves separately? and randomly sample only one half to avoid pseudo-replicates)
# Not taking care of the centromere (yet, it is also one hypothesis we want to test with a different approach)
# Arbitrary sampling of the ratio (why only 1 bin divided by 19?) -> sensitivity analysis of the sampling, compare mean periphery ratio for 1/19 with 2/18, 3/17...
# Gradient of ratios
#----------------------------------------------------------------------------#
# ASSESSING SENSITIVITY TO SAMPLING IN RATIO - HOW MANY BINS? ----
#----------------------------------------------------------------------------#
# Compute a gradient ratio for all bin number
gradient_ratio = data.frame(set = rep(chromosomes$set, 20), chromosome = rep(chromosomes$chromosome, 20), bin = rep(1:20, each = nrow(chromosomes)), ratio = NA)
for (i in 1:nrow(gradient_ratio)) {
  subs = subset(df, set == gradient_ratio$set[i] & chromosome == gradient_ratio$chromosome[i])
  # gradient_ratio$ratio[i] = (sum(subs$rec.rate[c(1:gradient_ratio$bin[i])], na.rm = TRUE)/sum(subs$rec.rate, na.rm = TRUE))
  # gradient_ratio[i, j] = (mean(subs$rec.rate[c(1:j)], na.rm = TRUE)/mean(subs$rec.rate[c(j:20)], na.rm = TRUE))/mean(subs$rec.rate, na.rm = TRUE) # Standardized ratio by the mean recombination rate
  gradient_ratio$ratio[i] = (mean(subs$rec.rate[c(1:gradient_ratio$bin[i])], na.rm = TRUE)/mean(subs$rec.rate, na.rm = TRUE))
}

# If recombination higher in distal regions, then log(ratio) would be > 0
plot(x = gradient_ratio$bin, y = log(gradient_ratio$ratio))
abline(h = 0)
# We observed that log(ratio) are distributed both above and below 0
# No universal feature, even if there was correlation

gradient_ratio_plot = ggplot(data = gradient_ratio, aes(x = bin, y = log(ratio))) +
  geom_point(aes(), shape=1, colour = "darkGrey") +
  geom_smooth(aes(x = bin, y = log(ratio)), method = "loess", se = TRUE, level = 0.95) +
  geom_hline(yintercept = 0) +
  xlab("Number of bins sampled at chromosome end") + ylab("Recombination periphery bias (ratio)") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
gradient_ratio_plot

#----------------------------------------------------------------------------#
#   TESTING THE HYPOTHESIS OF A BIAS OF RECOMBINATION TOWARD DISTAL REGIONS OF THE CHROMOSOME ----
#----------------------------------------------------------------------------#
# Test the hypothesis by resampling
# Randomly reshuffling bins to compute a mean theoretical ratio under H0 (no bias toward the chromosome end) and compare with the mean ratio

# Building H0 bootstrapped distribution


# Testing observed mean ratio vs H0















#============================================================================#
#   Rec. rates ~ distance quantiles; chromosomes pooled per species ----
#============================================================================#

# COMPUTE DATA
# Recombination rate expressed as a function of the distance to the telomere for all species independently
# Load data
metadata.clean = read.table("data-cleaned/Genetic_maps_ressources.csv", header = TRUE, sep = ";")

# List of maps
list = system(paste("ls ", wd, "/output/recombination_maps/loess/100kbwind/", sep = ""), intern = TRUE)
list = list[grep("*.txt", list)]
list = gsub(".txt", "", list)
list
# Initialize an empty file
df_dist2telomere = data.frame(set = character(0), chromosome = character(0), rec.rate = numeric(0), dist2telomere = numeric(0))
write.table(df_dist2telomere, "output/dist2telomere/AllMaps.txt", row.names = FALSE, col.names = TRUE)

# After determining the center of the chromosome, the distance to the telomere is simply physical position for the first half
# and absolute(max - physical position) for the other half
# In other words, it is also the minimal value of distances taken from both tips of the chromosome
# map = list[2]
for(map in list) {
  cat(map, "\n")
  set = gsub("_chromosome[A-Za-z0-9]*", "", map)
  chr = gsub("[A-Za-z0-9_]*_chromosome", "", map)
  recombination_map = read.table(paste("output/recombination_maps/loess/100kbwind/", map, ".txt", sep = ""), header = TRUE)
  recombination_map$dist2telomere = apply(data.frame(recombination_map$phys, abs(recombination_map$phys-max(recombination_map$phys, na.rm = TRUE))), 1, min)
  # plot(recombination_map$dist2telomere, recombination_map$rec.rate, main = map)
  # Results
  df = data.frame(set = rep(set, nrow(recombination_map)), chromosome = rep(chr, nrow(recombination_map)),
                  rec.rate = recombination_map$rec.rate, dist2telomere = recombination_map$dist2telomere)
  # Load file to write
  df_dist2telomere = read.table("output/dist2telomere/AllMaps.txt", header = TRUE)
  df_dist2telomere$set = as.character(df_dist2telomere$set)
  df_dist2telomere$chromosome = as.character(df_dist2telomere$chromosome)
  # remove older estimates
  df_dist2telomere = df_dist2telomere[!(df_dist2telomere$set == set & df_dist2telomere$chromosome == chr),]
  # Append data
  df_dist2telomere = rbind(df_dist2telomere, df)
  write.table(df_dist2telomere, "output/dist2telomere/AllMaps.txt", row.names = FALSE, col.names = TRUE)
}


#----------------------------------------------------------------------------#
# Analyses -----
#----------------------------------------------------------------------------#
df_dist2telomere = read.table("output/dist2telomere/AllMaps.txt", header = TRUE)
# Group values in quantiles of distances
nquantiles = 20 # Number of quantiles

# Initialize an empty file
df = data.frame(set = character(0), rec.rate = numeric(0), upper = numeric(0), lower = numeric(0), dist2telomere = numeric(0))
# write.table(df, "output/dist2telomere/Distances_quantiles.txt", row.names = FALSE, col.names = TRUE)

# Prepare the dataset with species and chromosomes to gather
# List of species
list_species = unique(df_dist2telomere$set)
# Compute the quantiles of distances for each species and chromsosome
for (set in list_species) {
  df_sp = df_dist2telomere[which(df_dist2telomere$set == set),]
  # Boundaries of 20 quantiles
  # Chromosomes pooled per species
  quantiles = quantile(df_dist2telomere$dist2telomere[df_dist2telomere$set == set], probs = seq(0, 1, length.out = nquantiles+1))
  quantile.lower = quantiles[1:nquantiles]
  quantile.upper = quantiles[2:nquantiles+1]
  df = rbind(df, data.frame(set = rep(set, nquantiles), rec.rate = NA, upper = NA, lower = NA, quantile.lower = quantile.lower, quantile.upper = quantile.upper))
  # # Process each chromosome
  # list_chromosomes = unique(df_sp$chromosome)
  # for (chr in list_chromosomes) {
  #   # Boundaries of 20 quantiles
  #   quantiles = quantile(df_dist2telomere$dist2telomere[df_dist2telomere$set == set & df_dist2telomere$chromosome == chr], probs = seq(0, 1, length.out = 21))
  #   quantile.lower = quantiles[1:20]
  #   quantile.upper = quantiles[2:21]
  #   df = rbind(df, data.frame(set = rep(set, 20), chromosome = rep(chr, 20), rec.rate = NA, upper = NA, lower = NA, quantile.lower = quantile.lower, quantile.upper = quantile.upper))
  # }
  rm(df_sp)
}

# Compute the mean recombination rate in distance interval
# Take a vector of indexes in the df dataframe
# and return a vector of rec.rate
rec_in_distance_quantile = function(idx) {
  rec.rate.quantile = mean(df_dist2telomere$rec.rate[df_dist2telomere$set == df$set[idx]
                                                     & df_dist2telomere$dist2telomere >= df$quantile.lower[idx]
                                                     & df_dist2telomere$dist2telomere <= df$quantile.upper[idx]], na.rm = TRUE)
  return(rec.rate.quantile)
}

idx = c(1:nrow(df))
require(pbmcapply)
rec.rate = unlist(pbmclapply(idx, function(x) rec_in_distance_quantile(x), mc.cores =  7))
df$rec.rate = rec.rate

# for (i in 1:25) {
#   cat(i, "\n")
#   # TODO bootstrap mean and CI
#   df$rec.rate[i] = mean(df_dist2telomere$rec.rate[df_dist2telomere$set == df$set[i]
#                                                   & df_dist2telomere$dist2telomere >= df$quantile.lower[i]
#                                                   & df_dist2telomere$dist2telomere <= df$quantile.upper[i]], na.rm = TRUE)
# }


# Save
write.table(df, "output/dist2telomere/Distances_quantiles.txt", row.names = FALSE, col.names = TRUE)
df = read.table("output/dist2telomere/Distances_quantiles.txt", header = TRUE)

#----------------------------------------------------------------------------#
# Graphical representation ----
#----------------------------------------------------------------------------#
# Trim the dataset to reduce the number of species
df$set = as.character(df$set)
list_set = as.character(unique(df$set))
list_set = list_set[1:10]

df = df[which(df$set %in% list_set),]

# Colors
# A color per dataset
color.set = data.frame(set = list_set, color = as.character(polychrome(length(list_set))))
# show_col(as.character(polychrome(length(unique(species)))))
df = merge(df, color.set, by = "set")
# Plot
plot(log(rowMeans(cbind(df$quantile.lower, df$quantile.upper), na.rm = TRUE)), df$rec.rate, xlab = "Distance to the telomere (Mb; log scale)", ylab = "Recombination rate (cM/Mb)",
     bg = as.character(df$color), pch = 21)
# segments(plotdf$gene_count, plotdf$lowerCI, plotdf$gene_count, plotdf$upperCI, col = as.character(plotdf$color))
# Add legend to top right, outside plot region
legend("topright", legend = unique(df$set), title = "Species",
       pt.bg = unique(as.character(df$color)), pch = 21, cex = 0.6)

#----------------------------------------------------------------------------#
# Represent species per phylogenetic families ----
#----------------------------------------------------------------------------#
phylofamilies = read.table("data-cleaned/species_metadata.csv", header = TRUE, sep = ";")
df$family = NA
for (i in 1:nrow(df)) {
  if(grepl("Zea_mays", as.character(df$set[i]))) {
    df$family[i] = as.character(phylofamilies$family[which(phylofamilies$species == "Zea_mays")])
  } else {
    df$family[i] = as.character(phylofamilies$family[which(phylofamilies$species == gsub("_[A-Za-z0-9]*$", "", df$set[i]))])
  }
}

png("figures/Dist2telomeres.png", width = 1600, height = 1600)
par(mfrow = c(2,2))
# Trim the dataset to select Poaceae
df_poaceae = df[which(df$family == "Poaceae"),]
# df = df[which(df$set == "Arabidopsis_thaliana_Serin2017"),]
# Colors
# A color per dataset
list_set = as.character(unique(df_poaceae$set))
color.set = data.frame(set = list_set, color = as.character(polychrome(length(list_set))))
# show_col(as.character(polychrome(length(unique(species)))))
df_poaceae = merge(df_poaceae, color.set, by = "set")
# Plot
plot(rowMeans(cbind(df_poaceae$quantile.lower, df_poaceae$quantile.upper), na.rm = TRUE), df_poaceae$rec.rate,
     xlab = "Distance to the telomere (Mb, log scale)",
     ylab = "Recombination rate (cM/Mb)", main = "Poaceae", lwd = 1.5, cex = 1.5, cex.main = 1.5,
     bg = as.character(df_poaceae$color), pch = 21, log = "x", cex.axis = 1.5, cex.lab = 1.5)
for (x in split(df_poaceae, df_poaceae$color)) lines(rowMeans(cbind(x$quantile.lower, x$quantile.upper), na.rm = TRUE), x$rec.rate,
                                                     col=as.character(x$color[1]), lwd = 2)

# segments(plotdf$gene_count, plotdf$lowerCI, plotdf$gene_count, plotdf$upperCI, col = as.character(plotdf$color))
# Add legend to top right, outside plot region
legend("topright", legend = unique(df_poaceae$set), title = "Species",
       pt.bg = unique(as.character(df_poaceae$color)), pch = 21, cex = 1.5)

# Trim the dataset to select Fabaceae and Brassicaceae
df_fabaceaebrassicaceae = df[which(df$family == "Fabaceae" | df$family == "Brassicaceae"),]
# df = df[which(df$set == "Arabidopsis_thaliana_Serin2017"),]
# Colors
# A color per dataset
list_set = as.character(unique(df_fabaceaebrassicaceae$set))
color.set = data.frame(set = list_set, color = as.character(polychrome(length(list_set))))
# show_col(as.character(polychrome(length(unique(species)))))
df_fabaceaebrassicaceae = merge(df_fabaceaebrassicaceae, color.set, by = "set")
# Plot
plot(rowMeans(cbind(df_fabaceaebrassicaceae$quantile.lower, df_fabaceaebrassicaceae$quantile.upper), na.rm = TRUE), df_fabaceaebrassicaceae$rec.rate, xlab = "Distance to the telomere (Mb, log scale)",
     ylab = "Recombination rate (cM/Mb)", main = "Fabaceae and Brassicaceae", lwd = 1.5, cex = 1.5, cex.main = 1.5,
     bg = as.character(df_fabaceaebrassicaceae$color), pch = 21, log = "x", cex.axis = 1.5, cex.lab = 1.5)
for (x in split(df_fabaceaebrassicaceae, df_fabaceaebrassicaceae$color)) lines(rowMeans(cbind(x$quantile.lower, x$quantile.upper), na.rm = TRUE), x$rec.rate,
                                                                               col=as.character(x$color[1]), lwd = 2)

# segments(plotdf$gene_count, plotdf$lowerCI, plotdf$gene_count, plotdf$upperCI, col = as.character(plotdf$color))
# Add legend to top right, outside plot region
legend("topright", legend = unique(df_fabaceaebrassicaceae$set), title = "Species",
       pt.bg = unique(as.character(df_fabaceaebrassicaceae$color)), pch = 21, cex = 1.5)

# Trim the dataset to select Euphorbiaceae, Salicaceae, Cucurbitaceae and Rubiaceae
# Seem to have a bell distribution
df_EuphSaliCucurRub = df[which(df$family == "Euphorbiaceae" | df$family == "Salicaceae" | df$family == "Cucurbitaceae" | df$family == "Rubiaceae"),]
# df = df[which(df$set == "Arabidopsis_thaliana_Serin2017"),]
# Colors
# A color per dataset
list_set = as.character(unique(df_EuphSaliCucurRub$set))
color.set = data.frame(set = list_set, color = as.character(polychrome(length(list_set))))
# show_col(as.character(polychrome(length(unique(species)))))
df_EuphSaliCucurRub = merge(df_EuphSaliCucurRub, color.set, by = "set")
# Plot
plot(rowMeans(cbind(df_EuphSaliCucurRub$quantile.lower, df_EuphSaliCucurRub$quantile.upper), na.rm = TRUE), df_EuphSaliCucurRub$rec.rate, xlab = "Distance to the telomere (Mb, log scale)",
     ylab = "Recombination rate (cM/Mb)", main = "Euphorbiaceae, Salicaceae, Cucurbitaceae and Rubiaceae", lwd = 1.5, cex = 1.5, cex.main = 1.5,
     bg = as.character(df_EuphSaliCucurRub$color), pch = 21, log = "x", cex.axis = 1.5, cex.lab = 1.5)
for (x in split(df_EuphSaliCucurRub, df_EuphSaliCucurRub$color)) lines(rowMeans(cbind(x$quantile.lower, x$quantile.upper), na.rm = TRUE), x$rec.rate,
                                                                       col=as.character(x$color[1]), lwd = 2)


# segments(plotdf$gene_count, plotdf$lowerCI, plotdf$gene_count, plotdf$upperCI, col = as.character(plotdf$color))
# Add legend to top right, outside plot region
legend("topright", legend = unique(df_EuphSaliCucurRub$set), title = "Species",
       pt.bg = unique(as.character(df_EuphSaliCucurRub$color)), pch = 21, cex = 1.5)

# Trim the dataset to select remaining families
df_remainingfamilies = df[-which(df$family == "Fabaceae" | df$family == "Brassicaceae"  | df$family == "Poaceae" |
                                   df$family == "Euphorbiaceae" | df$family == "Salicaceae" | df$family == "Cucurbitaceae" | df$family == "Rubiaceae"),]
# df = df[which(df$set == "Arabidopsis_thaliana_Serin2017"),]
# Colors
# A color per dataset
list_set = as.character(unique(df_remainingfamilies$set))
color.set = data.frame(set = list_set, color = as.character(polychrome(length(list_set))))
# show_col(as.character(polychrome(length(unique(species)))))
df_remainingfamilies = merge(df_remainingfamilies, color.set, by = "set")
# Plot
plot(rowMeans(cbind(df_remainingfamilies$quantile.lower, df_remainingfamilies$quantile.upper), na.rm = TRUE), df_remainingfamilies$rec.rate, xlab = "Distance to the telomere (Mb, log scale)",
     ylab = "Recombination rate (cM/Mb)", main = "Other families", lwd = 1.5, cex = 1.5, cex.main = 1.5,
     bg = as.character(df_remainingfamilies$color), pch = 21, log = "x", cex.axis = 1.5, cex.lab = 1.5)
for (x in split(df_remainingfamilies, df_remainingfamilies$color)) lines(rowMeans(cbind(x$quantile.lower, x$quantile.upper), na.rm = TRUE), x$rec.rate,
                                                                         col=as.character(x$color[1]), lwd = 2)


# segments(plotdf$gene_count, plotdf$lowerCI, plotdf$gene_count, plotdf$upperCI, col = as.character(plotdf$color))
# Add legend to top right, outside plot region
legend("topright", legend = unique(df_remainingfamilies$set), title = "Species",
       pt.bg = unique(as.character(df_remainingfamilies$color)), pch = 21, cex = 1.5)
par(mfrow = c(1,1))
dev.off()


#----------------------------------------------------------------------------#
# One figure per dataset - Batch processing
list_set = as.character(unique(df$set))
for (set in list_set) {
  subset = df[which(df$set == set),]
  jpeg(paste("output/dist2telomere/figures/", set, ".jpeg", sep = ""), width = 900, height = 900)
  plot(rowMeans(cbind(subset$quantile.lower, subset$quantile.upper), na.rm = TRUE), subset$rec.rate, xlab = "Distance to the telomere (Mb)", ylab = "Recombination rate (cM/Mb)")
  dev.off()
}

#----------------------------------------------------------------------------#
# Next step, compute the kernel of distance to telomeres? How are recombination events spatially distributed in the genome?


#============================================================================#
#
#   CENTROMERES: IS RECOMBINATION DRIVEN BY CENTROMERES? ----
#
#============================================================================#
# About the centromere, we have two hypotheses:
# (1) recombination is suppressed/lower in the centromeric region
# (2) asymmetry of the recombination landscape is driven by the centromere position


# I have access to data on centromere position for a large part of our dataset
# However, the position is not oriented in the same way as the Marey map/recombination map
# And we have no indication to find the correct orientation

# hence, I estimated a putative centromere position based on marey maps
# If i assume that (1) is exact (as it seems to be in many species), 
# the putative centromere position is then the part of the Marey function interpolated on the Marey map (loess regression) that has the lowest derivate (i.e. lowest recombination rate)

# I assumed that this putative centromere and the physical centromere of the C.I. should be on the same side of the chromosome
# Therefore it was possible to orient chromosomes and place the centromere at the correct position
source("sources/MareyMap.R")

# Estimate centromere positions for all Marey maps ----
chromosome.stats = read.table(file = "tables/chromosome.stats.csv", header = TRUE, sep = ";")

chromosome.stats$centromere_position_estimated = NA
for (i in 1:nrow(chromosome.stats)) {
  cat("Dataset", chromosome.stats$set[i], "chromosome", chromosome.stats$chromosome[i], "\n")
  chromosome.stats$centromere_position_estimated[i] = estim_centromere_pos(chromosome.stats$set[i], chromosome.stats$chromosome[i])
}

# Orient the Centromeric Index based on the putative centromere position
chromosome.stats$centromeric_index_position_oriented = NA
for (i in 1:nrow(chromosome.stats)) {
  cat("Dataset", chromosome.stats$set[i], "chromosome", chromosome.stats$chromosome[i], "\n")
  # If the ratio with centromere position is lower than 0.5,
  # then the CI and the centromere inferred are oriented the same way
  centromere_position = chromosome.stats$centromere_position_estimated[i]
  if (centromere_position / max(chromosome.stats$phys.map.length[i]/scale, na.rm = TRUE) < 0.5) {
    CI = chromosome.stats$centromeric_index[i]*chromosome.stats$phys.map.length[i]/scale
  } else {
    CI = (1 - chromosome.stats$centromeric_index[i])*chromosome.stats$phys.map.length[i]/scale
  }
  chromosome.stats$centromeric_index_position_oriented[i] = CI
}


# Save putative centromere position
write.table(chromosome.stats, file = "tables/chromosome.stats.csv",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ";")

#============================================================================#
#   Analyses on the new oriented Centromeric index ----
#============================================================================#
chromosome.stats = read.table(file = "tables/chromosome.stats.csv", header = TRUE, sep = ";")

# Diagnostic plot - Is the method working well?
# Correlation between the Centromeric Index and our putative centromere position
plot(chromosome.stats$centromere_position_estimated,
     chromosome.stats$centromeric_index_position_oriented,
     xlab = "Putative centromere position",
     ylab = "Centromeric Index position (oriented)",
     log = "xy")
abline(0, 1)

cor.test(chromosome.stats$centromere_position_estimated,
         chromosome.stats$centromeric_index_position_oriented)

# Discuss errors
# - the more the C.I. is close to 0.5, the less dramatic errors will be
# - we use the putative position only for orienting the C.I. If error in interpolation, 1/2 chance to have nonetheless the good orientation

# Method works well, so we could use directly the estimated putative centromere position as new data
# And make analyses with it.
# Besides, centromere position seems more accurate for our estimates than centromerix indexes from the litterature
# maybe due to differences between the assembly sequence and the cytological chromosome size
# Nonetheless, this correlation shows robustness of the method, despite inherent noise in the data.

#-------------------------------------------#
# HYPOTHESIS ONE ----
# (1) recombination is suppressed/lower in the centromeric region ----
# The reasoning is a bit circular because we inferred a centromere position as
# the position with the lowest recombination rate (lowest derivate of the interpolated Marey function)
# hence we can distinguish if it is truly a centromere biologically or just the result of a statistical definition
# Yet, the local suppressor effect of centromeres has been assessed in many species and is not yet to demonstrate

# Finally, hypothesis 1 is more an assumption used for inferences than an hypothesis to test:
# However, the clear correlation between the centromerix index fo the litterature and our putative centromeres inferred
# seems to be a clear evidence that centromeres are indeed local regions of reduced recombination
# The only approximation that is really made seems reasonable: it is that putative centromeres and centromeric indexes are on the same side of the chromosome

# Illustration
# FIGURES - Distance to centromeres


#============================================================================#
#   Estimate Rec. rate ~ relative distance to centromere ----
#============================================================================#
chromosome.stats = read.table(file = "tables/chromosome.stats.csv", header = TRUE,
                              sep = ";", stringsAsFactors = FALSE)
# Use the centromere position given in the literature
# Measure of the ratio short arm/long arm
# Relative size of each arm

#-------------------------#
# NUMBER OF CHROMOSOMES WITH CENTROMERE POSITION
# Maps without metadata about the centromeric index
unique(chromosome.stats$set[is.na(chromosome.stats$centromeric_index_position_oriented)])
# Maps with a centromeric index
unique(chromosome.stats$set[!is.na(chromosome.stats$centromeric_index_position_oriented)])
length(unique(chromosome.stats$set[!is.na(chromosome.stats$centromeric_index_position_oriented)]))

# Compute distances to the centromere
# Split short arms and long arms (arm_length = short/long)
df_dist2centromere = data.frame(set = character(0), chromosome = character(0),
                                rec.rate = numeric(0), dist2centromere = numeric(0),
                                arm_length = factor(levels = c("short", "long")))

# Remove maps without centromere information
list_maps = data.frame(set = as.character(chromosome.stats$set[which(!is.na(chromosome.stats$centromeric_index))]),
                       chromosome = as.character(chromosome.stats$chromosome[which(!is.na(chromosome.stats$centromeric_index))]),
                       stringsAsFactors = FALSE)


for(i in 1:nrow(list_maps)) {
  cat(list_maps[i,1], list_maps[i,2], "\n")
  recombination_map = read.table(paste("output/recombination_maps/loess/100kbwind/", list_maps[i,1], "_chromosome", list_maps[i,2], ".txt", sep = ""), header = TRUE)
  recombination_map = recombination_map[order(recombination_map$phys),]
  # Physical position of the centromere
  totalchrsize = max(recombination_map$phys, na.rm = TRUE)
  centromeric_position = chromosome.stats$centromeric_index_position_oriented[which(chromosome.stats$set == list_maps[i,1] &
                                                                 chromosome.stats$chromosome == list_maps[i,2])]
  # Identify short/long arm
  recombination_map$arm_length = NA
  if (centromeric_position < totalchrsize/2) {
    recombination_map$arm_length[which(recombination_map$phys < centromeric_position)] = "short"
    recombination_map$arm_length[which(recombination_map$phys >= centromeric_position)] = "long"
    
  } else {
    if (centromeric_position >= totalchrsize/2) {
      recombination_map$arm_length[which(recombination_map$phys <= centromeric_position)] = "long"
      recombination_map$arm_length[which(recombination_map$phys > centromeric_position)] = "short"
    }
  }
  # Compute distances
  # Absolute value of the difference window position - centromere position
  recombination_map$dist2centromere = abs(centromeric_position - recombination_map$phys)
  # Results
  df = data.frame(set = rep(list_maps[i,1], nrow(recombination_map)), chromosome = rep(list_maps[i,2], nrow(recombination_map)),
                  rec.rate = recombination_map$rec.rate, dist2centromere = recombination_map$dist2centromere,
                  arm_length = recombination_map$arm_length)
  # remove older estimates
  df_dist2centromere = df_dist2centromere[!(df_dist2centromere$set == list_maps[i,1] & df_dist2centromere$chromosome == list_maps[i,2]),]
  # Append data
  df_dist2centromere = rbind(df_dist2centromere, df)
}
df_dist2centromere

write.table(df_dist2centromere, "output/dist2centromere/AllChromosomes_Dist2centromeres.txt", row.names = FALSE, col.names = TRUE)


#============================================================================#
#   Raw data - Rec. rate ~ distance to centromere ----
#============================================================================#
df_dist2centromere = read.table("output/dist2centromere/AllChromosomes_Dist2centromeres.txt", header = TRUE)

# Remove NA and transform null values
dist2centromere = subset(df_dist2centromere, (!is.na(df_dist2centromere$rec.rate)))
# dist2centromere$rec.rate = dist2centromere$rec.rate + 0.0001
# dist2centromere$dist2centromere = dist2centromere$dist2centromere + 0.0001
# Remove null values
dist2centromere = dist2centromere[(dist2centromere$rec.rate > 0.001),]
dist2centromere = dist2centromere[(dist2centromere$dist2centromere != 0),]
# Species names
dist2centromere = merge(dist2centromere, unique(chromosome.stats[,c(1, 24)]), by = "set", all.x = TRUE)

# Distribution of data
hist(dist2centromere$dist2centromere, breaks = 20, xlab = "Distance to centromere (Mb)", main = "")
hist(log(dist2centromere$dist2centromere), breaks = 20, xlab = "Distance to centromere (log-transformed)", main = "")

hist(dist2centromere$rec.rate, breaks = 20, xlab = "Recombination rate (cM/Mb)", main = "")
hist(log(dist2centromere$rec.rate), breaks = 20, xlab = "Recombination rate (log-transformed)", main = "")


plot(dist2centromere$dist2centromere, dist2centromere$rec.rate,
     log = "xy",
     xlab = "Distance to centromere",
     ylab = "Recombination rate")

# A simple regression model
mod_dist2centromere = lm(log10(rec.rate) ~ log10(dist2centromere), data = dist2centromere)
abline(mod_dist2centromere, col = "Red", untf = FALSE)
# Linear Mixed Model
lmm_dist2centromere = lmer(log10(rec.rate) ~ log10(dist2centromere) + (log10(dist2centromere) | species), data = dist2centromere)
summary(lmm_dist2centromere)
# Significance of coefficients
anova(lmm_dist2centromere)

# Mixed regression line
abline(a = fixef(lmm_dist2centromere)[1], b = fixef(lmm_dist2centromere)[2],
       col = "Blue", unt = FALSE)

# Correlation
cor.test(dist2centromere$dist2centromere, dist2centromere$rec.rate, method = "spearman")

# Correlation is not very clear
# With LM and Spearman's correlation, recombination decreases with distance to the centromere
# While LMM show that recombination increases with distance to the centromere


#----------------------------------------------------------------------------#
df_dist2centromere = read.table("output/dist2centromere/AllChromosomes_Dist2centromeres.txt", header = TRUE)
# Compute the quantiles of distances
# Prepare the dataset with species and chromosomes to gather
# Subset short and long arms
nquantiles = 20
df = data.frame(set = character(0), chromosome = character(0),
                rec.rate = numeric(0), lower = numeric(0), upper = numeric(0),
                quantile.lower = numeric(0), quantile.upper = numeric(0),
                dist2centromere = numeric(0),
                arm_length = factor(levels = c("short", "long")))
# Compute the quantiles of distances for each species and chromosome
list_arms = data.frame(set = (rep(list_maps$set, each = 2)),
                       chromosome = rep(list_maps$chromosome, each = 2),
                       arm_length = rep(c("short", "long"), nrow(list_maps)),
                       stringsAsFactors = FALSE)

for (i in 1:nrow(list_arms)) {
  cat(list_arms[i,1], list_arms[i,2], list_arms[i,3], "\n")
  df_subset = df_dist2centromere[which(df_dist2centromere$set == list_arms[i,1] &
                                         df_dist2centromere$chromosome == list_arms[i,2] &
                                         df_dist2centromere$arm_length == list_arms[i,3]),]
  # Boundaries of 20 bins of equal distance
  bins = seq(from = 0, to = max(df_subset$dist2centromere, na.rm = TRUE), length.out = nquantiles+1)
  quantile.lower = bins[1:nquantiles]
  quantile.upper = bins[2:(nquantiles+1)]
  newdata = data.frame(set = rep(as.character(unique(df_subset$set)), each = nquantiles),
                       chromosome = rep(as.character(unique(df_subset$chromosome)), each = nquantiles),
                       rec.rate = NA, lower = NA, upper = NA,
                       quantile.lower = quantile.lower, quantile.upper = quantile.upper,
                       dist2centromere = (quantile.upper+quantile.lower)/2,
                       arm_length = rep(list_arms[i,3], nquantiles))
  df = rbind(df, newdata)
  rm(df_subset)
  rm(newdata)
}
df

#----------------------------------------------------------------------------#
# Compute the mean recombination rate in distance intervals
# Take a vector of indexes in the df dataframe
# and return a vector of rec.rate
# Standardized recombination rate among chromosomes (i.e. each chromosome standardized independently)
# for (i in 1:nrow(list_maps)) {
#   cat(list_maps$set[i], list_maps$chromosome[i], "\n")
#   df_dist2telomere$rec.rate[which(df_dist2telomere$set == list_maps$set[i] & df_dist2telomere$chromosome == list_maps$chromosome[i])] =
#     scale(df_dist2telomere$rec.rate[which(df_dist2telomere$set == list_maps$set[i] & df_dist2telomere$chromosome == list_maps$chromosome[i])])
# }
rec_in_distance_quantile_centromere = function(idx) {
  # subset = subset(df_dist2telomere, set == df$set[idx] & chromosome == df$chromosome[idx])
  # rec.rate.quantile = mean(subset$rec.rate, na.rm = TRUE)
  rec.rate.quantile = mean(df_dist2centromere$rec.rate[df_dist2centromere$set == df$set[idx]
                                                       & df_dist2centromere$chromosome == df$chromosome[idx]
                                                       & df_dist2centromere$dist2centromere >= df$quantile.lower[idx]
                                                       & df_dist2centromere$dist2centromere <= df$quantile.upper[idx]], na.rm = TRUE)
  return(rec.rate.quantile)
}

idx = c(1:nrow(df))
require(pbmcapply)
rec.rate = unlist(pbmclapply(idx, function(x) rec_in_distance_quantile_centromere(x), mc.cores =  7))
# hist(log(rec.rate))
df$rec.rate = rec.rate
hist(df$rec.rate)

# Scale recombination

# # Scale recombination rates per chromosome arms
# for (i in 1:nrow(list_arms)) {
#   cat(list_arms[i,1], list_arms[i,2], list_arms[i,3], "\n")
#   df$rec.rate[which(df$set == list_maps$set[i] &
#                       df$chromosome == list_maps$chromosome[i] &
#                       df$arm_length == list_arms$arm_length[i])] = as.numeric(scale(df$rec.rate[which(df$set == list_arms$set[i] &
#                                                                                                         df$chromosome == list_arms$chromosome[i] &
#                                                                                                         df$arm_length == list_arms$arm_length[i])]))
# }
df$rec.rate = as.numeric(scale(df$rec.rate))
hist(df$rec.rate)
mean(df$rec.rate, na.rm = TRUE)
var(df$rec.rate, na.rm = TRUE)

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
df$dist2centromere = rep(seq(0.025, 0.475, length.out = 20), nrow(df)/20)

#----------------------------------------------------------------------------#
write.table(df, "output/dist2centromere/DistancesRelativeScaled_bins.txt", row.names = FALSE, col.names = TRUE)

df = read.table("output/dist2centromere/DistancesRelativeScaled_bins.txt", header = TRUE)

summary(df)

#----------------------------------------------------------------------------#
# FIGURES ----
#----------------------------------------------------------------------------#
# Quantiles of rec. rate ~ distance to centromere for each chromosome
# one color per species (set)
DistancesRelative_bins = ggplot(data = df, aes(x = dist2centromere, y = rec.rate)) +
  geom_line(aes(group = interaction(set, chromosome, arm_length), colour = set)) +
  geom_point(aes(colour = set)) +
  xlab("Relative distance to centromere") + ylab("Recombination rate (cM/Mb)") +
  stat_smooth(aes(y = rec.rate), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Black", fill = "DarkGrey", alpha = 0.3) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
DistancesRelative_bins

# Unreadable, too much lines

#----------------------------------------------------------------------------#
# Split the dataset to identify groups/patterns ----
#----------------------------------------------------------------------------#
# Figures per classes of chromosome sizes?
# Indeed there seem to have co-existence of different patterns
hist(chromosome.stats$phys.map.length[which(!is.na(chromosome.stats$centromeric_index))], breaks = 60)
abline(v = 4e+08, col = "Red")
abline(v = 0.85e+08, col = "Red", lty = 2)
# Species with larger chromosomes
largerChr = chromosome.stats$set[which(chromosome.stats$phys.map.length > 4e+08)]
df_largerChr = subset(df, df$set %in% largerChr)
DistancesRelative_binsLargerChr = ggplot(data = df_largerChr, aes(x = dist2centromere, y = rec.rate)) +
  # geom_line(aes(group = interaction(set, chromosome), colour = set)) +
  geom_point(aes(colour = set, alpha = 0.3)) +
  stat_smooth(aes(y = rec.rate), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Black", fill = "DarkGrey", alpha = 0.3) +
  xlab("Relative distance to centromere") + ylab("Recombination rate (cM/Mb)") +
  xlim(0, 0.5) +
  # geom_abline(intercept = mod$coefficients[1], slope = mod$coefficients[2], linetype = "dashed") + # The linear regression line
  # geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
DistancesRelative_binsLargerChr

# Medium chromosomes
mediumChr = chromosome.stats$set[which(chromosome.stats$phys.map.length > 0.85e+08 & chromosome.stats$phys.map.length < 4e+08)]
df_mediumChr = subset(df, df$set %in% mediumChr)
DistancesRelative_binsMediumChr = ggplot(data = df_mediumChr, aes(x = dist2centromere, y = rec.rate)) +
  # geom_line(aes(group = interaction(set, chromosome), colour = set)) +
  geom_point(aes(colour = set, alpha = 0.3)) +
  stat_smooth(aes(y = rec.rate), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Black", fill = "DarkGrey", alpha = 0.3) +
  xlab("Relative distance to centromere") + ylab("Recombination rate (cM/Mb)") +
  xlim(0, 0.5) +
  # geom_abline(intercept = mod$coefficients[1], slope = mod$coefficients[2], linetype = "dashed") + # The linear regression line
  # geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
DistancesRelative_binsMediumChr
# And smaller chromosomes
smallerChr = chromosome.stats$set[which(chromosome.stats$phys.map.length < 0.85e+08)]
df_smallerChr = subset(df, (df$set %in% smallerChr))
DistancesRelative_binsSmallerChr = ggplot(data = df_smallerChr, aes(x = dist2centromere, y = rec.rate)) +
  # geom_line(aes(group = interaction(set, chromosome), colour = set)) +
  geom_point(aes(colour = set, alpha = 0.3)) +
  stat_smooth(aes(y = rec.rate), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Black", fill = "DarkGrey", alpha = 0.3) +
  xlab("Relative distance to centromere") + ylab("Recombination rate (cM/Mb)") +
  xlim(0, 0.5) +
  # geom_abline(intercept = mod$coefficients[1], slope = mod$coefficients[2], linetype = "dashed") + # The linear regression line
  # geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
DistancesRelative_binsSmallerChr

ggarrange(DistancesRelative_binsLargerChr, DistancesRelative_binsMediumChr, DistancesRelative_binsSmallerChr,
          ncol = 3)

#----------------------------------------------------------------------------#
# Chromosomes pooled per species ----
#----------------------------------------------------------------------------#
df_pooled = read.table("output/dist2centromere/DistancesRelativeScaled_bins.txt", header = TRUE)
# Quantiles of relative distances are pooled together
df_pooled$species = gsub("_[A-Za-z0-9]*$", "", df_pooled$set)
df_pooled$species = gsub("_MaizeGDBConsensus", "", df_pooled$species)
df_pooled$species = gsub("_", " ", df_pooled$species)

df_pooled = aggregate(df_pooled$rec.rate~df_pooled$dist2centromere+df_pooled$species+df_pooled$arm_length, FUN=mean)
colnames(df_pooled) = c("dist2centromere", "species", "arm_length", "rec.rate")

# LMM (model selection, model minimizing AIC)
hist(df_pooled$rec.rate)

lmer.model = lmer(rec.rate ~ dist2centromere + I(dist2centromere^2) + (1|species), data = df_pooled)
summary(lmer.model) # Quadratic coefficient is not significant
lmer.model = lmer(rec.rate ~ dist2centromere + (dist2centromere|species), data = df_pooled)
summary(lmer.model) # Model with random slopes

# Diagnostic plots
# Residuals distribution
library(car)
qqPlot(residuals(lmer.model))
# Residuals as a function of fitted
plot(fitted(lmer.model), residuals(lmer.model))
abline(h = 0)
# Observed~predicted qqplot
plot(predict(lmer.model), df_pooled$rec.rate, xlab = "Predicted",
     ylab = "Observed")
abline(a = 0, b = 1)
par(mfrow = c(1,1))

# Confidence interval of the coefficients was computed by bootstrap
# Parametric bootstrap
lmer.confint = confint(lmer.model, level = 0.95,
                       method = "boot", nsim = 1000, boot.type = "norm")

# Quadratic function and C.I.
quadrafun = function(x, c) {fixef(lmer.model)[[1]] + fixef(lmer.model)[[2]]*x}
lowerquadrafun = function(x, c) {lmer.confint[[5]] + lmer.confint[[6]]*x}
upperquadrafun = function(x, c) {lmer.confint[[11]] + lmer.confint[[12]]*x}
# The C.I. values
df.new = data.frame(dist2centromere = df_pooled$dist2centromere[1:20])
df.new$lwr.pred = lowerquadrafun(df.new$dist2centromere)
df.new$upr.pred = upperquadrafun(df.new$dist2centromere)

DistancesRelative_pooledChromosomes = ggplot(data = df_pooled, aes(x = dist2centromere, y = rec.rate)) +
  # geom_line(aes(group = species, colour = species)) +
  geom_point(aes(colour = species), alpha = 0.4) +
  xlab("Relative distance to centromere") + ylab("Recombination rate (cM/Mb)") +
  # stat_smooth(aes(y = rec.rate), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Black", fill = "DarkGrey", alpha = 0.3) +
  stat_function(fun = quadrafun, colour = "Black", size = 1) +
  geom_ribbon(data = df.new, aes(x = dist2centromere, ymin = lwr.pred, ymax = upr.pred), alpha = 0.4, inherit.aes = F, fill = "Black") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
DistancesRelative_pooledChromosomes

#----------------------------------------------------------------------------#
# Distribution of correlations between distance to telomere & recombination ----
#----------------------------------------------------------------------------#
# One correlation coefficient per chromosome
df_dist2centromere = read.table("output/dist2centromere/AllChromosomes_Dist2centromeres.txt", header = TRUE)
library(plyr)
# Trim for NA values
df_dist2centromere_trimmed = subset(df_dist2centromere,
                                    !(is.na(df_dist2centromere$rec.rate) | is.na(df_dist2centromere$dist2centromere)))
cor_dist2centromere = ddply(df_dist2centromere_trimmed, c("set", "chromosome"),
                            function(x) cor(x$rec.rate, x$dist2centromere, method = "spearman"))
colnames(cor_dist2centromere)[3] = "correlation"

hist(cor_dist2centromere$correlation, freq = TRUE,
     breaks = 20, xlab = "Correlation", main = "")





#-------------------------------------------#
# HYPOTHESIS TWO ----
# (2) asymmetry of the recombination landscape is driven by the centromere position ----

# If one CO mandatory per chromosom arm, then genetic length of short arm and long arm
# should be roughly equals
# It means the ratio short arm genetic length/total genetic length should be equal to 0.5
# Remove maps without centromere information
list_maps = data.frame(set = as.character(chromosome.stats$set[which(!is.na(chromosome.stats$centromeric_index))]),
                       chromosome = as.character(chromosome.stats$chromosome[which(!is.na(chromosome.stats$centromeric_index))]),
                       stringsAsFactors = FALSE)
short_arm_ratio = data.frame(set = character(nrow(list_maps)), chromosome = character(nrow(list_maps)))
short_arm_ratio$short_arm_genlength = NA
short_arm_ratio$total_genlength = NA
short_arm_ratio$short_arm_ratio = NA
for(i in 1:nrow(list_maps)) {
  cat(list_maps[i,1], list_maps[i,2], "\n")
  marey_map = read.table(paste("data-cleaned/marey_maps/", list_maps$set[i], ".txt", sep = ""), header = TRUE)
  marey_map = subset(marey_map, marey_map$map == list_maps$chromosome[i])
  marey_map$phys = marey_map$phys/1000000
  short_arm_ratio$set[i] = list_maps$set[i]
  short_arm_ratio$chromosome[i] = list_maps$chromosome[i]
  # Total genetic length
  short_arm_ratio$total_genlength[i] = max(marey_map$gen, na.rm = TRUE)
  # Centromeric breakpoint
  totalchrsize = max(marey_map$phys, na.rm = TRUE)
  centromeric_position = chromosome.stats$centromeric_index_position_oriented[which(chromosome.stats$set == list_maps[i,1] &
                                                                                      chromosome.stats$chromosome == list_maps[i,2])]
  # Identify short/long arm
  if (centromeric_position <= totalchrsize/2) {
    short_arm_ratio$short_arm_genlength[i] = max(marey_map$gen[marey_map$phys <= centromeric_position], na.rm = TRUE)
  } else {
    if (centromeric_position > totalchrsize/2) {
      short_arm_ratio$short_arm_genlength[i] = short_arm_ratio$total_genlength[i] - max(marey_map$gen[marey_map$phys <= centromeric_position], na.rm = TRUE)
    }
  }
}
short_arm_ratio$short_arm_ratio = short_arm_ratio$short_arm_genlength/short_arm_ratio$total_genlength

# CHROMOSOME LEVEL
# Distribution of the ratio per chromosome
hist(short_arm_ratio$short_arm_ratio, breaks = 40, main = "", xlab = "Short Arm Ratio")

# SPECIES LEVEL
# At a species level, is it a common pattern? i.e. the mean ratio is different to 0.5
# Test the hypothesis by randomly resampling chromosomes to estimate the mean ratio among all species
# short arm genetic length/total genetic length
# H0: The ratio short arm/total length is equal to 0.5 if there is
# the same amount of recombination on each chromosome arm
nboot = 1000
boot = numeric(nboot)
for (i in 1:nboot) {
  boot[i] = mean(sample(short_arm_ratio$short_arm_ratio, replace = TRUE), na.rm = TRUE)
}
mean(boot)
quantile(boot, 0.025)
quantile(boot, 0.975)

hist(short_arm_ratio$short_arm_ratio, breaks = 40, main = "", xlab = "Short Arm Ratio",
     ylim = c(0, 55))
abline(v = 0.5, lw = 2)
lines(density(boot), col = "red", lw = 2)
abline(v = mean(boot), col = "red", lw = 2)

# Illustration of the hypothesis,
# FIGURES - Return to the broken stick
# Re-plot the Broken Stick with all chromosomes oriented in the same way and a red dot to the centromere position

# See broken stick scripts





























#     ----..................................DEPRECATED ----

#============================================================================#
#
#   CENTROMERES: IS RECOMBINATION DRIVEN BY CENTROMERES? ----
#
#============================================================================#

# The hypothesis is that centromeres are interacting with recombination
# and suppressing/lowering recombination rates around them


#============================================================================#
#   Estimate Rec. rate ~ relative distance to centromere ----
#============================================================================#
chromosome.stats = read.table(file = "tables/chromosome.stats.csv", header = TRUE,
                              sep = ";", stringsAsFactors = FALSE)
# Use the centromere position given in the literature
# Measure of the ratio short arm/long arm
# Relative size of each arm

#-------------------------#
# NUMBER OF CHROMOSOMES WITH CENTROMERE POSITION
# Maps without metadata about the centromeric index
unique(chromosome.stats$set[is.na(chromosome.stats$centromeric_index)])
# Maps with a centromeric index
unique(chromosome.stats$set[!is.na(chromosome.stats$centromeric_index)])
length(unique(chromosome.stats$set[!is.na(chromosome.stats$centromeric_index)]))

# Compute distances to the centromere
# Split short arms and long arms (arm_length = short/long)
df_dist2centromere = data.frame(set = character(0), chromosome = character(0),
                                rec.rate = numeric(0), dist2centromere = numeric(0),
                                arm_length = factor(levels = c("short", "long")))

# Remove maps without centromere information
list_maps = data.frame(set = as.character(chromosome.stats$set[which(!is.na(chromosome.stats$centromeric_index))]),
                       chromosome = as.character(chromosome.stats$chromosome[which(!is.na(chromosome.stats$centromeric_index))]),
                       stringsAsFactors = FALSE)


for(i in 1:nrow(list_maps)) {
  cat(list_maps[i,1], list_maps[i,2], "\n")
  recombination_map = read.table(paste("output/recombination_maps/loess/100kbwind/", list_maps[i,1], "_chromosome", list_maps[i,2], ".txt", sep = ""), header = TRUE)
  recombination_map = recombination_map[order(recombination_map$phys),]
  # Physical position of the centromere
  totalchrsize = max(recombination_map$phys, na.rm = TRUE)
  centromeric_index = chromosome.stats$centromeric_index[which(chromosome.stats$set == list_maps[i,1] &
                                           chromosome.stats$chromosome == list_maps[i,2])]
  centromeric_position = totalchrsize*centromeric_index
  # Identify short/long arm
  recombination_map$arm_length = NA
  if (centromeric_index < 0.5) {
    recombination_map$arm_length[which(recombination_map$phys < centromeric_position)] = "short"
    recombination_map$arm_length[which(recombination_map$phys >= centromeric_position)] = "long"
    
  } else {
    if (centromeric_index >= 0.5) {
      recombination_map$arm_length[which(recombination_map$phys <= centromeric_position)] = "long"
      recombination_map$arm_length[which(recombination_map$phys > centromeric_position)] = "short"
    }
  }
  # Compute distances
  # Absolute value of the difference window position - centromere position
  recombination_map$dist2centromere = abs(centromeric_position - recombination_map$phys)
  # Results
  df = data.frame(set = rep(list_maps[i,1], nrow(recombination_map)), chromosome = rep(list_maps[i,2], nrow(recombination_map)),
                  rec.rate = recombination_map$rec.rate, dist2centromere = recombination_map$dist2centromere,
                  arm_length = recombination_map$arm_length)
  # remove older estimates
  df_dist2centromere = df_dist2centromere[!(df_dist2centromere$set == list_maps[i,1] & df_dist2centromere$chromosome == list_maps[i,2]),]
  # Append data
  df_dist2centromere = rbind(df_dist2centromere, df)
}
df_dist2centromere

write.table(df_dist2centromere, "output/dist2centromere/AllChromosomes_Dist2centromeres.txt", row.names = FALSE, col.names = TRUE)


#============================================================================#
#   Raw data - Rec. rate ~ distance to centromere ----
#============================================================================#
df_dist2centromere = read.table("output/dist2centromere/AllChromosomes_Dist2centromeres.txt", header = TRUE)

# Remove NA and transform null values
dist2centromere = subset(df_dist2centromere, (!is.na(df_dist2centromere$rec.rate)))
# dist2centromere$rec.rate = dist2centromere$rec.rate + 0.0001
# dist2centromere$dist2centromere = dist2centromere$dist2centromere + 0.0001
# Remove null values
dist2centromere = dist2centromere[(dist2centromere$rec.rate != 0),]
dist2centromere = dist2centromere[(dist2centromere$dist2centromere != 0),]
# Species names
dist2centromere = merge(dist2centromere, unique(chromosome.stats[,c(1, 24)]), by = "set", all.x = TRUE)

# Distribution of data
hist(dist2centromere$dist2centromere, breaks = 20, xlab = "Distance to centromere (Mb)", main = "")
hist(log(dist2centromere$dist2centromere), breaks = 20, xlab = "Distance to centroomere (log-transformed)", main = "")

hist(dist2centromere$rec.rate, breaks = 20, xlab = "Recombination rate (cM/Mb)", main = "")
hist(log(dist2centromere$rec.rate), breaks = 20, xlab = "Recombination rate (log-transformed)", main = "")


plot(log(dist2centromere$dist2centromere), log(dist2centromere$rec.rate),
     xlab = "Distance to centromere (log-transformed)",
     ylab = "Recombination rate (log-transformed)")

# A simple regression model
mod_dist2centromere = lm(log(rec.rate) ~ log(dist2centromere), data = dist2centromere)
abline(mod_dist2centromere, col = "Red")
# Linear Mixed Model
lmm_dist2centromere = lmer(log(rec.rate) ~ log(dist2centromere) + (log(dist2centromere) | species), data = dist2centromere)
summary(lmm_dist2centromere)
# Significance of coefficients
anova(lmm_dist2centromere)

# Mixed regression line
abline(a = fixef(lmm_dist2centromere)[1], b = fixef(lmm_dist2centromere)[2], col = "Blue")

# Correlation
cor.test(log(dist2centromere$dist2centromere), log(dist2centromere$rec.rate), method = "spearman")




#----------------------------------------------------------------------------#
df_dist2centromere = read.table("output/dist2centromere/AllChromosomes_Dist2centromeres.txt", header = TRUE)
# Compute the quantiles of distances
# Prepare the dataset with species and chromosomes to gather
# Subset short and long arms
nquantiles = 20
df = data.frame(set = character(0), chromosome = character(0),
                rec.rate = numeric(0), lower = numeric(0), upper = numeric(0),
                quantile.lower = numeric(0), quantile.upper = numeric(0),
                dist2centromere = numeric(0),
                arm_length = factor(levels = c("short", "long")))
# Compute the quantiles of distances for each species and chromosome
list_arms = data.frame(set = (rep(list_maps$set, each = 2)),
                       chromosome = rep(list_maps$chromosome, each = 2),
                       arm_length = rep(c("short", "long"), nrow(list_maps)),
                       stringsAsFactors = FALSE)

for (i in 1:nrow(list_arms)) {
  cat(list_arms[i,1], list_arms[i,2], list_arms[i,3], "\n")
  df_subset = df_dist2centromere[which(df_dist2centromere$set == list_arms[i,1] &
                                         df_dist2centromere$chromosome == list_arms[i,2] &
                                         df_dist2centromere$arm_length == list_arms[i,3]),]
  # Boundaries of 20 bins of equal distance
  bins = seq(from = 0, to = max(df_subset$dist2centromere, na.rm = TRUE), length.out = nquantiles+1)
  quantile.lower = bins[1:nquantiles]
  quantile.upper = bins[2:(nquantiles+1)]
  newdata = data.frame(set = rep(as.character(unique(df_subset$set)), each = nquantiles),
                       chromosome = rep(as.character(unique(df_subset$chromosome)), each = nquantiles),
                       rec.rate = NA, lower = NA, upper = NA,
                       quantile.lower = quantile.lower, quantile.upper = quantile.upper,
                       dist2centromere = (quantile.upper+quantile.lower)/2,
                       arm_length = rep(list_arms[i,3], nquantiles))
  df = rbind(df, newdata)
  rm(df_subset)
  rm(newdata)
}
df

#----------------------------------------------------------------------------#
# Compute the mean recombination rate in distance intervals
# Take a vector of indexes in the df dataframe
# and return a vector of rec.rate
# Standardized recombination rate among chromosomes (i.e. each chromosome standardized independently)
# for (i in 1:nrow(list_maps)) {
#   cat(list_maps$set[i], list_maps$chromosome[i], "\n")
#   df_dist2telomere$rec.rate[which(df_dist2telomere$set == list_maps$set[i] & df_dist2telomere$chromosome == list_maps$chromosome[i])] =
#     scale(df_dist2telomere$rec.rate[which(df_dist2telomere$set == list_maps$set[i] & df_dist2telomere$chromosome == list_maps$chromosome[i])])
# }
rec_in_distance_quantile_centromere = function(idx) {
  # subset = subset(df_dist2telomere, set == df$set[idx] & chromosome == df$chromosome[idx])
  # rec.rate.quantile = mean(subset$rec.rate, na.rm = TRUE)
  rec.rate.quantile = mean(df_dist2centromere$rec.rate[df_dist2centromere$set == df$set[idx]
                                                     & df_dist2centromere$chromosome == df$chromosome[idx]
                                                     & df_dist2centromere$dist2centromere >= df$quantile.lower[idx]
                                                     & df_dist2centromere$dist2centromere <= df$quantile.upper[idx]], na.rm = TRUE)
  return(rec.rate.quantile)
}

idx = c(1:nrow(df))
require(pbmcapply)
rec.rate = unlist(pbmclapply(idx, function(x) rec_in_distance_quantile_centromere(x), mc.cores =  7))
# hist(log(rec.rate))
df$rec.rate = rec.rate
hist(df$rec.rate)

# Scale recombination rates per chromosome arms
for (i in 1:nrow(list_arms)) {
  cat(list_arms[i,1], list_arms[i,2], list_arms[i,3], "\n")
  df$rec.rate[which(df$set == list_maps$set[i] &
                      df$chromosome == list_maps$chromosome[i] &
                      df$arm_length == list_arms$arm_length[i])] = as.numeric(scale(df$rec.rate[which(df$set == list_arms$set[i] &
                        df$chromosome == list_arms$chromosome[i] &
                        df$arm_length == list_arms$arm_length[i])]))
}
# df$rec.rate = as.numeric(scale(df$rec.rate))
hist(df$rec.rate)
mean(df$rec.rate, na.rm = TRUE)
var(df$rec.rate, na.rm = TRUE)

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
df$dist2centromere = rep(seq(0.025, 0.475, length.out = 20), nrow(df)/20)

#----------------------------------------------------------------------------#
write.table(df, "output/dist2centromere/DistancesRelativeScaled_bins.txt", row.names = FALSE, col.names = TRUE)

df = read.table("output/dist2centromere/DistancesRelativeScaled_bins.txt", header = TRUE)

summary(df)

#----------------------------------------------------------------------------#
# FIGURES ----
#----------------------------------------------------------------------------#
# Quantiles of rec. rate ~ distance to centromere for each chromosome
# one color per species (set)
DistancesRelative_bins = ggplot(data = df, aes(x = dist2centromere, y = rec.rate)) +
  geom_line(aes(group = interaction(set, chromosome, arm_length), colour = set)) +
  geom_point(aes(colour = set)) +
  xlab("Relative distance to centromere") + ylab("Recombination rate (cM/Mb)") +
  stat_smooth(aes(y = rec.rate), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Black", fill = "DarkGrey", alpha = 0.3) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
DistancesRelative_bins

# Unreadable, too much lines

#----------------------------------------------------------------------------#
# Split the dataset to identify groups/patterns ----
#----------------------------------------------------------------------------#
# Figures per classes of chromosome sizes?
# Indeed there seem to have co-existence of different patterns
hist(chromosome.stats$phys.map.length[which(!is.na(chromosome.stats$centromeric_index))], breaks = 60)
abline(v = 4e+08, col = "Red")
abline(v = 0.85e+08, col = "Red", lty = 2)
# Species with larger chromosomes
largerChr = chromosome.stats$set[which(chromosome.stats$phys.map.length > 4e+08)]
df_largerChr = subset(df, df$set %in% largerChr)
DistancesRelative_binsLargerChr = ggplot(data = df_largerChr, aes(x = dist2centromere, y = rec.rate)) +
  # geom_line(aes(group = interaction(set, chromosome), colour = set)) +
  geom_point(aes(colour = set, alpha = 0.3)) +
  stat_smooth(aes(y = rec.rate), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Black", fill = "DarkGrey", alpha = 0.3) +
  xlab("Relative distance to centromere") + ylab("Recombination rate (cM/Mb)") +
  xlim(0, 0.5) +
  # geom_abline(intercept = mod$coefficients[1], slope = mod$coefficients[2], linetype = "dashed") + # The linear regression line
  # geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
DistancesRelative_binsLargerChr

# Medium chromosomes
mediumChr = chromosome.stats$set[which(chromosome.stats$phys.map.length > 0.85e+08 & chromosome.stats$phys.map.length < 4e+08)]
df_mediumChr = subset(df, df$set %in% mediumChr)
DistancesRelative_binsMediumChr = ggplot(data = df_mediumChr, aes(x = dist2centromere, y = rec.rate)) +
  # geom_line(aes(group = interaction(set, chromosome), colour = set)) +
  geom_point(aes(colour = set, alpha = 0.3)) +
  stat_smooth(aes(y = rec.rate), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Black", fill = "DarkGrey", alpha = 0.3) +
  xlab("Relative distance to centromere") + ylab("Recombination rate (cM/Mb)") +
  xlim(0, 0.5) +
  # geom_abline(intercept = mod$coefficients[1], slope = mod$coefficients[2], linetype = "dashed") + # The linear regression line
  # geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
DistancesRelative_binsMediumChr
# And smaller chromosomes
smallerChr = chromosome.stats$set[which(chromosome.stats$phys.map.length < 0.85e+08)]
df_smallerChr = subset(df, (df$set %in% smallerChr))
DistancesRelative_binsSmallerChr = ggplot(data = df_smallerChr, aes(x = dist2centromere, y = rec.rate)) +
  # geom_line(aes(group = interaction(set, chromosome), colour = set)) +
  geom_point(aes(colour = set, alpha = 0.3)) +
  stat_smooth(aes(y = rec.rate), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Black", fill = "DarkGrey", alpha = 0.3) +
  xlab("Relative distance to centromere") + ylab("Recombination rate (cM/Mb)") +
  xlim(0, 0.5) +
  # geom_abline(intercept = mod$coefficients[1], slope = mod$coefficients[2], linetype = "dashed") + # The linear regression line
  # geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
DistancesRelative_binsSmallerChr

ggarrange(DistancesRelative_binsLargerChr, DistancesRelative_binsMediumChr, DistancesRelative_binsSmallerChr,
          ncol = 3)

#----------------------------------------------------------------------------#
# Chromosomes pooled per species ----
#----------------------------------------------------------------------------#
df_pooled = read.table("output/dist2centromere/DistancesRelativeScaled_bins.txt", header = TRUE)
# Quantiles of relative distances are pooled together
df_pooled$species = gsub("_[A-Za-z0-9]*$", "", df_pooled$set)
df_pooled$species = gsub("_MaizeGDBConsensus", "", df_pooled$species)
df_pooled$species = gsub("_", " ", df_pooled$species)

df_pooled = aggregate(df_pooled$rec.rate~df_pooled$dist2centromere+df_pooled$species+df_pooled$arm_length, FUN=mean)
colnames(df_pooled) = c("dist2centromere", "species", "arm_length", "rec.rate")

# Quadratic Phylogenetic Regression
# Load the phylogenetic tree estimated in 'Phylogeny.R'
load("output/phylogeny/tree.Rdata")
tree
tree$tip.label
# Drop unsampled tips
library(ape)
tree = drop.tip(tree, tip = tree$tip.label[!(tree$tip.label %in% df_pooled$species)])
tree

# PGLMM (model selection, model minimizing AIC)
hist(df_pooled$rec.rate)
library(phyr)
# pglmm.model = pglmm(rec.rate ~ dist2telomere + (1|species), data = df_pooled, family = 'gaussian',
#                     cov_ranef = list(species = tree))
pglmm.model = pglmm(rec.rate ~ dist2centromere + I(dist2centromere^2) + (1|species), data = df_pooled, family = 'gaussian',
                    cov_ranef = list(species = tree))
summary(pglmm.model)

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
pglmm.model = pglmm(rec.rate ~ dist2centromere + I(dist2centromere^2) + (1|species), data = df_pooled, family = 'gaussian',
                    cov_ranef = list(species = tree), bayes = TRUE)
summary(pglmm.model)
# Quadratic function and C.I.
quadrafun = function(x, c) fixef(pglmm.model)[[1]][1] + fixef(pglmm.model)[[1]][2]*x + fixef(pglmm.model)[[1]][3]*x^2
lowerquadrafun = function(x, c) fixef(pglmm.model)[[2]][1] + fixef(pglmm.model)[[2]][2]*x + fixef(pglmm.model)[[2]][3]*x^2
upperquadrafun = function(x, c) fixef(pglmm.model)[[3]][1] + fixef(pglmm.model)[[3]][2]*x + fixef(pglmm.model)[[3]][3]*x^2
# The C.I. values
df.new = data.frame(dist2centromere = df_pooled$dist2centromere[1:20])
df.new$lwr.pred = lowerquadrafun(df.new$dist2centromere)
df.new$upr.pred = upperquadrafun(df.new$dist2centromere)

DistancesRelative_pooledChromosomes = ggplot(data = df_pooled, aes(x = dist2centromere, y = rec.rate)) +
  # geom_line(aes(group = species, colour = species)) +
  geom_point(aes(colour = species), alpha = 0.4) +
  xlab("Relative distance to centromere") + ylab("Recombination rate (cM/Mb)") +
  # stat_smooth(aes(y = rec.rate), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Black", fill = "DarkGrey", alpha = 0.3) +
  stat_function(fun = quadrafun, colour = "Black", size = 1) +
  geom_ribbon(data = df.new, aes(x = dist2centromere, ymin = lwr.pred, ymax = upr.pred), alpha = 0.4, inherit.aes = F, fill = "Black") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
DistancesRelative_pooledChromosomes

#----------------------------------------------------------------------------#
# Distribution of correlations between distance to telomere & recombination ----
#----------------------------------------------------------------------------#
# One correlation coefficient per chromosome
df_dist2centromere = read.table("output/dist2centromere/AllChromosomes_Dist2centromeres.txt", header = TRUE)
library(plyr)
# Trim for NA values
df_dist2centromere_trimmed = subset(df_dist2centromere,
                                  !(is.na(df_dist2centromere$rec.rate) | is.na(df_dist2centromere$dist2centromere)))
cor_dist2centromere = ddply(df_dist2centromere_trimmed, c("set", "chromosome"),
                          function(x) cor(x$rec.rate, x$dist2centromere, method = "spearman"))
colnames(cor_dist2centromere)[3] = "correlation"

hist(cor_dist2centromere$correlation, freq = TRUE,
     breaks = 20, xlab = "Correlation", main = "")


#============================================================================#
#
#   DISENTANGLE EFFECTS OF CENTROMERE AND TELOMERE ----
#
#============================================================================#

dist2telomere = read.table("output/dist2telomere/AllChromosomes_Dist2telomeres.txt", header = TRUE)
dist2centromere = read.table("output/dist2centromere/AllChromosomes_Dist2centromeres.txt", header = TRUE)
#============================================================================#
# Joint estimates of dist2telomere and dist2centromere ----
#============================================================================#
#-------------------------#
# NUMBER OF CHROMOSOMES WITH CENTROMERE POSITION
# Maps without metadata about the centromeric index
unique(chromosome.stats$set[is.na(chromosome.stats$centromeric_index)])
# Maps with a centromeric index
unique(chromosome.stats$set[!is.na(chromosome.stats$centromeric_index)])
length(unique(chromosome.stats$set[!is.na(chromosome.stats$centromeric_index)]))

# Compute distances to the centromere
# Split short arms and long arms (arm_length = short/long)
df_dist2 = data.frame(set = character(0), chromosome = character(0),
                                rec.rate = numeric(0), dist2telomere = numeric(0),
                                dist2centromere = numeric(0),
                                arm_length = factor(levels = c("short", "long")))

# Remove maps without centromere information
list_maps = data.frame(set = as.character(chromosome.stats$set[which(!is.na(chromosome.stats$centromeric_index))]),
                       chromosome = as.character(chromosome.stats$chromosome[which(!is.na(chromosome.stats$centromeric_index))]),
                       stringsAsFactors = FALSE)


for(i in 1:nrow(list_maps)) {
  cat(list_maps[i,1], list_maps[i,2], "\n")
  recombination_map = read.table(paste("output/recombination_maps/loess/100kbwind/", list_maps[i,1], "_chromosome", list_maps[i,2], ".txt", sep = ""), header = TRUE)
  recombination_map = recombination_map[order(recombination_map$phys),]
  
  # Estimate distance to telomere
  recombination_map$dist2telomere = apply(data.frame(recombination_map$phys, abs(recombination_map$phys-totalchrsize)), 1, min)
  
  # Physical position of the centromere
  totalchrsize = max(recombination_map$phys, na.rm = TRUE)
  centromeric_index = chromosome.stats$centromeric_index[which(chromosome.stats$set == list_maps[i,1] &
                                                                 chromosome.stats$chromosome == list_maps[i,2])]
  centromeric_position = totalchrsize*centromeric_index
  # Identify short/long arm
  recombination_map$arm_length = NA
  if (centromeric_index < 0.5) {
    recombination_map$arm_length[which(recombination_map$phys < centromeric_position)] = "short"
    recombination_map$arm_length[which(recombination_map$phys >= centromeric_position)] = "long"
    
  } else {
    if (centromeric_index >= 0.5) {
      recombination_map$arm_length[which(recombination_map$phys <= centromeric_position)] = "long"
      recombination_map$arm_length[which(recombination_map$phys > centromeric_position)] = "short"
    }
  }
  # Compute distances
  # Absolute value of the difference window position - centromere position
  recombination_map$dist2centromere = abs(centromeric_position - recombination_map$phys)
  # Results
  df = data.frame(set = rep(list_maps[i,1], nrow(recombination_map)), chromosome = rep(list_maps[i,2], nrow(recombination_map)),
                  rec.rate = recombination_map$rec.rate, dist2telomere = recombination_map$dist2telomere,
                  dist2centromere = recombination_map$dist2centromere,
                  arm_length = recombination_map$arm_length)
  # remove older estimates
  df_dist2 = df_dist2[!(df_dist2$set == list_maps[i,1] & df_dist2$chromosome == list_maps[i,2]),]
  # Append data
  df_dist2 = rbind(df_dist2, df)
}

write.table(df_dist2, "output/AllChromosomes_Dist2.txt", row.names = FALSE, col.names = TRUE)



#============================================================================#
# Which effect is stronger between centromeres and telomeres? ----
#============================================================================#
dist2 = read.table("output/AllChromosomes_Dist2.txt", header = TRUE)
# Remove NA and transform null values
dist2 = subset(dist2, (!is.na(dist2$rec.rate)))
# dist2centromere$rec.rate = dist2centromere$rec.rate + 0.0001
# dist2centromere$dist2centromere = dist2centromere$dist2centromere + 0.0001
# Remove null values
dist2 = dist2[(dist2$rec.rate != 0),]
dist2 = dist2[(dist2$dist2telomere != 0),]
dist2 = dist2[(dist2$dist2centromere != 0),]
# Species names
dist2 = merge(dist2, unique(chromosome.stats[,c(1, 24)]), by = "set", all.x = TRUE)

# Full set model with both variables and decomposition of variance
lmm_dist2 = lmer(log(rec.rate) ~ log(dist2telomere) + log(dist2centromere) + 
                                     (1 | species), data = dist2)
summary(lmm_dist2)
# Both effects are significant in the full set model

library(partR2)
R2_positive = partR2(lmm_dist2,
                     partvars = c("log(dist2telomere)", "log(dist2centromere)"), data = dist2,
                     nboot = 1000, CI = 0.95)
summary(R2_positive)
# Save R2 decomposition of variance
save(R2_positive, file = "output/models/lmm_dist2.Rda")

load(file = "output/models/lmm_dist2.Rda")
summary(R2_positive)# Part R2 is the variance explained by each individual predictor, by iterative removal of each predictor
forestplot(R2_positive, type = "R2", line_size = 0.7, text_size = 14, point_size = 3)
# Check also inclusive R2 (gives information about the structure of variance explained; square of the structure coefficient, i.e. the correlation between predictors)
# and beta weights (standardized regression slopes)
forestplot(R2_positive, type = "IR2", line_size = 0.7, text_size = 14, point_size = 3)
forestplot(R2_positive, type = "BW", line_size = 0.7, text_size = 14, point_size = 3)
# Structure coefficients show the correlation between individual predictor variables
forestplot(R2_positive, type = "SC", line_size = 0.7, text_size = 14, point_size = 3)





#============================================================================#
# Patterns are different between metacentric and acrocentric chromosomes? ----
#============================================================================#
df = read.table("output/dist2centromere/DistancesRelativeScaled_bins.txt", header = TRUE)

hist(chromosome.stats$centromeric_index, breaks = 20)
# Few acro-telo-centric chromosomes compared to meta-centric ones.
sum(chromosome.stats$centromeric_index > 0.3, na.rm = TRUE)
sum(chromosome.stats$centromeric_index < 0.3, na.rm = TRUE)

# Species with larger chromosomes
metacentric = chromosome.stats$set[which(chromosome.stats$centromeric_index >= 0.3)]
df_metacentric = subset(df, df$set %in% metacentric)
df_metacentric = subset(df_metacentric, df_metacentric$arm_length == "long")
DistancesRelative_metacentric = ggplot(data = df_metacentric, aes(x = dist2centromere, y = rec.rate)) +
  # geom_line(aes(group = interaction(set, chromosome), colour = set)) +
  geom_point(aes(colour = set, alpha = 0.3)) +
  stat_smooth(aes(y = rec.rate), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Black", fill = "DarkGrey", alpha = 0.3) +
  xlab("Relative distance to centromere") + ylab("Recombination rate (cM/Mb)") +
  xlim(0, 0.5) +
  ylim(0, 25) +
  # geom_abline(intercept = mod$coefficients[1], slope = mod$coefficients[2], linetype = "dashed") + # The linear regression line
  # geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
DistancesRelative_metacentric

# Acro-telo-centric chromosomes
acrotelocentric = chromosome.stats$set[which(chromosome.stats$centromeric_index < 0.3)]
df_acrotelocentric = subset(df, df$set %in% acrotelocentric)
df_acrotelocentric = subset(df_acrotelocentric, df_acrotelocentric$arm_length == "long")
DistancesRelative_acrotelocentric = ggplot(data = df_acrotelocentric, aes(x = dist2centromere, y = rec.rate)) +
  # geom_line(aes(group = interaction(set, chromosome), colour = set)) +
  geom_point(aes(colour = set, alpha = 0.3)) +
  stat_smooth(aes(y = rec.rate), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Black", fill = "DarkGrey", alpha = 0.3) +
  xlab("Relative distance to centromere") + ylab("Recombination rate (cM/Mb)") +
  xlim(0, 0.5) +
  ylim(0, 25) +
  # geom_abline(intercept = mod$coefficients[1], slope = mod$coefficients[2], linetype = "dashed") + # The linear regression line
  # geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
DistancesRelative_acrotelocentric

ggarrange(DistancesRelative_metacentric, DistancesRelative_acrotelocentric, ncol = 2)





































#============================================================================#
# DEPRECATED ----
#============================================================================#



#============================================================================#
# Finding centromeres ----
#============================================================================#

# There exist many ways to decribe a centromere and get its localization
# Centromere linked proteins may not be in the centromere, leading to false positions
# But some histone variants are characetistics of centromeric regions
# e.g. CENH3 or CENP-A for wheat (Sun et al. 2019)



#-----------------------------------------------------------------------------#
# Find in annotations (gff file, if exists) ----
#-----------------------------------------------------------------------------#

for (species in list_species) {
  print(species)
  list_accession = system(paste("ls ", wd, "/data/Genome/", tolower(species), sep = ""), intern = TRUE)
  for (accession in list_accession) {
    # Count the number of chromosomes
    list.chr = system(paste("ls ", wd, "/data/Genome/", tolower(species), "/", accession, "/gff3 | grep 'chromosome.[0-9]*.gff3'", sep = ""), intern = TRUE)
    name.chr = gsub(".gff3", "", gsub("chromosome.", "", regmatches(list.chr, gregexpr("chromosome.[0-9]*.gff3", list.chr))))
    rm(list.chr)
    # Then retrieve centromere of each chromosome
    for (chromosome in name.chr) {
      # Open gff files and find positions in a shell script for speed
      res = system(paste("sh ", wd, "/sources/getCentromerePosition.sh ", tolower(species), " ", accession, " ", chromosome, sep = ""), intern = TRUE)
      print(res)
      if (length(res) > 0) {
        link = c()
        id = c()
        chr = c()
        start = c()
        end = c()
        for (e in 1:length(res)) {
          # print(e)
          link = c(link, strsplit(res, " ")[[e]][1])
          id = c(id, strsplit(res, " ")[[e]][5])
          chr = c(chr, strsplit(res, " ")[[e]][2])
          start = c(start, strsplit(res, " ")[[e]][3])
          end = c(end, strsplit(res, " ")[[e]][4])
        }
        df = data.frame(link = link, species = species, accession = accession, id = id,
                        chr = chr, start = start, end = end)
        # df = data.frame(link = strsplit(res, " ")[[1]][1], species = species, accession = accession, id = strsplit(res, " ")[[1]][5],
        #                 chr = strsplit(res, " ")[[1]][2], start = strsplit(res, " ")[[1]][3], end = strsplit(res, " ")[[1]][4])
        # Append new centromere positions to genome metadata, if non existing
        tmp = read.table(file = "data/Genome/genome_metadata.txt", header = TRUE, sep = ";")
        df = df[!(df$id %in% tmp$id) | is.na(df$id),]
        # print(df)
        type = rep("centromere", nrow(df))
        tmp = rbind(tmp, cbind(df, type))
        write.table(tmp, file = "data/Genome/genome_metadata.txt", sep = ";", row.names = F, col.names = T, quote = FALSE)
        rm(type)
        rm(tmp)
        rm(df)
      }
      
    }
  }
}

genome_metadata = read.table(file = "data/Genome/genome_metadata.txt", header = TRUE, sep = ";")
genome_metadata = genome_metadata[!(duplicated(genome_metadata)),]
write.table(genome_metadata, file = "data/Genome/genome_metadata.txt", sep = ";", row.names = F, col.names = T, quote = FALSE)

#-----------------------------------------------------------------------------#
# Find in proteins with Entrez & BLAST ----
#-----------------------------------------------------------------------------#

# species = "Arabidopsis_thaliana"
# accession = "GCA_902460315.1"

for (species in list_species) {
  print(species)
  list_accession = system(paste("ls ", wd, "/data/Genome/", tolower(species), sep = ""), intern = TRUE)
  for (accession in list_accession) {
    # get centromere proteins available at NCBI, for all species of interest
    # species identified by species name (ORGNANISM) in the genome DB
    # entrez_dbs()
    # entrez_db_summary("protein")
    # entrez_db_searchable("protein")
    search = c()
    r_search = entrez_search(db = "protein", term = paste(gsub("_", " ", species), "[ORGN] AND kinetochore[TITL]", sep = ""))
    search = c(search, r_search$ids)
    r_search = entrez_search(db = "protein", term = paste(gsub("_", " ", species), "[ORGN] AND centromer[TITL]", sep = ""))
    search = c(search, r_search$ids)
    r_links = entrez_link(dbfrom = 'protein', id = search, db='all')
    search = unique(search)
    
    # A DB with eight columns
    # link
    # Species
    # Genome accession
    # ID of the sequence
    # chromosome
    # Sequence (nuc)
    # Start position
    # End position
    size = length(r_links$links$protein_nuccore)
    df = data.frame(link = r_links$links$protein_nuccore, species = rep(species, size), accession = rep(accession, size), id = rep(NA, size), chr = rep(NA, size),
                    start = rep(NA, size), end = rep(NA, size))
    
    for (s in 1:size) {
      cat("Retrieving element", s, "of", size, "\n")
      # Retrieve the centromere sequence in nucleotides
      seq = entrez_fetch(db="nuccore", id = df$link[s], rettyp = "fasta")
      # Remove the first line of fasta (id line)
      id = gsub(">", "", strsplit(unlist(seq), "\n")[[1]][1])
      print(id)
      seq = strsplit(unlist(seq), "\n")[[1]][-1]
      seq = as.character(paste(seq, sep = "", collapse = ""))
      
      # Run blast on a local db in a repertory 'species/accession/fasta/'
      ## load a BLAST database (replace db with the location + name of the BLAST DB)
      ls = list.files(path = paste(wd, "/data/Genome/", tolower(species),"/" ,accession,"/fasta/", sep = ""))
      dbname = gsub(".nin", "", ls[grepl(".nin", ls)])
      #dbname = gsub(".gz", "",list.files(path = paste(wd, "/data/Genome/", tolower(species),"/" ,accession,"/fasta/", sep = ""))) # Older way to find the db name, not working with Ensembl db
      blastDB = blast(db = paste(wd, "/data/Genome/", tolower(species),"/" ,accession,"/fasta/", dbname, sep = ""), type = "blastn")
      # blastDB
      query = DNAStringSet(x = seq, use.names = TRUE) # Convert the sequence tag as a XStringSet
      hit = predict(blastDB, query) # Search for the primer position in the genome DB
      # print(hit)
      
      # Recover in the BlastDB the Chromosome number from the SubjectID
      res = system(paste("blastdbcmd -db ", wd, "/data/Genome/", tolower(species), "/" , accession,"/fasta/", dbname, " -entry ", as.character(hit$SubjectID[1]), sep = ""), intern = TRUE)
      res[1]
      map = gsub("chromosome:*[ ]*[A-Za-z]*", "",regmatches(res[1], regexpr("chromosome:*[ ]*[A-Za-z]*[0-9]+", res[1])))
      print(map)
      
      # save in df
      df$id[s] = as.character(id)
      df$chr[s] = as.character(map)
      # df$seq[s] = as.character(seq)
      # Positions of the best hit
      df$start[s] = as.numeric(hit$S.start[1])
      df$end[s] = as.numeric(hit$S.end[1])
    }
    rm(size)
    
    # Remove lines non DNA (e.g. mRNA)
    df = df[-grep("mRNA", df$id),]
    # Remove duplicates
    df = df[!duplicated(df$chr),]
    df
    
    # Append new entrez links to genome metadata
    tmp = read.csv(file = "data/Genome/genome_metadata.txt", header = TRUE, sep = ";")
    df = df[!(df$link %in% tmp$link),]
    type = rep("centromere", nrow(df))
    tmp = rbind(tmp, cbind(df, type))
    tmp
    write.table(tmp, file = "data/Genome/genome_metadata.txt", sep = ";", row.names = F, col.names = T, quote = FALSE)
    rm(type)
    rm(tmp)
    rm(df)
  }
}


# # get centromere proteins available at NCBI, for all species of interest
# # species identified by species name (ORGNANISM) in the genome DB
# # entrez_dbs()
# # entrez_db_summary("protein")
# # entrez_db_searchable("protein")
# search = c()
# r_search = entrez_search(db = "protein", term = paste(gsub("_", " ", species), "[ORGN] AND kinetochore[TITL]", sep = ""))
# search = c(search, r_search$ids)
# r_search = entrez_search(db = "protein", term = paste(gsub("_", " ", species), "[ORGN] AND centromer[TITL]", sep = ""))
# search = c(search, r_search$ids)
# r_links = entrez_link(dbfrom = 'protein', id = search, db='all')
# search = unique(search)
# 
# # A DB with eight columns
# # link
# # Species
# # Genome accession
# # ID of the sequence
# # chromosome
# # Sequence (nuc)
# # Start position
# # End position
# size = length(r_links$links$protein_nuccore)
# df = data.frame(link = r_links$links$protein_nuccore, species = rep(species, size), accession = rep(accession, size), id = rep(NA, size), chr = rep(NA, size),
#                 start = rep(NA, size), end = rep(NA, size))
# 
# for (s in 1:size) {
#   cat("Retrieving element", s, "of", size, "\n")
#   # Retrieve the centromere sequence in nucleotides
#   seq = entrez_fetch(db="nuccore", id = df$link[s], rettyp = "fasta")
#   # Remove the first line of fasta (id line)
#   id = gsub(">", "", strsplit(unlist(seq), "\n")[[1]][1])
#   print(id)
#   seq = strsplit(unlist(seq), "\n")[[1]][-1]
#   seq = as.character(paste(seq, sep = "", collapse = ""))
#   
#   # Run blast on a local db in a repertory 'species/accession/fasta/'
#   ## load a BLAST database (replace db with the location + name of the BLAST DB)
#   ls = list.files(path = paste(wd, "/data/Genome/", tolower(species),"/" ,accession,"/fasta/", sep = ""))
#   dbname = gsub(".nin", "", ls[grepl(".nin", ls)])
#   #dbname = gsub(".gz", "",list.files(path = paste(wd, "/data/Genome/", tolower(species),"/" ,accession,"/fasta/", sep = ""))) # Older way to find the db name, not working with Ensembl db
#   blastDB = blast(db = paste(wd, "/data/Genome/", tolower(species),"/" ,accession,"/fasta/", dbname, sep = ""), type = "blastn")
#   # blastDB
#   query = DNAStringSet(x = seq, use.names = TRUE) # Convert the sequence tag as a XStringSet
#   hit = predict(blastDB, query) # Search for the primer position in the genome DB
#   # print(hit)
#   
#   # Recover in the BlastDB the Chromosome number from the SubjectID
#   res = system(paste("blastdbcmd -db ", wd, "/data/Genome/", tolower(species), "/" , accession,"/fasta/", dbname, " -entry ", as.character(hit$SubjectID[1]), sep = ""), intern = TRUE)
#   res[1]
#   map = gsub("chromosome:*[ ]*[A-Za-z]*", "",regmatches(res[1], regexpr("chromosome:*[ ]*[A-Za-z]*[0-9]+", res[1])))
#   print(map)
#   
#   # save in df
#   df$id[s] = as.character(id)
#   df$chr[s] = as.character(map)
#   # df$seq[s] = as.character(seq)
#   # Positions of the best hit
#   df$start[s] = as.numeric(hit$S.start[1])
#   df$end[s] = as.numeric(hit$S.end[1])
# }
# rm(size)
# 
# # Remove lines non DNA (e.g. mRNA)
# df = df[-grep("mRNA", df$id),]
# # Remove duplicates
# df = df[!duplicated(df$chr),]
# df
# 
# # Append new entrez links to genome metadata
# tmp = read.csv(file = "data/Genome/genome_metadata.csv", header = TRUE, sep = ";")
# df = df[!(df$link %in% tmp$link),]
# type = rep("centromere", nrow(df))
# tmp = rbind(tmp, cbind(df, type))
# tmp
# write.table(tmp, file = "data/Genome/genome_metadata.csv", sep = ";", row.names = F, col.names = T, quote = FALSE)
# rm(type)
# rm(tmp)
# rm(df)

#============================================================================#
# END ----
#============================================================================#
