########################################################################### #
#                    ECOBIO - PhD
#
#       Creating recombination landscapes with Marey Map
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

metadata = read.csv(file = paste(wd, "/data/Genetic_maps/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")
metadata.clean = read.csv(file = paste(wd, "/data-cleaned/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")
# Load the joint database with all maps
data.final = read.table(file = "data-cleaned/marey_maps/AllMaps.txt", header = TRUE, sep = "\t")
chromosome.stats = read.table(file = "tables/chromosome.stats.csv", header = TRUE, sep = ";")


#============================================================================#
# STATISTICS OF ESTIMATED RECOMBINATION MAPS ----
#============================================================================#
# Marey map selection step
# A function returning a data frame with summary statistics of selected maps
# AND saving a file with summary statistics
source("sources/MareyMap.R")
chromosome.stats = recombination.map.statistics(data.final, save.file = TRUE)

# Check congruence between chromosome.stats and the list of recombination maps in 'loess/'
chrlist = paste(chromosome.stats$set, "_chromosome", chromosome.stats$chromosome, sep = "")
nrow(chromosome.stats)
# List of recombinatin map files
list = system(paste("ls ", wd, "/output/recombination_maps/loess/100kbwind/", sep = ""), intern = TRUE)
list = list[grep("*.txt", list)]
list = gsub(".txt", "", list)
length(list)
list
list[which(!(as.character(list) %in% chrlist))]

#============================================================================#
# Exploratory analyses ----
#============================================================================#

# Number of dataset
length(unique(chromosome.stats$set))
unique(chromosome.stats$set)
# Number of species
length(unique(regmatches(as.character(chromosome.stats$set), gregexpr("^[A-Z][a-z]+_[a-z]+", as.character(chromosome.stats$set)))))
unlist(unique(regmatches(as.character(chromosome.stats$set), gregexpr("^[A-Z][a-z]+_[a-z]+", as.character(chromosome.stats$set)))))
# Number of chromosomes
nrow(chromosome.stats)

# Describe maps
summary(metadata.clean$n_progeny)
summary(as.numeric(table(chromosome.stats$species)))

summary(chromosome.stats$nb.markers)
hist(chromosome.stats$nb.markers, breaks = 40)
summary(chromosome.stats$density.markers.bp*1000000)
hist(chromosome.stats$density.markers.bp*1000000, breaks = 40)
sum(chromosome.stats$density.markers.bp*1000000 < 0.5)
summary(chromosome.stats$density.markers.cM)
hist(chromosome.stats$density.markers.cM, breaks = 40)
sum(chromosome.stats$density.markers.cM < 1)
summary(chromosome.stats$largest.gap.cM)
hist(chromosome.stats$largest.gap.cM, breaks = 40)

chromosome.stats$set[which(chromosome.stats$nb.markers < 40)]

summary(chromosome.stats$linkage.map.length.correctedHW)
summary(chromosome.stats$phys.map.length/1000000)


summary(chromosome.stats$mean.recrate)
# Mean chromosome recombination rate by dataset
boxplot(chromosome.stats$mean.recrate~chromosome.stats$set)
# SD of chromosome recombination rate by dataset
boxplot(chromosome.stats$sd.recrate~chromosome.stats$set)
# Coeeficient of Variation of chromosome recombination rate by dataset
boxplot(chromosome.stats$sd.recrate/chromosome.stats$mean.recrate~chromosome.stats$set)
summary(chromosome.stats$sd.recrate)
plot(density(chromosome.stats$sd.recrate/chromosome.stats$mean.recrate, na.rm = TRUE), col = "Red")


plot(chromosome.stats$mean.recrate, chromosome.stats$sd.recrate, ylim = c(0, 10), xlim = c(0, 10))
# sd increases with mean
# hence, choosing a glm model with family = Poisson

# Explore the relationship between the chromosome physical length (Mb) and the chromosome average recombination rate (cM/Mb)
hist(chromosome.stats$phys.map.length, breaks = 20)
hist(chromosome.stats$mean.recrate, breaks = 20)
plot(log(chromosome.stats$phys.map.length), chromosome.stats$mean.recrate)
plot(log(chromosome.stats$phys.map.length), chromosome.stats$sd.recrate, ylim = c(0,30))
# Heterogeneity: coefficient of variation
plot(log(chromosome.stats$phys.map.length), chromosome.stats$cv.recrate)
cor.test(log(chromosome.stats$phys.map.length), chromosome.stats$cv.recrate)
# Without the outlier species
chromosome.stats.filtered = chromosome.stats[-which(chromosome.stats$set == "Cenchrus_americanus_Pucher2017"),]
plot(log(chromosome.stats.filtered$phys.map.length), chromosome.stats.filtered$cv.recrate)
cor.test(log(chromosome.stats.filtered$phys.map.length), chromosome.stats.filtered$cv.recrate)

# Coefficient of variation outliers
chromosome.stats$set[which(chromosome.stats$cv.recrate > 2.5)]
chromosome.stats$sd.recrate[which(chromosome.stats$cv.recrate > 2.5)]
chromosome.stats$mean.recrate[which(chromosome.stats$cv.recrate > 2.5)]
# Outliers are not deviating from ~1-2 because their sd is close to 0
# Mean not closed to 0 except for Cenchrus
# Estimator seems unbiased, but effect of sd close to zero may be important for Cenchrus americanus chromosome 2


# Log ~ log shows a linear relationship
plot(log(chromosome.stats$phys.map.length), log(chromosome.stats$mean.recrate))
mod = lm(log(mean.recrate) ~ log(phys.map.length), data = chromosome.stats)
abline(mod)
summary(mod)
cor.test(log(chromosome.stats$phys.map.length), log(chromosome.stats$mean.recrate))

# Show the linear relationship between chromosome wide rate and the mean estimated recombination rate as a cross-validation procedure
plot(chromosome.stats$chrwide.rate, chromosome.stats$mean.recrate)
abline(0, 1)
# Same mean: no bias in interpolation, but may be noise produced (variance)
# Rejection of maps too far from the QQline
# i.e. points that are out of the confidence interval around the QQline (arbitrary interval of +- 0.2)
# See interval calibration below
interval = 0.5
# Plotting the interval
quantile_interval = data.frame(x = round(seq(0, 10, 0.0001), digits = 4), y.upper = seq(0, 10, 0.0001) + interval, y.lower = seq(0, 10, 0.0001) - interval)
lines(quantile_interval$x, quantile_interval$y.upper, col = "Red")
lines(quantile_interval$x, quantile_interval$y.lower, col = "Red")
# test if a map is outside the interval
df_test = data.frame(x = round(chromosome.stats$chrwide.rate , digits = 4), y = chromosome.stats$mean.recrate) # A data frame of map coordinates to test
idx_reject = logical(nrow(df_test)) # the index of maps rejected
for (i in 1:nrow(df_test)) {
  idx_reject[i] = (df_test$y[i] > quantile_interval$y.upper[which(quantile_interval$x == df_test$x[i])]) | (df_test$y[i] < quantile_interval$y.lower[quantile_interval$x == df_test$x[i]])
}
sum(idx_reject)

idx_reject_Mean = which(idx_reject)
# plot rejected maps in red
points(chromosome.stats$chrwide.rate[idx_reject], chromosome.stats$mean.recrate[idx_reject], col = "Red", bg = "Red", pch = 16)
# Interval calibration by exploring distributions
plot(density(abs(chromosome.stats$mean.recrate-chromosome.stats$chrwide.rate)), col = "Red")
plot(density(abs(chromosome.stats$mean.recrate-chromosome.stats$chrwide.rate)), xlim = c(0, 2), col = "Red")
abline(v = 1, col = "Red", lty = 2)
# We see a peak around 0.1-0.2



# List of rejected maps
rejected_maps = unique(c(idx_reject_Mean, idx_reject_SD))
# Save in tables
chromosome.stats = read.table("tables/chromosome.stats.csv", header = TRUE, sep = ";")
chromosome.stats$valid[rejected_maps] = FALSE
write.table(chromosome.stats, file = "tables/chromosome.stats.csv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ";")


# Explore the relationship between the haploid chromosome number and the chromosome average recombination rate (cM/Mb)
# A vector of HCN retrieved in 'metadata'
HCN = as.character(chromosome.stats$set)
for (i in 1:length(HCN)) {
  HCN[i]= metadata$hcn[which(as.character(metadata$id) == as.character(HCN[i]))]
}
HCN
boxplot((chromosome.stats$mean.recrate)~HCN)

# Explore the relationship between the ploidy level and the chromosome average recombination rate (cM/Mb)
# A vector of ploidy  retrieved in 'metadata'
ploidy = as.character(chromosome.stats$set)
for (i in 1:length(ploidy)) {
  ploidy[i]= metadata$ploidy[which(as.character(metadata$id) == as.character(ploidy[i]))]
}
ploidy
boxplot(chromosome.stats$mean.recrate~ploidy)

# Explore the relationship between the total number of chromosomes and the chromosome average recombination rate (cM/Mb)
plot(as.numeric(ploidy)*as.numeric(HCN), chromosome.stats$mean.recrate)

#============================================================================#
# Centromeric index ----
#============================================================================#

# Maps without metadata about the centromeric index
unique(chromosome.stats$set[is.na(chromosome.stats$centromeric_index)])
# Maps with a centromeric index
unique(chromosome.stats$set[!is.na(chromosome.stats$centromeric_index)])
length(unique(chromosome.stats$set[!is.na(chromosome.stats$centromeric_index)]))

#============================================================================#
# Control if there is an effect of map quality on estimates ----
#============================================================================#
#----------------------------------------------------------------------------#
# Number of markers ----
#----------------------------------------------------------------------------#
qqnorm(log(chromosome.stats$mean.recrate))
qqnorm(chromosome.stats$nb.markers)
# Some outliers, maps with a very large number of markers (i.e. Panicum hallii)
qqnorm(log(chromosome.stats$mean.recrate[which(chromosome.stats$nb.markers < 2000)]))
qqnorm(chromosome.stats$nb.markers[which(chromosome.stats$nb.markers < 2000)])
# Test the correlation
cor.test(chromosome.stats$nb.markers[which(chromosome.stats$nb.markers < 2000)],
         log(chromosome.stats$mean.recrate[which(chromosome.stats$nb.markers < 2000)]), method = "spearman")
# Plot the relationship
plot(chromosome.stats$nb.markers[which(chromosome.stats$nb.markers < 2000)],
     log(chromosome.stats$mean.recrate[which(chromosome.stats$nb.markers < 2000)]),
     xlab = "Number of markers", ylab = "Mean recombination rate")

#----------------------------------------------------------------------------#
# Marker density (cM) ----
#----------------------------------------------------------------------------#
qqnorm(chromosome.stats$density.markers.cM)
# Some outliers, maps with a very large number of markers (i.e. Panicum hallii)
qqnorm(chromosome.stats$density.markers.cM[which(chromosome.stats$nb.markers < 2000)])
# Test the correlation
cor.test(chromosome.stats$density.markers.cM[which(chromosome.stats$nb.markers < 2000)],
         log(chromosome.stats$mean.recrate[which(chromosome.stats$nb.markers < 2000)]), method = "spearman")
# Plot the relationship
plot(chromosome.stats$density.markers.cM[which(chromosome.stats$nb.markers < 2000)],
     log(chromosome.stats$mean.recrate[which(chromosome.stats$nb.markers < 2000)]),
     xlab = "Marker density (cM)", ylab = "Mean recombination rate")

#----------------------------------------------------------------------------#
# Marker density (bp) ----
#----------------------------------------------------------------------------#
qqnorm(chromosome.stats$density.markers.bp)
# Some outliers, maps with a very large number of markers (i.e. Panicum hallii)
qqnorm(chromosome.stats$density.markers.bp[which(chromosome.stats$nb.markers < 2000)])
# Test the correlation
cor.test(chromosome.stats$density.markers.bp[which(chromosome.stats$nb.markers < 2000)],
         log(chromosome.stats$mean.recrate[which(chromosome.stats$nb.markers < 2000)]), method = "spearman")
# Plot the relationship
plot(chromosome.stats$density.markers.bp[which(chromosome.stats$nb.markers < 2000)],
     log(chromosome.stats$mean.recrate[which(chromosome.stats$nb.markers < 2000)]),
     xlab = "Marker density (bp)", ylab = "Mean recombination rate")

#----------------------------------------------------------------------------#
# Number of progeny ----
#----------------------------------------------------------------------------#
# Data per species
# Aggregate
species.stats = aggregate(chromosome.stats, by = list(chromosome.stats$species), mean)
colnames(species.stats)[1] = "species"
species.stats$nb_progeny = NA
for (i in 1:nrow(species.stats)) {
  species.stats$nb_progeny[i] = metadata.clean$n_progeny[which(gsub("_", " ", as.character(metadata.clean$species)) == species.stats$species[i])][1]
}
# Data normality check
qqnorm(species.stats$nb_progeny)
cor.test(species.stats$nb_progeny, log(species.stats$mean.recrate), method = "spearman")
# Plot the relationship
plot(species.stats$nb_progeny, log(species.stats$mean.recrate),
     xlab = "Progeny number", ylab = "Mean recombination rate")

# A small tendency of larger progeny number to decrease the recombination rate
# Lower number of descendants inflates genetic distances?

#============================================================================#
# Assess correlations between summary statistics and genome characteristics ----
# with Spearman correlation and Linear model
#============================================================================#

#----------------------------------------------------------------------------#
# Mean recombination rate ~ Chromosome size ----
#----------------------------------------------------------------------------#
# Log ~ log shows a linear relationship
plot(log(chromosome.stats$phys.map.length/1000000), log(chromosome.stats$mean.recrate))
mod = lm(log(mean.recrate) ~ log(phys.map.length/1000000), data = chromosome.stats)
abline(mod)
summary(mod)
AIC(mod)
cor.test(log(chromosome.stats$phys.map.length/1000000), log(chromosome.stats$mean.recrate), method = "spearman")

#----------------------------------------------------------------------------#
# Recombination rate coefficient of variation ~ Chromosome size ----
#----------------------------------------------------------------------------#
plot(log(chromosome.stats$phys.map.length/1000000), chromosome.stats$cv.recrate)
mod = lm(cv.recrate ~ log(phys.map.length/1000000), data = chromosome.stats)
abline(mod)
summary(mod)
AIC(mod)
cor.test(log(chromosome.stats$phys.map.length/1000000), chromosome.stats$cv.recrate, method = "spearman")

hist(chromosome.stats$cv.recrate)
hist(log(chromosome.stats$cv.recrate))

plot(log(chromosome.stats$phys.map.length/1000000), log(chromosome.stats$cv.recrate))
mod = lm(log(cv.recrate) ~ log(phys.map.length/1000000), data = chromosome.stats)
abline(mod)
summary(mod)
AIC(mod)
cor.test(log(chromosome.stats$phys.map.length/1000000), log(chromosome.stats$cv.recrate), method = "spearman")

#----------------------------------------------------------------------------#
# Gini coefficient ~ Chromosome size ----
#----------------------------------------------------------------------------#
plot(log(chromosome.stats$phys.map.length/1000000), chromosome.stats$gini)
mod = lm(gini ~ log(phys.map.length/1000000), data = chromosome.stats)
abline(mod)
summary(mod)
AIC(mod)
cor.test(log(chromosome.stats$phys.map.length/1000000), chromosome.stats$gini, method = "spearman")

hist(chromosome.stats$gini)
hist(log(chromosome.stats$gini))

plot(log(chromosome.stats$phys.map.length/1000000), log(chromosome.stats$gini))
mod = lm(log(gini) ~ log(phys.map.length/1000000), data = chromosome.stats)
abline(mod)
summary(mod)
AIC(mod)
cor.test(log(chromosome.stats$phys.map.length/1000000), log(chromosome.stats$gini), method = "spearman")

#----------------------------------------------------------------------------#
# Bias toward periphery - Skewness ~ Chromosome size ----
#----------------------------------------------------------------------------#
plot(log(chromosome.stats$phys.map.length/1000000), sqrt(chromosome.stats$peripherybias_ratio))
mod = lm(sqrt(peripherybias_ratio) ~ log(phys.map.length/1000000), data = chromosome.stats)
abline(mod)
summary(mod)
AIC(mod)
cor.test(log(chromosome.stats$phys.map.length/1000000), sqrt(chromosome.stats$peripherybias_ratio), method = "spearman")

#----------------------------------------------------------------------------#
# Mean recombination rate ~ Haploid Chromosome Number ----
#----------------------------------------------------------------------------#
# Data per species
# Aggregate
species.stats = aggregate(chromosome.stats, by = list(chromosome.stats$species), mean)
colnames(species.stats)[1] = "species"
species.stats$hcn = NA
for (i in 1:nrow(species.stats)) {
  species.stats$hcn[i] = metadata.clean$hcn[which(gsub("_", " ", as.character(metadata.clean$species)) == species.stats$species[i])][1]
}

plot(species.stats$hcn, log(species.stats$mean.recrate))
mod = lm(log(mean.recrate) ~ species.stats$hcn, data = species.stats)
abline(mod)
summary(mod)
cor.test(species.stats$hcn, log(species.stats$mean.recrate), method = "spearman")

#----------------------------------------------------------------------------#
# Cv recombination rate ~ Haploid Chromosome Number ----
#----------------------------------------------------------------------------#
# Data per species
plot(species.stats$hcn, species.stats$cv.recrate)
mod = lm(cv.recrate ~ species.stats$hcn, data = species.stats)
abline(mod)
summary(mod)
cor.test(species.stats$hcn, species.stats$cv.recrate, method = "spearman")

#----------------------------------------------------------------------------#
# Gini ~ Haploid Chromosome Number ----
#----------------------------------------------------------------------------#
# Data per species
plot(species.stats$hcn, species.stats$gini)
mod = lm(gini ~ species.stats$hcn, data = species.stats)
abline(mod)
summary(mod)
cor.test(species.stats$hcn, species.stats$gini, method = "spearman")

#----------------------------------------------------------------------------#
# Periphery-bias ratio ~ Haploid Chromosome Number ----
#----------------------------------------------------------------------------#
# Data per species
plot(species.stats$hcn, log(species.stats$peripherybias_ratio))
mod = lm(log(peripherybias_ratio) ~ species.stats$hcn, data = species.stats)
abline(mod)
summary(mod)
cor.test(species.stats$hcn, log(species.stats$peripherybias_ratio), method = "spearman")



#============================================================================#
# LINEAR MIXED MODEL ----
# Assess correlations between summary statistics and genome characteristics ----
# while controlling for phylogeny by fitting a Mixed Linear Model
#============================================================================#

lmm.data = cbind(data.frame(Species = gsub("_", " ", unlist(lapply(chromosome.stats$set, function(x){metadata.clean$species[which(as.character(metadata.clean$id) == as.character(x))]})))),
      chromosome.stats)
#----------------------------------------------------------------------------#
# Mean recombination rate ~ Chromosome size ----
#----------------------------------------------------------------------------#
lmer.model = lmer(log(mean.recrate)~log(phys.map.length/1000000) + (1|Species), data = lmm.data)
lmer.model = lmer(log(mean.recrate)~log(phys.map.length/1000000) + (log(phys.map.length/1000000)|Species), data = lmm.data)
summary(lmer.model)
anova(lmer.model)
plot(lmer.model)
AIC(lmer.model)
# Diagnostic plots
qqPlot(residuals(lmer.model))
# Variance explained by random effects
randomeff = as.data.frame(VarCorr(lmer.model,comp="Variance"))
randomeff$vcov[1]/sum(randomeff$vcov)
# Plot the fitted regression line over data
plot(log(chromosome.stats$phys.map.length/1000000), log(chromosome.stats$mean.recrate))
abline(a = fixef(lmer.model)[1], b = fixef(lmer.model)[2])

#----------------------------------------------------------------------------#
# Recombination rate coefficient of variation ~ Chromosome size ----
#----------------------------------------------------------------------------#
lmer.model = lmer(cv.recrate~log(phys.map.length/1000000) + (1|Species), data = lmm.data)
# lmer.model = lmer(log(cv.recrate)~log(phys.map.length/1000000) + (1|Species), data = lmm.data)
summary(lmer.model)
plot(lmer.model)
AIC(lmer.model)
# Variance explained by random effects
randomeff = as.data.frame(VarCorr(lmer.model,comp="Variance"))
randomeff$vcov[1]/sum(randomeff$vcov)
# Plot the fitted regression line over data
plot(log(chromosome.stats$phys.map.length/1000000), chromosome.stats$cv.recrate)
abline(a = fixef(lmer.model)[1], b = fixef(lmer.model)[2])

#----------------------------------------------------------------------------#
# Gini coefficient ~ Chromosome size ----
#----------------------------------------------------------------------------#
lmer.model = lmer(gini~log(phys.map.length/1000000) + (1|Species), data = lmm.data)
summary(lmer.model)
(lmer.model)
plot(lmer.model)
AIC(lmer.model)
# Variance explained by random effects
randomeff = as.data.frame(VarCorr(lmer.model,comp="Variance"))
randomeff$vcov
randomeff$vcov[1]/sum(randomeff$vcov)
# Plot the fitted regression line over data
plot(log(chromosome.stats$phys.map.length/1000000), chromosome.stats$gini)
abline(a = fixef(lmer.model)[1], b = fixef(lmer.model)[2])

# Yet, even if lmer take into account the autocorrelated effect of multiple chromosomes of the same species,
# it does not deal with phylogenetic autocorrelation
# To do this we would need to adjust the model with a covariance matrix of phylogenetic distances

#----------------------------------------------------------------------------#
# Bias toward periphery - Skewness ~ Chromosome size ----
#----------------------------------------------------------------------------#
lmer.model = lmer(sqrt(peripherybias_ratio)~log(phys.map.length/1000000) + (1|Species), data = lmm.data)
summary(lmer.model)
(lmer.model)
plot(lmer.model)
AIC(lmer.model)
# Variance explained by random effects
randomeff = as.data.frame(VarCorr(lmer.model,comp="Variance"))
randomeff$vcov
randomeff$vcov[1]/sum(randomeff$vcov)
# Plot the fitted regression line over data
plot(log(chromosome.stats$phys.map.length/1000000), chromosome.stats$peripherybias_ratio)
abline(a = fixef(lmer.model)[1], b = fixef(lmer.model)[2])


#============================================================================#
# PGLMM ----
# Assess correlations between summary statistics and genome characteristics  ----
# while controlling for phylogeny by fitting a Phylogenetic Generalized Linear Model
#============================================================================#
# In a basic linear regression, the response variable (y) is modelled as the product of the estimated
# coefficients (Beta) and the explanatory variables (X), plus residual variation (error): y = Beta.X + errors. If the
# cases in the model are related taxa, then the values in y and X are no longer independent: each
# value will be more similar to some values (closer relatives) than others (distant relatives).
# The pgls function addresses this problem by incorporating the covariance between taxa into
# the calculation of estimated coefficients: this is a generalized least squares (GLS) model. The
# covariance matrix (V), showing the expected covariance between each pair of tips is calculated
# using the branch lengths of a phylogeny showing how the taxa are related.

# chromosome.stats = recombination.map.statistics(data.final, save.file = TRUE)
# metadata.clean = read.csv(file = paste(wd, "/data-cleaned/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")

#----------------------------------------------------------------------------#
# Build the phylogeny ----
#----------------------------------------------------------------------------#
# library(rotl)
# Load the phylogenetic tree estimated in 'Phylogeny.R'
load("output/phylogeny/tree.Rdata")
tree$tip.label

# Initial covariance matrix
# The initial covariance matrix (V) taken from an ultrametric phylogeny assumes a Brownian
# model of evolution: that variation between tips accumulates along all branches of the tree at a
# rate proportional to the length of the branches.
library(phyr)
# vcv_tree = vcv2(tree)
# vcv_tree = vcv2(tree, corr = TRUE) # return a correlation matrix

#----------------------------------------------------------------------------#
# Build the traits dataset ----
#----------------------------------------------------------------------------#
# gsub("_", " ", unlist(lapply(chromosome.stats$set, function(x){metadata.clean$species[which(as.character(metadata.clean$id) == as.character(x))]})))
# phylogenetic_traits = cbind(data.frame(Species = gsub("_", " ", unlist(lapply(chromosome.stats$set, function(x){metadata.clean$species[which(as.character(metadata.clean$id) == as.character(x))]})))))
phylogenetic_traits = cbind(data.frame(Species = gsub("_", " ", unlist(lapply(chromosome.stats$set, function(x){metadata.clean$species[which(as.character(metadata.clean$id) == as.character(x))]})))),
                            chromosome.stats)

#----------------------------------------------------------------------------#
# Test for the phylogenetic signal ----
#----------------------------------------------------------------------------#
# When applied without predictor (independent) variables, PGLMM gives a test for phylogenetic signal.
pglmm.model = pglmm(log(mean.recrate) ~ 1 + (1|Species__), data = phylogenetic_traits, family = 'gaussian',
                    cov_ranef = list(Species = tree), bayes = TRUE)
summary(pglmm.model)
# Species random effect significantly different of 0 (variance explained C.I. = 0.09505  0.64216)
# And phylogenetic signal too was significant (variance explained C.I. = 0.11707  0.72653)
# Residuals had a very low variance (0.03261)

# Look for correlations between variables under phylogenetic autocorrelation
cor_phylo(~mean.recrate+cv.recrate+gini, species=unique(phylogenetic_traits$Species), phy = tree, data = phylogenetic_traits)


#----------------------------------------------------------------------------#
# Mean recombination rate ~ Chromosome size ----
#----------------------------------------------------------------------------#
pglmm.model = pglmm(log(mean.recrate)~log(phys.map.length/1000000) + (1|Species), data = phylogenetic_traits, family = 'gaussian',
                    cov_ranef = list(Species = tree))
pglmm.model = pglmm(log(mean.recrate)~log(phys.map.length/1000000) + (1|Species__), data = phylogenetic_traits, family = 'gaussian',
                    cov_ranef = list(Species = tree))
# pglmm.model = pglmm(log(mean.recrate)~log(phys.map.length/1000000) + (0|Species__) + (log(phys.map.length/1000000)|Species__),
#                     data = phylogenetic_traits, family = 'gaussian',
#                     cov_ranef = list(Species = tree))
# Model fail to converge, not enough individuals in each species for random slopes

summary(pglmm.model)
# Diagnostic plots
# Residuals distribution
qqPlot(residuals(pglmm.model))
# Residuals as a function of fitted
plot(fitted(pglmm.model), residuals(pglmm.model))
abline(h = 0)
# Observed~predicted qqplot
plot(unlist(pglmm_predicted_values(pglmm.model)), pglmm.model$Y)
abline(a = 0, b = 1)

# Control the influence of map quality on the statistical relationship
# by fitting number of markers as a function of residuals
control.mod = lm(phylogenetic_traits$nb.markers ~ residuals(pglmm.model))
summary(control.mod)
cor.test(phylogenetic_traits$nb.markers, residuals(pglmm.model))
cor.test(phylogenetic_traits$density.markers.bp, residuals(pglmm.model))
cor.test(phylogenetic_traits$density.markers.bp, residuals(pglmm.model))
# NS coefficient; QUlaity of the data did not influenced the relationship

mod = lm(log(mean.recrate) ~ log(phys.map.length/1000000), data = chromosome.stats)
# Simulate the expected regression line under the assumption of one CO per chromosome (i.e. a genetic length of 50cM)
expectedline = 100/(chromosome.stats$phys.map.length/1000000)
# Plot the fitted regression line over data
plot(log(chromosome.stats$phys.map.length/1000000), log(chromosome.stats$mean.recrate))
abline(a = fixef(pglmm.model)[[1]][1], b = fixef(pglmm.model)[[1]][2])
abline(mod, lty = 2)
lines(log(chromosome.stats$phys.map.length/1000000), log(expectedline), col = "red")

#----------------------------------------------------------------------------#
# Investigating the strange species effect that reduces the slope of the mixed model

# Look for species sorted by ascending mean recombination rate
phylogenetic_traits$Species[order(phylogenetic_traits$mean.recrate)]

# Problem was that some species didn't had enough chromosomes to check for intra-specific random slope
# Hence, subset only species with at least 5 chromosomes to check this effect
(list_species = as.character(unique(phylogenetic_traits$Species)[table(phylogenetic_traits$Species) > 4]))
phylogenetic_traits_subset = subset(phylogenetic_traits, Species %in% list_species)
phylogenetic_traits_subset$Species = as.character(phylogenetic_traits_subset$Species)
phylogenetic_traits_subset$chrsize = log(phylogenetic_traits_subset$phys.map.length/1000000)
phylogenetic_traits_subset$meanreclog = log(phylogenetic_traits_subset$mean.recrate)
# Drop tips in the phylogeny
tree_subset = drop.tip(tree, as.character(unique(phylogenetic_traits$Species)[table(phylogenetic_traits$Species) < 5]))
pglmm.model = pglmm(meanreclog~chrsize + (1|Species__) + (chrsize|Species__),
                    data = phylogenetic_traits_subset, family = 'gaussian',
                    cov_ranef = list(Species = tree_subset))
summary(pglmm.model)
fixef(pglmm.model)
ranef(pglmm.model)
# Diagnostic plots
# Residuals distribution
qqPlot(residuals(pglmm.model))
# Residuals as a function of fitted
plot(fitted(pglmm.model), residuals(pglmm.model))
abline(h = 0)
# Observed~predicted qqplot
plot(unlist(pglmm_predicted_values(pglmm.model)), pglmm.model$Y)
abline(a = 0, b = 1)

# PGLMM with phyr does not allow the estimation of random slopes coefficients
# Hence switch to lmer model that is very close to the PGLMM model
library(lme4)
library(sjPlot)
library(sjmisc)
mixed.model = lmer(log(mean.recrate)~log(phys.map.length/1000000) + (log(phys.map.length/1000000)|Species),
                    data = phylogenetic_traits_subset)
summary(mixed.model)
AIC(mixed.model)
plot(mixed.model)
# Diagnostic plots
# Residuals distribution
qqPlot(residuals(mixed.model))
# Residuals as a function of fitted
plot(fitted(mixed.model), residuals(mixed.model))
abline(h = 0)
# Observed~predicted qqplot
plot(predict(mixed.model), log(phylogenetic_traits_subset$mean.recrate))
abline(a = 0, b = 1)

# Plotting random effects
(ranint = coef(mixed.model)$Species[,1]) # Random intercepts
(ranslope = coef(mixed.model)$Species[,2]) # Random slopes

# Plotting fixed (blue) & random (green) effects
plot_model(mixed.model)
tab_model(mixed.model)
plot(log(chromosome.stats$phys.map.length/1000000), log(chromosome.stats$mean.recrate))
abline(a = fixef(pglmm.model)[[1]][1], b = fixef(pglmm.model)[[1]][2])
mod = lm(log(mean.recrate) ~ log(phys.map.length/1000000), data = chromosome.stats)
abline(mod, lty = 2)
abline(fixef(mixed.model)[1], fixef(mixed.model)[2], col = "Blue")
# Boundaries of random segments (mean of x values ± e^0.5)
boundaries = data.frame(species = rownames(coef(mixed.model)$Species), upper = NA, lower = NA)
boundaries
chromosome.stats$species = as.character(chromosome.stats$species)
boundaries$species = as.character(boundaries$species)
for (i in 1:nrow(boundaries)) {
  boundaries$upper[i] = log(mean(chromosome.stats$phys.map.length[which(chromosome.stats$species == boundaries$species[i])], na.rm = TRUE)/1000000) + 0.5
  boundaries$lower[i] = log(mean(chromosome.stats$phys.map.length[which(chromosome.stats$species == boundaries$species[i])], na.rm = TRUE)/1000000) - 0.5
}
boundaries
for (i in 1:length(ranint)) {
  clip(boundaries$lower[i], boundaries$upper[i], -5, 5)
  abline(a = ranint[i], b = ranslope[i], col = "red")
}
# Correlations with the random slope?
randomtable = coef(mixed.model)$Species
randomtable$species = row.names(randomtable)
randomtable$chrsize = NA
randomtable$meanrecrate = NA
for (i in 1:nrow(randomtable)) {
  randomtable$chrsize[i] = mean(chromosome.stats$phys.map.length[which(chromosome.stats$species == randomtable$species[i])], na.rm = TRUE)
  randomtable$meanrecrate[i] = mean(chromosome.stats$mean.recrate[which(chromosome.stats$species == randomtable$species[i])], na.rm = TRUE)
}
colnames(randomtable) = c("randomintercept", "randomslope", "species", "chrsize", "meanrecrate")
plot(log(randomtable$chrsize/1000000), randomtable$randomslope, xlab = "Chromosome size", ylab = "Species random slope")
cor.test(log(randomtable$chrsize/1000000), randomtable$randomslope, method = "spearman")

plot(log(randomtable$meanrecrate), randomtable$randomslope, xlab = "Species mean recombination rate", ylab = "Species random slope")
cor.test(log(randomtable$meanrecrate), randomtable$randomslope, method = "spearman")
# Strong correlation between the random slopes and the chromosome size/mean recombination rate
# Larger chromosomes have higher random slopes (i.e. stronger effect of chromosome size)
# Could it be an effect of the log scale?

# The specific correlation between mean recombination rate and chromosome size
species_correlation  = function(sp) {
  p = cor.test(chromosomes$mean.recrate[which(chromosomes$species == sp)], chromosomes$phys.map.length[which(chromosomes$species == sp)], method = "spearman")$p.value
  stat = cor.test(chromosomes$mean.recrate[which(chromosomes$species == sp)], chromosomes$phys.map.length[which(chromosomes$species == sp)], method = "spearman")$estimate
  return(data.frame(species = sp, mean = mean(chromosomes$mean.recrate[which(chromosomes$species == sp)], na.rm = TRUE),
                    chrsize = mean(chromosomes$phys.map.length[which(chromosomes$species == sp)], na.rm = TRUE),
                    correlation = stat, pvalue = p))
}
(cor_mean_chrsize = bind_rows(lapply(unique(phylogenetic_traits_subset$Species), function(x) species_correlation(x))))
# Significant correlations at 5%
(cor_mean_chrsize$significant = cor_mean_chrsize$pvalue < 0.05)
cor_mean_chrsize
plot(log(cor_mean_chrsize$chrsize/1000000), abs(cor_mean_chrsize$correlation))
boxplot(cor_mean_chrsize$mean~cor_mean_chrsize$significant)
boxplot(cor_mean_chrsize$chrsize~cor_mean_chrsize$significant)


# SIMULATION APPROACH

plot(phylogenetic_traits$phys.map.length, phylogenetic_traits$mean.recrate)
# Simulate random slopes of the original correlation before transformation
# Ten species with ten points each
n = 10 # number of chromosomes sampled
N = 10 # Number of species sampled
df_simu = data.frame(sp = rep(1:N, each = n))
df_simu$x = exp(jitter(df_simu$sp)) # Increasing exp
hist(df_simu$x)
df_simu$y = exp(jitter((n+1)-df_simu$sp)) # Decreasing exp
hist(df_simu$y)
plot(df_simu$x, df_simu$y)
plot(log(df_simu$x), log(df_simu$y))
mod = lm(log(y)~log(x), data = df_simu)
abline(mod, col = "blue")
mixed.model = lmer(log(y)~log(x) + (x|sp), data = df_simu)
abline(fixef(mixed.model)[1], fixef(mixed.model)[2], col = "red")

# Go further by simulating mean rec rate for a given physical length
# Divide 50cM by the physical length
# Same method as in Haenel et al. 2018 for a theoretical distribution
df_simu2 = data.frame(sp = phylogenetic_traits$Species, x = phylogenetic_traits$phys.map.length/1000000)
df_simu2$y = 50/df_simu2$x
plot(df_simu2$x, df_simu2$y)
hist(df_simu2$x)
hist(df_simu2$y)
plot(log(df_simu2$x), log(df_simu2$y))
mod = lm(log(y)~log(x), data = df_simu2)
abline(mod, col = "blue")
mixed.model = lmer(log(y)~log(x) + (log(x)|sp), data = df_simu2)
mixed.model
abline(fixef(mixed.model)[1], fixef(mixed.model)[2], col = "red")
pglmm.model = pglmm(log(y)~log(x) + (1| sp__), data = df_simu2, family = 'gaussian',
                    cov_ranef = list(sp = tree))
summary(pglmm.model)
abline(a = fixef(pglmm.model)[[1]][1], b = fixef(pglmm.model)[[1]][2], col = "green")

# Now test what happens if effect on recombination is stronger for larger chromosomes than smaller ones
# Hence effect is proportional to chromosome size
n = 10 # number of chromosomes sampled
N = 10 # Number of species sampled
df_simu3 = data.frame(sp = rep(1:N, each = n))
df_simu3$x = exp(jitter(df_simu3$sp)) # Increasing exp
hist(df_simu3$x)
# df_simu3$y = exp((n+1) - jitter(df_simu3$sp, df_simu3$sp)) # Decreasing exp
df_simu3$y = 1/exp(df_simu3$sp*(jitter(df_simu3$sp, factor = df_simu3$sp)/N)) # Decreasing exp
hist(df_simu3$y)
plot(log(df_simu3$x), log(df_simu3$y))
mod = lm(log(y)~log(x), data = df_simu3)
abline(mod, col = "blue")
mixed.model = lmer(log(y)~log(x) + (log(x)|sp), data = df_simu3)
mixed.model
abline(fixef(mixed.model)[1], fixef(mixed.model)[2], col = "red")
# Not very satisfying...
# A new approach, respecting the actual sampling
df_simu4 = data.frame(sp = phylogenetic_traits$Species, x = phylogenetic_traits$phys.map.length/1000000)
# Simulate the genetic map length with variance scaled by the physical genetic length
# Larger variance in genetic size when physical size increases
df_simu4$y = rnorm(nrow(df_simu4), mean = 100, sd = 20*(1+df_simu4$x/max(df_simu4$x)))/df_simu4$x
plot(df_simu4$x, df_simu4$y)
hist(df_simu4$x)
hist(df_simu4$y)
plot(log(df_simu4$x), log(df_simu4$y))
mod = lm(log(y)~log(x), data = df_simu4)
abline(mod, col = "blue")
mixed.model = lmer(log(y)~log(x) + (log(x)|sp), data = df_simu4)
mixed.model
abline(fixef(mixed.model)[1], fixef(mixed.model)[2], col = "red")
# Even testing phylogenetic signal with the actual sampling
pglmm.model = pglmm(log(y)~log(x) + (1| sp__), data = df_simu4, family = 'gaussian',
                    cov_ranef = list(sp = tree))
summary(pglmm.model)
# No phylogenetic signal observed when simulating data
# Hence actual phylogenetic signal in observed data must be either
# biological
# or sampling bias
abline(a = fixef(pglmm.model)[[1]][1], b = fixef(pglmm.model)[[1]][2], col = "green")





#----------------------------------------------------------------------------#
# Recombination rate coefficient of variation ~ Chromosome size ----
#----------------------------------------------------------------------------#
pglmm.model = pglmm(log(cv.recrate)~log(phys.map.length/1000000) + (1| Species), data = phylogenetic_traits, family = 'gaussian',
                    cov_ranef = list(Species = tree))
summary(pglmm.model)
pglmm.model = pglmm(log(cv.recrate)~log(phys.map.length/1000000) + (1| Species__), data = phylogenetic_traits, family = 'gaussian',
                    cov_ranef = list(Species = tree))
summary(pglmm.model)

plot(log(phylogenetic_traits$cv.recrate)~log(phylogenetic_traits$phys.map.length/1000000))
abline(a = fixef(pglmm.model)[[1]][1], b = fixef(pglmm.model)[[1]][2])

# Residuals distribution
qqPlot(residuals(pglmm.model))
# Residuals as a function of fitted
plot(fitted(pglmm.model), residuals(pglmm.model))
abline(h = 0)
# Observed~predicted qqplot
plot(unlist(pglmm_predicted_values(pglmm.model)), pglmm.model$Y)
abline(a = 0, b = 1)

# Control the influence of map quality on the statistical relationship
# by fitting number of markers as a function of residuals
control.mod = lm(phylogenetic_traits$nb.markers ~ residuals(pglmm.model))
summary(control.mod)
cor.test(phylogenetic_traits$nb.markers, residuals(pglmm.model))
cor.test(phylogenetic_traits$density.markers.bp, residuals(pglmm.model))
cor.test(phylogenetic_traits$density.markers.bp, residuals(pglmm.model))
# NS coefficient; QUlaity of the data did not influenced the relationship

# fixef(pglmm.model) # Fixed effects coefficients
# Plot the fitted regression line over data
plot(log(chromosome.stats$phys.map.length/1000000), chromosome.stats$cv.recrate)
abline(a = fixef(pglmm.model)[[1]][1], b = fixef(pglmm.model)[[1]][2])

mod = lm(cv.recrate ~ log(phys.map.length/1000000), data = chromosome.stats)
# Plot the fitted regression line over data
plot(log(chromosome.stats$phys.map.length/1000000), chromosome.stats$cv.recrate)
abline(a = fixef(pglmm.model)[[1]][1], b = fixef(pglmm.model)[[1]][2])
abline(mod, lty = 2)
# The expected relationship is sqrt(var)/mean, where var = mean
# expectedline = 100/(chromosome.stats$phys.map.length/1000000) # The expected mean recombination rate
# expectedline = sqrt(expectedline)/expectedline # Hence the expected coefficient of variation under the assumption of a Poisson process
# lines(log(chromosome.stats$phys.map.length/1000000), log(expectedline), col = "red")

# Simulate form data
# expectedline = (sqrt(chromosome.stats$mean.recrate)/chromosome.stats$mean.recrate)
# modexpected = lm(expectedline~log(chromosome.stats$phys.map.length/1000000))
# abline(modexpected, col = "red")
# pglmm.model.expected = pglmm((sqrt(mean.recrate)/mean.recrate)~log(phys.map.length/1000000) + (1| Species__), data = phylogenetic_traits, family = 'gaussian',
#                     cov_ranef = list(Species = tree))
# abline(a = fixef(pglmm.model.expected)[[1]][1], b = fixef(pglmm.model.expected)[[1]][2], col = "Red")


#----------------------------------------------------------------------------#
# Gini coefficient ~ Chromosome size ----
#----------------------------------------------------------------------------#
pglmm.model = pglmm(gini~log(phys.map.length/1000000) + (1|Species), data = phylogenetic_traits, family = 'gaussian',
                    cov_ranef = list(Species = tree))
summary(pglmm.model)
pglmm.model = pglmm(gini~log(phys.map.length/1000000) + (1|Species__), data = phylogenetic_traits, family = 'gaussian',
                    cov_ranef = list(Species = tree))
summary(pglmm.model)
# Residuals distribution
qqPlot(residuals(pglmm.model))
# Residuals as a function of fitted
plot(fitted(pglmm.model), residuals(pglmm.model))
abline(h = 0)
# Observed~predicted qqplot
plot(unlist(pglmm_predicted_values(pglmm.model)), pglmm.model$Y)
abline(a = 0, b = 1)

# Control the influence of map quality on the statistical relationship
# by fitting number of markers as a function of residuals
control.mod = lm(phylogenetic_traits$nb.markers ~ residuals(pglmm.model))
summary(control.mod)
cor.test(phylogenetic_traits$nb.markers, residuals(pglmm.model))
cor.test(phylogenetic_traits$density.markers.bp, residuals(pglmm.model))
cor.test(phylogenetic_traits$density.markers.bp, residuals(pglmm.model))
# NS coefficient; QUlaity of the data did not influenced the relationship

# fixef(pglmm.model) # Fixed effects coefficients
# Plot the fitted regression line over data
plot(log(chromosome.stats$phys.map.length/1000000), chromosome.stats$gini)
abline(a = fixef(pglmm.model)[[1]][1], b = fixef(pglmm.model)[[1]][2])

mod = lm(gini ~ log(phys.map.length/1000000), data = chromosome.stats)
# Plot the fitted regression line over data
plot(log(chromosome.stats$phys.map.length/1000000), chromosome.stats$gini)
abline(a = fixef(pglmm.model)[[1]][1], b = fixef(pglmm.model)[[1]][2])
abline(mod, lty = 2)
# The expected regression of the Gini coefficient
abline(h = 0, col = "red")

# Finally, there seems to be an effect of chromosome size only for  the mean recombination rate
# the variance of recombination rates seems more species-specific than chromosome size-specific

# Hence, is there correlation between the CV or Gini coefficient and the phylogeny?
# Or is it just a sampling bias? Poaceae are also the biggest genomes in the dataset,
# hence there seems to be a phylogenetic bias with the association Poaceae/genome size
# Represent both relationships: Spearman correlation (i.e. linear relationship) AND pglmm on the same plot

# Sensitivity analyses with sensiPhy
# https://cran.r-project.org/web/packages/sensiPhy/vignettes/sensiPhy_vignette.html

#----------------------------------------------------------------------------#
# Bias toward periphery - Skewness ~ Chromosome size ----
#----------------------------------------------------------------------------#
subset = phylogenetic_traits[!is.na(phylogenetic_traits$peripherybias_ratio),]
subset$peripherybias_ratio = subset$peripherybias_ratio + 0.0001
# subset = subset[-c(556, 309, 139, 566),]
pglmm.model = pglmm(sqrt(peripherybias_ratio)~log(phys.map.length/1000000) + (1|Species), data = subset, family = 'gaussian',
                    cov_ranef = list(Species = tree))
# pglmm.model = pglmm(sqrt(peripherybias_ratio)~log(phys.map.length/1000000) + (1|Species__), data = phylogenetic_traits, family = 'gaussian',
#                     cov_ranef = list(Species = tree))
summary(pglmm.model)
# Residuals distribution
qqPlot(residuals(pglmm.model))
# Residuals as a function of fitted
plot(fitted(pglmm.model), residuals(pglmm.model))
abline(h = 0)
# Observed~predicted qqplot
plot(unlist(pglmm_predicted_values(pglmm.model)), pglmm.model$Y)
abline(a = 0, b = 1)

# Control the influence of map quality on the statistical relationship
# by fitting number of markers as a function of residuals
control.mod = lm(subset$nb.markers ~ residuals(pglmm.model))
summary(control.mod)
cor.test(subset$nb.markers, residuals(pglmm.model))
cor.test(subset$density.markers.bp, residuals(pglmm.model))
cor.test(subset$density.markers.bp, residuals(pglmm.model))
# NS coefficient; QUlaity of the data did not influenced the relationship

# fixef(pglmm.model) # Fixed effects coefficients
# Plot the fitted regression line over data
plot(log(chromosome.stats$phys.map.length/1000000), sqrt(chromosome.stats$peripherybias_ratio))
abline(a = fixef(pglmm.model)[[1]][1], b = fixef(pglmm.model)[[1]][2])




#----------------------------------------------------------------------------#
# Variance in the broken stick - Skewness ~ Chromosome size ----
#----------------------------------------------------------------------------#

# Data is not normally distributed
# transform with logarithms
hist(log10(phylogenetic_traits$brokenstick_pvariance), breaks = 30)
# Test a model with LMER
lmer.model = lmer(log10(brokenstick_pvariance)~log(phys.map.length) + (1|species), data = lmm.data)
summary(lmer.model)

# Residuals distribution
qqPlot(residuals(lmer.model))
# Residuals as a function of fitted
plot(fitted(lmer.model), residuals(lmer.model))
abline(h = 0)
# Observed~predicted qqplot
plot(lmm.data$phys.map.length[which(!is.na(lmm.data$brokenstick_pvariance))], unlist(predict(lmer.model)))
abline(a = 0, b = 1)

# Model is not working
# Besides, there is no significant effect of chromosome on the variance of the broken stick
# The metric may not be adequate


# # PGLMM
hist(log10(phylogenetic_traits$brokenstick_pvariance), breaks = 30)
pglmm.model = pglmm(log10(brokenstick_pvariance)~log10(phys.map.length) + (1|Species__), data = phylogenetic_traits, family = 'gaussian',
                    cov_ranef = list(Species = tree))
# pglmm.model = pglmm(sqrt(peripherybias_ratio)~log(phys.map.length/1000000) + (1|Species__), data = phylogenetic_traits, family = 'gaussian',
#                     cov_ranef = list(Species = tree))
summary(pglmm.model)
# Residuals distribution
qqPlot(residuals(pglmm.model))
# Residuals as a function of fitted
plot(fitted(pglmm.model), residuals(pglmm.model))
abline(h = 0)
# Observed~predicted qqplot
plot(unlist(pglmm_predicted_values(pglmm.model)), pglmm.model$Y)
abline(a = 0, b = 1)

# fixef(pglmm.model) # Fixed effects coefficients
# Plot the fitted regression line over data
plot(log(chromosome.stats$phys.map.length/1000000), log(chromosome.stats$brokenstick_pvariance))
abline(a = fixef(pglmm.model)[[1]][1], b = fixef(pglmm.model)[[1]][2])





#----------------------------------------------------------------------------#
# Control the influence of map quality ----
#----------------------------------------------------------------------------#




#============================================================================#
# END ----
#============================================================================#