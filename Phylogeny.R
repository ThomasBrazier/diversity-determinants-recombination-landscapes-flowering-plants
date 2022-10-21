############################################################################ #
#                    ECOBIO - PhD
#
#       Phylogeny of recombination landscapes
#
############################################################################ #
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
library(rotl) # API for the iTOL
# BiocManager::install("treeio")
library(treeio)
# BiocManager::install("ggtree")
library(ggtree)
library(MareyMap)
library(stringr) 
library(tibble)
library(ggnewscale)
library(phyr)
library(lme4)
library(car)
library(dplyr)
library(ape)

# Custom functions
# source("sources/PhysMap.R") # Get the physical positions of a marker given a list of marker names
# source("sources/MareyMap.R") # Consolidate a MareyMap input file
source("sources/core_pipeline.R") # Basic functions but essential for the pipeline (e.g. saving logs)


#============================================================================#
# Loading variables & objects ----
#============================================================================#

# Get the directory of the file & set working directory
wd=dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

#============================================================================#
#       Building the phylogenetic tree data ----
#============================================================================#

# Import dataset
metadata = read.csv(paste(wd, "/data/Genetic_maps/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")
metadata.clean = read.csv(paste(wd, "/data-cleaned/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")
# Load genome database
genome = read.table(paste(wd, "/data/Genome/Genome_ressources.csv", sep = ""), header = TRUE, sep = ";")
# Load the joint database with all maps
data.final = read.table(file = "data-cleaned/marey_maps/AllMaps.txt", header = TRUE, sep = "\t")
chromosome.stats = read.table("tables/chromosome.stats.csv", header = TRUE, sep = ";")

# Import the big_seed_plant_trees megaphylogeny
# library(ape)
tree = read.tree(file = "data-cleaned/phylogeny/ALLMB.tre")
# Reduce the tree to the species list we had
# List of all validated genetic maps
list.tree = unique(as.character(metadata.clean$species))
(list.tree = gsub("Quercus_sp", "Quercus_robur", list.tree))
# Retrieve species, including varieties
(tips = unlist(lapply(list.tree, function(x) tree$tip.label[grep(x, tree$tip.label)])))
# If the species exists, remove all duplicated varieties
for (species in tree$tip.label[(tree$tip.label %in% list.tree)]) {
  cat(species, "\n")
  tips[grep(species, tips)] = species
}
(tips = unique(tips))
# Keep only one varieties per species (the first one, arbitrary)
for (species in list.tree[!(list.tree %in% tree$tip.label[(tree$tip.label %in% list.tree)])]) {
  tips[grep(species, tips)] = tips[grep(species, tips)][1]
}
(tips = unique(tips))
# Drop all other tips
tree = drop.tip(tree, tree$tip.label[!(tree$tip.label %in% tips)])
tree$tip.label
# Now we have 57 species
# Corrections for names misspelled
# Remove variety names
(tree$tip.label = gsub("_subsp._[A-Za-z]*$", "", tree$tip.label))
(tree$tip.label = gsub("_var._[A-Za-z]*$", "", tree$tip.label))
(tree$tip.label = gsub("_robur", "_sp", tree$tip.label))

# Remove underscores
(tree$tip.label = gsub("_", " ", tree$tip.label))

# Save the tree
save(tree, file = "output/phylogeny/tree.Rdata")
load("output/phylogeny/tree.Rdata")



#============================================================================#
# Retrieve the iTOL phylogenetic tree for species in our dataset
# Matching names in iTOL
(resolved_names = tnrs_match_names(list.tree))
resolved_names$unique_name = resolved_names$search_string

save(resolved_names, file = "output/phylogeny/resolved_names_reducedtree.Rdata")
load("output/phylogeny/resolved_names_reducedtree.Rdata")

tree = tol_induced_subtree(ott_ids = ott_id(resolved_names))

tips = tree$tip.label
tips = gsub("_ott[0-9]*", "", tips)
tips = gsub("_\\(species_in_domain_Eukaryota\\)", "", tips)
tree$tip.label = gsub("_", " ", tips)
tree$tip.label[which(tree$tip.label == "Quercus robur")] = "Quercus sp"
tree$tip.label[which(tree$tip.label == "Oryza rufipogon")] = "Oryza nivara"


#============================================================================#
# PHYLOGENETIC TREE ----
#============================================================================#


#============================================================================#
#       Making a circular phylogenetic tree ----
#============================================================================#
#http://yulab-smu.top/treedata-book/chapter4.html
#https://bioconductor.riken.jp/packages/3.4/bioc/vignettes/ggtree/inst/doc/treeAnnotation.html

# Load the phylogenetic tree estimated in 'Phylogeny.R'
load("output/phylogeny/tree.Rdata")
tree$tip.label
# Add annotations
# e.g. mean recombination rate, mean chromosome size
annotation = data.frame(species = tree$tip.label)
annotation$REC = NA
annotation$CHRSIZE = NA
for (i in 1:nrow(annotation)) {
  annotation$REC[i] = log10(mean(chromosome.stats$mean.recrate[which(chromosome.stats$species == annotation$species[i])], na.rm = TRUE))
  annotation$CHRSIZE[i] = log10(mean(chromosome.stats$phys.map.length[which(chromosome.stats$species == annotation$species[i])], na.rm = TRUE)/1000000)
}
rownames(annotation) = annotation$species
annotation = as.data.frame(annotation[,-1])
colnames(annotation) = c("REC", "CHRSIZE")

# Make the tree
p1 = ggtree(tree, aes(),layout = 'fan', open.angle = 30) +
  geom_tiplab(aes(fontface = "italic"), size = 5, offset = 20)
p1

p2 = gheatmap(p1, annotation[,"REC", drop = FALSE], offset=0, width=.05,
         colnames_angle=90, colnames_offset_y = 61, font.size = 4) +
  scale_fill_viridis_c(option = "D", name = "REC")
p2
p3 = p2 + new_scale_fill()
p3 = gheatmap(p3, annotation[,"CHRSIZE", drop = FALSE], offset=8, width=.05,
              colnames_angle=90, colnames_offset_y = 61, font.size = 4) +
  scale_fill_viridis_c(option = "C", name = "CHRSIZE", direction = -1)
p3

ggsave("figures/phylogeny.svg", plot = p3,
       device="svg", dpi = 320, units = "cm", width = 50, height = 50)


#============================================================================#
#       Circular phylogenetic tree with all phylogeny ----
#============================================================================#





ggsave("figures/phylogeny_complete.svg", plot = p3,
       device="svg", dpi = 320, units = "cm", width = 50, height = 50)









#============================================================================#
# EXPLICIT TEST FOR A PHYLOGENETIC SIGNAL ----
# Assess correlations between summary statistics and genome characteristics
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

metadata = read.csv(file = paste(wd, "/data/Genetic_maps/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")
metadata.clean = read.csv(file = paste(wd, "/data-cleaned/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")
# Load the joint database with all maps
data.final = read.table(file = "data-cleaned/marey_maps/AllMaps.txt", header = TRUE, sep = "\t")
chromosome.stats = read.table(file = "tables/chromosome.stats.csv", header = TRUE, sep = ";")


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

#============================================================================#
# (1) ASSESS IN CHROMOSOME POINTS ----
#============================================================================#

#----------------------------------------------------------------------------#
# Build the chromosome traits dataset ----
#----------------------------------------------------------------------------#
# gsub("_", " ", unlist(lapply(chromosome.stats$set, function(x){metadata.clean$species[which(as.character(metadata.clean$id) == as.character(x))]})))
# phylogenetic_traits = cbind(data.frame(Species = gsub("_", " ", unlist(lapply(chromosome.stats$set, function(x){metadata.clean$species[which(as.character(metadata.clean$id) == as.character(x))]})))))
phylogenetic_traits = cbind(data.frame(Species = gsub("_", " ", unlist(lapply(chromosome.stats$set, function(x){metadata.clean$species[which(as.character(metadata.clean$id) == as.character(x))]})))),
                            chromosome.stats)

#----------------------------------------------------------------------------#
# Chomosomes - Test for the phylogenetic signal ----
#----------------------------------------------------------------------------#
# 
# When applied without predictor (independent) variables, PGLMM gives a test for phylogenetic signal.
pglmm.model = pglmm(log(mean.recrate) ~ 1 + (1|Species__), data = phylogenetic_traits, family = 'gaussian',
                    cov_ranef = list(Species = tree), bayes = TRUE)
summary(pglmm.model)
plot_bayes(pglmm.model, sort = TRUE)

# Species random effect significantly different of 0 (variance explained C.I. = 0.09505  0.64216)
# And phylogenetic signal too was significant (variance explained C.I. = 0.11707  0.72653)
# Residuals had a very low variance (0.03261)
# Hence phylogeny should be considered

# ?cor_phylo
# Look for correlations between variables under phylogenetic autocorrelation
phylocor = cor_phylo(~ log(mean.recrate) + cv.recrate + gini, species = unique(phylogenetic_traits$Species), phy = tree,
          data = phylogenetic_traits, max_iter = 10000)
phylocor
phylocor$corrs # Correlation matrix
phylocor$d # The strength of stabilizing selection
phylocor$B # The coefficient and its significance

#----------------------------------------------------------------------------#
# Chomosomes - Is LMM a good alternative/approximation of the PGLMM? ----
#----------------------------------------------------------------------------#
lmm.data = cbind(data.frame(Species = gsub("_", " ", unlist(lapply(chromosome.stats$set, function(x){metadata.clean$species[which(as.character(metadata.clean$id) == as.character(x))]})))),
                 chromosome.stats)
#----------------------------------------------------------------------------#
# Mean recombination rate ~ Chromosome size ----
#----------------------------------------------------------------------------#
lmer.model = lmer(log(mean.recrate) ~ 1 + (1|Species), data = lmm.data)
summary(lmer.model)
AIC(lmer.model)
BIC(lmer.model)
# Compare with AIC of PGLMM
pglmm.model = pglmm(log(mean.recrate) ~ 1 + (1|Species__), data = phylogenetic_traits, family = 'gaussian',
                    cov_ranef = list(Species = tree))
summary(pglmm.model)
pglmm.model$AIC
pglmm.model$BIC
# The PGLMM minimizes AIC, hence a better likelihood/goodness of fit
# But, after correction for the number of parameters, BIC is better for LMM than PGLMM
# Suggesting that LMM is a most parcimonious model

# But is there common predictions/fit fot both models?
lmer.model = lmer(log(mean.recrate) ~ log(phys.map.length/1000000) + (1|Species), data = lmm.data)
pglmm.model = pglmm(log(mean.recrate) ~ log(phys.map.length/1000000) + (1|Species__), data = phylogenetic_traits, family = 'gaussian',
                    cov_ranef = list(Species = tree))
# Diagnostic plots
par(mfrow = c(1, 2))
# Compare residuals - Present LMER and PGLMM side by side
qqPlot(residuals(lmer.model))
qqPlot(residuals(pglmm.model))
# Residuals are very similar

# Residuals as a function of fitted
plot(fitted(lmer.model), residuals(lmer.model))
abline(h = 0)
plot(fitted(pglmm.model), residuals(pglmm.model))
abline(h = 0)
# Once again, results are very similar

# Observed~predicted qqplot
plot(predict(lmer.model), log(lmm.data$mean.recrate))
abline(a = 0, b = 1)
plot(unlist(pglmm_predicted_values(pglmm.model)), pglmm.model$Y)
abline(a = 0, b = 1)
# Once again, results are very similar

par(mfrow = c(1, 1))
# Plot the fitted regression line over data
plot(log(chromosome.stats$phys.map.length/1000000), log(chromosome.stats$mean.recrate))
abline(a = fixef(lmer.model)[1], b = fixef(lmer.model)[2], col = "blue")
abline(a = fixef(pglmm.model)[[1]][1], b = fixef(pglmm.model)[[1]][2], col = "red", lty = "dashed")
# Regression lines are very similar

# We can conclude that, besides PGLMM being more accurate with the phylogenetic structure in data,
# LMER are a good approximation/alternative to PGLMM models that can be computationnally restrictiv
# or intensive, due to the large number of parameters


#============================================================================#
# (2) ASSESS IN RECOMBINATION WINDOWS  ----
#============================================================================#

#----------------------------------------------------------------------------#
# Build the windows traits dataset ----
#----------------------------------------------------------------------------#
# Load the phylogenetic tree estimated in 'Phylogeny.R'
load("output/phylogeny/tree.Rdata")
tree$tip.label
# Import data for ALL chromosomes & species
# Repoduce to all species in dataset
# List of all maps
list = system(paste("ls ", wd, "/data-cleaned/genome/gene_count/", sep = ""), intern = TRUE)
readall = function(x) {
  data = read.table(file = gzfile(paste(wd, "/data-cleaned/genome/gene_count/", x, sep = "")),
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  gendens = read.table(file = gzfile(paste(wd, "/data-cleaned/genome/gene_density/", x, sep = "")),
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  data$gene_density = gendens$gene_density
  data$coding_density = gendens$coding_density
  data$species = gsub("_", " ", rep(gsub("_[A-Za-z0-9]+_chromosome[A-Za-z0-9]*.txt.gz", "", x), nrow(data)))
  data$species = gsub(" MaizeGDBConsensus", "", data$species)
  return(data)
}
df = bind_rows(lapply(list, function(x) readall(x)))
summary(df)
unique(df$species)
tree$tip.label
# Drop all tips without species (43/56 have been conserved in dataset)
tree = drop.tip(tree, tree$tip.label[!(tree$tip.label %in% unique(df$species))])
sort(tree$tip.label)

cbind(unique(df$species), sort(tree$tip.label))

# Remove all NA values
# df = subset(df, !is.na(df$rec.rate))
# hist(df$rec.rate)
# hist(log(df$rec.rate + 0.0001))
# Remove null values otherwise Bayesian INLA does not work
df = subset(df, (df$rec.rate != 0))
hist(log(df$rec.rate), breaks = 20)

#----------------------------------------------------------------------------#
# Windows - Test for the phylogenetic signal ----
#----------------------------------------------------------------------------#
# 
# pglmm.model = pglmm(log(rec.rate) ~ 1 + (1|species__), data = df, family = 'gaussian',
#                     cov_ranef = list(species = tree), bayes = FALSE)
# When applied without predictor (independent) variables, PGLMM gives a test for phylogenetic signal.
# Apply this method on Linux to compute without errors
pglmm.model = pglmm(log(rec.rate) ~ 1 + (1|species__), data = df, family = 'gaussian',
                    cov_ranef = list(species = tree), bayes = TRUE, verbose = TRUE,
                    prior = "pc.prior.auto")
save(pglmm.model, file = "output/phylogeny/pglmm.model_bayes.Rdata", compress = TRUE) # Computation is intensive
load("output/phylogeny/pglmm.model_bayes.Rdata")
summary(pglmm.model)
plot_bayes(pglmm.model, sort = TRUE)

# Species random effect significantly different of 0 (variance explained C.I. = 0.09505  0.64216)
# And phylogenetic signal too was significant (variance explained C.I. = 0.11707  0.72653)
# Residuals had a very low variance (0.03261)
# Hence phylogeny should be considered

# Look for correlations between variables under phylogenetic autocorrelation
df_noNA = subset(df, !is.na(df$gene_count))
phylocor = cor_phylo(~log(rec.rate) + gene_count, species = unique(df$species), phy = tree,
                     data = df_noNA, max_iter = 10000, boot = 1000)
phylocor
phylocor$corrs # Correlation matrix
phylocor$d # The strength of stabilizing selection
phylocor$B # The coefficient and its significance

# Bootstrap CI
boot_ci(phylocor)

#----------------------------------------------------------------------------#
# Windows - Is LMM a good alternative/approximation of the PGLMM? ----
#----------------------------------------------------------------------------#

#----------------------------------------------------------------------------#
# Mean recombination rate ~ Chromosome size ----
#----------------------------------------------------------------------------#
# Subset only 50,000 values to allow computation of PGLMM (otherwise, if sample too large, matrix is uncomputable)
n_subset = 50000
df_subset = df[sample(1:nrow(df), n_subset),]
lmer.model = lmer(log(rec.rate) ~ 1 + (1|species), data = df_subset)
save(pglmm.model, file = "output/phylogeny/lmer.model_windowsdata.Rdata", compress = TRUE) # Computation is intensive
load("output/phylogeny/lmer.model_windowsdata.Rdata")
summary(lmer.model)
AIC(lmer.model)
BIC(lmer.model)
# Compare with AIC of PGLMM
# PGLMM Likelihood is too computationally intensive for macOs Laptop - Try on Linux server
# TODO Compute in Linux ----
pglmm.model = pglmm(log(rec.rate) ~ 1 + (1|species__), data = df_subset, family = 'gaussian',
                    cov_ranef = list(species = tree))
save(pglmm.model, file = "output/phylogeny/pglmm.model_bayes_windowsdata.Rdata", compress =  TRUE) # Computation is intensive
load("output/phylogeny/pglmm.model_bayes_windowsdata.Rdata")
summary(pglmm.model)
pglmm.model$AIC
pglmm.model$BIC
# The PGLMM minimizes AIC, hence a better likelihood/goodness of fit
# But, after correction for the number of parameters, BIC is better for LMM than PGLMM
# Suggesting that LMM is a most parcimonious model

# But is there common predictions/fit fot both models?
lmer.model = lmer(log(rec.rate) ~ gene_count + (1|species), data = df)
pglmm.model = pglmm(log(rec.rate) ~ gene_count + (1|species__), data = df,
                    family = 'gaussian', cov_ranef = list(species = tree))
# Diagnostic plots
par(mfrow = c(1, 2))
# Compare residuals - Present LMER and PGLMM side by side
qqPlot(residuals(lmer.model))
qqPlot(residuals(pglmm.model))
# Residuals are very similar

# Residuals as a function of fitted
plot(fitted(lmer.model), residuals(lmer.model))
abline(h = 0)
plot(fitted(pglmm.model), residuals(pglmm.model))
abline(h = 0)
# Once again, results are very similar

# Observed~predicted qqplot
plot(predict(lmer.model), log(lmm.data$mean.recrate))
abline(a = 0, b = 1)
plot(unlist(pglmm_predicted_values(pglmm.model)), pglmm.model$Y)
abline(a = 0, b = 1)
# Once again, results are very similar

par(mfrow = c(1, 1))
# Plot the fitted regression line over data
plot(log(chromosome.stats$phys.map.length/1000000), log(chromosome.stats$mean.recrate))
abline(a = fixef(lmer.model)[1], b = fixef(lmer.model)[2], col = "blue")
abline(a = fixef(pglmm.model)[[1]][1], b = fixef(pglmm.model)[[1]][2], col = "red", lty = "dashed")
# Regression lines are very similar

# We can conclude that, besides PGLMM being more accurate with the phylogenetic structure in data,
# LMER are a good approximation/alternative to PGLMM models that can be computationnally restrictiv
# or intensive, due to the large number of parameters



#============================================================================#
# END ----
#============================================================================#
