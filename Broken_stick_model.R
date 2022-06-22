########################################################################### #
#                    ECOBIO - PhD
#
#       Broken stick model
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
# GRAPHIC LIBRARIES
library(ggplot2)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
library(ggplotify)

#============================================================================#
# Loading variables & objects ----
#============================================================================#

# Get the directory of the file & set working directory
filedir=dirname(rstudioapi::getSourceEditorContext()$path)
wd=filedir
setwd(wd)

# Custom functions
source("sources/PhysMap.R") # Get the physical positions of a marker given a list of marker names
source("sources/MareyMap.R") # Consolidate a MareyMap input file
source("sources/core_pipeline.R") # Basic functions but essential for the pipeline (e.g. saving logs)
source("sources/stats.R") # Statistical functions implemented

#============================================================================#
# Describing the data with the broken stick model ----
#============================================================================#
# Broken-stick model (White & Hill 2020)
# K segments of equal physical length (Mb) but different genetic length (cM)
# A model with three segments seemed sufficient to describe most patterns of recombination rate/most landscapes (broad-scale variation)
# Show heterogeneity telomere/center
# Show assymetry

# Estimate the proportions of a broken stick model for a single dataset (all chromosomes)
# i.e. proportion of total physical length in k segments of equal genetic length along the chromosome

#-------------------------------------------------------#
# Plots are messy with k = {4, 5, 6}, no clear structure to observe
# Increase k to recover finer structure? Increase to k = 10
# Divide genetic map in 10 segments of equal size
# Color gradient of segments is the relative length to expected proportion
# eq: expected - p
# If > 0, shorter segment than expected, hence high recombination
# Otherwise if < 0, longer segment than expected, hence low recombination
# Now compute for every map the ten proportions
source("sources/MareyMap.R") # Consolidate a MareyMap input file

# List of all maps
list = system(paste("ls ", wd, "/data-cleaned/marey_maps/", sep = ""), intern = TRUE)
list = list[grep("*.txt", list)]
list = list[!(list == "AllMaps.txt")]
# Remove .txt  in filename to get set name properly
list = gsub(".txt", "", list)
# list = list[-c(5, 17)]
# list = list[-c(1:2)]
# process all maps in the list
# lapply(list, batch_brokenstick)
lapply(list, batch_brokenstick2)

# DONE modify batch_brokenstick to get 10 genomic segements of equal size

# Plot it in a STRUCTURE like barplot with chromosome/species
# Add a bar at each null recombination point, i.e. presumed centromere
#-------------------------------------#
### Bar plot representation ###
# k = 10 and fill colour gradient as a function of the departure to the exprected proportion 1/k
# brokenstick = read.table("output/brokenstick/brokenstick_k10.txt", header = TRUE, sep = "\t")
brokenstick = read.table("output/brokenstick/brokenstick_k10_v2.txt", header = TRUE, sep = "\t")
brokenstick = brokenstick[!(is.na(brokenstick$p1) | is.na(brokenstick$p2) | is.na(brokenstick$p3) | is.na(brokenstick$p4) | is.na(brokenstick$p5) | is.na(brokenstick$p6) | is.na(brokenstick$p7) | is.na(brokenstick$p8) | is.na(brokenstick$p9) | is.na(brokenstick$p10)),]
brokenstick = melt(brokenstick)
brokenstick = cbind(paste(brokenstick$set, "_", brokenstick$chromosome, sep =""), brokenstick)

# Get species name
brokenstick$set = gsub("_", " ",regmatches(as.character(brokenstick$set), gregexpr("^[A-Z][a-z]+_[a-z]+", as.character(brokenstick$set))))
# Format columns in proper way
colnames(brokenstick)=c("sample", "set","chromosome","segment","proportion.length")

# Set a vector of gradient color
# Value of color is simply expected - p, the departure from the expected proportion 1/k
k = 10
brokenstick$color = (1/k) - brokenstick$proportion.length
# Besides, estimates the ratio observed/expected
# brokenstick$ratio = brokenstick$proportion.length/(1/k)

# Besides, estimates the ratio expected/observed (longer than expected will have lower relative recombination rate)
brokenstick$ratio = brokenstick$proportion.length/(1/k)
save(brokenstick, file = "output/brokenstick/brokenstick10_v2.Rda")

var(brokenstick$proportion.length)
var(brokenstick$color)
var(brokenstick$ratio)
# An alternative measure of the variance : the variance of the ratio observed/expected
# Ratio = proportion length / (1/10) 
# So the variance of the ratio is one hundred times the variance of the proportion length



# brokenstick$set = factor(brokenstick$set, levels = c("Arabidopsis thaliana", "Brachypodium distachyon", "Brassica napus",
#                                                        "Citrullus lanatus", "Cucumis sativus",
#                                                        "Malus domestica", "Manihot esculenta", "Oryza sativa", "Populus deltoides", "Populus simonii", 
#                                                        "Prunus mume", "Sesamum indicum", "Theobroma cacao", "Vitis vinifera", "Cucumis melo",
#                                                        "Setaria italica","Coffea canephora",
#                                                        "Arachis duranensis", "Lupinus angustifolus",
#                                                        "Capsicum annuum", "Cenchrus americanus", "Glycine max","Gossypium hirsutum",
#                                                        "Gossypium raimondii", "Phaseolus vulgaris","Solanum lycopersicum",
#                                                        "Solanum tuberosum", "Sorghum bicolor", "Hordeum_vulgare",
#                                                        "Triticum aestivum", "Triticum diccocoides", "Zea mays"
# ))

# TODO CAH classification on broken sitcks to group species, possibly two groups
# Species ordered by Hierarchical clustering
# load("output/classification/h1.clust.Rda")
# h1
# h1$order
# h1$labels[h1$order]
# ngroupes = 6
# plot(h1)
# rect.hclust(h1,k=ngroupes)
# # Clustering in groups
# groupes.h1 = as.factor(cutree(h1,k=ngroupes))
# # Individuals in groups
# names(groupes.h1[which(groupes.h1 == 1)])
# names(groupes.h1[which(groupes.h1 == 2)])
# names(groupes.h1[which(groupes.h1 == 3)])
# names(groupes.h1[which(groupes.h1 == 4)])
# names(groupes.h1[which(groupes.h1 == 5)])
# names(groupes.h1[which(groupes.h1 == 6)])
# # Reorder groups
# brokenstick$set = as.factor(brokenstick$set)
# brokenstick_levels = c(names(groupes.h1[which(groupes.h1 == 1)]),
#                        names(groupes.h1[which(groupes.h1 == 5)]),
#                        names(groupes.h1[which(groupes.h1 == 3)]),
#                        names(groupes.h1[which(groupes.h1 == 2)]),
#                        names(groupes.h1[which(groupes.h1 == 4)]),
#                        names(groupes.h1[which(groupes.h1 == 6)]))
# brokenstick_levels = brokenstick_levels[-1] # remove Aegilops
# brokenstick$set = factor(brokenstick$set, levels = brokenstick_levels)

# Not working very well

# # Reorder by increasing variance in proportions
# brokenspecies = brokenstick[,c(2,5)]
# brokenspecies$set = as.character(brokenspecies$set)
# brokenspecies = aggregate(brokenspecies, by = list(brokenspecies$set), sd)
# brokenspecies = brokenspecies[,-c(2)]
# colnames(brokenspecies) = c("set", "var.proportion")
# # Order by ascending variance in proportions
# brokenspecies = brokenspecies[order(brokenspecies$var.proportion),]
# # Reorder groups
# brokenstick$set = as.factor(brokenstick$set)
# brokenstick$set = factor(brokenstick$set, levels = brokenspecies$set)
# save(brokenstick, file = "output/brokenstick/brokenstick10.Rda")

######### GGPLOT
# Plotting the barplot
load(file = "output/brokenstick/brokenstick10.Rda")

brokenStick_barplot = ggplot(data = brokenstick, aes(x=sample, y=proportion.length, fill = color))+
  geom_bar(stat='identity', width = 1) +
  # scale_fill_manual(values = color) +
  scale_fill_viridis_c() +
  # scale_fill_gradient2() +
  facet_grid(~set, scales = "free", space="free_x") +
  labs(x="Chromosome", y="Proportion of\ntotal physical length", fill="Segment") +
  theme(axis.line = element_blank(),
        # axis.line.x = element_blank(), # No x axis
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=24,hjust = 0.5),
        axis.title.x = element_text(color="black", size=24),
        axis.title.y = element_text(color="black", size=24),
        axis.text=element_text(size=24, colour="black"),
        axis.text.x=element_blank(), # No samples names
        axis.ticks.x=element_blank(), # No x axis
        strip.text=element_text(size=18, colour="black", angle = 90),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=24),
        legend.title=element_text(size=24))
brokenStick_barplot
# ggsave("Figures/brokenStick_barplot_k10_gradient.png",
#        device="png",dpi=320,units="cm",width=70,height=30)

brokenStick_barplot = ggplot(data = brokenstick, aes(x=sample, y=proportion.length, fill = log10(ratio)))+
  geom_bar(stat='identity', width = 1) +
  # scale_fill_manual(values = color) +
  scale_fill_viridis_c(breaks = c(-1, 0, 1), labels = c("-1", "0", "1"), direction = -1, limits = c(-1, 1),
                       values = c(0,0.7,0.8,1), option = "D") +
  # scale_fill_gradient2() +
  facet_grid(~set, scales = "free", space="free_x") +
  labs(x="Chromosome", y="Proportion of\ntotal physical length", fill="Segment") +
  theme(axis.line = element_blank(),
        # axis.line.x = element_blank(), # No x axis
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=24,hjust = 0.5),
        axis.title.x = element_text(color="black", size=24),
        axis.title.y = element_text(color="black", size=24),
        axis.text=element_text(size=24, colour="black"),
        axis.text.x=element_blank(), # No samples names
        axis.ticks.x=element_blank(), # No x axis
        strip.text=element_text(size=18, colour="black", angle = 90),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=24),
        legend.title=element_text(size=24))
brokenStick_barplot


#============================================================================#
# Variance in segment proportions - An additional measure of evenness ----
#============================================================================#
load(file = "output/brokenstick/brokenstick10.Rda")
brokenstick_pvariance = data.frame(sample = unique(brokenstick$sample), species = NA, chromosome = NA,
                                   var.proportion = NA, var.ratio = NA)
for (i in 1:nrow(brokenstick_pvariance)) {
  brokenstick_pvariance$species[i] = as.character(unique(brokenstick$set[which(brokenstick$sample == brokenstick_pvariance$sample[i])]))
  brokenstick_pvariance$chromosome[i] = as.character(unique(brokenstick$chromosome[which(brokenstick$sample == brokenstick_pvariance$sample[i])]))
  brokenstick_pvariance$var.proportion[i] = var(brokenstick$proportion.length[which(brokenstick$sample == brokenstick_pvariance$sample[i])])
  brokenstick_pvariance$var.ratio[i] = var(brokenstick$ratio[which(brokenstick$sample == brokenstick_pvariance$sample[i])])
}
write.table(brokenstick_pvariance, file = "tables/brokenstick_pvariance.txt", col.names = TRUE, row.names = FALSE, sep = ";", quote = FALSE)



#============================================================================#
# Variance in inter-marker distances - An extended measure of the broken stick, at a smaller scale ----
#============================================================================#

# I needed a measure of the variance based on untransformed data
# One measure per chromosome, finally written in chromosome.stats
# Divide genetic distances between markers by their physical distance
# Giving the recombination rate between adjacent markers
# Data is untransformed, yet we can expect an effect of unbalanced sampling
# Dependence to sampling density: we need to check if residuals are correlated to marker density before validating the metric

# A function that computes the inter-marker variance in recombination rates for a given Marey map (i.e. one chromosome)
variance_intermarker = function(set = "", chromosome  = "") {
  # Import the map
  marey_map = read.table(paste("data-cleaned/marey_maps/", set, ".txt", sep = ""), header = TRUE, sep = "\t")
  # Subset chromosome
  marey_map = subset(marey_map, as.character(marey_map$map) == as.character(chromosome))
  # Convert physical distances to Mb
  marey_map$phys = marey_map$phys/1000000
  
  # Order the dataset in ascending genetic distances
  library(dplyr)
  library(tidyr)
  marey_map = arrange(marey_map, gen)

  # Compute inter-marker genetic distance divided by physical distance
  intermarker_recrate = (marey_map$gen[-1] - marey_map$gen[-nrow(marey_map)])/
    (marey_map$phys[-1] - marey_map$phys[-nrow(marey_map)])
  # Errors in ording markers: some markers have negative physical distances
  # Transform negative values to NA
  intermarker_recrate[intermarker_recrate < 0] = NA
  # Some markers are strictly at the same physical position, producing "Inf" values
  # Transform "Inf" to NA
  intermarker_recrate[is.infinite(intermarker_recrate)] = NA
  
  # Estimate the variance
  intermarker_variance = var(intermarker_recrate, na.rm = TRUE)
  
  return(intermarker_variance)
}

# test on Arabidospis
data = "Arabidopsis_thaliana_Serin2017"
chr = 1
variance_intermarker(data, chr)
variance_intermarker(set, chromosome)

# Re-iterates over all chromosomes
library(pbmcapply)
varbetweenmarkers = unlist(pbmclapply(1:nrow(chromosome.stats), function(x) variance_intermarker(chromosome.stats$set[x],
                                                                  chromosome.stats$chromosome[x])))
hist(varbetweenmarkers)
hist(log10(varbetweenmarkers))
plot(chromosome.stats$phys.map.length/1000000, varbetweenmarkers, log = "xy")


# We see no notable effect in the variance
# Indeed, the variance depends on the mean recombination rate that is also highly variable between maps
# Divide the variance by the mean recombination rate to get the coefficient of variation
# A function that computes the inter-marker coefficient of variation in recombination rates for a given Marey map (i.e. one chromosome)
cvariation_intermarker = function(set = "", chromosome  = "") {
  # Import the map
  marey_map = read.table(paste("data-cleaned/marey_maps/", set, ".txt", sep = ""), header = TRUE, sep = "\t")
  # Subset chromosome
  marey_map = subset(marey_map, as.character(marey_map$map) == as.character(chromosome))
  # Convert physical distances to Mb
  marey_map$phys = marey_map$phys/1000000
  
  # Order the dataset in ascending genetic distances
  library(dplyr)
  library(tidyr)
  marey_map = arrange(marey_map, gen)
  
  # Compute inter-marker genetic distance divided by physical distance
  intermarker_recrate = (marey_map$gen[-1] - marey_map$gen[-nrow(marey_map)])/
    (marey_map$phys[-1] - marey_map$phys[-nrow(marey_map)])
  # Errors in ording markers: some markers have negative physical distances
  # Transform negative values to NA
  intermarker_recrate[intermarker_recrate < 0] = NA
  # Some markers are strictly at the same physical position, producing "Inf" values
  # Transform "Inf" to NA
  intermarker_recrate[is.infinite(intermarker_recrate)] = NA
  
  # Estimate the variance
  intermarker_variance = var(intermarker_recrate, na.rm = TRUE)/mean(intermarker_recrate, na.rm = TRUE)
  
  return(intermarker_variance)
}

# test on Arabidospis
data = "Arabidopsis_thaliana_Serin2017"
chr = 1
cvariation_intermarker(data, chr)
cvariation_intermarker(set, chromosome)

# Re-iterates over all chromosomes
library(pbmcapply)
cvbetweenmarkers = unlist(pbmclapply(1:nrow(chromosome.stats), function(x) cvariation_intermarker(chromosome.stats$set[x],
                                                                                                 chromosome.stats$chromosome[x])))
hist(cvbetweenmarkers)
hist(log10(cvbetweenmarkers))
plot(chromosome.stats$phys.map.length/1000000, cvbetweenmarkers, log = "xy")


# Store the result in chromosome.stats
chromosome.stats = read.table(file = "tables/chromosome.stats.csv", header = TRUE, sep = ";")
chromosome.stats$varbetweenmarkers = varbetweenmarkers
chromosome.stats$cvbetweenmarkers = cvbetweenmarkers
write.table(chromosome.stats, file = "tables/chromosome.stats.csv", sep = ";", col.names = TRUE, row.names = FALSE, quote = FALSE)















#============================================================================#
# Test the randomness of each recombination landscape by simulating the null hypothesis ----
#============================================================================#
# Recompute with 20 bins
# List of all maps
# list = system(paste("ls ", wd, "/data-cleaned/marey_maps/", sep = ""), intern = TRUE)
# list = list[grep("*.txt", list)]
# list = list[!(list == "AllMaps.txt")]
# # Remove .txt  in filename to get set name properly
# list = gsub(".txt", "", list)
# # list = list[-c(5, 17)]
# list = list[-c(1:2)]
# # process all maps in the list
# source("sources/MareyMap.R") # Consolidate a MareyMap input file
# lapply(list, function(x) batch_brokenstick(x, k = 20))
# 
# brokenstick = read.table("output/brokenstick/brokenstick_k20.txt", header = TRUE, sep = "\t")
# brokenstick = brokenstick[!(is.na(brokenstick$p1) | is.na(brokenstick$p2) | is.na(brokenstick$p3) | is.na(brokenstick$p4) | is.na(brokenstick$p5) | is.na(brokenstick$p6) | is.na(brokenstick$p7) | is.na(brokenstick$p8) | is.na(brokenstick$p9) | is.na(brokenstick$p10)),]
# brokenstick = melt(brokenstick)
# brokenstick = cbind(paste(brokenstick$set, "_", brokenstick$chromosome, sep =""), brokenstick)
# 
# # Get species name
# brokenstick$set = gsub("_", " ",regmatches(as.character(brokenstick$set), gregexpr("^[A-Z][a-z]+_[a-z]+", as.character(brokenstick$set))))
# # Format columns in proper way
# colnames(brokenstick)=c("sample", "set","chromosome","segment","proportion.length")
# 
# # Set a vector of gradient color
# # Value of color is simply expected - p, the departure from the expected proportion 1/k
# k = 20
# brokenstick$color = (1/k) - brokenstick$proportion.length
# # Besides, estimates the ratio observed/expected
# brokenstick$ratio = brokenstick$proportion.length/(1/k)
# save(brokenstick, file = "output/brokenstick/brokenstick20.Rda")
# 
# 
# load(file = "output/brokenstick/brokenstick20.Rda")
load(file = "output/brokenstick/brokenstick10.Rda")

brokenstick_varproportion = data.frame(sample = unique(brokenstick$sample), species = NA, chromosome = NA, var.proportion = NA)
for (i in 1:nrow(brokenstick_varproportion)) {
  brokenstick_varproportion$species[i] = as.character(unique(brokenstick$set[which(brokenstick$sample == brokenstick_varproportion$sample[i])]))
  brokenstick_varproportion$chromosome[i] = as.character(unique(brokenstick$chromosome[which(brokenstick$sample == brokenstick_varproportion$sample[i])]))
  brokenstick_varproportion$var.proportion[i] = var(brokenstick$proportion.length[which(brokenstick$sample == brokenstick_varproportion$sample[i])])
}
hist(brokenstick_varproportion$var.proportion, breaks = 40)


i = 1
i = 644
# Draw k segment relative size with mean = 0.1 and variance = variance in data
k = 10
rnorm(n = k, mean = 1/k, sd = sqrt(brokenstick_varproportion$var.proportion[i]))
draws = rnorm(n = k, mean = 1/k, sd = sqrt(brokenstick_varproportion$var.proportion[i]))
sum(draws)
# Relative sizes - sum equal to 1
draws = draws/sum(draws)
draws
sum(draws)

(obs = brokenstick$proportion.length[which(brokenstick$sample == brokenstick_varproportion$sample[i])])
# Estimate parameters: sd
beta = sd(brokenstick$proportion.length)

# Bootstrap the null hypothesis
nboot =1000
boot = matrix(NA, ncol = k, nrow = nboot)
for (b in 1:nboot) {
  draws = rnorm(n = k, mean = 1/k, sd = beta)
  # Relative sizes - sum equal to 1
  draws = draws/sum(draws)
  boot[b,] = draws
}
# Now we have 'nboot' simulated recombination landscapes
# Compute the theoretical distribution of the variance for mean = 1/k
# bootvar = apply(boot, 1, sd)
# plot(density(bootvar))
# abline(v = beta, col = "red")

bootmean = colMeans(boot, na.rm = TRUE)
quantile.up = apply(boot, 2, function(x) quantile(x, 0.975))
quantile.low = apply(boot, 2, function(x) quantile(x, 0.025))
# Now we have a bootstrapped simulated recombination landscape
bootmean
obs
# quantile.up
# quantile.low

# We compared the two landscapes (i.e. distributions) with a chisq test
?chisq.test
chisq.test(bootmean, obs)
chisq.test(obs, p = rep(0.1, 10))

# Chisq test works with counts and not proportions
# Compare the theoretical simulated proportions to the observed proportions
# See supplementary in @bauerIntraspecificVariationRecombination2013 for the method to compare recombination landscapes with a chisq test
chisq.test(round(rep(bootmean, 2), digits = 3), round(rep(obs, 2), digits = 3))
chisq.test(round(rep(bootmean, 2), digits = 3)*1000, round(rep(obs, 2), digits = 3)*1000)
# Same results with proportions and count data
# Problem of significance is with sample size, only 10 bins is not enough to assess significant differences
# Increase the number of bins in the brokenstick to conclude

# Validity of the chisq test
# Check standardized residuals

#-------------------------------------------------------#
# Work with physical length instead of proportion
# Length follow a normal distribution of mean 0.1*totallength and variance = variance in the data
load(file = "output/brokenstick/brokenstick10.Rda")

i = 661
k = 10
brokenstick$sample[i]
(obs = brokenstick$proportion.length[which(brokenstick$sample == brokenstick$sample[i])])
# Convert to physical length
physlength = chromosome.stats$phys.map.length[which(paste(chromosome.stats$set, chromosome.stats$chromosome, sep = "_") == brokenstick$sample[i])]/1000000
(obs = obs*physlength)
hist(obs)
# Now simulate with bootstrap a theoretical landscape
# Bootstrap the null hypothesis
nboot =1000
boot = matrix(NA, ncol = k, nrow = nboot)
for (b in 1:nboot) {
  draws = rnorm(n = k, mean = physlength/k, sd = sd(obs))
  boot[b,] = draws
}
# Now we have 'nboot' simulated recombination landscapes
# Compute the theoretical distribution of the variance for mean = 1/k
bootvar = apply(boot, 1, sd)
plot(density(bootvar))
abline(v = sd(obs), col = "red")


(exp = colMeans(boot, na.rm = TRUE))
(quantile.up = apply(boot, 2, function(x) quantile(x, 0.975)))
(quantile.low = apply(boot, 2, function(x) quantile(x, 0.025)))
obs
# Compare with a chisq test
chisq.test(obs/sum(obs), exp/sum(exp))

#-------------------------------------------------------#
# In Bauer et al. 2013
# "Each associated variance can be
# obtained by considering that the number of crossovers in a meiosis follows a Poisson
# distribution with mean and variance given by the estimated genetic length in Morgan."

# Genetic length of bins
i = 1
(bingeneticlength = (chromosome.stats$linkage.map.length.correctedHW[which(paste(chromosome.stats$set, chromosome.stats$chromosome, sep = "_") == brokenstick$sample[i])])/k)
# Hence the mean number of COs in the bin and its associated variance

# The simulated landscape of COs number
rpois(n = 10, lambda = bingeneticlength)




#-------------------------------------------------------#
# Bootstrap observed values with replacement to construct H0
# The absolute difference between the theoretical value of 0.1 for segment size and the observed value
i = 1
i = 660
brokenstick$sample[i]
# Draw k segment relative size with mean = 0.1 and variance = variance in data
k = 10
(obs = brokenstick$proportion.length[which(brokenstick$sample == brokenstick$sample[i])])
# Bootstrap the null hypothesis
nboot =1000
boot = matrix(NA, ncol = 1, nrow = nboot)
for (b in 1:nboot) {
  boot[b,] = abs(sample(obs, 1, replace = TRUE) - 0.1)
}
boot
# hist(boot, breaks = 40)
plot(density(boot), xlim = c(0, 0.5))
abline(v = 0)
abline(v = quantile(boot, 0.025), lty = 2)
abline(v = quantile(boot, 0.975), lty = 2)
abline(v = mean(obs-0.1), col = "red")
# Significance testing
# Significant difference when 0 is not comprised in the 95% C.I.
lower = quantile(boot, 0.025)
upper = quantile(boot, 0.975)
(0 < lower | 0 > upper)
# Hence the segment size is significantly different from 0

#-------------------------------------------------------#

i = 1
i = 660
brokenstick$sample[i]
# Draw k segment relative size with mean = 0.1 and variance = variance in data
(obs = brokenstick$proportion.length[which(brokenstick$sample == brokenstick$sample[i])])
# Simulate the chisq stat with 1,000 random draws
# Compare to the observed stat
(stat = sum((obs-0.1)^2/0.1))
pchisq(stat, df = 10) # p-value?

exp = rchisq(n = 100000, df = 10)
plot(density(exp))
abline(v = quantile(exp, 0.025), lty = 2)
abline(v = quantile(exp, 0.975), lty = 2)
abline(v = stat, col = "red")












#============================================================================#
# DEPRECATED ----
#============================================================================#

#-------------------------------------------------------#
# Adopt an approach that maximizes segmentation entropy/contrasts between segments

#============================================================================#
# Finding the segment length that describe the best the data ----
#============================================================================#
# set = "Arabidopsis_thaliana_Serin2017"
# map = read.table(paste("data-cleaned/marey_maps/", set, ".txt", sep =""), header = TRUE)
# chrmap = map[map$map == 1,]

# Rationale
# library(segmented)
# lin.mod <- lm(phys~gen, data = chrmap)
# segmented.mod <- segmented(lin.mod, seg.Z = ~gen, npsi = 2)
# summary(segmented.mod)

# Now compute for every map the three proportions
batch_brokenstick = function(set = "") {
  source("sources/MareyMap.R") # Consolidate a MareyMap input file
  print(set)
  df = read.table("tables/brokenstick/brokenstick_segmented.txt", header = TRUE, sep = "\t")
  df = df[which(df$set != set),]
  res = brokenstick(set, method = "segmented")
  res = cbind(rep(set, length(set)), res)
  colnames(res) = c("set", "chromosome", "p1", "p2", "p3")
  # Problem of factor level
  res$chromosome = as.character(res$chromosome)
  df$chromosome = as.character(df$chromosome)
  # Add results to data & save
  df = rbind(df, res)
  write.table(df, file = "tables/brokenstick/brokenstick_segmented.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
}
# List of all maps
list = system(paste("ls ", wd, "/data-cleaned/marey_maps/", sep = ""), intern = TRUE)
list = list[grep("*.txt", list)]
list = list[!(list == "AllMaps.txt")]
# Remove .txt  in filename to get set name properly
list = gsub(".txt", "", list)
# list = list[-c(5, 17)]
list = list[-1]
# process all maps in the list
lapply(list, batch_brokenstick)


# Plot it in a STRUCTURE like barplot with chromosome/species
# Add a bar at each null recombination point, i.e. presumed centromere
#-------------------------------------#
### Bar plot representation ###
brokenstick = read.table("tables/brokenstick/brokenstick_segmented.txt", header = TRUE, sep = "\t")
# Set a vector of colors for each segment
colors = brewer.pal(6, "Set1")

brokenstick = melt(brokenstick)

brokenstick = cbind(paste(brokenstick$set, "_", brokenstick$chromosome, sep =""), brokenstick)

# Get species name
brokenstick$set = gsub("_", " ",regmatches(as.character(brokenstick$set), gregexpr("^[A-Z][a-z]+_[a-z]+", as.character(brokenstick$set))))

colnames(brokenstick)=c("sample", "set","chromosome","segment","proportion.length")

######### GGPLOT
# Plotting the barplot
brokenStick_barplot = ggplot(data = brokenstick, aes(x=sample, y=proportion.length, fill=segment))+
  geom_bar(stat='identity', width = 1) +
  scale_fill_manual(values = colors) +
  facet_grid(~set, scales = "free", space="free_x") +
  labs(x="Chromosome", y="Proportion of\ntotal physical length", fill="Segment") +
  theme(axis.line = element_blank(),
        # axis.line.x = element_blank(), # No x axis
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=28,hjust = 0.5),
        axis.title.x = element_text(color="black", size=24),
        axis.title.y = element_text(color="black", size=24),
        axis.text=element_text(size=28, colour="black"),
        axis.text.x=element_blank(), # No samples names
        axis.ticks.x=element_blank(), # No x axis
        strip.text=element_text(size=18, colour="black", angle = 90),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28))
brokenStick_barplot
ggsave("Figures/brokenStick_barplot_segmented.png",
       device="png",dpi=320,units="cm",width=70,height=20)



#============================================================================#
# Statistical test of significance - Distribution is not random ----
#============================================================================#
metadata.clean = read.csv(paste(wd, "/data/Genetic_maps/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")

# A little bit of reflexion to find a statistical test of significance
# Simulate a random model (equal size of all segments with n markers)
nmarkers = 10000
rbinom(10, nmarkers/10, 0.1)/(nmarkers/10)

# Random model
nboot = 1000
boot = numeric(nboot)
for (i in 1:nboot) {
  boot[i] = sd(rbinom(10, nmarkers/10, 0.1)/(nmarkers/10))
}
mean(boot)

# Now compare with data
# First with Arabidopsis_thaliana_Serin2017	chromosome 1 ()
# Observed heterogeneity (sd)
obs_sd = sd(brokenstick[1,3:12])

# Expected heterogeneity
nmarkers = sum(data.final$set == "Arabidopsis_thaliana_Serin2017" & data.final$map == 1)
# Random model
nboot = 1000
boot = numeric(nboot)
for (i in 1:nboot) {
  boot[i] = sd(rbinom(10, round(nmarkers/10), 0.1)/(round(nmarkers)/10))
}
mean(boot)
exp_sd = mean(boot)
exp_sd_lower = quantile(boot, 0.025)
exp_sd_upper = quantile(boot, 0.975)

# Is Obs_sd within 95% CI of expected sd?
(obs_sd >= exp_sd_lower & obs_sd <= exp_sd_upper)
hist(boot)
abline(v = obs_sd, col = "Red")

#============================================================================#
# Finding the number of segments that describe the best the data ----
#============================================================================#
# Finding structure
# Decompose heterogeneity in segments maximizing intra-segment homogeneity while maximing inter-segment heterogeneity





#============================================================================#
# OLDER STUFF ----
#============================================================================#

# Load a map: Arabidopsis and Zea to compare contrasting patterns
# set = "Arabidopsis_thaliana_Serin2017"
# set = "Zea_mays_IBM_MaizeSNP50"
# set = "Aegilops_speltoides_Zhang2019"
# brokenstick(set)

# # Now compute for every map the three proportions
# batch_brokenstick = function(set = "") {
#   source("sources/MareyMap.R") # Consolidate a MareyMap input file
#   print(set)
#   df = read.table("tables/brokenstick/brokenstick.txt", header = TRUE, sep = "\t")
#   df = df[which(df$set != set),]
#   res = brokenstick(set, method = "strict")
#   res = cbind(rep(set, length(set)), res)
#   colnames(res) = c("set", "chromosome", "p1", "p2", "p3")
#   # Problem of factor level
#   res$chromosome = as.character(res$chromosome)
#   df$chromosome = as.character(df$chromosome)
#   # Add results to data & save
#   df = rbind(df, res)
#   write.table(df, file = "tables/brokenstick/brokenstick.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
# }
# # List of all maps
# list = system(paste("ls ", wd, "/data-cleaned/marey_maps/", sep = ""), intern = TRUE)
# list = list[grep("*.txt", list)]
# list = list[!(list == "AllMaps.txt")]
# # Remove .txt  in filename to get set name properly
# list = gsub(".txt", "", list)
# # list = list[-c(5, 17)]
# # Remove Aegilops
# list = list[-1]
# # process all maps in the list
# lapply(list, batch_brokenstick)
# 
# 
# # Plot it in a STRUCTURE like barplot with chromosome/species
# # Add a bar at each null recombination point, i.e. presumed centromere
# #-------------------------------------#
# ### Bar plot representation ###
# brokenstick = read.table("tables/brokenstick/brokenstick.txt", header = TRUE, sep = "\t")
# brokenstick = brokenstick[!(is.na(brokenstick$p1) | is.na(brokenstick$p2) | is.na(brokenstick$p3)),]
# # Set a vector of colors for each segment
# colors = brewer.pal(6, "Set1")
# 
# brokenstick = melt(brokenstick)
# 
# brokenstick = cbind(paste(brokenstick$set, "_", brokenstick$chromosome, sep =""), brokenstick)
# 
# # Get species name
# brokenstick$set = gsub("_", " ",regmatches(as.character(brokenstick$set), gregexpr("^[A-Z][a-z]+_[a-z]+", as.character(brokenstick$set))))
# colnames(brokenstick)=c("sample", "set","chromosome","segment","proportion.length")
# 
# brokenstick$set = factor(brokenstick$set, levels = c("Arabidopsis thaliana", "Brachypodium distachyon", "Brassica napus",
#                                                        "Citrullus lanatus", "Cucumis sativus",
#                                                        "Malus domestica", "Manihot esculenta", "Oryza sativa", "Populus deltoides", "Populus simonii", 
#                                                        "Prunus mume", "Sesamum indicum", "Theobroma cacao", "Vitis vinifera", "Cucumis melo",
#                                                        "Setaria italica","Coffea canephora",
#                                                        "Arachis duranensis", "Lupinus angustifolus",
#                                                        "Capsicum annuum", "Cenchrus americanus", "Glycine max","Gossypium hirsutum",
#                                                        "Gossypium raimondii", "Phaseolus vulgaris","Solanum lycopersicum",
#                                                        "Solanum tuberosum", "Sorghum bicolor", "Hordeum_vulgare",
#                                                        "Triticum aestivum", "Triticum diccocoides", "Zea mays"
# ))
# # Add classification to cluster similar species
# # Manual clustering
# # brokenstick$set = factor(brokenstick$set, levels = c("Arabidopsis thaliana", "Malus domestica",
# #                                                        "Manihot esculenta", "Oryza sativa",
# #                                                        "Prunus mume", "Sesamum indicum",
# #                                                        
# #                                                        ))
# # factor(temp$type, levels=c("T","F","P"), labels=c("T","F","P"))
# 
# ######### GGPLOT
# # Plotting the barplot
# k = 3
# brokenStick_barplot = ggplot(data = brokenstick, aes(x=sample, y=proportion.length, fill=segment))+
#   geom_bar(stat='identity', width = 1) +
#   geom_hline(yintercept = c(1/k, 2/k)) +
#   scale_fill_manual(values = colors) +
#   facet_grid(~set, scales = "free", space="free_x") +
#   labs(x="Chromosome", y="Proportion of\ntotal physical length", fill="Segment") +
#   theme(axis.line = element_blank(),
#         # axis.line.x = element_blank(), # No x axis
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
#         plot.subtitle = element_text(color="black",size=28,hjust = 0.5),
#         axis.title.x = element_text(color="black", size=24),
#         axis.title.y = element_text(color="black", size=24),
#         axis.text=element_text(size=28, colour="black"),
#         axis.text.x=element_blank(), # No samples names
#         axis.ticks.x=element_blank(), # No x axis
#         strip.text=element_text(size=18, colour="black", angle = 90),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.text=element_text(size=28),
#         legend.title=element_text(size=28))
# brokenStick_barplot
# ggsave("Figures/brokenStick_barplot.png",
#        device="png",dpi=320,units="cm",width=70,height=20)


# Add classification to cluster similar species
# CAH clustering

# #-------------------------------------------------------#
# # Problem is this method create artefacts with exagerately long segments
# # Increase k to recover finer structure?
# # Now compute for every map the kth proportions
# batch_brokenstick = function(set = "") {
#   source("sources/MareyMap.R") # Consolidate a MareyMap input file
#   print(set)
#   df = read.table("tables/brokenstick/brokenstick_k5.txt", header = TRUE, sep = "\t")
#   df = df[which(df$set != set),]
#   res = brokenstick(set, k = 5, method = "strict")
#   resnames = colnames(res)
#   res = cbind(rep(set, length(set)), res)
#   colnames(res) = c("set", resnames)
#   # Problem of factor level
#   res$chromosome = as.character(res$chromosome)
#   df$chromosome = as.character(df$chromosome)
#   # Add results to data & save
#   df = rbind(df, res)
#   write.table(df, file = "tables/brokenstick/brokenstick_k5.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
# }
# # List of all maps
# list = system(paste("ls ", wd, "/data-cleaned/marey_maps/", sep = ""), intern = TRUE)
# list = list[grep("*.txt", list)]
# list = list[!(list == "AllMaps.txt")]
# # Remove .txt  in filename to get set name properly
# list = gsub(".txt", "", list)
# # list = list[-c(5, 17)]
# list = list[-1]
# # process all maps in the list
# lapply(list, batch_brokenstick)
# 
# 
# # Plot it in a STRUCTURE like barplot with chromosome/species
# # Add a bar at each null recombination point, i.e. presumed centromere
# #-------------------------------------#
# ### Bar plot representation ###
# brokenstick = read.table("tables/brokenstick/brokenstick_k5.txt", header = TRUE, sep = "\t")
# # Set a vector of colors for each segment
# colors = brewer.pal(6, "Set1")
# 
# brokenstick = melt(brokenstick)
# 
# brokenstick = cbind(paste(brokenstick$set, "_", brokenstick$chromosome, sep =""), brokenstick)
# 
# # Get species name
# brokenstick$set = gsub("_", " ",regmatches(as.character(brokenstick$set), gregexpr("^[A-Z][a-z]+_[a-z]+", as.character(brokenstick$set))))
# 
# colnames(brokenstick)=c("sample", "set","chromosome","segment","proportion.length")
# 
# ######### GGPLOT
# # Plotting the barplot
# brokenStick_barplot = ggplot(data = brokenstick, aes(x=sample, y=proportion.length, fill=segment))+
#   geom_bar(stat='identity', width = 1) +
#   scale_fill_manual(values = colors) +
#   facet_grid(~set, scales = "free", space="free_x") +
#   labs(x="Chromosome", y="Proportion of\ntotal physical length", fill="Segment") +
#   theme(axis.line = element_blank(),
#         # axis.line.x = element_blank(), # No x axis
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
#         plot.subtitle = element_text(color="black",size=28,hjust = 0.5),
#         axis.title.x = element_text(color="black", size=24),
#         axis.title.y = element_text(color="black", size=24),
#         axis.text=element_text(size=28, colour="black"),
#         axis.text.x=element_blank(), # No samples names
#         axis.ticks.x=element_blank(), # No x axis
#         strip.text=element_text(size=18, colour="black", angle = 90),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.text=element_text(size=28),
#         legend.title=element_text(size=28))
# brokenStick_barplot
# ggsave("Figures/brokenStick_barplot_k5.png",
#        device="png",dpi=320,units="cm",width=70,height=20)
# 





#============================================================================#
# END  ----
#============================================================================#