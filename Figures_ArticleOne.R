############################################################################
#                    ECOBIO - PhD
#
#       Generating figures for articles & presentation
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
library(MareyMap)
library(stringr)
library(dplyr) 
library(reshape2)
library(ggdendro)
# library(ggpubr)
library(egg)
library(RColorBrewer)
library(ggplotify)
library(pals)
library(cowplot)
library(grid)
library(dendextend)
library(phyr)
library(ggtree)

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

# Loading data
chromosome.stats = read.table("tables/chromosome.stats.csv", header = TRUE, sep = ";")
metadata = read.csv(file = paste(wd, "/data/Genetic_maps/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")
metadata.clean = read.csv(file = paste(wd, "/data-cleaned/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")


#============================================================================
# Quantitative assessment - Correlation between the chromosome wide recombination rate & the mean estimated recombination rate
# Show the linear relationship between chromosome wide rate and the mean estimated recombination rate as a cross-validation procedure
#============================================================================
png(filename = "figures/article_one/correlation_chrwide_estimations.png", units = "cm", width = 15, height = 10,
    res = 150)
plot(chromosome.stats$chrwide.rate, chromosome.stats$mean.recrate, xlab = "Chromosome-wide recombination rate",
     ylab = "Mean estimated recombination rate", cex = 1.2, cex.lab = 1.2, cex.axis = 1.2)
abline(a = 0, b = 1, col = "Black")
dev.off()

#############################################################################
# MAREY STATISTICS
#############################################################################
# Assess correlations between summary statistics and genome characteristics
# while controlling for phylogeny by fitting a Phylogenetic Generalized Linear Model
#----------------------------------------------------------------------------
# Build the phylogeny
#----------------------------------------------------------------------------
library(phyr)
# Load the phylogenetic tree estimated in 'Phylogeny.R'
load("output/phylogeny/tree.Rdata")

#----------------------------------------------------------------------------
# Build the traits dataset
#----------------------------------------------------------------------------
phylogenetic_traits = cbind(data.frame(Species = gsub("_", " ", unlist(lapply(chromosome.stats$set, function(x){metadata.clean$species[which(as.character(metadata.clean$id) == as.character(x))]})))),
                            chromosome.stats)
phylogenetic_traits$peripherybias_ratio = phylogenetic_traits$peripherybias_ratio + 0.0001
# three figures side by side

#============================================================================
# Figure. Mean recombination rate (log scale) as a function of physical length (log scale) for each chromosome
#============================================================================
df = data.frame(phys.map.length = phylogenetic_traits$phys.map.length, mean.recrate = phylogenetic_traits$mean.recrate, species = phylogenetic_traits$Species)
ncolors = nlevels(as.factor(phylogenetic_traits$Species))
div_col_pals = brewer.pal.info[brewer.pal.info$category == 'div',]
color.species = unlist(mapply(brewer.pal, div_col_pals$maxcolors, rownames(div_col_pals)))[1:ncolors]

#----------------------------------------------------------------------------
# Mean recombination rate ~ Chromosome size
#----------------------------------------------------------------------------
# Linear regression fit
mod = lm(log(chromosome.stats$mean.recrate) ~ log(phys.map.length/1000000), data = chromosome.stats)
# PGLMM
pglmm.model = pglmm(log(mean.recrate)~log(phys.map.length/1000000) + (1|Species), data = phylogenetic_traits, family = 'gaussian',
                  cov_ranef = list(Species = tree))

# Simulate the expected regression line under the assumption of one CO per chromosome (i.e. a genetic length of 50cM)
expectedline = 100/(phylogenetic_traits$phys.map.length/1000000)

meanrec_chrsize = ggplot(data = phylogenetic_traits, aes(x = log(phys.map.length/1000000), y = log(mean.recrate))) +
  geom_point(aes(fill = Species), pch=21) +
  scale_fill_manual(values = color.species) +
  # coord_trans(x = "log10", y = "log10") +
  scale_x_continuous(breaks = log(c(5, 20, 100, 200, 500, 1000)), labels = c(5, 20, 100, 200, 500, 1000), limits = log(c(3, 1000))) +
  scale_y_continuous(breaks = log(c(0.2, 0.5, 2, 4, 8, 16)), labels = as.character(c(0.2, 0.5, 2, 4, 8, 16)), limits = log(c(0.1, 16))) +
  xlab("Chromosome size (Mb, log scale)") + ylab("Mean recombination rate\n(cM/Mb, log scale)") +
  # geom_smooth(method = 'lm', formula = log(expectedline) ~ log(phylogenetic_traits$phys.map.length/1000000), colour = "Red") + # The expected regression line
  geom_line(aes(x = log(phys.map.length/1000000), y = log(50/(phylogenetic_traits$phys.map.length/1000000))), colour = "Red", linetype = "dashed") + # The expected regression line
  geom_line(aes(x = log(phys.map.length/1000000), y = log(100/(phylogenetic_traits$phys.map.length/1000000))), colour = "Red") + # The expected regression line
  geom_line(aes(x = log(phys.map.length/1000000), y = log(150/(phylogenetic_traits$phys.map.length/1000000))), colour = "Red", linetype = "dashed") + # The expected regression line
  geom_abline(intercept = mod$coefficients[1], slope = mod$coefficients[2], linetype = "dashed", size = 1.1) + # The linear regression line
  geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2], size = 1.1) + # The PGLMM regression line
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
meanrec_chrsize
ggsave(meanrec_chrsize,filename = "figures/article_one/correlations_meanrecrate~chromosomesize.png",
       device="png",dpi=320,units="cm",width=12,height=10)

#----------------------------------------------------------------------------
# Plot random slopes per species
(list_species = as.character(unique(phylogenetic_traits$Species)[table(phylogenetic_traits$Species) > 4]))
phylogenetic_traits_subset = subset(phylogenetic_traits, Species %in% list_species)
phylogenetic_traits_subset$Species = as.character(phylogenetic_traits_subset$Species)
phylogenetic_traits_subset$chrsize = log(phylogenetic_traits_subset$phys.map.length/1000000)
phylogenetic_traits_subset$meanreclog = log(phylogenetic_traits_subset$mean.recrate)
# Drop tips in the phylogeny
tree_subset = drop.tip(tree, as.character(unique(phylogenetic_traits$Species)[table(phylogenetic_traits$Species) < 5]))
# PGLMM with phyr does not allow the estimation of random slopes coefficients
# Hence switch to lmer model that is very close to the PGLMM model
library(lme4)
library(sjPlot)
library(sjmisc)
mixed.model = lmer(log(mean.recrate)~log(phys.map.length/1000000) + (log(phys.map.length/1000000)|Species),
                   data = phylogenetic_traits_subset)
# Boundaries of random segments (mean of x values ± e^0.5)
boundaries = data.frame(species = rownames(coef(mixed.model)$Species), upper = NA, lower = NA, color = NA)
chromosome.stats$species = as.character(chromosome.stats$species)
boundaries$species = as.character(boundaries$species)
for (i in 1:nrow(boundaries)) {
  boundaries$upper[i] = log(mean(chromosome.stats$phys.map.length[which(chromosome.stats$species == boundaries$species[i])], na.rm = TRUE)/1000000) + 0.5
  boundaries$lower[i] = log(mean(chromosome.stats$phys.map.length[which(chromosome.stats$species == boundaries$species[i])], na.rm = TRUE)/1000000) - 0.5
  # Add the same colors as in previous ggplot
  boundaries$color[i] = color.species[which(levels(as.factor(phylogenetic_traits$species)) == boundaries$species[i])]
}

meanrec_chrsize_randomslopes = ggplot(data = phylogenetic_traits, aes(x = log(phys.map.length/1000000), y = log(mean.recrate))) +
  geom_point(aes(fill = Species), pch = 21, alpha = 0.3) +
  scale_fill_manual(values = color.species) +
  # coord_trans(x = "log10", y = "log10") +
  scale_x_continuous(breaks = log(c(5, 20, 100, 200, 500, 1000)), labels = c(5, 20, 100, 200, 500, 1000), limits = log(c(3, 1000))) +
  scale_y_continuous(breaks = log(c(0.2, 0.5, 2, 4, 8, 16)), labels = as.character(c(0.2, 0.5, 2, 4, 8, 16)), limits = log(c(0.1, 16))) +
  xlab("Chromosome size (Mb, log scale)") + ylab("Mean recombination rate\n(cM/Mb, log scale)") +
  geom_abline(intercept = mod$coefficients[1], slope = mod$coefficients[2], linetype = "dashed", size = 1.1) + # The linear regression line
  geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2], size = 1.1) + # The PGLMM regression line
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
# Plotting random effects
(ranint = coef(mixed.model)$Species[,1]) # Random intercepts
(ranslope = coef(mixed.model)$Species[,2]) # Random slopes
# Add random slopes
for (i in 1:length(ranint)) {
  meanrec_chrsize_randomslopes = meanrec_chrsize_randomslopes +
    geom_segment(x = boundaries$lower[i], xend = boundaries$upper[i],
                 y = ranint[i] + ranslope[i]*boundaries$lower[i],
                 yend = ranint[i] + ranslope[i]*boundaries$upper[i],
                 colour = boundaries$color[i])
}
ggsave(meanrec_chrsize_randomslopes,filename = "figures/article_one/correlations_meanrecrate_randomslopes~chromosomesize.png",
       device="png",dpi=320,units="cm",width=12,height=10)

# Two figures in one
correlation_2plots_intensity = ggarrange(meanrec_chrsize, meanrec_chrsize_randomslopes, ncol = 2, nrow = 1, labels = c("(a)", "(b)"))
ggsave(correlation_2plots_intensity,filename = "figures/article_one/correlations_2plots_intensity~chromosomesize.png",
       device="png",dpi=320,units="cm",width=36,height=12)


#============================================================================
# Figure. Correlation of random slopes for Mean recombination rate (log scale) as a function of physical length (log scale)
#============================================================================



#============================================================================
# Figure. Coefficient of variation recombination rate (log scale) as a function of physical length (log scale) for each chromosome
#============================================================================
#----------------------------------------------------------------------------
# Recombination rate coefficient of variation ~ Chromosome size
#----------------------------------------------------------------------------
# Linear regression fit
mod = lm(chromosome.stats$cv.recrate ~ log(phys.map.length/1000000), data = chromosome.stats)
# PGLMM
pglmm.model = pglmm(cv.recrate~log(phys.map.length/1000000) + (1| Species__), data = phylogenetic_traits, family = 'gaussian',
                    cov_ranef = list(Species = tree))
# Smoothing enveloppe of the fitted regression line?
# The expected PGLMM regression line with Poisson distribution of COs (sqrt(mean)/mean)
pglmm.model.expected = pglmm((sqrt(mean.recrate)/mean.recrate)~log(phys.map.length/1000000) + (1| Species__), data = phylogenetic_traits, family = 'gaussian',
                             cov_ranef = list(Species = tree))

cvrec_chrsize = ggplot(data = phylogenetic_traits, aes(x = log(phys.map.length/1000000), y = cv.recrate)) +
  geom_point(aes(fill = Species), pch=21) +
  scale_fill_manual(values = color.species) +
  # coord_trans(x = "log10", y = "log10") +
  scale_x_continuous(breaks = log(c(5, 20, 100, 200, 500, 1000)), labels = c(5, 20, 100, 200, 500, 1000), limits = log(c(3, 1000))) +
  scale_y_continuous(breaks = c(1, 2, 3, 4, 5), labels = as.character(c(1, 2, 3, 4, 5))) +
  xlab("Chromosome size (Mb, log scale)") + ylab("Coefficient of variation") +
  geom_abline(intercept = mod$coefficients[1], slope = mod$coefficients[2], linetype = "dashed") + # The linear regression line
  geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
  # geom_abline(intercept = fixef(pglmm.model.expected)[[1]][1], slope = fixef(pglmm.model.expected)[[1]][2], colour = "Red") + # The expected PGLMM regression line
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
cvrec_chrsize
ggsave(cvrec_chrsize,filename = "figures/article_one/correlations_cvrecrate~chromosomesize.png",
       device="png",dpi=320,units="cm",width=12,height=10)

#============================================================================
# Figure. Gini recombination rate (log scale) as a function of physical length (log scale) for each chromosome
#============================================================================
#----------------------------------------------------------------------------
# Gini coefficient ~ Chromosome size
#----------------------------------------------------------------------------
# Linear regression fit
mod = lm(chromosome.stats$gini ~ log(phys.map.length/1000000), data = chromosome.stats)
# PGLMM
pglmm.model = pglmm(gini~log(phys.map.length/1000000) + (1| Species), data = phylogenetic_traits, family = 'gaussian',
                    cov_ranef = list(Species = tree))
# Smoothing enveloppe of the fitted regression line?

gini_chrsize = ggplot(data = phylogenetic_traits, aes(x = log(phys.map.length/1000000), y = gini)) +
  geom_point(aes(fill = Species), pch=21) +
  scale_fill_manual(values = color.species) +
  # coord_trans(x = "log10", y = "log10") +
  scale_x_continuous(breaks = log(c(5, 20, 100, 200, 500, 1000)), labels = c(5, 20, 100, 200, 500, 1000), limits = log(c(3, 1000))) +
  # scale_y_continuous(breaks = c(1, 2, 3, 4, 5), labels = as.character(c(1, 2, 3, 4, 5))) +
  xlab("Chromosome size (Mb, log scale)") + ylab("Gini coefficient") +
  geom_abline(intercept = mod$coefficients[1], slope = mod$coefficients[2], linetype = "dashed") + # The linear regression line
  geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
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
gini_chrsize
ggsave(gini_chrsize,filename = "figures/article_one/correlations_gini~chromosomesize.png",
       device="png",dpi=320,units="cm",width=12,height=10)

# Two figures in one
correlation_2plots_heterogeneity = ggarrange(cvrec_chrsize, gini_chrsize, ncol = 2, nrow = 1, labels = c("(a)", "(b)"))
ggsave(correlation_2plots_heterogeneity,filename = "figures/article_one/correlations_2plots_heterogeneity~chromosomesize.png",
       device="png",dpi=320,units="cm",width=36,height=12)

#============================================================================
# Figure. Periphery-bias ratio ~ physical length (log scale) for each chromosome
#============================================================================
chromosomes = read.table("tables/chromosome.stats.csv", header = TRUE, sep = ";")
# Remove NAs
chromosomes = chromosomes[which(!is.na(chromosomes$peripherybias_ratio)),]
chromosomes$peripherybias_ratio = chromosomes$peripherybias_ratio + 0.0001

mod = lm(sqrt(chromosomes$peripherybias_ratio)~log(chromosomes$phys.map.length/1000000))
pglmm.model = pglmm(sqrt(peripherybias_ratio)~log(phys.map.length/1000000) + (1|species), data = chromosomes, family = 'gaussian',
                    cov_ranef = list(species = tree))

peripherybias_chrsize = ggplot(data = chromosomes, aes(x = log(phys.map.length/1000000), y = sqrt(peripherybias_ratio))) +
  geom_point(aes(fill = species), pch=21) +
  scale_fill_manual(values = color.species) +
  # coord_trans(x = "log10", y = "log10") +
  scale_x_continuous(breaks = sqrt(c(5, 20, 100, 200, 500, 1000)), labels = c(5, 20, 100, 200, 500, 1000), limits = log(c(3, 1000))) +
  # scale_y_continuous(breaks = c(1, 2, 3, 4, 5), labels = as.character(c(1, 2, 3, 4, 5))) +
  xlab("Chromosome size (Mb, log scale)") + ylab("Periphery bias ratio (square root scale)") +
  geom_abline(intercept = mod$coefficients[1], slope = mod$coefficients[2], linetype = "dashed") + # The linear regression line
  geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=8),
        axis.title.y = element_text(color="black", size=8),
        axis.text=element_text(size=8, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
peripherybias_chrsize
ggsave(peripherybias_chrsize,filename = "figures/article_one/correlations_peripherybias~chromosomesize.png",
       device="png",dpi=320,units="cm",width=12,height=6)

#----------------------------------------------------------------------------
# ASSESSING SENSITIVITY TO SAMPLING IN RATIO - HOW MANY BINS?
#----------------------------------------------------------------------------
# Compute a gradient ratio for all bin number
df = read.table("output/dist2telomere/DistancesRelative_bins.txt", header = TRUE)
gradient_ratio = data.frame(set = rep(chromosomes$set, 20), chromosome = rep(chromosomes$chromosome, 20), bin = rep(1:20, each = nrow(chromosomes)), ratio = NA)
for (i in 1:nrow(gradient_ratio)) {
  subs = subset(df, set == gradient_ratio$set[i] & chromosome == gradient_ratio$chromosome[i])
  # gradient_ratio$ratio[i] = (sum(subs$rec.rate[c(1:gradient_ratio$bin[i])], na.rm = TRUE)/sum(subs$rec.rate, na.rm = TRUE))
  # gradient_ratio[i, j] = (mean(subs$rec.rate[c(1:j)], na.rm = TRUE)/mean(subs$rec.rate[c(j:20)], na.rm = TRUE))/mean(subs$rec.rate, na.rm = TRUE) # Standardized ratio by the mean recombination rate
  gradient_ratio$ratio[i] = (mean(subs$rec.rate[c(1:gradient_ratio$bin[i])], na.rm = TRUE)/mean(subs$rec.rate, na.rm = TRUE))
}

# Bootstrap C.I. for loess regression
# library(spatialEco)
# loess.boot = loess.boot(gradient_ratio$bin, gradient_ratio$ratio, nreps = 1000, confidence = 0.95)
# loess.boot
# loess.boot$fit
# plot(loess.boot)
# save(loess.boot, file = "output/gradient_peripherybiasratio_loessboot.Rda")
load("output/gradient_peripherybiasratio_loessboot.Rda")
loess.boot$fit

gradient_ratio_plot = ggplot() +
  geom_point(data = gradient_ratio, aes(x = bin, y = log(ratio)), shape = 1, colour = "darkGrey") +
  # geom_smooth(aes(x = bin, y = log(ratio)), method = "loess", se = TRUE, level = 0.95) +
  geom_ribbon(data = loess.boot$fit, aes(x = x, ymin = log(low.lim), ymax = log(up.lim)), fill = "grey70") +
  geom_line(data = loess.boot$fit, aes(x = x, y = log(y.fit))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab("Number of bins sampled at chromosome end") + ylab("Periphery-bias ratio\n(square root scale)") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        axis.title.x = element_text(color="black", size=9),
        axis.title.y = element_text(color="black", size=9),
        axis.text=element_text(size=9, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.height = unit(2,"line"),
        legend.key.width = unit(5,"line"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        legend.position='none')
gradient_ratio_plot
ggsave(gradient_ratio_plot,filename = "figures/article_one/peripherybias_sensitivity.png",
       device="png",dpi=320,units="cm",width=12,height=6)


#----------------------------------------------------------------------------
# DIAGNOSTIC PLOTS
#----------------------------------------------------------------------------
load("output/phylogeny/tree.Rdata")

#----------------------------------------------------------------------------
# Mean recombination rate ~ Chromosome size
#----------------------------------------------------------------------------
pglmm.model = pglmm(log(mean.recrate)~log(phys.map.length/1000000) + (1|Species), data = phylogenetic_traits, family = 'gaussian',
                    cov_ranef = list(Species = tree))
# Diagnostic plots
png(filename = "figures/article_one/diagnosticplots_pglmm_meanrecrate~chrsize.png", units = "cm", width = 20, height = 8, res = 150)
par(mfrow = c(1,3))
# Residuals distribution
qqPlot(residuals(pglmm.model), xlab = "Theoretical Normal quantiles", ylab = "Residuals")
title("(a)", adj = 0)
# Residuals as a function of fitted
plot(fitted(pglmm.model), residuals(pglmm.model), xlab = "Fitted", ylab = "Residuals")
abline(h = 0, col = "red")
title("(b)", adj = 0)
# Observed~predicted qqplot
plot(unlist(pglmm.predicted.values(pglmm.model)), pglmm.model$Y, xlab = "Observed", ylab = "Predicted")
abline(a = 0, b = 1, col = "red")
title("(c)", adj = 0)
par(mfrow = c(1,1))
dev.off()

#----------------------------------------------------------------------------
# Recombination rate coefficient of variation ~ Chromosome size
#----------------------------------------------------------------------------
pglmm.model = pglmm(cv.recrate~log(phys.map.length/1000000) + (1| Species__), data = phylogenetic_traits, family = 'gaussian',
                    cov_ranef = list(Species = tree))
# Diagnostic plots
png(filename = "figures/article_one/diagnosticplots_pglmm_cvrecrate~chrsize.png", units = "cm", width = 20, height = 8, res = 150)
par(mfrow = c(1,3))
# Residuals distribution
qqPlot(residuals(pglmm.model), xlab = "Theoretical Normal quantiles", ylab = "Residuals")
title("(a)", adj = 0)
# Residuals as a function of fitted
plot(fitted(pglmm.model), residuals(pglmm.model), xlab = "Fitted", ylab = "Residuals")
abline(h = 0, col = "red")
title("(b)", adj = 0)
# Observed~predicted qqplot
plot(unlist(pglmm.predicted.values(pglmm.model)), pglmm.model$Y, xlab = "Observed", ylab = "Predicted")
abline(a = 0, b = 1, col = "red")
title("(c)", adj = 0)
par(mfrow = c(1,1))
dev.off()

#----------------------------------------------------------------------------
# Gini coefficient ~ Chromosome size
#----------------------------------------------------------------------------
pglmm.model = pglmm(gini~log(phys.map.length/1000000) + (1|Species), data = phylogenetic_traits, family = 'gaussian',
                    cov_ranef = list(Species = tree))
# Diagnostic plots
png(filename = "figures/article_one/diagnosticplots_pglmm_gini~chrsize.png", units = "cm", width = 20, height = 8, res = 150)
par(mfrow = c(1,3))
# Residuals distribution
qqPlot(residuals(pglmm.model), xlab = "Theoretical Normal quantiles", ylab = "Residuals")
title("(a)", adj = 0)
# Residuals as a function of fitted
plot(fitted(pglmm.model), residuals(pglmm.model), xlab = "Fitted", ylab = "Residuals")
abline(h = 0, col = "red")
title("(b)", adj = 0)
# Observed~predicted qqplot
plot(unlist(pglmm.predicted.values(pglmm.model)), pglmm.model$Y, xlab = "Observed", ylab = "Predicted")
abline(a = 0, b = 1, col = "red")
title("(c)", adj = 0)
par(mfrow = c(1,1))
dev.off()

#----------------------------------------------------------------------------
# Bias toward periphery - Skewness ~ Chromosome size
#----------------------------------------------------------------------------
pglmm.model = pglmm(sqrt(peripherybias_ratio)~log(phys.map.length/1000000) + (1|Species), data = phylogenetic_traits, family = 'gaussian',
                    cov_ranef = list(Species = tree))
# Diagnostic plots
png(filename = "figures/article_one/diagnosticplots_pglmm_peripherybiasratio~chrsize.png", units = "cm", width = 20, height = 8, res = 150)
par(mfrow = c(1,3))
# Residuals distribution
qqPlot(residuals(pglmm.model), xlab = "Theoretical Normal quantiles", ylab = "Residuals")
title("(a)", adj = 0)
# Residuals as a function of fitted
plot(fitted(pglmm.model), residuals(pglmm.model), xlab = "Fitted", ylab = "Residuals")
abline(h = 0, col = "red")
title("(b)", adj = 0)
# Observed~predicted qqplot
plot(unlist(pglmm.predicted.values(pglmm.model)), pglmm.model$Y, xlab = "Observed", ylab = "Predicted")
abline(a = 0, b = 1, col = "red")
title("(c)", adj = 0)
par(mfrow = c(1,1))
dev.off()






#############################################################################
#     BROKEN STICK
#############################################################################
#============================================================================
# Figures. M&M Broken stick conceptual model
#============================================================================
df = data.frame(distribution = rep(c("Uniform", "Distal", "Proximal"), each = 3),
                segment = rep(c("p1","p2","p3"), times = 3),
                proportion = c(1/3, 1/3, 1/3, 0.2, 0.6, 0.2, 0.4, 0.2, 0.4),
                color = rep(c("black", "grey", "black"), times = 3))
df$distribution = factor(df$distribution, levels = c("Uniform", "Distal", "Proximal"))
df$segment = as.factor(df$segment)
df$rec = 0.3-df$proportion
k3 = ggplot(data = df, aes(x = proportion, y = distribution, group = segment, colour = color, fill = rec))+
  geom_bar(stat='identity', width = 1) +
  facet_grid(rows=df$distribution, scales = "free", space="fixed") +
  scale_colour_manual(values = c("black", "black")) +
  scale_fill_gradient(high = "grey80", low ="grey20") +
  labs(x="Relative physical length (K=3)", y="", fill="Relative recombination rate") +
  theme(axis.line = element_blank(),
        # axis.line.x = element_blank(), # No x axis
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=12, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=12,hjust = 0.5),
        axis.title.x = element_text(color="black", size=12),
        axis.title.y = element_text(color="black", size=12),
        axis.text=element_blank(),
        axis.text.x=element_blank(), # No samples names
        axis.ticks=element_blank(), # No x axis
        strip.text.y=element_text(size=12, colour="black", angle = 0),
        legend.key = element_rect(fill = "white", size = 1),
        strip.background=element_rect(fill="white"),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(size=12, angle = 0),
        legend.title=element_text(size=12),
        legend.position = "none")
k3
df = data.frame(distribution = rep(c("Uniform", "Distal", "Proximal"), each = 10),
                 segment = rep(c("p1","p2","p3", "p4", "p5", "p6", "p7", "p8", "p9", "p10"), times = 3),
                 proportion = c(rep(0.1, times = 10),
                                0.05, 0.05, 0.08, 0.1, 0.22, 0.22, 0.1, 0.08, 0.05, 0.05,
                                0.2, 0.2, 0.05, 0.03, 0.02, 0.02, 0.03, 0.05, 0.2, 0.2),
                 color = rep(c("black"), times = 10*3))
df$distribution = factor(df$distribution, levels = c("Uniform", "Distal", "Proximal"))
df$segment = factor(df$segment, levels = c("p1","p2","p3", "p4", "p5", "p6", "p7", "p8", "p9", "p10"))
df$rec = 0.1-df$proportion
k10 = ggplot(data = df, aes(x = proportion, y = distribution, group = segment, colour = color, fill = rec))+
  geom_bar(stat='identity', width = 1) +
  facet_grid(rows=df$distribution, scales = "free", space="fixed") +
  scale_colour_manual(values = c("black")) +
  scale_fill_gradient(high = "grey80", low ="grey20") +
  labs(x="Relative physical length (K=10)", y="", fill="Relative recombination rate") +
  theme(axis.line = element_blank(),
        # axis.line.x = element_blank(), # No x axis
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=12, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=12,hjust = 0.5),
        axis.title.x = element_text(color="black", size=12),
        axis.title.y = element_text(color="black", size=12),
        axis.text=element_blank(),
        axis.text.x=element_blank(), # No samples names
        axis.ticks=element_blank(), # No x axis
        strip.text.y=element_text(size=12, colour="black", angle = 0),
        legend.key = element_rect(fill = "white", size = 1),
        strip.background=element_rect(fill="white"),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(size=12, angle = 0),
        legend.title=element_text(size=12),
        legend.position = "none")
k10
ggarrange(k3, k10, nrow = 2, ncol = 1)
ggsave(filename = "figures/article_one/BrokenStick_concept.png", unit = "cm", dpi = 320, width = 15, height = 5)




#============================================================================
# Figures. Broken stick barplot solo
#============================================================================

#----------------------------------------------------------------------------
#       Making a vertical brokenstick
#----------------------------------------------------------------------------
# ordered by the phylogeny
load("output/phylogeny/tree_reduced.Rdata")
tree$tip.label
load(file = "output/brokenstick/brokenstick10.Rda")

brokenspecies = brokenstick[,c(2,5)]
brokenspecies$set = as.character(brokenspecies$set)
brokenspecies = aggregate(brokenspecies, by = list(brokenspecies$set), sd)
brokenspecies = brokenspecies[,-c(2)]
colnames(brokenspecies) = c("set", "var.proportion")
# Order by ascending variance in proportions
brokenspecies = brokenspecies[order(brokenspecies$var.proportion),]
# Reorder groups
brokenstick$set = as.factor(brokenstick$set)
brokenstick$set = factor(brokenstick$set, levels = brokenspecies$set)

# Plotting the barplot
p1 = ggplot(data = brokenstick, aes(x = sample, y = proportion.length, fill = color))+
  geom_bar(stat='identity', width = 1) +
  # scale_fill_manual(values = color) +
  scale_fill_viridis_c(breaks = c(-0.5, 0), labels = c("-0.5", "0"), option = "viridis") +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  # scale_fill_gradient2() +
  facet_grid(~set, scales = "free", space="fixed") +
  labs(x="Chromosome", y="\nRelative physical length\n", fill="Relative recombination rate") +
  theme(axis.line = element_blank(),
        # axis.line.x = element_blank(), # No x axis
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="black", angle = 0),
        axis.text.x=element_blank(), # No samples names
        axis.ticks.x=element_blank(), # No x axis
        strip.text=element_text(size=18, colour="black", angle = 90),
        legend.key = element_rect(fill = "white", size = 1),
        strip.background=element_rect(fill="white"),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(size=16, angle = 0),
        legend.title=element_text(size=16),
        legend.position = "bottom")
# p1
ggsave(filename = "figures/article_one/BrokenStick_solo.png", unit = "cm", dpi = 320, width = 60, height = 20)


#============================================================================
# Figures. Broken stick barplot annotated with other summary statistics
#============================================================================

#----------------------------------------------------------------------------
#       Making a vertical brokenstick
#----------------------------------------------------------------------------
# ordered by the phylogeny
load("output/phylogeny/tree.Rdata")
tree$tip.label
load(file = "output/brokenstick/brokenstick10.Rda")

brokenspecies = brokenstick[,c(2,5)]
brokenspecies$set = as.character(brokenspecies$set)
brokenspecies = aggregate(brokenspecies, by = list(brokenspecies$set), sd)
brokenspecies = brokenspecies[,-c(2)]
colnames(brokenspecies) = c("set", "var.proportion")
# Order by ascending variance in proportions
brokenspecies = brokenspecies[order(brokenspecies$var.proportion),]
# Reorder groups
brokenstick$set = as.factor(brokenstick$set)
brokenstick$set = factor(brokenstick$set, levels = brokenspecies$set)

# Plotting the barplot
p1 = ggplot(data = brokenstick, aes(x = sample, y = proportion.length, fill = color))+
  geom_bar(stat='identity', width = 1) +
  # scale_fill_manual(values = color) +
  scale_fill_viridis_c(breaks = c(-0.5, 0), labels = c("-0.5", "0"), option = "viridis") +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  # scale_fill_gradient2() +
  facet_grid(~set, scales = "free", space="fixed") +
  labs(x="", y="\nRelative physical length\n", fill="") +
  theme(axis.line = element_blank(),
        # axis.line.x = element_blank(), # No x axis
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16, angle = -90),
        axis.text=element_text(size=16, colour="black", angle = 0),
        axis.text.x=element_blank(), # No samples names
        axis.ticks.x=element_blank(), # No x axis
        strip.text=element_text(size=18, colour="black", angle = -90),
        legend.key = element_rect(fill = "white", size = 1),
        strip.background=element_rect(fill="white"),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(size=16, angle = 0),
        legend.title=element_text(size=16),
        legend.position = "right",
        legend.direction = "horizontal")
# p1
#----------------------------------------------------------------------------
# Adding covariates
#----------------------------------------------------------------------------
covariates = subset(brokenstick, brokenstick$segment == "p1")
covariates$proportion.length = rep(1, nrow(covariates))
covariates = covariates[,-c(4, 6)]
covariates$chr_size = NA
chromosome.stats$set = as.character(chromosome.stats$set)
chromosome.stats$chromosome = as.character(chromosome.stats$chromosome)
for (i in 1:nrow(covariates)) {
  covariates$chr_size[i] = chromosome.stats$phys.map.length[chromosome.stats$set == gsub("_[0-9A-Za-z]*$", "",covariates$sample[i]) & chromosome.stats$chromosome == covariates$chromosome[i]]/1000000
}
covariates$avgrecrate = NA
for (i in 1:nrow(covariates)) {
  covariates$avgrecrate[i] = chromosome.stats$mean.recrate[chromosome.stats$set == gsub("_[0-9A-Za-z]*$", "",covariates$sample[i]) & chromosome.stats$chromosome == covariates$chromosome[i]]
}
# rbarintra = read.table("tables/df_geneticshuffling_physical.csv", header = TRUE)
# brokenstick$rbarintra = NA
# for (i in 1:nrow(brokenstick)) {
#   brokenstick$rbarintra[i] = rbarintra$geneticshuffling[rbarintra$set == gsub("_[0-9A-Za-z]*$", "",brockenstick$sample[i])]
# }
p2 = ggplot(data = covariates, aes(x = sample, y = proportion.length, fill = chr_size)) +
  geom_bar(stat='identity', width = 1) +
  # scale_fill_manual(values = color) +
  scale_fill_viridis_c(breaks = c(20, 50, 500), labels = c("20", "50", "500"), trans = "log", option = "inferno", direction = -1) +
  # scale_fill_gradient2() +
  facet_grid(~set, scales = "free", space="fixed") +
  labs(x = "", y = "\n", fill = "") +
  theme(axis.line = element_blank(),
        # axis.line.x = element_blank(), # No x axis
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="white"),
        axis.text.x=element_blank(), # No samples names
        axis.ticks=element_blank(), # No axis ticks
        strip.text=element_blank(),
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(size=16, angle = 0),
        legend.title=element_text(size=16),
        legend.position = "right",
        legend.direction = "horizontal",
        plot.margin = margin(0, 0, 0, 0, "cm"))
# p2
p3 = ggplot(data = covariates, aes(x = sample, y = proportion.length, fill = avgrecrate)) +
  geom_bar(stat='identity', width = 1) +
  scale_fill_viridis_c(breaks = c(0, 1, 10), labels = c("0", "1", "10"), trans = "log", option = "inferno") +
  facet_grid(~set, scales = "free", space="fixed") +
  labs(x="Chromosome", y="\n", fill = "") +
  theme(axis.line = element_blank(),
        # axis.line.x = element_blank(), # No x axis
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="white"),
        axis.text.x=element_blank(), # No samples names
        axis.ticks=element_blank(), # No axis ticks
        strip.text=element_blank(),
        legend.position = "right",
        legend.direction = "horizontal",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(size=16, angle = 0),
        legend.title=element_text(size=16),
        plot.margin = margin(0, 0, 0, 0, "cm"))
p3

# Aligning complex ggplots
# gridExtra::grid.arrange provides no way to align the panels of individual plots. While this is achievable with low-level gtable functions, it often requires substantial effort on a case-by-case basis. The egg package introduces a general strategy for such layout manipulations, with the following steps:
# decompose each plot into a 3x3 layout, where the central cell corresponds to the core panels, surrounded by axes, legends, etc.
# set the core width and height to a fixed dimension
# align the individual 3x3 gtables using rbind/cbind
library(egg)
g1 = ggplotGrob(p1)
g2 = ggplotGrob(p2)
g3 = ggplotGrob(p3)

fg1 = gtable_frame(g1, debug = TRUE,  height=unit(1,'null'))
fg2 = gtable_frame(g2, debug = TRUE,  height=unit(1,'null'))
fg3 = gtable_frame(g3, debug = TRUE,  height=unit(1,'null'))

png("figures/article_one/BrokenStick_propvar.png", width = 2000, height = 800)
egg::ggarrange(p1, p2, p3, ncol = 1, nrow = 3, heights = c(10,1,1),
          labels = c("(a)", "(b)", "(c)"))
dev.off()

#############################################################################
#     PHYLOGENY - STABILITY AMONG PHYLOGENY
#############################################################################
#============================================================================
# Figures. Phylogeny annotated with summary statistics
#============================================================================




#############################################################################
#     GENE COUNT
#############################################################################
#============================================================================
# Rec rate ~ gene count, for each species, 2X2 plot
#============================================================================
meanrecrate_genecount = read.table("output/gene_count/meanrecrate_genecount.txt", header = TRUE)
meanrecrate_genecount$species = gsub("_", " ", meanrecrate_genecount$species)
# FIGURES
# Divide the figure in 4 plots for visualization (44 species)
# Divide in four subsets
d1 = meanrecrate_genecount[which(meanrecrate_genecount$species %in% unique(meanrecrate_genecount$species)[1:10]),]
d2 = meanrecrate_genecount[which(meanrecrate_genecount$species %in% unique(meanrecrate_genecount$species)[11:22]),]
d3 = meanrecrate_genecount[which(meanrecrate_genecount$species %in% unique(meanrecrate_genecount$species)[23:33]),]
d4 = meanrecrate_genecount[which(meanrecrate_genecount$species %in% unique(meanrecrate_genecount$species)[34:43]),]
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
textsize = 12
p1 = ggplot(data = d1, aes(x = gene_count, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(-3, 3)) +
  labs(x = "Gene count", y = "Standardized recombination rate", colour = "Species") +
  guides(color = guide_legend(override.aes = list(linetype = 0))) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=textsize, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=textsize,hjust = 0.5),
        axis.title.x = element_text(color="black", size=textsize),
        axis.title.y = element_text(color="black", size=textsize),
        axis.text=element_text(size=textsize, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=textsize),
        legend.title=element_text(size=textsize))
p2 = ggplot(data = d2, aes(x = gene_count, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(-3, 3)) +
  labs(x = "Gene count", y = "Standardized recombination rate", colour = "Species") +
  guides(color = guide_legend(override.aes = list(linetype = 0))) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=textsize, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=textsize,hjust = 0.5),
        axis.title.x = element_text(color="black", size=textsize),
        axis.title.y = element_text(color="black", size=textsize),
        axis.text=element_text(size=textsize, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=textsize),
        legend.title=element_text(size=textsize))
p3 = ggplot(data = d3, aes(x = gene_count, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(-3, 3)) +
  labs(x = "Gene count", y = "Standardized recombination rate", colour = "Species") +
  guides(color = guide_legend(override.aes = list(linetype = 0))) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=textsize, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=textsize,hjust = 0.5),
        axis.title.x = element_text(color="black", size=textsize),
        axis.title.y = element_text(color="black", size=textsize),
        axis.text=element_text(size=textsize, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=textsize),
        legend.title=element_text(size=textsize))
p4 = ggplot(data = d4, aes(x = gene_count, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(-3, 3)) +
  labs(x = "Gene count", y = "Standardized recombination rate", colour = "Species") +
  guides(color = guide_legend(override.aes = list(linetype = 0))) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=textsize, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=textsize,hjust = 0.5),
        axis.title.x = element_text(color="black", size=textsize),
        axis.title.y = element_text(color="black", size=textsize),
        axis.text=element_text(size=textsize, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=textsize),
        legend.title=element_text(size=textsize))
p = ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
ggsave("figures/article_one/meanrecrate~genecount_allspecies.png",
       plot = p, device = "png", units = "cm", width = 40, height = 20)

#============================================================================
# Rec rate ~ gene count, all species pooled (dark grey transparent dots) + quadratic phylogenetic regression with C.I.
#============================================================================
meanrecrate_genecount = read.table("output/gene_count/meanrecrate_genecount.txt", header = TRUE)
meanrecrate_genecount$species = gsub("_", " ", meanrecrate_genecount$species)
# Subset only for a gene count <= 20
meanrecrate_genecount = meanrecrate_genecount[which(meanrecrate_genecount$gene_count < 21),]
# Quadratic Phylogenetic Regression
# Load the phylogenetic tree estimated in 'Phylogeny.R'
load("output/phylogeny/tree.Rdata")
# Drop unsampled tips
library(ape)
tree = drop.tip(tree, tip = tree$tip.label[!(tree$tip.label %in% meanrecrate_genecount$species)])
# PGLMM
# Bayesian estimates for confidence interval (Bayesian INLA)
pglmm.model = pglmm(mean_rec ~ gene_count + I(gene_count^2) + (1|species__), data = meanrecrate_genecount, family = 'gaussian',
                    cov_ranef = list(species = tree), bayes = TRUE)
# Quadratic function and C.I.
quadrafun = function(x, c) fixef(pglmm.model)[[1]][1] + fixef(pglmm.model)[[1]][2]*x + fixef(pglmm.model)[[1]][3]*x^2
lowerquadrafun = function(x, c) fixef(pglmm.model)[[2]][1] + fixef(pglmm.model)[[2]][2]*x + fixef(pglmm.model)[[2]][3]*x^2
upperquadrafun = function(x, c) fixef(pglmm.model)[[3]][1] + fixef(pglmm.model)[[3]][2]*x + fixef(pglmm.model)[[3]][3]*x^2
# The C.I. values
df.new = data.frame(gene_count = 0:20)
df.new$lwr.pred = lowerquadrafun(df.new$gene_count)
df.new$upr.pred = upperquadrafun(df.new$gene_count)

# FIGURES
textsize = 16

p = ggplot(data = meanrecrate_genecount, aes(x = gene_count, y = mean_rec)) +
  # geom_point(colour = "Darkgrey", alpha = 0.5) +
  geom_point(aes(colour = species, alpha = 0.5)) +
  scale_x_continuous(limits = c(0, 20)) +
  # scale_y_continuous(limits = c(0, 6)) +
  labs(x = "Gene count", y = "Standardized recombination rate") +
  # geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
  stat_function(fun = quadrafun, colour = "Black", size = 1) +
  geom_ribbon(data = df.new, aes(x = gene_count, ymin = lwr.pred, ymax = upr.pred), alpha = 0.3, inherit.aes = F, fill = "Black") +
  # stat_smooth(aes(y = mean_rec), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Black", fill = "Black", alpha = 0.4) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=textsize, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=textsize,hjust = 0.5),
        axis.title.x = element_text(color="black", size=textsize),
        axis.title.y = element_text(color="black", size=textsize),
        axis.text=element_text(size=textsize, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "none")
p
ggsave("figures/article_one/meanrecrate~genecount_quadraticregression.png",
       plot = p, device = "png", units = "cm", width = 20, height = 15)

#============================================================================
# Rec rate ~ gene count, all species pooled (dark grey transparent dots), + bootstrapped mean ± C.I. (pooled species)
#============================================================================
# Compute mean and C.I. for each gene count
df = data.frame(gene_count = c(0:20), mean_rec = NA, lower_rec = NA, upper_rec = NA)

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

# FIGURES
p = ggplot(data = meanrecrate_genecount, aes(x = gene_count, y = mean_rec)) +
  geom_point(colour = "Darkgrey", alpha = 0.5) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0, 6)) +
  labs(x = "Gene count", y = "Recombination rate (cM/Mb, scaled)") +
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
ggsave("figures/article_one/meanrecrate~genecount_bootstrap.png",
       plot = p, device = "png", units = "cm", width = 20, height = 15)

#----------------------------------------------------------------------------
# DIAGNOSTIC PLOTS
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Mean recombination rate ~ Chromosome size
#----------------------------------------------------------------------------
mod = lm(mean_rec ~ gene_count + I(gene_count^2), data = meanrecrate_genecount)
# mod = glm(mean_rec ~ gene_count + I(gene_count^2), data = meanrecrate_genecount, family = "quasipoisson")
# mod = glm(mean_rec ~ as.factor(gene_count), data = meanrecrate_genecount, family = "quasipoisson")
summary(mod)
# Diagnostic plots
png(filename = "figures/article_one/diagnosticplots_lm_meanrecrate~genecount.png", units = "cm", width = 20, height = 8, res = 150)
par(mfrow = c(1,3))
# Residuals distribution
qqPlot(residuals(mod), xlab = "Theoretical Normal quantiles", ylab = "Residuals")
title("(a)", adj = 0)
# Residuals as a function of fitted
plot(fitted(mod), residuals(mod), xlab = "Fitted", ylab = "Residuals")
abline(h = 0, col = "red")
title("(b)", adj = 0)
# Observed~predicted qqplot
# plot(mod$model$mean_rec, mod$coefficients[1]  + mod$coefficients[2]*mod$model$gene_count + mod$coefficients[3]*(mod$model$gene_count)^2, xlab = "Observed", ylab = "Predicted")
plot(mod$model$mean_rec, predict(mod), xlab = "Observed", ylab = "Predicted")
abline(a = 0, b = 1, col = "red")
title("(c)", adj = 0)
par(mfrow = c(1,1))
dev.off()

#----------------------------------------------------------------------------
# PGLMM Mean recombination rate ~ Chromosome size
#----------------------------------------------------------------------------
pglmm.model = pglmm(mean_rec ~ gene_count + I(gene_count^2) + (1|species__), data = meanrecrate_genecount, family = 'gaussian',
                    cov_ranef = list(species = tree))
# Diagnostic plots
png(filename = "figures/article_one/diagnosticplots_pglmm_meanrecrate~genecount.png", units = "cm", width = 20, height = 8, res = 150)
par(mfrow = c(1,3))
# Residuals distribution
qqPlot(residuals(pglmm.model), xlab = "Theoretical Normal quantiles", ylab = "Residuals")
title("(a)", adj = 0)
# Residuals as a function of fitted
plot(fitted(pglmm.model), residuals(pglmm.model), xlab = "Fitted", ylab = "Residuals")
abline(h = 0, col = "red")
title("(b)", adj = 0)
# Observed~predicted qqplot
plot(pglmm.model$Y, unlist(pglmm.predicted.values(pglmm.model)), xlab = "Observed", ylab = "Predicted")
abline(a = 0, b = 1, col = "red")
title("(c)", adj = 0)
par(mfrow = c(1,1))
dev.off()



#############################################################################
#     GENE DENSITY
#############################################################################
#============================================================================
# Rec rate ~ gene density, for each species, 2X2 plot
#============================================================================
meanrecrate_genedensity = read.table("output/gene_density/meanrecrate_genedensity.txt", header = TRUE)
meanrecrate_genedensity$species = gsub("_", " ", meanrecrate_genedensity$species)
unique(meanrecrate_genedensity$species)
meanrecrate_genedensity$gene_density = rowMeans(cbind(meanrecrate_genedensity$quantile.lower, meanrecrate_genedensity$quantile.upper), na.rm = TRUE)
# FIGURES
# Divide the figure in 4 plots for visualization (44 species)
# Divide in four subsets
d1 = meanrecrate_genedensity[which(meanrecrate_genedensity$species %in% unique(meanrecrate_genedensity$species)[1:10]),]
d2 = meanrecrate_genedensity[which(meanrecrate_genedensity$species %in% unique(meanrecrate_genedensity$species)[11:22]),]
d3 = meanrecrate_genedensity[which(meanrecrate_genedensity$species %in% unique(meanrecrate_genedensity$species)[23:33]),]
d4 = meanrecrate_genedensity[which(meanrecrate_genedensity$species %in% unique(meanrecrate_genedensity$species)[34:43]),]

textsize = 12
p1 = ggplot(data = d1, aes(x = gene_density, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(-3, 3)) +
  labs(x = "Gene density", y = "Standardized recombination rate", colour = "Species") +
  guides(color = guide_legend(override.aes = list(linetype = 0))) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=textsize, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=textsize,hjust = 0.5),
        axis.title.x = element_text(color="black", size=textsize),
        axis.title.y = element_text(color="black", size=textsize),
        axis.text=element_text(size=textsize, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=textsize),
        legend.title=element_text(size=textsize))
p2 = ggplot(data = d2, aes(x = gene_density, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(-3, 3)) +
  labs(x = "Gene density", y = "Standardized recombination rate", colour = "Species") +
  guides(color = guide_legend(override.aes = list(linetype = 0))) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=textsize, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=textsize,hjust = 0.5),
        axis.title.x = element_text(color="black", size=textsize),
        axis.title.y = element_text(color="black", size=textsize),
        axis.text=element_text(size=textsize, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=textsize),
        legend.title=element_text(size=textsize))
p3 = ggplot(data = d3, aes(x = gene_density, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(-3, 3)) +
  labs(x = "Gene density", y = "Standardized recombination rate", colour = "Species") +
  guides(color = guide_legend(override.aes = list(linetype = 0))) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=textsize, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=textsize,hjust = 0.5),
        axis.title.x = element_text(color="black", size=textsize),
        axis.title.y = element_text(color="black", size=textsize),
        axis.text=element_text(size=textsize, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=textsize),
        legend.title=element_text(size=textsize))
p4 = ggplot(data = d4, aes(x = gene_density, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(-3, 3)) +
  labs(x = "Gene density", y = "Standardized recombination rate", colour = "Species") +
  guides(color = guide_legend(override.aes = list(linetype = 0))) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=textsize, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=textsize,hjust = 0.5),
        axis.title.x = element_text(color="black", size=textsize),
        axis.title.y = element_text(color="black", size=textsize),
        axis.text=element_text(size=textsize, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=textsize),
        legend.title=element_text(size=textsize))
p = ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
ggsave("figures/article_one/meanrecrate~genedensity_allspecies.png",
       plot = p, device = "png", units = "cm", width = 40, height = 20)

#============================================================================
# Rec rate ~ gene count, all species pooled (dark grey transparent dots) + quadratic phylogenetic regression with C.I.
#============================================================================
meanrecrate_genedensity = read.table("output/gene_density/meanrecrate_genedensity.txt", header = TRUE)
meanrecrate_genedensity$species = gsub("_", " ", meanrecrate_genedensity$species)
unique(meanrecrate_genedensity$species)
meanrecrate_genedensity$gene_density = rowMeans(cbind(meanrecrate_genedensity$quantile.lower, meanrecrate_genedensity$quantile.upper), na.rm = TRUE)
# Quadratic Phylogenetic Regression
# Load the phylogenetic tree estimated in 'Phylogeny.R'
load("output/phylogeny/tree.Rdata")
# Drop unsampled tips
library(ape)
tree = drop.tip(tree, tip = tree$tip.label[!(tree$tip.label %in% meanrecrate_genedensity$species)])
# PGLMM
# Bayesian estimates for confidence interval (Bayesian INLA)
pglmm.model = pglmm(mean_rec ~ gene_density + I(gene_density^2) + (1|species__), data = meanrecrate_genedensity, family = 'gaussian',
                    cov_ranef = list(species = tree), bayes = TRUE)
# Quadratic function and C.I.
quadrafun = function(x, c) fixef(pglmm.model)[[1]][1] + fixef(pglmm.model)[[1]][2]*x + fixef(pglmm.model)[[1]][3]*x^2
lowerquadrafun = function(x, c) fixef(pglmm.model)[[2]][1] + fixef(pglmm.model)[[2]][2]*x + fixef(pglmm.model)[[2]][3]*x^2
upperquadrafun = function(x, c) fixef(pglmm.model)[[3]][1] + fixef(pglmm.model)[[3]][2]*x + fixef(pglmm.model)[[3]][3]*x^2
# The C.I. values
nquantiles = 20
df.new = data.frame(gene_density = seq(0, 1, length.out = nquantiles+1))
df.new$lwr.pred = lowerquadrafun(df.new$gene_density)
df.new$upr.pred = upperquadrafun(df.new$gene_density)

# FIGURES
textsize = 16
p = ggplot(data = meanrecrate_genedensity, aes(x = gene_density, y = mean_rec)) +
  geom_point(colour = "Darkgrey", alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1)) +
  # scale_y_continuous(limits = c(0, 6)) +
  labs(x = "Gene density", y = "Standardized recombination rate") +
  # geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
  stat_function(fun = quadrafun, colour = "Black", size = 1) +
  geom_ribbon(data = df.new, aes(x = gene_density, ymin = lwr.pred, ymax = upr.pred), alpha = 0.3, inherit.aes = F, fill = "Black") +
  # stat_smooth(aes(y = mean_rec), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Black", fill = "Black", alpha = 0.4) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=textsize, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=textsize,hjust = 0.5),
        axis.title.x = element_text(color="black", size=textsize),
        axis.title.y = element_text(color="black", size=textsize),
        axis.text=element_text(size=textsize, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90))
p
ggsave("figures/article_one/meanrecrate~genedensity_quadraticregression.png",
       plot = p, device = "png", units = "cm", width = 20, height = 15)

#----------------------------------------------------------------------------
# DIAGNOSTIC PLOTS
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Mean recombination rate ~ Chromosome size
#----------------------------------------------------------------------------
mod = lm(mean_rec ~ gene_density + I(gene_density^2), data = meanrecrate_genedensity)
# mod = glm(mean_rec ~ gene_count + I(gene_count^2), data = meanrecrate_genecount, family = "quasipoisson")
# mod = glm(mean_rec ~ as.factor(gene_count), data = meanrecrate_genecount, family = "quasipoisson")
summary(mod)
# Diagnostic plots
png(filename = "figures/article_one/diagnosticplots_lm_meanrecrate~genedensity.png", units = "cm", width = 20, height = 8, res = 150)
par(mfrow = c(1,3))
# Residuals distribution
qqPlot(residuals(mod), xlab = "Theoretical Normal quantiles", ylab = "Residuals")
title("(a)", adj = 0)
# Residuals as a function of fitted
plot(fitted(mod), residuals(mod), xlab = "Fitted", ylab = "Residuals")
abline(h = 0, col = "red")
title("(b)", adj = 0)
# Observed~predicted qqplot
# plot(mod$model$mean_rec, mod$coefficients[1]  + mod$coefficients[2]*mod$model$gene_count + mod$coefficients[3]*(mod$model$gene_count)^2, xlab = "Observed", ylab = "Predicted")
plot(mod$model$mean_rec, predict(mod), xlab = "Observed", ylab = "Predicted")
abline(a = 0, b = 1, col = "red")
title("(c)", adj = 0)
par(mfrow = c(1,1))
dev.off()

#----------------------------------------------------------------------------
# PGLMM Mean recombination rate ~ Chromosome size
#----------------------------------------------------------------------------
pglmm.model = pglmm(mean_rec ~ gene_density + I(gene_density^2) + (1|species__), data = meanrecrate_genedensity, family = 'gaussian',
                    cov_ranef = list(species = tree))
# Diagnostic plots
png(filename = "figures/article_one/diagnosticplots_pglmm_meanrecrate~genedensity.png", units = "cm", width = 20, height = 8, res = 150)
par(mfrow = c(1,3))
# Residuals distribution
qqPlot(residuals(pglmm.model), xlab = "Theoretical Normal quantiles", ylab = "Residuals")
title("(a)", adj = 0)
# Residuals as a function of fitted
plot(fitted(pglmm.model), residuals(pglmm.model), xlab = "Fitted", ylab = "Residuals")
abline(h = 0, col = "red")
title("(b)", adj = 0)
# Observed~predicted qqplot
plot(pglmm.model$Y, unlist(pglmm.predicted.values(pglmm.model)), xlab = "Observed", ylab = "Predicted")
abline(a = 0, b = 1, col = "red")
title("(c)", adj = 0)
par(mfrow = c(1,1))
dev.off()


#############################################################################
#     GENE DENSITY
#############################################################################
#============================================================================
# Rec rate ~ gene density, for each species, 2X2 plot
#============================================================================
meanrecrate_genedensity = read.table("output/gene_density/meanrecrate_genedensity_divbygenecount.txt", header = TRUE)
meanrecrate_genedensity$species = gsub("_", " ", meanrecrate_genedensity$species)
unique(meanrecrate_genedensity$species)
meanrecrate_genedensity$gene_density = rowMeans(cbind(meanrecrate_genedensity$quantile.lower, meanrecrate_genedensity$quantile.upper), na.rm = TRUE)
# FIGURES
# Divide the figure in 4 plots for visualization (44 species)
# Divide in four subsets
d1 = meanrecrate_genedensity[which(meanrecrate_genedensity$species %in% unique(meanrecrate_genedensity$species)[1:10]),]
d2 = meanrecrate_genedensity[which(meanrecrate_genedensity$species %in% unique(meanrecrate_genedensity$species)[11:22]),]
d3 = meanrecrate_genedensity[which(meanrecrate_genedensity$species %in% unique(meanrecrate_genedensity$species)[23:33]),]
d4 = meanrecrate_genedensity[which(meanrecrate_genedensity$species %in% unique(meanrecrate_genedensity$species)[34:43]),]

textsize = 12
p1 = ggplot(data = d1, aes(x = gene_density, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(-3, 3)) +
  labs(x = "Gene density", y = "Standardized recombination rate", colour = "Species") +
  guides(color = guide_legend(override.aes = list(linetype = 0))) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=textsize, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=textsize,hjust = 0.5),
        axis.title.x = element_text(color="black", size=textsize),
        axis.title.y = element_text(color="black", size=textsize),
        axis.text=element_text(size=textsize, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=textsize),
        legend.title=element_text(size=textsize))
p2 = ggplot(data = d2, aes(x = gene_density, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(-3, 3)) +
  labs(x = "Gene density", y = "Standardized recombination rate", colour = "Species") +
  guides(color = guide_legend(override.aes = list(linetype = 0))) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=textsize, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=textsize,hjust = 0.5),
        axis.title.x = element_text(color="black", size=textsize),
        axis.title.y = element_text(color="black", size=textsize),
        axis.text=element_text(size=textsize, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=textsize),
        legend.title=element_text(size=textsize))
p3 = ggplot(data = d3, aes(x = gene_density, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(-3, 3)) +
  labs(x = "Gene density", y = "Standardized recombination rate", colour = "Species") +
  guides(color = guide_legend(override.aes = list(linetype = 0))) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=textsize, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=textsize,hjust = 0.5),
        axis.title.x = element_text(color="black", size=textsize),
        axis.title.y = element_text(color="black", size=textsize),
        axis.text=element_text(size=textsize, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=textsize),
        legend.title=element_text(size=textsize))
p4 = ggplot(data = d4, aes(x = gene_density, y = mean_rec, colour = species)) +
  geom_line() +
  geom_point() +
  geom_pointrange(aes(ymin = lower_rec, ymax = upper_rec)) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(-3, 3)) +
  labs(x = "Gene density", y = "Standardized recombination rate", colour = "Species") +
  guides(color = guide_legend(override.aes = list(linetype = 0))) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=textsize, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=textsize,hjust = 0.5),
        axis.title.x = element_text(color="black", size=textsize),
        axis.title.y = element_text(color="black", size=textsize),
        axis.text=element_text(size=textsize, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=textsize),
        legend.title=element_text(size=textsize))
p = ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
ggsave("figures/article_one/meanrecrate~genedensity_divbygenecount_allspecies.png",
       plot = p, device = "png", units = "cm", width = 40, height = 20)

#============================================================================
# Rec rate ~ gene count, all species pooled (dark grey transparent dots) + quadratic phylogenetic regression with C.I.
#============================================================================
meanrecrate_genedensity = read.table("output/gene_density/meanrecrate_genedensity_divbygenecount.txt", header = TRUE)
meanrecrate_genedensity$species = gsub("_", " ", meanrecrate_genedensity$species)
unique(meanrecrate_genedensity$species)
meanrecrate_genedensity$gene_density = rowMeans(cbind(meanrecrate_genedensity$quantile.lower, meanrecrate_genedensity$quantile.upper), na.rm = TRUE)
# Quadratic Phylogenetic Regression
# Load the phylogenetic tree estimated in 'Phylogeny.R'
load("output/phylogeny/tree.Rdata")
# Drop unsampled tips
library(ape)
tree = drop.tip(tree, tip = tree$tip.label[!(tree$tip.label %in% meanrecrate_genedensity$species)])
# PGLMM
# Bayesian estimates for confidence interval (Bayesian INLA)
pglmm.model2 = pglmm(mean_rec ~ gene_density + I(gene_density^2) + (1|species__), data = meanrecrate_genedensity, family = 'gaussian',
                    cov_ranef = list(species = tree), bayes = TRUE)
# Quadratic function and C.I.
quadrafun2 = function(x, c) fixef(pglmm.model2)[[1]][1] + fixef(pglmm.model2)[[1]][2]*x + fixef(pglmm.model2)[[1]][3]*x^2
lowerquadrafun2 = function(x, c) fixef(pglmm.model2)[[2]][1] + fixef(pglmm.model2)[[2]][2]*x + fixef(pglmm.model2)[[2]][3]*x^2
upperquadrafun2 = function(x, c) fixef(pglmm.model2)[[3]][1] + fixef(pglmm.model2)[[3]][2]*x + fixef(pglmm.model2)[[3]][3]*x^2
# The C.I. values
nquantiles = 20
df.new = data.frame(gene_density = seq(0, 1, length.out = nquantiles+1))
df.new$lwr.pred = lowerquadrafun2(df.new$gene_density)
df.new$upr.pred = upperquadrafun2(df.new$gene_density)

# FIGURES
textsize = 16
p_prime = ggplot(data = meanrecrate_genedensity, aes(x = gene_density, y = mean_rec)) +
  geom_point(colour = "Darkgrey", alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1)) +
  # scale_y_continuous(limits = c(0, 6)) +
  labs(x = "Gene density divided by gene count", y = "Standardized recombination rate") +
  # geom_abline(intercept = fixef(pglmm.model)[[1]][1], slope = fixef(pglmm.model)[[1]][2]) + # The PGLMM regression line
  stat_function(fun = quadrafun2, colour = "Black", size = 1) +
  geom_ribbon(data = df.new, aes(x = gene_density, ymin = lwr.pred, ymax = upr.pred), alpha = 0.3, inherit.aes = F, fill = "Black") +
  # stat_smooth(aes(y = mean_rec), method = "lm", formula = y ~ x + I(x^2), size = 1, colour = "Black", fill = "Black", alpha = 0.4) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=textsize, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=textsize,hjust = 0.5),
        axis.title.x = element_text(color="black", size=textsize),
        axis.title.y = element_text(color="black", size=textsize),
        axis.text=element_text(size=textsize, colour="black"),
        strip.text=element_text(size=12, colour="black", angle = 90))
p_prime
ggsave("figures/article_one/meanrecrate~genedensity_divbygenecount_quadraticregression.png",
       plot = p_prime, device = "png", units = "cm", width = 20, height = 15)

#----------------------------------------------------------------------------
# DIAGNOSTIC PLOTS
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# Mean recombination rate ~ Chromosome size
#----------------------------------------------------------------------------
mod = lm(mean_rec ~ gene_density + I(gene_density^2), data = meanrecrate_genedensity)
# mod = glm(mean_rec ~ gene_count + I(gene_count^2), data = meanrecrate_genecount, family = "quasipoisson")
# mod = glm(mean_rec ~ as.factor(gene_count), data = meanrecrate_genecount, family = "quasipoisson")
summary(mod)
# Diagnostic plots
png(filename = "figures/article_one/diagnosticplots_lm_meanrecrate~genedensity_divbygenecount.png", units = "cm", width = 20, height = 8, res = 150)
par(mfrow = c(1,3))
# Residuals distribution
qqPlot(residuals(mod), xlab = "Theoretical Normal quantiles", ylab = "Residuals")
title("(a)", adj = 0)
# Residuals as a function of fitted
plot(fitted(mod), residuals(mod), xlab = "Fitted", ylab = "Residuals")
abline(h = 0, col = "red")
title("(b)", adj = 0)
# Observed~predicted qqplot
# plot(mod$model$mean_rec, mod$coefficients[1]  + mod$coefficients[2]*mod$model$gene_count + mod$coefficients[3]*(mod$model$gene_count)^2, xlab = "Observed", ylab = "Predicted")
plot(mod$model$mean_rec, predict(mod), xlab = "Observed", ylab = "Predicted")
abline(a = 0, b = 1, col = "red")
title("(c)", adj = 0)
par(mfrow = c(1,1))
dev.off()

#----------------------------------------------------------------------------
# PGLMM Mean recombination rate ~ Chromosome size
#----------------------------------------------------------------------------
pglmm.model = pglmm(mean_rec ~ gene_density + I(gene_density^2) + (1|species__), data = meanrecrate_genedensity, family = 'gaussian',
                    cov_ranef = list(species = tree))
# Diagnostic plots
png(filename = "figures/article_one/diagnosticplots_pglmm_meanrecrate~genedensity_divbygenecount.png", units = "cm", width = 20, height = 8, res = 150)
par(mfrow = c(1,3))
# Residuals distribution
qqPlot(residuals(pglmm.model), xlab = "Theoretical Normal quantiles", ylab = "Residuals")
title("(a)", adj = 0)
# Residuals as a function of fitted
plot(fitted(pglmm.model), residuals(pglmm.model), xlab = "Fitted", ylab = "Residuals")
abline(h = 0, col = "red")
title("(b)", adj = 0)
# Observed~predicted qqplot
plot(pglmm.model$Y, unlist(pglmm.predicted.values(pglmm.model)), xlab = "Observed", ylab = "Predicted")
abline(a = 0, b = 1, col = "red")
title("(c)", adj = 0)
par(mfrow = c(1,1))
dev.off()


#----------------------------------------------------------------------------
# COMAPRE WITH/WITHOUT DIVIDING BY GENE COUNT
#----------------------------------------------------------------------------
p_compare = ggarrange(p, p_prime, ncol = 2, nrow = 1, labels = c("(a)", "(b)"))
ggsave("figures/article_one/meanrecrate~genedensity_quadraticregression_compare.png",
       plot = p_compare, device = "png", units = "cm", width = 30, height = 12)







#############################################################################
# PCA & CLASSIFICATION
#############################################################################
#============================================================================
# PCA - CHROMOSOMES
#============================================================================
load(file = "output/classification/data.pca.chromosomes.Rda")
#####################################
# # PLOTS
# # AXES 1&2
# Color gradient for chromosome size
list_chrsize = chromosomes$phys.map.length
# Plot of individuals 
plot.ind = fviz_pca_ind(pca, axes = c(1, 2), label = "none", # hide individual labels
                        col.ind = log(list_chrsize/1000000),
                        legend.title = "Mb") +
  scale_colour_viridis_c(option = "viridis", breaks = c(log(20), log(50), log(200), log(400)), labels = c(20, 50, 200, 400)) +
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
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        plot.margin = margin(0, 0, 0, 0, "cm"))
# Control variable colors using their contributions
plot.cor = fviz_pca_var(pca, axes = c(1, 2), col.var="contrib",
                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                        repel = TRUE, legend.title = "Contrib.", labelsize = 4) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="white"),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        plot.margin = margin(0, 0, 0, 0, "cm"))

png("figures/article_one/PCA_Chromosomes_Axes1-2.png", units = "cm", width = 30, height = 10, res = 150)
ggarrange(plot.ind, plot.cor, ncol = 2, nrow = 1, widths = c(2,1),
          labels = c("(a)", "(b)"))
dev.off()


# AXES 1&3
list_chrsize = chromosomes$phys.map.length
# Plot of individuals 
plot.ind = fviz_pca_ind(pca, axes = c(1, 3), label = "none", # hide individual labels
                        col.ind = log(list_chrsize/1000000),
                        legend.title = "Mb") +
  scale_colour_viridis_c(option = "viridis", breaks = c(log(20), log(50), log(200), log(400)), labels = c(20, 50, 200, 400)) +
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
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        plot.margin = margin(0, 0, 0, 0, "cm"))
# Control variable colors using their contributions
plot.cor = fviz_pca_var(pca, axes = c(1, 3), col.var="contrib",
                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                        repel = TRUE, legend.title = "Contrib.", labelsize = 4) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="white"),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        plot.margin = margin(0, 0, 0, 0, "cm"))
png("figures/article_one/PCA_Chromosomes_Axes1-3.png", units = "cm", width = 30, height = 10, res = 150)
ggarrange(plot.ind, plot.cor, ncol = 2, nrow = 1, widths = c(2,1),
          labels = c("(a)", "(b)"))
dev.off()

#============================================================================
# PCA - POOL CHROMOSOMES AND AVERAGE SUMMARY STATISTICS PER SPECIES
#============================================================================
load(file = "output/classification/data.pca.species.Rda")

# # PLOTS
# # AXES 1&2
# Color gradient for chromosome size
list_sp = levels(chromosomes$species)
list_chrsize = NA
for (i in 1:length(list_sp)) {
  list_chrsize[i] = mean(chromosomes$phys.map.length[chromosomes$species == list_sp[i]], na.rm = TRUE)
}
# Plot of individuals 
plot.ind = fviz_pca_ind(pca, axes = c(1, 2), repel = FALSE, geom.ind = c("point", "text"), col.ind = log(list_chrsize/1000000),
                        legend.title = "Mb", labelsize = 4) +
  scale_colour_viridis_c(breaks = c(log(20), log(50), log(200), log(400)), labels = c(20, 50, 200, 400)) +
  xlim(-6, 6) +
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
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        plot.margin = margin(0, 0, 0, 0, "cm"))
# Control variable colors using their contributions
plot.cor = fviz_pca_var(pca, axes = c(1, 2), col.var="contrib",
                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                        repel = TRUE, legend.title = "Contrib.", labelsize = 4) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="white"),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        plot.margin = margin(0, 0, 0, 0, "cm"))

png("figures/article_one/PCA_Species_Axes1-2.png", units = "cm", width = 30, height = 10, res = 150)
ggarrange(plot.ind, plot.cor, ncol = 2, nrow = 1, widths = c(2,1),
          labels = c("(a)", "(b)"))
dev.off()


# # AXES 1&3
# Plot of individuals 
plot.ind = fviz_pca_ind(pca, axes = c(1, 3), repel = FALSE, geom.ind = c("point", "text"), col.ind = log(list_chrsize/1000000),
                        legend.title = "Mb", labelsize = 4) +
  scale_colour_viridis_c(breaks = c(log(20), log(50), log(200), log(400)), labels = c(20, 50, 200, 400)) +
  xlim(-6, 6) +
  ylim(-2, 2.5) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=18, colour="black"),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        plot.margin = margin(0, 0, 0, 0, "cm"))
# Control variable colors using their contributions
plot.cor = fviz_pca_var(pca, axes = c(1, 3), col.var="contrib",
                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                        repel = TRUE, legend.title = "Contrib.", labelsize = 4) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="white"),
        legend.position = "right",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(0.8,"cm"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        plot.margin = margin(0, 0, 0, 0, "cm"))

png("figures/article_one/PCA_Species_Axes1-3.png", units = "cm", width = 30, height = 10, res = 150)
ggarrange(plot.ind, plot.cor, ncol = 2, nrow = 1, widths = c(2,1),
          labels = c("(a)", "(b)"))
dev.off()

#============================================================================
# HIERARCHICAL CLASSIFICATION ON SPECIES
#============================================================================
png(filename = "figures/article_one/HierarchicalClustering_Species.png", units = "cm", width = 20, height = 10,
    res = 150)
load(file = "output/classification/h1.clust.Rda")
hclustering = as.dendrogram(h1)
list_sp = labels(hclustering)
list_chrsize = NA
for (i in 1:length(list_sp)) {
  list_chrsize[i] = mean(chromosomes$phys.map.length[chromosomes$species == list_sp[i]], na.rm = TRUE)
}
color.gradient = function(x, colors = c("lightgrey", "black"), colsteps = 100) {
  return(colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}
list_color = color.gradient(log(list_chrsize))
list_color = color.gradient(log(list_chrsize), colors = viridis(3))
# ggdendrogram(hclustering, rotate = FALSE, size = 2, tip.color = list_chrsize) + geom_text(colour = list_color)
library(dendextend)
par(mar=c(12,3,3,1))
dend = as.dendrogram(h1)
dend %>% dendextend::set("labels_col", value = list_color) %>% plot
dend %>% rect.dendrogram(k=6, border = "Darkgrey", lty = 5, lwd = 2)
# dend %>% rect.dendrogram(k=5, lty = 5, lwd = 0, which = c(1, 3, 5), col = rgb(0.1, 0.2, 0.4, 0.1) ) 
par(mar =  c(5.1, 4.1, 4.1, 2.1))
dev.off()

#############################################################################









































#############################################################################
#     DEPRECATED
#############################################################################
#============================================================================
# Figures. Broken stick barplot solo
#============================================================================
# #============================================================================
# # Figures. Vertical phylogeny with Broken stick barplot annotated with other summary statistics
# #============================================================================
# # ordered by the phylogeny
# load("output/phylogeny/tree_reduced.Rdata")
# tree$tip.label
# load(file = "output/brokenstick/brokenstick10.Rda")
# 
# # It seems the key is to look at the edge property. The tips are always the first nodes to be given an ID, which will simply correspond to the position in the tip.label vector.
# # First step is to filter out internal nodes from the the second column of the edge matrix:
# is_tip <- tree$edge[,2] <= length(tree$tip.label)
# ordered_tips <- tree$edge[is_tip, 2]
# # Then you can use this vector to extract the tips in the right order:
# tree$tip.label[ordered_tips]
# p = ggtree(tree)
# p
# 
# # Reorder by phylogenetic order
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
# 
# # Plotting the barplot
# p1 = ggplot(data = brokenstick, aes(x = sample, y = proportion.length, fill = color))+
#   geom_bar(stat='identity', width = 1) +
#   # scale_fill_manual(values = color) +
#   scale_fill_viridis_c() +
#   # scale_fill_gradient2() +
#   facet_grid(~set, scales = "free", space="fixed") +
#   labs(x="", y="Proportion of\ntotal physical length", fill="Segment") +
#   theme(axis.line = element_blank(),
#         # axis.line.x = element_blank(), # No x axis
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
#         plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
#         axis.title.x = element_text(color="black", size=16),
#         axis.title.y = element_text(color="black", size=16),
#         axis.text=element_text(size=16, colour="black"),
#         axis.text.x=element_blank(), # No samples names
#         axis.ticks.x=element_blank(), # No x axis
#         strip.text=element_text(size=18, colour="black", angle = 90),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.width=unit(0.8,"cm"),
#         legend.text=element_text(size=16),
#         legend.title=element_text(size=16),
#         legend.position = "right")
# p1
# #----------------------------------------------------------------------------
# # Adding covariates
# #----------------------------------------------------------------------------
# covariates = subset(brokenstick, brokenstick$segment == "p1")
# covariates$proportion.length = rep(1, nrow(covariates))
# covariates = covariates[,-c(4, 6)]
# covariates$chr_size = NA
# chromosome.stats$set = as.character(chromosome.stats$set)
# chromosome.stats$chromosome = as.character(chromosome.stats$chromosome)
# for (i in 1:nrow(covariates)) {
#   covariates$chr_size[i] = chromosome.stats$phys.map.length[chromosome.stats$set == gsub("_[0-9A-Za-z]*$", "",covariates$sample[i]) & chromosome.stats$chromosome == covariates$chromosome[i]]/1000000
# }
# covariates$avgrecrate = NA
# for (i in 1:nrow(covariates)) {
#   covariates$avgrecrate[i] = chromosome.stats$mean.recrate[chromosome.stats$set == gsub("_[0-9A-Za-z]*$", "",covariates$sample[i]) & chromosome.stats$chromosome == covariates$chromosome[i]]
# }
# # rbarintra = read.table("tables/df_geneticshuffling_physical.csv", header = TRUE)
# # brokenstick$rbarintra = NA
# # for (i in 1:nrow(brokenstick)) {
# #   brokenstick$rbarintra[i] = rbarintra$geneticshuffling[rbarintra$set == gsub("_[0-9A-Za-z]*$", "",brockenstick$sample[i])]
# # }
# my_breaks = c(20, 50, 500)
# p2 = ggplot(data = covariates, aes(x = sample, y = proportion.length, fill = chr_size)) +
#   geom_bar(stat='identity', width = 1) +
#   # scale_fill_manual(values = color) +
#   scale_fill_viridis_c(breaks = my_breaks, labels = my_breaks, trans = "log", option = "inferno", direction = -1) +
#   # scale_fill_gradient2() +
#   facet_grid(~set, scales = "free", space="fixed") +
#   labs(x = "", y = "\n", fill = "Chromosome size (Mb)") +
#   theme(axis.line = element_blank(),
#         # axis.line.x = element_blank(), # No x axis
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
#         plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
#         axis.title.x = element_text(color="black", size=16),
#         axis.title.y = element_text(color="black", size=16),
#         axis.text=element_text(size=16, colour="white"),
#         axis.text.x=element_blank(), # No samples names
#         axis.ticks=element_blank(), # No axis ticks
#         strip.text=element_blank(),
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.width=unit(0.8,"cm"),
#         legend.text=element_text(size=16),
#         legend.title=element_text(size=16),
#         legend.position = "right",
#         plot.margin = margin(0, 0, 0, 0, "cm"))
# # p2
# my_breaks = c(0, 1, 2, 5, 10, 15)
# p3 = ggplot(data = covariates, aes(x = sample, y = proportion.length, fill = avgrecrate)) +
#   geom_bar(stat='identity', width = 1) +
#   scale_fill_viridis_c(breaks = my_breaks, labels = my_breaks, trans = "log", option = "plasma") +
#   facet_grid(~set, scales = "free", space="fixed") +
#   labs(x="Chromosome", y="\n", fill = "Mean recombination rate (cM/Mb)") +
#   theme(axis.line = element_blank(),
#         # axis.line.x = element_blank(), # No x axis
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
#         plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
#         axis.title.x = element_text(color="black", size=16),
#         axis.title.y = element_text(color="black", size=16),
#         axis.text=element_text(size=16, colour="white"),
#         axis.text.x=element_blank(), # No samples names
#         axis.ticks=element_blank(), # No axis ticks
#         strip.text=element_blank(),
#         legend.position = "right",
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.width=unit(0.8,"cm"),
#         legend.text=element_text(size=16),
#         legend.title=element_text(size=16),
#         plot.margin = margin(0, 0, 0, 0, "cm"))
# # p3
# 
# # Aligning complex ggplots
# # gridExtra::grid.arrange provides no way to align the panels of individual plots. While this is achievable with low-level gtable functions, it often requires substantial effort on a case-by-case basis. The egg package introduces a general strategy for such layout manipulations, with the following steps:
# # decompose each plot into a 3x3 layout, where the central cell corresponds to the core panels, surrounded by axes, legends, etc.
# # set the core width and height to a fixed dimension
# # align the individual 3x3 gtables using rbind/cbind
# library(egg)
# g1 = ggplotGrob(p1)
# g2 = ggplotGrob(p2)
# g3 = ggplotGrob(p3)
# 
# fg1 = gtable_frame(g1, debug = TRUE,  height=unit(6,'null'))
# fg2 = gtable_frame(g2, debug = TRUE,  height=unit(1,'null'))
# fg3 = gtable_frame(g3, debug = TRUE,  height=unit(1,'null'))
# 
# png("figures/article_one/BrokenStick_phylogeny.png", width = 2000, height = 800)
# ggarrange(p1, p2, p3, ncol = 1, nrow = 3, heights = c(5,1,1),
#           labels = c("(a)", "(b)", "(c)"))
# dev.off()
# 
# 
# 
# 
# #----------------------------------------------------------------------------
# #       Annotating with a vertical phylogenetic tree
# #----------------------------------------------------------------------------
# # y position of phylogeny is given by the center of each species facet
# library(ggtree)
# 
# # p = ggtree(tree, layout = 'fan', branch.length = "none", open.angle = 180) +
# #   geom_tiplab(size = 3) +
# #   xlim(-40, 40) +
# #   theme(legend.position = "none")
# # p = rotate_tree(p, -90)
# # p
# p = ggtree(tree, layout = 'rectangular', branch.length = "none") +
#   geom_tiplab(size = 3) +
#   xlim(-40, 40) +
#   theme(legend.position = "none")
# p
# 
# 
# 
# # dendroPlot = ggplot() +
# #   geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) +
# #   labs(x="", y="\n", fill="") +
# #   theme(axis.line=element_blank(),
# #         axis.ticks=element_blank(),
# #         axis.text=element_blank(),
# #         axis.title=element_blank(),
# #         panel.background=element_rect(fill="white"),
# #         panel.grid=element_blank())
# 
# #----------------------------------------------------------------------------
# # Adding covariates
# #----------------------------------------------------------------------------
# # Chromosome size, mean recombination rate, periphery-bias ratio, genetic shuffling rate
# library(ggstance)
# p = facet_plot(p, panel = 'Brokenstick', data = brokenstick, 
#                  geom = geom_barh,
#                  mapping = aes(x = sample, y = set, fill = as.factor(proportionlength)), 
#                  stat='identity' ) 
# 
# #----------------------------------------------------------------------------
# # Adding brokenstick layer
# #----------------------------------------------------------------------------
# p
# 
# 
# 
# 
# 
# #----------------------------------------------------------------------------
# # Saving...
# #----------------------------------------------------------------------------
# ggsave("figures/article_one/BrokenStick_phylogeny.png",
#        device="png",units="cm", width=60, height=20)
# 
# #============================================================================
# 
# 
# 
# 
# 
# brockenstick = read.table("tables/brockenstick/brockenstick_k10.txt", header = TRUE, sep = "\t")
# brockenstick = brockenstick[!(is.na(brockenstick$p1) | is.na(brockenstick$p2) | is.na(brockenstick$p3) | is.na(brockenstick$p4) | is.na(brockenstick$p5) | is.na(brockenstick$p6) | is.na(brockenstick$p7) | is.na(brockenstick$p8) | is.na(brockenstick$p9) | is.na(brockenstick$p10)),]
# 
# brockenstick$set = as.character(brockenstick$set)
# brockenstick$chromosome = as.character(brockenstick$chromosome)
# # Add the information about chromosome size in Mb
# # Barplot under the structure plot
# # df_chrsize = data.frame(set = as.character(brockenstick$set), chromosome = brockenstick$chromosome)
# # df_chrsize$chr_size = NA
# # for (i in 1:nrow(df_chrsize)) {
# #   df_chrsize$chr_size[i] = chromosome.stats$phys.map.length[chromosome.stats$set == df_chrsize$set[i] & chromosome.stats$chromosome == df_chrsize$chromosome[i]]/1000000
# # }
# # df_chrsize$set = gsub("_", " ",regmatches(as.character(brockenstick$set), gregexpr("^[A-Z][a-z]+_[a-z]+", as.character(brockenstick$set))))
# # df_chrsize$p = 1
# 
# # Add the information about the average recombination rate in cM/Mb
# # Barplot under the structure plot
# # df_avgrecrate = data.frame(set = as.character(brockenstick$set), chromosome = brockenstick$chromosome)
# # df_avgrecrate$avgrecrate = NA
# # for (i in 1:nrow(df_avgrecrate)) {
# #   df_avgrecrate$avgrecrate[i] = chromosome.stats$mean.recrate[chromosome.stats$set == df_avgrecrate$set[i] & chromosome.stats$chromosome == df_avgrecrate$chromosome[i]]
# # }
# # df_avgrecrate$set = gsub("_", " ",regmatches(as.character(brockenstick$set), gregexpr("^[A-Z][a-z]+_[a-z]+", as.character(brockenstick$set))))
# # df_avgrecrate$p = 1
# 
# 
# # reformat data
# brockenstick = melt(brockenstick)
# brockenstick = cbind(paste(brockenstick$set, "_", brockenstick$chromosome, sep =""), brockenstick)
# 
# # Get species name
# brockenstick$set = gsub("_", " ",regmatches(as.character(brockenstick$set), gregexpr("^[A-Z][a-z]+_[a-z]+", as.character(brockenstick$set))))
# # Format columns in proper way
# colnames(brockenstick) = c("sample", "set","chromosome","segment","proportion.length")
# 
# brockenstick$chr_size = NA
# for (i in 1:nrow(brockenstick)) {
#   brockenstick$chr_size[i] = chromosome.stats$phys.map.length[chromosome.stats$set == gsub("_[0-9A-Za-z]*$", "",brockenstick$sample[i]) & chromosome.stats$chromosome == brockenstick$chromosome[i]]/1000000
# }
# brockenstick$avgrecrate = NA
# for (i in 1:nrow(brockenstick)) {
#   brockenstick$avgrecrate[i] = chromosome.stats$mean.recrate[chromosome.stats$set == gsub("_[0-9A-Za-z]*$", "",brockenstick$sample[i]) & chromosome.stats$chromosome == brockenstick$chromosome[i]]
# }
# rbarintra = read.table("tables/df_geneticshuffling_physical.csv", header = TRUE)
# brockenstick$rbarintra = NA
# for (i in 1:nrow(brockenstick)) {
#   brockenstick$rbarintra[i] = rbarintra$geneticshuffling[rbarintra$set == gsub("_[0-9A-Za-z]*$", "",brockenstick$sample[i])]
# }
# # Set a vector of gradient color
# # Value of color is simply expected - p, the departure from the expected proportion 1/k
# k = 10
# brockenstick$color = (1/k) - brockenstick$proportion.length
# # brockenstick$color[which(brockenstick$color > 0)] = rescale(brockenstick$color[which(brockenstick$color > 0)], to = c(0, 1), from = range(brockenstick$color[which(brockenstick$color > 0)], na.rm = TRUE))
# # brockenstick$color[which(brockenstick$color < 0)] = rescale(brockenstick$color[which(brockenstick$color < 0)], to = c(-1, 0), from = range(brockenstick$color[which(brockenstick$color < 0)], na.rm = TRUE))
# # brockenstick$color = log10(brockenstick$proportion.length/(1/k))
# 
# 
# #-------------------------------------------------
# # Add hierarchical Clustering for ordering species
# #-------------------------------------------------
# # Load 'h1', results of the clustering procedure in 'Classification.R'
# load(file = "output/classification/hclustPooledChromosomes.Rda")
# 
# # h1$hclust$order
# # Manual ordering
# # new_order = c(21,22,23,10,17,12, 29, 16, 33,3, 15, 9, 25, 5, 30, 31, 32,  6, 14, 27,  2, 34, 13, 20, 26, 19,  8, 11, 28)
# # new_order = c(13, 15, 19, 24, 14, 25, 30, 32, 27,  2,  5,  6, 29, 12, 26,  9,  3, 11, 8, 23, 22, 28, 17, 20,  7, 18,  4, 21, 10, 16, 31,  1)
# # h1$hclust = rotate(h1$hclust, new_order)
# 
# dendr = dendro_data(h1, type="rectangle") # convert for ggplot
# clust = cutree(h1, k = 4)                    # find 4 clusters
# clust.df = data.frame(label = names(clust), cluster = factor(clust))
# # dendr[["labels"]] has the labels, merge with clust.df based on label column
# dendr[["labels"]] <- merge(dendr[["labels"]],clust.df, by="label")
# 
# # # # plot the dendrogram; note use of color=cluster in geom_text(...)
# # ggplot() +
# #   geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) +
# #   geom_text(data=label(dendr), aes(x, y, label=label, hjust=0, color=cluster),
# #             size=3) +
# #   coord_flip() + scale_y_reverse(expand=c(0.2, 0)) +
# #   theme(axis.line.y=element_blank(),
# #         axis.ticks.y=element_blank(),
# #         axis.text.y=element_blank(),
# #         axis.title.y=element_blank(),
# #         panel.background=element_rect(fill="white"),
# #         panel.grid=element_blank())
# 
# brockenstick$set = factor(brockenstick$set, levels = gsub("_", " ", as.character(h1$labels[h1$order])))
# #-------------------------------------------------
# # Add classification to cluster similar species
# # Manual clustering
# # [1] "Arabidopsis thaliana"    "Brachypodium distachyon" "Brassica napus"          "Capsicum annuum"
# #  [5] "Cenchrus americanus"     "Citrullus lanatus"       "Cucumis melo"            "Cucumis sativus"
# #  [9] "Glycine max"             "Gossypium raimondii"     "Malus domestica"         "Manihot esculenta"      
# # [13] "Oryza sativa"            "Phaseolus vulgaris"      "Prunus mume"             "Sesamum indicum"
# # [17] "Setaria italica"         "Solanum lycopersicum"    "Solanum tuberosum"       "Sorghum bicolor"
# # [21] "Theobroma cacao"         "Triticum aestivum"       "Vitis vinifera"          "Zea mays"  
# brockenstick$set = factor(brockenstick$set, levels = c("Arabidopsis thaliana", "Brachypodium distachyon", "Brassica napus",
#                                                        "Citrullus lanatus", "Cucumis sativus",
#                                                        "Malus domestica", "Manihot esculenta", "Oryza sativa", "Populus deltoides",
#                                                        "Populus simonii",
#                                                        "Prunus mume", "Sesamum indicum", "Theobroma cacao", "Vitis vinifera", "Cucumis melo",
#                                                        "Setaria italica","Coffea canephora",
#                                                        "Lupinus angustifolius", "Glycine max","Gossypium hirsutum",
#                                                        "Gossypium raimondii", "Phaseolus vulgaris","Solanum lycopersicum",
#                                                        "Solanum tuberosum", "Sorghum bicolor",
#                                                        "Arachis duranensis", "Capsicum annuum", "Cenchrus americanus", "Zea mays",
#                                                        "Hordeum vulgare", "Triticum aestivum", "Triticum dicoccoides"
# ))
# df_chrsize$set = factor(df_chrsize$set, levels = c("Arabidopsis thaliana", "Brachypodium distachyon", "Brassica napus",
#                                                    "Citrullus lanatus", "Cucumis sativus",
#                                                    "Malus domestica", "Manihot esculenta", "Oryza sativa", "Populus deltoides",
#                                                    "Populus simonii",
#                                                    "Prunus mume", "Sesamum indicum", "Theobroma cacao", "Vitis vinifera", "Cucumis melo",
#                                                    "Setaria italica","Coffea canephora",
#                                                    "Lupinus angustifolius", "Glycine max","Gossypium hirsutum",
#                                                    "Gossypium raimondii", "Phaseolus vulgaris","Solanum lycopersicum",
#                                                    "Solanum tuberosum", "Sorghum bicolor",
#                                                    "Arachis duranensis", "Capsicum annuum", "Cenchrus americanus", "Zea mays",
#                                                    "Hordeum vulgare", "Triticum aestivum", "Triticum dicoccoides"
# ))
# # df_avgrecrate$set = factor(df_avgrecrate$set, levels = c("Arabidopsis thaliana", "Brachypodium distachyon", "Brassica napus",
# #                                                        "Citrullus lanatus", "Cucumis sativus",
# #                                                        "Malus domestica", "Manihot esculenta", "Oryza sativa", "Populus deltoides",
# #                                                        "Populus simonii",
# #                                                        "Prunus mume", "Sesamum indicum", "Theobroma cacao", "Vitis vinifera", "Cucumis melo",
# #                                                        "Setaria italica","Coffea canephora",
# #                                                        "Lupinus angustifolius", "Glycine max","Gossypium hirsutum",
# #                                                        "Gossypium raimondii", "Phaseolus vulgaris","Solanum lycopersicum",
# #                                                        "Solanum tuberosum", "Sorghum bicolor",
# #                                                        "Arachis duranensis", "Capsicum annuum", "Cenchrus americanus", "Zea mays",
# #                                                        "Hordeum vulgare", "Triticum aestivum", "Triticum dicoccoides"
# #                                                        ))
# 
# ######### GGPLOT
# # Plotting the barplot
# my_breaks = c(0.1, 0, -0.2, -0.4, -0.6, -0.8)
# BrockenStick_barplot = ggplot(data = brockenstick, aes(x=sample, y=proportion.length, fill = color))+
#   geom_bar(stat='identity', width = 1) +
#   # scale_fill_manual(values = color) +
#   scale_fill_viridis_c(breaks = my_breaks, labels = my_breaks) +
#   # scale_fill_gradient2() +
#   facet_grid(~set, scales = "free", space="fixed") +
#   labs(x="Chromosome", y="Proportion of\ntotal physical length", fill="Segment relative length") +
#   theme(axis.line = element_blank(),
#         # axis.line.x = element_blank(), # No x axis
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
#         plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
#         axis.title.x = element_text(color="black", size=24),
#         axis.title.y = element_text(color="black", size=24),
#         axis.text=element_text(size=24, colour="black"),
#         axis.text.x=element_blank(), # No samples names
#         axis.ticks.x=element_blank(), # No x axis
#         strip.text=element_text(size=24, colour="black", angle = 90),
#         legend.position = "bottom",
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.width=unit(2,"cm"),
#         legend.text=element_text(size=24),
#         legend.title=element_text(size=24))
# BrockenStick_barplot
# 
# 
# # plotting the dendrogram of species ordering and grouping
# # Breaks are plotted as a function of facets coordinates in BrockenStick_barplot
# # Facets of unequal size
# # Position is the mean chromosome number of the facet (center of the group) + position of the first chromosome in the group (in number of chromosomes before)
# # labels = h1$labels[h1$order]
# # facet_center = NA
# # for (i in 1:length(labels)) {
# #   toMatch = labels[1:i]
# #   facet_center[i] = sum((grepl(paste(toMatch,collapse="|"), brockenstick$sample) & brockenstick$segment == "p1"), na.rm = TRUE)
# # }
# # 
# dendroPlot = ggplot() +
#   geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) +
#   labs(x="", y="\n", fill="") +
#   theme(axis.line=element_blank(),
#         axis.ticks=element_blank(),
#         axis.text=element_blank(),
#         axis.title=element_blank(),
#         panel.background=element_rect(fill="white"),
#         panel.grid=element_blank())
# 
# # dendroPlot = ggdendrogram(h1, rotate = FALSE, size = 2) + theme_dendro()
# # dendroPlot
# ggarrange(dendroPlot, BrockenStick_barplot, nrow = 2, heights = c(5,20))
# 
# my_breaks = c(20, 50, 500)
# ChrSize_barplot = ggplot(data = brockenstick, aes(x=sample, y=proportion.length, fill = chr_size))+
#   geom_bar(stat='identity', width = 1) +
#   # scale_fill_manual(values = color) +
#   scale_fill_viridis_c(breaks = my_breaks, labels = my_breaks, trans = "log", option = "inferno", direction = -1) +
#   # scale_fill_gradient2() +
#   facet_grid(~set, scales = "free", space="fixed") +
#   labs(x="", y="\n", fill="Chromosome size (Mb)") +
#   theme(axis.line = element_blank(),
#         # axis.line.x = element_blank(), # No x axis
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
#         plot.subtitle = element_text(color="black",size=24,hjust = 0.5),
#         axis.title.x = element_text(color="black", size=24),
#         axis.title.y = element_text(color="black", size=24),
#         axis.text=element_text(size=24, colour="white"),
#         axis.text.x=element_blank(), # No samples names
#         axis.ticks=element_blank(), # No axis ticks
#         strip.text=element_blank(),
#         legend.position = "bottom",
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.width=unit(2,"cm"),
#         legend.text=element_text(size=24),
#         legend.title=element_text(size=24))
# # ChrSize_barplot
# 
# my_breaks = c(0, 1, 2, 5, 10, 15)
# AvgRecRate_barplot = ggplot(data = brockenstick, aes(x=sample, y=proportion.length, fill = avgrecrate))+
#   geom_bar(stat='identity', width = 1) +
#   scale_fill_viridis_c(breaks = my_breaks, labels = my_breaks, trans = "log", option = "plasma") +
#   facet_grid(~set, scales = "free", space="fixed") +
#   labs(x="", y="\n", fill="Mean recombination rate (cM/Mb)") +
#   theme(axis.line = element_blank(),
#         # axis.line.x = element_blank(), # No x axis
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
#         plot.subtitle = element_text(color="black",size=24,hjust = 0.5),
#         axis.title.x = element_text(color="black", size=24),
#         axis.title.y = element_text(color="black", size=24),
#         axis.text=element_text(size=24, colour="white"),
#         axis.text.x=element_blank(), # No samples names
#         axis.ticks=element_blank(), # No axis ticks
#         strip.text=element_blank(),
#         legend.position = "bottom",
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.width=unit(2,"cm"),
#         legend.text=element_text(size=24),
#         legend.title=element_text(size=24))
# # AvgRecRate_barplot
# 
# my_breaks = c(0, 0.02, 0.04, 0.06, 0.08 ,0.1)
# rbar_barplot = ggplot(data = brockenstick, aes(x=sample, y=proportion.length, fill = rbarintra))+
#   geom_bar(stat='identity', width = 1) +
#   scale_fill_viridis_c(option = "viridis") +
#   facet_grid(~set, scales = "free", space = "fixed") +
#   labs(x="", y="\n", fill="Intra chromosomal genetic shuffling") +
#   theme(axis.line = element_blank(),
#         # axis.line.x = element_blank(), # No x axis
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
#         plot.subtitle = element_text(color="black",size=24,hjust = 0.5),
#         axis.title.x = element_text(color="black", size=24),
#         axis.title.y = element_text(color="black", size=24),
#         axis.text=element_text(size=24, colour="white"),
#         axis.text.x=element_blank(), # No samples names
#         axis.ticks=element_blank(), # No axis ticks
#         strip.text=element_blank(),
#         legend.position = "bottom",
#         legend.key = element_rect(fill = "white", size = 1),
#         legend.key.width=unit(2,"cm"),
#         legend.text=element_text(size=24),
#         legend.title=element_text(size=24))
# # rbar_barplot
# 
# plot = plot_grid(BrockenStick_barplot, ChrSize_barplot, AvgRecRate_barplot, rbar_barplot,
#                  nrow = 4, rel_heights = c(4, 1, 1, 1))
# 
# #
# png("figures/BrockenStick_barplot_clusters.png", width = 2000, height = 1600)
# grid.newpage()
# print(dendroPlot, vp = viewport(x = 0.525, y = 0.9, width = 1.03, height = 0.2))
# print(plot, vp = viewport(x = 0.5, y = 0.4, width = 1, height = 0.8))
# dev.off()
# # 
# 
# # library("cowplot")
# # plot_grid(dendroPlot, BrockenStick_barplot, ChrSize_barplot, AvgRecRate_barplot, rbar_barplot,
# #           align = "v", axis = "r",
# #           nrow = 5, rel_heights = c(1, 4, 1, 1, 1))
# png("figures/BrockenStick_barplot_clusters.png", width = 2000, height = 1600)
# grid.newpage()
# print(plot_grid(BrockenStick_barplot, ChrSize_barplot, AvgRecRate_barplot, rbar_barplot,
#                 nrow = 4, rel_heights = c(4, 1, 1, 1)))
# dev.off()
# # ggarrange(dendroPlot, BrockenStick_barplot, nrow = 2)
# 
# 
# 













































#############################################################################
# 
#############################################################################


#============================================================================
# Figures. Annotated recombination maps
#============================================================================

# Recombination maps annotated with regions of significant deviation from the mean (high/low recombination)
# Estimated by resampling
# Clean directories from previous results
list_dirs = c("figures/recombination_maps/loess/1Mbwind", "figures/recombination_maps/loess/100kbwind", "figures/recombination_maps/loess/bins",
              "figures/recombination_maps/loess/pointwise",
              "figures/recombination_maps/smooth.spline/1Mbwind", "figures/recombination_maps/smooth.spline/100kbwind", "figures/recombination_maps/smooth.spline/bins",
              "figures/recombination_maps/smooth.spline/pointwise")
unlink(list_dirs, recursive = TRUE)
sapply(list_dirs, function(x) dir.create(x, recursive = TRUE))

# A data frame with maps coordinates
# and regions of high/low significance indicated with horizontal lines






#============================================================================
# Figure. Marey plot (all species & chromosomes, color = species)
#============================================================================
# Scaled distances to {0,1} (= position/total length) for both genetic and physical distances
# Load the joint database with all maps
data.final = read.table(file = "data-cleaned/Marey_maps/AllMaps.txt", header = TRUE, sep = "\t")

df = data.frame(species = as.character(data.final$set), map = as.character(data.final$map),
                gen = as.numeric(data.final$gen), phys = as.numeric(data.final$phys))
# Get species name
df$species = as.character(df$species)
df$species = gsub("_", " ", unlist(regmatches(as.character(df$species), gregexpr("^[A-Z][a-z]+_[a-z]+", as.character(df$species)))))

# Normalized distances (divided by total length)
# Compute for each species and chromosome
for (i in 1:nrow(df)) {
  print(i)
  df$gen[i] = df$gen[i]/max(df$gen[which(df$species == df$species[i] & df$map == df$map[i])], na.rm = TRUE)
  df$phys[i] = df$phys[i]/max(df$phys[which(df$species == df$species[i] & df$map == df$map[i])], na.rm = TRUE)
}
# Remove NA data
df = df[!is.na(df$phys),]
# Define groups
df$group = paste(df$species, df$map)

write.table(df, file = "data-cleaned/Marey_plot.scaled.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

df = read.table(file = "data-cleaned/Marey_plot.scaled.txt", header = TRUE, sep = "\t")

# Plot only selected species
# df = df[(df$species %in% c("Arabidopsis thaliana", "Malus domestica", "Triticum aestivum", "Oryza sativa")),]
# df = df[df$species == "Arabidopsis thaliana",]
# df = df[df$species == "Malus domestica",]
# Plot one chromosome per species
df = df[(df$map == 2) | (df$map == "2A") | (df$map == "2B") | (df$map == "A02"),]

Marey_dots = ggplot(data = df, aes(x = phys, y = gen, group = group, colour = species)) +
  geom_point() +
  xlim(0, 1) +
  ylim(0, 1) +
  labs(x="Physical position (Mb)", y="Genetic position (cM)", fill="Species") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=28,hjust = 0.5),
        axis.title.x = element_text(color="black", size=28),
        axis.title.y = element_text(color="black", size=28),
        # axis.text=element_text(size=28, colour="black"),
        axis.text=element_blank(),
        strip.text=element_text(size=28, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28)
  )
Marey_dots
ggsave("figures/Marey_dots.png",
       device="png",dpi=320,units="cm",width=50,height=40)


Marey_plot = ggplot(data = df, aes(x = phys, y = gen, group = group, colour = species)) +
  # Add the regression line with loess method: local regression fitting
  geom_smooth(method = "loess", se = FALSE, na.rm = TRUE, span = 0.2) +
  xlim(0, 1) +
  ylim(0, 1) +
  labs(x="Physical position (Mb)", y="Genetic position (cM)", fill="Species") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=28,hjust = 0.5),
        axis.title.x = element_text(color="black", size=28),
        axis.title.y = element_text(color="black", size=28),
        # axis.text=element_text(size=28, colour="black"),
        axis.text=element_blank(),
        strip.text=element_text(size=28, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28)
  )
Marey_plot
ggsave("figures/Marey_plot.png",
       device="png",dpi=320,units="cm",width=50,height=40)

rm(data.final)
rm(df)

df = read.table(file = "data-cleaned/Marey_plot.scaled.txt", header = TRUE, sep = "\t")
# Plot one chromosome per species
df = df[(df$map == 2) | (df$map == "2A") | (df$map == "2B") | (df$map == "A02"),]
# Compute gamma parameters for each map
par = data.frame(species = rep(NA, length(unique(df$group))), map = rep(NA, length(unique(df$group))), group = unique(df$group), shape1 = rep(NA, length(unique(df$group))), shape2 = rep(NA, length(unique(df$group))))

# Vector of colours
# install.packages("viridis")  # Install
library("viridis")           # Load
colours = viridis(n = length(unique(df$species)))
# Init plot with neutral (abline)
neutral = data.frame(x = c(0,1), y = c(0,1))
Marey_plot = ggplot(data = neutral, aes(x, y)) +
  # stat_function(fun = pgamma, args = list(shape = par$shape, rate = par$rate)) +
  geom_line() +
  xlim(0, 1) +
  ylim(0, 1) +
  # geom_point(aes(group = map)) +
  labs(x="Physical position (Mb)", y="Genetic position (cM)", fill="Species") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=28,hjust = 0.5),
        axis.title.x = element_text(color="black", size=28),
        axis.title.y = element_text(color="black", size=28),
        # axis.text=element_text(size=28, colour="black"),
        axis.text=element_blank(),
        strip.text=element_text(size=28, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28)
  )
# For each map, compute the gamma parameters (shape and rate) estimated with fitdistr
for (i in 1:nrow(par)) {
  print(as.character(par$group[i]))
  # species name
  sp = as.character(unique(df$species[df$group == par$group[i]]))
  par$species[i] = sp
  # chromosome name
  chr = unique(df$map[df$group == par$group[i]])
  par$map[i] = as.character(chr)
    # Estimates
  #Remove 0 and 1 exact values (not estimated with fitdist)
  dat = df[df$species == sp & df$map == chr,]
  tmp = as.numeric(df$gen[df$species == sp & df$map == chr])
  tmp = tmp[!(tmp == 0)]
  tmp = tmp[!(tmp == 1)]
  Marey_plot = Marey_plot + geom_point(data = dat, aes(x = gen, y = phys), colour = colours[i])
  # Fitting
  fit = fitdist(tmp, distr = "beta")
  par$shape1[i] = as.numeric(fit$estimate[1])
  par$shape2[i] = as.numeric(fit$estimate[2])
  # Now calculate points on the cdf
  # cdf = pgamma(x, meanlog = as.numeric(fit$estimate[1]), sdlog = as.numeric(fit$estimate[2]))
  # Shown plotted here
  # lines(x,cdf)
  Marey_plot = Marey_plot + stat_function(fun = pbeta, args = list(shape1 = par$shape1[i], shape2 = par$shape2[i]), colour = colours[i])
}
# Generate points of the ECDF
Marey_plot = Marey_plot + geom_line(data = neutral, aes(x, y), colour = "red", size = 1.2)
Marey_plot
ggsave("figures/Marey_plot-fitdist.png",
       device="png",dpi=320,units="cm",width=40,height=40)


write.table(par, file = "data-cleaned/Marey_maps/Marey_plot.fitparameters.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
par = read.table(file = "data-cleaned/Marey_maps/Marey_plot.fitparameters.txt", header = TRUE, sep = "\t")


df = read.table(file = "data-cleaned/Marey_maps/Marey_plot.scaled.txt", header = TRUE, sep = "\t")
df = df[(df$species %in% c("Arabidopsis thaliana")),]

par = data.frame(species = rep(NA, length(unique(df$group))), map = rep(NA, length(unique(df$group))), group = unique(df$group), shape1 = rep(NA, length(unique(df$group))), shape2 = rep(NA, length(unique(df$group))))

# Init plot
neutral = data.frame(x = c(0,1), y = c(0,1))
Marey_plot = ggplot(data = neutral, aes(x, y)) +
  # stat_function(fun = pgamma, args = list(shape = par$shape, rate = par$rate)) +
  geom_line() +
  xlim(0, 1) +
  ylim(0, 1) +
  # geom_point(aes(group = map)) +
  labs(x="Physical position (Mb)", y="Genetic position (cM)", fill="Species") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=28,hjust = 0.5),
        axis.title.x = element_text(color="black", size=28),
        axis.title.y = element_text(color="black", size=28),
        # axis.text=element_text(size=28, colour="black"),
        axis.text=element_blank(),
        strip.text=element_text(size=28, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28)
  )
# For each map, compute the gamma parameters (shape and rate) estimated with fitdistr
for (i in 1:nrow(par)) {
  print(as.character(par$group[i]))
  # species name
  sp = as.character(unique(df$species[df$group == par$group[i]]))
  par$species[i] = sp
  # chromosome name
  chr = unique(df$map[df$group == par$group[i]])
  par$map[i] = as.character(chr)
  # Estimates
  #Remove 0 and 1 exact values (not estimated with fitdist)
  tmp = as.numeric(df$gen[df$species == sp & df$map == chr])
  tmp = tmp[!(tmp == 0)]
  tmp = tmp[!(tmp == 1)]
  fit = fitdist(tmp, distr = "beta")
  par$shape1[i] = as.numeric(fit$estimate[1])
  par$shape2[i] = as.numeric(fit$estimate[2])
  # Now calculate points on the cdf
  # cdf = pgamma(x, shape = as.numeric(fit$estimate[1]), rate = as.numeric(fit$estimate[2]))
  # Shown plotted here
  # lines(x,cdf)
  Marey_plot = Marey_plot + stat_function(fun = pbeta, args = list(shape1 = par$shape1[i], shape2 = par$shape2[i]), colour = "black")
}
# Generate points of the ECDF
Marey_plot = Marey_plot + geom_line(data = neutral, aes(x, y), colour = "red", size = 1.2)
Marey_plot
ggsave("figures/Marey_plot-fitdist_Arabidopsis.png",
       device="png",dpi=320,units="cm",width=40,height=40)

Marey_plot = ggplot(data = par, aes(x = x, group = group, colour = species)) +
  stat_function(fun = pgamma, args = list(shape = par$shape, rate = par$rate)) +
  xlim(0, 1) +
  ylim(0, 1) +
  # geom_point(aes(group = map)) +
  labs(x="Physical position (Mb)", y="Genetic position (cM)", fill="Species") +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=28,hjust = 0.5),
        axis.title.x = element_text(color="black", size=28),
        axis.title.y = element_text(color="black", size=28),
        # axis.text=element_text(size=28, colour="black"),
        axis.text=element_blank(),
        strip.text=element_text(size=28, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28)
        )
Marey_plot
ggsave("figures/Marey_plot.png",
       device="png",dpi=320,units="cm",width=50,height=40)

rm(data.final)
rm(df)



#============================================================================
# Figure. Recombination map in barplot, all species & chromosomes
#============================================================================
# Data: recombination maps, loess, 1000 bins of equal size
# All species & chromosomes
# Tidy dataframe
filenames = list.files(path = "output/recombination_maps/loess/bins", pattern = ".txt", full.names=TRUE)
recombinationbins = data.frame(set = as.character(), species = as.character(), map = as.character(), bin = as.numeric(),
                               phys = as.numeric(), rec.rate = as.numeric())
for (file  in filenames) {
  data = read.table(file, header = TRUE)
  set = gsub("_chromosome[A-Za-z0-9]*.txt", "", gsub("output/recombination_maps/loess/bins/", "", file))
  species = as.character(metadata$species[which(metadata$id == set)])
  # Compute the scaled recombination rate (i.e. divide by max rec rate)
  df = data.frame(set = rep(set, 999),
                  species = rep(species, 999),
                  map = rep(gsub(".txt", "", gsub("output/recombination_maps/loess/bins/[A-Za-z0-9[:punct:]]*_chromosome", "", file)), 999),
                  bin = c(1:999),
                  phys = data$phys,
                  rec.rate = data$rec.rate/max(data$rec.rate, na.rm = TRUE))
  recombinationbins = rbind(recombinationbins,
        df)
}
write.table(recombinationbins, "output/recombination_maps/loess/AllMapsBins.txt", row.names = FALSE, col.names = TRUE)
recombinationbins = read.table("output/recombination_maps/loess/AllMapsBins.txt", header = TRUE)

# Deal with infinite values in log scale
recombinationbins$rec.rate = recombinationbins$rec.rate + 0.0001

# Pool recombination rates in 100 or 10 bins instead of 1000?
# Scale effect? To see global patterns instead of smaller variations


######### GGPLOT
# Plotting the barplot
# my_breaks = c(0.1, 0, -0.2, -0.4, -0.6, -0.8)
recombinationbins_barplot = ggplot(data = recombinationbins, aes(x = map, y = bin/(bin*1000), fill = rec.rate))+
  geom_bar(stat='identity', width = 1) +
  # scale_fill_manual(values = color) +
  # scale_fill_viridis_c(breaks = my_breaks, labels = my_breaks) +
  scale_fill_viridis_c(option = "viridis") +
  # scale_fill_gradient2() +
  facet_grid(~set, scales = "free", space="free_x") +
  labs(x="Chromosome", y="Physical position (bin)", fill="Recombination rate (cM/Mb)") +
  theme(axis.line = element_blank(),
        # axis.line.x = element_blank(), # No x axis
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=16),
        axis.title.y = element_text(color="black", size=16),
        axis.text=element_text(size=16, colour="black"),
        axis.text.x=element_blank(), # No samples names
        axis.ticks.x=element_blank(), # No x axis
        strip.text=element_text(size=12, colour="black", angle = 90),
        legend.position = "top",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=14),
        legend.title=element_text(size=16))
recombinationbins_barplot


#============================================================================
# Figure. Recombination map with genome architecture in heat map
#============================================================================
# "Heat maps below the curves of recombination rates indicate gene density (low for cold colors and high for hot colors)." (see Bauer et al. 2013)







#============================================================================
# Figure. Brockenstick with explanatory variables
#============================================================================
brockenstick = read.table("tables/brockenstick/brockenstick_k10.txt", header = TRUE, sep = "\t")
brockenstick = brockenstick[!(is.na(brockenstick$p1) | is.na(brockenstick$p2) | is.na(brockenstick$p3) | is.na(brockenstick$p4) | is.na(brockenstick$p5) | is.na(brockenstick$p6) | is.na(brockenstick$p7) | is.na(brockenstick$p8) | is.na(brockenstick$p9) | is.na(brockenstick$p10)),]

brockenstick$set = as.character(brockenstick$set)
brockenstick$chromosome = as.character(brockenstick$chromosome)
# Add the information about chromosome size in Mb
# Barplot under the structure plot
# df_chrsize = data.frame(set = as.character(brockenstick$set), chromosome = brockenstick$chromosome)
# df_chrsize$chr_size = NA
# for (i in 1:nrow(df_chrsize)) {
#   df_chrsize$chr_size[i] = chromosome.stats$phys.map.length[chromosome.stats$set == df_chrsize$set[i] & chromosome.stats$chromosome == df_chrsize$chromosome[i]]/1000000
# }
# df_chrsize$set = gsub("_", " ",regmatches(as.character(brockenstick$set), gregexpr("^[A-Z][a-z]+_[a-z]+", as.character(brockenstick$set))))
# df_chrsize$p = 1

# Add the information about the average recombination rate in cM/Mb
# Barplot under the structure plot
# df_avgrecrate = data.frame(set = as.character(brockenstick$set), chromosome = brockenstick$chromosome)
# df_avgrecrate$avgrecrate = NA
# for (i in 1:nrow(df_avgrecrate)) {
#   df_avgrecrate$avgrecrate[i] = chromosome.stats$mean.recrate[chromosome.stats$set == df_avgrecrate$set[i] & chromosome.stats$chromosome == df_avgrecrate$chromosome[i]]
# }
# df_avgrecrate$set = gsub("_", " ",regmatches(as.character(brockenstick$set), gregexpr("^[A-Z][a-z]+_[a-z]+", as.character(brockenstick$set))))
# df_avgrecrate$p = 1


# reformat data
brockenstick = melt(brockenstick)
brockenstick = cbind(paste(brockenstick$set, "_", brockenstick$chromosome, sep =""), brockenstick)

# Get species name
brockenstick$set = gsub("_", " ",regmatches(as.character(brockenstick$set), gregexpr("^[A-Z][a-z]+_[a-z]+", as.character(brockenstick$set))))
# Format columns in proper way
colnames(brockenstick) = c("sample", "set","chromosome","segment","proportion.length")

brockenstick$chr_size = NA
for (i in 1:nrow(brockenstick)) {
  brockenstick$chr_size[i] = chromosome.stats$phys.map.length[chromosome.stats$set == gsub("_[0-9A-Za-z]*$", "",brockenstick$sample[i]) & chromosome.stats$chromosome == brockenstick$chromosome[i]]/1000000
}
brockenstick$avgrecrate = NA
for (i in 1:nrow(brockenstick)) {
  brockenstick$avgrecrate[i] = chromosome.stats$mean.recrate[chromosome.stats$set == gsub("_[0-9A-Za-z]*$", "",brockenstick$sample[i]) & chromosome.stats$chromosome == brockenstick$chromosome[i]]
}
rbarintra = read.table("tables/df_geneticshuffling_physical.csv", header = TRUE)
brockenstick$rbarintra = NA
for (i in 1:nrow(brockenstick)) {
  brockenstick$rbarintra[i] = rbarintra$geneticshuffling[rbarintra$set == gsub("_[0-9A-Za-z]*$", "",brockenstick$sample[i])]
}
# Set a vector of gradient color
# Value of color is simply expected - p, the departure from the expected proportion 1/k
k = 10
brockenstick$color = (1/k) - brockenstick$proportion.length
# brockenstick$color[which(brockenstick$color > 0)] = rescale(brockenstick$color[which(brockenstick$color > 0)], to = c(0, 1), from = range(brockenstick$color[which(brockenstick$color > 0)], na.rm = TRUE))
# brockenstick$color[which(brockenstick$color < 0)] = rescale(brockenstick$color[which(brockenstick$color < 0)], to = c(-1, 0), from = range(brockenstick$color[which(brockenstick$color < 0)], na.rm = TRUE))
# brockenstick$color = log10(brockenstick$proportion.length/(1/k))


#-------------------------------------------------
# Add hierarchical Clustering for ordering species
#-------------------------------------------------
# Load 'h1', results of the clustering procedure in 'Classification.R'
load(file = "output/classification/hclustPooledChromosomes.Rda")

# h1$hclust$order
# Manual ordering
# new_order = c(21,22,23,10,17,12, 29, 16, 33,3, 15, 9, 25, 5, 30, 31, 32,  6, 14, 27,  2, 34, 13, 20, 26, 19,  8, 11, 28)
# new_order = c(13, 15, 19, 24, 14, 25, 30, 32, 27,  2,  5,  6, 29, 12, 26,  9,  3, 11, 8, 23, 22, 28, 17, 20,  7, 18,  4, 21, 10, 16, 31,  1)
# h1$hclust = rotate(h1$hclust, new_order)

dendr = dendro_data(h1, type="rectangle") # convert for ggplot
clust = cutree(h1, k = 4)                    # find 4 clusters
clust.df = data.frame(label = names(clust), cluster = factor(clust))
# dendr[["labels"]] has the labels, merge with clust.df based on label column
dendr[["labels"]] <- merge(dendr[["labels"]],clust.df, by="label")

# # # plot the dendrogram; note use of color=cluster in geom_text(...)
# ggplot() +
#   geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) +
#   geom_text(data=label(dendr), aes(x, y, label=label, hjust=0, color=cluster),
#             size=3) +
#   coord_flip() + scale_y_reverse(expand=c(0.2, 0)) +
#   theme(axis.line.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.title.y=element_blank(),
#         panel.background=element_rect(fill="white"),
#         panel.grid=element_blank())

brockenstick$set = factor(brockenstick$set, levels = gsub("_", " ", as.character(h1$labels[h1$order])))
#-------------------------------------------------
# Add classification to cluster similar species
# Manual clustering
# [1] "Arabidopsis thaliana"    "Brachypodium distachyon" "Brassica napus"          "Capsicum annuum"
#  [5] "Cenchrus americanus"     "Citrullus lanatus"       "Cucumis melo"            "Cucumis sativus"
#  [9] "Glycine max"             "Gossypium raimondii"     "Malus domestica"         "Manihot esculenta"      
# [13] "Oryza sativa"            "Phaseolus vulgaris"      "Prunus mume"             "Sesamum indicum"
# [17] "Setaria italica"         "Solanum lycopersicum"    "Solanum tuberosum"       "Sorghum bicolor"
# [21] "Theobroma cacao"         "Triticum aestivum"       "Vitis vinifera"          "Zea mays"  
brockenstick$set = factor(brockenstick$set, levels = c("Arabidopsis thaliana", "Brachypodium distachyon", "Brassica napus",
                                                       "Citrullus lanatus", "Cucumis sativus",
                                                       "Malus domestica", "Manihot esculenta", "Oryza sativa", "Populus deltoides",
                                                       "Populus simonii",
                                                       "Prunus mume", "Sesamum indicum", "Theobroma cacao", "Vitis vinifera", "Cucumis melo",
                                                       "Setaria italica","Coffea canephora",
                                                       "Lupinus angustifolius", "Glycine max","Gossypium hirsutum",
                                                       "Gossypium raimondii", "Phaseolus vulgaris","Solanum lycopersicum",
                                                       "Solanum tuberosum", "Sorghum bicolor",
                                                       "Arachis duranensis", "Capsicum annuum", "Cenchrus americanus", "Zea mays",
                                                       "Hordeum vulgare", "Triticum aestivum", "Triticum dicoccoides"
))
df_chrsize$set = factor(df_chrsize$set, levels = c("Arabidopsis thaliana", "Brachypodium distachyon", "Brassica napus",
                                                       "Citrullus lanatus", "Cucumis sativus",
                                                       "Malus domestica", "Manihot esculenta", "Oryza sativa", "Populus deltoides",
                                                       "Populus simonii",
                                                       "Prunus mume", "Sesamum indicum", "Theobroma cacao", "Vitis vinifera", "Cucumis melo",
                                                       "Setaria italica","Coffea canephora",
                                                       "Lupinus angustifolius", "Glycine max","Gossypium hirsutum",
                                                       "Gossypium raimondii", "Phaseolus vulgaris","Solanum lycopersicum",
                                                       "Solanum tuberosum", "Sorghum bicolor",
                                                       "Arachis duranensis", "Capsicum annuum", "Cenchrus americanus", "Zea mays",
                                                       "Hordeum vulgare", "Triticum aestivum", "Triticum dicoccoides"
                                                       ))
# df_avgrecrate$set = factor(df_avgrecrate$set, levels = c("Arabidopsis thaliana", "Brachypodium distachyon", "Brassica napus",
#                                                        "Citrullus lanatus", "Cucumis sativus",
#                                                        "Malus domestica", "Manihot esculenta", "Oryza sativa", "Populus deltoides",
#                                                        "Populus simonii",
#                                                        "Prunus mume", "Sesamum indicum", "Theobroma cacao", "Vitis vinifera", "Cucumis melo",
#                                                        "Setaria italica","Coffea canephora",
#                                                        "Lupinus angustifolius", "Glycine max","Gossypium hirsutum",
#                                                        "Gossypium raimondii", "Phaseolus vulgaris","Solanum lycopersicum",
#                                                        "Solanum tuberosum", "Sorghum bicolor",
#                                                        "Arachis duranensis", "Capsicum annuum", "Cenchrus americanus", "Zea mays",
#                                                        "Hordeum vulgare", "Triticum aestivum", "Triticum dicoccoides"
#                                                        ))

######### GGPLOT
# Plotting the barplot
my_breaks = c(0.1, 0, -0.2, -0.4, -0.6, -0.8)
BrockenStick_barplot = ggplot(data = brockenstick, aes(x=sample, y=proportion.length, fill = color))+
  geom_bar(stat='identity', width = 1) +
  # scale_fill_manual(values = color) +
  scale_fill_viridis_c(breaks = my_breaks, labels = my_breaks) +
  # scale_fill_gradient2() +
  facet_grid(~set, scales = "free", space="fixed") +
  labs(x="Chromosome", y="Proportion of\ntotal physical length", fill="Segment relative length") +
  theme(axis.line = element_blank(),
        # axis.line.x = element_blank(), # No x axis
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=16, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=16,hjust = 0.5),
        axis.title.x = element_text(color="black", size=24),
        axis.title.y = element_text(color="black", size=24),
        axis.text=element_text(size=24, colour="black"),
        axis.text.x=element_blank(), # No samples names
        axis.ticks.x=element_blank(), # No x axis
        strip.text=element_text(size=24, colour="black", angle = 90),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=24),
        legend.title=element_text(size=24))
BrockenStick_barplot


# plotting the dendrogram of species ordering and grouping
# Breaks are plotted as a function of facets coordinates in BrockenStick_barplot
# Facets of unequal size
# Position is the mean chromosome number of the facet (center of the group) + position of the first chromosome in the group (in number of chromosomes before)
# labels = h1$labels[h1$order]
# facet_center = NA
# for (i in 1:length(labels)) {
#   toMatch = labels[1:i]
#   facet_center[i] = sum((grepl(paste(toMatch,collapse="|"), brockenstick$sample) & brockenstick$segment == "p1"), na.rm = TRUE)
# }
# 
dendroPlot = ggplot() +
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) +
  labs(x="", y="\n", fill="") +
  theme(axis.line=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.title=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())

# dendroPlot = ggdendrogram(h1, rotate = FALSE, size = 2) + theme_dendro()
# dendroPlot
ggarrange(dendroPlot, BrockenStick_barplot, nrow = 2, heights = c(5,20))

my_breaks = c(20, 50, 500)
ChrSize_barplot = ggplot(data = brockenstick, aes(x=sample, y=proportion.length, fill = chr_size))+
  geom_bar(stat='identity', width = 1) +
  # scale_fill_manual(values = color) +
  scale_fill_viridis_c(breaks = my_breaks, labels = my_breaks, trans = "log", option = "inferno", direction = -1) +
  # scale_fill_gradient2() +
  facet_grid(~set, scales = "free", space="fixed") +
  labs(x="", y="\n", fill="Chromosome size (Mb)") +
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
        axis.text=element_text(size=24, colour="white"),
        axis.text.x=element_blank(), # No samples names
        axis.ticks=element_blank(), # No axis ticks
        strip.text=element_blank(),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=24),
        legend.title=element_text(size=24))
# ChrSize_barplot

my_breaks = c(0, 1, 2, 5, 10, 15)
AvgRecRate_barplot = ggplot(data = brockenstick, aes(x=sample, y=proportion.length, fill = avgrecrate))+
  geom_bar(stat='identity', width = 1) +
  scale_fill_viridis_c(breaks = my_breaks, labels = my_breaks, trans = "log", option = "plasma") +
  facet_grid(~set, scales = "free", space="fixed") +
  labs(x="", y="\n", fill="Mean recombination rate (cM/Mb)") +
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
        axis.text=element_text(size=24, colour="white"),
        axis.text.x=element_blank(), # No samples names
        axis.ticks=element_blank(), # No axis ticks
        strip.text=element_blank(),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=24),
        legend.title=element_text(size=24))
# AvgRecRate_barplot

my_breaks = c(0, 0.02, 0.04, 0.06, 0.08 ,0.1)
rbar_barplot = ggplot(data = brockenstick, aes(x=sample, y=proportion.length, fill = rbarintra))+
  geom_bar(stat='identity', width = 1) +
  scale_fill_viridis_c(option = "viridis") +
  facet_grid(~set, scales = "free", space = "fixed") +
  labs(x="", y="\n", fill="Intra chromosomal genetic shuffling") +
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
        axis.text=element_text(size=24, colour="white"),
        axis.text.x=element_blank(), # No samples names
        axis.ticks=element_blank(), # No axis ticks
        strip.text=element_blank(),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white", size = 1),
        legend.key.width=unit(2,"cm"),
        legend.text=element_text(size=24),
        legend.title=element_text(size=24))
# rbar_barplot

plot = plot_grid(BrockenStick_barplot, ChrSize_barplot, AvgRecRate_barplot, rbar_barplot,
                 nrow = 4, rel_heights = c(4, 1, 1, 1))

#
png("figures/BrockenStick_barplot_clusters.png", width = 2000, height = 1600)
grid.newpage()
print(dendroPlot, vp = viewport(x = 0.525, y = 0.9, width = 1.03, height = 0.2))
print(plot, vp = viewport(x = 0.5, y = 0.4, width = 1, height = 0.8))
dev.off()
# 

# library("cowplot")
# plot_grid(dendroPlot, BrockenStick_barplot, ChrSize_barplot, AvgRecRate_barplot, rbar_barplot,
#           align = "v", axis = "r",
#           nrow = 5, rel_heights = c(1, 4, 1, 1, 1))
png("figures/BrockenStick_barplot_clusters.png", width = 2000, height = 1600)
grid.newpage()
print(plot_grid(BrockenStick_barplot, ChrSize_barplot, AvgRecRate_barplot, rbar_barplot,
          nrow = 4, rel_heights = c(4, 1, 1, 1)))
dev.off()
# ggarrange(dendroPlot, BrockenStick_barplot, nrow = 2)



#============================================================================
# Figure. Lorenz curves of recombination
#============================================================================


#============================================================================
# End
#============================================================================