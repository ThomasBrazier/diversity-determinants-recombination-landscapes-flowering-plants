############################################################################ #
#                    ECOBIO - PhD
#
#       Fit a distributionot the Marey plots - Estimates parameters
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
gc() #free up memory and report the memory usage.

# Loading packages
library(rstudioapi)
library(ggplot2)
library(MareyMap)
library(stringr) 
library(fitdistrplus)
library(optimx)
library(xlsx)

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
source("sources/ModelFit.R")

# df = read.table(file = "data-cleaned/Marey_plot.scaled.txt", header = TRUE, sep = "\t")
# Load the joint database with all maps
data.final = read.table(file = "data-cleaned/marey_maps/AllMaps.txt", header = TRUE, sep = "\t")

chromosome.stats = read.table(file = "tables/chromosome.stats.csv", header = TRUE, sep = ";")

#============================================================================#
# Gamma Telomere model of recombination (Sylvain & Thomas) ----
#============================================================================#

# Equation of the PDF in mathematica:
# p g[\[Alpha]1, \[Beta]1, z] + (1 - p) g[\[Alpha]2, \[Beta]2, L - z]

# Two Gamma distributions with two parameters (Alpha and Beta) with O being the start in telomeres
# \[Alpha] is the shape parameter of the distribution
# \[Alpha] = 1 is exponential distribution
# With \[Alpha] >>1 recombination is distributed away from telomeres
# p the coefficient of assymetry; p is the probability of recombination in left part of the chromosome
# (1-p) is the probability of recombination in right part of the chromosome
# Chromosome of size L
# Conditions were
# cond = {\[Alpha]1 > 0, \[Alpha]2 > 0, \[Beta]1 > 0, \[Beta]2 > 0, 0 < p < 1, L > 0, 0 <= x <= L}
# hence a PDF with six parameters: \[Alpha]1, \[Beta]1, \[Alpha]2, \[Beta]2, p, L
Alpha1 = 1
Alpha2 = 1
Beta1 = 1
Beta2 = 1
p = 0.5
L = 25

# PDF
plot(Gamma_Telomere_PDF(Alpha1, Alpha2, Beta1, Beta2, p, L), type = 'l')

# CDF
plot(Gamma_Telomere_CDF(Alpha1, Alpha2, Beta1, Beta2, p, L), type = 'l')

# CDF scaled on [0; 1], hence the genetic distance must be scaled (divided by total genetic length)

# Explore the parameter space manually
# Conditions:
# alpha >= 1
plot(Gamma_Telomere_CDF(Alpha1, Alpha2, 2, 2, p, L), type = 'l')
plot(Gamma_Telomere_CDF(Alpha1, Alpha2, 4, 4, p, L), type = 'l')
plot(Gamma_Telomere_CDF(Alpha1, Alpha2, 6, 6, p, L), type = 'l')

plot(Gamma_Telomere_CDF(1, 1, 2, 2, p, L), type = 'l')
plot(Gamma_Telomere_CDF(2, 2, 4, 4, p, L), type = 'l')
plot(Gamma_Telomere_CDF(0.5, 0.5, 6, 6, p, L), type = 'l')

plot(Gamma_Telomere_CDF(Alpha1, Alpha2, 2, 2, p, 20), type = 'l')
plot(Gamma_Telomere_CDF(Alpha1, Alpha2, 2, 2, p, 40), type = 'l')
plot(Gamma_Telomere_CDF(Alpha1, Alpha2, 2, 2, p, 60), type = 'l')

plot(Gamma_Telomere_CDF(20, 20, 20, 20, p, L), type = 'l')
plot(Gamma_Telomere_CDF(2, 2, 4, 4, p, L), type = 'l')

#============================================================================#
# Predictions with the Gamma Telomere Model that fit hypotheses ----
#============================================================================#

# Chromosome size L
# Chromosome size seems an important predictor of the mean recombination rate
# and the recombination landscape global pattern
# Compare two chromosomes with change in chromosome size while other parameter were fixed
# Both Gamma distributions set to exponential distribution (Alpha = 1)
# Symetry of the landscape (p = O.5)

# Represent the two chromosomes on the same scale
L = 5
p = 0.5
Alpha1 = 1
Alpha2 = 1
plot(Gamma_Telomere_CDF(Alpha1, Alpha2, Beta1, Beta2, p, L)$x/L,
     Gamma_Telomere_CDF(Alpha1, Alpha2, Beta1, Beta2, p, L)$cdf,
     xlab = "Relative physical position", xlim = c(0,1), ylim = c(0, 1),
     ylab = "CDF", col = "Black", type = "l")
L = 10
lines(Gamma_Telomere_CDF(Alpha1, Alpha2, Beta1, Beta2, p, L)$x/L,
     Gamma_Telomere_CDF(Alpha1, Alpha2, Beta1, Beta2, p, L)$cdf,
     xlab = "Relative physical position",
     ylab = "CDF", col = "Blue")
L = 20
lines(Gamma_Telomere_CDF(Alpha1, Alpha2, Beta1, Beta2, p, L)$x/L,
     Gamma_Telomere_CDF(Alpha1, Alpha2, Beta1, Beta2, p, L)$cdf,
     xlab = "Relative physical position",
     ylab = "CDF", col = "Red")
legend(x = "bottomright", legend = c("L = 5", "L = 10", "L = 20"),
       col = c("Black", "Blue", "Red"), lty = 1)

#----------------------------------------------------------------------------#
# Assymetry, p
L = 5
p = 0.5
Alpha1 = 1
Alpha2 = 1
plot(Gamma_Telomere_CDF(Alpha1, Alpha2, Beta1, Beta2, p, L)$x/L,
     Gamma_Telomere_CDF(Alpha1, Alpha2, Beta1, Beta2, p, L)$cdf,
     xlab = "Relative physical position", xlim = c(0,1), ylim = c(0, 1),
     ylab = "CDF", col = "Black", type = "l")
p = 0.25
lines(Gamma_Telomere_CDF(Alpha1, Alpha2, Beta1, Beta2, p, L)$x/L,
      Gamma_Telomere_CDF(Alpha1, Alpha2, Beta1, Beta2, p, L)$cdf,
      xlab = "Relative physical position",
      ylab = "CDF", col = "Blue")
p = 0.75
lines(Gamma_Telomere_CDF(Alpha1, Alpha2, Beta1, Beta2, p, L)$x/L,
      Gamma_Telomere_CDF(Alpha1, Alpha2, Beta1, Beta2, p, L)$cdf,
      xlab = "Relative physical position",
      ylab = "CDF", col = "Red")
legend(x = "bottomright", legend = c("p = 0.5", "p = 0.25", "p = 0.75"),
       col = c("Black", "Blue", "Red"), lty = 1)


#----------------------------------------------------------------------------#
# Shape parameter, Alpha
# Departure from 1 is departure from the exponential distribution
L = 5
p = 0.5
Alpha1 = 1
Alpha2 = 1
plot(Gamma_Telomere_CDF(Alpha1, Alpha2, Beta1, Beta2, p, L)$x/L,
     Gamma_Telomere_CDF(Alpha1, Alpha2, Beta1, Beta2, p, L)$cdf,
     xlab = "Relative physical position", xlim = c(0,1), ylim = c(0, 1),
     ylab = "CDF", col = "Black", type = "l")
Alpha1 = 2
Alpha2 = 2
lines(Gamma_Telomere_CDF(Alpha1, Alpha2, Beta1, Beta2, p, L)$x/L,
      Gamma_Telomere_CDF(Alpha1, Alpha2, Beta1, Beta2, p, L)$cdf,
      xlab = "Relative physical position",
      ylab = "CDF", col = "Blue")
Alpha1 = 4
Alpha2 = 4
lines(Gamma_Telomere_CDF(Alpha1, Alpha2, Beta1, Beta2, p, L)$x/L,
      Gamma_Telomere_CDF(Alpha1, Alpha2, Beta1, Beta2, p, L)$cdf,
      xlab = "Relative physical position",
      ylab = "CDF", col = "Red")
legend(x = "bottomright", legend = c("Alpha = 1", "Alpha = 2", "Alpha = 4"),
       col = c("Black", "Blue", "Red"), lty = 1)

# With Alpha > 4, CDF does not converge on 1; y limits are 0.95 (Alpha = 2) and 0.75 (Alpha = 4)
# Not working with Alpha < 1; infinite values

#============================================================================#
# Fitting the Gamma Telomere Model on Marey empirical data ----
#============================================================================#
# data.final$set = as.character(data.final$set)
# data.final$map = as.character(data.final$map)

# 
# source("sources/ModelFit.R")
# # Three test dataset
# set = "Arabidopsis_thaliana_Serin2017"
# set= "Zea_mays_IBM_MaizeSNP50"
# map = "2"
# X = as.numeric(as.character(data.final$phys[data.final$set == set & data.final$map == map]))/1000000
# Y = data.final$gen[data.final$set == set & data.final$map == map]
# plot(X, Y)
# 
# # X = df$phys[df$group == "Triticum aestivum 1B"]
# # Y = df$gen[df$group == "Triticum aestivum 1B"]
# 
# # X = df$phys[df$group == "Malus domestica 12"]
# # Y = df$gen[df$group == "Malus domestica 12"]
# 
# # # Sort X/Y by X (ascending) to have a continuous ascending function from 0 to L
# Y = Y[order(X)]
# X = X[order(X)]
# 
# # TODOne Work only on a scaled dataset
# Ys = Y/max(Y, na.rm = TRUE)
# removing outliers 0 values
# Y[Y == 0] = 0.0001
# Y[Y == 1] = 0.9999


# Calibrating starting parameters
# The model proved to be very sensitive to starting values, with multiple local minima
# Alpha1 = 2
# Alpha2 = 2
# Beta1 = 2
# Beta2 = 2
# p = 0.5
# L = max(X)

# plot(X, Ys, ylim  = c(0, 1))
# lines(Gamma_Telomere_ECDF(Alpha1, Alpha2, Beta1, Beta2, p, L, X))

# Beware, computation can be long
# fit.optim = optim(par = c(Alpha1, Alpha2, Beta1, Beta2, p), fitcdf, method = "L-BFGS-B", lower = 0.01, control = list(trace = 1, maxit = 100))
# gc()

# Careful to the boundaries: 0.1 to 10
# Besides, fitting is very sensitive to starting parameters
  # TODO Iterative procedure with parameters randomly drawed from a prior distribution?
# fit.optim = optim(par = c(Alpha1, Alpha2, Beta1, Beta2, p), fitcdf, method = "L-BFGS-B", lower = 0.1, upper = 10, control = list(trace = 1, maxit = 1000))

# Method "L-BFGS-B" that can handle box (or bounds) constraints.
# fit.optim = optim(par = c(Alpha1, Alpha2, Beta1, Beta2, p), fitcdf, method = "L-BFGS-B", lower = c(0.1, 0.1, 0.1, 0.1, 0.49), upper = c(10, 10, 10, 10, 0.51),
#                   control = list(trace = 1, maxit = 1000))
# fit.optim = optimx(par = c(Alpha1, Alpha2, Beta1, Beta2, p), fitcdf, method = "L-BFGS-B", hessian = TRUE, lower = c(0.1, 0.1, 0.1, 0.1, 0.49), upper = c(10, 10, 10, 10, 0.51),
                  # control = list(trace = 1, maxit = 1000))
# nlm with boundaries
# fit.optim = nlminb(c(Alpha1, Alpha2, Beta1, Beta2, p), fitcdf, lower = c(0.1, 0.1, 0.1, 0.1, 0.49), upper = c(10, 10, 10, 10, 0.51))

# Parameters:
# fit.optim$par
# fit.optim
# Plot of the observed ECDF with the fitted beta CDF (parameters fitted on the ECDF)
# plot(X, Ys)
# lines(Gamma_Telomere_CDF(fit.optim$p1,  fit.optim$p2, fit.optim$p3, fit.optim$p4,
                          # fit.optim$p5, L), col = "Green")


# Explore parameter space to restrain priors of the Gamma distribution to more fittable ranges
    # Done in Mathematica with interactive graphical representations



# Whatever the method used, multiple local minima can't be disentangle from varying starting values
# Multiple starts from randomized initializations, take the best-of-N results.
  # Look for the posterior mode with a random draw of multiple starting values?
  # N (default = 1,000) iterations and selection of the set of parameters minimizing the value
  # Explore the entire parameter space to search for the best local minima

# Save the five parameters and the associated criterion
# A df of 6 columns and N lines
# N = 100
# iterated_parameters = data.frame(p1 = numeric(N), p2 = numeric(N), p3 = numeric(N), p4 = numeric(N),
#                                   p5 = numeric(N), crit = numeric(N))
# Search parameters among N random set of starting values
# for (i in 1:N) {
#   print(i)
#   parameters = c(runif(5, min = c(1, 1, 0.1, 0.1, 0.49), max = c(10, 10, 10, 10, 0.51)))
#   fit.optim = optimx(par = parameters, fitcdf, method = "L-BFGS-B", hessian = TRUE, lower = c(1, 1, 0.1, 0.1, 0.3),
#                      upper = c(10, 10, 10, 10, 0.7),
#                      control = list(trace = 1, maxit = 1000))
#   iterated_parameters$p1[i] = fit.optim$p1
#   iterated_parameters$p2[i] = fit.optim$p2
#   iterated_parameters$p3[i] = fit.optim$p3
#   iterated_parameters$p4[i] = fit.optim$p4
#   iterated_parameters$p5[i] = fit.optim$p5
#   iterated_parameters$crit[i] = fit.optim$value
# }
# Explore "posterior" parameters distributions to see if there is values stuck to boundaries ("priors" too restricted)
# hist(iterated_parameters$p1, breaks = 20)
# hist(iterated_parameters$p2, breaks = 20)
# hist(iterated_parameters$p3, breaks = 20)
# hist(iterated_parameters$p4, breaks = 20)
# hist(iterated_parameters$p5, breaks = 20)
# 
# source("sources/ModelFit.R")
# # Iterative fitting to find the best local minima out of N random sets of starting parameters
# N = 100
# iterated_parameters = iterative_fitcdf()
# # Select the best out of N runs
# best_parameter =  iterated_parameters[which.min(iterated_parameters$crit),]
# # best_parameter = data.frame(p1 = 2, p2 = 2, p3 = 10, p4 = 10, p5 = 0.5)
# best_parameter
# # Plot of the observed ECDF with the fitted CDF (parameters fitted on the ECDF)
# plot(X, Ys)
# lines(Gamma_Telomere_ECDF(best_parameter$p1,  best_parameter$p2, best_parameter$p3, best_parameter$p4,
#                          best_parameter$p5, L, X), col = "Red")
# lines(Gamma_Telomere_ECDF(best_parameter$p1,  best_parameter$p2, best_parameter$p3, best_parameter$p4,
#                           best_parameter$p5, L, seq(0, L, by = 0.1)), col = "Green")
# lines(Gamma_Telomere_CDF(best_parameter$p1,  best_parameter$p2, best_parameter$p3, best_parameter$p4,
#                           best_parameter$p5, L), col = "Blue")

# CDF = Gamma_Telomere_CDF(best_parameter$p1,  best_parameter$p2, best_parameter$p3, best_parameter$p4,
#                          best_parameter$p5, L)
# 
# ECDF = Gamma_Telomere_ECDF(best_parameter$p1,  best_parameter$p2, best_parameter$p3, best_parameter$p4,
#                            best_parameter$p5, L, X)

# source("sources/ModelFit.R")
# best_parameter = data.frame(p1 = 2, p2 = 2, p3 = 10, p4 = 10, p5 = 0.5)
# plot(Gamma_Telomere_CDF(best_parameter$p1,  best_parameter$p2, best_parameter$p3, best_parameter$p4,
#                           best_parameter$p5, L), col = "Red")

# Plot of the observed ECDF with the fitted beta CDF (parameters fitted on the PDF)
# plot(X, Y)
# x = seq(0, L, by = 0.1)
# lines(Gamma_Telomere_CDF(fit.optim$par[1],  fit.optim$par[2], fit.optim$par[3], fit.optim$par[4],
#                          fit.optim[5]$par), add=T, col = "Green")





# Fitting procedure seems to have anormal behavior when trying to fit p, or on not so simple Marey maps (e.g. Z. mays)
# Big influence of starting values, may be multiple local minima
# Hence model selection and parameter sensitivity/selection
# Or directly Bayesian estimates

# Hierarchical procedure?
  # Estimate parameters with fixed p, then estimate p with fixed parameters and see if it significantly improve the fit
# Draw the distribution of the parameters estimated after 10,000 random draws of initial values
  # Iterative process: replace starting values by new fitted values if convergence criterion < previous convergence criterion


# Is it possible to fit to unscaled data? i.e. Y predictions not in the range 0-1
# plot(X, Y)
# ECDF = Gamma_Telomere_ECDF(fit.optim$par[1],  fit.optim$par[2], fit.optim$par[3], fit.optim$par[4],
#                     fit.optim$par[5], L, X)
# lines(ECDF$x, ECDF$cdf*max(Y, na.rm = TRUE), col = "Green")



#-----------------------------------------------------#
# Process all maps and save parameters inferred  ----
#-----------------------------------------------------#
source("sources/ModelFit.R")
# Import data
data.final = read.table(file = "data-cleaned/marey_maps/AllMaps.txt", header = TRUE, sep = "\t")
data.final$set = as.character(data.final$set)
data.final$map = as.character(data.final$map)

# Initialize the data frame which will store results
sets = unique(paste(data.final$set, data.final$map, sep = ":"))
nsets = length(sets)
parameters_GammaModel = data.frame(set = character(nsets), map = character(nsets), L = NA, a1 = NA, a2 = NA, b1 = NA,
                                   b2 = NA, p = NA, crit = NA)
parameters_GammaModel$set = gsub(":[A-Za-z0-9]*", "", sets)
parameters_GammaModel$map = gsub("[_A-Za-z0-9]*:", "", sets)

# Batch fitting to the whole dataset
for (j in 1:nrow(parameters_GammaModel)) {
  set = as.character(parameters_GammaModel$set[j])
  chromosome = parameters_GammaModel$map[j]
  cat(as.character(parameters_GammaModel$set[j]), "chromosome", as.character(parameters_GammaModel$map[j]), "\n")
  # Import data
  data.final = read.table(file = "data-cleaned/marey_maps/AllMaps.txt", header = TRUE, sep = "\t")
  data.final$set = as.character(data.final$set)
  data.final$map = as.character(data.final$map)
  
  X = as.numeric(as.character(data.final$phys[data.final$set == set & data.final$map == chromosome]))/1000000
  Y = data.final$gen[data.final$set == set & data.final$map == chromosome] 
  # Sort X/Y by X (ascending) to have a continuous ascending function from 0 to L
  Y = Y[order(X)]
  X = X[order(X)]
  # Work only on a scaled dataset for the recombination
  Ys = Y/max(Y, na.rm = TRUE)
  # Standard starting parameter
  L = max(X)
  
  results = fit_Mareymap(as.character(parameters_GammaModel$set[j]),
                     as.character(parameters_GammaModel$map[j]), niter = 100,
                     X, Ys, L)
  
  # Plot of the observed ECDF with the fitted CDF (parameters fitted on the ECDF)
  plot(X, Ys, main = paste(set, "chromosome", chromosome, sep = " "))
  lines(Gamma_Telomere_ECDF(results$param$a1,  results$param$a2, results$param$b1, results$param$b2,
                            results$param$p, L, X), col = "Red", lwd = 2)
  lines(Gamma_Telomere_ECDF(results$param$a1,  results$param$a2, results$param$b1, results$param$b2,
                            results$param$p, L, seq(0, L, by = 0.1)),
        col = "Green", lty = "dashed", lwd = 2)
  # Add the putative centromere position, p
  abline(v = results$param$p*results$L)
  
  # Save results in the dataset
  parameters_GammaModel$L[j] = results$L
  parameters_GammaModel$a1[j] = results$param$a1
  parameters_GammaModel$a2[j] = results$param$a2
  parameters_GammaModel$b1[j] = results$param$b1
  parameters_GammaModel$b2[j] = results$param$b2
  parameters_GammaModel$p[j] = results$param$p
  parameters_GammaModel$crit[j] = results$param$crit
}

write.table(parameters_GammaModel, file = "tables/parameters_GammaModel.txt", col.names = TRUE, row.names = FALSE,
            sep = "\t", quote = FALSE)



#============================================================================#
# Parameters estimates ----
#============================================================================#
parameters_GammaModel = read.table("tables/parameters_GammaModel.txt", header = TRUE, sep = "\t")

# Explore parameters distributions
plot(density(parameters_GammaModel$p1, na.rm = TRUE), main = "", lwd = 2)
plot(density(parameters_GammaModel$p2, na.rm = TRUE), main = "", lwd = 2)
plot(density(parameters_GammaModel$p3, na.rm = TRUE), main = "", lwd = 2)
plot(density(parameters_GammaModel$p4, na.rm = TRUE), main = "", lwd = 2)
plot(density(parameters_GammaModel$p5, na.rm = TRUE), main = "", lwd = 2)
# Folded distribution of assymmetry parameter
# i.e. p = 1-p if p>0.5, to get values between 0 and 0.5 (proportion of the shorter arm)
parameters_GammaModel$p5[which(parameters_GammaModel$p5 > 0.5)] = 1 - parameters_GammaModel$p5[which(parameters_GammaModel$p5 > 0.5)]
plot(density(parameters_GammaModel$p5, na.rm = TRUE), main = "", lwd = 2)
hist(parameters_GammaModel$p5, main = "", lwd = 2, breaks = 50)
# Length of chromosomes (not estimated)
plot(density(parameters_GammaModel$L, na.rm = TRUE), main = "", lwd = 2)

# FIGURES
# Export descriptive figures
png("figures/GammaModel_alpha1.png", width = 600, height = 400)
# Alpha 1
plot(density(parameters_GammaModel$p1, na.rm = TRUE), main = "", lwd = 2, cex = 2)
dev.off()
# Beta 1
png("figures/GammaModel_beta1.png", width = 600, height = 400)
plot(density(parameters_GammaModel$p3, na.rm = TRUE), main = "", lwd = 2)
dev.off()
# Folded p
png("figures/GammaModel_p.png", width = 600, height = 400)
plot(density(parameters_GammaModel$p5, na.rm = TRUE), main = "", lwd = 2)
dev.off()


# Figures of parameters distributions per phylogenetic groups
# Import phylogenetic families
phylofamilies = read.table("data-cleaned/species_metadata.csv", header = TRUE, sep = ";")
parameters_GammaModel$family = NA
for (i in 1:nrow(parameters_GammaModel)) {
  if(grepl("Zea_mays", as.character(parameters_GammaModel$set[i]))) {
    parameters_GammaModel$family[i] = as.character(phylofamilies$family[which(phylofamilies$species == "Zea_mays")])
  } else {
    parameters_GammaModel$family[i] = as.character(phylofamilies$family[which(phylofamilies$species == gsub("_[A-Za-z0-9]*$", "", parameters_GammaModel$set[i]))])
  }
}
# Figures with species grouped in families
boxplot(parameters_GammaModel$p1~parameters_GammaModel$family)
boxplot(parameters_GammaModel$p2~parameters_GammaModel$family)
boxplot(parameters_GammaModel$p3~parameters_GammaModel$family)
boxplot(parameters_GammaModel$p4~parameters_GammaModel$family)
boxplot(parameters_GammaModel$p5~parameters_GammaModel$family)
boxplot(parameters_GammaModel$L~parameters_GammaModel$family)


# Batch figures of Marey maps with fitted distributions
# One figure per chromosome
for (i in 1:nrow(parameters_GammaModel)) {
  cat(as.character(parameters_GammaModel$set[i]), "chromosome", parameters_GammaModel$map[i], "\n")
  # Import data
  set = parameters_GammaModel$set[i]
  map = parameters_GammaModel$map[i]
  X = as.numeric(as.character(data.final$phy[data.final$set == set & data.final$map == map]))/1000000
  Y = data.final$gen[data.final$set == set & data.final$map == map] 
  # Sort X/Y by X (ascending) to have a continuous ascending function from 0 to L
  Y = Y[order(X)]
  X = X[order(X)]
  # Work only on a scaled dataset for the recombination
  Ys = Y/max(Y, na.rm = TRUE)
  # Standard sarting parameter
  L = max(X)
  
  if (!is.na(parameters_GammaModel$p1[i])) {
    # Plot of the observed ECDF with the fitted CDF (parameters fitted on the ECDF)
    png(paste("output/model_gamma_telomere/figures_fit/Gammatelomere_",set, "_", map,"_fit.png", sep =""), width = 600, height = 600)
    plot(X, Ys, main = paste(parameters_GammaModel$set[i], "chromosome", parameters_GammaModel$map[i], sep = " "))
    # lines(Gamma_Telomere_ECDF(parameters_GammaModel$p1[i],  parameters_GammaModel$p2[i], parameters_GammaModel$p3[i], parameters_GammaModel$p4[i],
    #                           parameters_GammaModel$p5[i], L, X), col = "Red")
    lines(Gamma_Telomere_ECDF(parameters_GammaModel$p1[i],  parameters_GammaModel$p2[i], parameters_GammaModel$p3[i], parameters_GammaModel$p4[i],
                              parameters_GammaModel$p5[i], L, seq(0, L, by = 0.1)), lwd = 4, col = "Green")
    dev.off()
  } else {
    cat("Non evaluated parameters (NA)\n")
  }
}

# Explore relationships between model parameters and recombination rates
chromosome.stats = read.table("tables/chromosome.stats.csv", header = TRUE)
parameters_GammaModel = read.table("tables/parameters_GammaModel.txt", header = TRUE, sep = "\t")
chromosome.stats$set = as.character(chromosome.stats$set)
colnames(chromosome.stats)[2] = "map"
parameters_GammaModel$set = as.character(parameters_GammaModel$set)

df = merge(chromosome.stats, parameters_GammaModel)

plot(df$phys.map.length, df$p1, main = "", log = "x", xlab = "Chromosome size (Mb, log-scaled)", pch = 19)
plot(df$phys.map.length, df$p2, main = "", log = "x", xlab = "Chromosome size (Mb, log-scaled)", pch = 19)
plot(df$phys.map.length, df$p3, main = "", log = "x", xlab = "Chromosome size (Mb, log-scaled)", pch = 19)
plot(df$phys.map.length, df$p4, main = "", log = "x", xlab = "Chromosome size (Mb, log-scaled)", pch = 19)
plot(df$phys.map.length, df$p5, main = "", log = "x", xlab = "Chromosome size (Mb, log-scaled)", pch = 19)
# None relationship convincing among mean rec rate, cv rec rate...
# Only the physical size of the chromosome seems to correlate with parameters of the model
# Shape might really be driven by the chromosome size
plot(df$mean.recrate, df$p1, main = "", pch = 19)
plot(df$mean.recrate, df$p2, main = "", pch = 19)
plot(df$mean.recrate, df$p3, main = "", pch = 19)
plot(df$mean.recrate, df$p4, main = "", pch = 19)
plot(df$mean.recrate, df$p5, main = "", pch = 19)

plot(df$sd.recrate/df$mean.recrate, df$p1, main = "", pch = 19)
plot(df$sd.recrate/df$mean.recrate, df$p2, main = "", pch = 19)
plot(df$sd.recrate/df$mean.recrate, df$p3, main = "", pch = 19)
plot(df$sd.recrate/df$mean.recrate, df$p4, main = "", pch = 19)
plot(df$sd.recrate/df$mean.recrate, df$p5, main = "", pch = 19)



#============================================================================#
# Testing hypotheses ----
#============================================================================#

# Compare the p parameter estimated from Marey maps to the centromere position from litterature
karyotypes = read.xlsx("data-cleaned/Karyotypes.xlsx", sheetIndex = 1)
# karyotypes$sp = paste(karyotypes$Genus, karyotypes$Species, sep = "_")
# Use the centromeric index as a measure of the centromere position
karyotypes$Centromeric_Index = as.numeric(as.character(karyotypes$Centromeric_Index))
plot(density(karyotypes$Centromeric_Index, na.rm = TRUE), main = "Density of Centromeric Index (black)\n and parameter p (red)")
lines(density(parameters_GammaModel$p5, na.rm = TRUE), col = "red")
# Compare estimated parameter with 'reference' centromere position
# Merging datasets
parameters_GammaModel = read.table("tables/parameters_GammaModel.txt", header = TRUE, sep = "\t")
# Folded distribution of assymmetry parameter
# i.e. p = 1-p if p>0.5, to get values between 0 and 0.5 (proportion of the shorter arm)
parameters_GammaModel$p5[which(parameters_GammaModel$p5 > 0.5)] = 1 - parameters_GammaModel$p5[which(parameters_GammaModel$p5 > 0.5)]

parameters_GammaModel$set = as.character(parameters_GammaModel$set)
parameters_GammaModel$map = as.character(parameters_GammaModel$map)
parameters_GammaModel$id = NA
for (i in 1:nrow(parameters_GammaModel)) {
  parameters_GammaModel$id[i] = paste(unlist(strsplit(parameters_GammaModel$set[i], split = "_"))[1],
                                      unlist(strsplit(parameters_GammaModel$set[i], split = "_"))[2],
                                      parameters_GammaModel$map[i],
                                      sep = "_")
}

karyotypes$id = paste(karyotypes$Genus, karyotypes$Species, karyotypes$Chromosome, sep = "_")

df = merge(parameters_GammaModel, karyotypes, by = "id")

# Compare differences in a Q-Q plot and density
plot(df$Centromeric_Index, df$p5, xlim = c(0, 0.5), ylim = c(0, 0.5),
     xlab = "Centromeric Index", ylab = "Parameter p")
abline(a = 0, b = 1)
df$diff = df$Centromeric_Index - df$p5
plot(density(df$diff, na.rm = TRUE))
# Not working very well, at least...

#============================================================================#




























































#---------------------------#
# Fitting a distribution on Marey plots
#---------------------------#

# Beta CDF fitted to the Marey plot with scaled values to {0,1} for both distances
# In order to compare all maps of different length, they were scaled to {0, 1}

# The marey plot is a CDF 
# Since data is defined  on the interval {0,1} and the PDF looks like a U-shape distribution,
# we can assume a Beta distribution with two parameters

# Likelihood can not be computed on a CDF (only PDF have likelihood functions),
# hence we used a Least Squares optimization method with optim()

# Yet, the Beta distribution did not fit well to the Marey plot
# We could have better results with a Weighted Least Squares optimization

# Alternatively, the Beta parameters could be fitted directly on the recombination rate estimated from curve fitting
# The recombination rate is the derivative of the Marey plot curve (i.e. the PDF)
# It would provide a cross-validation for parameters estimates

# Three test dataset
X = df$phys[df$group == "Arabidopsis thaliana 2"] 
Y = df$gen[df$group == "Arabidopsis thaliana 2"]

X = df$phys[df$group == "Triticum aestivum 1B"]
Y = df$gen[df$group == "Triticum aestivum 1B"]

X = df$phys[df$group == "Malus domestica 12"] 
Y = df$gen[df$group == "Malus domestica 12"]

library(extraDistr)
# removing outliers 0 values
# Y[Y == 0] = 0.0001
# Y[Y == 1] = 0.9999
distcdf = function(x) {
  shape1 = x[1]
  shape2 = x[2]
  ncp = x[3]

  fitted = pbeta(X, shape1, shape2)
  # Noncentral beta distribution
  # fitted = pbeta(X, shape1, shape2, ncp)
  # Dirichlet
  # fitted = pdirichlet(X, shape1)
  sum((Y-fitted)^2, na.rm = TRUE) # Criterion: sum of squared differences
  
  # mod = lm(Y~pbeta(X, shape1, shape2, ncp)) # Criterion: sum of squared differences
  # plot(X, abs(residuals(mod)))
  # Weights are a function of absolute values of residuals 
  # w = abs(residuals(mod))
  # Minimize the weight of values with lower residual
  # sum(((Y-fitted) * (1+w))^2, na.rm = TRUE) # Criterion: sum of squared differences
  
  # Minimize SSR Sum of Squared Residuals
  # sum(mod$residuals^2)
  
  # try some other criterion
  # https://stats.stackexchange.com/questions/122708/fitting-parametric-cdf-to-ecdf
  # minimizing the Anderson-Darling statistic (precision-weighted MSE)
  # as.numeric(ad.test(Y, "pbeta", X, shape1, shape2)$statistic)
  # minimizing the Kolmogorov Smirnoff statistic
  # as.numeric(ks.test(Y, pbeta(X, shape1, shape2))$statistic)
  
  
   # sum(abs(Y-pbeta(X, shape1, shape2, ncp))) # Criterion: sum of absolute values differences
}

# https://stats.stackexchange.com/questions/80985/how-do-i-assign-weight-to-each-observation-of-the-data-in-weighted-least-square

# library(logitnorm)
# 
# distcdf = function(x) {
#   shape1 = x[1]
#   shape2 = x[2]
#   
#   # Noncentral beta distribution
#   fitted = plogitnorm(X, shape1, shape2)
#   sum((Y-fitted)^2, na.rm = TRUE) # Criterion: sum of squared differences
# }

# fit.nlm = nlm(distcdf, p=c(1,1,1))$estimate
# # fit.nlm
# fit.optim = optim(par = c(0.1,0.1), distcdf, control = list(trace = 1, maxit = 1000))$par
# fit.optim = optim(par = c(1,1,1), distcdf)$par
fit.optim = optim(par = c(0.5, 0.5), distcdf, method = "L-BFGS-B", lower = 0.01, control = list(trace = 1, maxit = 10000))$par
fit.optim = optim(par = c(0.5, 0.5, 1), distcdf, method = "L-BFGS-B", lower = 0.01, control = list(trace = 1, maxit = 10000))$par

# Same parameter values found for different strating parameters, hence assessing a global maximum
# Avoiding local maxima

# TO DO iterative selection of the best parameter (1,000 iterations) by using random start parameters

# Parameters:
fit.optim

# Plot of the observed ECDF with the fitted beta CDF (parameters fitted on the PDF)
plot(X,Y, xlim = c(0, 1), ylim = c(0, 1))
x = seq(0, 1, by = 0.01)
curve(pbeta(x, fit.optim[1],  fit.optim[2]), add=T, col = "Green")
curve(pbeta(x, fit.optim[1],  fit.optim[2], ncp = fit.optim[3]), add=T, col = "Blue")
# curve(pbeta(x, shape1 = fit.optim[1], shape2 = fit.optim[2], ncp = fit.optim[3]), add=T, col = "Red")

# Comparing means to assess the fit
mean(Y)
mean(pbeta(x, shape1 = fit.optim[1], shape2 = fit.optim[2], ncp = fit.optim[3]))




#---------------------------#
# Fitting a distribution on recombination maps
#---------------------------#








































#---------------------------#
# Various tests of implementation, none working better than above...
#---------------------------#


# Fitting quantiles





# Fitting a polynom to the monotonic function (CDF) to compute the derivative in bins of 0.1% of the total length
# Then, fit a Beta to the derivative (i.e. the PDF)



# Maximum goodness of fit
mge.fit = mgedist(Y, "beta", gof = "CvM", start = list(shape1 = 0.5, shape2 = 0.5))$estimate
plot(X,Y)
x = seq(0.01, 0.99, by = 0.01)
# est = fit.nlm
# curve(pbeta(q = x, shape1 = est[1], shape2 = est[2]), add=T, col = "Blue")
est = mge.fit
curve(pbeta(x, shape1 = est[1], shape2 = est[2]), add=T, col = "Red")
est = fit.optim
curve(pbeta(x, shape1 = est[1], shape2 = est[2]), add=T, col = "Blue")
abline(0, 1)







# Setting phys. distances points as the center of each interval
dX = rowMeans(embed(X,2)) # centers the X values for plotting

# Computing derivative of the genetic distances for the interval in these points -> the PDF
dY = diff(Y)/diff(X) # the derivative of Y

# Fitting the beta distribution on the PDF
hist(dY)
# removing outliers and 0 derivatives
dX = dX[dY != 0]
dY = dY[dY != 0]
dX = dX[dY < 40]
dY = dY[dY < 40]

plot(dX, dY, ylim = c(0,100))

fit = fitdistr(dY, "beta", start = list(shape1 = 0.5, shape2 = 0.5))

# Plot of the observed ECDF with the fitted beta CDF (parameters fitted on the PDF)
x = seq(0, 1, by = 0.01)
pdf = dbeta(x, shape1 = fit$estimate[1], shape2 = fit$estimate[2])
pdf = dbeta(x, shape1 = 0.5, shape2 = 0.5)
lines(x, pdf)

X = X[Y!=0]
Y = Y[Y!=0]
X = X[Y!=1]
Y = Y[Y!=1]

# fit = fitdistr(as.numeric(df$gen+0.001), densfun = "lnorm")
tmp = df$gen[df$group == "Arabidopsis thaliana 2"]
?pbeta
CDF = df$gen[df$group == "Arabidopsis thaliana 2"]
breaks = df$phys[df$group == "Arabidopsis thaliana 2"]
# Now generate points on the cdf
x = seq(0, 1, by = 0.01)
# The Probability Density function (PDF) is the derivative of CDF.
PDF = diff(CDF)/diff(breaks)
hist(PDF)
descdist(PDF, discrete = FALSE)



plot(Y, X)
fit = fitdist(as.numeric(Y), method = "mme", "beta")
# cdf = pgamma(x, shape = 2.2840497, rate = 4.3303002)
# Now calculate points on the cdf
x = seq(0, 1, by = 0.01)
cdf = dbeta(x, shape1 = fit$estimate[1], shape2 = fit$estimate[2])
# Shown plotted here
lines(x,cdf)
# Now calculate points on the cdf
# cdf = pgamma(x, shape = 1.83317062, rate = 3.53284703)
# Shown plotted here
# lines(x,cdf)

# library(MASS)
# X = X[-171]
# Y = Y[-171]
# gammafit <- fitdistr(X, "beta", start = list(shape1 = 1, shape2 = 2))   
# #    shape        rate   
# #    17.552961    2.902459 
# #    ( 5.366214) ( 0.900112)
# plot(X, Y, col = "red")
# lines(X, pbeta(X, gammafit$estimate[1], gammafit$estimate[2]))



# Custom fit for CDF
# customFit <- function(x, data) {
#   d.data <- rev(cumsum(dnorm(1:length(data), x[1], x[2]))) * max(data)
#   SS <- sum((d.data - data)^2)
#   return(SS)
# }
# 
# fit.optim <- optim(c(5, 8), customFit, data = boot.mean)
# 
# plot(boot.mean)
# lines(rev(cumsum(dnorm(1:length(boot.mean), 
#                        fit.optim$par[1], fit.optim$par[2]))) * max(boot.mean), 
#       col = "red")


# Second method: fit with nlm (non-linear minimization) and squared differences
X = df$phys[df$group == "Arabidopsis thaliana 2"]
Y = df$gen[df$group == "Arabidopsis thaliana 2"]

X = df$phys[df$group == "Triticum aestivum 1B"]
Y = df$gen[df$group == "Triticum aestivum 1B"]
X = X[!(Y == 0)]
Y = Y[!(Y == 0)]
X = X[!(Y == 1)]
Y = Y[!(Y == 1)]

logit = function(x) {
  log(x/(1-x))
}

distcdf = function(x) {
  shape1 = x[1]
  shape2 = x[2]
  # mod = lm(Y~pbeta(X, 0.5, 0.5)) # Criterion: sum of squared differences
  # plot(fitted(mod),residuals(mod))
  # w = 
  fitted = pbeta(X, shape1, shape2)
  sum((Y-fitted)^2, na.rm = TRUE) # Criterion: sum of squared differences
  # sum(abs(Y-pbeta(X, shape1, shape2))) # Criterion: sum of absolute values differences
}

# https://stats.stackexchange.com/questions/80985/how-do-i-assign-weight-to-each-observation-of-the-data-in-weighted-least-square

# fit.nlm = nlm(distcdf, c(1,1))$estimate
# fit.optim = optim(par = c(0.1,0.1), distcdf, control = list(trace = 1, maxit = 1000))$par
# fit.optim = optim(par = c(1,1), distcdf)$par
# fit.optim = optim(par = c(0.5,0.5), distcdf)$par
fit.optim = optim(par = c(0.3,0.2), distcdf)$par
# Same parameter values found for different strating parameters, hence assessing a global maximum
# Avoiding local maxima

# fit.nlm
fit.optim

plot(X,Y)
x = seq(0.01, 0.99, by = 0.01)
# est = fit.nlm
# curve(pbeta(q = x, shape1 = est[1], shape2 = est[2]), add=T, col = "Blue")
est = fit.optim
curve(pbeta(q = x, shape1 = est[1], shape2 = est[2]), add=T, col = "Red")






# function needed for visualization purposes
# Fitting a logistic for sigmoid curves
sigmoid = function(params, x) {
  params[1] / (1 + exp(-params[2] * (x - params[3])))
}

# x = 1:53
# y = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0.1,0.18,0.18,0.18,0.33,0.33,0.33,0.33,0.41,
#       0.41,0.41,0.41,0.41,0.41,0.5,0.5,0.5,0.5,0.68,0.58,0.58,0.68,0.83,0.83,0.83,
#       0.74,0.74,0.74,0.83,0.83,0.9,0.9,0.9,1,1,1,1,1,1,1)

distcdf = function(x) {
  a = x[1]
  b = x[2]
  c = x[3]
  sum((Y - (a/(1+exp(-b*(x-c)))))^2) # Criterion: sum of squared differences
}
# fitting code
fitmodel = nlm(distcdf, c(1,1,1))$estimate

# fitmodel <- nls(y~a/(1 + exp(-b * (x-c))), start=list(a=1,b=.5,c=25))
# fitmodel <- nls(Y~a/(1 + exp(-b * (X-c))), start=list(a=0.5,b=.5,c=2))

# visualization code
# get the coefficients using the coef function
plot(X,Y)

params=c(a=fitmodel[1], b=fitmodel[2], c=fitmodel[3])
params
# params=c(a=0.5,b=.5,c=2)
x = seq(0.01, 0.99, by = 0.01)
p <- sigmoid(params, x)
p = params[1]/(1 + exp(-params[2] * (x - params[3])))
lines(p~x)




# https://stats.stackexchange.com/questions/112614/determining-beta-distribution-parameters-alpha-and-beta-from-two-arbitrary
#

x = df$phys[df$group == "Triticum aestivum 1B"]
y = df$gen[df$group == "Triticum aestivum 1B"]
# Logistic transformation of the fitted Beta CDF.
#
f.beta <- function(alpha, beta, x) {
  p <- pbeta(x, alpha, beta)
  log(p/(1-p)) # Transform into logit log(1/(1-p))
}
#
# Logistic transformation of the Beta empirical CDF.
#
y.p = log(y/(1-y))
#
# Sums of squares.
#
delta <- function(fit, actual) sum((fit-actual)^2)
#
# The objective function handles the transformed parameters `theta` and
# uses `f.beta` and `delta` to fit the values and measure their discrepancies.
#
objective <- function(theta, x, obs, ...) {
  ab <- exp(theta) # Parameters are the *logs* of alpha and beta
  fit <- f.beta(ab[1], ab[2], x, ...)
  return (delta(fit, obs))
}
#
# Solve two problems.
#
# x.p <- f.beta(alpha, beta, x)        # The correct values of the p_i
start <- log(c(0.5, 0.5))
sol <- nlm(objective, start, x=x, obs = y)
parms <- exp(sol$estimate)           # Estimates of alpha and beta
parms
plot(y~x, type="l", lwd=2)
curve(pbeta(x, parms[1], parms[2]), add=TRUE, col="Red")




# Four steps:
# 1) Model/function choice: hypothesize families of distributions;
# 2) Estimate parameters;
# 3) Evaluate quality of fit;
# 4) Goodness of fit statistical tests.

# Two methods

# Least Squares Estimation  (LSE)
# Less efficient and accurate, but works with CDF
# https://fr.mathworks.com/help/stats/examples/fitting-a-univariate-distribution-using-cumulative-probabilities.html

# If variance is not constant, Weighted Least Squares (https://en.wikipedia.org/wiki/Weighted_least_squares)
# Besides, some nonlinear regression problems can be moved to a linear domain by a suitable transformation of the model formulation. 

# Maximum Likelihood Estimation (MSE)
# More efficient but requires the PDF
# For more information on maximum likelihood parameter estimation, read @myungTutorialMaximumLikelihood2003
# MLE must be performed on PDF and not CDF
# CDF does not allow to asses the peak of likelihood: https://stats.stackexchange.com/questions/322375/why-do-we-need-density-in-estimation-and-cumulative-distribution-in-transformati
# https://stats.stackexchange.com/questions/220783/is-it-better-to-look-at-ecdfs-than-pdfs-when-exploring-empirical-sample-distribu?rq=1
# https://stats.stackexchange.com/questions/327625/why-maximum-likelihood-estimation-use-the-product-of-pdfs-rather-than-cdfs?rq=1


















df = read.table(file = "data-cleaned/Marey_maps/Marey_plot.norm.txt", header = TRUE, sep = "\t")
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
  tmp = as.numeric(df$gen[df$species == sp & df$map == chr])
  tmp = tmp[!(tmp == 0)]
  tmp = tmp[!(tmp == 1)]
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





df = read.table(file = "data-cleaned/Marey_maps/Marey_plot.norm.txt", header = TRUE, sep = "\t")
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



# https://stats.stackexchange.com/questions/206073/looking-for-function-to-fit-sigmoid-like-curve