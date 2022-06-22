############################################################################
#                    ECOBIO - PhD
#
#       Creating recombination landscapes with Marey Map
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
library(myTAI)

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


set = "Coffea_arabica_Crouzillat2020.txt"

#============================================================================
# Map quality control & Data cleaning
#============================================================================

#-------------------------------------
# Diagnostic plots
# Genetic distance (cM) as a function of physical position (bp)
marey.diagplot(set = set, dir = "data/Marey_maps/", output = "output/marey_maps/diagnostic_plots/", display = TRUE)


#-------------------------------------
# Check & correct bad orientation of linkage groups
# Markers oriented in the direction of a reference map if necessary
# A reference genome with markers mapped onto it seems inevitable

# FULL MAP CORRECTION
# Apply correction (automatic check not implemented yet)
# Set = name of a map
# id_chr = index of chromosomes to flip
new_map = check_map_orientation(set, id_chr = c("I_E"))
# Plot the data to validate the procedure
marey.diagplot(set = new_map)
# Save the corrected map once validated
write.table(new_map, file = paste("data/Marey_maps/", set, sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
save.log(msg = paste("Check and correct bad orientation of linkage groups: Map saved.", sep = " "))


#-------------------------------------
# Filter out outliers automatically
# Selection of valid markers

######### STEP 1 ######### 

# A rough filtering of clear outliers, i.e. points clearly aberrant
#-------------------------------------
# Filter out outliers manually
# Visually check for outliers in the Marey map plot
source("sources/MareyMap.R")
# Plot the data
marey.diagplot(set = set, dir = "data/Marey_maps/", save.plot = FALSE, display = TRUE)
# Apply correction (automatic check not implemented yet)
# Set = name of a map
# id_chr = index of chromosomes to flip
chr = "K_R"
new_map = outliers.selection(set , chr)
# Plot the data to validate the procedure
marey.diagplot(set = new_map, save.plot = FALSE, display = TRUE)
# Save the corrected map once validated
write.table(new_map, file = paste("data/Marey_maps/", set, sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
save.log(msg = paste("Manually removed some outliers in a chromosome: Map saved.", sep = " "))

# Discard a complete chromosome
chr = "H_E"
new_map = read.table(file = paste("data/Marey_maps/", set, sep = ""), header = TRUE)
new_map$vld[new_map$map == chr] = FALSE
write.table(new_map, file = paste("data/Marey_maps/", set, sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
save.log(msg = paste("Chromosome", chr,"in", set,"discarded: Map saved.", sep = " "))


# Reprint the MArey maps with outliers in red
marey.diagplot(set = set, save.plot = TRUE, display = TRUE)


#============================================================================
# Estimate the local recombination rate
#============================================================================

# Estimate the recombination rate at many positions
# (for instance all the genes of a genome)
# This can be done by up-loading a text file including all the positions.
# a text file (txt extension) containing at least a “map” column and a “phys” column
# indicating respectively the map and the physical position of each gene.

#----------------------------------------------------------------------
# LOESS method
#----------------------------------------------------------------------
# ?loess # As in MareyMap

# Par ailleurs, dans la mesure où elle repose sur des régressions par les moindres carrés,
# la régression locale bénéficie aussi de la plupart des outils liés à ces méthodes de régression, notamment la théorie de calcul des incertitudes de prédiction et de calibrage.
# Beaucoup d'autres tests et procédures utilisés pour valider les modèles par les moindres carrés peuvent également être étendus aux modèles de régression locale. 
# https://fr.wikipedia.org/wiki/R%C3%A9gression_locale
# https://en.wikipedia.org/wiki/Local_regression
#source("sources/MareyMap.R")

# A modified version of recombination.map, with ylim = c(0,100), to take into account high recombination rates
recombination.map.coffea = function(set ="", data = data.final, chr = "all", method = c("loess", "smooth.spline"), K = 5, boot = 1000,...) { # Optional arguments passed to loess, smooth.spline and pther fitting methods
  # Physical distances are in bp, yet analyses are in Mb
  # Convert bp to Mb
  scale = 1000000
  
  if (method == "veller") {
    boot = 0
  }
  set.name = set
  if (chr == "all") { # Chromosomes to map are given in a list of names of "all" to evaluate every chromosome in the dataset
    chr.list = unique(data$map[data$set == set.name])
  } else {
    chr.list = chr
  }
  # For each chosen chromosome in a dataset
  for (c in chr.list) {
    cat("Treating chromosome ", c, "...\n")
    # Regression is made on a whole chromosome
    chr.nb = c # Chromosome number
    data$set = as.character(data$set)
    data$map = as.character(data$map)
    chr.data = subset(data, map == chr.nb & set == set.name) # Sample only distances of the chromosome
    chr.data$gen = as.numeric(as.character(chr.data$gen))
    chr.data$phys = as.numeric(as.character(chr.data$phys))
    # Convert physical distances in Mb
    chr.data$phys = chr.data$phys/scale
    
    # Remove NA distances
    chr.data = chr.data[!is.na(chr.data$gen),]
    chr.data = chr.data[!is.na(chr.data$phys),]
    
    # span = 0.1
    degree = 2
    # library(spatialEco)
    # fit.loess.bis = loess.ci(chr.data$gen, chr.data$phys, p = 0.95, plot = TRUE, span = 0.1, degree = 2) # Another implementation giving same results
    
    #-----------------------------------------
    # NO BOOTSTRAPS  
    if (boot == 0) {
      pointwise = sort(chr.data$phys)
      # warning("Confidence intervals not estimated (no bootstrap).")
      # Bin recombination map
      # Partitioning of the physical map (bp): 1,000 bins of equal size
      phys.bins = seq(0, max(chr.data$phys, na.rm = TRUE), len = 1000) # A sequence of bins coordinates
      
      # 1Mb windows recombination map
      # Partitioning of the physical map in windows of 1Mb for further comparison with genomic content
      phys.1Mbwind = c(seq(0, max(chr.data$phys, na.rm = TRUE), by = 1), max(chr.data$phys, na.rm = TRUE)) # A sequence of windows coordinates
      
      # 100kb windows recombination map
      # Partitioning of the physical map in windows of 100kb for further comparison with genomic content
      phys.100kbwind = c(seq(0, max(chr.data$phys, na.rm = TRUE), by = 0.1), max(chr.data$phys, na.rm = TRUE)) # A sequence of windows coordinates
      
      #-----------------------------------------
      # LOESS
      #-----------------------------------------
      if (method == "loess") {
        #---------------------------------
        # Calibrate span
        span = calibrate.span(x = chr.data, K, from = 0.3, to = 0.5)
        
        #---------------------------------
        # Fit the interpolation curve on the data
        fit = fit.loess(x = chr.data, span = span, degree = degree)
        # Save the recombination map
        # Pointwise estimates (same as physical positions observed)
        predicted.pointwise = predict(fit, newdata = pointwise) # Predict the local recombination rate at the exact physical positions of the Marey map
        
        # Bin recombination map
        # Partitioning of the physical map (bp): 1,000 bins of equal size
        predicted.physbins = predict(fit, newdata = phys.bins) # Predict the local recombination rate at the exact physical positions of the Marey map
        
        # 1Mb windows recombination map
        # Partitioning of the physical map in windows of 1Mb for further comparison with genomic content
        predicted.phys.1Mbwind = predict(fit, newdata = phys.1Mbwind) # Predict the local recombination rate at the exact physical positions of the Marey map
        # 100kb windows recombination map
        predicted.phys.100kbwind = predict(fit, newdata = phys.100kbwind) # Predict the local recombination rate at the exact physical positions of the Marey map
        
        rec.pointwise = diff(predicted.pointwise)/diff(chr.data$phys)
        rec.bins = diff(predicted.physbins)/diff(phys.bins)
        rec.1Mbwind = diff(predicted.phys.1Mbwind)/diff(phys.1Mbwind)
        rec.100kbwind = diff(predicted.phys.100kbwind)/diff(phys.100kbwind)
        
      } else {
        #-----------------------------------------
        # CUBIC SMOOTH SPLINE
        #-----------------------------------------
        if (method == "smooth.spline") {
          #---------------------------------
          # Calibrate span
          spar = calibrate.spar(data = chr.data, K)
          
          #---------------------------------
          # Fit the interpolation curve on the data
          fit = fit.smooth.spline(x = chr.data, spar = spar)
          
          # Pointwise estimates (same as physical positions observed)
          predicted.pointwise = predict(fit, x = pointwise)$y # Predict the local recombination rate at the exact physical positions of the Marey map
          
          # Bin recombination map
          # Partitioning of the physical map (bp): 1,000 bins of equal size
          predicted.physbins = predict(fit, x = phys.bins)$y # Predict the local recombination rate at the exact physical positions of the Marey map
          
          # 1Mb windows recombination map
          # Partitioning of the physical map in windows of 1Mb for further comparison with genomic content
          predicted.phys.1Mbwind = predict(fit, x = phys.1Mbwind)$y # Predict the local recombination rate at the exact physical positions of the Marey map
          # 100kb windows recombination map
          predicted.phys.100kbwind = predict(fit, x = phys.100kbwind)$y # Predict the local recombination rate at the exact physical positions of the Marey map
          
          rec.pointwise = diff(predicted.pointwise)/diff(chr.data$phys)
          rec.bins = diff(predicted.physbins)/diff(phys.bins)
          rec.1Mbwind = diff(predicted.phys.1Mbwind)/diff(phys.1Mbwind)
          rec.100kbwind = diff(predicted.phys.100kbwind)/diff(phys.100kbwind)
          
        }  # END OF SMOOTH SPLINE PROCEDURE
        else {
          #---------------------------------
          # Veller method: interpolate 1,000 evenly distributed physical positions
          #---------------------------------
          if (method == "veller") {
            #---------------------------------
            # Calibrate span
            span = calibrate.span(x = chr.data, K)
            # span = 0.3
            #---------------------------------
            # Fit the interpolation curve on data
            fit = fit.loess(x = chr.data, span = span, degree = degree)
            # Partitioning of the physical map (bp): 1,000 bins of equal size
            predicted.physbins = predict(fit, newdata = phys.bins) # Predict the genetic distance at the exact physical positions of the Marey map
            
          } else {
            warning("Method must be loess, smooth.spline or veller.")
          }
        } # END OF VELLER METHOD
      }
      
      if (method == "loess" | method == "smooth.spline") {
        #-----------------------------------------
        # DERIVATIVES
        # The local recombination rate is the derivative of the polynom (i.e. fitted loess) in a point
        #------------------------------
        # Pointwise recombination map
        dY = rec.pointwise # the derivative of LOESS function
        dX = rowMeans(embed(pointwise, 2)) # centers the X values for plotting
        # Return a data.frame of four columns for each method -> final plotting df
        df.pointwise = data.frame(phys = dX, rec.rate = dY)
        write.table(df.pointwise, file = paste("output/recombination_maps/", method,"/pointwise/", set.name, "_chromosome", chr.nb, ".txt", sep = ""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
        
        # Save graphical predictions
        png(paste("output/recombination_maps_fig/", method, "/pointwise/", set.name,"_chr", chr.nb,".png",sep = ""), width = 1200, height = 600)
        plot(dX, dY, type="l", ylim = c(0, 30), main = paste("Recombination landscape (", method, ") of chromosome ", chr.nb, " in ", set.name, sep = ""), xlab = "Physical position (Mb)", ylab = "Recombination rate (cM/Mb)") #check
        abline(h = 0, col = "Red", lty = 2) # Print a threshold for rec. rate = 0
        dev.off()
        rm(dX, dY)
        rm(rec.pointwise)
        
        #------------------------------
        # Bin recombination map
        dY = rec.bins # the derivative of LOESS function
        dX = rowMeans(embed(phys.bins, 2)) # centers the X values for plotting
        # Return a data.frame of four columns for each method -> final plotting df
        df.bins = data.frame(phys = dX, rec.rate = dY)
        write.table(df.bins, file = paste("output/recombination_maps/", method,"/bins/", set.name, "_chromosome", chr.nb, ".txt", sep = ""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
        
        # Save graphical predictions
        png(paste("output/recombination_maps_fig/", method, "/bins/", set.name,"_chr", chr.nb,".png",sep = ""), width = 1200, height = 600)
        plot(dX, dY, type="l", ylim = c(0, 30), main = paste("Recombination landscape (", method, ") of chromosome ", chr.nb, " in ", set.name, sep = ""), xlab = "Physical position (Mb)", ylab = "Recombination rate (cM/Mb)") #check
        abline(h = 0, col = "Red", lty = 2) # Print a threshold for rec. rate = 0
        dev.off()
        rm(dX, dY)
        rm(phys.bins)
        rm(rec.bins)
        
        #------------------------------
        # 1Mb windows recombination map
        dY = rec.1Mbwind # the derivative of LOESS function
        dX = phys.1Mbwind = rowMeans(embed(phys.1Mbwind,2)) # Take the center of each window for interpolation
        # Return a data.frame of four columns for each method -> final plotting df
        df.1Mbwind = data.frame(phys = dX, rec.rate = dY)
        write.table(df.1Mbwind, file = paste("output/recombination_maps/", method,"/1Mbwind/", set.name, "_chromosome", chr.nb, ".txt", sep = ""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
        
        # Save graphical predictions
        png(paste("output/recombination_maps_fig/", method, "/1Mbwind/", set.name,"_chr", chr.nb,".png",sep = ""), width = 1200, height = 600)
        plot(dX, dY, type="l", ylim = c(0, 30), main = paste("Recombination landscape (", method, ") of chromosome ", chr.nb, " in ", set.name, sep = ""), xlab = "Physical position (Mb)", ylab = "Recombination rate (cM/Mb)") #check
        abline(h = 0, col = "Red", lty = 2) # Print a threshold for rec. rate = 0
        dev.off()
        rm(dX, dY)
        rm(phys.1Mbwind)
        rm(rec.1Mbwind)
        
        #------------------------------
        # 100kb windows recombination map
        dY = rec.100kbwind # the derivative of LOESS function
        dX = rowMeans(embed(phys.100kbwind,2)) # Take the center of each window for interpolation
        # Return a data.frame of four columns for each method -> final plotting df
        df.100kbwind = data.frame(phys = dX, rec.rate = dY)
        write.table(df.100kbwind, file = paste("output/recombination_maps/", method,"/100kbwind/", set.name, "_chromosome", chr.nb, ".txt", sep = ""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
        
        # Save graphical predictions
        png(paste("output/recombination_maps_fig/", method, "/100kbwind/", set.name,"_chr", chr.nb,".png",sep = ""), width = 1200, height = 600)
        plot(dX, dY, type="l", ylim = c(0, 30), main = paste("Recombination landscape (", method, ") of chromosome ", chr.nb, " in ", set.name, sep = ""), xlab = "Physical position (Mb)", ylab = "Recombination rate (cM/Mb)") #check
        abline(h = 0, col = "Red", lty = 2) # Print a threshold for rec. rate = 0
        dev.off()
        rm(dX, dY)
        rm(phys.100kbwind)
        rm(rec.100kbwind)
      } else {
        if (method == "veller") {
          # A Marey map of 1,000 evenly spaced markers
          df.tmp = read.table("data-cleaned/veller/GeneticDistances.txt", header = TRUE)
          # Return a data.frame
          df.veller = data.frame(set = rep(set.name, 1000), map = rep(chr.nb, 1000), mkr = paste("chr", chr.nb,"_",1:1000, sep = ""), phys = phys.bins, gen = predicted.physbins)
          # Append new data to replace
          df.tmp = df.tmp[which(!(df.tmp$set == set.name & df.tmp$map == chr.nb)),]
          df.tmp = rbind(df.tmp, df.veller)
          write.table(df.tmp, file = "data-cleaned/veller/GeneticDistances.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
          
        }
      }
      
      # END OF NO BOOTSTRAP PROCEDURE
    } else {
      if (boot > 0 & is.numeric(boot)) {
        pointwise = sort(chr.data$phys)
        pointwise.boot = matrix(NA, nrow = length(chr.data$phys)-1, ncol = boot) # Predict the local recombination rate at the exact physical positions of the Marey map
        # Bin recombination map
        # Partitioning of the physical map (bp): 1,000 bins of equal size
        phys.bins = seq(0, max(chr.data$phys, na.rm = TRUE), len = 1000) # A sequence of bins coordinates
        physbins.boot = matrix(NA, nrow = length(phys.bins)-1, ncol = boot) # Predict the local recombination rate at the exact physical positions of the Marey map
        
        # 1Mb windows recombination map
        # Partitioning of the physical map in windows of 1Mb for further comparison with genomic content
        phys.1Mbwind = c(seq(0, max(chr.data$phys, na.rm = TRUE), by = 1), max(chr.data$phys, na.rm = TRUE)) # A sequence of windows coordinates
        # phys.1Mbwind = rowMeans(embed(phys.1Mbwind,2)) # Take the center of each window for interpolation
        phys.1Mbwind.boot = matrix(NA, nrow = length(phys.1Mbwind)-1, ncol = boot)
        
        # 100kb windows recombination map
        # Partitioning of the physical map in windows of 100kb for further comparison with genomic content
        phys.100kbwind = c(seq(0, max(chr.data$phys, na.rm = TRUE), by = 0.1), max(chr.data$phys, na.rm = TRUE)) # A sequence of windows coordinates
        # phys.100kbwind = rowMeans(embed(phys.100kbwind,2)) # Take the center of each window for interpolation
        phys.100kbwind.boot = matrix(NA, nrow = length(phys.100kbwind)-1, ncol = boot)
        
        #-----------------------------------------
        # LOESS
        #-----------------------------------------
        if (method == "loess") {
          #---------------------------------
          # Calibrate span
          span = calibrate.span(x = chr.data, K)
          # if (is.na(span)) {span = 0.3}
          # span = 0.3
          #---------------------------------
          # Bootstrapping
          # Bootstrap computation can be very time demanding
          # TODO Code optimization of bootstrap, estimate 'spar' parameter only once, outside the loop
          pb = txtProgressBar(min = 1, max = boot, initial = 1)
          for (b in 1:boot) {
            # print(b)
            setTxtProgressBar(pb, b)
            # Resample data
            resampled = chr.data[sample(1:nrow(chr.data), replace = TRUE),]
            
            # Fit the interpolation curve on the resampled data
            fit = fit.loess(x = resampled, span = span, degree = degree)
            # Return a list list(pointwise = predicted.pointwise, physbins = predicted.physbins, phys.1Mbwind = predicted.phys.1Mbwind)
            # Store results for each map
            
            # if (b == 1) {
            #   # Init a matrix of temporary bootstraps results
            #   pointwise.boot = matrix(NA, nrow = length(predicted.pointwise), ncol = boot)
            #   physbins.boot = matrix(NA, nrow = length(predicted.physbinsmatrix), ncol = boot)
            #   phys.1Mbwind.boot = matrix(NA, nrow = length(predicted.phys.1Mbwind), ncol = boot)
            # }
            # Save the recombination map
            # Pointwise estimates (same as physical positions observed)
            predicted.pointwise = predict(fit, newdata = pointwise) # Predict the local recombination rate at the exact physical positions of the Marey map
            
            # Bin recombination map
            # Partitioning of the physical map (bp): 1,000 bins of equal size
            predicted.physbins = predict(fit, newdata = phys.bins) # Predict the local recombination rate at the exact physical positions of the Marey map
            
            # 1Mb windows recombination map
            # Partitioning of the physical map in windows of 1Mb for further comparison with genomic content
            predicted.phys.1Mbwind = predict(fit, newdata = phys.1Mbwind) # Predict the local recombination rate at the exact physical positions of the Marey map
            # 100kb windows recombination map
            predicted.phys.100kbwind = predict(fit, newdata = phys.100kbwind) # Predict the local recombination rate at the exact physical positions of the Marey map
            
            pointwise.boot[,b] = diff(predicted.pointwise)/diff(chr.data$phys)
            physbins.boot[,b] = diff(predicted.physbins)/diff(phys.bins)
            phys.1Mbwind.boot[,b] = diff(predicted.phys.1Mbwind)/diff(phys.1Mbwind)
            phys.100kbwind.boot[,b] = diff(predicted.phys.100kbwind)/diff(phys.100kbwind)
          }
          Sys.sleep(1)
          close(pb)
        } else {
          #-----------------------------------------
          # CUBIC SMOOTH SPLINE
          #-----------------------------------------
          if (method == "smooth.spline") {
            #---------------------------------
            # Calibrate span
            spar = calibrate.spar(data = chr.data, K)
            
            #---------------------------------
            # Bootstrapping
            # Bootstrap computation can be very time demanding
            pb = txtProgressBar(min = 1, max = boot, initial = 1)
            for (b in 1:boot) {
              # print(b)
              setTxtProgressBar(pb, b)
              # Resample data
              resampled = chr.data[sample(1:nrow(chr.data), replace = TRUE),]
              
              # Fit the interpolation curve on the resampled data
              fit = fit.smooth.spline(x = resampled, spar = spar)
              # Return a list list(pointwise = predicted.pointwise, physbins = predicted.physbins, phys.1Mbwind = predicted.phys.1Mbwind)
              # Store results for each map
              
              # if (b == 1) {
              #   # Init a matrix of temporary bootstraps results
              #   pointwise.boot = matrix(NA, nrow = length(predicted.pointwise), ncol = boot)
              #   physbins.boot = matrix(NA, nrow = length(predicted.physbinsmatrix), ncol = boot)
              #   phys.1Mbwind.boot = matrix(NA, nrow = length(predicted.phys.1Mbwind), ncol = boot)
              # }
              # Save the recombination map
              # Pointwise estimates (same as physical positions observed)
              predicted.pointwise = predict(fit, x = pointwise)$y # Predict the local recombination rate at the exact physical positions of the Marey map
              
              # Bin recombination map
              # Partitioning of the physical map (bp): 1,000 bins of equal size
              predicted.physbins = predict(fit, x = phys.bins)$y # Predict the local recombination rate at the exact physical positions of the Marey map
              
              # 1Mb windows recombination map
              # Partitioning of the physical map in windows of 1Mb for further comparison with genomic content
              predicted.phys.1Mbwind = predict(fit, x = phys.1Mbwind)$y # Predict the local recombination rate at the exact physical positions of the Marey map
              # 100kb windows recombination map
              predicted.phys.100kbwind = predict(fit, x = phys.100kbwind)$y # Predict the local recombination rate at the exact physical positions of the Marey map
              
              pointwise.boot[,b] = diff(predicted.pointwise)/diff(chr.data$phys)
              physbins.boot[,b] = diff(predicted.physbins)/diff(phys.bins)
              phys.1Mbwind.boot[,b] = diff(predicted.phys.1Mbwind)/diff(phys.1Mbwind)
              phys.100kbwind.boot[,b] = diff(predicted.phys.100kbwind)/diff(phys.100kbwind)
            }
            Sys.sleep(1)
            close(pb) 
            # # 1,000 markers in the Marey map
            # # Partitioning of the physical map (bp): 1,000 bins of equal size
            # phys.bins = seq(0, max(chr.data$phys, na.rm = TRUE), len = 1000) # A sequence of bins coordinates
            # 
            # # Fit the interpolation curve on the resampled data
            # fit = fit.smooth.spline(x = resampled, spar = spar)
            # # Partitioning of the physical map (bp): 1,000 bins of equal size
            # predicted.physbins = predict(fit, x = phys.bins)$y # Predict the local recombination rate at the exact physical positions of the Marey map
            
          } else {
            warning("Method must be loess or smooth.spline.")
          }
          # END OF SMOOTH SPLINE PROCEDURE
        }
        
        #-----------------------------------------
        # DERIVATIVES
        # The local recombination rate is the derivative of the polynom (i.e. fitted loess) in a point
        #------------------------------
        upperCI = function(x) {
          quantile(x, 0.975, na.rm = TRUE)
        }
        lowerCI = function(x) {
          quantile(x, 0.025, na.rm = TRUE)
        }
        # Pointwise recombination map
        dY = rowMeans(pointwise.boot, na.rm = TRUE) # the derivative of LOESS function
        # Compute the 95% CI of the derivative, apply on rows
        dY.upper = apply(pointwise.boot, MARGIN = 1, FUN = upperCI)
        dY.lower = apply(pointwise.boot, MARGIN = 1, FUN = lowerCI)
        dX = rowMeans(embed(pointwise, 2)) # centers the X values for plotting
        # dX = rowMeans(embed(chr.data$phys,2)) # centers the X values for plotting
        # Negative estimates are impossible
        # Negative values are forced to zero
        dY[dY < 0] = 0
        dY.upper[dY.upper < 0] = 0
        dY.lower[dY.lower < 0] = 0
        
        # Return a data.frame of four columns for each method -> final plotting df
        df.pointwise = data.frame(phys = dX, rec.rate = dY, upper = dY.upper, lower = dY.lower)
        
        write.table(df.pointwise, file = paste("output/recombination_maps/", method,"/pointwise/", set.name, "_chromosome", chr.nb, ".txt", sep = ""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
        
        # Save graphical predictions
        png(paste("output/recombination_maps_fig/", method, "/pointwise/", set.name,"_chr", chr.nb,".png",sep = ""), width = 1200, height = 600)
        plot(dX, dY, type="l", ylim = c(0, 60), main = paste("Recombination landscape (", method, ") of chromosome ", chr.nb, " in ", set.name, sep = ""), xlab = "Physical position (Mb)", ylab = "Recombination rate (cM/Mb)") #check
        lines(dX, dY.upper, type = "l", lty = 2, col = "Grey")
        lines(dX, dY.lower, type = "l", lty = 2, col = "Grey")
        abline(h = 0, col = "Red", lty = 2) # Print a threshold for rec. rate = 0
        dev.off()
        rm(dX, dY, dY.upper, dY.lower)
        rm(pointwise.boot)
        
        #------------------------------
        # Bin recombination map
        dY = rowMeans(physbins.boot, na.rm = TRUE) # the derivative of LOESS function
        # Compute the 95% CI of the derivative, apply on rows
        dY.upper = apply(physbins.boot, MARGIN = 1, FUN = upperCI)
        dY.lower = apply(physbins.boot, MARGIN = 1, FUN = lowerCI)
        dX = rowMeans(embed(phys.bins,2)) # centers the X values for plotting
        # Negative estimates are impossible
        # Negative values are forced to zero
        dY[dY < 0] = 0
        dY.upper[dY.upper < 0] = 0
        dY.lower[dY.lower < 0] = 0
        
        # Return a data.frame of four columns for each method -> final plotting df
        df.bins = data.frame(phys = dX, rec.rate = dY, upper = dY.upper, lower = dY.lower)
        
        write.table(df.bins, file = paste("output/recombination_maps/", method,"/bins/", set.name, "_chromosome", chr.nb, ".txt", sep = ""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
        
        # Save graphical predictions
        png(paste("output/recombination_maps_fig/", method, "/bins/", set.name,"_chr", chr.nb,".png",sep = ""), width = 1200, height = 600)
        plot(dX, dY, type="l", ylim = c(0, 100), main = paste("Recombination landscape (", method, ") of chromosome ", chr.nb, " in ", set.name, sep = ""), xlab = "Physical position (Mb)", ylab = "Recombination rate (cM/Mb)") #check
        lines(dX, dY.upper, type = "l", lty = 2, col = "Grey")
        lines(dX, dY.lower, type = "l", lty = 2, col = "Grey")
        abline(h = 0, col = "Red", lty = 2) # Print a threshold for rec. rate = 0
        dev.off()
        rm(dX, dY, dY.upper, dY.lower)
        rm(physbins.boot)
        
        #------------------------------
        # 1Mb windows recombination map
        dY = rowMeans(phys.1Mbwind.boot, na.rm = TRUE) # the derivative of LOESS function
        # Compute the 95% CI of the derivative, apply on rows
        dY.upper = apply(phys.1Mbwind.boot, MARGIN = 1, FUN = upperCI)
        dY.lower = apply(phys.1Mbwind.boot, MARGIN = 1, FUN = lowerCI)
        dX = rowMeans(embed(phys.1Mbwind,2)) # centers the X values for plotting
        # Negative estimates are impossible
        # Negative values are forced to zero
        dY[dY < 0] = 0
        dY.upper[dY.upper < 0] = 0
        dY.lower[dY.lower < 0] = 0
        
        # Return a data.frame of four columns for each method -> final plotting df
        df.1Mbwind = data.frame(phys = dX, rec.rate = dY, upper = dY.upper, lower = dY.lower)
        
        write.table(df.1Mbwind, file = paste("output/recombination_maps/", method,"/1Mbwind/", set.name, "_chromosome", chr.nb, ".txt", sep = ""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
        
        # Save graphical predictions
        png(paste("output/recombination_maps_fig/", method, "/1Mbwind/", set.name,"_chr", chr.nb,".png",sep = ""), width = 1200, height = 600)
        plot(dX, dY, type="l", ylim = c(0, 100), main = paste("Recombination landscape (", method, ") of chromosome ", chr.nb, " in ", set.name, sep = ""), xlab = "Physical position (Mb)", ylab = "Recombination rate (cM/Mb)") #check
        lines(dX, dY.upper, type = "l", lty = 2, col = "Grey")
        lines(dX, dY.lower, type = "l", lty = 2, col = "Grey")
        abline(h = 0, col = "Red", lty = 2) # Print a threshold for rec. rate = 0
        dev.off()
        rm(dX, dY, dY.upper, dY.lower)
        rm(phys.1Mbwind.boot)
        
        #------------------------------
        # 100kb windows recombination map
        dY = rowMeans(phys.100kbwind.boot, na.rm = TRUE) # the derivative of LOESS function
        # Compute the 95% CI of the derivative, apply on rows
        dY.upper = apply(phys.100kbwind.boot, MARGIN = 1, FUN = upperCI)
        dY.lower = apply(phys.100kbwind.boot, MARGIN = 1, FUN = lowerCI)
        dX = rowMeans(embed(phys.100kbwind,2)) # centers the X values for plotting
        # Negative estimates are impossible
        # Negative values are forced to zero
        dY[dY < 0] = 0
        dY.upper[dY.upper < 0] = 0
        dY.lower[dY.lower < 0] = 0
        
        # Return a data.frame of four columns for each method -> final plotting df
        df.100kbwind = data.frame(phys = dX, rec.rate = dY, upper = dY.upper, lower = dY.lower)
        
        write.table(df.100kbwind, file = paste("output/recombination_maps/", method,"/100kbwind/", set.name, "_chromosome", chr.nb, ".txt", sep = ""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
        
        # Save graphical predictions
        png(paste("output/recombination_maps_fig/", method, "/100kbwind/", set.name,"_chr", chr.nb,".png",sep = ""), width = 1200, height = 600)
        plot(dX, dY, type="l", ylim = c(0, 100), main = paste("Recombination landscape (", method, ") of chromosome ", chr.nb, " in ", set.name, sep = ""), xlab = "Physical position (Mb)", ylab = "Recombination rate (cM/Mb)") #check
        lines(dX, dY.upper, type = "l", lty = 2, col = "Grey")
        lines(dX, dY.lower, type = "l", lty = 2, col = "Grey")
        abline(h = 0, col = "Red", lty = 2) # Print a threshold for rec. rate = 0
        dev.off()
        rm(dX, dY, dY.upper, dY.lower)
        rm(phys.100kbwind.boot)
        # END OF BOOTSTRAP PROCEDURE
      } else {
        warning("Boot value must be 0 (no bootstrap) or positive integer (number of bootstrap).")
      }
    }
    
    
    # The span paramater is a smoothing parameter
    # the one parameter that need to be adjusted in order to avoid overfitting of the Marey map, larger values of span give smoother estimates
    # "Fitting is done locally. That is, for the fit at point x, the fit is made using points in a neighbourhood of x,
    # weighted by their distance from x (with differences in ‘parametric’ variables being ignored when computing the distance). 
    # The size of the neighbourhood is controlled by α (set by span or enp.target).
    # For α < 1, the neighbourhood includes proportion α of the points, and these have tricubic weighting (proportional to (1 - (dist/maxdist)^3)^3).
    # For α > 1, all points are used, with the ‘maximum distance’ assumed to be α^(1/p) times the actual maximum distance for p explanatory variables."
    
    # Compare predicted to observed in order to assess estimates robustness (RMSE?)
    # Can be used for model selection and/or parameter optimization
    
    #-----------------------------------------
    # PARAMETER OPTIMIZATION
    # Optim
    # define function that returns the SSE
    # calcSSE = function(x){
    #   loessMod = try(loess(gen ~ phys, data = chr.data, span=x, degree = 2), silent=T)
    #   res = try(loessMod$residuals, silent=T)
    #   if(class(res)!="try-error"){
    #     if((sum(res, na.rm=T) > 0)){
    #       sse = sum(res^2)  
    #     }
    #   }else{
    #     sse = 99999
    #   }
    #   return(sse)
    # }
    # # Run optim to find span that gives min criterion (SSE), starting at 0.1
    # optim(par=c(0.1, 0.5), calcSSE, method="SANN")
    # Yet, optimizing the SSE had little sense because it didn't prevented us from overfitting
    
    # New approach, based on data... cross-validation
    
    # Systematic K-fold cross-validation
    # 2 criterion: (see Lee & Cox 2010)
    # Least Square CV, sensitive to outliers
    # Absolute CV, less sensitive to gross outliers
    
    # Automated parameter selection
    # https://www.r-bloggers.com/automated-parameter-selection-for-loess-regression
    # https://cran.r-project.org/web/packages/fANCOVA/index.html
    # Leave-one-out robust cross validation
    
    # Integrated Squared Error (ISE), decribed in Lee & Cox 2010
    # Integrated Absolute Error (IAE)
  }
}

# For a chosen dataset
set.name = "Coffea_arabica_Crouzillat2020"
data = read.table(paste("data/Marey_maps/", set.name, ".txt", sep = ""), header = TRUE)
recombination.map.coffea(set = set.name, data = data, chr = "all", method = "loess", K = 3, boot = 1000)

set.name = "Coffea_canephora_Crouzillat2020"
data = read.table(paste("data-cleaned/marey_maps/", set.name, ".txt", sep = ""), header = TRUE)
recombination.map(set = set.name, data = data, chr = "all", method = "loess", K = 3, boot = 1000)

#----------------------------------------------------------------------
# Combined representation of recombination ladnscapes for Coffea robusta/Coffea arabica
#----------------------------------------------------------------------

# Extract R chromosomes in coffea arabica
# and get corresponding chromosomes in coffea robusta
set.name = "Coffea_arabica_Crouzillat2020"
data = read.table(paste("data/Marey_maps/", set.name, ".txt", sep = ""), header = TRUE)
commonChromosomes = gsub("_R", "",unique(data$map[which(data$vld == TRUE)])[grep("_R", unique(data$map[which(data$vld == TRUE)]))])

wd = "/Users/tbrazier/Documents/Coffea_Crouzillat/"
# Save graphical predictions
for (chr.nb in commonChromosomes) {
  set.name = "Coffea_arabica_Crouzillat2020"
  df.arabica = read.table(file = paste(wd, "recombination_maps/", set.name, "_chromosome", chr.nb, "_R.txt", sep = ""), sep = "\t", header = TRUE)
  set.name = "Coffea_canephora_Crouzillat2020"
  df.robusta = read.table(file = paste(wd, "recombination_maps/", set.name, "_chromosome", chr.nb, ".txt", sep = ""), sep = "\t", header = TRUE)
  
  png(paste(wd, "combined_recombination_landscapes/Combined_Coffea_chr", chr.nb,".png",sep = ""), width = 1200, height = 600)
  plot(df.arabica$phys, df.arabica$rec.rate, type="l", ylim = c(0, 100), xlim = c(0 , max(c(df.arabica$phys, df.robusta$phys), na.rm = TRUE)),
       main = paste("Coffea, chromosome ", chr.nb, sep = ""), xlab = "Physical position (Mb)", ylab = "Recombination rate (cM/Mb)") #check
  lines(df.arabica$phys, df.arabica$upper, type = "l", lty = 2, col = "Grey")
  lines(df.arabica$phys, df.arabica$lower, type = "l", lty = 2, col = "Grey")
  abline(h = 0, col = "Black", lty = 2) # Print a threshold for rec. rate = 0
  
  lines(df.robusta$phys, df.robusta$rec.rate, type="l", col = "Red")
  lines(df.robusta$phys, df.robusta$upper, type = "l", lty = 2, col = "Red")
  lines(df.robusta$phys, df.robusta$lower, type = "l", lty = 2, col = "Red")
  dev.off()
}
# Huge differences in chromosome sizes between robusta and arabica. Robusta is almost two times bigger than arabica.

# Scaled version on the physical distances (physical distances from 0 to 1)
for (chr.nb in commonChromosomes) {
  set.name = "Coffea_arabica_Crouzillat2020"
  df.arabica = read.table(file = paste(wd, "recombination_maps/", set.name, "_chromosome", chr.nb, "_R.txt", sep = ""), sep = "\t", header = TRUE)
  set.name = "Coffea_canephora_Crouzillat2020"
  df.robusta = read.table(file = paste(wd, "recombination_maps/", set.name, "_chromosome", chr.nb, ".txt", sep = ""), sep = "\t", header = TRUE)
  
  # Scaling physical distances
  df.arabica$phys = df.arabica$phys/max(df.arabica$phys, na.rm = TRUE)
  df.robusta$phys = df.robusta$phys/max(df.robusta$phys, na.rm = TRUE)
  
  png(paste(wd, "combined_recombination_landscapes_scaled/Combined_Coffea_chr", chr.nb,"_scaled.png",sep = ""), width = 1200, height = 600)
  plot(df.arabica$phys, df.arabica$rec.rate, type="l", ylim = c(0, 100), xlim = c(0 , max(c(df.arabica$phys, df.robusta$phys), na.rm = TRUE)),
       main = paste("Coffea, chromosome ", chr.nb, sep = ""), xlab = "Physical position (Mb)", ylab = "Recombination rate (cM/Mb)") #check
  lines(df.arabica$phys, df.arabica$upper, type = "l", lty = 2, col = "Grey")
  lines(df.arabica$phys, df.arabica$lower, type = "l", lty = 2, col = "Grey")
  abline(h = 0, col = "Black", lty = 2) # Print a threshold for rec. rate = 0
  
  lines(df.robusta$phys, df.robusta$rec.rate, type="l", col = "Red")
  lines(df.robusta$phys, df.robusta$upper, type = "l", lty = 2, col = "Red")
  lines(df.robusta$phys, df.robusta$lower, type = "l", lty = 2, col = "Red")
  dev.off()
}


# Add the syntenic markers information
# Common markers
set.name = "Coffea_arabica_Crouzillat2020"
marey.arabica = read.table(file = paste(wd, "marey_maps/", set.name, ".txt", sep = ""), sep = "\t", header = TRUE)
set.name = "Coffea_canephora_Crouzillat2020"
marey.robusta = read.table(file = paste(wd, "marey_maps/", set.name, ".txt", sep = ""), sep = "\t", header = TRUE)
# Common markers have the same name in Robusta and Arabica, yet they have a suffix (A) or (B) in arabica to remove
marey.arabica$mkr = as.character(marey.arabica$mkr)
marey.robusta$mkr = as.character(marey.robusta$mkr)

commonMarkers = marey.arabica$mkr[which(gsub("\\([A-B]+\\)$", "", marey.arabica$mkr) %in% marey.robusta$mkr)]

write.table(commonMarkers, file = paste(wd, "commonMarkers.txt", sep = ""), col.names = F, row.names = F, quote = F)




for (chr.nb in commonChromosomes) {
  set.name = "Coffea_arabica_Crouzillat2020"
  df.arabica = read.table(file = paste(wd, "recombination_maps/", set.name, "_chromosome", chr.nb, "_R.txt", sep = ""), sep = "\t", header = TRUE)
  marey.arabica = read.table(file = paste(wd, "marey_maps/", set.name, ".txt", sep = ""), sep = "\t", header = TRUE)
  set.name = "Coffea_canephora_Crouzillat2020"
  df.robusta = read.table(file = paste(wd, "recombination_maps/", set.name, "_chromosome", chr.nb, ".txt", sep = ""), sep = "\t", header = TRUE)
  marey.robusta = read.table(file = paste(wd, "marey_maps/", set.name, ".txt", sep = ""), sep = "\t", header = TRUE)
  
  # Which markers are common in Arabica and Robusta?
  # Some markers have the same name in Robusta and Arabica, yet they have a suffix (A) or (B) in arabica to remove
  marey.arabica$mkr = as.character(marey.arabica$mkr)
  marey.robusta$mkr = as.character(marey.robusta$mkr)
  # Subset to marker on the chromosome
  marey.arabica = subset(marey.arabica, map == paste(chr.nb, "_R", sep = ""))
  marey.robusta = subset(marey.robusta, map == chr.nb)
  
  commonMarkers = marey.arabica$mkr[which(gsub("\\([A-B]+\\)$", "", marey.arabica$mkr) %in% marey.robusta$mkr)]
  #commonMarkers = unique(marey.arabica$mkr[grep("\\([A-B]+\\)$", marey.arabica$mkr)])
  #marey.robusta$mkr
  
  # Extract x-y positions of common markers in both maps
  # First physical positions in Marey maps
  commonArabica = data.frame(mkr = commonMarkers, x = NA, y = NA)
  for (i in 1:nrow(commonArabica)) {
    commonArabica$x[i] = (marey.arabica$phys[marey.arabica$mkr == commonArabica$mkr[i]])/1000000
  }
  
  commonRobusta = data.frame(mkr = gsub("\\([A-B]+\\)$", "", commonMarkers), x = NA, y = NA)
  for (i in 1:nrow(commonRobusta)) {
    if (commonRobusta$mkr[i] %in% marey.robusta$mkr) {
      commonRobusta$x[i] = (marey.robusta$phys[marey.robusta$mkr == commonRobusta$mkr[i]])/1000000
    }
  }
  
  # How many common markers?
  sum(!is.na(commonRobusta$x))
  
  # Find the window in the recombination landscape
  # i.e. the physical window position in the recombination landscape immediately after the marker position
  for (i in 1:nrow(commonArabica)) {
    commonArabica$y[i] = df.arabica$rec.rate[min(which(df.arabica$phys > commonArabica$x[i]), na.rm = TRUE)]
  }
  for (i in 1:nrow(commonRobusta)) {
    commonRobusta$y[i] = df.robusta$rec.rate[min(which(df.robusta$phys > commonRobusta$x[i]), na.rm = TRUE)]
  }
  
  png(paste(wd, "combined_recombination_landscapes/Combined_Coffea_chr", chr.nb,".png",sep = ""), width = 1200, height = 600)
  # Arabica map
  plot(df.arabica$phys, df.arabica$rec.rate, type="l", ylim = c(0, 100), xlim = c(0 , max(c(df.arabica$phys, df.robusta$phys), na.rm = TRUE)),
       main = paste("Coffea, chromosome ", chr.nb, sep = ""), xlab = "Physical position (Mb)", ylab = "Recombination rate (cM/Mb)") #check
  lines(df.arabica$phys, df.arabica$upper, type = "l", lty = 2, col = "Grey")
  lines(df.arabica$phys, df.arabica$lower, type = "l", lty = 2, col = "Grey")
  abline(h = 0, col = "Black", lty = 2) # Print a threshold for rec. rate = 0
  # Robusta map
  lines(df.robusta$phys, df.robusta$rec.rate, type="l", col = "Red")
  lines(df.robusta$phys, df.robusta$upper, type = "l", lty = 2, col = "Red")
  lines(df.robusta$phys, df.robusta$lower, type = "l", lty = 2, col = "Red")
  
  # Draw arrows between common markers
  arrows(x0 = commonArabica$x, y0 = commonArabica$y, x1 = commonRobusta$x, y1 = commonRobusta$y, col = "Blue", code = 0)
  
  dev.off()
}


# Scaled version on the physical distances (physical distances from 0 to 1)
for (chr.nb in commonChromosomes) {
  set.name = "Coffea_arabica_Crouzillat2020"
  df.arabica = read.table(file = paste(wd, "recombination_maps/", set.name, "_chromosome", chr.nb, "_R.txt", sep = ""), sep = "\t", header = TRUE)
  marey.arabica = read.table(file = paste(wd, "marey_maps/", set.name, ".txt", sep = ""), sep = "\t", header = TRUE)
  set.name = "Coffea_canephora_Crouzillat2020"
  df.robusta = read.table(file = paste(wd, "recombination_maps/", set.name, "_chromosome", chr.nb, ".txt", sep = ""), sep = "\t", header = TRUE)
  marey.robusta = read.table(file = paste(wd, "marey_maps/", set.name, ".txt", sep = ""), sep = "\t", header = TRUE)
  
  # Scaling physical distances
  df.arabica$phys = df.arabica$phys/max(df.arabica$phys, na.rm = TRUE)
  df.robusta$phys = df.robusta$phys/max(df.robusta$phys, na.rm = TRUE)
  marey.arabica$phys = marey.arabica$phys/max(marey.arabica$phys, na.rm = TRUE)
  marey.robusta$phys = marey.robusta$phys/max(marey.robusta$phys, na.rm = TRUE)
  
  # Common markers have the same name in Robusta and Arabica, yet they have a suffix (A) or (B) in arabica to remove
  marey.arabica$mkr = as.character(marey.arabica$mkr)
  marey.robusta$mkr = as.character(marey.robusta$mkr)
  # Subset to marker on the chromosome
  marey.arabica = subset(marey.arabica, map == paste(chr.nb, "_R", sep = ""))
  marey.robusta = subset(marey.robusta, map == chr.nb)
  
  commonMarkers = marey.arabica$mkr[which(gsub("\\([A-B]+\\)$", "", marey.arabica$mkr) %in% marey.robusta$mkr)]
  #commonMarkers = unique(marey.arabica$mkr[grep("\\([A-B]+\\)$", marey.arabica$mkr)])
  #marey.robusta$mkr
  
  # Extract x-y positions of common markers in both maps
  # First physical positions in Marey maps
  commonArabica = data.frame(mkr = commonMarkers, x = NA, y = NA)
  for (i in 1:nrow(commonArabica)) {
    commonArabica$x[i] = (marey.arabica$phys[marey.arabica$mkr == commonArabica$mkr[i]])
  }
  
  commonRobusta = data.frame(mkr = gsub("\\([A-B]+\\)$", "", commonMarkers), x = NA, y = NA)
  for (i in 1:nrow(commonRobusta)) {
    if (commonRobusta$mkr[i] %in% marey.robusta$mkr) {
      commonRobusta$x[i] = (marey.robusta$phys[marey.robusta$mkr == commonRobusta$mkr[i]])
    }
  }
  
  # How many common markers?
  sum(!is.na(commonRobusta$x))
  
  # Find the window in the recombination landscape
  # i.e. the physical window position in the recombination landscape immediately after the marker position
  for (i in 1:nrow(commonArabica)) {
    commonArabica$y[i] = df.arabica$rec.rate[min(which(df.arabica$phys > commonArabica$x[i]), na.rm = TRUE)]
  }
  for (i in 1:nrow(commonRobusta)) {
    commonRobusta$y[i] = df.robusta$rec.rate[min(which(df.robusta$phys > commonRobusta$x[i]), na.rm = TRUE)]
  }
  
  
  png(paste(wd, "combined_recombination_landscapes_scaled/Combined_Coffea_chr", chr.nb,"_scaled.png",sep = ""), width = 1200, height = 600)
  plot(df.arabica$phys, df.arabica$rec.rate, type="l", ylim = c(0, 100), xlim = c(0 , max(c(df.arabica$phys, df.robusta$phys), na.rm = TRUE)),
       main = paste("Coffea, chromosome ", chr.nb, sep = ""), xlab = "Physical position (Mb)", ylab = "Recombination rate (cM/Mb)") #check
  lines(df.arabica$phys, df.arabica$upper, type = "l", lty = 2, col = "Grey")
  lines(df.arabica$phys, df.arabica$lower, type = "l", lty = 2, col = "Grey")
  abline(h = 0, col = "Black", lty = 2) # Print a threshold for rec. rate = 0
  
  lines(df.robusta$phys, df.robusta$rec.rate, type="l", col = "Red")
  lines(df.robusta$phys, df.robusta$upper, type = "l", lty = 2, col = "Red")
  lines(df.robusta$phys, df.robusta$lower, type = "l", lty = 2, col = "Red")
  # Draw arrows between common markers
  arrows(x0 = commonArabica$x, y0 = commonArabica$y, x1 = commonRobusta$x, y1 = commonRobusta$y, col = "Blue", code = 0)
  
  dev.off()
}





#----------------------------------------------------------------------
# Sliding window
#----------------------------------------------------------------------
# "sliding window, which is the simplest and most widely used method. The idea is just sliding a window along a chromosome and getting the local 
# estimate with the slope of the best line fit to the data in the local window. Parameters are window size and shift." Rezvoy et al. 2007

# May be the worst interpolation method, so no implementation yet

#----------------------------------------------------------------------
# Cubic spline
#----------------------------------------------------------------------
# "cubic splines, which is probably the best method to estimate recombination rate with Marey map approach (Berloff et al., 2002; Yu et al., 2001)." Rezvoy et al. 2007
# ?smooth.spline # As in MareyMap
# Another implementation is the function qsreg for robust smoothing spline from the package fields (Oh, et al., 2004),
# For a chosen dataset
#recombination.map(set = set, chr = "all", method = "smooth.spline", K = 5, boot = 1000)

# Cubic smooth splines seems to be very sensitive to missing data and produces some artefacts
# Lower confidence in results
# Estimates are smoother
# But may be used as a confirmation of the global pattern for LOESS

# Below are methods of interest to improve esimates, but they have to be tested

#============================================================================
# Validate the recombination map
#============================================================================
# A new round of validation for recombination maps
# (1) Maps must not present regions of inconsistent recombination rate (i.e. over 30 cM.Mb-1)
# (2) CI must be under a threshold of quality
# (3) loess and smooth.spline must show convergence
# (4) Maps must not show a high discrepancy between the chromosome-wide recombiantion rate and the estimated mean recombination rate

# Otherwise, maps are discarded in the original data and outputs/figures are removed
discard.map() # Remove data in the original dataset

# Then analyses are runned again to update the results


#----------------------------------------------------------------------
# RMSE method
#----------------------------------------------------------------------





#============================================================================
# Significance testing of the recombination map
#============================================================================
# Which part of the recombination landscape have a recombination rate significantly higher than the average recombination rate?





#============================================================================
# Save the recombination map
#============================================================================

#----------------------------------------------------------------------
# Saving in the directory '/output/MareyMaps'
# Maps can be saved to R data files (rda, Rda, rdata or Rdata) or to text files (txt).


#----------------------------------------------------------------------
# Save in a txt file the estimated local recombination rate for each marker position



#----------------------------------------------------------------------
# Indicate the centromere position
# BLAST the centromeric motif, since centromeres are not necessarily at the center of the chromosome (e.g. metacentric, acrocentric or telocentric)

#============================================================================
# END
#============================================================================