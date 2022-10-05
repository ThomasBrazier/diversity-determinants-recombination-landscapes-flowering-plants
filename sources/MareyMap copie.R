# Work with Marey maps

#============================================================================
# Check orientation of the map
#============================================================================
# Check & correct bad orientation of linkage groups
# Markers oriented in the direction of a reference map if necessary
# A reference genome with markers mapped onto it seems inevitable for an exhaustive check
# set = name of the dataset, ONE map
# chr = chromosome number to flip entirely, automatic check if chr = c()
# dir = path to the dataset

# Yet, the Marey map curve can only be ascending
check_map_orientation = function(set = "", id_chr = c(), dir = "data/Marey_maps/", auto = FALSE) {
    marey = read.table(paste(dir , set, sep = ""), header = TRUE)
    # The function gen~phys must be monotonously increasing
    # In case you search to identify segments going the wrong direction in a chromosome
    if (auto == TRUE) {
      new_map = data.frame()
      for (chr in unique(marey$map)) {
        # Find groups of markers with descending genetic position
        # (i.e. genetic and physical maps not oriented in the same direction)
        flip_group = c()
        #---------------------------------------------------
        # Exact approach
        #---------------------------------------------------
        # "Three markers interspaced at least one cM is considered the minimum requisite
        # to ensure the same orientation of the map.
        # For each interval between successive markers, compute if curve is increasing/decreasing
        # A segment is a vector of more than intervals of the same direction
        # and subsequently oriented in the direction of a reference map if supplied" (Mapfuser, van Muijen et al. 2017)
        # subset = marey[marey$map == chr,]
        # # remove NA values
        # na_phys = subset[is.na(subset$phys),]
        # subset = subset[!is.na(subset$phys),]
        # # sort subset by physical position
        # subset = subset[order(subset$phys),]
        # for (n in 1:(nrow(subset))) {
        #   print(n)
        #   # If (n+1) - n < 0 for two consecutive markers, then genetic positions are in the wrong direction
        #   # Gather indexes of markers, until a new group of markers going in the true direction is found,
        #   # Once two good markers in a row are encountered, if the group has more than one marker, flip indexed markers
        #   if ((subset$gen[n+1] - subset$gen[n]) < 0) {
        #     flip_group = append(flip_group, n)
        #   } else {
        #     if ((subset$gen[n+2] - subset$gen[n+1]) > 0 & (length(flip_group) > 1)) {
        #       # Flip genetic distances within the interval of markers (ascending order)
        #       subset$gen[min(flip_group):max(flip_group)] = abs(subset$gen[min(flip_group):max(flip_group)] - max(subset$gen[min(flip_group):max(flip_group)], na.rm = TRUE))
        #       # Re-init indexes
        #       flip_group = c()
        #     }
        #   }
        #   # At the end of the chromosome, flip the last group remaining
        #   if (length(flip_group) > 1) {
        #     # Flip genetic distances within these groups of markers (ascending order)
        #     subset$gen[flip_group] = abs(subset$gen[flip_group] - max(subset$gen[flip_group], na.rm = TRUE))
        #   }
        #   # New positions and NA values
        #   new_map = rbind(new_map, subset, na_phys)
        #   rm(subset)
        # }
        #---------------------------------------------------
        # Likelihood approach
        #---------------------------------------------------
        require(segmentr)
        require(tidyr)
        require(tibble)
        require(dplyr)
        require(lubridate)
        require(magrittr)
        require(purrr)
        # Identify monotonous segments with 'segmentr'
        # Segmentr is a package that implements a handful of algorithms to segment a given data set,
        # by finding the change points that maximize the collective likelihood of the segments according to an arbitrary likelihood function
        lm_likelihood = function (data) {
          fit = t(data) %>%
            as_tibble() %>%
            rowid_to_column() %>%
            gather(station, temperature, -rowid) %>%
            with(lm(temperature ~ rowid))
          
          -mean(fit$residuals ^ 2)
        }
        # Penalize the likelihood to avoid small segments
        penalized_likelihood = auto_penalize(marey, lm_likelihood)
        seg = segment(data = marey,
          likelihood = penalized_likelihood,
          algorithm = "hierarchical"
        )
        # Tuning the likelihood to obtain segments fitting to the data
        boundaries = c()
        # Identify which segments are decreasing
        flip = c()
        # Flip decreasing segments
        
        save.log(msg = paste("check_map_orientation: Found ", length(flip)," segments that have been flipped", sep=""))

    }
      } else { # in case of chromosome ids to flip entirely (i.e. auto = FALSE)
        new_map = marey[-which(marey$map %in% unique(marey$map)[id_chr]),] # Unchanged chromosome are saved in new_map before procedure
      for (chr in unique(marey$map)[id_chr]) {
        subset = marey[marey$map == chr,]
        # remove NA values
        na_phys = subset[is.na(subset$phys),]
        subset = subset[!is.na(subset$phys),]
        # sort subset by physical position
        subset = subset[order(subset$phys),]
        # Flip the entire chromosome genetic distances
        subset$gen = abs(subset$gen - max(subset$gen))
        # New positions and NA values
        new_map = rbind(new_map, subset, na_phys)
        rm(subset)
      }

      }
    #
    
    # Saving the call to the function in a log file for reproducibility purposes
    save.log()
    # Return the map woth corrections applied as a data frame
    return(new_map)
}



#============================================================================
# Filter out outliers automatically
#============================================================================
# Selection of valid markers

# How to filter outliers?
# 1/ Fit a curve and filter markers away from the curve (threshold):
# bootstrap the dist. of the mean departure from the curve, exclude points far away than 95% (2.5% each side)
# Advantage, selection is curve specific
# 2/ Sample the subset of markers minimizing the RMSE (resampling approach)


#============================================================================
# Filter out chromosomes with insufficient number of markers,
#============================================================================
# based on an interval threshold or number of markers specified by the user
# mean interval distance



#============================================================================
# Diagnostic plots
#============================================================================
# Genetic distance (cM) as a function of physical position (bp)

# set (character) is a vector of names of Marey map files to plot, or a data.frame object
# dir (character) is the path of datasets
# output (character) is the path where to save plots
# save = TRUE if you want to save your plot as a png file
# display = TRUE if you want to print the plot in Rstudio

marey.diagplot = function(set = "", dir = "data/Marey_maps/", output = "output/Marey_maps/diagnostic_plots/", save.plot = TRUE, display = TRUE) {
  require(ggplot2)
  if (is.data.frame(set)) {
    marey = set
    marey = marey[!(is.na(marey$phys) | is.na(marey$gen) | is.na(marey$map)),]
    diagnostic_plot = ggplot(data=marey, aes(x=phys/1000000, y=gen, col = vld)) +
      geom_point() +
      scale_colour_manual(values=c("TRUE" = "dodgerblue4", "FALSE" = "firebrick4")) +
      facet_wrap(~as.factor(map), nrow=5) +
      labs(x="Physical position (Mb)", y="Genetic distance (cM)") +
      theme(axis.line = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(color="black", size=22, face="bold.italic",hjust = 0.5),
            plot.subtitle = element_text(color="black",size=22,hjust = 0.5),
            axis.title.x = element_text(color="black", size=22),
            axis.title.y = element_text(color="black", size=22),
            axis.text=element_text(size=22, colour="black"),
            # axis.text.x=element_blank(),
            strip.text=element_text(size=18, colour="black"),
            legend.position = "none")
    plot(diagnostic_plot)
    rm(diagnostic_plot)
    rm(marey)
  } else {
    if (is.vector(set)) {
      for (i in set) {
        marey = read.table(paste(dir , i, sep = ""), header = TRUE, sep = "\t")
        print(i)
        marey = marey[!(is.na(marey$phys) | is.na(marey$gen) | is.na(marey$map)),]
        marey$phys = as.numeric(marey$phys)
        marey$gen = as.numeric(marey$gen)
        
        diagnostic_plot = ggplot(data=marey, aes(x=phys/1000000, y=gen, col = vld)) +
          geom_point() +
          scale_colour_manual(values=c("TRUE" = "dodgerblue4", "FALSE" = "firebrick4")) +
          facet_wrap(~as.factor(map), nrow=5) +
          labs(x="Physical position (Mb)", y="Genetic distance (cM)") +
          theme(axis.line = element_line(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.border = element_blank(),
                panel.background = element_blank(),
                plot.title = element_text(color="black", size=22, face="bold.italic",hjust = 0.5),
                plot.subtitle = element_text(color="black",size=22,hjust = 0.5),
                axis.title.x = element_text(color="black", size=22),
                axis.title.y = element_text(color="black", size=22),
                axis.text=element_text(size=22, colour="black"),
                # axis.text.x=element_blank(),
                strip.text=element_text(size=18, colour="black"),
                legend.position = "none")
        
        if (save.plot == TRUE) {
          ggsave(paste(output , gsub(".txt", "", i), ".png", sep = ""),
                 device="png",dpi=320,units="cm",width=60,height=60)
        }
        
        if (display == TRUE) {
          plot(diagnostic_plot)  
        }
        
        rm(diagnostic_plot)
        rm(marey)
      }
    } else {
      warning("'set' must be a vector or a data frame")
    }
  }
  
}



#============================================================================
# Map a physical map onto a genetic map
#============================================================================
# Consolidate a MareyMap input file
# i.e. Physical and genetic positons of corresponding markers merged into one Marey map input file

# Access to the markers databse for the species and retrieve physical positions of the genetic map markers
# Otherwise, fail to recover markers position and give a warning
consolidateMareyMap = function(set = "",  markerDB = "data/Physical_maps/Markers/", input = "data/Genetic_maps/", output = "data/Marey_maps/") {
  require(gdata)
  # Saving the call to the function in a log file for reproducibility purposes
  save.log()
  for (i in set) {
    cat("Processing", i, "\n")
    # Take a list of genetic map as input
    gen_map = read.table(paste(input, i, sep = ""), header = TRUE)
    # Get species name
    sp = regmatches(as.character(i), gregexpr("^[A-Z][a-z]+_[a-z]+", as.character(i)))
    # Import the marker database
    phys_map = read.csv(paste(markerDB, sp, ".csv", sep = ""), header = TRUE, sep = ";")
    # Reduce MarkerDB to the species
    # phys_map = phys_map[phys_map$species == sp,] #### Deprecated
    
    #-------------------------------------
    cat("Matching physical positions...\n")
    # Matching physical positions with marker name
    match = as.character(gen_map$Marker) # List of markers in the genetic map
    map = rep(NA, nrow(gen_map)) # The vector of chromosome numbers associated to each marker in genetic map
    no_match = c() # Index of markers without matching physical position
    
    # For each marker:
    for (j in 1:length(match)) {
      if (match[j] %in% as.character(phys_map$mkr)) {
        map[j] = phys_map$map[which(match[j] == as.character(phys_map$mkr))] # The chromosome associated to the position
        match[j] = as.character(phys_map$phys[which(match[j] == as.character(phys_map$mkr))]) # Marker ID is replaced by corresponding physical position
        # If the marker physical position is NA, then vld = FALSE
        if (is.na(match[j])) {no_match = c(no_match, j)}
      } else {
        no_match = c(no_match, j) # The marker is not in the physical map, no position
        match[j] = NA
      }
    }
    
    # Number of markers without position
    cat("Number of markers with no physical position retrieved:", length(no_match), "\n")
    save.log(msg = paste("Number of markers with no physical position retrieved: ", length(no_match), sep = ""))
    # Consolidate the Marey map
    marey = data.frame(set = sp, map = map, mkr = gen_map$Marker, phys = match,
                       gen = gen_map$cM, vld = TRUE)
    colnames(marey) = c("set", "map", "mkr", "phys", "gen", "vld")
    # Unvalidate markers with NA in physical positions
    marey$vld[is.na(marey$phys)] = FALSE
    
    
    
    #-------------------------------------
    # Map chromosome number for each marker
    # When the information is not available, use linkage group (LG) instead of chromosome number in genetic maps
    # Map the chromosome number afterward onto the Marey_map input file
    
    # Which genetic maps have LG numbers (i.e. LG[0-9]* in the column 'Chr')?
    
    
    #-------------------------------------
    # Save the Marey map input file
    write.table(marey, file = paste(output, set, sep = ""), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    # Memory clean
    rm(sp)
    rm(gen_map)
    rm(phys_map)
    rm(match)
    rm(no_match)
    rm(marey)
  }
}



#============================================================================
# Compute the number of markers per genetic map
#============================================================================
# Return a data frame with the map name and the number of markers associated
get.nbMkr = function(path = "data/Genetic_maps/") {
  metadata = read.csv(file = paste(path, "Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")
  set = metadata$id[which(!is.na(metadata$id))]
  df = data.frame(set = set, nbMkr = rep(NA, length(set)))
  
  for (i in 1:length(set)) {
    if (file.exists(paste(path, set[i], ".txt", sep = ""))) {
      map = read.csv(file = paste(path, set[i], ".txt", sep = ""), header = TRUE, sep = ";")
      df$set[i] = set[i]
      df$nbMkr[i] = sum(map$vld)
    } else {
      df$set[i] = set[i]
      df$nbMkr[i] = NA
    }
  }
  return(df)
}

#============================================================================
# Compute the linkage map length
#============================================================================
# Return a data frame with the map name and the number of markers associated
# The total linkage map length is the sum of max genetic distances of all chromosomes
get.maplength = function(path = "data/Genetic_maps/") {
  metadata = read.csv(file = paste(path, "Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")
  set = metadata$id[which(!is.na(metadata$id))]
  df = data.frame(set = set, linkage_map_length = rep(NA, length(set)))
  
  for (i in 1:length(set)) {
    if (file.exists(paste(path, set[i], ".txt", sep = ""))) {
      map = read.csv(file = paste(path, set[i], ".txt", sep = ""), header = TRUE, sep = "\t")
      df$set[i] = set[i]
      # The total linkage map length is the sum of max genetic distances of all chromosomes
      map = map[!is.na(map$Chr),]
      map = map[!is.na(map$cM),]
      map$Chr = as.factor(map$Chr)
      map$cM = as.numeric(map$cM)
      subset = split(map$cM, map$Chr)
      max.distances = lapply(subset, function(x) max(x, na.rm = TRUE))
      df$linkage_map_length[i] = sum(unlist(max.distances), na.rm = TRUE)
    } else {
      df$set[i] = set[i]
      df$linkage_map_length[i] = NA
    }
  }
  # Remove NA map length
  df = df[!is.na(df$linkage_map_length),]
  return(df)
}


#============================================================================
# Estimation of the local recombination rate (cM/Mb)
#============================================================================
recombination.map = function(set ="", data = data.final, scale = 1000000, chr = "all", method = c("loess", "smooth.spline"), K = 5, boot = 1000,...) { # Optional arguments passed to loess, smooth.spline and pther fitting methods
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
    chr.data = subset(data.final, map == chr.nb & set == set.name) # Sample only distances of the chromosome
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
    # BOOTSTRAPS  
    if (boot == 0) {
      
    } else {
      if (boot > 0 & is.numeric(boot)) {
        if (method == "loess") {
          
        } else {
          if (method == "smooth.spline") {
            
          } else {
            warning("Method must be loess or smooth.spline.")
          }
        }
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
    
    
    #-----------------------------------------
    # LOESS
    #-----------------------------------------
    if ((adj.span == "auto") & (method == "loess")) {
      span2test = seq(from = 0.1, to = 0.5, len = 9) # Values of parameter to compare
      # A cross-validation data frame with
      # value of span to test
      # physical positions
      # observed genetic position
      # Predicted genetic positions
      size = length(span2test)*nrow(chr.data) # Size of the data frame (number of rows)
      crossvalidation = data.frame(span = rep(span2test, each = nrow(chr.data)), phys = rep(chr.data$phys, length(span2test)), obs = rep(chr.data$gen, length(span2test)), pred = numeric(size))
      # First, subsetting the dataset in K clusters
      K = 5
      # Systematic K-fold selection, points are uniformly sampled along the map (a sequence of K-th nearest neighbors)
      # i.e., for K = 5, sample 1, 6, 11, 17... then 2, 7, 12, 18...
      subset = rep(seq(1,5), ceiling(nrow(chr.data)/5)) # Vector of clusters' number
      subset = subset[1:nrow(chr.data)] # Trim the vector to its maximal size
      
      if (method == "loess") {
        for (s in span2test) {
          # Run K times the cross-validation
          # Each time the cluster K is leaved out for fitting and used as validation
          # results.cv = data.frame(K = 1:K, ACV = numeric(K), LSCV = numeric(K))
          for (k in 1:K) {
            training = rbind( chr.data[1,], subset(chr.data, subset %in% seq(1,5)[-k]), chr.data[nrow(chr.data),]) # Only data out of cluster K + boundaries of the map
            validation = subset(chr.data, subset == k)
            
            fit.loess.cv = loess(gen~phys, training, span = s, degree = degree) # Same defaults parameters as in MareyMap
            summary(fit.loess.cv)
            # Predict the local recombination rate at the exact physical positions of the Marey map
            crossvalidation$pred[which(crossvalidation$span == s)[subset == k]] = predict(fit.loess.cv, newdata = validation$phys)
            # validation.loess = data.frame(pos = validation$phys, obs = validation$gen, pred = predict(fit.loess.cv, newdata = validation$phys)) # Predict the local recombination rate at the exact physical positions of the Marey map
            # validation.loess
            
            
          }
        }
        
        # Compute the criterion for each span
        # source("sources/stats.R")
        crit = data.frame(span = span2test, ACV = numeric(length(span2test)), LSCV = numeric(length(span2test)))
        for (s in span2test) {
          # Absolute Cross-validation
          # Compute the absolute value of the difference between the mean observed value and the predicted value in the left out dataset
          # lambda value to optimize (the span parameter)
          # Input is a data frame with three columns of physical position and observed/predicted values
          crit$ACV[crit$span == s] = ACV(data = crossvalidation[which(crossvalidation$span == s), -1])
          
          # Least squares Cross validation
          crit$LSCV[crit$span == s] = LSCV(data = crossvalidation[which(crossvalidation$span == s), -1])
        }
        
        # Automated selection of the span value which minimizes the chosen criterion
        span = crit$span[crit$ACV == min(crit$ACV, na.rm = TRUE)]
      } else {
        if (is.numeric(span)) {
          span = span
        } else {
          warning("The argument 'span' must be 'auto' or a numeric value between 0 and 1.\n")
        }
      }
      
      #-----------------------------------------
      # FITTING
      fit = loess(gen~phys, chr.data, span = span, degree = degree) # Same defaults parameters as in MareyMap
      # summary(fit)
      # plot(chr.data$gen~chr.data$phys)
      # lines(fit.loess, col = "Blue", lwd = 2)
      # predicted = data.frame(obs = chr.data$phys, pred = predict(fit, newdata = chr.data$phys, se = TRUE)$fit, se.pred = predict(fit, newdata = chr.data$phys, se = TRUE)$se.fit) # Predict the local recombination rate at the exact physical positions of the Marey map, with a s.e.
      # points(pred.loess, col = "Red")
      # Add a confidence interval
      # points(pred.loess$obs, pred.loess$pred+1.96*pred.loess$se.pred, col="Red",type="l", lwd=1, lty=2)
      # points(pred.loess$obs, pred.loess$pred-1.96*pred.loess$se.pred, col="Red",type="l", lwd=1, lty=2)
      # lines(fit.loess, col = "Blue", lwd = 2)
      
      # Save actions in log with parameters
      # print(set)
      save.log(path = "/output/", msg = paste("LOESS fitting on chromosome", chr.nb,"with span =", span, "(automated span optimization) and degree =", degree, ", with LOESS graphical predictions saved as figure", sep = " "), append = TRUE)
      
    }
    

    #-----------------------------------------
    # CUBIC SMOOTH SPLINE
    #-----------------------------------------
    if (method == "smooth.spline") {
      pointwise.boot = matrix(NA, nrow = length(chr.data$phys)-1, ncol = boot) # Predict the local recombination rate at the exact physical positions of the Marey map
      # Bin recombination map
      # Partitioning of the physical map (bp): 1,000 bins of equal size
      phys.bins = seq(0, max(chr.data$phys, na.rm = TRUE), len = 1000) # A sequence of bins coordinates
      physbins.boot = matrix(NA, nrow = length(phys.bins)-1, ncol = boot) # Predict the local recombination rate at the exact physical positions of the Marey map
      
      # 1Mb windows recombination map
      # Partitioning of the physical map in windows of 1Mb for further comparison with genomic content
      phys.1Mbwind = c(seq(0, max(chr.data$phys, na.rm = TRUE), by = 1), max(chr.data$phys, na.rm = TRUE)) # A sequence of windows coordinates
      phys.1Mbwind = rowMeans(embed(phys.1Mbwind,2)) # Take the center of each window for interpolation
      phys.1Mbwind.boot = matrix(NA, nrow = length(phys.1Mbwind)-1, ncol = boot)
      
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
        fit = fit.smooth.spline(chr.data = resampled, spar = spar, adj.spar = adj.spar)
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
        predicted.pointwise = predict(fit, x = chr.data$phys)$y # Predict the local recombination rate at the exact physical positions of the Marey map
        
        # Bin recombination map
        # Partitioning of the physical map (bp): 1,000 bins of equal size
        predicted.physbins = predict(fit, x = phys.bins)$y # Predict the local recombination rate at the exact physical positions of the Marey map
        
        # 1Mb windows recombination map
        # Partitioning of the physical map in windows of 1Mb for further comparison with genomic content
        predicted.phys.1Mbwind = predict(fit, x = phys.1Mbwind)$y # Predict the local recombination rate at the exact physical positions of the Marey map
        
        pointwise.boot[,b] = diff(predicted.pointwise)/diff(chr.data$phys)
        physbins.boot[,b] = diff(predicted.physbins)/diff(phys.bins)
        phys.1Mbwind.boot[,b] = diff(predicted.phys.1Mbwind)/diff(phys.1Mbwind)
      }
      Sys.sleep(1)
      close(pb)


      #------------------------------
      # Pointwise recombination map
      dY = rowMeans(pointwise.boot, na.rm = TRUE) # the derivative of LOESS function
      # Compute the 95% CI of the derivative, apply on rows
      upperCI = function(x) {
        quantile(x, 0.975)
      }
      dY.upper = apply(pointwise.boot, MARGIN = 1, FUN = upperCI)
      lowerCI = function(x) {
        quantile(x, 0.025)
      }
      dY.lower = apply(pointwise.boot, MARGIN = 1, FUN = lowerCI)
      dX = rowMeans(embed(chr.data$phys,2)) # centers the X values for plotting
      
      # Return a data.frame of four columns for each method -> final plotting df
      df.pointwise = data.frame(phys = dX, rec.rate = dY, upper = dY.upper, lower = dY.lower)
      write.table(df.pointwise, file = paste("output/recombination_maps/", method,"/bins/", set.name, "_chromosome", chr.nb, ".txt", sep = ""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
      
      # Save graphical predictions
      png(paste("figures/mareymap/", method, "/pointwise/", set.name,"_chr", chr.nb,".png",sep = ""), width = 1200, height = 600)
      plot(dX, dY, type="l", main = paste("Recombination landscape (", method, ") of chromosome ", chr.nb, " in ", set.name, sep = ""), xlab = "Physical position (Mb)", ylab = "Recombination rate (cM/Mb)") #check
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
      upperCI = function(x) {
        quantile(x, 0.975)
      }
      dY.upper = apply(physbins.boot, MARGIN = 1, FUN = upperCI)
      lowerCI = function(x) {
        quantile(x, 0.025)
      }
      dY.lower = apply(physbins.boot, MARGIN = 1, FUN = lowerCI)
      dX = rowMeans(embed(phys.bins,2)) # centers the X values for plotting
      
      # Return a data.frame of four columns for each method -> final plotting df
      df.bins = data.frame(phys = dX, rec.rate = dY, upper = dY.upper, lower = dY.lower)
      write.table(df.bins, file = paste("output/recombination_maps/", method,"/bins/", set.name, "_chromosome", chr.nb, ".txt", sep = ""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
      
      # Save graphical predictions
      png(paste("figures/mareymap/", method, "/bins/", set.name,"_chr", chr.nb,".png",sep = ""), width = 1200, height = 600)
      plot(dX, dY, type="l", main = paste("Recombination landscape (", method, ") of chromosome ", chr.nb, " in ", set.name, sep = ""), xlab = "Physical position (Mb)", ylab = "Recombination rate (cM/Mb)") #check
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
      upperCI = function(x) {
        quantile(x, 0.975)
      }
      dY.upper = apply(phys.1Mbwind.boot, MARGIN = 1, FUN = upperCI)
      lowerCI = function(x) {
        quantile(x, 0.025)
      }
      dY.lower = apply(phys.1Mbwind.boot, MARGIN = 1, FUN = lowerCI)
      dX = rowMeans(embed(phys.1Mbwind,2)) # centers the X values for plotting
      
      # Return a data.frame of four columns for each method -> final plotting df
      df.1Mbwind = data.frame(phys = dX, rec.rate = dY, upper = dY.upper, lower = dY.lower)
      write.table(df.1Mbwind, file = paste("output/recombination_maps/", method,"/bins/", set.name, "_chromosome", chr.nb, ".txt", sep = ""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
      
      # Save graphical predictions
      png(paste("figures/mareymap/", method, "/1Mbwind/", set.name,"_chr", chr.nb,".png",sep = ""), width = 1200, height = 600)
      plot(dX, dY, type="l", main = paste("Recombination landscape (", method, ") of chromosome ", chr.nb, " in ", set.name, sep = ""), xlab = "Physical position (Mb)", ylab = "Recombination rate (cM/Mb)") #check
      lines(dX, dY.upper, type = "l", lty = 2, col = "Grey")
      lines(dX, dY.lower, type = "l", lty = 2, col = "Grey")
      abline(h = 0, col = "Red", lty = 2) # Print a threshold for rec. rate = 0
      dev.off()
      rm(dX, dY, dY.upper, dY.lower)
      rm(phys.1Mbwind.boot)

            
    }
    
    #-----------------------------------------
    # DERIVATIVES
    # The local recombination rate is the derivative of the polynom (i.e. fitted loess) in a point
    # dY = diff(pred.loess$pred)/diff(pred.loess$obs) # the derivative of LOESS function
    # # Compute the 95% CI of the derivative
    # dY.upper = diff(pred.loess$pred+1.96*pred.loess$se.pred)/diff(pred.loess$obs)
    # dY.lower = diff(pred.loess$pred-1.96*pred.loess$se.pred)/diff(pred.loess$obs)
    # dX = rowMeans(embed(pred.loess$obs,2)) # centers the X values for plotting
    # plot(dX, dY, type="l", main = paste("Recombination landscape of chromosome", chr.nb, "in", set.name, sep = " "), xlab = "Physical position (bp)", ylab = "Recombination rate (cM)") #check
    # lines(dX, dY.upper, type = "l", lty = 2)
    # lines(dX, dY.lower, type = "l", lty = 2)
    # abline(h = 0, col = "Red", lty = 2) # Print a threshold for rec. rate = 0
    
    # dY = diff(fit.loess.bis$loess)/diff(chr.data$phys) # the derivative of LOESS function
    # # Compute the 95% CI of the derivative
    # dY.upper = diff(fit.loess.bis$uci)/diff(chr.data$phys)
    # dY.lower = diff(fit.loess.bis$lci)/diff(chr.data$phys)
    # dX = rowMeans(embed(chr.data$phys, 2)) # centers the X values for plotting
    # plot(dX, dY, type="l", main = paste("Recombination landscape of chromosome", chr.nb, "in", set.name, sep = " "), xlab = "Physical position (bp)", ylab = "Recombination rate (cM)") #check
    # lines(dX, dY.upper, type = "l", lty = 2)
    # lines(dX, dY.lower, type = "l", lty = 2)
    # abline(h = 0, col = "Red", lty = 2) # Print a threshold for rec. rate = 0
    
    
    # Save the recombination map
    # Pointwise estimates (same as physical positions observed)
    # if (method == "loess") {
    #   predicted.pointwise = data.frame(phys = chr.data$phys, pred = predict(fit, newdata = chr.data$phys, se = TRUE)$fit, se.pred = predict(fit, newdata = chr.data$phys, se = TRUE)$se.fit) # Predict the local recombination rate at the exact physical positions of the Marey map, with a s.e.
    # }
    # if (method == "smooth.spline") {
    #   predicted.pointwise = data.frame(phys = chr.data$phys, pred = predict(fit, x = chr.data$phys)$y, se.pred = rep(NA, length(chr.data$phys))) # Predict the local recombination rate at the exact physical positions of the Marey map
    # }
    # dY = diff(predicted.pointwise$pred)/diff(predicted.pointwise$phys) # the derivative of LOESS function
    # # Compute the 95% CI of the derivative
    # dY.upper = diff(predicted.pointwise$pred+1.96*predicted.pointwise$se.pred)/diff(predicted.pointwise$phys)
    # dY.lower = diff(predicted.pointwise$pred-1.96*predicted.pointwise$se.pred)/diff(predicted.pointwise$phys)
    # dX = rowMeans(embed(predicted.pointwise$phys,2)) # centers the X values for plotting
    # rec.pointwise = data.frame(phys = dX, rec.rate = dY, upper.ci = dY.upper,lower.ci = dY.lower)
    # 
    # # Save graphical predictions
    # png(paste("figures/mareymap/", method, "/pointwise/", set.name,"_chr", chr.nb,".png",sep = ""), width = 1200, height = 600)
    # plot(dX, dY, type="l", main = paste("Recombination landscape (", method, ") of chromosome ", chr.nb, " in ", set.name, sep = ""), xlab = "Physical position (Mb)", ylab = "Recombination rate (cM/Mb)") #check
    # lines(dX, dY.upper, type = "l", lty = 2)
    # lines(dX, dY.lower, type = "l", lty = 2)
    # abline(h = 0, col = "Red", lty = 2) # Print a threshold for rec. rate = 0
    # dev.off()
    # write.table(rec.pointwise, file = paste("output/recombination_maps/loess/pointwise/", set.name, "_chromosome", chr.nb, ".txt", sep = ""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
    # rm(predicted.pointwise)
    # rm(rec.pointwise)
    # 
    # # Bin recombination map
    # # Partitioning of the physical map (bp): 1,000 bins of equal size
    # phys.bins = seq(0, max(chr.data$phys, na.rm = TRUE), len = 1000) # A sequence of bins coordinates
    # if (method == "loess") {
    #   predicted.physbins = data.frame(phys = phys.bins, pred = predict(fit, newdata = phys.bins, se = TRUE)$fit, se.pred = predict(fit, newdata = phys.bins, se = TRUE)$se.fit) # Predict the local recombination rate at the exact physical positions of the Marey map, with a s.e.
    # }
    # if (method == "smooth.spline") {
    #   predicted.physbins = data.frame(phys = phys.bins, pred = predict(fit, x = phys.bins)$y, se.pred = rep(NA, length(phys.bins))) # Predict the local recombination rate at the exact physical positions of the Marey map
    # }
    # dY = diff(predicted.physbins$pred)/diff(predicted.physbins$phys) # the derivative of LOESS function
    # # Compute the 95% CI of the derivative
    # dY.upper = diff(predicted.physbins$pred+1.96*predicted.physbins$se.pred)/diff(predicted.physbins$phys)
    # dY.lower = diff(predicted.physbins$pred-1.96*predicted.physbins$se.pred)/diff(predicted.physbins$phys)
    # dX = rowMeans(embed(predicted.physbins$phys,2)) # centers the X values for plotting
    # rec.physbins = data.frame(phys = dX, rec.rate = dY, upper.ci = dY.upper,lower.ci = dY.lower)
    # write.table(rec.physbins, file = paste("output/recombination_maps/", method,"/bins/", set.name, "_chromosome", chr.nb, ".txt", sep = ""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
    # # plot(rec.physbins$rec.rate)
    # png(paste("figures/mareymap/",method, "/bins/", set.name,"_chr", chr.nb,"_bins.png",sep = ""), width = 800, height = 400)
    # plot(dX, dY, type="l", main = paste("Recombination landscape (", method, ") of chromosome ", chr.nb, " in ", set.name, sep = ""), xlab = "Physical position (Mb) in 1,000 bins of equal distances", ylab = "Recombination rate (cM/Mb)") #check
    # lines(dX, dY.upper, type = "l", lty = 2)
    # lines(dX, dY.lower, type = "l", lty = 2)
    # abline(h = 0, col = "Red", lty = 2) # Print a threshold for rec. rate = 0
    # dev.off()
    # rm(predicted.physbins)
    # rm(rec.physbins)
    # 
    # # 1Mb windows recombination map
    # # Partitioning of the physical map in windows of 1Mb for further comparison with genomic content
    # phys.1Mbwind = c(seq(0, max(chr.data$phys, na.rm = TRUE), by = 1), max(chr.data$phys, na.rm = TRUE)) # A sequence of windows coordinates
    # phys.1Mbwind = rowMeans(embed(phys.1Mbwind,2)) # Take the center of each window for interpolation
    # 
    # if (method == "loess") {
    #   predicted.phys.1Mbwind = data.frame(phys = phys.1Mbwind, pred = predict(fit, newdata = phys.1Mbwind, se = TRUE)$fit, se.pred = predict(fit, newdata = phys.1Mbwind, se = TRUE)$se.fit) # Predict the local recombination rate at the exact physical positions of the Marey map, with a s.e.
    # }
    # if (method == "smooth.spline") {
    #   predicted.phys.1Mbwind = data.frame(phys = phys.1Mbwind, pred = predict(fit, x = phys.1Mbwind)$y, se.pred = rep(NA, length(phys.1Mbwind))) # Predict the local recombination rate at the exact physical positions of the Marey map
    # }
    # dY = diff(predicted.phys.1Mbwind$pred)/diff(predicted.phys.1Mbwind$phys) # the derivative of LOESS function + conversion in Mb (x 1000)
    # # Compute the 95% CI of the derivative
    # dY.upper = diff(predicted.phys.1Mbwind$pred+1.96*predicted.phys.1Mbwind$se.pred)/diff(predicted.phys.1Mbwind$phys)
    # dY.lower = diff(predicted.phys.1Mbwind$pred-1.96*predicted.phys.1Mbwind$se.pred)/diff(predicted.phys.1Mbwind$phys)
    # dX = rowMeans(embed(predicted.phys.1Mbwind$phys, 2)) # centers the X values for plotting
    # rec.phys.1Mbwind = data.frame(phys = dX, rec.rate = dY, upper.ci = dY.upper,lower.ci = dY.lower)
    # write.table(rec.phys.1Mbwind, file = paste("output/recombination_maps/", method,"/1Mbwind/", set.name, "_chromosome", chr.nb, ".txt", sep = ""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
    # # plot(rec.phys.1Mbwind$rec.rate)
    # png(paste("figures/mareymap/",method, "/1Mbwind/", set.name,"_chr", chr.nb,"_1Mbwind.png",sep = ""), width = 800, height = 400)
    # plot(dX, dY, type="l", main = paste("Recombination landscape (", method, ") of chromosome ", chr.nb, " in ", set.name, sep = ""), xlab = "Physical position (Mb)", ylab = "Recombination rate (cM/Mb)") #check
    # lines(dX, dY.upper, type = "l", lty = 2)
    # lines(dX, dY.lower, type = "l", lty = 2)
    # abline(h = 0, col = "Red", lty = 2) # Print a threshold for rec. rate = 0
    # dev.off()
    # rm(predicted.phys.1Mbwind)
    # rm(rec.phys.1Mbwind)
    
    # 100kb windows recombination map
    # Partitioning of the physical map in windows of 1Mb for further comparison with genomic content
    # phys.100kbwind = c(seq(0, max(chr.data$phys, na.rm = TRUE), by = 100), max(chr.data$phys, na.rm = TRUE)) # A sequence of windows coordinates
    # phys.100kbwind = as.numeric(rowMeans(embed(phys.100kbwind,2))) # Take the center of each window for interpolation
    # 
    # if (method == "loess") {
    #   #Integer overflow on large genomes (e.g. Gossypium chr 5), due to s.e. calculation
    #   # .Machine$integer.max
    #   # in a first time, do not calculate se for 100kbwind
    #   predicted.phys.100kbwind = data.frame(phys = phys.100kbwind, pred = predict(fit, newdata = as.numeric(phys.100kbwind), se = FALSE), se.pred = rep(NA, length(phys.100kbwind))) # Predict the local recombination rate at the exact physical positions of the Marey map, with a s.e.
    # }
    # if (method == "smooth.spline") {
    #   predicted.phys.100kbwind = data.frame(phys = phys.100kbwind, pred = predict(fit, x = phys.100kbwind)$y, se.pred = rep(NA, length(phys.100kbwind))) # Predict the local recombination rate at the exact physical positions of the Marey map
    # }
    # dY = diff(predicted.phys.100kbwind$pred)/diff(predicted.phys.100kbwind$phys)*1000 # the derivative of LOESS function + conversion in Mb (x 1000)
    # # Compute the 95% CI of the derivative
    # dY.upper = diff(predicted.phys.100kbwind$pred+1.96*predicted.phys.100kbwind$se.pred)/diff(predicted.phys.100kbwind$phys)*1000
    # dY.lower = diff(predicted.phys.100kbwind$pred-1.96*predicted.phys.100kbwind$se.pred)/diff(predicted.phys.100kbwind$phys)*1000
    # dX = rowMeans(embed(predicted.phys.100kbwind$phys,2)) # centers the X values for plotting
    # rec.phys.100kbwind = data.frame(phys = dX, rec.rate = dY, upper.ci = dY.upper,lower.ci = dY.lower)
    # write.table(rec.phys.100kbwind, file = paste("output/recombination_maps/", method,"/100kbwind/", set.name, "_chromosome", chr.nb, ".txt", sep = ""), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
    # # plot(rec.phys.100kbwind$rec.rate)
    # png(paste("figures/mareymap/",method, "/100kbwind/", set.name,"_chr", chr.nb,"_100kbwind.png",sep = ""), width = 800, height = 400)
    # plot(dX, dY, type="l", main = paste("Recombination landscape (", method, ") of chromosome ", chr.nb, " in ", set.name, sep = ""), xlab = "Physical position (bp) in windows of 100kb", ylab = "Recombination rate (cM/Mb)") #check
    # lines(dX, dY.upper, type = "l", lty = 2)
    # lines(dX, dY.lower, type = "l", lty = 2)
    # abline(h = 0, col = "Red", lty = 2) # Print a threshold for rec. rate = 0
    # dev.off()
    # rm(predicted.phys.100kbwind)
    # rm(rec.phys.100kbwind)
  }
}


# the fitting function for Smooth spline used in recombination.map
# Return a vector of fitted values of same size as original dataset (same number of markers)
fit.smooth.spline = function(chr.data = chr.data, spar = spar, adj.spar = adj.spar) {
  source("sources/stats.R")
  x = chr.data$phys
  y = chr.data$gen
  if (adj.spar == "auto") {
    spar2test = seq(from = 0.1, to = 0.5, len = 9) # Values of parameter to compare
    # A cross-validation data frame with
    # value of span to test
    # physical positions
    # observed genetic position
    # Predicted genetic positions
    size = length(spar2test)*nrow(chr.data) # Size of the data frame (number of rows)
    crossvalidation = data.frame(spar = rep(spar2test, each = nrow(chr.data)), phys = rep(chr.data$phys, length(spar2test)), obs = rep(chr.data$gen, length(spar2test)), pred = numeric(size))
    # First, subsetting the dataset in K clusters
    K = 5
    # Systematic K-fold selection, points are uniformly sampled along the map (a sequence of K-th nearest neighbors)
    # i.e., for K = 5, sample 1, 6, 11, 17... then 2, 7, 12, 18...
    subset = rep(seq(1,5), ceiling(nrow(chr.data)/5)) # Vector of clusters' number
    subset = subset[1:nrow(chr.data)] # Trim the vector to its maximal size
    
    for (s in spar2test) {
      # Run K times the cross-validation
      # Each time the cluster K is leaved out for fitting and used as validation
      # results.cv = data.frame(K = 1:K, ACV = numeric(K), LSCV = numeric(K))
      for (k in 1:K) {
        training = rbind(chr.data[1,], subset(chr.data, subset %in% seq(1,5)[-k]), chr.data[nrow(chr.data),]) # Only data out of cluster K + boundaries of the map
        validation = subset(chr.data, subset == k)
        
        fit.cv = smooth.spline(x = training$phys, y = training$gen, spar = s) # Same defaults parameters as in MareyMap
        summary(fit.cv)
        # Predict the local recombination rate at the exact physical positions of the Marey map
        crossvalidation$pred[which(crossvalidation$spar == s)[subset == k]] = predict(fit.cv, x = validation$phys)$y
      }
      
      # Compute the criterion for each spar value
      # source("sources/stats.R")
      crit = data.frame(spar = spar2test, ACV = numeric(length(spar2test)), LSCV = numeric(length(spar2test)))
      for (s in spar2test) {
        # Absolute Cross-validation
        # Compute the absolute value of the difference between the mean observed value and the predicted value in the left out dataset
        # lambda value to optimize (the span parameter)
        # Input is a data frame with three columns of physical position and observed/predicted values
        crit$ACV[crit$spar == s] = ACV(data = crossvalidation[which(crossvalidation$spar == s), -1])
        
        # Least squares Cross validation
        crit$LSCV[crit$spar == s] = LSCV(data = crossvalidation[which(crossvalidation$spar == s), -1])
      }
      
      # Automated selection of the span value which minimizes the chosen criterion
      spar = crit$spar[crit$ACV == min(crit$ACV, na.rm = TRUE)]
    }
  } else {
    if (is.numeric(spar)) {
      spar = spar
    } else {
      warning("The argument 'spar' must be 'auto' or a numeric value between 0 and 1.\n")
    }
  }
  
  #-----------------------------------------
  # FITTING
  fit = smooth.spline(x = chr.data$phys, y = chr.data$gen, spar = spar) # Same defaults parameters as in MareyMap
  # Save actions in log with parameters
  # print(set)
  save.log(path = "/output/", msg = paste("Cubic smooth spline fitting on chromosome", unique(chr.data$map), "with graphical predictions saved as figure", sep = " "), append = TRUE)
  
  return(fit)
}
