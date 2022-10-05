# Work with Marey maps

#============================================================================#
# Check orientation of the map ----
#============================================================================#
# Check & correct bad orientation of linkage groups
# Markers oriented in the direction of a reference map if necessary
# A reference genome with markers mapped onto it seems inevitable for an exhaustive check
# set = name of the dataset, ONE map
# chr = chromosome number to flip entirely, automatic check if chr = c()
# dir = path to the dataset

# Yet, the Marey map curve can only be ascending
check_map_orientation = function(set = "", id_chr = c(), dir = "data/Marey_maps/", auto = FALSE) {
    marey = read.table(paste(dir , set, sep = ""), header = TRUE)
    marey$map = as.character(marey$map)
    marey$gen = as.numeric(marey$gen)
    # The function gen~phys must be monotonously increasing
    # In case you search to identify segments going the wrong direction in a chromosome
    if (auto == TRUE) {
      new_map = data.frame()
      for (chr in unique(marey$map)) {
        # Find groups of markers with descending genetic position
        # (i.e. genetic and physical maps not oriented in the same direction)
        flip_group = c()
        #---------------------------------------------------#
        # Exact approach
        #---------------------------------------------------#
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
        #---------------------------------------------------#
        # Likelihood approach
        #---------------------------------------------------#
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
        new_map = marey[-which(marey$map %in% id_chr),] # Unchanged chromosome are saved in new_map before procedure
      for (chr in id_chr) {
        subset = marey[marey$map == chr,]
        subset$gen = as.numeric(subset$gen)
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

#============================================================================#
# Flip orientation of a segment by manual selection ----
#============================================================================#
# Some segments on the Marey map have the wrong orientation
# This function allows the manual selection of the segment
# that is flipped in the good orientation

# USAGE
# flip.segment(x = data, set ="", chr = "")

# ARGUMENTS
# set =  the name of the dataset
# chr = the name of the chromosome

# VALUES
# Return the same data frame as input, with genetic positions within the segment flipped
flip.segment = function(set ="", chr = "") {
  # x = a data frame of form marey map
  x = read.table(paste("data/Marey_maps/", set, sep=""), header = TRUE)
  x$map = as.character(x$map)
  tmp = x[as.character(x$map) == chr,]
  # Select points in Marey plot
  idx = points.selection(tmp)
  # A vector of genetic distances values to work on
  vec = x$gen[row.names(x) %in% as.character(idx)]
  # Flip the selection selected segment
  new_gendist = abs(vec - max(vec, na.rm = TRUE)) + min(vec, na.rm = TRUE)
  # Reaffect values to the map
  x$gen[row.names(x) %in% as.character(idx)] = new_gendist
  return(x)
}


# EXAMPLE
# set = "Raphanus_sativus_Luo2019.txt"
# chr = "1"
# marey.diagplot(set = set, dir = "data/Marey_maps/", save.plot = FALSE, display = TRUE)
# new_map = flip.segment(set, chr)
# marey.diagplot(set = new_map, save.plot = FALSE)


#============================================================================#
# Filter out outliers automatically
#============================================================================#
# Selection of valid markers

# How to filter outliers?
# 1/ Fit a curve and filter markers away from the curve (threshold):
# bootstrap the dist. of the mean departure from the curve, exclude points far away than 95% (2.5% each side)
# Advantage, selection is curve specific
# 2/ Sample the subset of markers minimizing the RMSE (resampling approach)



#============================================================================#
# Filter out outliers automatically ----
#============================================================================#



#============================================================================#
# Filter out outliers manually ----
#============================================================================#
# Exclude some markers from further analyses by a manual selection in a Marey plot (graphical selection)
# Draw a selection around markers to exclude

# USAGE
# outliers.selection(x = data, set ="", chr = "")

# ARGUMENTS
# x = a data frame of form marey map
# set =  the name of the dataset
# chr = the name of the chromosome

# VALUES
# Return the same data frame as input, with excluded markers flagged with vld = FALSE

outliers.selection = function(set ="", chr = "") {
  # x = a data frame of form marey map
  x = read.table(paste("data/Marey_maps/", set, sep=""), header = TRUE)
  x$phys = as.numeric(as.character(x$phys))
  tmp = x[as.character(x$map) == chr,]
  # Select points in Marey plot
  idx = points.selection(tmp)
  # Flag selected points as FALSE
  x$vld[as.character(row.names(x)) %in% idx] = FALSE
  return(x)
}

# EXAMPLE
# tmp = outliers.selection(x = data.final, set = "Arabidopsis_thaliana_Serin2017", chr = "1")

#============================================================================#
# Points selection
#============================================================================#
# Plot a Marey map and draw a selection around points you want to select

# USAGE
# points.selection(x = data)

# ARGUMENTS
# x = a data frame of form marey map

# VALUES
# Return a vector of row names for selected points

points.selection = function(x) {
  require(gatepoints)
  df = x[,c(4:5)]
  df = df[!is.na(df$phys),]
  # Display the Marey map
  X11()
  plot(df, col = "red")
    # Draw the selection around points
  selectedPoints = fhs(df, mark = TRUE)
  dev.off()
  # Return a vector of row names for selected points
  points.names = as.character(selectedPoints)
  return(points.names)
}



#============================================================================#
# Discard a chromosome
#============================================================================#


discard.map = function(set ="", chr = "") {
  # x = a data frame of form marey map
  x = read.table(paste("data/Marey_maps/", set, sep=""), header = TRUE)
  x$vld[x$set == gsub(".txt", "", set) & x$map == chr] = FALSE
  write.table(x, file = paste("data/Marey_maps/", set, sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
  save.log()
}

# EXAMPLE
# discard.map(set, chr)

#============================================================================#
# Filter out chromosomes with insufficient number of markers,
#============================================================================#
# based on an interval threshold or number of markers specified by the user
# mean interval distance



#============================================================================#
# Diagnostic plots
#============================================================================#
# Genetic distance (cM) as a function of physical position (bp)

# set (character) is a vector of names of Marey map files to plot, or a data.frame object
# dir (character) is the path of datasets
# output (character) is the path where to save plots
# save = TRUE if you want to save your plot as a png file
# display = TRUE if you want to print the plot in Rstudio

marey.diagplot = function(set = "", dir = "data/Marey_maps/", output = "output/marey_maps/diagnostic_plots/", save.plot = TRUE, display = TRUE) {
  require(ggplot2)
  if (is.data.frame(set)) {
    marey = set
    marey$map = as.character(marey$map)
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
        marey = read.table(paste(dir, i, sep = ""), header = TRUE, sep = "\t")
        print(i)
        marey$map = as.character(marey$map)
        marey = marey[!(is.na(marey$phys) | is.na(marey$gen) | is.na(marey$map)),]
        marey$phys = as.numeric(as.character(marey$phys))
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



#============================================================================#
# Map a physical map onto a genetic map
#============================================================================#
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
    phys_map$map = as.character(phys_map$map)
    gen_map$Chr = as.character(gen_map$Chr)
    #-----------------------------------------#
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
    marey = data.frame(set = rep(gsub(".txt", "", set), length(map)), map = map, mkr = gen_map$Marker, phys = match,
                       gen = gen_map$cM, vld = TRUE)
    colnames(marey) = c("set", "map", "mkr", "phys", "gen", "vld")
    # Unvalidate markers with NA in physical positions
    marey$vld[is.na(marey$phys)] = FALSE
    
    
    
    #-----------------------------------------#
    # Map chromosome number for each marker
    # When the information is not available, use linkage group (LG) instead of chromosome number in genetic maps
    # Map the chromosome number afterward onto the Marey_map input file
    
    # Which genetic maps have LG numbers (i.e. LG[0-9]* in the column 'Chr')?
    
    
    #-----------------------------------------#
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



#============================================================================#
# Compute the number of markers per genetic map
#============================================================================#
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

#============================================================================#
# Compute the linkage map length
#============================================================================#
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


#============================================================================#
# Estimation of the local recombination rate (cM/Mb)
#============================================================================#
recombination.map = function(set ="", data = data.final, chr = "all", method = c("loess", "smooth.spline"), K = 5, boot = 1000, parallel = TRUE, span = NA,...) { # Optional arguments passed to loess, smooth.spline and pther fitting methods
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
    
    #-----------------------------------------#
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

      #-----------------------------------------#
      # LOESS
      #-----------------------------------------#
      if (method == "loess") {
        #-------------------------------------#
        # Calibrate span
        if (is.na(span)) {
          span = calibrate.span(x = chr.data, K, from = 0.2, to = 0.5, parallel = parallel)
        } else {
          span = span
        }
        
        #-------------------------------------#
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
        #-----------------------------------------#
        # CUBIC SMOOTH SPLINE
        #-----------------------------------------#
        if (method == "smooth.spline") {
          #-------------------------------------#
          # Calibrate span
          spar = calibrate.spar(data = chr.data, K)
          
          #-------------------------------------#
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
          #-------------------------------------#
          # Veller method: interpolate 1,000 evenly distributed physical positions
          #-------------------------------------#
          if (method == "veller") {
            #-------------------------------------#
            # Calibrate span
            span = calibrate.span(x = chr.data, K)
            # span = 0.3
            #-------------------------------------#
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
        #-----------------------------------------#
        # DERIVATIVES
        # The local recombination rate is the derivative of the polynom (i.e. fitted loess) in a point
        #------------------------------#
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
        
        #------------------------------#
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
        
        #------------------------------#
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
        
        #------------------------------#
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
        
        #-----------------------------------------#
        # LOESS
        #-----------------------------------------#
        if (method == "loess") {
          #-------------------------------------#
          # Calibrate span
          if (is.na(span)) {
            span = calibrate.span(x = chr.data, K, from = 0.2, to = 0.5, parallel = parallel)
          } else {
            span = span
          }
          # if (is.na(span)) {span = 0.3}
          # span = 0.3
          #-------------------------------------#
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
          #-----------------------------------------#
          # CUBIC SMOOTH SPLINE
          #-----------------------------------------#
          if (method == "smooth.spline") {
            #-------------------------------------#
            # Calibrate span
            spar = calibrate.spar(data = chr.data, K)
            
            #-------------------------------------#
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
        
        #-----------------------------------------#
        # DERIVATIVES
        # The local recombination rate is the derivative of the polynom (i.e. fitted loess) in a point
        #------------------------------#
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
        plot(dX, dY, type="l", ylim = c(0, 30), main = paste("Recombination landscape (", method, ") of chromosome ", chr.nb, " in ", set.name, sep = ""), xlab = "Physical position (Mb)", ylab = "Recombination rate (cM/Mb)") #check
        lines(dX, dY.upper, type = "l", lty = 2, col = "Grey")
        lines(dX, dY.lower, type = "l", lty = 2, col = "Grey")
        abline(h = 0, col = "Red", lty = 2) # Print a threshold for rec. rate = 0
        dev.off()
        rm(dX, dY, dY.upper, dY.lower)
        rm(pointwise.boot)
        
        #------------------------------#
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
        plot(dX, dY, type="l", ylim = c(0, 30), main = paste("Recombination landscape (", method, ") of chromosome ", chr.nb, " in ", set.name, sep = ""), xlab = "Physical position (Mb)", ylab = "Recombination rate (cM/Mb)") #check
        lines(dX, dY.upper, type = "l", lty = 2, col = "Grey")
        lines(dX, dY.lower, type = "l", lty = 2, col = "Grey")
        abline(h = 0, col = "Red", lty = 2) # Print a threshold for rec. rate = 0
        dev.off()
        rm(dX, dY, dY.upper, dY.lower)
        rm(physbins.boot)
        
        #------------------------------#
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
        plot(dX, dY, type="l", ylim = c(0, 30), main = paste("Recombination landscape (", method, ") of chromosome ", chr.nb, " in ", set.name, sep = ""), xlab = "Physical position (Mb)", ylab = "Recombination rate (cM/Mb)") #check
        lines(dX, dY.upper, type = "l", lty = 2, col = "Grey")
        lines(dX, dY.lower, type = "l", lty = 2, col = "Grey")
        abline(h = 0, col = "Red", lty = 2) # Print a threshold for rec. rate = 0
        dev.off()
        rm(dX, dY, dY.upper, dY.lower)
        rm(phys.1Mbwind.boot)
        
        #------------------------------#
        # 100kb windows recombination map
        # dY = rowMeans(phys.100kbwind.boot, na.rm = TRUE) # the derivative of LOESS function
        dY = apply(phys.100kbwind.boot, MARGIN = 1, function(x) mean(x, na.rm = TRUE)) # the derivative of LOESS function
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
        plot(dX, dY, type="l", ylim = c(0, 30), main = paste("Recombination landscape (", method, ") of chromosome ", chr.nb, " in ", set.name, sep = ""), xlab = "Physical position (Mb)", ylab = "Recombination rate (cM/Mb)") #check
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
    
    #-----------------------------------------#
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


#-----------------------------------------#
# LOESS span parameter calibration
#-----------------------------------------#
# Automatic calibration of the span parameter
# TODO Be careful with this function, that seems to converge too often to the same lower bounds of tested values.
# I need more work on it to refine the process and to be more stringent to reject overfitting
# Actually, overfitting is not always reduced if initial values of span include values less than 0.2
# i.e. in the 'span2test' list
# The same remark is valid for the 'spar' parameter in smooth spline interpolation (i.e. 'calibrate.spar' function)
# See https://stats.stackexchange.com/questions/2002/how-do-i-decide-what-span-to-use-in-loess-regression-in-r
calibrate.span = function(x, K = 5, from = 0.2, to = 0.5, parallel = TRUE) {
  source("sources/stats.R")
  
  # For very large datasets (>10,000 markers per chromosome), need to sample a smaller fraction (10%)
  if (nrow(x) > 10000) {
    cat("Dataset too large (>10,000 markers). Subsampling 10%...\n")
    x = x[sample(1:nrow(x), size = nrow(x)/10, replace = FALSE),]
  }
  
  # The resampling function to parallelize
  cv.resampling = function() {
    # Resample nresampling times in training and validation sets of proportion (k-1)/k and 1/k respectively
    # Sample randomly training and validation
    idx = sample(1:nrow(x), size = nrow(x)*((K-1)/K), replace = FALSE) # Index of individuals to keep
    training = x[idx,]
    validation = x[-idx,]
    # Fit the curve on training
    fit.loess.cv = loess(gen~phys, training, span = s) # Same defaults parameters as in MareyMap
    summary(fit.loess.cv)
    # Predict the local recombination rate at the exact physical positions of the Marey map on validation
    predict = predict(fit.loess.cv, newdata = validation$phys, se = TRUE)
    validation$pred = predict$fit
    # Compute the criterion for validation data
    # Input is a data frame with three columns of physical position and observed/predicted values
    df = validation[,c(4,5,7)]
    colnames(df) = c("phys", "obs", "pred")
    # ACV.boot[boot] = ACV(data = df)
    # LSCV.boot[boot] = LSCV(data = df)
    # SECV.boot[boot] = sd(predict$se.fit, na.rm = TRUE)
    mean_diff = mean((df$obs - df$pred)^2, na.rm = TRUE)
    rm(df)
    return(mean_diff)
    # END OF RESAMPLING
  }

  span2test = seq(from = from, to = to, by = 0.05) # Values of parameter to compare
  # If number of individuals too small, minimal span to test need to be increased or K decreased
  if (nrow(x) < 90) {
    span2test = seq(from = max(0.3, from), to = to, by = 0.05)
  }
  # A cross-validation data frame with
  # value of span to test and associated criterion after resampling
  size = length(span2test)*nrow(x) # Size of the data frame (number of rows)
  crossvalidation = data.frame(span = span2test, ACV = NA, LSCV = NA, MSE = NA)
  
  # crossvalidation = data.frame(span = rep(span2test, each = nrow(x)), phys = rep(x$phys, length(span2test)), obs = rep(x$gen, length(span2test)), pred = numeric(size))
  # First, subsetting the dataset in K clusters
  # Systematic K-fold selection, points are uniformly sampled along the map (a sequence of K-th nearest neighbors)
  # i.e., for K = 5, sample 1, 6, 11, 17... then 2, 7, 12, 18...
  # subset = rep(seq(1,5), ceiling(nrow(x)/5)) # Vector of clusters' number
  # subset = subset[1:nrow(x)] # Trim the vector to its maximal size
  
  # Systematic K-fold cross validation
  # for (s in span2test) {
  #   # Run K times the cross-validation
  #   # Each time the cluster K is leaved out for fitting and used as validation
  #   # results.cv = data.frame(K = 1:K, ACV = numeric(K), LSCV = numeric(K))
  #   for (k in 1:K) {
  #     training = rbind(x[1,], subset(x, subset %in% seq(1,5)[-k]), x[nrow(x),]) # Only data out of cluster K + boundaries of the map
  #     validation = subset(x, subset == k)
  #     
  #     fit.loess.cv = loess(gen~phys, training, span = s) # Same defaults parameters as in MareyMap
  #     summary(fit.loess.cv)
  #     # Predict the local recombination rate at the exact physical positions of the Marey map
  #     crossvalidation$pred[which(crossvalidation$span == s)[subset == k]] = predict(fit.loess.cv, newdata = validation$phys)
  #     # validation.loess = data.frame(pos = validation$phys, obs = validation$gen, pred = predict(fit.loess.cv, newdata = validation$phys)) # Predict the local recombination rate at the exact physical positions of the Marey map
  #     # validation.loess
  #   }
  # }
  # Use random resampling cross-validation to improve calibration and decrease sensitivity to missing data/outliers
  for (s in span2test) {
    # Run K times the cross-validation
    # Each time the cluster K is leaved out for fitting and used as validation
    # results.cv = data.frame(K = 1:K, ACV = numeric(K), LSCV = numeric(K))
    nresampling = 1000
    
    # Initiate bootstrap vectors of criterions
    # ACV.boot = numeric(nresampling) # Absolute differences
    # LSCV.boot = numeric(nresampling) # Least Squares differences
    # SECV.boot = numeric(nresampling) # Mean Standard Error
    MSE.boot = numeric(nresampling) # Mean Squared Error
    
    # Parallelize bootstraps for performance issues
    # Yet mclapply can have big memory issues on large dataset
    cat("Calibrating span", s, "...\n")
    if (parallel == TRUE) {
      # require(parallel)
      require(pbmcapply)
      MSE.boot = unlist(pbmclapply(X = MSE.boot, function(x) cv.resampling(), mc.cores = detectCores(all.tests = FALSE, logical = TRUE)-1))
    } else {
      MSE.boot = unlist(lapply(X = MSE.boot, function(x) cv.resampling()))
    }
    
        # # if (parallel == TRUE) {
        # if (nrow(x) < 10000) {
        #   require(parallel)
        #   MSE.boot = unlist(mclapply(X = MSE.boot, function(x) cv.resampling(), mc.cores = 7))
        # } else {
        #   require(parallel)
        #   MSE.boot = unlist(mclapply(X = MSE.boot, function(x) cv.resampling(), mc.cores = 7))
        #   # require(parallel)
        #   # Reduce the sampling to prevent memory overflow
        #   MSE.boot = numeric(10)
        #   # MSE.boot = unlist(mclapply(X = MSE.boot, function(x) cv.resampling(), mc.cores = 7))
        #   # for (i in 1:length(MSE.boot)) {
        #   #   cat("Boot", i, "...")
        #   #     MSE.boot[i] = cv.resampling()
        #   # }
        #   MSE.boot = unlist(lapply(X = MSE.boot, function(x) cv.resampling()))
        # }

    # # Resample nresampling times in training and validation sets of proportion (k-1)/k and 1/k respectively
    # for (boot in 1:nresampling) {
    #   # Sample randomly training and validation
    #   idx = sample(1:nrow(x), size = nrow(x)*((K-1)/K), replace = FALSE) # Index of individuals to keep
    #   training = x[idx,]
    #   validation = x[-idx,]
    #   # Fit the curve on training
    #   fit.loess.cv = loess(gen~phys, training, span = s) # Same defaults parameters as in MareyMap
    #   summary(fit.loess.cv)
    #   # Predict the local recombination rate at the exact physical positions of the Marey map on validation
    #   predict = predict(fit.loess.cv, newdata = validation$phys, se = TRUE)
    #   validation$pred = predict$fit
    #   # Compute the criterion for validation data
    #   # Input is a data frame with three columns of physical position and observed/predicted values
    #   df = validation[,c(4,5,7)]
    #   colnames(df) = c("phys", "obs", "pred")
    #   # ACV.boot[boot] = ACV(data = df)
    #   # LSCV.boot[boot] = LSCV(data = df)
    #   # SECV.boot[boot] = sd(predict$se.fit, na.rm = TRUE)
    #   MSE.boot[boot] = mean((df$obs - df$pred)^2, na.rm = TRUE)
    #   rm(df)
    #   # END OF RESAMPLING
    # }
    # crossvalidation$ACV[which(crossvalidation$span == s)] = mean(ACV.boot)
    # crossvalidation$LSCV[which(crossvalidation$span == s)] = mean(LSCV.boot)
    # crossvalidation$SECV[which(crossvalidation$span == s)] = mean(SECV.boot)
    crossvalidation$MSE[which(crossvalidation$span == s)] = mean(MSE.boot, na.rm = TRUE)
    # END OF LOOP SPAN2TEST
  }  
  
  # # Compute the criterion for each span
  # # source("sources/stats.R")
  # crit = data.frame(span = span2test, ACV = numeric(length(span2test)), LSCV = numeric(length(span2test)))
  # for (s in span2test) {
  #   # Absolute Cross-validation
  #   # Compute the absolute value of the difference between the mean observed value and the predicted value in the left out dataset
  #   # lambda value to optimize (the span parameter)
  #   # Input is a data frame with three columns of physical position and observed/predicted values
  #   crit$ACV[crit$span == s] = ACV(data = crossvalidation[which(crossvalidation$span == s), -1])
  #   
  #   # Least squares Cross validation
  #   crit$LSCV[crit$span == s] = LSCV(data = crossvalidation[which(crossvalidation$span == s), -1])
  # }
  
  # Automated selection of the span value which minimizes the chosen criterion
  # Here MSE, the Mean Squared Error
  span = crossvalidation$span[crossvalidation$MSE == min(crossvalidation$MSE, na.rm = TRUE)]
  
  # Save the span parameter and the criterion in output file "output/recombination_maps/span.txt"
  output_span = read.table("output/recombination_maps/span.txt", header = TRUE)
  idx = which(unique(x$set) == output_span$set & unique(x$map) == output_span$map)
  # Check first if the span parameter is already in the dataframe
  if (length(idx) > 0) {
    output_span$span[idx] = span
    output_span$crit[idx] = min(crossvalidation$MSE, na.rm = TRUE)
    output_span$K[idx] = K
  } else { # Add a new line if the span parameter for this map have never been saved
    output_span = rbind(output_span,
                        data.frame(set =  unique(x$set),
                                   map =  unique(x$map),
                                   span = span,
                                   crit = min(crossvalidation$MSE, na.rm = TRUE),
                                   K = K))
  }
  write.table(output_span, file = "output/recombination_maps/span.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  return(span)
}

#-----------------------------------------#
# LOESS fitting function
#-----------------------------------------#
fit.loess = function(x = chr.data, span = span, degree = degree, K = K) {
  
  #-----------------------------------------#
  # FITTING
  fit = loess(gen~phys, x, span = span) # Same defaults parameters as in MareyMap
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
  #save.log(path = "/output/", msg = paste("LOESS fitting on chromosome", unique(x$map), ", with LOESS graphical predictions saved as figure", sep = " "), append = TRUE)
  
  return(fit)
}

#-----------------------------------------#
# CUBIC SMOOTH SPLINE spar parameter calibration
#-----------------------------------------#
calibrate.spar = function(data, K = 5) {
  source("sources/stats.R")
  # For very large datasets (>10,000 markers per chromosome), need to sample a smaller fraction (10%)
  if (nrow(data) > 10000) {
    cat("Dataset too large (>10,000 markers). Subsampling 10%...\n")
    data = data[sample(1:nrow(data), size = nrow(data)/10, replace = FALSE),]
  }
  
  x = data$phys
  y = data$gen
  

  
  #------------------------------------------#
  # Parameter Spar calibration
  #------------------------------------------#
  spar2test = seq(from = 0.2, to = 0.5, by = 0.05) # Values of parameter to compare
  # A cross-validation data frame with
  # value of span to test
  # physical positions
  # observed genetic position
  # Predicted genetic positions
  size = length(spar2test)*nrow(data) # Size of the data frame (number of rows)
  crossvalidation = data.frame(spar = rep(spar2test, each = nrow(data)), phys = rep(data$phys, length(spar2test)), obs = rep(data$gen, length(spar2test)), pred = numeric(size))
  # First, subsetting the dataset in K clusters
  # K = 5
  # Systematic K-fold selection, points are uniformly sampled along the map (a sequence of K-th nearest neighbors)
  # i.e., for K = 5, sample 1, 6, 11, 17... then 2, 7, 12, 18...
  subset = rep(seq(1,5), ceiling(nrow(data)/5)) # Vector of clusters' number
  subset = subset[1:nrow(data)] # Trim the vector to its maximal size
  
  for (s in spar2test) {
    # Run K times the cross-validation
    # Each time the cluster K is leaved out for fitting and used as validation
    # results.cv = data.frame(K = 1:K, ACV = numeric(K), LSCV = numeric(K))
    for (k in 1:K) {
      training = rbind(data[1,], subset(data, subset %in% seq(1,5)[-k]), data[nrow(data),]) # Only data out of cluster K + boundaries of the map
      validation = subset(data, subset == k)
      
      fit.cv = smooth.spline(x = training$phys, y = training$gen, spar = s) # Same defaults parameters as in MareyMap
      summary(fit.cv)
      # Predict the local recombination rate at the exact physical positions of the Marey map
      crossvalidation$pred[which(crossvalidation$spar == s)[subset == k]] = predict(fit.cv, x = validation$phys)$y
    }
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
  
  return(spar)
}

#-----------------------------------------#
# CUBIC SMOOTH SPLINE fitting function
#-----------------------------------------#
# the fitting function for Smooth spline used in recombination.map
# Return a vector of fitted values of same size as original dataset (same number of markers)
fit.smooth.spline = function(x = chr.data, spar = spar, K = K) {
  
  #-----------------------------------------#
  # FITTING
  fit = smooth.spline(x = x$phys, y = x$gen, spar = spar) # Same defaults parameters as in MareyMap
  # Save actions in log with parameters
  # print(set)
  # save.log(path = "/output/", msg = paste("Cubic smooth spline fitting on chromosome", unique(x$map), "with graphical predictions saved as figure", sep = " "), append = TRUE)
  
  return(fit)
}


#----------------------------------------------------------------------#
# Natural splines method of Berloff (2002) with penalized max likelihood estimation
#----------------------------------------------------------------------#
# Custom implementation of the penalized max likelihood estimation with natural splines developped by Berloff et al. (2002)

smooth.spline.berloff = function()  {
  
}


#============================================================================#
# Estimation of the genetic shuffling in terms of gene distance
#============================================================================#
GeneDistancesShuffling = function(set ="", chr = "all", K = 3,...) { # Optional arguments passed to loess, smooth.spline and other fitting methods
  cat("##################################################################\nDataset", set, "\n")
  # Gene distances are in cumulated number of genes
  data = read.table(gzfile(paste("data-cleaned/genome/gene_distance/", set, ".txt.gz", sep = "")), header = TRUE,
                    stringsAsFactors = FALSE)
  
  if (chr == "all") { # Chromosomes to map are given in a list of names of "all" to evaluate every chromosome in the dataset
    chr.list = unique(data$map)
  } else {
    chr.list = chr
  }
  
  # Sometimes NA or unwanted scaffold in data as chromosomes
  chr.list = chr.list[!is.na(chr.list)]
  chr.list = chr.list[!grepl("scaffold", chr.list)]
  
  # Discard manually some chromosomes were loess is not working
  # No explanation yet...
  if (set == "Helianthus_annuus_Talukder2014") {chr.list = chr.list[chr.list != 7]}
  
  # For each chosen chromosome in a dataset
  for (c in chr.list) {
    cat("Treating chromosome ", c, "...\n")
    # Regression is made on a whole chromosome
    chr.nb = c # Chromosome number
    data$set = as.character(data$set)
    data$map = as.character(data$map)
    chr.data = subset(data, map == chr.nb & set == set) # Sample only distances of the chromosome
    chr.data$gen = as.numeric(as.character(chr.data$gen))
    # Physical positions are replaced by gene position
    chr.data$phys = as.numeric(as.character(chr.data$gene_distance))

    # Remove NA distances
    chr.data = chr.data[!is.na(chr.data$gen),]
    chr.data = chr.data[!is.na(chr.data$phys),]
    
    degree = 2
    
    # If the gene count failed and there is no gene annotated
    if ((max(chr.data$phys, na.rm = TRUE) == 0) | sum(!is.na(chr.data$phys)) == 0) {
      warning("Gene count failed: no genes have been annotated and counted.")
    } else {
      # Partitioning of the physical map (bp): 1,000 bins of equal size
      phys.bins = seq(0, max(chr.data$phys, na.rm = TRUE), len = 1000) # A sequence of bins coordinates
      #-------------------------------------#
      # Calibrate span
      # span = calibrate.span(x = chr.data, K)
      # # span = 0.3
      # span = span[-which(is.na(span))]
      
      # Re-use the already inferred span
      span = read.table("output/recombination_maps/span.txt", header = TRUE)
      span = span$span[which(span$set == set & span$map == c)]
      if (length(span) == 0) {span = 0.4}
      if (is.na(span)) {span = 0.4}
      
      #-------------------------------------#
      out = try(
        {# Fit the interpolation curve on data
          fit = fit.loess(x = chr.data, span = span, degree = degree)      # Partitioning of the physical map (bp): 1,000 bins of equal size
          predicted.physbins = predict(fit, newdata = phys.bins) # Predict the genetic distance at the exact physical positions of the Marey map
        }
      )

      # A Marey map of 1,000 evenly spaced markers
      df.tmp = read.table("output/veller/GeneDistances.txt", header = TRUE, stringsAsFactors = FALSE)
      # Return a data.frame
      df.veller = data.frame(set = rep(set, 1000), map = rep(chr.nb, 1000), mkr = paste("chr", chr.nb,"_",1:1000, sep = ""), phys = phys.bins, gen = predicted.physbins)
      # Append new data to replace
      df.tmp = df.tmp[which(!(df.tmp$set == set & df.tmp$map == chr.nb)),]
      df.tmp = rbind(df.tmp, df.veller)
      write.table(df.tmp, file = "output/veller/GeneDistances.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
    }

    
  }
  }




#============================================================================#
# BROCKEN STICK MODEL
#============================================================================#

# Estimate the proportions of a broken stick model for a single dataset (all chromosomes)
# i.e. proportion of total physical length in k segments of equal genetic length along the chromosome

# ARGUMENTS
# set = "", name of the dataset to import in the 'data-cleaned/marey_maps' directory
# k = 3, number of segments to which the map has to be decomposed
# method = c("strict", "segmented"). With method "strict", segements are strictly of equal genetic size
# with method "segmented", break points in the piecewise-linear relationshiop are estimated with a likelihood approach implemented in the R package 'segmented'

# VALUE
# Return a data frame with k columns (k1, k2,... kn), with proportions of total physical length in k segments
# with chromosomes in lines

# Load a map: Arabidopsis and Zea to compare contrasting patterns
# set = "Arabidopsis_thaliana_Serin2017"
# set = "Zea_mays_IBM_MaizeSNP50"

# DETAILS
# Broken-stick model (White & Hill 2020)
# K segments of equal physical length (Mb) but different genetic length (cM)
# A model with three segments seemed sufficient to describe most patterns of recombination rate/most landscapes (broad-scale variation)
# Show heterogeneity telomere/center
# Show assymetry

brokenstick = function(set = "", k = 3, method = "strict") {
  map = read.table(paste("data-cleaned/marey_maps/", set, ".txt", sep =""), header = TRUE)
  map$map = as.character(map$map)
  map$mkr = as.character(map$mkr)
  map$phys = as.numeric(as.character(map$phys))
  # List of chromosomes to process
  listchr = unique(map$map)
  # Remove NA chromosomes
  listchr = listchr[!is.na(listchr)]
  # Creat empty dataset
  stick_proportions = as.data.frame(matrix(NA, nrow = length(listchr), ncol = k))
  colnames(stick_proportions) = paste("p", c(1:k), sep ="")
  
  # Proportions for each chromosome
  count = 0
  for (chr in listchr) {
    print(chr)
    # Select the chromosome
    chrmap = map[map$map == chr,]
    # Remove NAs in genetic distances
    chrmap = subset(chrmap, !is.na(chrmap$gen))
    count = count + 1
    
    if (method == "strict") {
      # Segment in k segments of equal genetic size (cM)
      # Size of segment is total genetic size divided by k
      segment_size = round(max(chrmap$gen, na.rm = TRUE)/k)
      
      # A vector of segments physical size pi (proportion of total length)
      p = numeric(k)
      for (i in 1:k) {
        # print(max(which(chrmap$gen < segment_size*i)))
        # p[i] = chrmap$phys[max(which(chrmap$gen < segment_size*i))]
        p[i] = max(chrmap$phys[which(chrmap$gen < segment_size*i)], na.rm = TRUE)
        # p[i] = max(chrmap$gen[chrmap$gen < segment_size*i])
        # p[i] = chrmap$phys[max(which(chrmap$gen < quantile(chrmap$gen, i/k)))]
        # p[i] = quantile(chrmap$gen, i/k)
        
        if (i > 1) {
          p[i] = p[i] - sum(p[i-1:(i-1)])
        }
      }
      # END OF METHOD STRICT
    }
    if (method == "segmented") {
      library(segmented)
      res = c()
      lin.mod <- lm(phys~gen, data = chrmap)
      # Estimated breakpoints in genetic distances
      res = try(segmented.mod <- segmented(lin.mod, seg.Z = ~gen, npsi=(k-1), control = seg.control(n.boot = 100, fix.npsi=TRUE)))
      # In cases of segmentation failure
      if (class(res)[1] == "try-error") {
        # Failure because not enough data
        breakpoints = c(NA, NA, NA)
      } else {
        breakpoints = c(summary(segmented.mod)$psi[,2], max(chrmap$gen, na.rm = TRUE))
      }
      
      for (i in 1:k) {
        p[i] = max(chrmap$phys[chrmap$gen < breakpoints[i]])

        if (i > 1) {
          p[i] = p[i] - sum(p[i-1:(i-1)])
        }
      }
      # END OF METHOD SEGMENTED
    }
    
    # p is a proportion of total physical length
    p = p/sum(p)
    
    stick_proportions[count,] = rbind(p)
    # END OF CHROMOSOME LOOP
  }
  stick_proportions = cbind(listchr, stick_proportions)
  colnames(stick_proportions)[1] = "chromosome"
  stick_proportions
  
  return(stick_proportions)
}

batch_brokenstick = function(set = "", k = 10) {
  print(set)
  if (!file.exists(paste("output/brokenstick/brokenstick_k", k,".txt", sep = ""))) {
    if (k == 10) {
      df = data.frame(set = character(),	chromosome = character(),	p1 = numeric(),	p2 = numeric(),	p3 = numeric(),
                      p4 = numeric(),	p5 = numeric(),	p6 = numeric(),	p7 = numeric(),	p8 = numeric(),	p9 = numeric(),	p10 = numeric())
    }
    if (k == 20) {
      df = data.frame(set = character(),	chromosome = character(),	p1 = numeric(),	p2 = numeric(),	p3 = numeric(),
                      p4 = numeric(),	p5 = numeric(),	p6 = numeric(),	p7 = numeric(),	p8 = numeric(),	p9 = numeric(),	p10 = numeric(),
                      p11 = numeric(),	p12 = numeric(),	p13 = numeric(),
                      p14 = numeric(),	p15 = numeric(),	p16 = numeric(),	p17 = numeric(),	p18 = numeric(),
                      p19 = numeric(),	p20 = numeric())
    }
    
    write.table(df, file = paste("output/brokenstick/brokenstick_k", k,".txt", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  }
  # First clear final table from older results
  df = read.table(paste("output/brokenstick/brokenstick_k", k,".txt", sep = ""), header = TRUE, sep = "\t")
  df = df[which(df$set != set),]
  # Estimate proportions of the brokenstick
  res = brokenstick(set, k = k, method = "strict")
  resnames = colnames(res)
  res = cbind(rep(set, length(set)), res)
  colnames(res) = c("set", resnames)
  # Problem of factor level
  res$chromosome = as.character(res$chromosome)
  df$chromosome = as.character(df$chromosome)
  # Add results to data & save
  df = rbind(df, res)
  write.table(df, file = paste("output/brokenstick/brokenstick_k", k,".txt", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
}

# Changed the way Brokensticks were computed
# K segments of equal genomic size (Mb) and evaluate the genetic relative size
brokenstick2 = function(set = "", k = 3, method = "strict") {
  map = read.table(paste("data-cleaned/marey_maps/", set, ".txt", sep =""), header = TRUE)
  map$map = as.character(map$map)
  map$mkr = as.character(map$mkr)
  map$phys = as.numeric(as.character(map$phys))
  # List of chromosomes to process
  listchr = unique(map$map)
  # Remove NA chromosomes
  listchr = listchr[!is.na(listchr)]
  # Creat empty dataset
  stick_proportions = as.data.frame(matrix(NA, nrow = length(listchr), ncol = k))
  colnames(stick_proportions) = paste("p", c(1:k), sep ="")
  
  # Proportions for each chromosome
  count = 0
  for (chr in listchr) {
    print(chr)
    # Select the chromosome
    chrmap = map[map$map == chr,]
    # Remove NAs in genetic distances
    chrmap = subset(chrmap, !is.na(chrmap$gen))
    count = count + 1
    
    if (method == "strict") {
      # Segment in k segments of equal genomic size (Mb)
      # Size of segment is total genomic size divided by k
      segment_size = round(max(chrmap$phys, na.rm = TRUE)/k)
      
      # A vector of segments genetic size pi (proportion of total length)
      p = numeric(k)
      for (i in 1:k) {
        # print(max(which(chrmap$gen < segment_size*i)))
        # p[i] = chrmap$phys[max(which(chrmap$gen < segment_size*i))]
        p[i] = max(chrmap$gen[which(chrmap$phys < segment_size*i)], na.rm = TRUE)
        # p[i] = max(chrmap$gen[chrmap$gen < segment_size*i])
        # p[i] = chrmap$phys[max(which(chrmap$gen < quantile(chrmap$gen, i/k)))]
        # p[i] = quantile(chrmap$gen, i/k)
        
        if (i > 1) {
          p[i] = p[i] - sum(p[i-1:(i-1)])
        }
      }
      # END OF METHOD STRICT
    }
    if (method == "segmented") {
      library(segmented)
      res = c()
      lin.mod <- lm(phys~gen, data = chrmap)
      # Estimated breakpoints in genomic distances
      res = try(segmented.mod <- segmented(lin.mod, seg.Z = ~phys, npsi=(k-1), control = seg.control(n.boot = 100, fix.npsi=TRUE)))
      # In cases of segmentation failure
      if (class(res)[1] == "try-error") {
        # Failure because not enough data
        breakpoints = c(NA, NA, NA)
      } else {
        breakpoints = c(summary(segmented.mod)$psi[,2], max(chrmap$phys, na.rm = TRUE))
      }
      
      for (i in 1:k) {
        p[i] = max(chrmap$gen[chrmap$phys < breakpoints[i]])
        
        if (i > 1) {
          p[i] = p[i] - sum(p[i-1:(i-1)])
        }
      }
      # END OF METHOD SEGMENTED
    }
    
    # p is a proportion of total physical length
    p = p/sum(p)
    
    stick_proportions[count,] = rbind(p)
    # END OF CHROMOSOME LOOP
  }
  stick_proportions = cbind(listchr, stick_proportions)
  colnames(stick_proportions)[1] = "chromosome"
  stick_proportions
  
  return(stick_proportions)
}

batch_brokenstick2 = function(set = "", k = 10) {
  print(set)
  if (!file.exists(paste("output/brokenstick/brokenstick_k", k,"_v2.txt", sep = ""))) {
    if (k == 10) {
      df = data.frame(set = character(),	chromosome = character(),	p1 = numeric(),	p2 = numeric(),	p3 = numeric(),
                      p4 = numeric(),	p5 = numeric(),	p6 = numeric(),	p7 = numeric(),	p8 = numeric(),	p9 = numeric(),	p10 = numeric())
    }
    if (k == 20) {
      df = data.frame(set = character(),	chromosome = character(),	p1 = numeric(),	p2 = numeric(),	p3 = numeric(),
                      p4 = numeric(),	p5 = numeric(),	p6 = numeric(),	p7 = numeric(),	p8 = numeric(),	p9 = numeric(),	p10 = numeric(),
                      p11 = numeric(),	p12 = numeric(),	p13 = numeric(),
                      p14 = numeric(),	p15 = numeric(),	p16 = numeric(),	p17 = numeric(),	p18 = numeric(),
                      p19 = numeric(),	p20 = numeric())
    }

    write.table(df, file = paste("output/brokenstick/brokenstick_k", k,"_v2.txt", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  }
  # First clear final table from older results
  df = read.table(paste("output/brokenstick/brokenstick_k", k,"_v2.txt", sep = ""), header = TRUE, sep = "\t")
  df = df[which(df$set != set),]
  # Estimate proportions of the brokenstick
  res = brokenstick2(set, k = k, method = "strict")
  resnames = colnames(res)
  res = cbind(rep(set, length(set)), res)
  colnames(res) = c("set", resnames)
  # Problem of factor level
  res$chromosome = as.character(res$chromosome)
  df$chromosome = as.character(df$chromosome)
  # Add results to data & save
  df = rbind(df, res)
  write.table(df, file = paste("output/brokenstick/brokenstick_k", k,"_v2.txt", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
}



#============================================================================#
# UPDATE CLEANED MAREY MAPS AND ASSOCIATED METADATA
#============================================================================#

update_datacleaned = function() {
  # Return a data frame with the name of the map and the flag 'valid' associated, after filtering
  metadata = read.table(file = paste(wd, "/data/Genetic_maps/Genetic_maps_ressources.csv", sep = ""),
                      header = TRUE, sep = ";")
  vld.maps = c() # init the vector of maps selected and validated
  
  # List of all maps
  list = system(paste("ls ", wd, "/data/Marey_maps/", sep = ""), intern = TRUE)
  list = list[grep("*.txt", list)]
  length(list)
  list
  # Gossypium_hirsutum_Jia2016 removed, sept 22th 2020, too much overfitting
  # Replicates removed: "Oryza_sativa_Jiang2017", "Vitis_vinifera_Zhu2018", "Zea_mays_IBM_MaizeSNP50",
  # "Vitis_vinifera_Tello2019", "Citrullus_lanatus_Branham2017", "Camelina_sativa_Singh2015", "Glycine_max_SoyBase"
  # "Prunus_avium_Shirasawa2017", "Setaria_italica_Ni2017",
  vld.maps = sort(unique(c(vld.maps, "Oryza_sativa_DeLeon2016", "Sesamum_indicum_Wang2016",
                           "Triticum_aestivum_GutierrezGonzalez2019", "Arabidopsis_thaliana_Serin2017",
                           "Solanum_tuberosum_Endelman2016", "Solanum_lycopersicum_Gonda2018", "Brachypodium_distachyon_Huo2011",
                           "Citrullus_lanatus_Ren2015", "Gossypium_raimondii_Wang2013",
                           "Malus_domestica_DiPierro2016", "Phaseolus_vulgaris_Song2015", "Prunus_mume_Zhang2015",
                           "Theobroma_cacao_Royaert2016", "Zea_mays_MaizeGDBConsensus_v4",
                           "Cucumis_sativus_Zhu2016", "Manihot_esculenta_ICGMC2015", "Sorghum_bicolor_Zou2012",
                           "Cucumis_melo_Pereira2018",
                           "Hordeum_vulgare_MunozAmatriain2014", "Lupinus_angustifolius_Zhou2017",
                           "Cenchrus_americanus_Pucher2017",
                           "Brassica_napus_Yang2017", "Prunus_persica_Verde2017",
                           "Capsicum_annuum_Han2016", "Aegilops_speltoides_Zhang2019", "Vitis_vinifera_Brault2020", "Coffea_canephora_Crouzillat2020",
                           "Triticum_dicoccoides_Jorgensen2017", "Triticum_urartu_Ling2018", "Arachis_duranensis_Bertioli2016", "Panicum_hallii_Lovell2018",
                           "Arachis_hypogaea_Zhuang2019", "Boechera_stricta_Lee2020", "Momordica_charantia_Mastumura2020",
                           "Vigna_unguiculata_Lonardi2019", "Elaeis_guineensis_Yaakub2020", "Brassica_rapa_CodyMarkelz2017",
                           "Cucurbita_pepo_MonteroPau2017", "Dioscorea_alata_Cormier2019", "Draba_nivalis_Birkeland2020",
                           "Quercus_sp_Bodenes2016", "Gossypium_hirsutum_Zhang2019", "Setaria_italica_Bennetzen2012",
                           "Raphanus_sativus_Mun2015", "Lupinus_albus_Ksiazkiewicz2017", "Mangifera_indica_Luo2016",
                           "Nelumbo_nucifera_Gui2018", "Oryza_nivara_Ma2016",
                           "Aegilops_speltoides_Zhang2020", "Camelina_sativa_King2019", "Capsella_rubella_Slotte2013",
                           "Helianthus_annuus_Talukder2014", "Glycine_max_Watanabe2017", 
                           "Citrus_sinensis_Huang2018", "Eucalyptus_grandis_Bertholome2015", "Juglans_regia_Luo2015",
                           "Cucurbita_maxima_Wang2020", "Camellia_sinensis_Xu2018")))
  
  # Deprecated maps
  # "Phoenix_dactylifera_Mathew2014"
  # 
  # Map filtering per chromosome
  # by Number of markers: reject (vld = FALSE) chromosomes with less than the required number of markers
  
  # A subset of maps both automatically & manually validated with the flag 'valid = TRUE'
  metadata$valid = FALSE # Reset the flag before updating
  for (i in 1:length(vld.maps)) {
    metadata$valid[which(as.character(metadata$id) == vld.maps[i])] = TRUE
  }
  
  # Saving the updated database: TRUE/FALSE when map used/not used
  write.table(metadata, file = paste(wd, "/data/Genetic_maps/Genetic_maps_ressources.csv", sep = ""),
              col.names = TRUE, row.names = FALSE, sep = ";")
  # rm(metadata)
  #------------------------#
  # Clean directory first
  system("rm data-cleaned/marey_maps/*")
  # Maps with the flag 'valid = TRUE' are then copied into data-cleaned for records
  # but only vld markers are copied
  for (i in vld.maps) {
    print(i)
    df = read.table(paste("data/Marey_maps/", i, ".txt", sep = ""), header = TRUE, sep = "\t")
    df = df[df$vld == TRUE,]
    write.table(df, file = paste("data-cleaned/marey_maps/", i, ".txt", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  }
  # And then all maps are concatened in a single final dataset
  # starting by dataset number one
  df = read.table(paste("data/Marey_maps/", vld.maps[1], ".txt", sep = ""), header = TRUE, sep = "\t")
  df$mkr = as.character(df$mkr)
  df$map = as.character(df$map)
  for (i in vld.maps[-1]) {
    print(i)
    tmp = read.table(paste("data/Marey_maps/", i, ".txt", sep = ""), header = TRUE, sep = "\t")
    tmp$map = as.character(tmp$map)
    df = rbind(df, tmp)
  }
  # Remove FALSE markers
  df = df[df$vld == TRUE,]
  # Remove NA distances & chromosomes
  df = df[!is.na(df$phys),]
  df = df[!is.na(df$gen),]
  df = df[!is.na(df$map),]
  
  write.table(df, file = "data-cleaned/marey_maps/AllMaps.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  allMaps = read.table(file = "data-cleaned/marey_maps/AllMaps.txt", header = TRUE, sep = "\t")
  
  
  # Batch copy in shell
  #file.copy(from = paste("data/Marey_maps/", vld.maps, ".txt", sep = ""), to = paste("data-cleaned/Marey_maps/", vld.maps, ".txt", sep = ""), overwrite = TRUE)
  #------------------------#
  # Diagnostic plots
  # Genetic distance (cM) as a function of physical position (bp)
  # Make a diagnostic plot of all maps in 'data-cleaned/Marey_maps/'
  # List of all maps
  list = system(paste("ls ", wd, "/data-cleaned/marey_maps/", sep = ""), intern = TRUE)
  list = list[grep("*.txt", list)]
  length(list)
  list = list[!(list == "AllMaps.txt")]
  list
  # Save diagnostic plots as one figure per dataset in 'output/marey_maps/diagnostic_plots/'
  # Valid markers are in blue and unvalid markers are in red
  # list = list[14]
  # Clean directory first
  system("rm output/marey_maps/diagnostic_plots-cleaned/*")
  marey.diagplot(set = list, dir = "data-cleaned/marey_maps/", output = "output/marey_maps/diagnostic_plots-cleaned/", save.plot = TRUE, display = TRUE)
  rm(list)
  
  #------------------------#
  # Saving the updated & cleaned database (only clean marey maps)
  metadata.clean = read.table(paste(wd, "/data/Genetic_maps/Genetic_maps_ressources.csv", sep = ""),
                            header = TRUE, sep = ";")
  metadata.clean = metadata.clean[which(metadata$valid == TRUE),]
  write.table(metadata.clean, file = paste(wd, "/data-cleaned/Genetic_maps_ressources.csv", sep = ""),
              col.names = TRUE, row.names = FALSE, sep = ";")
  
  # Display the list of maps already cleaned
  list_cleaned = system(paste("ls ", wd, "/data-cleaned/marey_maps/", sep = ""), intern = TRUE)
  (list_cleaned = list_cleaned[grep("*.txt", list_cleaned)])
  # How many maps?
  length(list_cleaned)
}


#============================================================================#
# STATISTICS OF SELECTED MAPS
#============================================================================#
# Marey map selection step
# Summary statistics to trim the dataset
# A function returning a data frame with summary statistics of selected maps
map.statistics = function(data.final) {
  # Subset valid markers
  data.final = subset(data.final, vld == TRUE)
  data.final$map = as.character(data.final$map)
  data.final$mkr = as.character(data.final$mkr)
  data.final$phys = as.numeric(as.character(data.final$phys))
  data.final$gen = as.numeric(data.final$gen)

  # Preallocate memory for the dataframe
  list.maps = unique(paste(data.final$set, data.final$map, sep = "-")) # List of all maps, i.e. chromosomes
  nb.maps = length(list.maps) # The total number of maps, i.e. chromosomes
  chromosome.stats = data.frame(set = as.character(gsub("-[0-9A-Za-z]+", "", list.maps)), chromosome = as.character(str_extract(list.maps, "[0-9A-Za-z]+$")), valid = rep(TRUE, nb.maps), linkage.map.length = numeric(nb.maps),
                                linkage.map.length.correctedCh = numeric(nb.maps), linkage.map.length.correctedHW = numeric(nb.maps), phys.map.length = numeric(nb.maps), nb.markers = numeric(nb.maps),
                                density.markers.cM = numeric(nb.maps), density.markers.bp = numeric(nb.maps), marker.interval.cM = numeric(nb.maps), largest.gap.cM = numeric(nb.maps), marker.interval.bp = numeric(nb.maps),
                                chrwide.rate = numeric(nb.maps), map.coverage = numeric(nb.maps), rejected.rate = numeric(nb.maps)
  )
  # Then, estimate summary statistics for each map
  data.final$set = as.character(data.final$set)
  chromosome.stats$set = as.character(chromosome.stats$set)
  chromosome.stats$chromosome = as.character(chromosome.stats$chromosome)
  for (i in 1:nrow(chromosome.stats)) {
    chromosome.stats$linkage.map.length[i] = max(data.final$gen[which(as.character(data.final$set) == as.character(chromosome.stats$set[i]) & as.character(data.final$map) == as.character(chromosome.stats$chromosome[i]))], na.rm = TRUE)
    # (1) correction of Chakravarti et al. (1991): multiplies the linkage group length by (m + 1)/(m - 1), where m is the number of framework markers on each group
    m = length(data.final$phys[which(data.final$set == chromosome.stats$set[i] & data.final$map == chromosome.stats$chromosome[i])])
    chromosome.stats$linkage.map.length.correctedCh[i] = max(data.final$gen[which(data.final$set == chromosome.stats$set[i] & data.final$map == chromosome.stats$chromosome[i])], na.rm = TRUE) * ((m + 1)/(m - 1))
    # (2) correction 2 of Hall & Willis (2005): add 2s to the length of each linkage group, s being the average marker spacing in the group
    tmp = data.final$gen[which(data.final$set == chromosome.stats$set[i] & data.final$map == chromosome.stats$chromosome[i])]
    chromosome.stats$linkage.map.length.correctedHW[i] = max(data.final$gen[which(data.final$set == chromosome.stats$set[i] & data.final$map == chromosome.stats$chromosome[i])], na.rm = TRUE) + 2*abs(mean(tmp[-1]-head(tmp, -1), na.rm = TRUE))
    
    chromosome.stats$phys.map.length[i] = max(as.numeric(data.final$phys[which(data.final$set == chromosome.stats$set[i] & data.final$map == chromosome.stats$chromosome[i])]), na.rm = TRUE)
    
    chromosome.stats$nb.markers[i] = length(data.final$phys[which(data.final$set == chromosome.stats$set[i] & data.final$map == chromosome.stats$chromosome[i])])
    chromosome.stats$density.markers.cM[i] = chromosome.stats$nb.markers[i]/chromosome.stats$linkage.map.length[i]
    chromosome.stats$density.markers.bp[i] = chromosome.stats$nb.markers[i]/(chromosome.stats$phys.map.length[i])
    tmp = data.final$gen[which(data.final$set == chromosome.stats$set[i] & data.final$map == chromosome.stats$chromosome[i])]
    # Ordering markers
    tmp = sort(tmp)
    chromosome.stats$marker.interval.cM[i] = mean(abs(tmp[-1]-head(tmp, -1)), na.rm = TRUE)
    chromosome.stats$largest.gap.cM[i] = max(abs(tmp[-1]-head(tmp,-1)), na.rm = TRUE)
    tmp = data.final$phys[which(data.final$set == chromosome.stats$set[i] & data.final$map == chromosome.stats$chromosome[i])]
    tmp = sort(tmp)
    chromosome.stats$marker.interval.bp[i] = mean(abs(tmp[-1]-head(tmp, -1)), na.rm = TRUE)
    chromosome.stats$chrwide.rate[i] = chromosome.stats$linkage.map.length.correctedHW[i]/(chromosome.stats$phys.map.length[i]/1000000)
    chromosome.stats$map.coverage[i] = 1 - exp(-2 * chromosome.stats$nb.markers[i]/chromosome.stats$linkage.map.length[i])
    rm(tmp)
  }
  
  # Estimating the rejected rate
  for (i in 1:nrow(chromosome.stats)) {
    # The rejected rate is the number of markers with flag = FALSE divided by the total number of markers
    # Need to assess it in /data/, since cleaned maps do not have FALSE markers
    # Open original map in /data/
    if (file.exists(paste("data/Marey_maps/", chromosome.stats$set[i], ".txt", sep = ""))) {
      tmp.map = read.table(file = paste("data/Marey_maps/", chromosome.stats$set[i], ".txt", sep = ""), header = TRUE, sep ="\t",
                           stringsAsFactors = FALSE)
      tmp.map$map = as.character(tmp.map$map)
      tmp.map$phys = as.numeric(as.character(tmp.map$phys))
      tmp.map$gen = as.numeric(as.character(tmp.map$gen))
      # Subset only non-NA values
      tmp.map = subset(tmp.map, !is.na(map))
      # Select chromosome
      tmp.map = tmp.map[tmp.map$map == as.character(chromosome.stats$chromosome[i]),]
      # Count number of vld == FALSE and divide by total number of markers
      chromosome.stats$rejected.rate[i] = sum(tmp.map$vld == FALSE)/nrow(tmp.map)
    } else {
      chromosome.stats$rejected.rate[i] = NA
    }
  }
  
  return(chromosome.stats)
}


#============================================================================##
# STATISTICS OF ESTIMATED RECOMBINATION MAPS ----
#============================================================================##
# Marey map selection step
# A function returning a data frame with summary statistics of selected maps
# AND saving a file with summary statistics
recombination.map.statistics = function(data = data.final, save.file = TRUE) {
  require(readxl)
  # PER CHROMOSOME SUMMARY STATISTICS #
  #-----------------------------------------##
  # Describe first the Marey maps
  # Linkage map length (cM), the y axis
  # The corrected map length, using the methods of Chakravarti et al. (1991) for variation in marker density and Hall and Willis (2005), for undected events of recombination in distal regions
  #       (1) correction of Chakravarti et al. (1991): multiplies the linkage group length by (m + 1)/(m - 1), where m is the number of framework markers on each group
  #       (2) correction 2 of Hall & Willis (2005): add 2s to the length of each linkage group, s being the average marker spacing in the group
  # Physical map length (bp), the x axis
  # Number of markers
  # The density in markers (number of markers/cM & number of markers/bp)
  # Average marker interval (cM)
  # Largest gap (cM)
  # Average marker interval (bp)
  # C-values (pg)
  # Map coverage (Chakravarti et al. (1991), Hall and Willis (2005))
  #         "map coverage c. proportion c of the genome that is within distance d cM of a marker, assuming random distribution of markers,
  #         was estimated using c = 1 - exp(-2dn/L), where n is the number of markers and L is the estimated genome length"
  #   We estimated the map coverage for a distance d of 1 cM
  # Chromosome-wide recombination rate (cM/Mb) estimated by dividing the genetic total length by the physical total length
  
  # Parameters of the Marey map: fitdistR
  # Recombination kernel
  
  # Satistics of the estimated recombination landscape
  # Mean of the local recombination rate (cM/Mb, estimated)
  # Median of the local recombination rate (cM/Mb, estimated)
  # SD of the local recombination rate (cM/Mb, estimated)
  
  # Crossover interference
  # Gamma Srinkling model (Copenhaver et al. 2002, Bauer et al. 2013) -> parameters nu and p
  
  # rbar, a new metrics for measuring genome-wide shuffling developped by Veller et al. (2019) taking into account the genetic position of the CO in the genetic shuffling rate
  # "This measure, rbar, is the probability that a randomly chosen pair of loci shuffles their alleles in meiosis."
  # rbar can be decomposed in two components of intra-chr and inter-chr shuffling
  
  # And statistics of the genome architecture
  # Centromere physical position (bp)
  # Type of chromosome (centromeric classification)
  # Mean GC content (%)
  # Nonsynonymous heterozygosity in coding sequences
  # dN/dS
  
  # + Metadata
  # Species
  # Sex
  # Wild/Domesticated
  # Hybrid
  
  
  # Statistics of map quality
  # Map coverage
  # Rejected markers rate
  
  cat("Adding linkage map statistics...\n")
  chromosome.stats = map.statistics(data)
  
  # # Preallocate memory for the dataframe
  # list.maps = unique(paste(data.final$set, data.final$map, sep = "-")) # List of all maps, i.e. chromosomes
  # nb.maps = length(list.maps) # The total number of maps, i.e. chromosomes
  # chromosome.stats = data.frame(set = as.character(gsub("-[0-9A-Za-z]+", "", list.maps)), chromosome = as.character(str_extract(list.maps, "[0-9A-Za-z]+$")), valid = rep(TRUE, nb.maps), linkage.map.length = numeric(nb.maps),
  #                               linkage.map.length.correctedCh = numeric(nb.maps), linkage.map.length.correctedHW = numeric(nb.maps), phys.map.length = numeric(nb.maps), nb.markers = numeric(nb.maps),
  #                               density.markers.cM = numeric(nb.maps), density.markers.bp = numeric(nb.maps), marker.interval.cM = numeric(nb.maps), largest.gap.cM = numeric(nb.maps), marker.interval.bp = numeric(nb.maps),
  #                               chrwide.rate = numeric(nb.maps), map.coverage = numeric(nb.maps), rejected.rate = numeric(nb.maps)
  # )
  
  
  # Save the template file and then re-import for each operations...
  # write.table(chromosome.stats, file = "tables/chromosome.stats.csv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  # chromosome.stats = read.table("tables/chromosome.stats.csv", header = TRUE)
  
  # # Then, estimate summary statistics for each map
  # data.final$set = as.character(data.final$set)
  # chromosome.stats$set = as.character(chromosome.stats$set)
  # data.final$map = as.character(data.final$map)
  # data.final$mkr = as.character(data.final$mkr)
  # data.final$phys = as.numeric(as.character(data.final$phys))
  # chromosome.stats$chromosome = as.character(chromosome.stats$chromosome)
  # for (i in 1:nrow(chromosome.stats)) {
  #   chromosome.stats$linkage.map.length[i] = max(data.final$gen[which(as.character(data.final$set) == as.character(chromosome.stats$set[i]) & as.character(data.final$map) == as.character(chromosome.stats$chromosome[i]))], na.rm = TRUE)
  #   # (1) correction of Chakravarti et al. (1991): multiplies the linkage group length by (m + 1)/(m - 1), where m is the number of framework markers on each group
  #   m = length(data.final$phys[which(data.final$set == chromosome.stats$set[i] & data.final$map == chromosome.stats$chromosome[i])])
  #   chromosome.stats$linkage.map.length.correctedCh[i] = max(data.final$gen[which(data.final$set == chromosome.stats$set[i] & data.final$map == chromosome.stats$chromosome[i])], na.rm = TRUE) * ((m + 1)/(m - 1))
  #   # (2) correction 2 of Hall & Willis (2005): add 2s to the length of each linkage group, s being the average marker spacing in the group
  #   tmp = data.final$gen[which(data.final$set == chromosome.stats$set[i] & data.final$map == chromosome.stats$chromosome[i])]
  #   chromosome.stats$linkage.map.length.correctedHW[i] = max(data.final$gen[which(data.final$set == chromosome.stats$set[i] & data.final$map == chromosome.stats$chromosome[i])], na.rm = TRUE) + 2*abs(mean(tmp[-1]-head(tmp, -1), na.rm = TRUE))
  #   chromosome.stats$phys.map.length[i] = max(data.final$phys[which(data.final$set == chromosome.stats$set[i] & data.final$map == chromosome.stats$chromosome[i])], na.rm = TRUE)
  #   chromosome.stats$nb.markers[i] = length(data.final$phys[which(data.final$set == chromosome.stats$set[i] & data.final$map == chromosome.stats$chromosome[i])])
  #   chromosome.stats$density.markers.cM[i] = chromosome.stats$nb.markers[i]/chromosome.stats$linkage.map.length[i]
  #   chromosome.stats$density.markers.bp[i] = chromosome.stats$nb.markers[i]/(chromosome.stats$phys.map.length[i])
  #   tmp = data.final$gen[which(data.final$set == chromosome.stats$set[i] & data.final$map == chromosome.stats$chromosome[i])]
  #   tmp = sort(tmp)
  #   chromosome.stats$marker.interval.cM[i] = mean(abs(tmp[-1]-head(tmp, -1)), na.rm = TRUE)
  #   chromosome.stats$largest.gap.cM[i] = max(abs(tmp[-1]-head(tmp,-1)), na.rm = TRUE)
  #   tmp = data.final$phys[which(data.final$set == chromosome.stats$set[i] & data.final$map == chromosome.stats$chromosome[i])]
  #   tmp = sort(tmp)
  #   chromosome.stats$marker.interval.bp[i] = mean(abs(tmp[-1]-head(tmp, -1)), na.rm = TRUE)
  #   chromosome.stats$chrwide.rate[i] = chromosome.stats$linkage.map.length[i]/(chromosome.stats$phys.map.length[i]/1000000)
  #   chromosome.stats$map.coverage[i] = 1 - exp(-2 * chromosome.stats$nb.markers[i]/chromosome.stats$linkage.map.length[i])
  # }
  
  # # Estimating the rejected rate
  # for (i in 1:nrow(chromosome.stats)) {
  #   # The rejected rate is the number of markers with flag = FALSE divided by the total number of markers
  #   # Need to assess it in /data/, since cleaned maps do not have FALSE markers
  #   # Open original map in /data/
  #   tmp.map = read.table(file = paste("data/Marey_maps/", chromosome.stats$set[i], ".txt", sep = ""), header = TRUE, sep ="\t")
  #   # Select chromosome
  #   tmp.map = tmp.map[tmp.map$map == as.character(chromosome.stats$chromosome[i]),]
  #   # Count number of vld == FALSE and divide by total number of markers
  #   chromosome.stats$rejected.rate[i] = sum(tmp.map$vld == FALSE)/nrow(tmp.map)
  # }
  
  #-----------------------------------------#
  # How to estimate/quantify the quality of a map? Quality score?
  
  # Density: number of markers per chromosome size OR linkage length -> retain only high density maps
  # Compute density per chromosome
  # Average marker interval (cM)
  
  #-----------------------------------------#
  # Satistics of the estimated recombination landscape
  # Mean of the local recombination rate (cM/Mb, estimated)
  # Median of the local recombination rate (cM/Mb, estimated)
  # SD of the local recombination rate (cM/Mb, estimated)
  cat("Adding recombination landscapes statistics...\n")
  # Initialise new columns with col names
  chromosome.stats$mean.recrate = NA
  chromosome.stats$median.recrate = NA
  chromosome.stats$sd.recrate = NA
  chromosome.stats$max = NA
  chromosome.stats$cv.recrate = NA
  
  chromosome.stats$set = as.character(chromosome.stats$set)
  chromosome.stats$chromosome = as.character(chromosome.stats$chromosome)
  
  # Which recombination maps are considered as reference for statistics
  ref.path = "/output/recombination_maps/loess/100kbwind/"
  # Compute statistics per chromosome
  for (i in 1:nrow(chromosome.stats)) {
    print(i)
    # print(paste(chromosome.stats$set[i], "_chromosome", chromosome.stats$chromosome[i], sep = ""))
    # Import reference map
    if (file.exists(file = paste(wd, ref.path, chromosome.stats$set[i], "_chromosome", chromosome.stats$chromosome[i], ".txt", sep = ""))) {
      ref.map = read.table(file = paste(wd, ref.path, chromosome.stats$set[i], "_chromosome", chromosome.stats$chromosome[i], ".txt", sep = ""), sep = "\t", header = TRUE)
      # Compute mean
      chromosome.stats$mean.recrate[i] = mean(ref.map$rec.rate, na.rm = TRUE)
      # Compute median
      chromosome.stats$median.recrate[i] = median(ref.map$rec.rate, na.rm = TRUE)
      # Compute SD
      chromosome.stats$sd.recrate[i] = sd(ref.map$rec.rate, na.rm = TRUE)
      # Maximum recombination rate
      chromosome.stats$max[i] = max(ref.map$rec.rate, na.rm = TRUE)
      # Coefficient of variation
      chromosome.stats$cv.recrate[i] = chromosome.stats$sd.recrate[i]/chromosome.stats$mean.recrate[i]
    } else {
      cat(chromosome.stats$set[i], "_chromosome", chromosome.stats$chromosome[i], " don't exist.\n", sep = "")
    }
  }
  
  # The Gini Index
  chromosome.stats$gini = NA
  for (i in 1:nrow(chromosome.stats)) {
    if (file.exists(file = paste(wd, ref.path, chromosome.stats$set[i], "_chromosome", chromosome.stats$chromosome[i], ".txt", sep = ""))) {
      # Compute the Gini index on the map
      chromosome.stats$gini[i] = gini.coefficient(set = chromosome.stats$set[i], map = chromosome.stats$chromosome[i])
    } else {
      cat(chromosome.stats$set[i], "_chromosome", chromosome.stats$chromosome[i], " don't exist.\n", sep = "")
    }
  }
  
  # The span parameter estimated by cross-validation is to be saved in 'chromosome.stats'
  cat("Adding calibrated span parameter...\n")
  span = read.table("output/recombination_maps/span.txt", header = TRUE)
  chromosome.stats$span = NA
  for (i in 1:nrow(chromosome.stats)) {
    if (length(span$span[which(span$set == chromosome.stats$set[i] & span$map == chromosome.stats$chromosome[i])]) > 0) {
      chromosome.stats$span[i] = span$span[which(span$set == chromosome.stats$set[i] & span$map == chromosome.stats$chromosome[i])]
    }
  }
  
  # Add the Centromeric Index to 'chromosome.stats'
  cat("Adding centromerix indexes...\n")
  # require(readxl)
  karyotypes = read_xlsx("data-cleaned/Karyotypes.xlsx", sheet = 1)
  # karyotypes$sp = paste(karyotypes$Genus, karyotypes$Species, sep = "_")
  # Use the centromeric index as a measure of the centromere position
  karyotypes$Centromeric_Index = as.numeric(as.character(karyotypes$Centromeric_Index))
  karyotypes$species = paste(karyotypes$Genus, karyotypes$Species, sep = " ")
  karyotypes$Chromosome = as.character(karyotypes$Chromosome)
  chromosome.stats$chromosome = as.character(chromosome.stats$chromosome)
  chromosome.stats$species = gsub("_", " ", regmatches(chromosome.stats$set, regexpr("^[A-Za-z]+_[A-Za-z]+", chromosome.stats$set, perl=TRUE)))
  chromosome.stats$centromeric_index = NA
  for (i in 1:nrow(chromosome.stats)) {
    idx = which(karyotypes$species == chromosome.stats$species[i] & karyotypes$Chromosome == chromosome.stats$chromosome[i])
    if (length(idx) > 0) {
      chromosome.stats$centromeric_index[i] = round(karyotypes$Centromeric_Index[idx][1], digits = 2)
    } else {
      chromosome.stats$centromeric_index[i] = NA
    }
  }
  
  # summary(chromosome.stats$centromeric_index)
  # sum(!is.na(chromosome.stats$centromeric_index))
  # sum(is.na(chromosome.stats$centromeric_index))

  #============================================================================#
  # Periphery-bias ratio
  #============================================================================#
  cat("Adding periphery-bias statistics...\n")
  chromosome.stats$peripherybias_ratio = NA
  chromosome.stats$peripherybias_mode = NA
  df = read.table("output/dist2telomere/DistancesRelative_bins.txt", header = TRUE)
  for (i in 1:nrow(chromosome.stats)) {
    subs = subset(df, set == chromosome.stats$set[i] & chromosome == chromosome.stats$chromosome[i])
    # The ratio of end/rest of the chromosome
    N = 2
    chromosome.stats$peripherybias_ratio[i] = (mean(subs$rec.rate[c(1:N)], na.rm = TRUE)/mean(subs$rec.rate, na.rm = TRUE))
    chromosome.stats$peripherybias_mode[i] = which(subs$rec.rate == max(subs$rec.rate, na.rm = TRUE)) # The bin with the maximal value of recombination
  }
  # Estimate the median point, where half the recombination already occured
  # chromosome.stats$peripherybias_median = NA
  # Cut the half chromosome in two parts of equal genetic length
  
  #============================================================================#
  # Variance in segment proportions - An additional measure of evenness
  #============================================================================#
  cat("Adding variance in brokenstick proportions...\n")
  load(file = "output/brokenstick/brokenstick10.Rda")
  chromosome.stats$brokenstick_pvariance = NA
  for (i in 1:nrow(chromosome.stats)) {
    sample = paste(chromosome.stats$set[i], chromosome.stats$chromosome[i], sep = "_")
    chromosome.stats$brokenstick_pvariance[i] = var(brokenstick$proportion.length[which(brokenstick$sample == sample)], na.rm = TRUE)
  }
  

  #============================================================================#
  # Variance in inter-marker recombination rates
  #============================================================================#
  # cat("Adding variance in inter-marker recombination rates...\n")
  # library(pbmcapply)
  # varbetweenmarkers = unlist(pbmclapply(1:nrow(chromosome.stats), function(x) variance_intermarker(chromosome.stats$set[x],
  #                                                                                                  chromosome.stats$chromosome[x])))
  # cvbetweenmarkers = unlist(pbmclapply(1:nrow(chromosome.stats), function(x) cvariation_intermarker(chromosome.stats$set[x],
  #                                                                                                   chromosome.stats$chromosome[x])))
  # chromosome.stats$varbetweenmarkers = varbetweenmarkers
  # chromosome.stats$cvbetweenmarkers = cvbetweenmarkers

  #============================================================================#
  # Number of genes in the chromosome
  #============================================================================#
  chromosome.stats$genecount = NA
  pb = txtProgressBar(min = 0, max = nrow(chromosome.stats), style = 3)
  for (i in 1:nrow(chromosome.stats)) {
    if (file.exists(paste("data-cleaned/genome/gene_positions/", gsub(" ", "_", chromosome.stats$species[i]),
                          "_gene_positions.txt.gz" , sep = ""))) {
      genepos = read.table(gzfile(paste("data-cleaned/genome/gene_positions/", gsub(" ", "_", chromosome.stats$species[i]),
                                        "_gene_positions.txt.gz" , sep = "")), header = TRUE)
      
      # Deal with exceptions in chromosome naming
      if (chromosome.stats$species[i] == "Brassica rapa" | chromosome.stats$species[i] == "Brassica napus" |
          chromosome.stats$species[i] == "Arachis hypogaea") {
        # Rename chromosomes A1 to A01
        genepos$chromosome = gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", genepos$chromosome, perl = TRUE)
      }
      
      if (chromosome.stats$species[i] == "Theobroma cacao") {
        genepos$chromosome = gsub("CM001879.1", "1", genepos$chromosome)
        genepos$chromosome = gsub("CM001880.1", "2", genepos$chromosome)
        genepos$chromosome = gsub("CM001881.1", "3", genepos$chromosome)
        genepos$chromosome = gsub("CM001882.1", "4", genepos$chromosome)
        genepos$chromosome = gsub("CM001883.1", "5", genepos$chromosome)
        genepos$chromosome = gsub("CM001884.1", "6", genepos$chromosome)
        genepos$chromosome = gsub("CM001885.1", "7", genepos$chromosome)
        genepos$chromosome = gsub("CM001886.1", "8", genepos$chromosome)
        genepos$chromosome = gsub("CM001887.1", "9", genepos$chromosome)
        genepos$chromosome = gsub("CM001888.1", "10", genepos$chromosome)
      }
      
      # Count the number of unique IDs
      # chromosome.stats$genenumber[i] = length(unique(genepos$id[genepos$type == "gene" & genepos$chromosome == chromosome.stats$chromosome[i]]))
      # Problem, sometimes there is no ID
      chromosome.stats$genecount[i] = sum(genepos$type == "gene" & genepos$chromosome == chromosome.stats$chromosome[i])
      setTxtProgressBar(pb, i)
    }
  }
  close(pb)
  
  chromosome.stats$genecount[which(chromosome.stats$genecount == 0)] = NA
  #============================================================================#
  # Centromere position
  cat("Estimating the centromere position and orienting the centromeric index...\n")
  chromosome.stats$centromere_position_estimated = NA
  pb = txtProgressBar(min = 0, max = nrow(chromosome.stats), style = 3)
  for (i in 1:nrow(chromosome.stats)) {
    chromosome.stats$centromere_position_estimated[i] = estim_centromere_pos(chromosome.stats$set[i], chromosome.stats$chromosome[i])
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Orient the Centromeric Index based on the putative centromere position
  scale = 1000000
  chromosome.stats$centromeric_index_position_oriented = NA
  for (i in 1:nrow(chromosome.stats)) {
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
  
  #============================================================================#
  # Save in tables
  if (save.file == TRUE) {
    cat("File 'chromosome.stats.csv' saved.\n")
    write.table(chromosome.stats, file = "tables/chromosome.stats.csv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ";")
    chromosome.stats = read.table("tables/chromosome.stats.csv", header = TRUE, sep =";")
  }
  
  
  return(chromosome.stats)
}


#============================================================================#
# GINI COEFFICIENT
#============================================================================#
# Estimate the Gini coefficient of a given Marey map
# Where map is a chromosome name
# the Gini coefficient expresses the unevenness/inequality of recombination
# Calculated on the distribution of recombination rates
# The corrected Gini coefficient is between 0 and 1. Values close to 0 are the most equal.
gini.coefficient = function(set = "", map = "") {
  require(reldist)
  # Import reference map
  ref.map = read.table(file = paste(wd, "/output/recombination_maps/loess/100kbwind/", set, "_chromosome",map, ".txt", sep = ""), sep = "\t", header = TRUE)
  # Compute the Gini index on the map
  # Twice the area under the Lorenz curve (i.e. proportion of recombination ~ proportion of sequence)
  # see xueReducedMeioticRecombination2020
  # G = 2A = 1-2B
  # Remove NA values
  gini = gini(ref.map$rec.rate[!is.na(ref.map$rec.rate)])
  
  return(gini)
}

#============================================================================#
# INTER-MARKER VARIANCE
#============================================================================#
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



#============================================================================#
# INTER-MARKER COEFFICIENT OF VARIATION
#============================================================================#
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


#============================================================================#
#   Orient the Centromeric Index on Marey maps
#============================================================================#
# Estimate the position of the putative centromere on the Marey map
# Use the interpolated loess function with span already inferred in doing the recombination map
# and identify the plateau region, where the monotonic function stops increasing or has the lowest slope
# i.e. the region with the lowest recombination
# Marey maps were cut in 100 bins (N=100) of equal physical size, to calculate the derivate of each bin
# Large scale smoothing, since we want to detect the most important region of low recombination
# rather than local peaks 

# estim_centromere_pos(set, chromosome)
# Function to estimate a putative centromere position on a Marey map
# Return the putative centromere position in physical position (Mb, same orientation as the reference genome)
estim_centromere_pos = function(set, chromosome, N = 100) {
  source("sources/MareyMap.R")
  # Physical distances are in bp, yet analyses are in Mb
  # Convert bp to Mb
  scale = 1000000
  
  # Load the Marey map
  data = read.table(file = paste("data-cleaned/marey_maps/", set,".txt", sep = ""), header = TRUE, sep = "\t")
  data$set = as.character(data$set)
  data$map = as.character(data$map)
  chr.data = data[which(data$map == as.character(chromosome) & data$set == set),] # Sample only distances of the chromosome
  # chr.data = subset(data, (map == as.character(chromosome) & set == set)) # Sample only distances of the chromosome
  chr.data$gen = as.numeric(as.character(chr.data$gen))
  chr.data$phys = as.numeric(as.character(chr.data$phys))
  # Convert physical distances in Mb
  chr.data$phys = chr.data$phys/scale
  # Remove NA distances
  chr.data = chr.data[!is.na(chr.data$gen),]
  chr.data = chr.data[!is.na(chr.data$phys),]
  
  # Interpolate with loess regression
  degree = 2
  
  # Span parameter already calibrated
  span_table = read.table("output/recombination_maps/span.txt", header = TRUE)
  span = span_table$span[which(span_table$set == set & span_table$map == chromosome)]
  
  # Partitioning of the physical map (bp): N bins of equal size
  phys.bins = seq(0, max(chr.data$phys, na.rm = TRUE), len = N) # A sequence of bins coordinates
  
  # Fit the interpolation curve on the data
  fit = fit.loess(x = chr.data, span = span, degree = degree)
  # Bin loess map
  predicted.physbins = predict(fit, newdata = phys.bins) # Predict the local recombination rate at the exact physical positions of the Marey map
  
  derivates = diff(predicted.physbins)/diff(phys.bins)
  
  # We have a problem with incomplete maps (large missing data in one or both tips)
  # and maps with strong decrease in recombination rates in telomeric regions (seems to be a biological process)
  # Hence, when we know that the chromosome is not telocentric (i.e. C.I. > 0.1)
  # we can use this prior information to remove telomeric regions from the derivates
  
  # Check if C.I. is < 0.1; if this is the case, no correction; otherwise, remove derivates on the tips of the chromosome to avoid errors
  # Remove 10% on both sides
  # If no prior information on C.I. (missing data), then remove 10% on both sides because most chromosome are metacentric in plants
  chromosome.stats = read.table(file = "tables/chromosome.stats.csv", header = TRUE, sep = ";")
  centromeric_index = chromosome.stats$centromeric_index[which(chromosome.stats$chromosome == as.character(chromosome) & chromosome.stats$set == set)]
  if (is.na(centromeric_index) | centromeric_index > 0.1) {
    derivates[1:(0.1*N)] = NA
    derivates[(length(derivates) - 0.1*N):length(derivates)] = NA
  }
  
  # The approximate point position of the centromere is the middle of the bin
  centromere_position = sum(phys.bins[(which.min(derivates) - 1):which.min(derivates)])/2
  
  # Add the orientation corrected centromeric index too
  # If the ratio with centromere position is lower than 0.5,
  # then the CI and the centromere inferred are oriented the same way
  centromere_proportion = centromere_position / max(chromosome.stats$phys.map.length[which(chromosome.stats$chromosome == as.character(chromosome) & chromosome.stats$set == set)]/scale, na.rm = TRUE)
  if (centromere_proportion < 0.5) {
    CI = chromosome.stats$centromeric_index[which(chromosome.stats$chromosome == as.character(chromosome) & chromosome.stats$set == set)]*chromosome.stats$phys.map.length[which(chromosome.stats$chromosome == as.character(chromosome) & chromosome.stats$set == set)]/scale
  } else {
    CI = (1 - chromosome.stats$centromeric_index[which(chromosome.stats$chromosome == as.character(chromosome) & chromosome.stats$set == set)])*chromosome.stats$phys.map.length[which(chromosome.stats$chromosome == as.character(chromosome) & chromosome.stats$set == set)]/scale
  }
  
  # Plot diagnostic figure
  marey_map = read.table(paste("output/recombination_maps/loess/100kbwind/", set,"_chromosome", chromosome, ".txt", sep = ""),
                         header = TRUE)
  chromosome.stats = read.table(file = "tables/chromosome.stats.csv", header = TRUE, sep = ";")
  png(paste("output/centromere_position/diagnostic_plots/", set,"_chr", chromosome,".png",sep = ""),
      width = 1200, height = 600)
  plot(marey_map$phys, marey_map$rec.rate, type="l", ylim = c(0, 30), main = paste("Recombination landscape of chromosome ", chromosome, " in ", set, sep = ""), xlab = "Physical position (Mb)", ylab = "Recombination rate (cM/Mb)") #check
  lines(marey_map$phys, marey_map$upper, type = "l", lty = 2, col = "Grey")
  lines(marey_map$phys, marey_map$lower, type = "l", lty = 2, col = "Grey")
  abline(h = 0, col = "Red", lty = 2) # Print a threshold for rec. rate = 0
  # Add the centromere position
  abline(v = centromere_position, col = "Black", lty = 1, lwd = 2)
  # Add the Centromeric Index oriented
  abline(v = CI, col = "Black", lty = 4, lwd = 2)
  dev.off()
  
  # End of function
  return(centromere_position)
}
