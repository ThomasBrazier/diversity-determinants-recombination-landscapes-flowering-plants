#============================================================================#
# Estimate the measure of genetic shuffling as described by Veller et al. (2019)  ----
#============================================================================#
R_intra = function(data, set = "", chromosome = "") {
  metadata.clean = read.csv(file = paste(wd, "/data-cleaned/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")
  # data = read.table("data-cleaned/veller/GeneticDistances.txt", header = TRUE)
  d = data[which(data$set == set & data$map == chromosome),]
  # Chromosome length
  chrlength = max(data$phys[which(data$set == set & data$map == chromosome)], na.rm = TRUE)
  # Total genome length
  totalchrlength = 0
  listchr = unique(data$map[which(data$set == set)])
  for (chr in listchr) {
    totalchrlength = totalchrlength + max(data$phys[which(data$set == set & data$map == chr)], na.rm = TRUE)
  }
  # Compute inter-loci distances for each locus pair (i,j) in a matrix
  dij = matrix(NA, nrow = nrow(d), ncol = nrow(d))
  for (i in 1:nrow(d)) {
    for (j in 1:nrow(d)) {
      dij[i,j] = d$gen[j] - d$gen[i]
    }
  }
  # 0 distances in the diagonal passed to NA
  dij[which(dij == 0)] = NA
  # hist(dij)
  
  # Convert distances back to recombinant fractions
  if (is.na(metadata.clean$mapping_function[which(metadata.clean$id == set)])) {
    # reversed Morgan, converting Morgans to recombinant fraction
    cat("Converting Morgans to recombinant fraction\n")
    rij = dij/100
    # cat(mean(rij, na.rm = TRUE))
  } else {
    if (metadata.clean$mapping_function[which(metadata.clean$id == set)] == "haldane") {
      # Reversed Haldane, converting cM to recombinant fraction
      cat("Reversed Haldane, converting cM to recombinant fraction\n")
      rij = 0.5*(1-exp(-2*dij/100))
      # cat(mean(rij, na.rm = TRUE))
    } else {
      if (metadata.clean$mapping_function[which(metadata.clean$id == set)] == "kosambi") {
        # Reversed Kosambi, converting cM to recombinant fraction
        cat("Reversed Kosambi, converting cM to recombinant fraction\n")
        rij = 0.5*tanh(2*dij/100)
        # cat(mean(rij, na.rm = TRUE))
      }
    }
  }
  
  # Sum of rij for i<j, i.e. the upper triangular matrix rij
  rij[lower.tri(rij)] = NA
  # hist(rij)
  # Intra-chromosomal component is the averaged rate of shuffling rij for each locus pair (i,j)
  # lambda = nrow(d)*(nrow(d)-1)/2 # loci x (loci-1)/2
  lambda = sum(!is.na(rij)) # is the total number of pairs of loci, i.e. number of values not NA
  r_intra = (sum(rij, na.rm = TRUE)/(lambda*(lambda-1)*(1/2)))*(chrlength/totalchrlength)^2
  return(r_intra)
}

# A function that compute the intrachromosomal component of the Rbar measure from a given dataset and returns its value
Rbar_intra = function(data, set) {
  list_chr = unique(data$map[which(data$set == set)])
  rbarintra = 0
  for (chr in list_chr) {
    tmp = R_intra(data, set, chr)
    # print(tmp)
    rbarintra = rbarintra + tmp
  }
  return(rbarintra)
}

#============================================================================#
# Genetic shuffling rate of chromosomes ----
#============================================================================#
# Reversed mapping function: from the genetic distances to the recombination fraction
# Reversed Kosambi from the formula given in Veller et al. 2019
# r_ij = r(d_ij) = 1/2(tanh(2d_ij))

# A function that compute the genetic shuffling rate per chromosomes,
# According to equation 10 in Veller et al. (2019)
# for a given set and chromosome, and return its value
R_chr = function(data, set = "", chromosome = "") {
  cat(as.character(set), "chromosome", chromosome, "\n")
  # data = read.table("data-cleaned/veller/GeneticDistances.txt", header = TRUE)
  d = data[which(data$set == set & data$map == chromosome),]
  metadata.clean = read.table(file = paste(wd, "/data-cleaned/Genetic_maps_ressources.csv", sep = ""),
                              header = TRUE, sep = ";", stringsAsFactors = FALSE)

  # Compute inter-loci distances for each locus pair (i,j) in a matrix
  dij = matrix(NA, nrow = nrow(d), ncol = nrow(d))
  for (i in 1:nrow(d)) {
    for (j in 1:nrow(d)) {
      dij[i,j] = d$gen[j] - d$gen[i]
    }
  }
  # 0 distances in the diagonal passed to NA
  dij[which(dij == 0)] = NA
  # Sum of dij for i<j, i.e. the upper triangular matrix dij
  dij[lower.tri(dij)] = NA
  
  # Convert distances back to recombinant fractions
  if (is.na(metadata.clean$mapping_function[which(metadata.clean$id == set)]) | metadata.clean$mapping_function[which(metadata.clean$id == set)] == "bin") {
    # reversed Morgan, converting Morgans to recombinant fraction
    cat("Converting Morgans to recombinant fraction\n")
    rij = dij/100
    # cat(mean(rij, na.rm = TRUE))
  } else {
    if (metadata.clean$mapping_function[which(metadata.clean$id == set)] == "haldane") {
      # Reversed Haldane, converting cM to recombinant fraction
      cat("Reversed Haldane, converting cM to recombinant fraction\n")
      rij = 0.5*(1-exp(-2*dij/100))
      # cat(mean(rij, na.rm = TRUE))
    } else {
      if (metadata.clean$mapping_function[which(metadata.clean$id == set)] == "kosambi") {
        # Reversed Kosambi, converting cM to recombinant fraction
        cat("Reversed Kosambi, converting cM to recombinant fraction\n")
        rij = 0.5*tanh(2*dij/100)
        # cat(mean(rij, na.rm = TRUE))
      }
    }
  }
  
  # Intra-chromosomal component is the averaged rate of shuffling rij for each locus pair (i,j)
  #lambda = nrow(d)*(nrow(d)-1)/2 # loci x (loci-1)/2
  lambda = nrow(d) # is the total number of loci
  r_chr = (sum(rij, na.rm = TRUE)/(lambda*(lambda-1)/2))
  # r_intra = (sum(rij, na.rm = TRUE)/(lambda*(lambda-1)*(1/2)))*(chrlength/totalchrlength)^2
  return(r_chr)
}

# I had a doubt on the way I implemented the equation(10) of Veller et al. 2019
# Hence I re-implemented the chromosome shuffling rate directly from the Matlab script at disposition
R_chr_fromMatlab = function(data, set = "", chromosome = "") {
  cat(as.character(set), "chromosome", chromosome, "\n")
  # data = read.table("data-cleaned/veller/GeneticDistances.txt", header = TRUE)
  metadata.clean = read.table(file = paste(wd, "/data-cleaned/Genetic_maps_ressources.csv", sep = ""), header = TRUE, sep = ";")
  
  # This loads an n x 6 data frame of inter-SNP map length data where,
  # n = 1000 pseudo-markers
  # and in row i
  # 1 - set
  # 2 - chromosome
  # 3 - marker name
  # 4 - marker physical position
  # 5 - marker genetic position
  # 6 - marker estimated map distance (cM) between the SNP in row i and that in row i-1 (if same chromosome).
  d = data[which(data$set == set & data$map == chromosome),]
  # if (is.na(d$gen[1])) {d$gen[1] = 0}
  d$dist = c(0, diff(d$gen))

  loci = nrow(d) # Number of evenly spaced pseudo-loci
  
  # Distances between all locus pairs
  pairwise_distances = as.matrix(dist(d$gen))
  dim(pairwise_distances)
  # Sum of rij for i<j, i.e. the upper triangular matrix rij
  pairwise_distances[lower.tri(pairwise_distances)] = NA
  # 0 distances in the diagonal passed to NA
  
  pairwise_distances[which(pairwise_distances == 0)] = NA
  
  # Convert distances back to recombinant fractions
  if (is.na(metadata.clean$mapping_function[which(metadata.clean$id == set)]) | metadata.clean$mapping_function[which(metadata.clean$id == set)] == "bin") {
    # reversed Morgan, converting Morgans to recombinant fraction
    cat("Converting Morgans to recombinant fraction\n")
    pairwise_distances = pairwise_distances/100
    # cat(mean(rij, na.rm = TRUE))
  } else {
    if (metadata.clean$mapping_function[which(metadata.clean$id == set)] == "haldane") {
      # Reversed Haldane, converting cM to recombinant fraction
      cat("Reversed Haldane, converting cM to recombinant fraction\n")
      pairwise_distances = 0.5*(1-exp(-2*pairwise_distances/100))
      # cat(mean(rij, na.rm = TRUE))
    } else {
      if (metadata.clean$mapping_function[which(metadata.clean$id == set)] == "kosambi") {
        # Reversed Kosambi, converting cM to recombinant fraction
        cat("Reversed Kosambi, converting cM to recombinant fraction\n")
        pairwise_distances = 0.5*tanh(2*pairwise_distances/100)
        # cat(mean(rij, na.rm = TRUE))
      }
    }
  }
  
  intra_sum = sum(pairwise_distances, na.rm = TRUE) # track the running total of the inra-chromosomal contribution to rbar as we move from locus pair to locus pair

  count_pairs = loci*(loci - 1)/2 # Will keep track of how many pairs we've counted -> loci x (loci-1)/2

  intra_contrib = (intra_sum/count_pairs)
  
  return(intra_contrib)
}



# A function that compute the rbar genetic shuffling rate for a species and returns its value
R_chr_mean = function(data, set) {
  cat(as.character(set), "\n")
  list_chr = unique(data$map[which(data$set == set)])
  r_chr = 0
  for (chr in list_chr) {
    tmp = R_chr(data, set, chr)
    # print(tmp)
    r_chr = r_chr + tmp
  }
  r_chr_mean = r_chr/length(list_chr)
  return(r_chr_mean)
}


#============================================================================#
# Estimation of the local recombination rate (cM/Mb) ----
#============================================================================#
recombinationMap_Veller = function(set ="", data = data.final, chr = "all", K = 5, parallel = TRUE,...) { # Optional arguments passed to loess, smooth.spline and pther fitting methods
  source("sources/MareyMap.R") 
  # Physical distances are in bp, yet analyses are in Mb
  # Convert bp to Mb
  scale = 1000000
  cat("Loading data\n")
  data.final = read.table(file = "data-cleaned/marey_maps/AllMaps.txt", header = TRUE, sep = "\t")
  
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
    data.final$set = as.character(data.final$set)
    data.final$map = as.character(data.final$map)
    chr.data = subset(data.final, map == chr.nb & set == set.name) # Sample only distances of the chromosome
    chr.data$gen = as.numeric(as.character(chr.data$gen))
    chr.data$phys = as.numeric(as.character(chr.data$phys))
    # Convert physical distances in Mb
    chr.data$phys = chr.data$phys/scale
    
    # Remove NA distances
    chr.data = chr.data[!is.na(chr.data$gen),]
    chr.data = chr.data[!is.na(chr.data$phys),]

    degree = 2

    pointwise = sort(chr.data$phys)
    phys.bins = seq(0, max(chr.data$phys, na.rm = TRUE), len = 1000) # A sequence of bins coordinates
    #---------------------------------#
    # Veller method: interpolate 1,000 evenly distributed physical positions
    #---------------------------------#
    # Calibrate span
    # span = calibrate.span(x = chr.data, K)
    # span = calibrate.span(x = chr.data, K, from = 0.2, to = 0.5, parallel = parallel)
    # Use the previously calibrated span parameter
    cat("Loading span\n")
    span = NA
    span_df = read.table("output/recombination_maps/span.txt", header = TRUE)
    if (length(span_df$span[which(span_df$set == set & span_df$map == chr.nb)]) > 0) {
      span = span_df$span[which(span_df$set == set & span_df$map == chr.nb)]  
    }
    if(is.na(span)) {span = 0.25} # default value
    
    # Fit the interpolation curve on data
    cat("Fitting the interpolated curve\n")
    fit = fit.loess(x = chr.data, span = span, degree = degree)
    # Partitioning of the physical map (bp): 1,000 bins of equal size
    predicted.physbins = predict(fit, newdata = phys.bins) # Predict the genetic distance at the exact physical positions of the Marey map
    # A Marey map of 1,000 evenly spaced markers
    df.tmp = read.table("output/veller/GeneticDistances.txt", header = TRUE)
    # Return a data.frame
    df.veller = data.frame(set = rep(set.name, 1000), map = rep(chr.nb, 1000), mkr = paste("chr", chr.nb,"_",1:1000, sep = ""), phys = phys.bins, gen = predicted.physbins)
    # Append new data to replace
    df.tmp = df.tmp[which(!(df.tmp$set == set.name & df.tmp$map == chr.nb)),]
    df.tmp = rbind(df.tmp, df.veller)
    write.table(df.tmp, file = "output/veller/GeneticDistances.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
  }
}


