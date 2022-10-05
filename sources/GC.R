#============================================================================#
# ESTIMATE GC CONTENT AT A GENE LEVEL ----
#============================================================================#
# GFF-based files
# with supplementary columns
# col. 10 = exon rank (if pertinent)
# col. 11 = gene id
# col. 12 = recombination rate
# col. 13-16= GC/GC1/GC2/GC3

# TODO optim, use shell/python
# Look at https://gist.github.com/darencard/9497e151882c3ff366335040e20b6714
# for a method to extract CDS sequences based on parent features

# Produce one file for each dataset, saved in the "gc_genes" directory
GC_content = function(map = "", ncores = 8, parallel = TRUE) {
  cat("Analysing ", map, ".\n", sep = "")
  
  metadata.clean = read.table("data-cleaned/Genetic_maps_ressources.csv", header = TRUE, sep = ";")
  species = as.character(metadata.clean$species[which(metadata.clean$id == map)])
  accession = as.character(metadata.clean$accession[which(metadata.clean$id == map)])
  # Packages
  require(ape) # dependency to read GFF files with read.gff()
  require(seqinr)
  require(rBLAST)
  require(Biostrings)
  require(parallel)
  require(MASS)
  require(pbmcapply)
  
  # Detect the number of cores to use in parallel tasks
  # Either the max number specified by user of actual number of cores in hardware - 1
  numcores = detectCores()
  ncores = min(numcores - 1, ncores)
  
  if (dir.exists(paste(wd, "data/Genome", tolower(species), accession, "gff3",sep = "/"))) {
    #--------------------------------#
    # DATA
    cat("Read data.\n")
    # Find the path to data
    path = paste(wd, "data/Genome", tolower(species), accession, "gff3",sep = "/")
    # path
    # Retrieve list of chromosomes
    # gff files are saved in .gff.gz with usually the pattern 'chromosome.[A-Za-z0-9]+.' in the directory /gff3
    # Get the list of gff.gz files in the path and keep only chromosomes
    list = system(paste("ls ", path, sep = ""), intern = TRUE)
    list = list[grep('chromosome.[A-Za-z]*[0-9]*', list)]
    list = list[grep('.gff[3]*.gz', list)]
    # Yet, if no individual files for chromosomes, select the 'genomic.gff.gz'
    if (length(list) == 0) {
      list = system(paste("ls ", path, sep = ""), intern = TRUE)
      list = list[grep('genomic.gff[3]*.gz', list)]
    }
    # Eventually, if none of above selections returned a valid gff file, try to get all gff files in the directory
    if (length(list) == 0) {
      list = system(paste("ls ", path, sep = ""), intern = TRUE)
      list = list[grep('.gff[3]*.gz', list)]
    }
    # Finally, if no gff file, then check for gbff in 'Genbank' (this will be the case for most NCBI files)
    # if (length(list) == 0) {
    #   # New path
    #   path = paste(wd, "data/Genome", tolower(species), accession, "genbank",sep = "/")
    #   list = system(paste("ls ", path, sep = ""), intern = TRUE)
    #   list = list[grep('.gbff*.gz', list)]
    # }
    # list
    # Concatenate all GFF files in a single 'data' data frame
    data = data.frame(seqid = character(0), source  = character(0), type = character(0), start  = character(0),
                      end = character(0), score = character(0), strand = character(0), phase = character(0),
                      attributes = character(0))
    # Make paths
    paths = paste(path, list, sep = "/")
    # Read all chromosomes at once
    # gb = readGenBank(list) # too long
    data = lapply(paths, function(x){read.gff(file = x)})
    data = do.call("rbind", data)
    
    #--------------------------------#
    # Remove data that is not chromosome, gene, CDS or exon ----
    data = subset(data, (data$type == "chromosome" | data$type == "gene" | data$type == "exon" | data$type == "CDS"))

    #--------------------------------#
    cat("Retrieve feature metadata...\n")
    # e.g. exon rank, gene id...
    data$rank = NA
    data$id = NA
    data$recombination_rate = NA
    data$gc = NA
    data$gc1 = NA
    data$gc2 = NA
    data$gc3 = NA
    
    data$id = gsub(";", "",
                   gsub("^[A-Za-z_]+=[A-Za-z]+[:|\\-]", "",
                        regmatches(data$attributes,
                                   gregexpr("^[A-Za-z_]+=[A-Za-z0-9|_\\-:.]*;", data$attributes, perl = TRUE)), perl = TRUE))
    
    # Retrieve queried informations and format output
    # GENE IDs
    # data$id[which(data$type == "gene")] = gsub(";", "", gsub("ID=gene[:-]", "",
    #                                                          regmatches(data$attributes[which(data$type == "gene")],
    #                                                                     gregexpr("ID=gene:[A-Za-z0-9_]+;",
    #                                                                              data$attributes[which(data$type == "gene")]))))
    # # CDS IDs
    # data$id[which(data$type == "CDS")] = gsub(";", "", gsub("ID=CDS[:-]", "",
    #                                                         regmatches(data$attributes[which(data$type == "CDS")],
    #                                                                    gregexpr("ID=CDS:[A-Za-z0-9\\.]+;",
    #                                                                             data$attributes[which(data$type == "CDS")]))))
    

    
    #--------------------------------#
    # FILTER DATASET: KEEP ONLY CHROMOSOMES ----
    #--------------------------------#
    # Remove scaffolds, mitochondrial and other non chromosomal parts of the genome
    # All seqid that are not chromosome names, i.e. not matching "[A-Z]{0-1}[0-9]+"
    # unique(data$seqid)
    data = subset(data, grepl("[A-Z]{0-1}[0-9]+", data$seqid))
    
    # EXON IDs and RANK
    # Exception for particular case
    # For some species, exon_id is not working
    # if (accession == "GCA_004115385.1") {
    #   data$id[which(data$type == "exon")] = gsub(";", "",
    #                                              regmatches(data$attributes[which(data$type == "exon")],
    #                                                         gregexpr("HF[A-Za-z0-9]+-[A-Za-z0-9]+-[A-Za-z0-9]+;",
    #                                                                  data$attributes[which(data$type == "exon")])))
    # } else { # General case
    #   data$id[which(data$type == "exon")] = gsub(";", "", gsub("exon_id=", "",
    #                                                            regmatches(data$attributes[which(data$type == "exon")],
    #                                                                       gregexpr("exon_id=[A-Za-z0-9[:punct:]]+;",
    #                                                                                data$attributes[which(data$type == "exon")]))))
    # }
    # EXON RANK
    # data$rank[which(data$type == "exon")] = as.character(regmatches(data$id[which(data$type == "exon")],
    #                                                                      gregexpr("[0-9]+$", data$id[which(data$type == "exon")])))
    
    # Search exon rank, from the most accessible to the most difficult annotation way...
    # If exon is directly encoded in attributes, then retrieve it directly
    if (sum(grepl("exon", unique(data$type))) > 0) {
      exon_exists = TRUE
    } else {
      exon_exists = FALSE
    }
    
    if (exon_exists == FALSE) {
      cat("Annotation file does not contain exons.\n")
    } else {
      if (sum(grepl("exon=", data$attributes[which(data$type == "exon")]))/sum(data$type == "exon") > 0.95 | sum(grepl("rank=", data$attributes[which(data$type == "exon")]))/sum(data$type == "exon") > 0.95) {
        data$rank[which(data$type == "exon")] = as.character(gsub("rank=", "", regmatches(data$attributes[which(data$type == "exon")],
                                                                                               gregexpr("rank=[0-9]+", data$attributes[which(data$type == "exon")]))))
      } else {
        # If most exons are encoded in the ID (threshold of 95% matching the condition), with pattern like 'id-rank', then retrieve rank form id
        # exon rank is the number at the end of ID
        if (sum(grepl("[.\\-][0-9]+$", data$id[which(data$type == "exon")], perl = TRUE))/sum(data$type == "exon") > 0.95) {
          data$rank[which(data$type == "exon")] = as.character(regmatches(data$id[which(data$type == "exon")],
                                                                               gregexpr("[0-9]+$", data$id[which(data$type == "exon")], perl = TRUE)))
        }
      }
    }
    
    # if (accession == "GCA_004115385.1") {
    #   data$rank[which(data$type == "exon")] = gsub(";", "",
    #                                                     gsub("HF[A-Za-z0-9]+-[A-Za-z0-9]+-", "",
    #                                                          regmatches(data$attributes[which(data$type == "exon")],
    #                                                                     gregexpr("HF[A-Za-z0-9]+-[A-Za-z0-9]+-[A-Za-z0-9]+;",
    #                                                                              data$attributes[which(data$type == "exon")]))))
    # } else { # General case
    #   # data$rank[which(data$type == "exon")] = as.character(gsub(";rank=", "", regmatches(data$attributes[which(data$type == "exon")],
    #   #                                                                                         gregexpr(";rank=[0-9]+", data$attributes[which(data$type == "exon")]))))
    #   data$rank[which(data$type == "exon")] = as.character(regmatches(data$id[which(data$type == "exon")],
    #                                                                                           gregexpr("[0-9]+$", data$id[which(data$type == "exon")])))
    # }
    # 
    # Yet sometimes id is still not recovered properly or not indicated in gff
    if (length(which(as.character(data$id) == "character(0)")) > 0) {
      data$id[which(as.character(data$id) == "character(0)")] = NA
      warning("Errors in IDs: NA produced.")
    }
    
    
    #--------------------------------#
    # Retrieve CDS rank ----
    #--------------------------------#
    cat("Retrieve CDS rank...\n")
    # CDS must be numbered by their position from 5' to 3' in the gene
    # Independently of exon rank, since all gff does not have exons annotated
    # For each feature (i.e. row), number if a CDS, otherwise return actual value
    # A function that takes a data frame and a row number x as argument
    
    # Alternative variant (splicing) must be taken into account
    # Use the parent information in Attributes if available
    
    # TODO Optimize: time too long with splicing variants version
    # Subset data as soon as possible
    
    get_CDS_rank = function(data, x) {
      # Is it a CDS?
      if (data$type[x] == "CDS") {
        # Subset data as soon as possible
        # # First at a +/-10kb scale
        # subset_gene = subset(data, data$start >= (data$start[x] - 10000) &
        #                        data$end <= (data$end[idx_gene] + 10000) & data$seqid == data$seqid[idx_gene])
        # Get gene coordinates
        idx_gene = which(data$type == "gene" & data$start <= data$start[x] & data$end >= data$end[x]
                         & data$seqid == data$seqid[x])
        # If more than one gene, select the one that has the same id in attributes
        # Yet if no attributes to make the link, select the gene index just before the row x 
        if (length(idx_gene) > 1 & length(which(unlist(lapply(idx_gene, function(i) grepl(data$id[i], data$attributes[x]))))) > 0) {
          idx_gene = idx_gene[which(unlist(lapply(idx_gene, function(i) grepl(data$id[i], data$attributes[x]))))]
        } else {
          # Select only idx before x, then the max idx position, which is the closest to x
          idx_gene = idx_gene[which.min(idx_gene[which((idx_gene - x) < 0)])]
        }
        subset_gene = subset(data, data$start >= data$start[idx_gene] &
                               data$end <= data$end[idx_gene] & data$seqid == data$seqid[idx_gene])
        # if gene id is present in attributes of exons/CDS (useful with multiple splicing variants)
        if (sum(grepl(data$id[idx_gene], subset_gene$attributes)) > 1) {
          subset_gene = subset(subset_gene, grepl(data$id[idx_gene], subset_gene$attributes))
        }
        # Associate CDS to the parent exon if available
        # Otherwise manually count CDS rank (beware, don't take into account splicing variants)
        if (grepl("[Pp]arent=", data$attributes[x])) {
          # Identify the parent information (i.e. splicing variant)
          parent = as.character(gsub(";$", "", gsub(";[Pp]arent=", "", regmatches(data$attributes[x], gregexpr(";[Pp]arent=[A-Za-z0-9:\\-_.]*;", data$attributes[x], perl = TRUE)))))
          # How many CDS before (strand +) or after (strand -)?
          if (data$strand[idx_gene] == "+") {
            cds_rank = length(which(subset_gene$type == "CDS" & subset_gene$start >= data$start[idx_gene] &
                                      subset_gene$end <= data$end[x] & subset_gene$seqid == data$seqid[x] &
                                      grepl(parent, subset_gene$attributes)))
          } else {
            cds_rank = length(which(subset_gene$type == "CDS" & subset_gene$start >= data$start[x] &
                                      subset_gene$end <= data$end[idx_gene] & subset_gene$seqid == data$seqid[x] &
                                      grepl(parent, subset_gene$attributes)))
          }
        } else {
          # How many CDS before (strand +) or after (strand -)?
          if (data$strand[idx_gene] == "+") {
            cds_rank = length(which(subset_gene$type == "CDS" & subset_gene$start >= data$start[idx_gene] & subset_gene$end <= data$end[x] & subset_gene$seqid == data$seqid[x]))
          } else {
            cds_rank = length(which(subset_gene$type == "CDS" & subset_gene$start >= data$start[x] & subset_gene$end <= data$end[idx_gene] & subset_gene$seqid == data$seqid[x]))
          }
        }
      } else {
          cds_rank = data$rank[x]
        }
      return(cds_rank)
      # END OF FUNCTION
    }
    
    # DEBUG ONLY
    # for (i in 1:nrow(data)) {
    #   cat(i, "\n")
    #   get_CDS_rank(data, i)
    # }

    if (parallel == TRUE) {
      # CDS_ranks = unlist(pbmclapply(X = 1:100, function(x) get_CDS_rank(data, x), mc.cores = ncores))
      CDS_ranks = unlist(pbmclapply(X = 1:nrow(data), function(x) get_CDS_rank(data, x), mc.cores = ncores))
    } else {
      CDS_ranks = unlist(pbmclapply(X = 1:nrow(data), function(x) get_CDS_rank(data, x), mc.cores = 1))
    }
    # data$rank[1:100] = CDS_ranks
    data$rank = CDS_ranks
    rm(CDS_ranks)
    # Free unused memory
    gc()
    
    # I tried to optim memory usage with foreach, but same memory overflow problem...
    # library(foreach)
    # library(doMC)
    # registerDoMC(ncores)
    # iterations = nrow(data)
    # pb = txtProgressBar(max = iterations, style = 3)
    # progress = function(n) setTxtProgressBar(pb, n)
    # opts = list(progress = progress)
    # CDS_ranks = foreach(i = 1:iterations, .options.snow = opts) %dopar% {
    #   get_CDS_rank(data, i)
    # }
    # close(pb)
    # CDS_ranks = unlist(CDS_ranks)
    
    # # The list of genes
    # list_genes = as.character(unique(data$id[which(data$type == "CDS")]))
    # # For each genes, number CDS b their order of appearance in the gene
    # # And it depends on the direction of the reading frame (strand +/-)
    # cat("Retrieve rank of CDS.\n")
    # # pb = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
    # #                        total = length(list_genes),
    # #                        complete = "=",   # Completion bar character
    # #                        incomplete = "-", # Incomplete bar character
    # #                        current = ">",    # Current bar character
    # #                        clear = FALSE,    # If TRUE, clears the bar when finish
    # #                        width = 100)      # Width of the progress bar
    # pb = txtProgressBar(min = 0, max = length(list_genes), initial = 0)
    # for (i in 1:length(list_genes)) {
    #   setTxtProgressBar(pb, i)
    #   # pb$tick()
    #   # cat(i, "\n")
    #   # Get CDS of a gene i
    #   idx = which(grepl(list_genes[i], data$id) & data$type == "CDS")
    #   cds_data = data[idx,]
    #   cds_data$idx = idx
    #   # Order them by position and strand
    #   if (unique(cds_data$strand) == "+") {
    #     cds_data = cds_data[order(cds_data$start, decreasing = FALSE),]
    #   } else {
    #     cds_data = cds_data[order(cds_data$start, decreasing = TRUE),]
    #   }
    #   # Save CDS ranks in data
    #   data$rank[cds_data$idx] = seq(1, nrow(cds_data))
    # }
    # close(pb)
    
    #--------------------------------#
    # CREATE INTRONS ----
    #--------------------------------#
    cat("Creating introns.\n")
    
    
    
    
    #--------------------------------#
    # ADD 100KB WINDOWS OF RECOMBINATION ----
    #--------------------------------#
    cat("Adding windows of recombination (100kb).\n")
    # Load recombination maps (Marey maps)
    # Read the recombination map of the associated dataset
    for (chr in unique(chromosome.stats$chromosome[which(chromosome.stats$set == map)])) {
      if (file.exists(paste("output/recombination_maps/loess/100kbwind/", map, "_chromosome", chr ,".txt", sep = ""))) {
        rec_map = read.table(paste("output/recombination_maps/loess/100kbwind/", map, "_chromosome", chr ,".txt", sep = ""),
                             header = TRUE, stringsAsFactors = FALSE)
        # index must begin at one
        new_data = data.frame(seqid = chr, source = "Marey_map", type = "100kbWindows", start = ((rec_map$phys - 0.05)*10^6 + 1),
                              end = ((rec_map$phys + 0.05)*10^6 + 1), score = NA, strand = "+", phase = NA,
                              attributes = NA, rank = NA,
                              id = NA, recombination_rate = rec_map$rec.rate, gc = NA, gc1 = NA, gc2 = NA, gc3 = NA)
        data = rbind(data, new_data)
      }
    }
    rm(new_data)
    rm(rec_map)

    # Add windows as features of type 'windows' with start-end coordinates
    
    
    
    
    #--------------------------------#
    # ADD FLANKING REGIONS ----
    #--------------------------------#
    cat("Adding flanking regions (5kb upstream/downstream the complete transcript).\n")
    # Flanking regions are associated to genes (i.e. transcripts),
    # thus they must have the exact same id as genes to perform correlations
    
    #--------------------------------#
    # TAKE CARE OF CHROMOSOME NAMES - TRANSLATE...  ---- 
    #--------------------------------#
    
    # TODONE Optimization
    # Use factors to translate chromosome names
    
    cat("Translate chromosome names...\n")
    # Retrieve factor levels
    data$seqid = as.factor(data$seqid)
    chr_levels = levels(data$seqid)
    save_levels = levels(data$seqid)
    
    # generic filtering
    # data$seqid = as.character(data$seqid)
    chr_levels = gsub("chr0", "", chr_levels)
    chr_levels = gsub("chr", "", chr_levels)
    chr_levels = gsub("Chr0", "", chr_levels)
    chr_levels = gsub("Chr", "", chr_levels)
    chr_levels = gsub("LG0", "", chr_levels)
    chr_levels = gsub("LG", "", chr_levels)
    # Solanum tuberosum
    chr_levels = gsub("ST4.03ch[0]*", "", chr_levels)
    # Gossypium hirsutum
    chr_levels = gsub("Ghir_", "", chr_levels)
    # Cucurbita maxima
    chr_levels = gsub("Cma_", "", chr_levels)
    
    # Converting from a translation table
    table_translate = read.table("data/chromosome_translate.csv", header = TRUE, sep = " ")
    # TODO optimisation
    if (parallel == FALSE) {
      cat("Parallel = FALSE\n")
      chrnames = mclapply(X = 1:length(chr_levels), function(x){if (chr_levels[x] %in% table_translate$annotation_name) {
        chr_levels[x] = unique(table_translate$dataset_name[which(table_translate$annotation_name == chr_levels[x])])
      } else {chr_levels[x] = chr_levels[x]}}, mc.cores = 1)
      chr_levels = unlist(cbind(chrnames))
    } else {
      chrnames = mclapply(X = 1:length(chr_levels), function(x){if (chr_levels[x] %in% table_translate$annotation_name) {
        chr_levels[x] = unique(table_translate$dataset_name[which(table_translate$annotation_name == chr_levels[x])])
      } else {chr_levels[x] = chr_levels[x]}})
      chr_levels = unlist(cbind(chrnames))
    }
    # Free unused memory
    gc()
    # Treating exceptions
    # Capsella rubella have pseudochromomes with scaffold names
    if (accession == "GCA_000695605.1") {
      chr_levels = gsub("scaffold[0]+", "", chr_levels)
      data = data[as.numeric(chr_levels) <= 9,]
    }
    # A. speltoides S genome is mapped on A. tauschii D genome
    if (species == "Aegilops_speltoides" & accession == "GCA_000347335.2") {
      chr_levels = gsub("D", "S", chr_levels)
    }
    # replace by the translated level
    data$seqid = factor(data$seqid, levels = save_levels, labels = chr_levels)
    rm(table_translate)

    
    #--------------------------------#
    # ESTIMATE GC CONTENT ----
    #--------------------------------#
    # GC content must be done before changing chromosome names, to keep reference between annotation and fasta files
    cat("Estimate GC content...\n")
    # Retrieve GC content for each feature
    # the get_gc() function estimates GC/GC1/GC2/GC3 contents in the features provided
    # Take a data frame of dataset name, accession, chromosomes and positions (start-end), 
    # and return a data frame of values (4 columns)
    
    # Beware that...
    # "For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame.
    # The phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the beginning
    # of this feature to reach the first base of the next codon.
    # In other words, a phase of "0" indicates that the next codon begins at the first base of the region described by the current line,
    # a phase of "1" indicates that the next codon begins at the second base of this region,
    # and a phase of "2" indicates that the codon begins at the third base of this region.
    # This is NOT to be confused with the frame, which is simply start modulo 3.
    # If there is no phase, put a "." (a period) in this field.
    # For forward strand features, phase is counted from the start field.
    # For reverse strand features, phase is counted from the end field." (http://seqanswers.com/forums/showthread.php?t=9500)
    # The phase is required for all CDS features.
    # data$seqid = gsub("^[Cc]hr", "", data$seqid)
    # data$seqid = gsub("^0", "", data$seqid)
    # data$seqid = gsub("^LG[0]*", "", data$seqid)
    pos = data.frame(dataset = map, accession = accession,
                     seqid = data$seqid,
                     start = data$start, end = data$end,
                     frame = data$phase, strand = data$strand, type = data$type)
    
    # A function to return the four GC proportions for each line in the gff
    get_gc = function(pos) {
      # Load fasta file
      # Names of fasta files susceptible to contain sequence data
      ls = list.files(path = paste(wd, "/data/Genome/", tolower(species),"/" ,accession,"/fasta/", sep = ""))
      ls = ls[grepl("fna.gz", ls) | grepl("fa.gz", ls) | grepl("fasta.gz", ls)]
      
      # If no fasta file found, i.e. no genomic data
      if (length(ls) == 0) {
        return("No fasta file for this dataset.")
      } else {  # Open the fasta file and concatenate them if multiple fasta files: make a single data for all chromosomes
        seq = readDNAStringSet(paste(wd, "/data/Genome/", tolower(species),"/" ,accession,"/fasta/", ls, sep = ""))
      }
      # Chromosome names
      # Either they are the same in fasta and gff, or they must be formatted prior to feature request
      chrnames = names(seq)
      translate_chrnames = function(chrnames) {
        for (i in 1:length(chrnames)) {
          # remove any literal naming (e.g. chr, Chr, chromosome)
          # generic filtering
          chrnames[i] = gsub("[Cc]hromosome", "", chrnames[i])
          chrnames[i] = gsub("[Cc]hr", "", chrnames[i])
          chrnames[i] = gsub("LG[0]{0-1}", "", chrnames[i])
          chrnames[i] = gsub("^0", "", chrnames[i])
          # Solanum tuberosum
          chrnames[i] = gsub("ST4.03ch[0]*", "", chrnames[i])
          # Gossypium hirsutum
          chrnames[i] = gsub("Ghir_", "", chrnames[i])
          # Cucurbita maxima
          chrnames[i] = gsub("Cma_", "", chrnames[i])
          # Keep the first element only (line begins by the chromosome name usually)
          chrnames[i] = strsplit(chrnames[i], split = " ")[[1]][1]
          
          # Deal directly this part with accession
          # if it is not a valid chromosome name, search a translation
          # Converting from a translation table
          if (accession %in% c("GCA_000001735.2", "GCA_000004075.2", "GCA_002240015.2", "GCA_002285895.2", "GCA_001433935.1")) {
            table_translate = read.table("data/chromosome_translate.csv", header = TRUE, sep = " ",
                                         stringsAsFactors = FALSE)
            if (chrnames[i] %in% table_translate$annotation_name){
              chrnames[i] = unique(table_translate$dataset_name[which(table_translate$annotation_name == chrnames[i])])
            }
          }
        }
        return(chrnames)
      }
      names(seq) = translate_chrnames(chrnames)
      
      # Estimate GC for the whole dataset
      if (parallel == FALSE) {
        cat("Parallel = FALSE\n")
        gc_feature = pbmclapply(X = 1:nrow(pos),
                                function(x) get_gc_feature(seq, chr = as.character(pos$seqid[x]), start = pos$start[x], end = pos$end[x],
                                                           type = pos$type[x], 
                                                           frame = as.character(pos$frame[x]), strand = as.character(pos$strand[x])),
                                mc.cores = 1)
        gc_feature = do.call(rbind, gc_feature)
      } else {
        gc_feature = pbmclapply(X = 1:nrow(pos),
                                function(x) get_gc_feature(seq, chr = as.character(pos$seqid[x]),
                                                           start = pos$start[x], end = pos$end[x],
                                                           type = as.character(pos$type[x]), 
                                                           frame = as.character(pos$frame[x]),
                                                           strand = as.character(pos$strand[x])),
                                mc.cores = ncores)
        gc_feature = do.call(rbind, gc_feature)
      }
      # DEBUG ONLY
      # for (i in 1:nrow(pos)) {
      #   cat(i, ":", i/nrow(pos)*100, "%\n")
      #   get_gc_feature(seq, chr = pos$seqid[i], start = pos$start[i], end = pos$end[i],
      #                  frame = pos$frame[i], strand = pos$strand[i])
      # }
      
      return(gc_feature)
    }
    
    # a small function to estimate GC proportions for a single feature with start-end positions
    get_gc_feature = function(seq, chr = "", start = "", end = "", frame = "", strand = "", type = "") {
      frame = as.numeric(frame)
      strand = as.character(strand)
      type = as.character(type)
      # Subset to chromosome
      seqwin = seq[max(which(names(seq) == chr))] 
      gc = NA
      gc1 = NA
      gc2 = NA
      gc3 = NA
      
      if (end < seqwin@ranges@width) { # If the feature is out of the boundaries of the reference genome
        # Subset to start-end positions
        seqwin = subseq(seqwin, start, end)
        # NAs produced in tails of the sequence
        seqwin[is.na(seqwin)] = "N"
        seqwin = s2c(tolower(as.character(seqwin))) # Convert to a vector of single characters
        # NA frames are set to 0 (non-CDS)
        frame[is.na(frame)] = 0

        # IF reverse strand (strand = -)
        if (!is.na(strand)) {
          if (strand == "-") {
            seqwin = rev(seqwin)
          }
          # Some sequences are very short... Force to NA if less than 10 nucleotides
          # Besides, some sequences are not the correct length with the reading frame
          if (length(seqwin) > 10 & (length(seqwin) - frame) >= 3) {
            # forcetoLower = FALSE require to convert sequence to lower case before, yet it reduces the computation time by 3
            # exact computation estimates GC content with ambiguous bases taken in account. Time is increased by 3.
            # Yet, it revealed differences between estimations with exact = FALSE
            gc = GC(seqwin, forceToLower = FALSE, exact = TRUE)
            if (type == "CDS" | type == "exon") { # GC1/2/3 make sense only for CDS/exons; save time
              gc1 = GC1(seqwin, frame = frame, forceToLower = FALSE, exact = TRUE)
              gc2 = GC2(seqwin, frame = frame, forceToLower = FALSE, exact = TRUE)
              gc3 = GC3(seqwin, frame = frame, forceToLower = FALSE, exact = TRUE)
            }
          }
        }
      }
      
      gc_feature = data.frame(gc = gc, gc1 = gc1, gc2 = gc2, gc3 = gc3)
      return(gc_feature)
    }
    # Estimate
    gc_content = get_gc(pos)
    data$gc = gc_content$gc
    data$gc1 = gc_content$gc1
    data$gc2 = gc_content$gc2
    data$gc3 = gc_content$gc3
    rm(gc_content)
    rm(pos)
    # Free unused memory
    gc()
    #--------------------------------#
    # CpG CONTENT ----
    #--------------------------------#
    cat("Estimate CpG (%) in 100kb windows, exons/CDS (5'-3') introns and flanking regions...\n")
    
    data$CpG = NA
    
    
    #--------------------------------#
    # RECOMBINATION RATES ----
    #--------------------------------#
    cat("Retrieve recombination rates...\n")
    # After changing chromosome names to the one used in Marey maps, retrieve recombination rates
    # Get the recombination rate where the feature start
    # Need only a start position and a chromosome name in argument, with a dataset name
    # and return the recombination rate of the 100kb window in which the start position is
    
    # TODO Optim: read.table only once
    
    get_recRate = function(pos, chr, rec_map) {
        rec_map = subset(rec_map, rec_map$chr == chr)
        # Query the local recombination rate of a feature by its physical position
        # i.e. windows in which start position is
        # And in any case, prevent for NA estimates out of boundaries
        if ((pos < min((rec_map$phys*10^6))) | (pos > max((rec_map$phys*10^6)))) {
          recrate = NA
        } else {
          # windows beginning just before the position
          recrate = rec_map$rec.rate[max(which(pos > (rec_map$phys*10^6)))]
        }
      return(recrate)
    }
    
    # Ignore 100kb windows: already have a recombination rate
    data_subset_one = subset(data, data$type != "100kbWindows")
    data_subset_two = subset(data, data$type == "100kbWindows")
    gc()
    # Load data once
    path = paste(wd, "output/recombination_maps/loess/100kbwind", sep = "/")
    list = system(paste("ls ", path, sep = ""), intern = TRUE)
    list = list[grep(map, list)]
    list = gsub(".txt", "", list)
    list = gsub("[A-Za-z0-9_]*_chromosome", "", list)
    # Read all chromosomes at once
    rec_map = lapply(list,
                     function(x){tmp = read.table(file = paste(wd, "/output/recombination_maps/loess/100kbwind/", map, "_chromosome", x, ".txt", sep = ""),
                                                  header = TRUE, stringsAsFactors = FALSE)
                     tmp$chr = x
                     return(tmp)})
    rec_map = do.call("rbind", rec_map)
    # Get recombination rates
    rec_rates = pbmclapply(X = 1:nrow(data_subset_one),
                           function(x) get_recRate(data_subset_one$start[x], data_subset_one$seqid[x], rec_map),
                           mc.cores = ncores)
    data_subset_one$recombination_rate = unlist(rec_rates)
    rm(rec_rates)
    # Free unused memory
    gc()
    # Reunite,datasets
    # data = rbind(data_subset_one, data_subset_two)
    
    cat("Saving files.\n")
    # Save in two different files
    write.table(data_subset_one, file = gzfile(paste("data-cleaned/genome/gc_genes/", map, "_gcgenes.txt.gz", sep = "")),
                row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    write.table(data_subset_two, file = gzfile(paste("data-cleaned/genome/gc_map/", map, "_gcmap.txt.gz", sep = "")),
                row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    
    cat("Finished.\n")
    # Free unused memory
    gc()
    return(data)
    # END OF 'IF DIR.EXIST'  
  } else {
    warning("No genome or annotation found !")
  }
  # Free unused memory
  gc()
  # END OF FUNCTION
}

