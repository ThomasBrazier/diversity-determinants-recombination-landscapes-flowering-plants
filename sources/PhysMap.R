# Work with physical maps
# With a given set of maps or all maps in a directory
require(ggplot2)

#============================================================================
# Update the marker DB
#============================================================================
# All markers stored in monospecific databases

#----------------------------------------------------------------------------
# Add new markers to the database and compute NEW UNKNOWN physical positions
# set = species name
# data = path to the marker directory
markerDB.update = function(set = "", data = "data/Physical_maps/Markers/") {
  # Import the handmade marker DB as a xls file
  require(gdata)
  mkr = read.xls (paste(data, set, ".xlsx", sep = ""), sheet = 1, header = TRUE)
  
  # Remove duplicates - only one version of each marker for the same reference genome
  
  # compute NEW UNKNOWN physical positions
  for (i in 1:nrow(mkr)) {
    if (is.na(mkr$phys[i])) {
      
    }
  }
  
  
  # Columns of the .xls file to add in .txt:
  # All columns other than 'phys'
  
  # Column of the .txt to read only: 'phys' (do not rewrite - results of computations/BLAST)
  # if 'phys' empty, compute position of the marker
  #     else (do not modify)
  
  # Merge the two dataframes
  
  # Save in '/data'. user must copy to '/data-cleaned' for a read-only use in further analyses
  write.table(mkr, file = paste(data, set, ".txt", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
} 

#----------------------------------------------------------------------------
# Recompute the physical position of a subset of markers (error management, updates)
markerDB.redo = function() {
  
}


#============================================================================
# Summary of the marker DB
#============================================================================
# Give some summary informations about number of species, list of species, number of markers... in the DB
markerDB.infos = function() {
  
}





############### 2 procedures exist

#-------------------------------------
#     1
# Get physical positions with a list of markers names
# The list of markers names corresponds to the first column 'Marker' in the genetic map file
# And will be the 'phys' (#4) in the Marey map input file

# If the column 'mkr_data' is TRUE in the 'Genetic_maps_ressources.xlxs' database:
# If the file exists in the 'Genetic_maps' directory but does not exist yet in the 'Marey_maps' directory, find the physical position of markers in genetic maps
# and consolidate a Marey map

# First test with 'Brassica_oleracea_Wang2012.txt'
#gen_map = read.table("data/Genetic_maps/Brassica_oleracea_Wang2012.txt", header = TRUE)
# Get the reference genome on which map the markers
#phys_positions(gen_map$Marker)

#MareyMap_input()

#rm(gen_map)


#-------------------------------------
#     2
# Get physical positions with a list of primer sequences
# If markers names are not provided or not assessed in public databases
# BLASTN on a reference genome with primer sequences in a FASTA file
# Genomes in a local database (formatdb)
# Get start positions for each sequence (marker)

# If the column 'primer_data' is TRUE in the 'Genetic_maps_ressources.xlxs' database:
# If the file exists in the 'Genetic_maps' directory but does not exist yet in the 'Marey_maps' directory, find the physical position of markers in genetic maps
# and consolidate a Marey map

# BLASTN the FASTA file of primers onto the reference genome

# Keep only best hits
# Forward primer FASTA
# Reverse primer FASTA
# Compare positions to confirme the best hit


# Once genetic and physical maps are retrieved...





#============================================================================
# Get the sequences of a marker given a list of marker names
#============================================================================
# Then, sequences will be blasted against the reference genome
# Use Entrez API utilities to gete the sequence
# species = "Nelumbo_nucifera"
# mkrDB = read.csv(paste("data/Physical_maps/Markers/", species, ".csv", sep = ""), sep = ";")

mkr.names = function(species = "") {
  # Load the file
  mkrDB = read.csv(paste("data/Physical_maps/Markers/", species, ".csv", sep = ""), sep = ";")
  # If marker name but no sequence or primer associated, process...
  idx = which(is.na(mkrDB$forward_primer) & is.na(mkrDB$reverse_primer) & is.na(mkrDB$sequence_tag))
  cat(length(idx), "markers without nuleotide sequence. Checking NCBI...")
  for (i in idx) {
    mkrDB$sequence_tag[i]
  }
}

#============================================================================
# Get the physical positions of a marker given its position on a scaffold
#============================================================================
# 1/ Retrieve the scaffold if neccessary, to blast for scaffold position
# 2/ Calculate the SNP position from Scaffold start + SNP position on scaffold
# If update = FALSE, then only markers whithout a physical position are processed

# species = "Nelumbo_nucifera"
# mkrDB = read.csv(paste("data/Physical_maps/Markers/", species, ".csv", sep = ""), sep = ";")

blast.scaffolds = function(species = "", update = FALSE, tol = 95, submitter = "tbrazier") {
  library(rBLAST)
  library(Biostrings)
  library(rentrez)

  # Load the file
  mkrDB = read.csv(paste("data/Physical_maps/Markers/", species, ".csv", sep = ""), sep = ";")
  mkrDB$sequence_tag = as.character(mkrDB$sequence_tag)
  mkrDB$scaffold = as.character(mkrDB$scaffold)
  
  if (update == FALSE) {
    # If scaffold name but no physical position associated, process...
    idx = which(is.na(mkrDB$phys) & !is.na(mkrDB$scaffold))
    cat(length(idx), "scaffolds without assigned position. Checking NCBI...\n")
  } else {
    # If scaffold name but no physical position associated, process...
    idx = which(!is.na(mkrDB$scaffold))
    cat(length(idx), "scaffolds to assign. All scaffolds will be updated. Checking NCBI...\n")
  }
  
  # Retrieve scaffold sequences
  # Long process, since NCBI is called for each marker
  pb = txtProgressBar(min = 0, max = length(idx), style = 3)
  for (i in idx) {
    setTxtProgressBar(pb, which(idx == i))
    # Retrieve the fasta from NCBI
    id = NA
    getID = function(x) {
      return(entrez_search(mkrDB$scaffold[x], db ="nuccore")$id)
    }
    id = unlist(try(getID(i), silent = TRUE))
    # Prevent scaffolds not retunring results with NA
    if (length(id) > 0) {
      if (!is.na(id)) {
        fasta = entrez_fetch(db = "nuccore", id, rettype = "fasta")
        # Remove fasta header and end of lines
        # Fasta header is first line toward the first \n encountered
        fasta = strsplit(fasta, "\n")
        fasta = paste(fasta[[1]][-1], collapse = "")
        mkrDB$sequence_tag[i] = as.character(fasta)
        rm(fasta)
      } else {
        mkrDB$sequence_tag[i] = NA
      }
    } else {
      mkrDB$sequence_tag[i] = NA
    }
  
  }
  Sys.sleep(1)
  close(pb)
  
  # Subset scaffold sequences
  # Take 150bp around the SNP position on the scaffold to reduce computational load
  for (i in idx) {
    mkrDB$sequence_tag[i] = substr(mkrDB$sequence_tag[i], start = mkrDB$SNPonScaffold[i] - 75,  stop = mkrDB$SNPonScaffold[i] + 75) 
  }
  
  # Blast to get scaffold position
  cat("Blast to get scaffold position\n")
  mkrDB$accession = as.character(mkrDB$accession)
  accession = unique(mkrDB$accession[idx])[1]
  mkrDB$submitter = as.character(mkrDB$submitter)
  mkrDB$date = as.character(mkrDB$date)
  # Run blast on a local db in a repertory 'species/accession/fasta/'
  ## load a BLAST database (replace db with the location + name of the BLAST DB)
  ls = list.files(path = paste(wd, "/data/Genome/", tolower(species),"/" ,accession,"/fasta/", sep = ""))
  dbname = gsub(".nin", "", ls[grepl(".nin", ls)])
  #dbname = gsub(".gz", "",list.files(path = paste(wd, "/data/Genome/", tolower(species),"/" ,accession,"/fasta/", sep = ""))) # Older way to find the db name, not working with Ensembl db
  blastDB = blast(db = paste(wd, "/data/Genome/", tolower(species),"/" ,accession,"/fasta/", dbname, sep = ""), type = "blastn")
  blastDB
  # For each marker independently
  pb = txtProgressBar(min = 0, max = length(idx), style = 3)
  for (i in idx) {
    map = NA
    # cat(i)
    setTxtProgressBar(pb, which(idx == i))
    tag = as.character(mkrDB$sequence_tag[i])
    tag = gsub("[A-Z]/[A-Z]", "N",tag) # Replacing any substitution site ([A/C type]) by a N nucleotide
    tag = gsub("_", "N",tag) # Replacing any substitution site (_ type) by a N nucleotide
    tag = gsub("\\[", "",tag)
    tag = gsub("\\]", "",tag)
    tag = gsub("/", "",tag)
    if (is.na(tag)) {
      map = NA
      phys = NA
      start = NA
      end = NA
    } else {
      ## query a sequence using BLAST
      seq = DNAStringSet(x=tag, use.names=TRUE) # Convert the sequence tag as a XStringSet
      # Short sequences need an optional task to adjust the expected value
      if (seq@ranges@width < 30) {
        hit = predict(blastDB, seq, BLAST_args = "-task blastn-short") # Search for the sequence tag position in the genome DB, with adjustement for short sequences
      } else {
        hit = predict(blastDB, seq) # Search for the sequence tag position in the genome DB
      }      
      #hit[1:5,]
      # print(hit, infor = TRUE)
      # Keep only the first hit, with more than 95% of identity (it is sequence tags, so we expect a perfect match)
      query = hit[1,]
      query = query[which(query$Perc.Ident > tol),]
      
      if (nrow(query) == 1) {
        id = as.character(query$SubjectID[1])
        # In case of Ensembl sequences, the naming standard is different. Subject id is already the chromosome number
        # If 'id' is '[0-9]+', (1) then it is the chr number of an Ensembl fasta sequence
        # Else, (2) if 'id' contains scaffold, it is also an Ensembl fasta but we're not interested in the position on the scaffold (irrelevant for a physcial map)
        # Otherwise, it is not an Ensembl fasta (3)
        id = gsub("lg", "", id)
        if (grepl("^[0-9]+$", id)) { # (1)
          map = id
        } else {
          if (grepl("scaffold", id)) { # (2)
            # Do nothing with scaffolds for the moment
          } else { # (3)
            # Recover in the BlastDB the Chromosome number from the SubjectID
            res = system(paste("blastdbcmd -db ", wd, "/data/Genome/", tolower(species), "/" , accession,"/fasta/", dbname, " -entry ", id, sep = ""), intern = TRUE)
            res[1]
            map = gsub("chromosome:*[ ]*[A-Za-z]*[[:punct:]]*", "", regmatches(res[1], regexpr("chromosome:*[ ]*[A-Za-z0-9]*[[:punct:]]*[A-Za-z]*[0-9]+[A-Za-z]*", res[1])))
          }
        }
        # Find the position of the sequence on the blasted chromosome
        if (length(map) == 0) { # If no chromosome number has been found
          map = NA
          phys = NA
          start = NA
          end = NA
        } else {
          # get positions
          phys = query$S.start[1] + mkrDB$SNPonScaffold[i]
          start = query$S.start[1]
          end = query$S.end[1]
        }
      } else {
        map = NA
        phys = NA
        start = NA
        end = NA
      }
    }
    
    # Saving results in the marker database
    mkrDB$map[i] = map
    # Physical position is scaffold position + SNP position on scaffold
    mkrDB$phys[i] = phys
    mkrDB$start[i] = start
    mkrDB$end[i] = end
    # Record submitter for the blast position and date
    mkrDB$submitter[i] = submitter
    mkrDB$date[i] = as.character(gsub(" ", "_", date()))
  } # END OF LOOP ON MARKERS
  Sys.sleep(1)
  close(pb)  

  # Save actions
  save.log()
  save.log(msg = paste(sum(!is.na(mkrDB$phys)), " physical positions retrieved and ", sum(is.na(mkrDB$phys))," NA." , sep = ""))
  
  return(mkrDB)

}

# esearch -db nucleotide -query "scaffold00123 AND Pyrus" | eftech -format fasta

#============================================================================
# Get the physical positions of a marker given a list of sequence tag or primers
#============================================================================
# Given a sequence tag or primers of a marker, the function retrieve the position of the marker
# An accession to a reference genome can be provided to narrow the search
# Given a marker database, the function retrieve the position of the markers
# An accession to a reference genome must be provided to narrow the search
# Return a marker database updated with positions
# tol: a tolerance parameter, the minimum percentage of identity with the reference genome accepted to accept a marker position
# If update.all is false, compute only for marker positions that are NA, otherwise, compute positions of all marker, hence erasing older positions (FALSE by default)
blast.markers = function(species = "", accession = "", db = mkrDB, tol = 95, submitter = "tbrazier", update.all = FALSE) {
  library(rBLAST)
  library(Biostrings)
  mkrDB$map = as.character(mkrDB$map)
  mkrDB$species = as.character(mkrDB$species)
  mkrDB$accession = as.character(mkrDB$accession)
  # Cleaning sequence tags in mkrDB, removing unwanted characters
  mkrDB$sequence_tag = gsub(" ", "", mkrDB$sequence_tag)
  mkrDB$sequence_tag = gsub("[*]", "", mkrDB$sequence_tag)
  mkrDB$sequence_tag = gsub("-", "", mkrDB$sequence_tag)
  mkrDB$sequence_tag = gsub("–", "", mkrDB$sequence_tag)
  mkrDB$sequence_tag = gsub("_", "", mkrDB$sequence_tag)
  mkrDB$sequence_tag = gsub("\t", "", mkrDB$sequence_tag)
  mkrDB$sequence_tag = gsub("\n", "", mkrDB$sequence_tag)
  mkrDB$sequence_tag = gsub("#", "", mkrDB$sequence_tag)
  mkrDB$sequence_tag = gsub("#", "", mkrDB$sequence_tag)
  mkrDB$sequence_tag = gsub("[Xx]", "N", mkrDB$sequence_tag)
  # mkrDB$sequence_tag = gsub("/", "", mkrDB$sequence_tag)
  # Run blast on a local db in a repertory 'species/accession/fasta/'
  ## load a BLAST database (replace db with the location + name of the BLAST DB)
  ls = list.files(path = paste(wd, "/data/Genome/", tolower(species),"/" ,accession,"/fasta/", sep = ""))
  dbname = gsub(".nin", "", ls[grepl(".nin", ls)])
  #dbname = gsub(".gz", "",list.files(path = paste(wd, "/data/Genome/", tolower(species),"/" ,accession,"/fasta/", sep = ""))) # Older way to find the db name, not working with Ensembl db
  blastDB = blast(db = paste(wd, "/data/Genome/", tolower(species),"/" ,accession,"/fasta/", dbname, sep = ""), type = "blastn")
  blastDB
  # For each marker independently
  pb = txtProgressBar(min = 0, max = nrow(mkrDB), style = 3)
  for (i in 1:nrow(mkrDB)) {
    map = NA
    setTxtProgressBar(pb, i)
    accession = mkrDB$accession[i]
    # Check first if the marker position need to be computed/updated
    if (is.na(mkrDB$phys[i]) | (!(is.na(mkrDB$phys[i])) & update.all == TRUE)) {
      if (!is.na(mkrDB$sequence_tag[i])) { # If your marker is identified by a sequence tag
        tag = as.character(mkrDB$sequence_tag[i])
        tag = gsub("[A-Z]/[A-Z]", "N",tag) # Replacing any substitution site ([A/C type]) by a N nucleotide
        tag = gsub("_", "N",tag) # Replacing any substitution site (_ type) by a N nucleotide
        tag = gsub("\\[", "",tag)
        tag = gsub("\\]", "",tag)
        tag = gsub("/", "",tag)
        ## query a sequence using BLAST
        seq = DNAStringSet(x=tag, use.names=TRUE) # Convert the sequence tag as a XStringSet
        # Short sequences need an optional task to adjust the expected value
        if (seq@ranges@width < 30) {
          hit = predict(blastDB, seq, BLAST_args = "-task blastn-short") # Search for the sequence tag position in the genome DB, with adjustement for short sequences
        } else {
          hit = predict(blastDB, seq) # Search for the sequence tag position in the genome DB
        }      
        #hit[1:5,]
        # print(hit, infor = TRUE)
        # Keep only the first hit, with more than 95% of identity (it is sequence tags, so we expect a perfect match)
        query = hit[1,]
        query = query[which(query$Perc.Ident > tol),]
        
        if (nrow(query) == 1) {
          id = as.character(query$SubjectID[1])
          # In case of Ensembl sequences, the naming standard is different. Subject id is already the chromosome number
          # If 'id' is '[0-9]+', (1) then it is the chr number of an Ensembl fasta sequence
          # Else, (2) if 'id' contains scaffold, it is also an Ensembl fasta but we're not interested in the position on the scaffold (irrelevant for a physcial map)
          # Otherwise, it is not an Ensembl fasta (3)
          id = gsub("lg", "", id)
          id = gsub("Chr", "", id)
          id = gsub("^chr[0-9]+LG", "", id)
          if (grepl("^[0-9]+$", id)) { # (1)
            map = id
          } else {
            if (grepl("scaffold", id)) { # (2)
              # Do nothing with scaffolds for the moment
            } else { # (3)
              # Recover in the BlastDB the Chromosome number from the SubjectID
              res = system(paste("blastdbcmd -db ", wd, "/data/Genome/", tolower(species), "/" , accession,"/fasta/", dbname, " -entry ", id, sep = ""), intern = TRUE)
              res[1]
              map = gsub("chromosome:*[ ]*[A-Za-z]*[[:punct:]]*", "", regmatches(res[1], regexpr("chromosome:*[ ]*[A-Za-z0-9]*[[:punct:]]*[A-Za-z]*[0-9]+[A-Za-z]*", res[1])))
            }
          }
          # Find the position of the sequence on the blasted chromosome
          if (length(map) == 0) { # If no chromosome number has been found
            map = NA
            phys = NA
            start = NA
            end = NA
          } else {
            # get positions
            phys = query$S.start[1]
            start = query$S.start[1]
            end = query$S.end[1]
          }
        } else {
          map = NA
          phys = NA
          start = NA
          end = NA
        }
        # Saving results in the marker database
        mkrDB$map[i] = map
        mkrDB$phys[i] = phys
        mkrDB$start[i] = start
        mkrDB$end[i] = end
        # Record submitter for the blast position and date
        mkrDB$submitter[i] = submitter
        mkrDB$date[i] = as.character(gsub(" ", "_", date()))
        
      } else {
        if (!is.na(mkrDB$forward_primer[i]) & !is.na(mkrDB$reverse_primer[i])) { # Otherwise, if your marker is identified by two primers
          # A good position is one with a minimal distance between the forward and reverse primers, delimitating short fragments
          # Hence, blast the forward an reverse, and check within the best hits which pair of positions have a minimal distance
          forward = as.character(mkrDB$forward_primer[i])
          reverse = as.character(mkrDB$reverse_primer[i])
          
          # Remove unattended characters
          forward = gsub(" *", "", forward)
          reverse = gsub(" *", "", reverse)
          
          ## query a sequence using BLAST
          seq.forward = DNAStringSet(x = forward, use.names = TRUE) # Convert the sequence tag as a XStringSet
          seq.reverse = DNAStringSet(x = reverse, use.names = TRUE) # Convert the sequence tag as a XStringSet
          
          # Search for short sequences need little adjustments of expected value
          # http://www.ncbi.nlm.nih.gov/BLAST/Why.shtml
          # You can adjust both the word size and the expect value on the standard BLAST pages to work with short sequences.
          # However, we do provide a BLAST page with these values preset to give optimum results with short sequences.
          # This page ("Search for short and nearly exact matches") is linked under the nucleotide BLAST section of the main BLAST page. The adjustments are described in the table below.
          # 
          # Program                             Word Size   Filter Setting  Expect Value
          # ------------------------------------------------------------------------------
          #   Standard Nucleotide BLAST              11         On (DUST)           10
          # Search for short/near exact matches     7         Off               1000
          
          
          hitf = predict(blastDB, seq.forward, BLAST_args = "-task blastn-short") # Search for the primer position in the genome DB
          hitr = predict(blastDB, seq.reverse, BLAST_args = "-task blastn-short") # Search for the primer position in the genome DB
          
          # Keep only hits with more than 99% (default value) of identity: primers must be perfect hits
          hitf = hitf[which(hitf$Perc.Ident > 99),]
          hitr = hitr[which(hitr$Perc.Ident > 99),]
          
          # Continue if both primers have at least one hit
          if (nrow(hitf) > 0 & nrow(hitr) > 0) {
            # Find the most probable dyad of positions
            # Primers are often positions separated by a short physical distance (identification of beginning and end of a short DNA fragment)
            # Too many hits, restrict to positions in the same chromosome as linkage map provided (assume correct linkage map)
            # Find the pair of positions for which positions are on the same chromosome as the linkage map
            # and dist(forward  - reverse) > 0 (conservative - true only if genome oriented in the same way as markers)
            # and dist(forward  - reverse) < 1000bp (targeting short fragments)
            
            # Chromosome number in the linkage map
            chr.linkage = as.character(mkrDB$map[i])
            # get a list of chromosome IDs in DB
            list.chr = c()
            IDs = unique(as.character(hitf$SubjectID))
            
            for (c in 1:length(IDs)) {
              # print(IDs[c])
              # In case of Ensembl sequences, the naming standard is different. Subject id is already the chromosome number
              # If 'id' is '[0-9]+', (1) then it is the chr number of an Ensembl fasta sequence
              # Else, (2) if 'id' contains scaffold, it is also an Ensembl fasta but we're not interested in the position on the scaffold (irrelevant for a physcial map)
              # Otherwise, it is not an Ensembl fasta (3)
              if (grepl("^[A-Z]*[0-9]+[A-Z]*$", IDs[c])) { # (1)
                map = IDs[c]
              } else {
                if (grepl("scaffold", IDs[c])) { # (2)
                  # Do nothing with scaffolds for the moment
                } else { # (3)
                  # Recover in the BlastDB the Chromosome number from the SubjectID
                  res = system(paste("blastdbcmd -db ", wd, "/data/Genome/", tolower(species), "/" , accession,"/fasta/", dbname, " -entry ", IDs[c], sep = ""), intern = TRUE)
                  res[1]
                  map = gsub("chromosome:*[ ]*[A-Za-z]*[[:punct:]]*", "", regmatches(res[1], regexpr("chromosome:*[ ]*[A-Za-z0-9]*[[:punct:]]*[A-Za-z]*[0-9]+[A-Za-z]*", res[1])))
                  # map = gsub("chromosome[ ]*[A-Za-z]*", "",regmatches(res[1], regexpr("chromosome[ ]*[A-Za-z]*.*[0-9]+", res[1])))
                  # map = gsub("[.]", "", map)
                }
              }
              list.chr = c(list.chr, map)
            }
            rm(map)

            # Retain IDs with the same chromosome number as in linkage map
            if (!is.na(chr.linkage)) {
              # Which ID is those of the chromosome in linkage map
              # The only SubjectID to retain is...
              mask = IDs[which(as.character(list.chr) == as.character(chr.linkage))]
            } else{ # But if no linkage map provided...
              # Else, first ID, those with the higher E value
              mask = IDs[1]
            }
            # Besides, if no corresponding chromosome names
            if (is.null(mask)) {
              mask = IDs[1]
            }

            
            # IDs are line number in the query results
            # List of all possible dyads (i.e. same chromosome)
            # Stored in a matrix
            res = matrix(NA, nrow = 0, ncol = 3)
            # Which rows for SubjectID?
            hitf = hitf[hitf$SubjectID == mask,]
            
            for (q in 1:nrow(hitf)) {
              # For each forward primer, associate all hits on the same chromosome and compute the physical distance
              # Count reverse hits...
              chr = hitf$SubjectID[q]
              match = hitr[as.character(hitr$SubjectID) == chr,]
              # If same chromosome as linkage and matching with reverse
              if (nrow(match) > 0) {
                for (m in 1:nrow(match)) {
                  # Add a new line
                  res = rbind(res, c(row.names(hitf[q,]), row.names(match[m,]), as.numeric(hitf$S.end[q] - match$S.start[m])))
                }
              }
              rm(chr)
              rm(match)
            }
            res = as.data.frame(res)
            res[,3] = as.numeric(levels(res[,3]))[res[,3]]
            res = res[!is.na(res[,3]),]
            # Keep only positive distances, assuming same orientation in genome and markers
            # sel = res[as.numeric(res[,3]) > 0,]
            # Transform all values in absolute distances
            sel = res
            sel[,3] = abs(sel[,3])
            # Keep fragment < 200
            sel = sel[as.numeric(sel[,3]) < 200,]
            # Select the shortest fragment
            if (nrow(sel) > 1) {
              sel = sel[which(sel[,3] == min(sel[,3])),]
            }
            sel
            
            if (nrow(sel) > 0) {
              # Get chromosome name from the subject ID
              id = as.character(hitf$SubjectID[1])
              # In case of Ensembl sequences, the naming standard is different. Subject id is already the chromosome number
              # If 'id' is '[0-9]+', (1) then it is the chr number of an Ensembl fasta sequence
              # Else, (2) if 'id' contains scaffold, it is also an Ensembl fasta but we're not interested in the position on the scaffold (irrelevant for a physcial map)
              # Otherwise, it is not an Ensembl fasta (3)
              if (grepl("^[A-Z]*[0-9]+[A-Z]*$", id)) { # (1)
                map = id
              } else {
                if (grepl("scaffold", id)) { # (2)
                  # Do nothing with scaffolds for the moment
                } else { # (3)
                  # Recover in the BlastDB the Chromosome number from the SubjectID
                  res = system(paste("blastdbcmd -db ", wd, "/data/Genome/", tolower(species), "/" , accession,"/fasta/", dbname, " -entry ", id, sep = ""), intern = TRUE)
                  res[1]
                  map = gsub("chromosome:*[ ]*[A-Za-z]*[[:punct:]]*", "", regmatches(res[1], regexpr("chromosome:*[ ]*[A-Za-z0-9]*[[:punct:]]*[A-Za-z]*[0-9]+[A-Za-z]*", res[1])))
                  # map = gsub("chromosome[ ]*[A-Za-z]*", "",regmatches(res[1], regexpr("chromosome[ ]*[A-Za-z]*.*[0-9]+", res[1])))
                  # map = gsub("[.]", "", map)
                }
              }
              
              # get positions
              map = hitf[which(row.names(hitf) == as.numeric(levels(droplevels(sel[1,1])))),]$SubjectID
              phys = hitf[which(row.names(hitf) == as.numeric(levels(droplevels(sel[1,1])))),]$S.start
              start = hitf[which(row.names(hitf) == as.numeric(levels(droplevels(sel[1,1])))),]$S.start
              end = hitr[which(row.names(hitr) == as.numeric(levels(droplevels(sel[1,2])))),]$S.end
            } else {
              map = NA
              phys = NA
              start = NA
              end = NA
            }
            
            # start = hitf[which(row.names(hitf) == as.numeric(levels(sel[1,1]))[sel[1,1]]),]$S.start
            # end = hitr[which(row.names(hitr) == as.numeric(levels(sel[1,2]))[sel[1,2]]),]$S.end
            # Saving results in the marker database
            mkrDB$map[i] = map
            mkrDB$phys[i] = phys
            mkrDB$start[i] = start
            mkrDB$end[i] = end
            # Record submitter for the blast position and date
            mkrDB$submitter[i] = submitter
            mkrDB$date[i] = as.character(gsub(" ", "_", date()))
          }
          
        }
      }
    }
    # end of the loop for each marker
  }
  
  mkrDB$map = as.character(mkrDB$map)
  Sys.sleep(1)
  close(pb)
  # Save actions
  save.log()
  save.log(msg = paste(sum(!is.na(mkrDB$phys)), " physical positions retrieved and ", sum(is.na(mkrDB$phys))," NA." , sep = ""))
  
  return(mkrDB)
}




#============================================================================
# Summary statistics of a map
#============================================================================




