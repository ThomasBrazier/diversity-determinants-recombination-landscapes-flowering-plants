#==========================================================
# Import & compile bibliographic database (DB)
#==========================================================
# data = path and filename of the xlsx file with references
# output = path and filename where to save the DB

genomeDB.compile = function(data = "data/Genome/Genome_ressources.xlsx", output = "data/Genome/GenomeDB.txt", forced = FALSE) {
  require(gdata)
  require(taxize)
  require(taxizedb)
  require(dplyr)
  require(rentrez)
  # Import the handmade bibliographic database in 'data'
  bib = read.xls (data, sheet = 1, header = TRUE)
  
  # We work at a chromosome level (full chromosome sequence); remove scaffolds and contigs
  bib = bib[bib$assembly_level == "chromosome",]
  # Clean the 'species_complete' names in a reduced name 'Genus_sp' without variety and subspecies names
  bib$species = as.character(regmatches(as.character(bib$species_complete), gregexpr("^[A-Z][a-z]+_[a-z]+", as.character(bib$species_complete))))
  # Particular case of crossings "*_x_*" -> get full name of the crossing
  bib$species[grepl("_x_", bib$species_complete)] = as.character(bib$species_complete[grepl("_x_", bib$species_complete)])
  # Remove crossings from dataset
  bib = bib[-(grepl("_x_", bib$species_complete)),]
  
  # 'id' is 'Genus_sp_accession'
  bib$id = paste(bib$species, bib$accession, sep = "_")
  id = gsub("_", " ", bib$species)
  
  #---------------------------------------------------
  # Get 'taxid' from species names when not completed
  for (i in 1:length(bib$species)) {
    print(i)
    # Unique ID in NCBI
    if (is.na(bib$species_taxid[i])) {
      uid = get_ids(names = bib$species[i], db = c('ncbi'))
      if (!is.na(uid$ncbi)) {
        bib$species_taxid[i] =  as.character(uid$ncbi)
        
        # print("taxon calling")
        # Classify plants in 'Bryophyta', 'Fern', 'Gymnosperm' or 'Angiosperm' in 'taxon'
        # hierarchial classification
        class = classification.local(as.character(bib$species_taxid[i]), db = 'ncbi', callopts = "&api_key=5f9324d393fefb5b55946444205a0ea24c08") # callopts = "api_key=5f9324d393fefb5b55946444205a0ea24c08"
        class = paste(class[[1]]$name, collapse = ";") # Format hierarchy as a character string
        # Search hierarchically which group of plant it is
        if (is.na(bib$taxon[i]) | !is.na(bib$species_taxid[i])) {
          if (grepl("Magnoliopsida", class) | grepl("Liliopsida", class)) {bib$taxon[i] = "angiosperm"} else {
            if (grepl("Acrogymnospermae", class)) {bib$taxon[i] = "gymnosperm"} else {
              if (grepl("Polypodiopsida", class)) {bib$taxon[i] = "fern"} else {
                if (grepl("Lycopodiopsida", class)) {bib$taxon[i] = "moss"}
              }
            }
          }
        }
        }
      }
  }
  bib$species_taxid = as.character(bib$species_taxid)
  
  #---------------------------------------------------
  bib$release_date = as.character(bib$release_date)
  bib$assembly_description = as.character(bib$assembly_description)
  bib$assembly_name = as.character(bib$assembly_name)
  bib$submitter = as.character(bib$submitter)
  # Retrieve metadata
  cat("Retrieve metadata...\n")
  for (i in 1:nrow(bib)) {
    print(i)
    # Retrieve metadata only if reference not validated; Update ref if forced = TRUE
    if (bib$vld[i]== FALSE | forced == TRUE) {
      # Get release dates from the accession
      id_search = entrez_search(db="assembly",
                                term=bib$accession[i],
                                retmax=1)
      if (length(id_search$ids) > 0) { # In case of empty list (no result from NCBI search)
        id_summ = entrez_summary(db="assembly",
                                 id=id_search$ids)
        bib$release_date[i] = as.character(gsub("/" , "-",gsub(pattern = " [0-9]+:[0-9]+", replacement = "", id_summ$asmreleasedate_genbank)))
        # Get submitter name or organization
        bib$submitter[i] = shQuote(id_summ$submitterorganization)
        # Get assembly name
        bib$assembly_name[i] = shQuote(id_summ$assemblyname)
        # Get assembly description from the accession, standardized from NCBI, if NA
        bib$assembly_description[i] = shQuote(id_summ$assemblydescription)
        rm(id_search)
        rm(id_summ)}
      
      # Bibliographic metadata: ref, year, title, DOI
      
      # Validate references fetched on ncbi
      bib$vld[i] = TRUE
    }
  }

  #---------------------------------------------------
  # Manual assessment of reference genomes that need to be dl in the dataset
  # Remove duplicates and selection of reference genomes for analyses
  bib$vld[(bib$taxon != "angiosperm") | is.na(bib$taxon)] = FALSE
  # Hierarchy of selection:
  # Keep first Ensembl, NCBI and Phytozome, with preference to the latest accession,
  # Then more specialised databases...
  
  
  
  #---------------------------------------------------
  # Reformating columns
  bib$assembly_level = as.character(bib$assembly_level)
  # Save 'title' as a character string with quotes, because of spaces
  bib$title = paste("\"",as.character(bib$title),"\"", sep ="")
  # Save the DB into a clean DB .txt file in 'data'
  write.table(bib, file = output, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}



#==========================================================
# Summary statistics
#==========================================================
# Number of genomes
# Number of unique species
# Number of genomes at chromosome level
# Number of unique species at chromosome level
# List of unique species at chromosome level

genomeDB.summary = function(data = "data/Genome/GenomeDB.txt") {
  # Import database
  bib = read.table(data, header = TRUE, sep = "\t")
  
  # Phylogeny
  cat("PHYLOGENY\n")
  print(ftable(bib$taxon))
  cat("\n")
  
  cat("ANGIOSPERMS ONLY\n\n")
  bib = bib[bib$taxon == "angiosperm" & bib$vld == TRUE,]
  # Number of genomes
  cat("Number of unique genomes:", length(unique(bib$accession)), "\n")
  # Number of unique species
  cat("Number of unique species:", length(unique(bib$species)), "\n")
  # Number of genomes at chromosome level
  # cat("Number of genomes at a chromosome level:", sum(as.character(bib$assembly_level) == "Chromosome"), "\n")
  # cat("Number of genomes at a chromosome level:", as.character(table(genDB$assembly_level)[[1]]), "\n")
  # Number of unique species at chromosome level
  # cat("Number of unique species with a genome at a chromosome level:", length(unique(bib$species[which(bib$assembly_level == "Chromosome")])), "\n")
  # Source of genomes at a chromosome level
  cat("\n")
  cat("Source of genomes\n")
  print(ftable(bib$source[which(bib$assembly_level == "chromosome")]))
  # List of unique species at chromosome level
  cat("\n")
  cat("List of species\n")
  print(sort(unique(bib$species[which(bib$assembly_level == "chromosome")])))
  
  # Number of accessions
  cat("\n")
  cat("Number of accessions\n")
  print(length(unique(bib$accession)))
  # Number of accessions
  cat("List of accessions\n")
  print(sort(unique(bib$accession)))
  
  # Number of valid genomes integrated in analyses
  cat("\n")
  cat("Number of valid genomes integrated in analyses:", sum(bib$vld, na.rm = TRUE), "\n")
  # List of valid genomes integrated in analyses
  cat("\n")
  cat("List of valid genomes integrated in analyses\n")
  print(sort(bib$id[bib$vld]))
}


#==========================================================
# Download genomes
#==========================================================

#----------------------------------------------------------
# List of genomes to download, by accessions
#----------------------------------------------------------
# Make a list of genomes to download in files for Ensembl, NCBI, Phytozome and custom requests to specific databases
mklist_genomes_dl = function() {
  
}


# #----------------------------------------------------------
# # NCBI
# #----------------------------------------------------------
# # accn = a list of accessions for genomes to download
# # output = path to the output directory
# genome_dl_ncbi = function(accn = c(), output = "/data/Genome/") {
#   for (g in accn) {
#     print(g)
#     # Get Fasta files from NCBI with a list of accession
#     cat("Fasta...\n")
#     system(paste("/Users/tbrazier/miniconda3/bin/esearch -db assembly -query ", as.character(g)," | /Users/tbrazier/miniconda3/bin/elink -target nucleotide -name assembly_nuccore_insdc | /Users/tbrazier/miniconda3/bin/efetch -format fasta > ", wd, output,"fasta/", as.character(g),".fa", sep = ""),
#            intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
#     system((paste("/Users/tbrazier/miniconda3/bin/esearch -db assembly -query ", g, sep = "")),
#            intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
#     # Get Genbank files from NCBI with a list of accession
#     cat("Genbank...\n")
#     # system(paste("/Users/tbrazier/miniconda3/bin/esearch -db assembly -query ", g," | /Users/tbrazier/miniconda3/bin/elink -target nucleotide -name assembly_nuccore_insdc | /Users/tbrazier/miniconda3/bin/efetch -format gb > ", wd, output,"genbank/", g,".gbk", sep = ""),
#     #        intern = FALSE)
#   }
# }


#==========================================================
# NCBI Taxonomy browser
#==========================================================
# https://github.com/ropensci/taxize
# https://github.com/ropensci/taxizedb
# install.packages("taxize")
# install.packages("taxizedb")
require("taxize")
# library("taxizedb")
source("sources/taxizedb_custom.R")
require("dplyr")
# Download NCBI database in a local cached version
# db_download_ncbi()



