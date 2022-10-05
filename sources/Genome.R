#----------------------------------------------------------------------------
# Chromosome size
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
# Retrieve feature positions in a GFF file
#----------------------------------------------------------------------------
# 
# Retrieve positions  of a specified feature in a chromosome sequence coded in .gff
# Return a data frame easy to use in R,
# with feature id, feature type, chromosome, start, end, type (e.g. family of transposon, +/- for genes)
# Execute in the directory of your genomes
# 
# ARGUMENTS
#   species = species names in format 'Genus_species'
#   accession = accession number of the genome
#   feature = feature type to retrieve in c("gene", "transposon", "exon", "CDS")
# 
# VALUES
# Return a data frame

feature.position = function(species = "", accession = "", feature = "") {
  require(ape) # dependency to read GFF files with read.gff()
  # read.delim("gfffile.gff", header=F, comment.char = ""#)
  
  # TODO optim -> parallel
  
  # A dataset as example
  # species = "Arabidopsis_thaliana"
  # accession = "GCA_902460315.1"
  
  if (dir.exists(paste(wd, "data/Genome", tolower(species), accession, "gff3",sep = "/"))) {
    #----------------------------------------------------------------------------
    # DATA
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
    cat("Read data.\n")
    data = lapply(paths, function(x){read.gff(file = x)})
    data = do.call("rbind", data)

    #----------------------------------------------------------------------------
    # PROCESSING DATA: SEARCH FEATURES...
    cat("Processing data: search features...\n")
    # Features names can have multiple synonyms
    if (feature == "gene") {
      query = c("gene", "SO:0000704")
    } else {
      if (feature == "exon") {
        query = c("exon", "coding_exon", "coding exon", "SO:0000195")
      } else {
        if (feature == "CDS") {
          query = c("CDS", "coding_sequence", "coding sequence", "SO:0000316")
        } else {
          warning("The searched feature must be 'gene', 'exon' or CDS'.")
        }
      }
    }
    # Searching the required features
    data = data[data$type %in% query,]
    
    #----------------------------------------------------------------------------
    # Retrieve queried informations and format output
    if (feature == "gene") {
      ids = gsub(";", "", gsub("ID=gene:", "", regmatches(data$attributes, gregexpr("ID=gene:[A-Za-z0-9_]+;", data$attributes))))
    } else {
      if (feature == "exon") {
        # Exception for particular case
        # For some species, exon_id is not working
        if (accession == "GCA_004115385.1") {
          ids = gsub(";", "", regmatches(data$attributes, gregexpr("HF[A-Za-z0-9]+-[A-Za-z0-9]+-[A-Za-z0-9]+;", data$attributes)))
        } else { # General case
          ids = gsub(";", "", gsub("exon_id=", "", regmatches(data$attributes, gregexpr("exon_id=[A-Za-z0-9[:punct:]]+;", data$attributes))))
          # TODO problem here
          if (identical(unique(ids), character(0))) {
            ids = gsub(";", "", gsub("ID=", "", regmatches(data$attributes, gregexpr("ID=[A-Za-z0-9[:punct:]]+;*", data$attributes))))
          }
        }
      } else {
        if (feature == "CDS") {
          ids = as.character(gsub(";", "", gsub("ID=CDS:", "", regmatches(data$attributes, gregexpr("ID=CDS:[A-Za-z0-9\\.]+;", data$attributes)))))
        } else {
          warning("The searched feature must be 'gene', 'exon' or CDS'.")
        }
      }
    }
    # Sometimes id is not recovered properly or not indicated in gff
    ids[which(as.character(ids) == "character(0)")] = NA
    # For exons, get the rank of the exon
    if (feature == "exon") {
      # Exception for particular case
      # For some species, exon_id is not working
      if (accession == "GCA_004115385.1") {
        rank = gsub(";", "", gsub("HF[A-Za-z0-9]+-[A-Za-z0-9]+-", "", regmatches(data$attributes, gregexpr("HF[A-Za-z0-9]+-[A-Za-z0-9]+-[A-Za-z0-9]+;", data$attributes))))
      } else { # General case
        rank = as.character(gsub(";rank=", "", regmatches(data$attributes, gregexpr(";rank=[0-9]+", data$attributes))))
      }
      featurepos = data.frame(id = ids, type = rep(feature, length(ids)),
                              chromosome = data$seqid, strand = data$strand, rank = rank,
                              start = data$start, end = data$end)
    } else {
      featurepos = data.frame(id = ids, type = rep(feature, length(ids)),
                              chromosome = data$seqid, strand = data$strand, phase = data$phase,
                              start = data$start, end = data$end)  
    }
    #----------------------------------------------------------------------------
    # Take care of chromosome names - Naming scheme must be respected
    # Malus domestica, for example, have a naming like 'CM014[0-9]+.1' for chromosomes
    # Table of translation for Malus domestica for 17 chromosomes
    featurepos$chromosome = as.character(featurepos$chromosome)
    featurepos$chromosome[featurepos$chromosome == "CM014049.1"] = "1"
    featurepos$chromosome[featurepos$chromosome == "CM014050.1"] = "2"
    featurepos$chromosome[featurepos$chromosome == "CM014051.1"] = "3"
    featurepos$chromosome[featurepos$chromosome == "CM014052.1"] = "4"
    featurepos$chromosome[featurepos$chromosome == "CM014053.1"] = "5"
    featurepos$chromosome[featurepos$chromosome == "CM014054.1"] = "6"
    featurepos$chromosome[featurepos$chromosome == "CM014055.1"] = "7"
    featurepos$chromosome[featurepos$chromosome == "CM014056.1"] = "8"
    featurepos$chromosome[featurepos$chromosome == "CM014057.1"] = "9"
    featurepos$chromosome[featurepos$chromosome == "CM014058.1"] = "10"
    featurepos$chromosome[featurepos$chromosome == "CM014059.1"] = "11"
    featurepos$chromosome[featurepos$chromosome == "CM014060.1"] = "12"
    featurepos$chromosome[featurepos$chromosome == "CM014061.1"] = "13"
    featurepos$chromosome[featurepos$chromosome == "CM014062.1"] = "14"
    featurepos$chromosome[featurepos$chromosome == "CM014063.1"] = "15"
    featurepos$chromosome[featurepos$chromosome == "CM014064.1"] = "16"
    featurepos$chromosome[featurepos$chromosome == "CM014065.1"] = "17"
    
    # For Theobroma cacao, chromosomes are in latin numerals
    # Convert to arabic numbers
    featurepos$chromosome[featurepos$chromosome == "I"] = "1"
    featurepos$chromosome[featurepos$chromosome == "II"] = "2"
    featurepos$chromosome[featurepos$chromosome == "III"] = "3"
    featurepos$chromosome[featurepos$chromosome == "IV"] = "4"
    featurepos$chromosome[featurepos$chromosome == "V"] = "5"
    featurepos$chromosome[featurepos$chromosome == "VI"] = "6"
    featurepos$chromosome[featurepos$chromosome == "VII"] = "7"
    featurepos$chromosome[featurepos$chromosome == "VIII"] = "8"
    featurepos$chromosome[featurepos$chromosome == "IX"] = "9"
    featurepos$chromosome[featurepos$chromosome == "X"] = "10"
    featurepos$chromosome[featurepos$chromosome == "XI"] = "11"
    featurepos$chromosome[featurepos$chromosome == "XII"] = "12"
    featurepos$chromosome[featurepos$chromosome == "XIII"] = "13"
    featurepos$chromosome[featurepos$chromosome == "XIV"] = "14"
    featurepos$chromosome[featurepos$chromosome == "XV"] = "15"
    featurepos$chromosome[featurepos$chromosome == "XVI"] = "16"
    featurepos$chromosome[featurepos$chromosome == "XVII"] = "17"
    featurepos$chromosome[featurepos$chromosome == "XVIII"] = "18"
    featurepos$chromosome[featurepos$chromosome == "XIX"] = "19"
    featurepos$chromosome[featurepos$chromosome == "XX"] = "20"
    
    # generic filtering
    featurepos$chromosome = gsub("chr0", "", featurepos$chromosome)
    featurepos$chromosome = gsub("chr", "", featurepos$chromosome)
    featurepos$chromosome = gsub("Chr0", "", featurepos$chromosome)
    featurepos$chromosome = gsub("Chr", "", featurepos$chromosome)
    featurepos$chromosome = gsub("LG0", "", featurepos$chromosome)
    featurepos$chromosome = gsub("LG", "", featurepos$chromosome)
    # Solanum tuberosum
    featurepos$chromosome = gsub("ST4.03ch[0]*", "", featurepos$chromosome)
    # Gossypium hirsutum
    featurepos$chromosome = gsub("Ghir_", "", featurepos$chromosome)
    # Cucurbita maxima
    featurepos$chromosome = gsub("Cma_", "", featurepos$chromosome)
    
    # Translate chromosome names for some accessions
    if (accession == "GCF_000309985.1") {
      featurepos$chromosome[featurepos$chromosome == "NC_024804.1"] = "A10"
      featurepos$chromosome[featurepos$chromosome == "NC_024795.1"] = "A1"
      featurepos$chromosome[featurepos$chromosome == "NC_024796.1"] = "A2"
      featurepos$chromosome[featurepos$chromosome == "NC_024797.1"] = "A3"
      featurepos$chromosome[featurepos$chromosome == "NC_024798.1"] = "A4"
      featurepos$chromosome[featurepos$chromosome == "NC_024799.1"] = "A5"
      featurepos$chromosome[featurepos$chromosome == "NC_024800.1"] = "A6"
      featurepos$chromosome[featurepos$chromosome == "NC_024801.1"] = "A7"
      featurepos$chromosome[featurepos$chromosome == "NC_024802.1"] = "A8"
      featurepos$chromosome[featurepos$chromosome == "NC_024803.1"] = "A9"
      # Translation table from ncbi
      # Chromosome A10	CM020897.1	=	NC_024804.2
      # Chromosome A01	CM020888.1	=	NC_024795.2
      # Chromosome A02	CM020889.1	=	NC_024796.2
      # Chromosome A03	CM020890.1	=	NC_024797.2
      # Chromosome A04	CM020891.1	=	NC_024798.2
      # Chromosome A05	CM020892.1	=	NC_024799.2
      # Chromosome A06	CM020893.1	=	NC_024800.2
      # Chromosome A07	CM020894.1	=	NC_024801.2
      # Chromosome A08	CM020895.1	=	NC_024802.2
      # Chromosome A09	CM020896.1	=	NC_024803.2
      # featurepos$chromosome = gsub("[A-Z]+_[0-9]+.", "", featurepos$chromosome)
    }
    # Citrus sinensis
    # Scaffold number is chromosome number
    # Only 9 chromosomes/linkage groups
    if (accession == "GCA_000695605.1") {
      featurepos$chromosome = gsub("scaffold[0]+", "", featurepos$chromosome)
      featurepos = featurepos[as.numeric(featurepos$chromosome) <= 9,]
    }
    # Capsella rubella
    # The eight first scaffolds are pseudo-chromosomes
    if (accession == "GCA_000375325.1") {
      featurepos$chromosome[featurepos$chromosome == "KB870805.1"] = "1"
      featurepos$chromosome[featurepos$chromosome == "KB870806.1"] = "2"
      featurepos$chromosome[featurepos$chromosome == "KB870807.1"] = "3"
      featurepos$chromosome[featurepos$chromosome == "KB870808.1"] = "4"
      featurepos$chromosome[featurepos$chromosome == "KB870809.1"] = "5"
      featurepos$chromosome[featurepos$chromosome == "KB870810.1"] = "6"
      featurepos$chromosome[featurepos$chromosome == "KB870811.1"] = "7"
      featurepos$chromosome[featurepos$chromosome == "KB870812.1"] = "8"
    }
    # Prunus persica
    if (accession == "GCA_000346465.2") {
      featurepos$chromosome[featurepos$chromosome == "CM007651.1"] = "1"
      featurepos$chromosome[featurepos$chromosome == "CM007652.1"] = "2"
      featurepos$chromosome[featurepos$chromosome == "CM007653.1"] = "3"
      featurepos$chromosome[featurepos$chromosome == "CM007654.1"] = "4"
      featurepos$chromosome[featurepos$chromosome == "CM007655.1"] = "5"
      featurepos$chromosome[featurepos$chromosome == "CM007656.1"] = "6"
      featurepos$chromosome[featurepos$chromosome == "CM007657.1"] = "7"
      featurepos$chromosome[featurepos$chromosome == "CM007658.1"] = "8"
    }
    # Arachis hypogaea
    if (accession == "GCA_003713155.1") {
      featurepos$chromosome[featurepos$chromosome == "CP030983.1"] = "A1"
      featurepos$chromosome[featurepos$chromosome == "CP030984.1"] = "A2"
      featurepos$chromosome[featurepos$chromosome == "CP030985.1"] = "A3"
      featurepos$chromosome[featurepos$chromosome == "CP030986.1"] = "A4"
      featurepos$chromosome[featurepos$chromosome == "CP030987.1"] = "A5"
      featurepos$chromosome[featurepos$chromosome == "CP030988.1"] = "A6"
      featurepos$chromosome[featurepos$chromosome == "CP030989.1"] = "A7"
      featurepos$chromosome[featurepos$chromosome == "CP030990.1"] = "A8"
      featurepos$chromosome[featurepos$chromosome == "CP030991.1"] = "A9"
      featurepos$chromosome[featurepos$chromosome == "CP030992.1"] = "A10"
      featurepos$chromosome[featurepos$chromosome == "CP030993.1"] = "B1"
      featurepos$chromosome[featurepos$chromosome == "CP030994.1"] = "B2"
      featurepos$chromosome[featurepos$chromosome == "CP030995.1"] = "B3"
      featurepos$chromosome[featurepos$chromosome == "CP030996.1"] = "B4"
      featurepos$chromosome[featurepos$chromosome == "CP030997.1"] = "B5"
      featurepos$chromosome[featurepos$chromosome == "CP030998.1"] = "B6"
      featurepos$chromosome[featurepos$chromosome == "CP030999.1"] = "B7"
      featurepos$chromosome[featurepos$chromosome == "CP031000.1"] = "B8"
      featurepos$chromosome[featurepos$chromosome == "CP031001.1"] = "B9"
      featurepos$chromosome[featurepos$chromosome == "CP031002.1"] = "B10"
    }
    if (accession == "GCA_002211085.2") {
      featurepos$chromosome[featurepos$chromosome == "CM008046.2"] = "1"
      featurepos$chromosome[featurepos$chromosome == "CM008047.2"] = "2"
      featurepos$chromosome[featurepos$chromosome == "CM008048.2"] = "3"
      featurepos$chromosome[featurepos$chromosome == "CM008049.2"] = "4"
      featurepos$chromosome[featurepos$chromosome == "CM008050.2"] = "5"
      featurepos$chromosome[featurepos$chromosome == "CM008051.2"] = "6"
      featurepos$chromosome[featurepos$chromosome == "CM008052.2"] = "7"
      featurepos$chromosome[featurepos$chromosome == "CM008053.2"] = "8"
      featurepos$chromosome[featurepos$chromosome == "CM008054.2"] = "9"
    }
    
    # Elaeis_guineensis
    if (accession == "GCF_000442705.1") {
      featurepos$chromosome[featurepos$chromosome == "NC_025993.1"] = "1"
      featurepos$chromosome[featurepos$chromosome == "NC_025994.1"] = "2"
      featurepos$chromosome[featurepos$chromosome == "NC_025995.1"] = "3"
      featurepos$chromosome[featurepos$chromosome == "NC_025996.1"] = "4"
      featurepos$chromosome[featurepos$chromosome == "NC_025997.1"] = "5"
      featurepos$chromosome[featurepos$chromosome == "NC_025998.1"] = "6"
      featurepos$chromosome[featurepos$chromosome == "NC_025999.1"] = "7"
      featurepos$chromosome[featurepos$chromosome == "NC_026000.1"] = "8"
      featurepos$chromosome[featurepos$chromosome == "NC_026001.1"] = "9"
      featurepos$chromosome[featurepos$chromosome == "NC_026002.1"] = "10"
      featurepos$chromosome[featurepos$chromosome == "NC_026003.1"] = "11"
      featurepos$chromosome[featurepos$chromosome == "NC_026004.1"] = "12"
      featurepos$chromosome[featurepos$chromosome == "NC_026005.1"] = "13"
      featurepos$chromosome[featurepos$chromosome == "NC_026006.1"] = "14"
      featurepos$chromosome[featurepos$chromosome == "NC_026007.1"] = "15"
      featurepos$chromosome[featurepos$chromosome == "NC_026008.1"] = "16"
    }
    # Vigna unguiculata
    if (accession == "GCF_004118075.1") {
      featurepos$chromosome[featurepos$chromosome == "NC_040279.1"] = "1"
      featurepos$chromosome[featurepos$chromosome == "NC_040280.1"] = "2"
      featurepos$chromosome[featurepos$chromosome == "NC_040281.1"] = "3"
      featurepos$chromosome[featurepos$chromosome == "NC_040282.1"] = "4"
      featurepos$chromosome[featurepos$chromosome == "NC_040283.1"] = "5"
      featurepos$chromosome[featurepos$chromosome == "NC_040284.1"] = "6"
      featurepos$chromosome[featurepos$chromosome == "NC_040285.1"] = "7"
      featurepos$chromosome[featurepos$chromosome == "NC_040286.1"] = "8"
      featurepos$chromosome[featurepos$chromosome == "NC_040287.1"] = "9"
      featurepos$chromosome[featurepos$chromosome == "NC_040288.1"] = "10"
      featurepos$chromosome[featurepos$chromosome == "NC_040289.1"] = "11"
    }
    # Panicum hallii
    if (accession == "GCA_009771035.1") {
      featurepos$chromosome[featurepos$chromosome == "CM019999.1"] = "1"
      featurepos$chromosome[featurepos$chromosome == "CM020000.1"] = "2"
      featurepos$chromosome[featurepos$chromosome == "CM020001.1"] = "3"
      featurepos$chromosome[featurepos$chromosome == "CM020002.1"] = "4"
      featurepos$chromosome[featurepos$chromosome == "CM020003.1"] = "5"
      featurepos$chromosome[featurepos$chromosome == "CM020004.1"] = "6"
      featurepos$chromosome[featurepos$chromosome == "CM020005.1"] = "7"
      featurepos$chromosome[featurepos$chromosome == "CM020006.1"] = "8"
      featurepos$chromosome[featurepos$chromosome == "CM020007.1"] = "9"
      featurepos$chromosome[featurepos$chromosome == "CM020008.1"] = "10"
      featurepos$chromosome[featurepos$chromosome == "CM020009.1"] = "11"
      featurepos$chromosome[featurepos$chromosome == "CM020010.1"] = "12"
      featurepos$chromosome[featurepos$chromosome == "CM020011.1"] = "13"
      featurepos$chromosome[featurepos$chromosome == "CM020012.1"] = "14"
      featurepos$chromosome[featurepos$chromosome == "CM020013.1"] = "15"
      featurepos$chromosome[featurepos$chromosome == "CM020014.1"] = "16"
      featurepos$chromosome[featurepos$chromosome == "CM020015.1"] = "17"
      featurepos$chromosome[featurepos$chromosome == "CM020016.1"] = "18"
      featurepos$chromosome[featurepos$chromosome == "CM020017.1"] = "19"
      featurepos$chromosome[featurepos$chromosome == "CM020018.1"] = "20"
      featurepos$chromosome[featurepos$chromosome == "CM020019.1"] = "21"
      featurepos$chromosome[featurepos$chromosome == "CM020020.1"] = "22"
      featurepos$chromosome[featurepos$chromosome == "CM020021.1"] = "23"
      featurepos$chromosome[featurepos$chromosome == "CM020022.1"] = "24"
      featurepos$chromosome[featurepos$chromosome == "CM020023.1"] = "25"
    }
    # Glycine max
    if (accession == "GCF_000004515.4") {
      featurepos$chromosome[featurepos$chromosome == "NC_016088.3"] = "1"
      featurepos$chromosome[featurepos$chromosome == "NC_016089.3"] = "2"
      featurepos$chromosome[featurepos$chromosome == "NC_016090.3"] = "3"
      featurepos$chromosome[featurepos$chromosome == "NC_016091.3"] = "4"
      featurepos$chromosome[featurepos$chromosome == "NC_038241.1"] = "5"
      featurepos$chromosome[featurepos$chromosome == "NC_038242.1"] = "6"
      featurepos$chromosome[featurepos$chromosome == "NC_038243.1"] = "7"
      featurepos$chromosome[featurepos$chromosome == "NC_038244.1"] = "8"
      featurepos$chromosome[featurepos$chromosome == "NC_038245.1"] = "9"
      featurepos$chromosome[featurepos$chromosome == "NC_038246.1"] = "10"
      featurepos$chromosome[featurepos$chromosome == "NC_038247.1"] = "11"
      featurepos$chromosome[featurepos$chromosome == "NC_038248.1"] = "12"
      featurepos$chromosome[featurepos$chromosome == "NC_038249.1"] = "13"
      featurepos$chromosome[featurepos$chromosome == "NC_038250.1"] = "14"
      featurepos$chromosome[featurepos$chromosome == "NC_038251.1"] = "15"
      featurepos$chromosome[featurepos$chromosome == "NC_038252.1"] = "16"
      featurepos$chromosome[featurepos$chromosome == "NC_038253.1"] = "17"
      featurepos$chromosome[featurepos$chromosome == "NC_038254.1"] = "18"
      featurepos$chromosome[featurepos$chromosome == "NC_038255.1"] = "19"
      featurepos$chromosome[featurepos$chromosome == "NC_038256.1"] = "20"
    }
    # Camelina sativa
    if (accession == "GCF_000633955.1") {
      featurepos$chromosome[featurepos$chromosome == "NC_025685.1"] = "1"
      featurepos$chromosome[featurepos$chromosome == "NC_025686.1"] = "2"
      featurepos$chromosome[featurepos$chromosome == "NC_025687.1"] = "3"
      featurepos$chromosome[featurepos$chromosome == "NC_025688.1"] = "4"
      featurepos$chromosome[featurepos$chromosome == "NC_025689.1"] = "5"
      featurepos$chromosome[featurepos$chromosome == "NC_025690.1"] = "6"
      featurepos$chromosome[featurepos$chromosome == "NC_025691.1"] = "7"
      featurepos$chromosome[featurepos$chromosome == "NC_025692.1"] = "8"
      featurepos$chromosome[featurepos$chromosome == "NC_025693.1"] = "9"
      featurepos$chromosome[featurepos$chromosome == "NC_025694.1"] = "10"
      featurepos$chromosome[featurepos$chromosome == "NC_025695.1"] = "11"
      featurepos$chromosome[featurepos$chromosome == "NC_025696.1"] = "12"
      featurepos$chromosome[featurepos$chromosome == "NC_025697.1"] = "13"
      featurepos$chromosome[featurepos$chromosome == "NC_025698.1"] = "14"
      featurepos$chromosome[featurepos$chromosome == "NC_025699.1"] = "15"
      featurepos$chromosome[featurepos$chromosome == "NC_025700.1"] = "16"
      featurepos$chromosome[featurepos$chromosome == "NC_025701.1"] = "17"
      featurepos$chromosome[featurepos$chromosome == "NC_025702.1"] = "18"
      featurepos$chromosome[featurepos$chromosome == "NC_025703.1"] = "19"
      featurepos$chromosome[featurepos$chromosome == "NC_025704.1"] = "20"
    }
    # Arachis duranensis
    if (accession == "GCF_000817695.2") {
      featurepos$chromosome[featurepos$chromosome == "NC_029772.2"] = "A1"
      featurepos$chromosome[featurepos$chromosome == "NC_029773.2"] = "A2"
      featurepos$chromosome[featurepos$chromosome == "NC_029774.2"] = "A3"
      featurepos$chromosome[featurepos$chromosome == "NC_029775.2"] = "A4"
      featurepos$chromosome[featurepos$chromosome == "NC_029776.2"] = "A5"
      featurepos$chromosome[featurepos$chromosome == "NC_029777.2"] = "A6"
      featurepos$chromosome[featurepos$chromosome == "NC_029778.2"] = "A7"
      featurepos$chromosome[featurepos$chromosome == "NC_029779.2"] = "A8"
      featurepos$chromosome[featurepos$chromosome == "NC_029780.2"] = "A9"
      featurepos$chromosome[featurepos$chromosome == "NC_029781.2"] = "A10"
    }
    # Setaria italica
    if (accession == "GCA_000263155.2") {
      featurepos$chromosome[featurepos$chromosome == "CM003528.1"] = "1"
      featurepos$chromosome[featurepos$chromosome == "CM003529.1"] = "2"
      featurepos$chromosome[featurepos$chromosome == "CM003530.1"] = "3"
      featurepos$chromosome[featurepos$chromosome == "CM003531.1"] = "4"
      featurepos$chromosome[featurepos$chromosome == "CM003532.1"] = "5"
      featurepos$chromosome[featurepos$chromosome == "CM003533.1"] = "6"
      featurepos$chromosome[featurepos$chromosome == "CM003534.1"] = "7"
      featurepos$chromosome[featurepos$chromosome == "CM003535.1"] = "8"
      featurepos$chromosome[featurepos$chromosome == "CM003536.1"] = "9"
    }
    # Camellia sinensis
    if (accession == "GCA_013676235.1") {
      featurepos$chromosome[featurepos$chromosome == "CM025496.1"] = "1"
      featurepos$chromosome[featurepos$chromosome == "CM025497.1"] = "2"
      featurepos$chromosome[featurepos$chromosome == "CM025498.1"] = "3"
      featurepos$chromosome[featurepos$chromosome == "CM025499.1"] = "4"
      featurepos$chromosome[featurepos$chromosome == "CM025500.1"] = "5"
      featurepos$chromosome[featurepos$chromosome == "CM025501.1"] = "6"
      featurepos$chromosome[featurepos$chromosome == "CM025502.1"] = "7"
      featurepos$chromosome[featurepos$chromosome == "CM025503.1"] = "8"
      featurepos$chromosome[featurepos$chromosome == "CM025504.1"] = "9"
      featurepos$chromosome[featurepos$chromosome == "CM025505.1"] = "10"
      featurepos$chromosome[featurepos$chromosome == "CM025506.1"] = "11"
      featurepos$chromosome[featurepos$chromosome == "CM025507.1"] = "12"
      featurepos$chromosome[featurepos$chromosome == "CM025508.1"] = "13"
      featurepos$chromosome[featurepos$chromosome == "CM025509.1"] = "14"
      featurepos$chromosome[featurepos$chromosome == "CM025510.1"] = "15"
    }
    
    # Brassica napus
    if (accession == "GCF_000686985.2") {
      featurepos$chromosome[featurepos$chromosome == "NC_027757.2"] = "A1"
      featurepos$chromosome[featurepos$chromosome == "NC_027758.2"] = "A2"
      featurepos$chromosome[featurepos$chromosome == "NC_027759.2"] = "A3"
      featurepos$chromosome[featurepos$chromosome == "NC_027760.2"] = "A4"
      featurepos$chromosome[featurepos$chromosome == "NC_027761.2"] = "A5"
      featurepos$chromosome[featurepos$chromosome == "NC_027762.2"] = "A6"
      featurepos$chromosome[featurepos$chromosome == "NC_027763.2"] = "A7"
      featurepos$chromosome[featurepos$chromosome == "NC_027764.2"] = "A8"
      featurepos$chromosome[featurepos$chromosome == "NC_027765.2"] = "A9"
      featurepos$chromosome[featurepos$chromosome == "NC_027766.2"] = "A10"
      featurepos$chromosome[featurepos$chromosome == "NC_027767.2"] = "C1"
      featurepos$chromosome[featurepos$chromosome == "NC_027768.2"] = "C2"
      featurepos$chromosome[featurepos$chromosome == "NC_027769.2"] = "C3"
      featurepos$chromosome[featurepos$chromosome == "NC_027770.2"] = "C4"
      featurepos$chromosome[featurepos$chromosome == "NC_027771.2"] = "C5"
      featurepos$chromosome[featurepos$chromosome == "NC_027772.2"] = "C6"
      featurepos$chromosome[featurepos$chromosome == "NC_027773.2"] = "C7"
      featurepos$chromosome[featurepos$chromosome == "NC_027774.2"] = "C8"
      featurepos$chromosome[featurepos$chromosome == "NC_027775.2"] = "C9"
    }
    # Prunus mume
    if (accession == "GCF_000346735.1") {
      featurepos$chromosome[featurepos$chromosome == "NC_024126.1"] = "1"
      featurepos$chromosome[featurepos$chromosome == "NC_024127.1"] = "2"
      featurepos$chromosome[featurepos$chromosome == "NC_024128.1"] = "3"
      featurepos$chromosome[featurepos$chromosome == "NC_024129.1"] = "4"
      featurepos$chromosome[featurepos$chromosome == "NC_024130.1"] = "5"
      featurepos$chromosome[featurepos$chromosome == "NC_024131.1"] = "6"
      featurepos$chromosome[featurepos$chromosome == "NC_024132.1"] = "7"
      featurepos$chromosome[featurepos$chromosome == "NC_024133.1"] = "8"
    }
    # Cucurbita pepo
    if (accession == "GCF_002806865.2") {
      featurepos$chromosome[featurepos$chromosome == "NC_036638.1"] = "1"
      featurepos$chromosome[featurepos$chromosome == "NC_036639.1"] = "2"
      featurepos$chromosome[featurepos$chromosome == "NC_036640.1"] = "3"
      featurepos$chromosome[featurepos$chromosome == "NC_036641.1"] = "4"
      featurepos$chromosome[featurepos$chromosome == "NC_036642.1"] = "5"
      featurepos$chromosome[featurepos$chromosome == "NC_036643.1"] = "6"
      featurepos$chromosome[featurepos$chromosome == "NC_036644.1"] = "7"
      featurepos$chromosome[featurepos$chromosome == "NC_036645.1"] = "8"
      featurepos$chromosome[featurepos$chromosome == "NC_036646.1"] = "9"
      featurepos$chromosome[featurepos$chromosome == "NC_036647.1"] = "10"
      featurepos$chromosome[featurepos$chromosome == "NC_036648.1"] = "11"
      featurepos$chromosome[featurepos$chromosome == "NC_036649.1"] = "12"
      featurepos$chromosome[featurepos$chromosome == "NC_036650.1"] = "13"
      featurepos$chromosome[featurepos$chromosome == "NC_036651.1"] = "14"
      featurepos$chromosome[featurepos$chromosome == "NC_036652.1"] = "15"
      featurepos$chromosome[featurepos$chromosome == "NC_036653.1"] = "16"
      featurepos$chromosome[featurepos$chromosome == "NC_036654.1"] = "17"
      featurepos$chromosome[featurepos$chromosome == "NC_036655.1"] = "18"
      featurepos$chromosome[featurepos$chromosome == "NC_036656.1"] = "19"
      featurepos$chromosome[featurepos$chromosome == "NC_036657.1"] = "20"
    }
    # Sesamum indicum
    if (accession == "GCF_000512975.1") {
      featurepos$chromosome[featurepos$chromosome == "NC_026145.1"] = "1"
      featurepos$chromosome[featurepos$chromosome == "NC_026146.1"] = "2"
      featurepos$chromosome[featurepos$chromosome == "NC_026147.1"] = "3"
      featurepos$chromosome[featurepos$chromosome == "NC_026148.1"] = "4"
      featurepos$chromosome[featurepos$chromosome == "NC_026149.1"] = "5"
      featurepos$chromosome[featurepos$chromosome == "NC_026150.1"] = "6"
      featurepos$chromosome[featurepos$chromosome == "NC_026151.1"] = "7"
      featurepos$chromosome[featurepos$chromosome == "NC_026152.1"] = "8"
      featurepos$chromosome[featurepos$chromosome == "NC_026153.1"] = "9"
      featurepos$chromosome[featurepos$chromosome == "NC_026154.1"] = "10"
      featurepos$chromosome[featurepos$chromosome == "NC_026155.1"] = "11"
      featurepos$chromosome[featurepos$chromosome == "NC_026156.1"] = "12"
      featurepos$chromosome[featurepos$chromosome == "NC_026157.1"] = "13"
      featurepos$chromosome[featurepos$chromosome == "NC_026158.1"] = "14"
      featurepos$chromosome[featurepos$chromosome == "NC_026159.1"] = "15"
      featurepos$chromosome[featurepos$chromosome == "NC_026160.1"] = "16"
    }
    

    # A. speltoides S genome is mapped on A. tauschii D genome
    if (species == "Aegilops_speltoides" & accession == "GCA_000347335.2") {
      featurepos$chromosome = gsub("D", "S", featurepos$chromosome)
    }
    
    
    #----------------------------------------------------------------------------
    # Return a data frame of feature positions
    return(featurepos)
  } else {
    warning("Accession is not available in your dataset.")
  }
  
}

#----------------------------------------------------------------------------
# Retrieve gene positions in GFF files
#----------------------------------------------------------------------------
# ARGUMENTS
#   species = species names in format 'Genus_species'
#   accession = accession number of the genome
#
# VALUES
# Return a data frame with queried informations for gene positions in chromosomes

gene.position = function(species = "", accession = "") {
  # The feature position function is used to retrieve all positions tagged as 'gene' in the GFF
  genpos = feature.position(species, accession, "gene")
  #----------------------------------------------------------------------------
  # Return a data frame of gene positions
  return(genpos)
}

# EXAMPLE
# species = "Arabidopsis_thaliana"
# accession = "GCA_902460315.1"
# res = gene.position(species, accession)

batch.gene.position = function(list_accessions) {
  for (acc in 1:length(list_accessions)) { # LOOP ON ACCESSIONS
    # Getting infos on dataset to treat
    accession = list_accessions[acc]
    species = list_species[acc]
    cat(species, accession, "\n")
    # Retrieving feature position, i.e. genes
    genpos = gene.position(species, accession)
    # Some data is not available or no features have been reported in annotation files (bad annotation)
    if (nrow(genpos) > 0) {
      # if (((genpos) != "Accession is not available in your dataset.")[1]) {
      # Add species and accession to the results
      genpos = cbind(data.frame(species = rep(species, nrow(genpos)), accession = rep(accession, nrow(genpos))),
                     genpos)
      #Load the dataset of the species if it exists, otherwise create a new one
      if (file.exists(paste("data-cleaned/genome/gene_positions/", species,"_gene_positions.txt.gz", sep = ""))) {
        # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
        dataset = read.table(gzfile(paste("data-cleaned/genome/gene_positions/", species,"_gene_positions.txt.gz", sep = "")),
                             header = TRUE, sep = "\t")
      } else {
        # Template of a new empty dataset
        dataset = data.frame(species = character(0), accession = character(0), id = character(0),
                             type = character(0), chromosome = character(0),
                             strand = character(0), phase = numeric(0), start = numeric(0), end = numeric(0))
        write.table(dataset, gzfile(paste("data-cleaned/genome/gene_positions/", species,"_gene_positions.txt.gz", sep = "")), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)      # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
        # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
        dataset = read.table(gzfile(paste("data-cleaned/genome/gene_positions/", species,"_gene_positions.txt.gz", sep = "")),
                             header = TRUE, sep = "\t")
      }
      
      # Remove older rows in the dataset
      dataset = dataset[!(dataset$accession == accession),]
      # Bind new rows
      dataset = rbind(dataset, genpos)
      # Init an empty dataset
      # dataset = data.frame(species = character(0), accession = character(0), id = character(0),
      #                      type = character(0), chromosome = character(0), strand = character(0),
      #                      start = character(0), end = character(0)) 
      # Write updated dataset in 'data-cleaned'
      write.table(dataset, gzfile(paste("data-cleaned/genome/gene_positions/", species,"_gene_positions.txt.gz", sep = "")), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)
      rm(dataset)
      # } else {
      #   warning("Accession is not available in your dataset.")
      # }
    } else {
      warning(paste(species, accession,"No feature found in data.", sep = " "))
    }
    # END OF ACCESSION LOOP
  }
}


#----------------------------------------------------------------------------
# Retrieve exons positions in GFF files
#----------------------------------------------------------------------------
# ARGUMENTS
#   species = species names in format 'Genus_species'
#   accession = accession number of the genome
#
# VALUES
# Return a data frame with queried informations for exons positions in chromosomes

exon.position = function(species = "", accession = "") {
  # The feature position function is used to retrieve all positions tagged as 'exon' in the GFF
  exonpos = feature.position(species, accession, "exon")
  #----------------------------------------------------------------------------
  # Return a data frame of exon positions
  return(exonpos)
}


#----------------------------------------------------------------------------
# Retrieve CDS positions in GFF files
#----------------------------------------------------------------------------
# ARGUMENTS
#   species = species names in format 'Genus_species'
#   accession = accession number of the genome
#
# VALUES
# Return a data frame with queried informations for CDS positions in chromosomes

cds.position = function(species = "", accession = "") {
  # The feature position function is used to retrieve all positions tagged as 'CDS' in the GFF
  cdspos = feature.position(species, accession, "CDS")
  # "For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame.
  # The phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the beginning of this feature to reach the first base of the next codon.
  # In other words, a phase of "0" indicates that the next codon begins at the first base of the region described by the current line, a phase of "1" indicates that the next codon begins at the second base of this region, and a phase of "2" indicates that the codon begins at the third base of this region.
  # This is NOT to be confused with the frame, which is simply start modulo 3. If there is no phase, put a "." (a period) in this field.
  # For forward strand features, phase is counted from the start field. For reverse strand features, phase is counted from the end field." (http://seqanswers.com/forums/showthread.php?t=9500)
  # The phase is required for all CDS features.
  #----------------------------------------------------------------------------
  # Return a data frame of CDS positions
  return(cdspos)
}


batch.cds.position = function(list_accessions) {
  for (acc in 1:length(list_accessions)) { # LOOP ON ACCESSIONS
    # Getting infos on dataset to treat
    accession = list_accessions[acc]
    species = list_species[acc]
    cat(species, accession, "\n")
    # Retrieving feature position, i.e. genes
    cdspos = cds.position(species, accession)
    # Some data is not available or no features have been reported in annotation files (bad annotation)
    if (nrow(cdspos) > 0) {
      # if (((genpos) != "Accession is not available in your dataset.")[1]) {
      # Add species and accession to the results
      cdspos = cbind(data.frame(species = rep(species, nrow(cdspos)), accession = rep(accession, nrow(cdspos))),
                     cdspos)
      #Load the dataset of the species if it exists, otherwise create a new one
      if (file.exists(paste("data-cleaned/genome/cds_positions/", species,"_cds_positions.txt.gz", sep = ""))) {
        # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
        dataset = read.table(gzfile(paste("data-cleaned/genome/cds_positions/", species,"_cds_positions.txt.gz", sep = "")),
                             header = TRUE, sep = "\t")
      } else {
        # Template of a new empty dataset
        dataset = data.frame(species = character(0), accession = character(0), id = character(0),
                             type = character(0), chromosome = character(0),
                             strand = character(0), phase = numeric(0), start = numeric(0), end = numeric(0))
        write.table(dataset, gzfile(paste("data-cleaned/genome/cds_positions/", species,"_cds_positions.txt.gz", sep = "")), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)      # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
        # Load the dataset of gene positions for the species (for all accessions of this species) in 'data-cleaned'
        dataset = read.table(gzfile(paste("data-cleaned/genome/cds_positions/", species,"_cds_positions.txt.gz", sep = "")),
                             header = TRUE, sep = "\t")
      }
      
      # Remove older rows in the dataset
      dataset = dataset[!(dataset$accession == accession),]
      # Bind new rows
      dataset = rbind(dataset, cdspos)
      # Init an empty dataset
      # dataset = data.frame(species = character(0), accession = character(0), id = character(0),
      #                      type = character(0), chromosome = character(0), strand = character(0),
      #                      start = character(0), end = character(0)) 
      # Write updated dataset in 'data-cleaned'
      write.table(dataset, gzfile(paste("data-cleaned/genome/cds_positions/", species,"_cds_positions.txt.gz", sep = "")), sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)
      rm(dataset)
      # } else {
      #   warning("Accession is not available in your dataset.")
      # }
    } else {
      warning(paste(species, accession,"No feature found in data.", sep = " "))
    }
    # END OF ACCESSION LOOP
  }
}

#============================================================================
#   GC content in recombination windows
#============================================================================
# Method 'map'
# Global composition in GC, in windows of 100kb and 1Mb of same coordinates as recombination maps
# Composition in 1rst, 2nd and 3rd codon position in the CDS contained by the 100kb window (i.e. CDS with start position in the window)
# In order to be explanatory variables for recombination landscapes, windows need to be the same as used for recombination rates
# Hence, importing each recombination maps to compute GC proportions within defined windows

# The R package 'seqinr' have functions to compute GC, GC1, GC2 and GC3 given a nucleic sequence
# https://www.rdocumentation.org/packages/seqinr/versions/3.6-1/topics/G%2BC%20Content
# Another way to implement GC proportions in sliding windows: https://web.stanford.edu/class/bios221/labs/biostrings/lab_1_biostrings.html

# Besides,
# Method 'gene'

# Method 'exon'

# USAGE
# gc.map(map = "", scale = 0.1, method = c("map", "gene", "exon"))

# ARGUMENTS
# map = the recombination map on which you want to estimate the local GC proportions
# scale in Mb (e.g. 0.1 for 100kb windows)
# ex = TRUE, the argument 'exact' given to GC function to estimate GC proportion
# method = c("map", "gene", "exon"); for 'map', compute statistics in windows of the recombination map
#         for 'gene', compute statistics for gene sequences (CDS for GC1, GC2, GC3)
#         for 'exon', compute statistics for exon sequences (CDS for GC1, GC2, GC3)
# ncors = 15, the maximum number of cores to use, overrided by the number of cores detected in the hardware
# VALUES
# Return the exact recombination map + statistics computed on the genome
# GC
# GC1
# GC2
# GC3
# GC1 first exon
# GC2 first exon
# GC3 first exon

gc.map = function(map = "", scale = 0.1, ex = TRUE, method = "map", ncores = 15) {
  require(ape)
  require(seqinr)
  require(rBLAST)
  require(Biostrings)
  # require(foreach)
  # require(doParallel)
  # require(doSNOW)
  require(parallel)
  require(MASS)
  require(pbmcapply)
  
  # TODO optim -> process either single map or a list of map
  # TODO optim -> parallel
  # TODO can use a vector of methods, multiple estimates in a single run
  
  cat("Processing ", as.character(map), "\n")
  # Load the db for the species, given a species name and accession number
  metadata.clean = read.table("data-cleaned/Genetic_maps_ressources.csv", header = TRUE, sep = ";")
  species = as.character(metadata.clean$species[metadata.clean$id == map])
  accession = as.character(metadata.clean$accession[metadata.clean$id == map])
  # Import list of genes, only for method = 'gene'
  if (method == "gene") {
    gene_positions = read.table(gzfile(paste("data-cleaned/genome/gene_positions/", species,"_gene_positions.txt.gz", sep = "")), header = TRUE, sep = "\t")
  }
  # Import list of CDS and exons
  CDS_positions = read.table(gzfile(paste("data-cleaned/genome/cds_positions/",species,"_cds_positions.txt.gz", sep = "")), header = TRUE, sep = "\t")
  exon_positions = read.table(gzfile(paste("data-cleaned/genome/exon_positions/", species,"_exon_positions.txt.gz", sep = "")), header = TRUE, sep = "\t")
  
  # Names of fasta files susceptible to contain sequence data
  ls = list.files(path = paste(wd, "/data/Genome/", tolower(species),"/" ,accession,"/fasta/", sep = ""))
  ls = ls[grepl("fna.gz", ls) | grepl("fa.gz", ls) | grepl("fasta.gz", ls)]

  # If no fasta file found, i.e. no genomic data
  if (length(ls) == 0) {
    return("No fasta file for this dataset.")
  } else {  # Open the fasta file and concatenate them if multiple fasta files: make a single data for all chromosomes
    seq = readDNAStringSet(paste(wd, "/data/Genome/", tolower(species),"/" ,accession,"/fasta/", ls, sep = ""))
  }
  
  # Retrieve sequence names in a vector
  # Sequences in DNAStringSet are...
  chrnames = names(seq)
  # Sometimes chromosome names are given in latin numerals
  # So, give it a translation before...
  source("sources/core_pipeline.R") # Basic functions but essential for the pipeline (e.g. saving logs)
  if (sum(grepl("[IVX]+", chrnames)) > 0) {
    chrnames = latin2arabic(chrnames)
  }
  # Reformat to chromosome name
  chrinseq = gsub("chromosome:*[ ]*[A-Za-z]*", "",regmatches(chrnames, regexpr("chromosome:*[ ]*[A-Za-z]*[0-9]+[A-Z]*", chrnames)))
  # Yet, chromosome name may also be in format ">[0-9]+"
  # Or, like in Zea mays, 'chromosome:' contain genome name, not chromosome name. Discard if chrinseq contains multiple occurences of a single name
  if (length(chrinseq) == 0 | (length(unique(chrinseq)) == 1 & length(chrinseq) > 5)) {
    chrinseq = gsub(" dna:chromosome", "", gsub(">", "",regmatches(chrnames, regexpr("^[>]*[A-Za-z]*[0-9]+[A-Za-z]* dna:chromosome", chrnames))))
  }
  # Exceptions in chromosome naming
  # Some species have different naming schemes, depending on sources
  # Translate chromosome names for some species to the common used naming in our dataset
  # Glycine max
  if (grepl("Glycine_max", species)) {
    # TODO Code optimization, very long process although it should be short
    translate_gly = data.frame(A1 = 5, A2 = 8, B1 = 11, B2 = 14, C1 = 4, C2 = 6, D1a = 1, D1b = 2,
                               D2 = 17, E = 15, F = 13, G = 18, H = 12, I = 20, J = 16, K = 9, L =19, M = 7, N = 3, O = 10)
    translation = function(x) {return(as.character(colnames(translate_gly[which(x == translate_gly)])))}
    chrinseq = lapply(chrinseq, translation)
    CDS_positions$chromosome = as.character(CDS_positions$chromosome)
    CDS_positions$chromosome = lapply(CDS_positions$chromosome, translation)
    rm(translate_gly)
    rm(translation)
  }

  # For each chromosome for which a recombination map exists
  list = system(paste("ls ", wd, "/output/recombination_maps/loess/100kbwind/", sep = ""), intern = TRUE)
  list = list[grep(map, list)]
  list = gsub(".txt", "", list)
  
  # Detect the number of cores to use in parallel tasks
  # Either the max number specified by user of actual number of cores in hardware - 1
  numcores = detectCores()
  ncores = min(numcores - 1, ncores)
  
  #-------------------------------------------
  # Declare functions to estimate GC
  # the getRecRate() function retrieve the recombination rate of the windows in which the feature is
  # Take the index of the features in argument and return a vector of values
  # data must be given, a gene_position or exon_position dataset
  getRecRate = function(g, data) {
    # Query the local recombination rate of a gene by its physical position
    # i.e. windows in which start position is
    # And in any case, prevent for NA estimates out of boundaries
    if ((data$start[g] < min((mapdata$phys*1000000))) | (data$start[g] > max((mapdata$phys*1000000)))) {
      rec.rate = NA
    } else {
      rec.rate = mapdata$rec.rate[max(which(data$start[g] > (mapdata$phys*1000000)))]
    }
    
    return(rec.rate)
  }
  # the getGC() function estimates GC/GC1/GC2/GC3 contents in the features provided
  # Take the index of the features in argument and return a vector of values
  # data must be given, a gene_position or exon_position dataset
  getGC = function(g, data) {
    # Estimates GC on the whole gene sequence
    # First retrieve the sequence with gene coordinates
    # MEAN GC, sequence of the complete gene sequence
    seqwin = seq_wholechr[data$start[g]:min(c(data$end[g], length(seq_wholechr)))]
    # NAs produced in tails of the sequence
    seqwin[is.na(seqwin)] = "N"
    seqwin = s2c(tolower(as.character(seqwin))) # Convert to a vector of single characters
    # Compute local GC proportion along the complete gene sequence
    if (length(seqwin) > 0) {
      # forcetoLower = FALSE require to convert sequence to lower case before, yet it reduces the computation time by 3
      # exact computation estimates GC content with ambiguous bases taken in account. Time is increased by 3.
      # Yet, it revealed differences between estimations with exact = FALSE
      gc = GC(seqwin, forceToLower = FALSE, exact = ex)
      # distGC = c(distGC, GC(seqwin, forceToLower = FALSE, exact = ex))
    } else {
      gc = NA
    }
    
    return(gc)
  }
  
  # Get the GC content in codon position i, i given by the 'position' argument
  getGCi = function(g, data, position = 1) {
    # GC IN CODONS
    # Trim sequences to CDS contained in map windows in a sequence object
    seqCDS = data.frame(start = CDS_positions$start[CDS_positions$accession == accession & CDS_positions$chromosome == chrinrecmaps],
                        end = CDS_positions$end[CDS_positions$accession == accession & CDS_positions$chromosome == chrinrecmaps],
                        strand = CDS_positions$strand[CDS_positions$accession == accession & CDS_positions$chromosome == chrinrecmaps],
                        phase = CDS_positions$phase[CDS_positions$accession == accession & CDS_positions$chromosome == chrinrecmaps])
    
    # Retain only CDS associated to the exon
    windowed_CDS = seqCDS[(seqCDS$start >= data$start[g]) & (seqCDS$start <= data$end[g]),]
    # Init vector of results: sum of GCi proportions for each sequence
    mean_gci = numeric(nrow(windowed_CDS))
    
    # cat("GC proportions in CDS")
    # If some CDS have been selected in the current window...
    if (nrow(windowed_CDS) > 0) {
      # LOOP ON CDS
      # For each CDS within the window
      for (c in 1:nrow(windowed_CDS)) {
        # print(c)
        # Retrieve the fasta sequence of the window; respects bounds of the sequence at the end
        # If subset is out of bounds of reference genome sequence, then add NA values
        if (data$start[g] < length(seq_wholechr)) {
          seqwin = seq_wholechr[windowed_CDS$start[c]:min(c(windowed_CDS$end[c], length(seq_wholechr)))]
          seqwin[is.na(seqwin)] = "N"
          seqwin = s2c(tolower(as.character(seqwin))) # Convert to a vector of single characters in lower case
          # IF reverse strand (strand = -)
          if (windowed_CDS$strand[c] == "-") {
            seqwin = rev(seqwin)
          }
          if ((length(seqwin) - windowed_CDS$phase[c]) >= 3) {
            mean_gci[c] = GCpos(seqwin, frame = windowed_CDS$phase[c],
                                forceToLower = FALSE, exact = ex, pos = position)
          }
        } else {
          mean_gci[c] = NA
        }
        # END OF CDS LOOP
      }
      
      # Compute mean GC proportions according to the position in codons on a list of CDS sequences
      gci = mean(mean_gci, na.rm = TRUE)
      
      # distGC3 = c(distGC3, mean_gc3)
    } else { # Sometimes there is no CDS recovered within the window
      gci = NA
    }
    rm(mean_gci)
    
    return(gci)
  }
  
  #-------------------------------------------
  if (method == "map") {
    for (chromosome in list) {
      # chromosome = list[1] # for debug purpose only
      # print(chromosome) # print for debug purpose only
      # Init vectors of results
      gc = numeric()
      gc1 = numeric()
      gc2 = numeric()
      gc3 = numeric()
      gc1ex1 = numeric()
      gc2ex1 = numeric()
      gc3ex1 = numeric()
      # distGC = numeric()
      # distGC3 = numeric()
      # Load the recombination map
      mapdata = read.table(file = paste("output/recombination_maps/loess/100kbwind/", chromosome, ".txt", sep = ""), header = TRUE, sep = "\t")
      # Retrieve the sequence for the given chromosome in the list
      chrinrecmaps = gsub("chromosome", "", regmatches(chromosome, regexpr("chromosome[A-Za-z]*[0-9]*[A-Za-z]*", chromosome)))
      cat("Chromosome ", chrinrecmaps, "\n")
      idx = which(chrinseq == chrinrecmaps)
      # Select the whole chromosome sequence for further analyses
      # In case of multiple fasta file, select arbitrary the first one
      seq_wholechr = seq[[idx[1]]]
      
      # Now, Create windows in bp
      # windows = seq(0, max(mapdata$phys), scale)
      # windows = c(windows, max(windows) + scale)
      # windows
      windows.start = (mapdata$phys - (scale/2))*1000000
      windows.stop = (mapdata$phys + (scale/2))*1000000
      
      #------------------------------
      # GC proportions in codons
      
      # Trim sequences to CDS contained in map windows in a sequence object
      seqCDS = data.frame(start = CDS_positions$start[CDS_positions$accession == accession & CDS_positions$chromosome == chrinrecmaps],
                          end = CDS_positions$end[CDS_positions$accession == accession & CDS_positions$chromosome == chrinrecmaps],
                          strand = CDS_positions$strand[CDS_positions$accession == accession & CDS_positions$chromosome == chrinrecmaps],
                          phase = CDS_positions$phase[CDS_positions$accession == accession & CDS_positions$chromosome == chrinrecmaps])
      
      # Gather the end position of CDS of first exons for each gene
      # exon_positions$end[exon_positions$rank == 1 & exon_positions$accession == accession & exon_positions$chromosome == chrinrecmaps]
      rank1exons = exon_positions$end[exon_positions$rank == 1 & exon_positions$accession == accession & exon_positions$chromosome == chrinrecmaps]
      
      # A LOOP ON EACH WINDOWS
      # For each 100kb window in the map
      pb = txtProgressBar(min = 0, max = length(windows.start), initial = 0)
      for (w in 1:length(windows.start)) {
        # print(w)
        setTxtProgressBar(pb, w)
        # cat("Average GC proportions")
        # Retrieve the fasta sequence of the window; respects bounds of the sequence at the end
        # e.g. reference genome slightly shorter than recombination map
        # If subset is out of bounds of reference genome sequence, then add NA values
        if (windows.start[w] < length(seq_wholechr)) {
          # MEAN GC, sequence of the complete window (100kb)
          seqwin = seq_wholechr[windows.start[w]:min(c(windows.stop[w], length(seq_wholechr)))]
          # NAs produced in tails of the sequence
          seqwin[is.na(seqwin)] = "N"
          seqwin = s2c(tolower(as.character(seqwin))) # Convert to a vector of single characters
          # Compute local GC proportion along the complete window
          if (length(seqwin) > 0) {
            # forcetoLower = FALSE require to convert sequence to lower case before, yet it reduces the computation time by 3
            # exact computation estimates GC content with ambiguous bases taken in account. Time is increased by 3.
            # Yet, it revealed differences between estimations with exact = FALSE
            gc[w] = GC(seqwin, forceToLower = FALSE, exact = ex)
            # distGC = c(distGC, GC(seqwin, forceToLower = FALSE, exact = ex))
          } else {
            gc[w] = NA
          }
          
          # GC IN CODONS
          # Retain only CDS within the window, discard CDS beginning in the previous window
          windowed_CDS = seqCDS[(seqCDS$start >= windows.start[w]) & (seqCDS$start <= windows.stop[w]),]
          # Init vector of results: sum of GCi proportions for each sequence
          mean_gc1 = numeric(nrow(windowed_CDS))
          mean_gc2 = numeric(nrow(windowed_CDS))
          mean_gc3 = numeric(nrow(windowed_CDS))
          
          # cat("GC proportions in CDS")
          # If some CDS have been selected in the current window...
          if (nrow(windowed_CDS) > 0) {
            # LOOP ON CDS
            # For each CDS within the window
            for (c in 1:nrow(windowed_CDS)) {
              # print(c)
              # Retrieve the fasta sequence of the window; respects bounds of the sequence at the end
              # If subset is out of bounds of reference genome sequence, then add NA values
              if (windowed_CDS$start[c] < length(seq_wholechr)) {
                seqwin = seq_wholechr[windowed_CDS$start[c]:min(c(windowed_CDS$end[c], length(seq_wholechr)))]
                seqwin[is.na(seqwin)] = "N"
                seqwin = s2c(tolower(as.character(seqwin))) # Convert to a vector of single characters in lower case
                # IF reverse strand (strand = -)
                if (windowed_CDS$strand[c] == "-") {
                  seqwin = rev(seqwin)
                }
                if ((length(seqwin) - windowed_CDS$phase[c]) >= 3) {
                  mean_gc1[c] = GC1(seqwin, frame = windowed_CDS$phase[c], forceToLower = FALSE, exact = ex)
                  mean_gc2[c] = GC2(seqwin, frame = windowed_CDS$phase[c], forceToLower = FALSE, exact = ex)
                  mean_gc3[c] = GC3(seqwin, frame = windowed_CDS$phase[c], forceToLower = FALSE, exact = ex)
                  # mean_gc1 = c(mean_gc1, GC1(seqwin, forceToLower = FALSE), exact = TRUE)
                  # mean_gc2 = c(mean_gc2, GC2(seqwin, forceToLower = FALSE), exact = TRUE)
                  # mean_gc3 = c(mean_gc3, GC3(seqwin, forceToLower = FALSE), exact = TRUE)
                }
              } else {
                mean_gc1[c] = NA
                mean_gc2[c] = NA
                mean_gc3[c] = NA
              }
              # END OF CDS LOOP
            }
            
            # Compute mean GC proportions according to the position in codons on a list of CDS sequences
            # mean(mean_gc1, na.rm = TRUE)
            # mean(mean_gc2, na.rm = TRUE)
            # mean(mean_gc3, na.rm = TRUE)
            gc1[w] = mean(mean_gc1, na.rm = TRUE)
            gc2[w] = mean(mean_gc2, na.rm = TRUE)
            gc3[w] = mean(mean_gc3, na.rm = TRUE)
            # distGC3 = c(distGC3, mean_gc3)
          } else { # Sometimes there is no CDS recovered within the window
            gc1[w] = NA
            gc2[w] = NA
            gc3[w] = NA
          }
          rm(mean_gc1, mean_gc2, mean_gc3)
          
          # GCi in first exon ONLY
          # Retain only CDS within the window
          windowed_CDS = seqCDS[(seqCDS$start >= windows.start[w]) & (seqCDS$start <= windows.stop[w]),]
          # Then, find exons with matching end position for the CDS
          # Keep only rank = 1
          # i.e. CDS of first exon
          windowed_CDS = windowed_CDS[windowed_CDS$end %in% rank1exons,]
          # Init vector of results: sum of GCi proportions for each sequence
          mean_gc1ex1 = numeric(nrow(windowed_CDS))
          mean_gc2ex1 = numeric(nrow(windowed_CDS))
          mean_gc3ex1 = numeric(nrow(windowed_CDS))
          
          if (nrow(windowed_CDS) > 0) {
            # For each CDS within the window
            for (c in 1:nrow(windowed_CDS)) {
              # print(c)
              # Retrieve the fasta sequence of the window; respects bounds of the sequence at the end
              # If subset is out of bounds of reference genome sequence, then add NA values
              if (windowed_CDS$start[c] < length(seq_wholechr)) {
                seqwin = seq_wholechr[windowed_CDS$start[c]:min(c(windowed_CDS$end[c], length(seq_wholechr)))]
                seqwin[is.na(seqwin)] = "N"
                seqwin = s2c(tolower(as.character(seqwin))) # Convert to a vector of single characters
                # IF reverse strand (strand = -)
                if (windowed_CDS$strand[c] == "-") {
                  seqwin = rev(seqwin)
                }
                if ((length(seqwin) - windowed_CDS$phase[c]) >= 3) {
                  # mean_gc1ex1 = c(mean_gc1ex1, GC1(seqwin, forceToLower = FALSE), exact = TRUE)
                  # mean_gc2ex1 = c(mean_gc2ex1, GC2(seqwin, forceToLower = FALSE), exact = TRUE)
                  # mean_gc3ex1 = c(mean_gc3ex1, GC3(seqwin, forceToLower = FALSE), exact = TRUE)
                  mean_gc1ex1[c] = GC1(seqwin, frame = windowed_CDS$phase[c], forceToLower = FALSE, exact = ex)
                  mean_gc2ex1[c] = GC2(seqwin, frame = windowed_CDS$phase[c], forceToLower = FALSE, exact = ex)
                  mean_gc3ex1[c] = GC3(seqwin, frame = windowed_CDS$phase[c], forceToLower = FALSE, exact = ex)
                }
              } else {
                mean_gc1ex1[c] = NA
                mean_gc2ex1[c] = NA
                mean_gc3ex1[c] = NA
              }
              # END OF EXON LOOP
            }
            
            # Compute mean GC proportions according to the position in codons on a list of CDS sequences
            # mean(mean_gc1ex1, na.rm = TRUE)
            # mean(mean_gc2ex1, na.rm = TRUE)
            # mean(mean_gc3ex1, na.rm = TRUE)
            gc1ex1[w] = mean(mean_gc1ex1, na.rm = TRUE)
            gc2ex1[w] = mean(mean_gc2ex1, na.rm = TRUE)
            gc3ex1[w] = mean(mean_gc3ex1, na.rm = TRUE)
            
          } else { # Sometimes there is no CDS recovered within the window
            gc1ex1[w] = NA
            gc2ex1[w] = NA
            gc3ex1[w] = NA
          }
          rm(mean_gc1ex1, mean_gc2ex1, mean_gc3ex1)
        } else { # Iterator w Out of the reference genome sequence
          gc[w] = NA
          gc1[w] = NA
          gc2[w] = NA
          gc3[w] = NA
          gc1ex1[w] = NA
          gc2ex1[w] = NA
          gc3ex1[w] = NA
        }
        
        
        
        
        # END OF WINDOWS LOOP - iterator w
      }
      close(pb)
      # DEBUG RESULTS
      # plot(gc_res)
      # plot(gc1_res)
      # plot(gc2_res)
      # plot(gc3_res)
      
      # Write directly mapped proportions
      # mapdata = cbind(mapdata, gc, gc1, gc2, gc3)
      mapdata = cbind(mapdata, gc, gc1, gc2, gc3, gc1ex1, gc2ex1, gc3ex1)
      
      # hist(distGC)
      # hist(distGC3)
      # mapdata = cbind(mapdata, gc)
      
      # File is saved in a txt data frame and gunzipped for efficient storage management
      # See https://rce-docs.hmdc.harvard.edu/faq/how-do-i-use-compressed-data-r
      write.table(mapdata, file = gzfile(paste("data-cleaned/genome/gc_maps/", chromosome, ".txt.gz", sep = "")),
                  quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
      
      # END OF CHROMOSOME LOOP
    }
  }
  
  #-------------------------------------------
  # COMPUTE GC CONTENT FOR GENES
                  # if (method == "gene") {
                  #   # Init empty result dataframe
                  #   gc4genes = data.frame(species = gene_positions$species, accession = gene_positions$accession, id = gene_positions$id,
                  #                         chromosome = gene_positions$chromosome, rec.rate = numeric(nrow(gene_positions)),
                  #                         gc = numeric(nrow(gene_positions)), gc1 = numeric(nrow(gene_positions)), gc2 = numeric(nrow(gene_positions)),
                  #                         gc3 = numeric(nrow(gene_positions)))
                  #   # Remove any line that is not in the chromosome list
                  #   listchr = gsub("_chromosome", "", regmatches(list, gregexpr("_chromosome[A-Za-z0-9[:punct:]]+$", list)))
                  #   gc4genes = gc4genes[which((gc4genes$chromosome %in% listchr)),]
                  #   
                  #   # LOOP ON EACH CHROMOSOME
                  #   for (chromosome in list) {
                  #     # Load the recombination map
                  #     mapdata = read.table(file = paste("output/recombination_maps/loess/100kbwind/", chromosome, ".txt", sep = ""), header = TRUE, sep = "\t")
                  #     # Work only on the current chromosome
                  #     chrname = gsub("_chromosome", "", regmatches(chromosome, gregexpr("_chromosome[A-Za-z0-9[:punct:]]+$", chromosome)))
                  #     # Retrieve the sequence for the given chromosome in the list
                  #     chrinrecmaps = gsub("chromosome", "", regmatches(chromosome, regexpr("chromosome[A-Za-z]*[0-9]*[A-Za-z]*", chromosome)))
                  #     cat("Chromosome ", chrinrecmaps, "\n")
                  #     idx = which(chrinseq == chrinrecmaps)
                  #     # Select the whole chromosome sequence for further analyses
                  #     # In case of multiple fasta file, select arbitrary the first one
                  #     seq_wholechr = seq[[idx[1]]]
                  #     
                  #     # LOOP ON GENES Compute statistics for each position
                  #     # Index of genes to process
                  #     idx.genes = which(gene_positions$chromosome == chrname)
                  #     pb = txtProgressBar(min = 0, max = length(idx.genes), initial = 0)
                  #     count = 0
                  #     for (g in idx.genes) {
                  #       count = count + 1 # Iterator for progress bar
                  #       setTxtProgressBar(pb, count)
                  #       # Query the local recombination rate of a gene by its physical position
                  #       # i.e. windows in which start position is
                  #       # And in any case, prevent for NA estimates out of boundaries
                  #       if ((gene_positions$start[g] < min((mapdata$phys*1000000))) | (gene_positions$start[g] > max((mapdata$phys*1000000)))) {
                  #         gc4genes$rec.rate[g] = NA
                  #       } else {
                  #         gc4genes$rec.rate[g] = mapdata$rec.rate[max(which(gene_positions$start[g] > (mapdata$phys*1000000)))]
                  #       }
                  #       
                  #       # Next estimates GC on the whole gene sequence
                  #       # First retrieve the sequence with gene coordinates
                  #       # MEAN GC, sequence of the complete gene sequence
                  #       seqwin = seq_wholechr[gene_positions$start[g]:min(c(gene_positions$end[g], length(seq_wholechr)))]
                  #       # NAs produced in tails of the sequence
                  #       seqwin[is.na(seqwin)] = "N"
                  #       seqwin = s2c(tolower(as.character(seqwin))) # Convert to a vector of single characters
                  #       # Compute local GC proportion along the complete gene sequence
                  #       if (length(seqwin) > 0) {
                  #         # forcetoLower = FALSE require to convert sequence to lower case before, yet it reduces the computation time by 3
                  #         # exact computation estimates GC content with ambiguous bases taken in account. Time is increased by 3.
                  #         # Yet, it revealed differences between estimations with exact = FALSE
                  #         gc4genes$gc[g] = GC(seqwin, forceToLower = FALSE, exact = ex)
                  #         # distGC = c(distGC, GC(seqwin, forceToLower = FALSE, exact = ex))
                  #       } else {
                  #         gc4genes$gc[g] = NA
                  #       }
                  # 
                  #       # GC IN CODONS
                  #       # Trim sequences to CDS contained in map windows in a sequence object
                  #       seqCDS = data.frame(start = CDS_positions$start[CDS_positions$accession == accession & CDS_positions$chromosome == chrinrecmaps],
                  #                           end = CDS_positions$end[CDS_positions$accession == accession & CDS_positions$chromosome == chrinrecmaps],
                  #                           strand = CDS_positions$strand[CDS_positions$accession == accession & CDS_positions$chromosome == chrinrecmaps],
                  #                           phase = CDS_positions$phase[CDS_positions$accession == accession & CDS_positions$chromosome == chrinrecmaps])
                  #       
                  #       # Retain only CDS within the window, discard CDS beginning in the previous window
                  #       windowed_CDS = seqCDS[(seqCDS$start >= gene_positions$start[g]) & (seqCDS$start <= gene_positions$end[g]),]
                  #       # Init vector of results: sum of GCi proportions for each sequence
                  #       mean_gc1 = numeric(nrow(windowed_CDS))
                  #       mean_gc2 = numeric(nrow(windowed_CDS))
                  #       mean_gc3 = numeric(nrow(windowed_CDS))
                  #       
                  #       # cat("GC proportions in CDS")
                  #       # If some CDS have been selected in the current window...
                  #       if (nrow(windowed_CDS) > 0) {
                  #         # LOOP ON CDS
                  #         # For each CDS within the window
                  #         for (c in 1:nrow(windowed_CDS)) {
                  #           # print(c)
                  #           # Retrieve the fasta sequence of the window; respects bounds of the sequence at the end
                  #           # If subset is out of bounds of reference genome sequence, then add NA values
                  #           if (gene_positions$start[g] < length(seq_wholechr)) {
                  #             seqwin = seq_wholechr[windowed_CDS$start[c]:min(c(windowed_CDS$end[c], length(seq_wholechr)))]
                  #             seqwin[is.na(seqwin)] = "N"
                  #             seqwin = s2c(tolower(as.character(seqwin))) # Convert to a vector of single characters in lower case
                  #             # IF reverse strand (strand = -)
                  #             if (windowed_CDS$strand[c] == "-") {
                  #               seqwin = rev(seqwin)
                  #             }
                  #             if ((length(seqwin) - windowed_CDS$phase[c]) >= 3) {
                  #               mean_gc1[c] = GC1(seqwin, frame = windowed_CDS$phase[c], forceToLower = FALSE, exact = ex)
                  #               mean_gc2[c] = GC2(seqwin, frame = windowed_CDS$phase[c], forceToLower = FALSE, exact = ex)
                  #               mean_gc3[c] = GC3(seqwin, frame = windowed_CDS$phase[c], forceToLower = FALSE, exact = ex)
                  #               # mean_gc1 = c(mean_gc1, GC1(seqwin, forceToLower = FALSE), exact = TRUE)
                  #               # mean_gc2 = c(mean_gc2, GC2(seqwin, forceToLower = FALSE), exact = TRUE)
                  #               # mean_gc3 = c(mean_gc3, GC3(seqwin, forceToLower = FALSE), exact = TRUE)
                  #             }
                  #           } else {
                  #             mean_gc1[c] = NA
                  #             mean_gc2[c] = NA
                  #             mean_gc3[c] = NA
                  #           }
                  #           # END OF CDS LOOP
                  #         }
                  #         
                  #         # Compute mean GC proportions according to the position in codons on a list of CDS sequences
                  #         # mean(mean_gc1, na.rm = TRUE)
                  #         # mean(mean_gc2, na.rm = TRUE)
                  #         # mean(mean_gc3, na.rm = TRUE)
                  #         gc4genes$gc1[g] = mean(mean_gc1, na.rm = TRUE)
                  #         gc4genes$gc2[g] = mean(mean_gc2, na.rm = TRUE)
                  #         gc4genes$gc3[g] = mean(mean_gc3, na.rm = TRUE)
                  #         # distGC3 = c(distGC3, mean_gc3)
                  #       } else { # Sometimes there is no CDS recovered within the window
                  #         gc4genes$gc1[g] = NA
                  #         gc4genes$gc2[g] = NA
                  #         gc4genes$gc3[g] = NA
                  #       }
                  #       rm(mean_gc1, mean_gc2, mean_gc3)
                  #       
                  #       # END OF LOOP ON EACH POSITION
                  #     }
                  #     close(pb)
                  # 
                  #   # END OF LOOP ON CHROMOSOMES
                  #   }
                  # 
                  #   # File is saved in a txt data frame and gunzipped for efficient storage management
                  #   # See https://rce-docs.hmdc.harvard.edu/faq/how-do-i-use-compressed-data-r
                  #   cat("Writing file...\n")
                  #   write.table(gc4genes, file = gzfile(paste("data-cleaned/genome/gc_genes/", map, ".txt.gz", sep = "")),
                  #               quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
                  # 
                  #   # END OF METHOD 'GENE'
                  # }
  if (method == "gene") {
    # Init empty result dataframe
    gc4genes = data.frame(species = gene_positions$species, accession = gene_positions$accession, id = gene_positions$id,
                          chromosome = gene_positions$chromosome, rec.rate = numeric(nrow(gene_positions)),
                          gc = numeric(nrow(gene_positions)), gc1 = numeric(nrow(gene_positions)), gc2 = numeric(nrow(gene_positions)),
                          gc3 = numeric(nrow(gene_positions)))
    # Remove any line that is not in the chromosome list
    listchr = gsub("_chromosome", "", regmatches(list, gregexpr("_chromosome[A-Za-z0-9[:punct:]]+$", list)))
    gc4genes = gc4genes[which((gc4genes$chromosome %in% listchr)),]
    
    # LOOP ON EACH CHROMOSOME
    for (chromosome in list) {
      # Load the recombination map
      mapdata = read.table(file = paste("output/recombination_maps/loess/100kbwind/", chromosome, ".txt", sep = ""), header = TRUE, sep = "\t")
      # Work only on the current chromosome
      chrname = gsub("_chromosome", "", regmatches(chromosome, gregexpr("_chromosome[A-Za-z0-9[:punct:]]+$", chromosome)))
      # Retrieve the sequence for the given chromosome in the list
      chrinrecmaps = gsub("chromosome", "", regmatches(chromosome, regexpr("chromosome[A-Za-z]*[0-9]*[A-Za-z]*", chromosome)))
      cat("Chromosome ", chrinrecmaps, "\n")
      idx = which(chrinseq == chrinrecmaps)
      # Select the whole chromosome sequence for further analyses
      # In case of multiple fasta file, select arbitrary the first one
      seq_wholechr = seq[[idx[1]]]
      
      # LOOP ON GENES Compute statistics for each position
      # Index of genes to process
      idx.genes = which(gene_positions$chromosome == chrname)
      
      cat("Recovering Recombination rates... ")
      gc4genes$rec.rate[idx.genes] = pbmclapply(idx.genes, function(X) getRecRate(X, CDS_positions), mc.cores = ncores)
      gc4genes$rec.rate = unlist(gc4genes$rec.rate)
      cat("Processing GC... ")
      gc4genes$gc[idx.genes] = pbmclapply(idx.genes, function(X) getGC(X, CDS_positions), mc.cores = ncores)
      gc4genes$gc = unlist(gc4genes$gc)
      cat("GC1... ")
      gc4genes$gc1[idx.genes] = pbmclapply(idx.genes, function(X) getGCi(X, CDS_positions, position = 1), mc.cores = ncores)
      gc4genes$gc1 = unlist(gc4genes$gc1)
      cat("GC2... ")
      gc4genes$gc2[idx.genes] = pbmclapply(idx.genes, function(X) getGCi(X, CDS_positions, position = 2), mc.cores = ncores)
      gc4genes$gc2 = unlist(gc4genes$gc2)
      cat("GC3... ")
      gc4genes$gc3[idx.genes] = pbmclapply(idx.genes, function(X) getGCi(X, CDS_positions, position = 3), mc.cores = ncores)
      gc4genes$gc3 = unlist(gc4genes$gc3)
      
      # END OF LOOP ON CHROMOSOMES
    }
    
    # File is saved in a txt data frame and gunzipped for efficient storage management
    # See https://rce-docs.hmdc.harvard.edu/faq/how-do-i-use-compressed-data-r
    cat("Writing file...\n")
    write.table(gc4genes, file = gzfile(paste("data-cleaned/genome/gc_genes/", map, ".txt.gz", sep = "")),
                quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
    
    # END OF METHOD 'GENE'
  }
  #-------------------------------------------
  if (method == "exon") {
    # Init empty result dataframe
    gc4exons = data.frame(species = exon_positions$species, accession = exon_positions$accession, id = exon_positions$id,
                          chromosome = exon_positions$chromosome, rank = exon_positions$rank, rec.rate = numeric(nrow(exon_positions)),
                          gc = numeric(nrow(exon_positions)), gc1 = numeric(nrow(exon_positions)), gc2 = numeric(nrow(exon_positions)),
                          gc3 = numeric(nrow(exon_positions)))
    # Remove any line that is not in the chromosome list
    listchr = gsub("_chromosome", "", regmatches(list, gregexpr("_chromosome[A-Za-z0-9[:punct:]]+$", list)))
    gc4exons = gc4exons[which((gc4exons$chromosome %in% listchr)),]
    
    # LOOP ON EACH CHROMOSOME
    for (chromosome in list) {
      # chromosome = list[1] # DEBUG
      # Load the recombination map
      mapdata = read.table(file = paste("output/recombination_maps/loess/100kbwind/", chromosome, ".txt", sep = ""), header = TRUE, sep = "\t")
      # Work only on the current chromosome
      chrname = gsub("_chromosome", "", regmatches(chromosome, gregexpr("_chromosome[A-Za-z0-9[:punct:]]+$", chromosome)))
      # Retrieve the sequence for the given chromosome in the list
      chrinrecmaps = gsub("chromosome", "", regmatches(chromosome, regexpr("chromosome[A-Za-z]*[0-9]*[A-Za-z]*", chromosome)))
      cat("Chromosome ", chrinrecmaps, "\n")
      idx = which(chrinseq == chrinrecmaps)
      # Select the whole chromosome sequence for further analyses
      # In case of multiple fasta file, select arbitrary the first one
      seq_wholechr = seq[[idx[1]]]
      
      # LOOP ON EXONS Compute statistics for each position
      # Index of exons to process
      # This loop can easily be parallelisable
      # since estimates are independant, yet iterator and progress bar is not possible anymore
      idx.exons = which(exon_positions$chromosome == chrname)
      
      cat("Recovering Recombination rates... ")
      gc4exons$rec.rate[idx.exons] = pbmclapply(idx.exons, function(X) getRecRate(X, exon_positions), mc.cores = ncores)
      gc4exons$rec.rate = unlist(gc4exons$rec.rate)
      cat("Processing GC... ")
      gc4exons$gc[idx.exons] = pbmclapply(idx.exons, function(X) getGC(X, exon_positions), mc.cores = ncores)
      gc4exons$gc = unlist(gc4exons$gc)
      cat("GC1... ")
      gc4exons$gc1[idx.exons] = pbmclapply(idx.exons, function(X) getGCi(X, exon_positions, position = 1), mc.cores = ncores)
      gc4exons$gc1 = unlist(gc4exons$gc1)
      cat("GC2... ")
      gc4exons$gc2[idx.exons] = pbmclapply(idx.exons, function(X) getGCi(X, exon_positions, position = 2), mc.cores = ncores)
      gc4exons$gc2 = unlist(gc4exons$gc2)
      cat("GC3... ")
      gc4exons$gc3[idx.exons] = pbmclapply(idx.exons, function(X) getGCi(X, exon_positions, position = 3), mc.cores = ncores)
      gc4exons$gc3 = unlist(gc4exons$gc3)
                  # # Foreach option not working properly, mclapply option is cleaner
                  # foreach (g=idx.exons) %dopar% {
                  #      # count = count + 1 # Iterator for progress bar
                  #       # setTxtProgressBar(pb, count)
                  #   cat(g/length(idx.exons), "\n")
                  #   # Query the local recombination rate of a gene by its physical position
                  #   # i.e. windows in which start position is
                  #   # And in any case, prevent for NA estimates out of boundaries
                  #   if ((exon_positions$start[g] < min((mapdata$phys*1000000))) | (exon_positions$start[g] > max((mapdata$phys*1000000)))) {
                  #     gc4exons$rec.rate[g] = NA
                  #   } else {
                  #     gc4exons$rec.rate[g] = mapdata$rec.rate[max(which(exon_positions$start[g] > (mapdata$phys*1000000)))]
                  #   }
                  #   
                  #   # Next estimates GC on the whole gene sequence
                  #   # First retrieve the sequence with gene coordinates
                  #   # MEAN GC, sequence of the complete gene sequence
                  #   seqwin = seq_wholechr[exon_positions$start[g]:min(c(exon_positions$end[g], length(seq_wholechr)))]
                  #   # NAs produced in tails of the sequence
                  #   seqwin[is.na(seqwin)] = "N"
                  #   seqwin = s2c(tolower(as.character(seqwin))) # Convert to a vector of single characters
                  #   # Compute local GC proportion along the complete gene sequence
                  #   if (length(seqwin) > 0) {
                  #     # forcetoLower = FALSE require to convert sequence to lower case before, yet it reduces the computation time by 3
                  #     # exact computation estimates GC content with ambiguous bases taken in account. Time is increased by 3.
                  #     # Yet, it revealed differences between estimations with exact = FALSE
                  #     gc4exons$gc[g] = GC(seqwin, forceToLower = FALSE, exact = ex)
                  #     # distGC = c(distGC, GC(seqwin, forceToLower = FALSE, exact = ex))
                  #   } else {
                  #     gc4exons$gc[g] = NA
                  #   }
                  #   
                  #   # GC IN CODONS
                  #   # Trim sequences to CDS contained in map windows in a sequence object
                  #   seqCDS = data.frame(start = CDS_positions$start[CDS_positions$accession == accession & CDS_positions$chromosome == chrinrecmaps],
                  #                       end = CDS_positions$end[CDS_positions$accession == accession & CDS_positions$chromosome == chrinrecmaps],
                  #                       strand = CDS_positions$strand[CDS_positions$accession == accession & CDS_positions$chromosome == chrinrecmaps],
                  #                       phase = CDS_positions$phase[CDS_positions$accession == accession & CDS_positions$chromosome == chrinrecmaps])
                  #   
                  #   # Retain only CDS associated to the exon
                  #   windowed_CDS = seqCDS[(seqCDS$start >= exon_positions$start[g]) & (seqCDS$start <= exon_positions$end[g]),]
                  #   # Init vector of results: sum of GCi proportions for each sequence
                  #   mean_gc1 = numeric(nrow(windowed_CDS))
                  #   mean_gc2 = numeric(nrow(windowed_CDS))
                  #   mean_gc3 = numeric(nrow(windowed_CDS))
                  #   
                  #   # cat("GC proportions in CDS")
                  #   # If some CDS have been selected in the current window...
                  #   if (nrow(windowed_CDS) > 0) {
                  #     # LOOP ON CDS
                  #     # For each CDS within the window
                  #     for (c in 1:nrow(windowed_CDS)) {
                  #       # print(c)
                  #       # Retrieve the fasta sequence of the window; respects bounds of the sequence at the end
                  #       # If subset is out of bounds of reference genome sequence, then add NA values
                  #       if (exon_positions$start[g] < length(seq_wholechr)) {
                  #         seqwin = seq_wholechr[windowed_CDS$start[c]:min(c(windowed_CDS$end[c], length(seq_wholechr)))]
                  #         seqwin[is.na(seqwin)] = "N"
                  #         seqwin = s2c(tolower(as.character(seqwin))) # Convert to a vector of single characters in lower case
                  #         # IF reverse strand (strand = -)
                  #         if (windowed_CDS$strand[c] == "-") {
                  #           seqwin = rev(seqwin)
                  #         }
                  #         if ((length(seqwin) - windowed_CDS$phase[c]) >= 3) {
                  #           mean_gc1[c] = GC1(seqwin, frame = windowed_CDS$phase[c], forceToLower = FALSE, exact = ex)
                  #           mean_gc2[c] = GC2(seqwin, frame = windowed_CDS$phase[c], forceToLower = FALSE, exact = ex)
                  #           mean_gc3[c] = GC3(seqwin, frame = windowed_CDS$phase[c], forceToLower = FALSE, exact = ex)
                  #           # mean_gc1 = c(mean_gc1, GC1(seqwin, forceToLower = FALSE), exact = TRUE)
                  #           # mean_gc2 = c(mean_gc2, GC2(seqwin, forceToLower = FALSE), exact = TRUE)
                  #           # mean_gc3 = c(mean_gc3, GC3(seqwin, forceToLower = FALSE), exact = TRUE)
                  #         }
                  #       } else {
                  #         mean_gc1[c] = NA
                  #         mean_gc2[c] = NA
                  #         mean_gc3[c] = NA
                  #       }
                  #       # END OF CDS LOOP
                  #     }
                  #     
                  #     # Compute mean GC proportions according to the position in codons on a list of CDS sequences
                  #     # mean(mean_gc1, na.rm = TRUE)
                  #     # mean(mean_gc2, na.rm = TRUE)
                  #     # mean(mean_gc3, na.rm = TRUE)
                  #     gc4exons$gc1[g] = mean(mean_gc1, na.rm = TRUE)
                  #     gc4exons$gc2[g] = mean(mean_gc2, na.rm = TRUE)
                  #     gc4exons$gc3[g] = mean(mean_gc3, na.rm = TRUE)
                  #     # distGC3 = c(distGC3, mean_gc3)
                  #   } else { # Sometimes there is no CDS recovered within the window
                  #     gc4exons$gc1[g] = NA
                  #     gc4exons$gc2[g] = NA
                  #     gc4exons$gc3[g] = NA
                  #   }
                  #   rm(mean_gc1, mean_gc2, mean_gc3)
                  #   
                  #   # END OF LOOP ON EACH POSITION
                  # }
      
      # END OF LOOP ON CHROMOSOMES
    }
    
    # File is saved in a txt data frame and gunzipped for efficient storage management
    # See https://rce-docs.hmdc.harvard.edu/faq/how-do-i-use-compressed-data-r
    cat("Writing file...\n")
    write.table(gc4exons, file = gzfile(paste("data-cleaned/genome/gc_exons/", map, ".txt.gz", sep = "")),
                quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
    
    
    
    # END OF METHOD 'EXON'  
  }
  

  return("Done.")
  # END OF FUNCTION
}

#------------------
# EXAMPLE
# map = "Arabidopsis_thaliana_Serin2017"
# gc.map(map)
# gcmap = read.table(file = paste("data-cleaned/genome/gc_maps/", map, "_chromosome4.txt", sep = ""), header = TRUE, sep = "\t")
# plot(gcmap$gc)
# plot(gcmap$gc1)
# plot(gcmap$gc2)
# plot(gcmap$gc3)





#============================================================================
#   Gene count in recombination windows
#============================================================================
# Count the number of genes in a window of the recombination map (i.e. genes with a start position in the window)

gene_count = function(set = "", scale = 0.1) {

  cat("Processing ", as.character(set), "\n")
  # Load the db for the species, given a species name and accession number
  metadata.clean = read.table("data-cleaned/Genetic_maps_ressources.csv", header = TRUE, sep = ";")
  species = as.character(metadata.clean$species[grep(set, metadata.clean$id)[1]])
  accession = as.character(metadata.clean$accession[grep(set, metadata.clean$id)[1]])
  # Import list of genes, only for method = 'gene'
  gene_positions = read.table(gzfile(paste("data-cleaned/genome/gene_positions/", species,"_gene_positions.txt.gz", sep = "")), header = TRUE, sep = "\t")
 
  # For each chromosome for which a recombination map exists
  list_chromosomes = system(paste("ls ", wd, "/output/recombination_maps/loess/100kbwind/", sep = ""), intern = TRUE)
  list_chromosomes = list_chromosomes[grep(set, list_chromosomes)]
  list_chromosomes = gsub(".txt", "", list_chromosomes)

  for (chromosome in list_chromosomes) {
    # chromosome = list_chromosomes[1] # Debug purpose only

    # Load the recombination map
    genecount_map = read.table(file = paste("output/recombination_maps/loess/100kbwind/", chromosome, ".txt", sep = ""), header = TRUE, sep = "\t")
    # Initialise the columns for new statistics to compute
    genecount_map$gene_count = NA # Count of genes in each windows
    # genedensity_map$gene_distance = NA # Cumulative number of genes, for genetic distance = 0 cM
    
    # A function that compute the gene density in a window
    # Take the line the map data and the gene_positions as arguments
    gene_count = function(data, genepos = gene_positions) {
      # print(data)
      # Estimate the window boundaries from the scale parameter, and multiply per 6 to convert Mb to bp
      # Gene density is the number of gene start positions in the window
      # data$gene_density = sum(gene_positions$start > (data$phys - scale/2)*1000000 & gene_positions$start < (data$phys + scale/2)*1000000)
      return(sum((genepos$start > as.numeric((data[1] - scale/2)*1000000)) & (genepos$start < as.numeric((data[1] + scale/2)*1000000)) & (genepos$chromosome == gsub(paste(set, "_chromosome", sep = ""), "", chromosome)), na.rm = TRUE))
      
    }
    
    # A LOOP ON EACH WINDOW
    # For each 100kb window in the map
    genecount_map$gene_count = apply(genecount_map, 1, function(x) gene_count(x, gene_positions))
    # END OF WINDOW LOOP  
    
    # # A function that compute the gene distance of the window
    # # It is the cumulative distance of all previous windows
    # # Take a line of map data and the genedensity_map as arguments
    # gene_distance = function(data, genedensity_map) {
    #   # Gene distance is the sum of all gene count in current and previous windows
    #   return(sum(genedensity_map$gene_density[genedensity_map$phys <= data[1]]))
    # }
    # 
    # # A LOOP ON EACH WINDOW
    # # For each 100kb window in the map
    # genedensity_map$gene_distance = apply(genedensity_map, 1, function(x) gene_distance(x, genedensity_map))
    # # END OF WINDOW LOOP  
    
    # File is saved in a txt data frame and gunzipped for efficient storage management
    # See https://rce-docs.hmdc.harvard.edu/faq/how-do-i-use-compressed-data-r
    cat("Writing file...\n")
    write.table(genecount_map, file = gzfile(paste("data-cleaned/genome/gene_count/", chromosome, ".txt.gz", sep = "")),
                  quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
      
    # END OF CHROMOSOME LOOP
    }
  
  return("Done.")
  # END OF FUNCTION
}


#============================================================================
#   Gene density in recombination windows
#============================================================================
# Estimate the proportion of gene sequence and CDS in each recombination window

gene_density = function(set = "", scale = 0.1) {
  cat("Processing ", as.character(set), "\n")
  # Load the db for the species, given a species name and accession number
  metadata.clean = read.table("data-cleaned/Genetic_maps_ressources.csv", header = TRUE, sep = ";")
  species = as.character(metadata.clean$species[grep(set, metadata.clean$id)[1]])
  accession = as.character(metadata.clean$accession[grep(set, metadata.clean$id)[1]])
  # Import list of genes
  gene_positions = read.table(gzfile(paste("data-cleaned/genome/gene_positions/", species,"_gene_positions.txt.gz", sep = "")), header = TRUE, sep = "\t")
  cds_positions = read.table(gzfile(paste("data-cleaned/genome/CDS_positions/", species,"_CDS_positions.txt.gz", sep = "")), header = TRUE, sep = "\t")
  
  # For each chromosome for which a recombination map exists
  list_chromosomes = system(paste("ls ", wd, "/output/recombination_maps/loess/100kbwind/", sep = ""), intern = TRUE)
  list_chromosomes = list_chromosomes[grep(set, list_chromosomes)]
  list_chromosomes = gsub(".txt", "", list_chromosomes)
  
  # A function that compute the gene density in a window
  # Take the line the map data and the gene_positions as arguments
  gene_density_estimate = function(data, genepos = gene_positions) {
    # Estimate the window boundaries from the scale parameter, and multiply per 6 to convert Mb to bp
    # Gene density is the proportion of the window that is within a gene
    # Subset genes within the window
    # Boundaries of the window
    lower = as.numeric(data[1] - scale/2)*1000000
    upper = as.numeric(data[1] + scale/2)*1000000
    if (is.na(lower) | is.na(upper)) {
      return(NA)
    }
    # i.e. genes on the chromosome and either start or end position is within the window
    subset = genepos[which(((genepos$start > lower & genepos$start < upper) |
                              (genepos$end > lower & genepos$end < upper))
                           & (genepos$chromosome == gsub(paste(set, "_chromosome", sep = ""), "", chromosome))),]
    # Cut extra-window sequences
    subset$start[subset$start < lower] = lower
    subset$end[subset$end > upper] = upper
    # We need to deal with overlapping genes
    # Avoid counting twice the same nucleotide
    # Merge genes if one of their boundaries is overlapping
    library(valr)
    colnames(subset)[5] = "chrom" # Group by chromosomes
    mergedsequences = as.data.frame(bed_merge(subset))
    # Estimate the proportion of gene sequence
    return(sum(mergedsequences$end - mergedsequences$start, na.rm = TRUE)/(1000000*scale))
  }
  # A function that compute the coding sequence density in a window
  # Take the line the map data and the cds_positions as arguments
  coding_density_estimate = function(data, cdspos = cds_positions) {
    # Estimate the window boundaries from the scale parameter, and multiply per 6 to convert Mb to bp
    # Gene density is the proportion of the window that is within a gene
    # Subset CDS within the window
    # Only take the first splicing variant
    # cdspos$id[grep(".[2-20]$", cdspos$id)]
    cdspos = cdspos[which(!grepl(".[2-9]$", cdspos$id)),]
    # Boundaries of the window
    lower = as.numeric(data[1] - scale/2)*1000000
    upper = as.numeric(data[1] + scale/2)*1000000
    if (is.na(lower) | is.na(upper)) {
      return(NA)
    }
    # i.e. genes on the chromosome and either start or end position is within the window
    subset = cdspos[which(((cdspos$start > lower & cdspos$start < upper) |
                             (cdspos$end > lower & cdspos$end < upper))
                          & (cdspos$chromosome == gsub(paste(set, "_chromosome", sep = ""), "", chromosome))),]
    # Cut extra-window sequences
    subset$start[subset$start < lower] = lower
    subset$end[subset$end > upper] = upper
    # We need to deal with overlapping genes
    # Avoid counting twice the same nucleotide
    # Merge genes if one of their boundaries is overlapping
    library(valr)
    colnames(subset)[5] = "chrom" # Group by chromosomes
    mergedsequences = as.data.frame(bed_merge(subset))
    # Estimate the proportion of gene sequence
    return(sum(mergedsequences$end - mergedsequences$start, na.rm = TRUE)/(1000000*scale))
  }
  
  
  # Parallelize chromosomes
  chromosome_genedensity = function(chromosome) {
    library(parallel)
    # library(foreach)
    # create cluster
    clust = detectCores()-1
    # Load the recombination map
    genedensity_map = read.table(file = paste("output/recombination_maps/loess/100kbwind/", chromosome, ".txt", sep = ""), header = TRUE, sep = "\t")
    # Initialise the columns for new statistics to compute
    genedensity_map$gene_density = NA # Proportion of gene sequence in each windows
    genedensity_map$coding_density = NA # Proportion of CDS in each windows
    # A LOOP ON EACH WINDOW
    # For each 100kb window in the map
    # genedensity_map$gene_density = unlist(foreach(i = 1:nrow(genedensity_map)) %do% gene_density_estimate(genedensity_map[i,], cds_positions))
    # genedensity_map$gene_density = parApply(cl = clust, genedensity_map[1:100,], 1, function(x) gene_density_estimate(x, gene_positions))
    genedensity_map$gene_density = unlist(mclapply(as.list(c(1:nrow(genedensity_map))), function(x) gene_density_estimate(genedensity_map[x,], gene_positions), mc.cores = clust))
    # END OF WINDOW LOOP
    
    # A LOOP ON EACH WINDOW
    # For each 100kb window in the map
    # genedensity_map$coding_density = unlist(foreach(i = 1:nrow(genedensity_map)) %do% coding_density_estimate(genedensity_map[i,], cds_positions))
    # genedensity_map$coding_density = parApply(cl = clust, genedensity_map, 1, function(x) coding_density_estimate(x, cds_positions), mc.cores = 7)
    genedensity_map$coding_density = unlist(mclapply(as.list(c(1:nrow(genedensity_map))), function(x) coding_density_estimate(genedensity_map[x,], cds_positions), mc.cores = clust))
    # END OF WINDOW LOOP 
    
    # File is saved in a txt data frame and gunzipped for efficient storage management
    # See https://rce-docs.hmdc.harvard.edu/faq/how-do-i-use-compressed-data-r
    cat("Writing file...", chromosome, "\n")
    write.table(genedensity_map, file = gzfile(paste("data-cleaned/genome/gene_density/", chromosome, ".txt.gz", sep = "")),
                quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
    
    # END OF CHROMOSOME LOOP
  }
  
  for (chromosome in list_chromosomes) {
    # cat(chromosome, "\n")
    chromosome_genedensity(chromosome)
  }
  # lapply(list_chromosomes, function(x) chromosome_genedensity(x))
  
  # for (chromosome in list_chromosomes) {
  #   # chromosome = list_chromosomes[1] # Debug purpose only
  #   
  #   # Load the recombination map
  #   genedensity_map = read.table(file = paste("output/recombination_maps/loess/100kbwind/", chromosome, ".txt", sep = ""), header = TRUE, sep = "\t")
  #   # Initialise the columns for new statistics to compute
  #   genedensity_map$gene_density = NA # Proportion of gene sequence in each windows
  #   genedensity_map$coding_density = NA # Proportion of CDS in each windows
  #   
  #   # A function that compute the gene density in a window
  #   # Take the line the map data and the gene_positions as arguments
  #   gene_density_estimate = function(data, genepos = gene_positions) {
  #     # Estimate the window boundaries from the scale parameter, and multiply per 6 to convert Mb to bp
  #     # Gene density is the proportion of the window that is within a gene
  #     # Subset genes within the window
  #     # Boundaries of the window
  #     lower = as.numeric(data[1] - scale/2)*1000000
  #     upper = as.numeric(data[1] + scale/2)*1000000
  #     if (is.na(lower) | is.na(upper)) {
  #       return(NA)
  #     }
  #     # i.e. genes on the chromosome and either start or end position is within the window
  #     subset = genepos[which(((genepos$start > lower & genepos$start < upper) |
  #                               (genepos$end > lower & genepos$end < upper))
  #                             & (genepos$chromosome == gsub(paste(set, "_chromosome", sep = ""), "", chromosome))),]
  #     # Cut extra-window sequences
  #     subset$start[subset$start < lower] = lower
  #     subset$end[subset$end > upper] = upper
  #     # We need to deal with overlapping genes
  #     # Avoid counting twice the same nucleotide
  #     # Merge genes if one of their boundaries is overlapping
  #     library(valr)
  #     colnames(subset)[5] = "chrom" # Group by chromosomes
  #     mergedsequences = as.data.frame(bed_merge(subset))
  #     # Estimate the proportion of gene sequence
  #     return(sum(mergedsequences$end - mergedsequences$start, na.rm = TRUE)/(1000000*scale))
  #   }
  #   # A LOOP ON EACH WINDOW
  #   # For each 100kb window in the map
  #   genedensity_map$gene_density = apply(genedensity_map, 1, function(x) gene_density_estimate(x, gene_positions))
  #   # END OF WINDOW LOOP
  #   
  #   # A function that compute the coding sequence density in a window
  #   # Take the line the map data and the cds_positions as arguments
  #   coding_density_estimate = function(data, cdspos = cds_positions) {
  #     # Estimate the window boundaries from the scale parameter, and multiply per 6 to convert Mb to bp
  #     # Gene density is the proportion of the window that is within a gene
  #     # Subset CDS within the window
  #     # Only take the first splicing variant
  #     # cdspos$id[grep(".[2-20]$", cdspos$id)]
  #     cdspos = cdspos[which(!grepl(".[2-9]$", cdspos$id)),]
  #     # Boundaries of the window
  #     lower = as.numeric(data[1] - scale/2)*1000000
  #     upper = as.numeric(data[1] + scale/2)*1000000
  #     if (is.na(lower) | is.na(upper)) {
  #       return(NA)
  #     }
  #     # i.e. genes on the chromosome and either start or end position is within the window
  #     subset = cdspos[which(((cdspos$start > lower & cdspos$start < upper) |
  #                               (cdspos$end > lower & cdspos$end < upper))
  #                            & (cdspos$chromosome == gsub(paste(set, "_chromosome", sep = ""), "", chromosome))),]
  #     # Cut extra-window sequences
  #     subset$start[subset$start < lower] = lower
  #     subset$end[subset$end > upper] = upper
  #     # We need to deal with overlapping genes
  #     # Avoid counting twice the same nucleotide
  #     # Merge genes if one of their boundaries is overlapping
  #     library(valr)
  #     colnames(subset)[5] = "chrom" # Group by chromosomes
  #     mergedsequences = as.data.frame(bed_merge(subset))
  #     # Estimate the proportion of gene sequence
  #     return(sum(mergedsequences$end - mergedsequences$start, na.rm = TRUE)/(1000000*scale))
  #   }
  #   
  #   # A LOOP ON EACH WINDOW
  #   # For each 100kb window in the map
  #   genedensity_map$coding_density = apply(genedensity_map, 1, function(x) coding_density_estimate(x, cds_positions))
  #   # END OF WINDOW LOOP 
  # 
  #   # File is saved in a txt data frame and gunzipped for efficient storage management
  #   # See https://rce-docs.hmdc.harvard.edu/faq/how-do-i-use-compressed-data-r
  #   # cat("Writing file...\n")
  #   write.table(genedensity_map, file = gzfile(paste("data-cleaned/genome/gene_density/", chromosome, ".txt.gz", sep = "")),
  #               quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
  #   
  #   # END OF CHROMOSOME LOOP
  # }
  
  return("Done.")
  # END OF FUNCTION
}



#============================================================================
#   Gene distances in Marey maps
#============================================================================
# Method 'map'
# Count the cumulative number of genes in a window of the recombination map (i.e. distance on gene number)

# USAGE
# gene_distance.map(map = "", scale = 0.1)

# ARGUMENTS
# map = the recombination map on which you want to estimate the local GC proportions
# scale in Mb (e.g. 0.1 for 100kb windows)

# VALUES
# Return the exact recombination map + statistics computed on the genome
# Gene distance

gene_distance.map = function(species = "", scale = 0.1) {
  library(pbmcapply)
  cat("=================================================\nProcessing ", as.character(species), "\n")
  # A function that compute the gene density in a window
  # Take the line the map data and the gene_positions as arguments
  gene_distance = function(data, genepos = gene_positions) {
    # print(data)
    # data = map[10,]
    # genepos$chromosome = as.character(genepos$chromosome)
    # Estimate the window boundaries from the scale parameter, and multiply per 6 to convert Mb to bp
    # Gene density is the number of gene start positions in the window
    # print(data[2])
    # print(as.character(data[2]) %in% as.character(genepos$chromosome))
    
    # print(sum((as.character(genepos$chromosome) == as.character(data[2])) & (genepos$start < as.numeric(data[4]))))
    return(sum((as.character(gene_positions$chromosome) == as.character(data$map)) & (gene_positions$start < as.numeric(data$phys)),
               na.rm = TRUE))
  }
  # DONE Bugfix: apply() works only with the last chromosome (tenth) FIXED
  
  
  # Load the db for the species, given a species name and accession number
  metadata.clean = read.table("data-cleaned/Genetic_maps_ressources.csv", header = TRUE, sep = ";",
                              stringsAsFactors = FALSE)
  # Many dataset can exists for a species
  list_dataset = as.character(metadata.clean$id[which(metadata.clean$species == species)])
  
  # Treat independantly each dataset
  for (s in list_dataset) {
    cat("Dataset ", as.character(s), "\n")
    gene_positions = read.table(gzfile(paste("data-cleaned/genome/gene_positions/", species,"_gene_positions.txt.gz",
                                             sep = "")), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    gene_positions$chromosome = as.character(gene_positions$chromosome)
    
    map = read.table(file = paste("data-cleaned/marey_maps/", s, ".txt", sep = ""), header = TRUE, sep = "\t")
    map$gene_distance = NA
    map$map = as.character(map$map)
    
    # A LOOP ON EACH WINDOW
    # For each 100kb window in the map
    map$gene_distance = unlist(pbmclapply(c(1:nrow(map)), function(x) gene_distance(map[x,], gene_positions)))
    # END OF WINDOW LOOP
    
    # File is saved in a txt data frame and gunzipped for efficient storage management
    # See https://rce-docs.hmdc.harvard.edu/faq/how-do-i-use-compressed-data-r
    cat("Writing file...\n")
    write.table(map, file = gzfile(paste("data-cleaned/genome/gene_distance/", s, ".txt.gz", sep = "")),
                quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
  }
  return("Done.")
  # END OF FUNCTION
}


#============================================================================
#   Number of exons per genes and associated statistics
#============================================================================
# A function that take map name and return & save a file containing the list of genes and the number of exons associated
# Plus additional informations such as recombination and GC3 (gene) + GC3 exon 1

genesize_nbexons = function(map) {
  # Load data
  gcexons = read.table(file = gzfile(paste("data-cleaned/genome/gc_exons/", map, ".txt.gz", sep = "")), header = TRUE, sep = "\t")
  
  # Remove alternative splicing variants for a subset of species with this information
  if (map %in% c("Arabidopsis_thaliana_Serin2017")) {
    # Index of splicing variants number 1
    idx = which(grepl(pattern = "\\.1\\.", gcexons$id))
  }
  if (map %in% c("Oryza_sativa_DeLeon2016")) {
    # Index of splicing variants number 1
    idx = which(grepl(pattern = "-01-", gcexons$id))
  }
  if (map %in% c("Zea_mays_MaizeGDBConsensus_v4", "Zea_mays_IBM_MaizeSNP50")) {
    # Index of splicing variants number 1
    idx = which(grepl(pattern = "_T001", gcexons$id))
  }
  # Remove variants from data
  gcexons = gcexons[idx,]
  
  # Make the list of genes
  # Function take the list of exons and save the list of genes in a new column, geneid
  genes_from_exons = function(data) {
    # Management of exceptions
    # Each genome dataset has its own pattern of naming exons
    if (map %in% c("Malus_domestica_DiPierro2016")) {
      pattern = "-RA-[0-9]+"
    }
    if (map %in% c("Arabidopsis_thaliana_Serin2017")) {
      pattern = "\\.[0-9]*\\.[A-Za-z0-9]*"
    }
    if (map %in% c("Oryza_sativa_DeLeon2016")) {
      pattern = "-[0-9]*-[A-Za-z0-9]*"
    }
    if (map %in% c("Zea_mays_MaizeGDBConsensus_v4", "Zea_mays_IBM_MaizeSNP50")) {
      pattern = "_T001.exon[A-Za-z0-9]*"
    }
    if (map %in% c("Brachypodium_distachyon_Huo2011", "Cucumis_sativus_Zhu2016", "Gossypium_raimondii_Wang2013",
                   "Phaseolus_vulgaris_Song2015", "Sorghum_bicolor_Zou2012", "Setaria_italica_Ni2017", 
                   "Triticum_aestivum_GutierrezGonzalez2019")) {
      pattern = "-[A-Za-z0-9]*$"
    }
    if (map %in% c("Solanum_lycopersicum_Gonda2018")) {
      pattern = "\\.[0-9]*\\.[0-9]*-E[0-9]*$"
    }
    
    # remove the pattern from exon id
    data$geneid = gsub(pattern, "", data$id)
    return(data)
  }
  
  gcexons = genes_from_exons(gcexons)
  
  # The list of genes (unique values in geneid)
  list_genes = unique(gcexons$geneid)
  list_genes = list_genes[!is.na(list_genes)]
  
  # For each gene retrieve the max rank of exons, i.e. the number of exons
  # For each gene, the splicing variant chosen is systematically the version 1
  # Take the name of a gene in argument (x) and return the number of exons
  cat("Processing number of exons")
  number_exons = function(x) {
    return(max(gcexons$rank[which(gcexons$geneid == x)]))
  }
  nb_exons = unlist(pbmclapply(list_genes, number_exons))
  
  # Assembling results in a data frame
  # Recombination rate of a gene is the mean of recombination rates of all exons
  # GC3 of a gene is the mean of GC3 of all exons
  # GC3ex1 is simply the GC3 of the first exon (rank == 1)
  cat("Processing recombination rate of genes")
  recombination_rate = function(x) {
    return(mean(gcexons$rec.rate[which(gcexons$geneid == x)], na.rm = TRUE))
  }
  rec.rate = unlist(pbmclapply(list_genes, recombination_rate))
  
  cat("Processing GC3 of genes")
  est_GC3 = function(x) {
    return(mean(gcexons$gc3[which(gcexons$geneid == x)], na.rm = TRUE))
  }
  GC3 = unlist(pbmclapply(list_genes, est_GC3))
  
  cat("Processing GC3 first exon")
  est_GC3ex1 = function(x) {
    if (length(which(gcexons$geneid == x & gcexons$rank == 1)) == 1) {
      ret = gcexons$gc3[which(gcexons$geneid == x & gcexons$rank == 1)]
    }
    else {
      ret = NA
    }
    return(ret)
  }
  GC3ex1 = unlist(pbmclapply(list_genes, est_GC3ex1))
  
  
  # Chromosome for genes
  get_chrname = function(x) {
    if (length(unique(gcexons$chromosome[which(gcexons$geneid == x)])) == 0) {
      ret = NA
    } else {ret = unique(gcexons$chromosome[which(gcexons$geneid == x)])}
    return(ret)
  }
  chr = unlist(pbmclapply(list_genes, get_chrname))
  
  cat("Saving...") 
  res = data.frame(gene_id = list_genes, chromosome = chr, number_exons = nb_exons, rec.rate = rec.rate, GC3 = GC3, GC3ex1 = GC3ex1)
  # Save the result
  write.table(res, file = gzfile(paste("data-cleaned/genome/number_exons/", map, "_exonnumber.txt.gz", sep = "")),
              quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
  return(res)
}