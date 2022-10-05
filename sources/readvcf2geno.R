
# Start of vcf.R ###############################################################

# getSamplesVCF ----------------------------------------------------------------
#' Get sample IDs from VCF file.
#' 
#' @param file A VCF file path.
#' 
#' @return Vector of the sample IDs in the given VCF file.
#'   
#' @export
#' @family VCF functions
#' @rdname getSamplesVCF
getSamplesVCF <- function(file) {
  stopifnot( isSingleString(file) )
  stopifnot( file.exists(file) )
  header <- VariantAnnotation::scanVcfHeader(file)
  return( VariantAnnotation::samples(header) )
}

# hasGenoQualVCF ---------------------------------------------------------------
#' Test if VCF file contains GQ scores.
#' 
#' @param file A VCF file path.
#' 
#' @return \code{TRUE} if the specified VCF file contains
#' genotype quality (GQ) scores; \code{FALSE} otherwise.
#'   
#' @keywords internal
#' @rdname hasGenoQualVCF
hasGenoQualVCF <- function(file) {
  stopifnot( isSingleString(file) )
  stopifnot( all( file.exists(file) ) )
  header <- VariantAnnotation::scanVcfHeader(file)
  return( 'GQ' %in% rownames( VariantAnnotation::geno(header) ) )
}

# readGenoVCF ------------------------------------------------------------------
#' Read genotype data from a VCF file.
#' 
#' This function reads SNP genotype data from one or more VCF files, and
#' returns these as an \pkg{R/qtl} \code{cross} \code{geno} object.
#' 
#' If no founder samples are specified, this function assigns enumerated
#' genotypes at each locus according to the observed raw SNP genotype. So
#' for example, if the raw SNP alleles at a given locus are \code{'A'} and
#' \code{'C'}, samples are assigned the genotypes \code{'1'} and \code{'2'},
#' respectively.
#' 
#' If founder samples are specified, this function assigns to each marker a
#' genotype symbol consisting of alleles that each correspond to a specific
#' founder.
#' 
#' If the \code{alleles} parameter is specified, this must be a mapping of
#' founder sample IDs to allele symbols (e.g.
#' \code{mapping( c(DBVPG6044 = 'W', Y12 = 'S') )}). If the \code{alleles}
#' parameter is not specified, allele symbols are taken from the letters of the
#' alphabet (i.e. \code{'A'}, \code{'B'} etc.).
#' 
#' @param ... Input VCF file paths.
#' @param samples Cross sample IDs.
#' @param founders Founder sample IDs.
#' @param alleles Mapping of founder sample IDs to founder allele symbols.
#' 
#' @return An \pkg{R/qtl} \code{cross} \code{geno} object.
#' 
#' @template section-geno-raw
#' @template section-geno-enum
#' @template section-geno-founder
#' 
#' @export
#' @family VCF functions
#' @rdname readGenoVCF
readGenoVCF <- function(..., samples, founders=NULL, alleles=NULL) {
  
  # TODO: optimise.
  
  sample.data <- readSnpsVCF(..., samples=samples, require.any=TRUE,
                             require.polymorphic=TRUE)
  
  if ( ! is.null(founders) ) {
    
    clashing <- intersect(samples, founders)
    if ( length(clashing) > 0 ) {
      stop("clashing sample IDs - ", toString(clashing), "'")
    }
    
    founder.data <- readSnpsVCF(..., samples=founders,
                                require.all=TRUE, require.polymorphic=TRUE)
    
  } else {
    
    founder.data <- NULL
  }
  
  geno <- makeGeno(sample.data, founder.data, alleles=alleles)
  
  return(geno)
}

# readSnpsVCF ------------------------------------------------------------------
#' Read raw SNP genotypes from VCF files.
#' 
#' This function reads SNP genotype data from one or more VCF files, and returns
#' these as a sequence of raw SNP genotypes - effectively variant base calls -
#' for each sample. If all relevant input VCF files contain genotype quality
#' (GQ) scores for the samples of interest, these are used - along with variant
#' quality scores - to calculate an error probability for each SNP genotype, and
#' the result is returned as an \code{array} with two slices - \code{'geno'} and
#' \code{'prob'} - containing raw SNP genotypes and their corresponding error
#' probabilities, respectively. Otherwise, SNP genotypes are returned as an
#' \code{array} with one slice, \code{'geno'}, which contains raw SNP genotype
#' values. In any case, the rows of each slice contain data for a given sample,
#' while the columns of each slice contain data for a given SNP locus.
#' 
#' @param ... Input VCF file paths.
#' @param samples Vector of samples for which SNP genotypes should be obtained.
#' If not specified, genotypes are returned for all available samples.
#' @param require.all Remove variants that are not completely genotyped
#' with respect to the given samples.
#' @param require.any Remove variants that do not have at least one genotype
#' call among the given samples.
#' @param require.polymorphic Remove variants that do not have at least two
#' different genotype calls among the given samples.
#' 
#' @return An \code{array} object containing raw genotype data, and if available,
#' genotype error probabilities.
#'   
#' @importFrom BiocGenerics start
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges CharacterList
#' @importFrom IRanges ranges
#' @keywords internal
#' @rdname readSnpsVCF
readSnpsVCF <- function(..., samples=NULL, require.all=FALSE, require.any=FALSE,
                        require.polymorphic=FALSE) {
  
  infiles <- unlist( list(...) )
  stopifnot( length(infiles) > 0 )
  stopifnot( isBOOL(require.all) )
  stopifnot( isBOOL(require.polymorphic) )
  
  # Set reference genome.
  # TODO: check contig info against reference genome.
  genome <- genomeOpt()
  
  # Get samples in each input VCF file.
  sample.list <- lapply(infiles, getSamplesVCF)
  
  # If samples specified, filter sample list, keep only those specified. 
  if ( ! is.null(samples) ) {
    
    stopifnot( is.character(samples) )
    stopifnot( length(samples) > 0 )
    
    for ( i in seq_along(infiles) ) {
      sample.list[[i]] <- sample.list[[i]][ sample.list[[i]] %in% samples ]
    }
    
    unfound <- samples[ ! samples %in% unlist(sample.list) ]
    if ( length(unfound) > 0 ) {
      stop("samples not found - '", toString(unfound), "'")
    }
  }
  
  samples <- unlist(sample.list)
  
  if ( length(samples) == 0 ) {
    stop("no samples found")
  }
  
  dup.samples <- samples[ duplicated(samples) ]
  if ( length(dup.samples) > 0 ) {
    stop("duplicate samples - '", toString(dup.samples), "'")
  }
  
  invalid.ids <- samples[ ! isValidID(samples) ]
  if ( length(invalid.ids) > 0 ) {
    stop("invalid sample IDs - '", toString(invalid.ids), "'")
  }
  
  # Keep only files relevant for the given samples.
  relevant <- lengths(sample.list) > 0
  sample.list <- sample.list[relevant]
  infiles <- infiles[relevant]
  
                # # Get mask indicating which relevant files contain genotype qualities.
                # gq.mask <- sapply(infiles, hasGenoQualVCF, USE.NAMES=FALSE)
                # 
                # # If all relevant files contain genotype qualities, 
                # # calculate SNP genotype quality scores..
                # if ( all(gq.mask) ) {
                #   reading.qualities <- TRUE
                #   slices <- c('geno', 'prob')
                # } else { # ..otherwise return only SNP genotypes.
                #   reading.qualities <- FALSE
                #   slices <- c('geno')
                #   if ( any(gq.mask) ) {
                #     warning("some VCF files lack GQ scores, reading genotypes only")
                #   }
                # }
  
  # Init list of SNP data.
  snps <- vector('list', length(infiles))
  
  # Init list for mapping genotype strings to allele indices.
  geno.memo <- list()
  
  # Read SNP genotypes from each input VCF file.
  for ( i in seq_along(infiles) ) {
    
    file.samples <- sample.list[[i]]
    infile <- infiles[[i]]
    
    # Set genotype fields to fetch.
    if (reading.qualities) {
      geno.fields <- c('GT', 'GQ')
    } else {
      geno.fields <- 'GT'
    }
    
    # Set VCF parameters.
    param <- VariantAnnotation::ScanVcfParam(fixed=c('ALT', 'QUAL'), 
                                             geno=geno.fields, samples=file.samples)
    
    # Read VCF.
    vcf <-  VariantAnnotation::readVcf(infile, genome=genome, param=param)
    stopifnot( length(vcf) > 0 )
    
    # Get variant data from VCF.
    var.seqs <- as.character( GenomeInfoDb::seqnames(vcf) ) # reference sequence
    var.pos <- BiocGenerics::start( IRanges::ranges(vcf) )  # position
    var.ref <- as.character( VariantAnnotation::ref(vcf) )  # reference allele
    var.alt <- lapply( IRanges::CharacterList(              # alternative alleles
      VariantAnnotation::alt(vcf) ), as.character)
    var.qual <- VariantAnnotation::qual(vcf)                # variant quality
    var.geno <- VariantAnnotation::geno(vcf)                # genotype info
    geno.matrix <- t( var.geno$GT )                         # genotype calls
    
    # Set indices of SNP variants.
    snp.indices <- which( VariantAnnotation::isSNV(vcf) )
    num.snps <- length(snp.indices)
    stopifnot( num.snps > 0 )
    
    # Filter genotype data to retain only SNP variants.
    geno.matrix <- geno.matrix[, snp.indices, drop=FALSE]
    
    # Combine REF and ALT alleles for each SNP.
    var.alleles <- lapply(snp.indices, function(j)
      c(var.ref[j], var.alt[[j]]))
    
    # Create dataframe with variant locus info.
    snp.loc <- data.frame(chr=var.seqs[snp.indices],
                          pos=var.pos[snp.indices])
    
    # Sanity check for ordered VCF records.
    if ( any( sapply( unique(var.seqs), function(var.seq)
      is.unsorted(snp.loc[snp.loc$chr == var.seq, 'pos']) ) ) ) {
      stop("unordered variants in file - '", infile, "'")
    }
    
    # Get default marker IDs for SNPs in this file.
    file.snps <- makeDefaultMarkerIDs( as.mapframe(snp.loc,
                                                   map.unit='bp') )
    
    # Check for multiple variants coinciding at same locus.
    # TODO: handle coinciding variants?
    if ( anyDuplicated(file.snps) ) {
      stop("coinciding variants in file - '", infile, "'")
    }
    
    colnames(geno.matrix) <- file.snps
    
    # Resolve raw genotypes for each variant as concatenated base calls.
    for ( j in 1:num.snps ) {
      
      # Get vector of alleles for this variant.
      variant.alleles <- var.alleles[[j]]
      
      # Get genotype strings for each sample in this variant record.
      geno.data <- unname( geno.matrix[, j] )
      
      # Get mask of genotypes that may have been called.
      # NB: this will misidentify diploid/polyploid missing values
      # as being genotyped, but such cases are handled below.
      called <- geno.data != const$vcf.missing.value
      
      # Set uncalled genotypes to missing value.
      geno.data[ ! called ] <- const$missing.value
      
      if ( any(called) ) {
        
        # Resolve previously unseen genotype strings
        # to their corresponding allele indices.
        for ( unique.call in unique(geno.data[called]) ) {
          if ( ! unique.call %in% names(geno.memo) ) {
            indicesC <- unlist( strsplit(unique.call, '[/|]') )
            indices0 <- suppressWarnings( as.integer(indicesC) )
            indices1 <- indices0 + 1
            geno.memo[[unique.call]] <- indices1
          }
        }
        
        # Get resolved allele indices for called genotypes.
        allele.list <- lapply(geno.data[called], function(geno.call)
          geno.memo[[geno.call]])
        
        # Check for consistent ploidy.
        ploidy <- unique( lengths(allele.list) )
        if ( length(ploidy) > 1 ) {
          k <- snp.indices[j]
          stop("mixed ploidy in variant at position ", var.pos[k],
               " of ", var.seqs[k], " in file - '", infile, "'")
        }
        
        # Create matrix of allele indices, set alleles (or missing values).
        allele.matrix <- matrix(unlist(allele.list), ncol=ploidy, byrow=TRUE)
        for ( a in seq_along(variant.alleles) ) {
          allele.matrix[ allele.matrix == a ] <- variant.alleles[a]
        }
        allele.matrix[ is.na(allele.matrix) ] <- const$missing.value
        
        # Concatenate alleles in each called genotype.
        geno.data[called] <- apply(allele.matrix, 1, paste0, collapse='')
      }
      
      geno.matrix[, j] <- geno.data
    }
    
    # Get symbols from genotype matrix.
    g.symbols <- unique( as.character(geno.matrix) )
    
    # Get genotype symbols.
    genotypes <- g.symbols[ g.symbols != const$missing.value ]
    
    # Get characters in genotype symbols.
    a.symbols <- unique( unlist( strsplit(genotypes, '') ) )
    
    # Get allele symbols.
    alleles <- a.symbols[ a.symbols != const$missing.value ]
    
    unknown <- alleles[ ! isRawAllele(alleles) ]
    if ( length(unknown) > 0 ) {
      stop("unknown raw SNP alleles (", toString(unknown),
           ") in file - '", infile, "'")
    }
    
    # If genotype qualities available, create  
    # object with quality-scaled genotypes..
    if (reading.qualities) {
      
      # Set vector of variant error probabilities.
      variant.error.probs <- 10 ^ ( -0.1 * var.qual[snp.indices] )
      
      # Set matrix of genotype error probabilities.
      genotype.error.probs <- 10 ^ ( -0.1 * var.geno$GQ[snp.indices,, drop=FALSE] )
      
      # Calculate error probability for each sample genotype
      # from converted variant/genotype quality scores.
      qual.matrix <- sapply( seq_along(variant.error.probs), function(i)
        variant.error.probs[i] + genotype.error.probs[i, ])
      
      # If qual matrix simplified, restore to matrix.
      if ( ! is.matrix(qual.matrix) ) {
        qual.matrix <- matrix(qual.matrix, ncol=length(variant.error.probs))
      }
      
      # Replace NA values with maximum error probability.
      qual.matrix[ is.na(qual.matrix) ] <- const$qual$prob$range[2]
      
      # Clamp quality scores within range of error probabilities.
      qual.matrix <- sapply(seq_along(file.snps), function(i)
        clamp(qual.matrix[, i], const$qual$prob$range))
      
      # If qual matrix simplified, restore to matrix.
      if ( ! is.matrix(qual.matrix) ) {
        qual.matrix <- matrix(qual.matrix, ncol=length(file.snps))
      }
      
      geno.data <- c(geno.matrix, qual.matrix)
      
    } else { # ..otherwise create object with only genotypes.
      
      geno.data <- geno.matrix
    }
    
    # Prep file variant data.
    data.names <- list(file.samples, file.snps, slices)
    data.shape <- lengths(data.names)
    
    # Set file variant data.
    snps[[i]] <- array(geno.data, dim=data.shape, dimnames=data.names)
  }
  
  # If multiple relevant files, combine SNP genotypes,
  # taking the union of all variants in input files..
  if ( length(infiles) > 1 ) {
    
    # Prep combined variant data.
    combined.samples <- sort( unique( unlist( lapply(snps, rownames) ) ) )
    combined.snps <- sort( unique( unlist( lapply(snps, colnames) ) ) )
    combined.names <- list(combined.samples, combined.snps, slices)
    combined.shape <- lengths(combined.names)
    
    # Init combined variant data.
    result <- array(const$missing.value, dim=combined.shape, dimnames=combined.names)
    
    # Set combined variant data from each input file.
    # NB: we previously checked for duplicate samples and coinciding
    # variants, so we can assume that there will be no conflicts.
    for ( i in seq_along(snps) ) {
      for ( slice in slices ) {
        for ( sample.id in rownames(snps[[i]]) ) {
          for ( snp.id in colnames(snps[[i]]) ) {
            result[sample.id, snp.id, slice] <- snps[[i]][sample.id, snp.id, slice]
          }
        }
      }
    }
    
  } else { # ..otherwise take single set of SNP genotypes.
    
    result <- snps[[1]]
  }
  
  # Get combined number of SNP variants.
  num.snps <- ncol(result)
  stopifnot( num.snps > 0 )
  
  # If filters specified, filter SNP variants.
  if ( require.all || require.any || require.polymorphic ) {
    
    mask <- rep(TRUE, num.snps)
    
    if (require.all) { # Remove incomplete variants.
      mask <- mask & sapply(1:num.snps, function(i)
        all(result[, i, 'geno'] != const$missing.value))
    }
    
    if (require.any) { # Remove variants without genotypes.
      mask <- mask & sapply(1:num.snps, function(i)
        any(result[, i, 'geno'] != const$missing.value))
    }
    
    if (require.polymorphic) { # Remove monomorphic variants.
      mask <- mask & sapply(1:num.snps, function(i) {
        geno.data <- result[, i, 'geno']
        g.symbols <- unique(geno.data)
        genotypes <- g.symbols[ g.symbols != const$missing.value ]
        return( length(genotypes) > 1 )
      })
    }
    
    # Apply filter mask.
    result <- result[, mask,, drop=FALSE]
  }
  
  return(result)
}

# End of vcf.R #################################################################


# Start of util.R ##############################################################

# addXInfo ---------------------------------------------------------------------
#' Add experiment info to the given matrix.
#' 
#' @param tab A \code{matrix} to which experiment info will be added.
#' @param col.index Index of the column in \code{tab} at which experiment
#' info will be inserted. Columns in the original \code{tab} matrix that
#' have an index equal to or greater than \code{col.index} will be moved
#' to the right of the inserted experiment info in the returned matrix.
#' @param xinfo A \code{matrix} of experiment info.
#' 
#' @return The input matrix with experiment info added at the specified column.
#' 
#' @keywords internal
#' @rdname addXInfo
addXInfo <- function(x, xinfo) {
  
  stopifnot( is.matrix(x) )
  stopifnot( all( dim(x) > 0 ) )
  stopifnot( typeof(x) == 'character' )
  stopifnot( colnames(x)[1] == 'File' )
  stopifnot( is.matrix(xinfo) )
  stopifnot( all( dim(xinfo) > 0 ) )
  stopifnot( typeof(x) == 'character' )
  stopif( anyDuplicated( rownames(xinfo) ) )
  stopif( any( colnames(xinfo) %in% colnames(x) ) )
  
  scanfiles <- as.character(x[, 'File'])
  
  index.list <- lapply( seq_along(scanfiles), function(i)
    which( rownames(xinfo) == scanfiles[i] ) )
  
  mxinfo <- matrix(NA_character_, nrow=nrow(x), ncol=ncol(xinfo),
                   dimnames=list(NULL, colnames(xinfo)))
  for ( i in seq_along(index.list) ) {
    if ( length(index.list[[i]]) == 1 ) {
      mxinfo[i, ] <- xinfo[ index.list[[i]], ]
    }
  }
  
  x <- cbind(x[, 1, drop=FALSE], mxinfo, x[, -1, drop=FALSE])
  
  return(x)
}

# allKwargs --------------------------------------------------------------------
#' Test if ellipsis arguments are all keyword arguments.
#' 
#' @param ... Ellipsis arguments.
#'     
#' @return \code{TRUE} if all ellipsis arguments are keyword arguments;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname allKwargs
allKwargs <- function(...) {
  args <- list(...)
  return( length(args) == 0 || hasNames(args) )
}

# allNA ------------------------------------------------------------------------
#' Test if all elements are \code{NA} values.
#' 
#' @param x Test vector.
#'     
#' @return \code{TRUE} if vector is of length zero or
#' contains only \code{NA} values; \code{FALSE} otherwise.
#' 
#' @template author-thomas-walsh
#' @template author-yue-hu
#' 
#' @keywords internal
#' @rdname allNA
allNA <- function(x) {
  stopifnot( is.vector(x) || is.factor(x) )
  return( all( is.na(x) ) )
}

# allWhite ---------------------------------------------------------------------
#' Test if vector is all whitespace.
#' 
#' @param x Character vector.
#' 
#' @return \code{TRUE} if character vector is of length zero or
#' contains only whitespace characters; \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname allWhite
allWhite <- function(x) {
  stopifnot( is.character(x) )
  return( all( grepl( '^[[:space:]]*$', x)  ) )
}

# anyKwargs --------------------------------------------------------------------
#' Test if any ellipsis arguments are keyword arguments.
#' 
#' @param ... Ellipsis arguments.
#' 
#' @return \code{TRUE} if any ellipsis arguments are keyword arguments;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname anyKwargs
anyKwargs <- function(...) {
  args <- list(...)
  return( length(args) > 0 && ! is.null( names(args) ) &&
            any( ! is.na(names(args)) & names(args) != '' ) )
}

# bstripBlankRows --------------------------------------------------------------
#' Strip blank rows from bottom of \code{data.frame}.
#' 
#' @param x A \code{data.frame} with columns of type \code{character}.
#' 
#' @return Input object in which bottommost blank rows (i.e. rows
#' containing no non-whitespace characters) have been stripped.
#' 
#' @keywords internal
#' @rdname bstripBlankRows
bstripBlankRows <- function(x) {
  stopifnot( is.data.frame(x) )
  stopifnot( all( sapply(x, class) == 'character' ) )
  while( allWhite( as.character( x[nrow(x), ]) ) ) {
    x <- x[-nrow(x),, drop=FALSE]
  }
  return(x)
}

# clamp ------------------------------------------------------------------------
#' Clamp numbers within a range.
#' 
#' @param n Numeric vector.
#' @param interval Numeric vector containing the minimum and maximum values of
#' the range, respectively. 
#'     
#' @return Input vector with values clamped within the specified range.
#' 
#' @keywords internal
#' @rdname clamp
clamp <- function(n, interval) {
  
  stopifnot( is.numeric(n) )
  stopifnot( is.numeric(interval) )
  stopifnot( length(interval) == 2 )
  stopif( anyNA(interval) )
  stopif( any( is.nan(interval) ) )
  stopifnot( diff(interval) >= 0 )
  
  n[ ! is.na(n) & n < interval[1] ] <- interval[1]
  n[ ! is.na(n) & n > interval[2] ] <- interval[2]
  
  return(n)
}

# coerceDataFrame --------------------------------------------------------------
#' Coerce \code{data.frame} columns to the specified classes.
#' 
#' @param x A \code{data.frame}.
#' @param classes Classes to which columns of \code{data.frame} are to be 
#' coerced. This should be a character vector containing class names, with 
#' the same number of elements as there are columns in the \code{data.frame}.
#'     
#' @return Coerced \code{data.frame}.
#' 
#' @keywords internal
#' @rdname coerceDataFrame
coerceDataFrame <- function(x, classes) {
  
  stopifnot( is.data.frame(x) )
  stopifnot( nrow(x) > 0 )
  stopifnot( ncol(x) > 0 )
  stopifnot( is.character(classes) )
  stopifnot( length(classes) == ncol(x) )
  
  coercions <- lapply(classes, getCoercionFromClassS3)
  
  stopif( any( is.null(coercions) ) )
  
  for ( i in getColIndices(x) ) {
    x[[i]] <- coercions[[i]](x[[i]])
  }
  
  return(x)
}

# deleteColumn -----------------------------------------------------------------
#' Delete column from object.
#' 
#' @param x A \code{data.frame} or \code{matrix}.
#' @param col.index Column index of column to remove.
#' @param col.name Column name of column to remove.
#'     
#' @return Input object with column removed.
#' 
#' @keywords internal
#' @rdname deleteColumn
deleteColumn <- function(x, col.index=NULL, col.name=NULL) {
  
  stopifnot( is.data.frame(x) || is.matrix(x) )
  
  if ( ! xor( is.null(col.index), is.null(col.name) ) ) {
    stop("deleteColumn takes either a column index or column name, but not both")
  }
  
  if ( ! is.null(col.index) ) {
    stopifnot( isSingleWholeNumber(col.index) )
    stopifnot( inRange( col.index, c( 1, ncol(x) ) ) )
  } else {
    stopifnot( isSingleString(col.name) )
    stopif( is.null( colnames(x) ) )
    stopif( anyDuplicated( colnames(x) ) )
    col.index <- which( colnames(x) == col.name )
  } 
  
  others <- otherattributes(x)
  
  # Delete column.
  x <- x[, -col.index, drop=FALSE]
  
  otherattributes(x) <- others
  
  return(x)
}

# dispatchFromClassS3 ----------------------------------------------------------
#' Dispatch method with respect to the given class vector.
#' 
#' @param generic Generic function name.
#' @param class.vector Class vector of an R object.
#' @param package Package from which available methods are taken.
#'     
#' @return Method suited to the specified class.
#' 
#' @importFrom utils lsf.str
#' @keywords internal
#' @rdname dispatchFromClassS3
dispatchFromClassS3 <- function(generic, class.vector, package) {
  
  stopifnot( isSingleString(generic) )
  stopifnot( is.character(class.vector) )
  stopifnot( length(class.vector) > 0 )
  stopifnot( isSingleString(package) )
  
  pkg <- paste0('package:', package)
  pattern <- paste0('^', generic, '[.]([^[:space:]]+)')
  x <- utils::lsf.str(pkg, pattern=pattern)
  
  m <- regexec(pattern, x)
  matches <- regmatches(x, m)
  
  functions <- sapply(matches, getElement, 1)
  suffixes <- sapply(matches, getElement, 2)
  
  default.index <- which( suffixes == 'default' )
  class.indices <- which( suffixes != 'default' )
  
  func.name <- NULL
  
  if ( length(class.indices) > 0 ) {
    
    classes <- suffixes[ class.indices ]
    
    class.match <- match(classes, class.vector)
    
    if ( length(class.match) > 0 && ! allNA(class.match) ) {
      
      first.match <- min(class.match, na.rm = TRUE)
      cls <- classes[ ! is.na(class.match) & class.match == first.match ]
      func.name <- functions[ suffixes == cls ]
    }
  }
  
  if ( is.null(func.name) ) {
    
    if ( length(default.index) > 0 ) {
      
      func.name <- functions[default.index]
      
    } else {
      
      stop("no matching class, no default")
    }
  }
  
  return( get(func.name) )
}

# ellipt -----------------------------------------------------------------------
#' Ellipt string to the given length.
#' 
#' @param s String to ellipt.
#' @param n Maximum number of characters in the output string. If the string is
#' ellipted, this includes the number of characters in the ellipses themselves.
#' @param left Indicates if the string should be
#' ellipted at the start ('left') of the string.
#' @param right Indicates if the string should be
#' ellipted at the end ('right') of the string.
#' 
#' @return Input string ellipted to the specified length.
#' 
#' @keywords internal
#' @rdname ellipt
ellipt <- function(s, n, left=FALSE, right=FALSE) {
  
  stopifnot( isSingleString(s) )
  stopifnot( isSinglePositiveNumber(n) )
  stopifnot( isBOOL(left) )
  stopifnot( isBOOL(right) )
  
  ellipsis = '...'
  
  if ( left && right ) {
    ellipses.nchar <- 2 * nchar(ellipsis)
  } else {
    ellipses.nchar <- nchar(ellipsis)
  }
  
  # Get length of input string that will
  # remain after it has been ellipted.
  ellipted.nchar <- n - ellipses.nchar
  
  stopifnot( ellipted.nchar > 0 )
  
  # Ellipt input string if longer than specified number of characters.
  if ( nchar(s) > n ) {
    
    if ( left && right ) {
      
      m <- nchar(s) %/% 2
      rm <- nchar(s) %% 2
      h <- ellipted.nchar %/% 2
      rh <- ellipted.nchar %% 2
      
      i <- (m + rm) - (h + rh) + 1
      j <- i + ellipted.nchar - 1
      
      s <- paste0(ellipsis, substr(s, i, j), ellipsis)
      
    } else if (left) {
      
      j <- nchar(s)
      i <- j - ellipted.nchar + 1
      
      s <- paste0(ellipsis, substr(s, i, j))
      
    } else if (right) {
      
      i <- 1
      j <- ellipted.nchar
      
      s <- paste0(substr(s, i, j), ellipsis)
      
    } else {
      
      h <- ellipted.nchar %/% 2
      rh <- ellipted.nchar %% 2
      
      i1 <- 1
      j1 <- h + rh
      j2 <- nchar(s)
      i2 <- j2 - h + 1
      
      s <- paste0(substr(s, i1, j1), ellipsis, substr(s, i2, j2))
    }
  }
  
  return(s)
}

# emptyArgs --------------------------------------------------------------------
#' Test if ellipsis arguments are empty.
#' 
#' @param ... Ellipsis arguments.
#' 
#' @return \code{TRUE} if there are no ellipsis arguments;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname emptyArgs
emptyArgs <- function(...) {
  return( length( list(...) ) == 0 )
}

# getCoercionFromClassS3 -------------------------------------------------------
#' Get coercion function for the given class.
#' 
#' @param class.vector Class vector of an R object.
#'     
#' @return Coercion function.
#' 
#' @importFrom utils as.roman
#' @keywords internal
#' @rdname getCoercionFromClassS3
getCoercionFromClassS3 <- function(class.vector) {
  
  stopifnot( is.character(class.vector) )
  stopifnot( length(class.vector) > 0 )
  
  known.coercions <- list(
    character=as.character,
    double=as.double,
    factor=as.factor,
    integer=as.integer,
    logical=as.logical,
    numeric=as.numeric,
    raw=as.raw,
    roman=utils::as.roman
  )
  
  coercion <- NULL
  
  for ( class.name in class.vector ) {
    
    if ( class.name %in% names(known.coercions) ) {
      coercion <- known.coercions[[class.name]]
      break
    }
  }
  
  return(coercion)
}

# getColIndices ----------------------------------------------------------------
#' Get column indices of object.
#' 
#' @param x An object with columns.
#' @param requested A character vector of column names, a logical vector of
#' length equal to the number of object columns, or a numeric vector of column
#' indices of the input object. If this is specified, only the requested column
#' indices are returned; otherwise, this function returns all column indices.
#' @param strict Option indicating that \code{requested}, if specified, must
#' be in the same order as the corresponding columns of the input object.
#' 
#' @return Integer vector of column indices.
#' 
#' @keywords internal
#' @rdname getColIndices
getColIndices <- function(x, requested=NULL, strict=FALSE) {
  
  num.cols <- ncol(x)
  
  if ( ! isSingleNonNegativeWholeNumber(num.cols) ) {
    stop("cannot get column indices - object does not have columns")
  }
  
  available <- if ( num.cols > 0 ) { 1:num.cols } else { integer() }
  names(available) <- colnames(x)
  
  resolved <- getIndices(available, requested=requested, strict=strict)
  indices <- unname(available[resolved])
  
  return(indices)
}

# getColMask -------------------------------------------------------------------
#' Get logical mask of object columns.
#' 
#' Get logical mask of the columns of an object, as constrained by the
#' \code{requested} parameter. If using this function without \code{requested}
#' constraints, consider using the command \code{rep(TRUE, ncol(x))} instead.
#' 
#' @param x An object with columns.
#' @param requested A character vector of column names, a logical vector of
#' length equal to the number of object columns, or a numeric vector of column
#' indices of the input object. If this parameter is not specified, the returned
#' logical mask is \code{TRUE} for all columns.
#' 
#' @return Logical mask of object columns.
#' 
#' @keywords internal
#' @rdname getColMask
getColMask <- function(x, requested=NULL) {
  available <- getColIndices(x)
  requested <- getColIndices(x, requested=requested)
  return( available %in% requested )
}

# getIdColIndex ----------------------------------------------------------------
#' Get sample ID column index.
#' 
#' Get the index of the sample ID column in the given object. For example,
#' in a \code{cross} object, this is the column of the phenotype
#' \code{data.frame} with the heading \code{'ID'} (case-insensitive). It can
#' also be used to get the ID column index directly from a \code{data.frame},
#' in which case the ID column is that which has heading \code{'ID'}
#' (also case-insensitive).
#' 
#' @param x An object that may contain a sample ID column.
#' 
#' @return Sample ID column index.
#' 
#' @export
#' @family cross object functions
#' @rdname getIdColIndex
getIdColIndex <- function(x) {
  UseMethod('getIdColIndex', x)
}

# getIdColIndex.cross ----------------------------------------------------------
#' @export
#' @rdname getIdColIndex
getIdColIndex.cross <- function(x) {
  
  id.col <- which( tolower( colnames(x$pheno) ) == 'id' )
  
  if ( length(id.col) == 0 ) {
    id.col <- NULL
  } else if ( length(id.col) > 1 ) {
    stop("multiple ID columns found")
  }
  
  return(id.col)
}

# getIdColIndex.data.frame -----------------------------------------------------
#' @export
#' @method getIdColIndex data.frame
#' @rdname getIdColIndex
getIdColIndex.data.frame <- function(x) {
  
  id.col <- which( tolower( colnames(x) ) == 'id' )
  
  if ( length(id.col) == 0 ) {
    id.col <- NULL
  } else if ( length(id.col) > 1 ) {
    stop("multiple ID columns found")
  }
  
  return(id.col)
}

# getIndices -------------------------------------------------------------------
#' Get indices of object elements.
#' 
#' Get the indices of an object, as constrained by the \code{requested}
#' parameter. If using this function without \code{requested} constraints,
#' consider using the faster primitive R function \code{seq_along} instead.
#' 
#' @param x An object with elements that are accessible by an index.
#' @param requested A character vector of names, a logical vector of the same
#' length as the input object, or a numeric vector containing indices of the
#' input object. If this parameter is not specified, all indices are returned.
#' @param strict Option indicating that \code{requested}, if specified, must
#' be in the same order as the corresponding elements of the input object.
#' 
#' @return Integer vector of object indices.
#' 
#' @keywords internal
#' @rdname getIndices
getIndices <- function(x, requested=NULL, strict=FALSE) {
  UseMethod('getIndices', x)
}

# getIndices.default -----------------------------------------------------------
#' @export
#' @rdname getIndices
getIndices.default <- function(x, requested=NULL, strict=FALSE) {
  
  stopifnot( isBOOL(strict) )
  
  object.length <- length(x)
  
  if ( ! isSingleNonNegativeWholeNumber(object.length) ) {
    stop("cannot get object indices - object does not have length")
  }
  
  indices <- seq_along(x)
  
  if ( ! is.null(requested) ) {
    
    if ( is.numeric(requested) ) {
      
      nonintegers <- requested[ ! isWholeNumber(requested) ]
      if ( length(nonintegers) > 0 ) {
        stop("requested indices are not integers - '", toString(nonintegers), "'")
      }
      
      exrange <- requested[ ! requested %in% indices ]
      if ( length(exrange) > 0 ) {
        stop("requested indices out of range - '", toString(exrange), "'")
      }
      
      indices <- indices[requested]
      
    } else if ( is.logical(requested) ) {
      
      if ( length(requested) != object.length ) {
        stop("cannot resolve indices by logical vector - length mismatch")
      }
      
      indices <- unname( which(requested) )
      
    } else if ( is.character(requested) ) {
      
      object.names <- names(x)
      
      if ( anyNA(requested) ) {
        stop("cannot resolve indices by name - requested names are incomplete")
      }
      
      if ( is.null(object.names) ) {
        stop("cannot resolve indices by name - no object names found")
      }
      
      if ( anyDuplicated(object.names) ) {
        stop("cannot resolve indices by name - duplicate object names found")
      }
      
      unfound <- requested[ ! requested %in% object.names ]
      if ( length(unfound) > 0 ) {
        stop("requested names not found - '", toString(unfound), "'")
      }
      
      indices <- match(requested, object.names)
      
    } else {
      
      stop("requested indices must be specified by index, logical mask, or name")
    }
    
    if (strict) { # NB: also ensures no duplicates
      if ( is.unsorted(indices, strictly=TRUE) ) {
        stop("requested indices not specified in strictly increasing order")
      }
    }
  }
  
  return(indices)
}

# getIndices.qtl ---------------------------------------------------------------
#' @export
#' @rdname getIndices
getIndices.qtl <- function(x, requested=NULL, strict=FALSE) {
  
  available <- seq(x$n.qtl)
  names(available) <- x$name
  
  resolved <- getIndices(available, requested=requested, strict=strict)
  indices <- unname(available[resolved])
  
  return(indices)
}

# getLodColIndex ---------------------------------------------------------------
#' Get LOD column index.
#' 
#' @param x A \code{scanone} or equivalent object,
#' or any object with LOD-column-associated elements.
#' @template param-lodcolumn
#' 
#' @return LOD column index.
#' 
#' @keywords internal
#' @rdname getLodColIndex
getLodColIndex <- function(x, lodcolumn=NULL) {
  
  lodcol.index <- getLodColIndices(x, lodcolumns=lodcolumn)
  
  if ( length(lodcol.index) > 1 ) {
    stop("object has multiple LOD columns - please choose one")
  } else if ( length(lodcol.index) == 0 ) {
    stop("no LOD column found")
  }
  
  return(lodcol.index)
}

# getLodColIndices -------------------------------------------------------------
#' Get LOD column indices.
#' 
#' @param x A \code{scanone} or equivalent object,
#' or any object with LOD-column-associated elements.
#' @template param-lodcolumns
#' @param strict Option indicating that the LOD columns, if specified, must
#' be in the same order as the corresponding elements of the input object.
#' 
#' @return Vector of LOD column indices.
#' 
#' @keywords internal
#' @rdname getLodColIndices
getLodColIndices <- function(x, lodcolumns=NULL, strict=FALSE) {
  UseMethod('getLodColIndices', x)
}

# getLodColIndices.mapframe ----------------------------------------------------
#' @rdname getLodColIndices
getLodColIndices.mapframe <- function(x, lodcolumns=NULL, strict=FALSE) {
  return( getDatColIndices(x, datcolumns=lodcolumns, strict=strict) )
}

# getLodColIndices.scanone -----------------------------------------------------
#' @rdname getLodColIndices
getLodColIndices.scanone <- function(x, lodcolumns=NULL, strict=FALSE) {
  return( getDatColIndices(x, datcolumns=lodcolumns, strict=strict) )
}

# getLodColIndices.scanonebins -------------------------------------------------
#' @rdname getLodColIndices
getLodColIndices.scanonebins <- function(x, lodcolumns=NULL, strict=FALSE) {
  available <- seq_len( dim(x)[3] )
  names(available) <- dimnames(x)[[3]]
  resolved <- getIndices(available, requested=lodcolumns, strict=strict)
  indices <- unname(available[resolved])
  return(indices)
}

# getLodColIndices.scanoneperm -------------------------------------------------
#' @rdname getLodColIndices
getLodColIndices.scanoneperm <- function(x, lodcolumns=NULL, strict=FALSE) {
  available <- getColIndices(x)
  names(available) <- colnames(x)
  resolved <- getIndices(available, requested=lodcolumns, strict=strict)
  indices <- unname(available[resolved])
  return(indices)
}

# getLodColIndices.scantwo -----------------------------------------------------
#' @rdname getLodColIndices
getLodColIndices.scantwo <- function(x, lodcolumns=NULL, strict=FALSE) {
  available <- if ( is.matrix(x$lod) ) { 1L } else { seq( dim(x$lod)[[3]] ) }
  names(available) <- attr(x, 'phenotypes')
  resolved <- getIndices(available, requested=lodcolumns, strict=strict)
  indices <- unname(available[resolved])
  return(indices)
}

# getLodColIndices.scantwoperm -------------------------------------------------
#' @rdname getLodColIndices
getLodColIndices.scantwoperm <- function(x, lodcolumns=NULL, strict=FALSE) {
  phenames.vectors <- lapply(unname(x), colnames)
  stopifnot( length( unique(phenames.vectors) ) == 1 )
  pheno.names <- phenames.vectors[[1]]
  available <- structure(seq_along(pheno.names), names=pheno.names)
  resolved <- getIndices(available, requested=lodcolumns, strict=strict)
  indices <- unname(available[resolved])
  return(indices)
}

# getLodColIndices.summary.scanonebins -----------------------------------------
#' @export
#' @method getLodColIndices summary.scanonebins
#' @rdname getLodColIndices
getLodColIndices.summary.scanonebins <- function(x, lodcolumns=NULL,
                                                 strict=FALSE) {
  available <- getColIndices(x)
  names(available) <- colnames(x)
  resolved <- getIndices(available, requested=lodcolumns, strict=strict)
  indices <- unname(available[resolved])
  return(indices)
}

# getLodColIndices.summary.scanoneperm -----------------------------------------
#' @export
#' @method getLodColIndices summary.scanoneperm
#' @rdname getLodColIndices
getLodColIndices.summary.scanoneperm <- function(x, lodcolumns=NULL,
                                                 strict=FALSE) {
  available <- getColIndices(x)
  names(available) <- colnames(x)
  resolved <- getIndices(available, requested=lodcolumns, strict=strict)
  indices <- unname(available[resolved])
  return(indices)
}

# getLodColIndices.summary.scantwoperm -----------------------------------------
#' @export
#' @method getLodColIndices summary.scantwoperm
#' @rdname getLodColIndices
getLodColIndices.summary.scantwoperm <- function(x, lodcolumns=NULL,
                                                 strict=FALSE) {
  
  # Get available LOD column indices.
  lodcolumn.counts <- lapply(unname(x), ncol)
  stopifnot( length( unique(lodcolumn.counts) ) == 1 )
  available <- seq(lodcolumn.counts[[1]])
  
  # Set LOD column names from phenotypes, if available.
  names(available) <-  getPhenotypes.summary.scantwoperm(x) # TODO: implement generic.
  
  # Resolve LOD column indices.
  resolved <- getIndices(available, requested=lodcolumns, strict=strict)
  indices <- unname(available[resolved])
  
  return(indices)
}

# getMask ----------------------------------------------------------------------
#' Get logical mask of object elements.
#' 
#' Get logical mask of the elements of an object, as constrained by the
#' \code{requested} parameter. If using this function without \code{requested}
#' constraints, consider using the command \code{rep(TRUE, length(x))} instead.
#' 
#' @param x An object that is subsettable by a logical mask.
#' @param requested A character vector of names, a logical vector of the same
#' length as the input object, or a numeric vector containing indices of the
#' input object. If this parameter is not specified, the returned logical mask
#' is \code{TRUE} for all elements.
#' 
#' @return Logical mask of object elements.
#' 
#' @keywords internal
#' @rdname getMask
getMask <- function(x, requested=NULL) {
  available <- getIndices(x)
  requested <- getIndices(x, requested=requested)
  return( available %in% requested )
}

# getMissingValueFromClassS3 ---------------------------------------------------
#' Get missing value for the given class.
#' 
#' @param class.vector Class vector of an R object.
#'     
#' @return Missing value.
#' 
#' @keywords internal
#' @rdname getMissingValueFromClassS3
getMissingValueFromClassS3 <- function(class.vector) {
  
  stopifnot( is.character(class.vector) )
  stopifnot( length(class.vector) > 0 )
  
  known.missing <- list(
    character = NA_character_,
    double    = NA_real_,
    factor    = NA_character_,
    integer   = NA_integer_,
    logical   = NA,
    numeric   = NA_real_
  )
  
  missing.value <- NULL
  
  for ( class.name in class.vector ) {
    
    if ( class.name %in% names(known.missing) ) {
      missing.value <- known.missing[[class.name]]
      break
    }
  }
  
  return(missing.value)
}

# getPhenoColIndices -----------------------------------------------------------
#' Get phenotype column indices.
#' 
#' Get the indices of the phenotype columns in the given object. In a
#' \code{data.frame}, phenotype columns are those columns that do not have
#' reserved phenotype headings (e.g. \code{'ID'}). In a \code{cross} object,
#' these are the corresponding columns of the phenotype \code{data.frame}.
#' 
#' @param x A \code{data.frame} containing phenotype data, or an \pkg{R/qtl}
#' \code{cross} object containing such a \code{data.frame}.
#' @param pheno.col Vector indicating the phenotype indices to return. This can
#' be an integer vector of phenotype column indices with respect to the columns
#' of the \code{data.frame} (which may or may not be within a \code{cross}
#' object), or a character vector that contains phenotype IDs or their
#' syntactically valid names. If none are specified, all phenotype column
#' indices are returned.
#' 
#' @return Phenotype column indices.
#' 
#' @export
#' @family cross object functions
#' @rdname getPhenoColIndices
getPhenoColIndices <- function(x, pheno.col=NULL) {
  UseMethod('getPhenoColIndices', x)
}

# getPhenoColIndices.cross -----------------------------------------------------
#' @export
#' @rdname getPhenoColIndices
getPhenoColIndices.cross <- function(x, pheno.col=NULL) {
  return( getPhenoColIndices(x$pheno, pheno.col=pheno.col) )
}

# getPhenoColIndices.data.frame ------------------------------------------------
#' @export
#' @method getPhenoColIndices data.frame
#' @rdname getPhenoColIndices
getPhenoColIndices.data.frame <- function(x, pheno.col=NULL) {
  
  stopif( anyNA( colnames(x) ) )
  
  pheno.names <- colnames(x)
  
  reserved.indices <- which( tolower(pheno.names) %in%
                               const$reserved.phenotypes )
  
  if ( ! is.null(pheno.col) ) {
    
    stopifnot( length(pheno.col) > 0 )
    stopif( anyNA(pheno.col) )
    
    if ( is.character(pheno.col) ) {
      
      indices <- vector('integer', length(pheno.col))
      
      for ( i in seq_along(pheno.col) ) {
        
        p <- pheno.col[i]
        
        pheno.index <- which( pheno.names == make.names(p) )
        
        if ( length(pheno.index) == 0 ) {
          stop("index not found for pheno.col - '", p, "'")
        }
        
        if ( length(pheno.index) > 1 ) {
          stop("multiple indices found for pheno.col - '", p, "'")
        }
        
        if ( pheno.index %in% reserved.indices ) {
          stop("pheno.col heading matches reserved phenotype - '", p, "'")
        }
        
        indices[i] <- pheno.index
      }
      
    } else if ( is.numeric(pheno.col) && all( isWholeNumber(pheno.col) ) ) {
      
      exrange <- pheno.col[ pheno.col < 1 | pheno.col > ncol(x) ]
      if ( length(exrange) > 0 ) {
        stop("pheno.col indices out of range - '", toString(exrange), "'")
      }
      
      reserved <- pheno.col[ pheno.col %in% reserved.indices ]
      if ( length(reserved) > 0 ) {
        stop("pheno.col indices point to reserved phenotype headings - '",
             toString(reserved), "'")
      }
      
      indices <- pheno.col
      
    } else {
      
      stop("pheno.col must contain phenotype indices or phenotype names")
    }
    
  } else {
    
    indices <- getColIndices(x)
    indices <- indices[ ! indices %in% reserved.indices ]
  }
  
  return(indices)
}

# getPhenotypes.summary.scantwoperm --------------------------------------------
#' Get phenotype names.
#' 
#' @param x Object that may contain phenotype names.
#' 
#' @return Character vector of phenotype names. Returns
#' \code{NULL} if the object does not contain phenotype names.
#' 
#' @export
#' @keywords internal
#' @rdname getPhenotypes.summary.scantwoperm
getPhenotypes.summary.scantwoperm <- function(x) {
  
  stopifnot( 'summary.scantwoperm' %in% class(x) ) # TODO: implement generic.
  
  # Assume consistent phenotype names.
  inconsistent <- FALSE
  
  # Get column names of each scantwo threshold matrix.
  phenames.vectors <- lapply(unname(x), colnames)
  
  # Get names of thresholds for which threshold matrix has no column names.
  unnamed.columns <- names(x)[ sapply(phenames.vectors, is.null) ]
  
  # If all threshold matrices have column names..
  if ( length(unnamed.columns) == 0 ) {
    
    # ..and if those column names are consistent, take as phenotype names..
    if ( length( unique(phenames.vectors) ) == 1 ) {
      phenotypes <- phenames.vectors[[1]]
    } else {
      inconsistent <- TRUE
    }
    
    # ..otherwise check that phenotype names are at least
    # consistently absent from all threshold matrices.
  } else if ( length(unnamed.columns) == length(x) ) {
    phenotypes <- NULL
  } else {
    inconsistent <- TRUE
  }
  
  if (inconsistent) {
    stop("inconsistent phenotype names in 'summary.scantwoperm' object")
  }
  
  return(phenotypes)
}

# getRowIndices ----------------------------------------------------------------
#' Get row indices of object.
#' 
#' @param x An object with rows.
#' @param requested A character vector of row names, a logical vector of length
#' equal to the number of object rows, or a numeric vector of row indices of the
#' input object. If this is specified, only the requested row indices are
#' returned; otherwise, this function returns all row indices.
#' @param strict Option indicating that \code{requested}, if specified, must
#' be in the same order as the corresponding rows of the input object.
#' 
#' @return Integer vector of row indices.
#' 
#' @keywords internal
#' @rdname getRowIndices
getRowIndices <- function(x, requested=NULL, strict=FALSE) {
  
  num.rows <- nrow(x)
  
  if ( ! isSingleNonNegativeWholeNumber(num.rows) ) {
    stop("cannot get row indices - object does not have rows")
  }
  
  available <- if ( num.rows > 0 ) { 1:num.rows } else { integer() }
  names(available) <- rownames(x)
  
  resolved <- getIndices(available, requested=requested, strict=strict)
  indices <- unname(available[resolved])
  
  return(indices)
}

# getRowMask -------------------------------------------------------------------
#' Get logical mask of object rows.
#' 
#' Get logical mask of the rows of an object, as constrained by the
#' \code{requested} parameter. If using this function without \code{requested}
#' constraints, consider using the command \code{rep(TRUE, nrow(x))} instead.
#' 
#' @param x An object with rows.
#' @param requested A character vector of row names, a logical vector of length
#' equal to the number of object rows, or a numeric vector of row indices of the
#' input object. If this parameter is not specified, the returned logical mask
#' is \code{TRUE} for all rows.
#' 
#' @return Logical mask of object rows.
#' 
#' @keywords internal
#' @rdname getRowMask
getRowMask <- function(x, requested=NULL) {
  available <- getRowIndices(x)
  requested <- getRowIndices(x, requested=requested)
  return( available %in% requested )
}

# getRunIndexList --------------------------------------------------------------
#' Get index list of successive runs in a vector.
#' 
#' @param x A vector.
#'     
#' @return List of integer vectors, each containing indices for a run of
#' repeated values in the input vector. Each list element takes its name
#' from the corresponding repeated value. Returns an empty list if the
#' input vector is of length zero.
#'    
#' @keywords internal
#' @rdname getRunIndexList
getRunIndexList <- function(x, na.rm=FALSE) {
  
  stopifnot( is.vector(x) )
  stopifnot( isBOOL(na.rm) )
  
  if ( length(x) > 0 ) {
    
    # Get run-length encoding of vector.
    runs <- rle(x)
    
    # Set run names from RLE values.
    run.names <- runs$values
    
    # Get number of runs in RLE.
    num.runs <- unique( lengths(runs) )
    
    # Get last index of each run.
    J <- cumsum(runs$lengths)
    
    # Get first index of each run.
    if ( num.runs > 1 ) {
      I <- c( 1, sapply(J[1:(length(J)-1)], function(j) j + 1) )
    } else {
      I <- 1
    }
    
    # Remove NA values, if specified.
    if (na.rm) {
      mask <- ! is.na(runs$values)
      run.names <- run.names[mask]
      I <- I[mask]
      J <- J[mask]
    }
    
    # Set index list from run index ranges.
    index.list <- mapply(function(i, j) i:j, I, J, SIMPLIFY=FALSE)
    
    # Set names of index list from run values.
    names(index.list) <- run.names
    
  } else {
    
    index.list <- list()
  }
  
  return(index.list)
}

# getRunIndices ----------------------------------------------------------------
#' Get indices of successive runs in a run-length encoding.
#' 
#' @param x An \code{rle} object.
#' 
#' @return Integer vector of all run indices in the given run-length encoding,
#' which can be used to index into the \code{lengths} and \code{values} of the
#' given \code{rle} object. Returns an empty integer vector if the run-length
#' encoding has zero runs.
#' 
#' @keywords internal
#' @rdname getRunIndices
getRunIndices <- function(x) {
  stopifnot( 'rle' %in% class(x) )
  num.runs <- union( length(x$lengths), length(x$values) )
  stopifnot( isSingleNonNegativeWholeNumber(num.runs) )
  return( if ( num.runs > 0 ) { 1:num.runs } else { integer() } )
}

# getXInfoFromFilenames --------------------------------------------------------
#' Get experiment info from filenames.
#' 
#' @param filenames Character vector of file names.
#' @param pattern Pattern for extracting experiment info from file names. This
#' must be a valid Perl regex with named capture groups. Neither the capture
#' groups nor the pattern itself are required to match any given filename, but
#' all capture groups must have a name, and that name cannot clash with other
#' names that might be used alongside the extracted info.
#' 
#' @return Character matrix in which each row contains the values of capture
#' groups for a given filename, and each column contains values of a given
#' capture group across filenames. Unmatched capture groups are represented by
#' \code{NA} values.
#' 
#' @keywords internal
#' @rdname getXInfoFromFilenames
getXInfoFromFilenames <- function(filenames, pattern) {
  
  stopifnot( is.character(filenames) )
  stopif( anyDuplicated(filenames) )
  stopifnot( isSingleString(pattern) )
  
  # Init scanfile info.
  xinfo <- NULL
  
  if ( length(filenames) > 0 ) {
    
    # Parse scan file names by the given pattern.
    parsed <- parseFilenames(filenames, pattern)
    
    # Check for capture-group names clashing with disallowed info tags.
    clashing <- colnames(parsed)[ colnames(parsed) %in%
                                    const$disallowed.xinfotags ]
    if ( length(clashing) > 0 ) {
      stop("info tags clash with Excel headings - '",
           toString(clashing), "'")
    }
    
    # Remove empty columns.
    parsed <- removeColsNA(parsed)
    
    if ( ncol(parsed) > 0 ) {
      xinfo <- parsed
    }
  }
  
  return(xinfo)
}

# getScanoneThresholdInfo ------------------------------------------------------
#' Get \code{scanone} threshold info from object.
#' 
#' @param x An object containing a single \code{scanone} threshold,
#' or an object associated with a single \code{scanone} threshold.
#' 
#' @return A named list containing the \code{scanone} threshold
#' info that was extracted from the input object.
#' 
#' @keywords internal
#' @rdname getScanoneThresholdInfo
getScanoneThresholdInfo <- function(x) {
  UseMethod('getScanoneThresholdInfo', x)
}

# getScanoneThresholdInfo.numeric ----------------------------------------------
#' @rdname getScanoneThresholdInfo
getScanoneThresholdInfo.numeric <- function(x) {
  stopifnot( isSingleNonNegativeNumber(x) )
  return( list( threshold=unname(x), alpha=NULL, fdr=NULL) )
}

# getScanoneThresholdInfo.qtlintervals -----------------------------------------
#' @rdname getScanoneThresholdInfo
getScanoneThresholdInfo.qtlintervals <- function(x) {
  stopifnot( 'threshold' %in% names( attributes(x) ) )
  return( list( threshold=unname( attr(x, 'threshold') ),
                alpha=unname( attr(x, 'alpha') ), fdr=unname( attr(x, 'fdr') ) ) )
}

# getScanoneThresholdInfo.summary.scanonebins ----------------------------------
#' @method getScanoneThresholdInfo summary.scanonebins
#' @rdname getScanoneThresholdInfo
getScanoneThresholdInfo.summary.scanonebins <- function(x) {
  stopifnot( ncol(x) == 1 )
  stopifnot( nrow(x) == 1 )
  threshold <- x[1, 1]
  alpha <- NULL
  fdr <- 0.01 * as.numeric( sub('%', '', rownames(x)[1]) )
  return( list( threshold=threshold, alpha=alpha, fdr=fdr) )
}

# getScanoneThresholdInfo.summary.scanoneperm ----------------------------------
#' @method getScanoneThresholdInfo summary.scanoneperm
#' @rdname getScanoneThresholdInfo
getScanoneThresholdInfo.summary.scanoneperm <- function(x) {
  stopifnot( ncol(x) == 1 )
  stopifnot( nrow(x) == 1 )
  threshold <- x[1, 1]
  alpha <- 0.01 * as.numeric( sub('%', '', rownames(x)[1]) )
  fdr <- NULL
  return( list( threshold=threshold, alpha=alpha, fdr=fdr) )
}

# getSeqinfo -------------------------------------------------------------------
#' Get \pkg{GenomeInfoDb} \code{Seqinfo} object.
#' 
#' @param genome Genome for which a \pkg{GenomeInfoDb} \code{Seqinfo} object
#' should be returned. If this is not specified, a \code{Seqinfo} object is
#' returned for the current genome.
#' 
#' @return A \code{Seqinfo} object for the given genome.
#'
#' @importFrom GenomeInfoDb Seqinfo
#' @keywords internal
#' @rdname getSeqinfo
getSeqinfo <- function(genome=NULL) {
  
  seqtab <- getSeqTable(genome)
  
  seqinfo <- GenomeInfoDb::Seqinfo(
    seqnames   = seqtab$seqids,
    seqlengths = seqtab$seqlengths,
    isCircular = seqtab$isCircular,
    genome     = seqtab$genome
  )
  
  return(seqinfo)
}

# getSeqTable ------------------------------------------------------------------
#' Get genome sequence info table.
#' 
#' @param genome Genome for which sequence info should be returned. If this
#' is not specified, sequence info is returned for the current genome.
#' 
#' @return A \code{data.frame} containing sequence info for the given genome.
#' This can be used to create, but is distinct from, a \pkg{GenomeInfoDb}
#' \code{Seqinfo} object.
#' 
#' @keywords internal
#' @rdname getSeqTable
getSeqTable <- function(genome=NULL) {
  prev.genome <- genomeOpt(genome)
  on.exit( genomeOpt(prev.genome) )
  seqtab <- const$seqtab[[ genomeOpt() ]]
  return(seqtab)
}

# getSpecialAttributeNames -----------------------------------------------------
#' Get special attributes for the given object.
#' 
#' @param x R object.
#'     
#' @return Vector of special attribute names.
#' 
#' @keywords internal
#' @rdname getSpecialAttributeNames
getSpecialAttributeNames <- function(x) {
  
  default.specials <- const$special.attributes[['default']]
  
  attrset.names <- names(const$special.attributes)
  
  class.specials <- NULL
  for ( class.name in class(x) ) {
    if ( class.name %in% attrset.names ) {
      class.specials <- const$special.attributes[[class.name]]
      break
    }
  }
  
  return( union(class.specials, default.specials) )
}

# hasNames ---------------------------------------------------------------------
#' Test if object has names.
#' 
#' @param x Test object.
#'     
#' @return \code{TRUE} if object has nonempty names with no \code{NA} values;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname hasNames
hasNames <- function(x) {
  return( ! is.null( names(x) ) && all( ! is.na(names(x)) & names(x) != '' ) )
}

# hasRownames ------------------------------------------------------------------
#' Test if object has rownames.
#' 
#' @param x A \code{data.frame} or \code{matrix}.
#'     
#' @return \code{TRUE} if input object has non-default row names;
#' \code{FALSE} otherwise.
#'    
#' @keywords internal
#' @rdname hasRownames
hasRownames <- function(x) {
  
  rowname.status <- FALSE
  
  if ( is.data.frame(x) ) {
    row.names <- attr(x, 'row.names')
    rowname.status <- ! is.null(row.names) && ! all( row.names == getRowIndices(x) )
  } else if ( ! is.null( rownames(x) ) ) {
    rowname.status <- TRUE
  } 
  
  return(rowname.status) 
}

# inferFormatFromFilename ------------------------------------------------------
#' Infer format from filename.
#' 
#' @param filename A file name.
#' 
#' @return File format inferred from extension of file name.
#' Returns \code{unknown} if format could not be inferred.
#' 
#' @keywords internal
#' @rdname inferFormatFromFilename
inferFormatFromFilename <- function(filename) {
  
  stopifnot( isSingleString(filename) )
  
  fmt <- 'unknown'
  
  ext <- tolower( tools::file_ext(filename) )
  
  indices <- which( sapply(const$ext, function(fmt.exts) ext %in% fmt.exts) )
  
  if ( length(indices) == 1 ) {
    fmt <- names(const$ext)[indices]
  } else {
    stop("ambiguous extension on file - '", filename, "'")
  }
  
  return(fmt)
}

# inferStepSize ----------------------------------------------------------------
#' Infer step size from step values.
#' 
#' @param steps Numeric vector of steps (i.e. differences between consecutive 
#' values). 
#' @param tol Tolerance for step equality.
#'     
#' @return Inferred step size. Returns \code{NULL} if step size cannot be 
#' inferred.
#' 
#' @keywords internal
#' @rdname inferStepSize
inferStepSize <- function(steps, tol=.Machine$double.eps^0.5) {
  
  stopifnot( is.numeric(steps) )
  stopifnot( all( steps >= 0 ) )
  stopifnot( isSingleNonNegativeNumber(tol) )
  
  # Get frequency table of steps.
  step.freq <- table(steps)
  
  # Get numeric step sizes.
  step.sizes <- as.numeric( names(step.freq) )
  
  # Get differences between step sizes.
  step.diff <- diff(step.sizes)
  
  # Group step-size values that are very similar.
  size.groups <- vector('list')
  i <- 1L
  while ( i <= length(step.freq) ) {
    
    j <- i
    
    while ( j < length(step.freq) && step.diff[j] < tol ) {
      j <- j + 1L
    }
    
    size.groups <- append( size.groups, list(i:j) )
    
    i <- j + 1L
  }
  
  # Merge similar step values.
  merged.freq <- integer( length=length(size.groups) )
  for ( i in seq_along(size.groups) ) {
    g <- unlist(size.groups[i])
    merged.freq[i] <- sum(step.freq[g])
    names(merged.freq)[i] <- names( sort(step.freq[g], decreasing=TRUE) )[1]
  }
  
  # Sort steps by decreasing frequency.
  sorted.freq <- sort(merged.freq, decreasing=TRUE)
  
  # Verify that most frequent step is in majority.
  if( length(sorted.freq) > 1 && sorted.freq[1] < sum(sorted.freq[2:length(sorted.freq)]) ) {
    return(NULL)
  }
  
  # Get most frequent step.
  step.size <- as.numeric( names(sorted.freq)[1] )
  
  return(step.size)
}

# inRange ----------------------------------------------------------------------
#' Test if numbers lie within a range.
#' 
#' @param n Numeric vector.
#' @param interval Numeric vector containing the minimum and maximum values of
#' the range, respectively. 
#'     
#' @return Logical vector indicating which elements of the input numeric vector
#' lie within the specified range.
#' 
#' @keywords internal
#' @rdname inRange
inRange <- function(n, interval) {
  
  stopifnot( is.numeric(n) )
  stopifnot( is.numeric(interval) )
  stopifnot( length(interval) == 2 )
  stopif( anyNA(interval) )
  stopif( any( is.nan(interval) ) )
  stopifnot( diff(interval) >= 0 )
  
  return( findInterval(n, interval, rightmost.closed=TRUE) == 1 )
}

# insertColumn -----------------------------------------------------------------
#' Insert column in an object.
#' 
#' @param x A \code{data.frame} or \code{matrix}.
#' @param col.index Column index of inserted column.
#' @param col.name Column name of inserted column. If a column name is 
#' specified, but the object has no existing column names, default names 
#' (e.g. \code{'COL05'}) will be assigned to the existing columns. If a column
#' name is not specified, but the object has existing column names, the 
#' inserted column will be assigned a default column name.
#' @param data Optional vector of data to insert in the new column. The length 
#' of this vector should be evenly divisible by the number of rows in the data 
#' frame.
#'     
#' @return Input object with column inserted.
#' 
#' @keywords internal
#' @rdname insertColumn
insertColumn <- function(x, col.index, col.name=NULL, data=NA) {
  
  stopifnot( is.data.frame(x) || is.matrix(x) )
  stopifnot( isSingleWholeNumber(col.index) )
  stopifnot( inRange(col.index, c(1, ncol(x) + 1)) )
  
  stopifnot( isWholeNumber( nrow(x) / length(data) ) )
  
  others <- otherattributes(x)
  
  # Get column names of input object.
  object.colnames <- colnames(x)
  
  # Get number of columns before and after column insertion.
  prev.ncol <- ncol(x)
  post.ncol <- prev.ncol + 1
  
  # If this is a matrix, insert new column to   
  # the right of the rightmost existing column.
  if ( is.matrix(x) ) {
    x <- cbind( x, rep( NA, nrow(x) ) )
  }
  
  # If new column index is within previously existing columns, 
  # nudge subsequent columns one column to the right.
  if ( col.index <= prev.ncol ) {
    x[, (col.index + 1):post.ncol] <- x[, col.index:prev.ncol]
  }
  
  # Set new column data.
  x[, col.index] <- data
  
  # Update column names, if appropriate.
  if ( ! is.null(col.name) && ! is.null(object.colnames) ) {
    colnames(x) <- append(object.colnames, col.name, after=(col.index - 1))
  }
  
  otherattributes(x) <- others
  
  return(x)
}

# isBOOL -----------------------------------------------------------------------
#' Test for a single logical value.
#' 
#' @param x Test object.
#'      
#' @return \code{TRUE} if object is a single logical value;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isBOOL
isBOOL <- function(x) {
  return( isTRUE(x) | isFALSE(x) )
}

# isDefaultMarkerID ------------------------------------------------------------
#' Test for default marker IDs.
#'  
#' @param ids Vector of locus IDs.
#'      
#' @return Logical vector indicating which elements are valid
#' IDs that follow the pattern of a default marker ID.
#' 
#' @template section-default-marker-ids
#' 
#' @export
#' @rdname isDefaultMarkerID
isDefaultMarkerID <- function(ids) {
  return( isValidID(ids) & grepl(const$pattern$default.marker.id, ids) )
}

# isDefaultQTLName -------------------------------------------------------------
#' Test for default QTL names.
#' 
#' @param ids Vector of item IDs.
#'      
#' @return Logical vector indicating which elements are valid
#' IDs that follow the pattern of a default QTL name.
#' 
#' @template section-qtl-names
#' 
#' @export
#' @rdname isDefaultQTLName
isDefaultQTLName <- function(ids) {
  return( isValidID(ids) & grepl(const$pattern$default.qtl.name, ids) )
}

# isEnumAllele -----------------------------------------------------------------
#' Test if symbol is a valid enumerated allele.
#' 
#' @param x Test object.
#' 
#' @return \code{TRUE} if symbol is a valid enumerated allele;
#' \code{FALSE} otherwise.
#' 
#' @template section-geno-enum
#' 
#' @keywords internal
#' @rdname isEnumAllele
isEnumAllele <- function(x) {
  return( x %in% const$enum.geno.charset )
}

# isEnumGenotype ---------------------------------------------------------------
#' Test if symbol is a valid enumerated genotype.
#' 
#' @param x Test object.
#' 
#' @return \code{TRUE} if symbol is a valid enumerated
#'  genotype; \code{FALSE} otherwise.
#' 
#' @template section-geno-enum
#' 
#' @keywords internal
#' @rdname isEnumGenotype
isEnumGenotype <- function(x) {
  return( x %in% const$enum.geno.charset )
}

# isFALSE ----------------------------------------------------------------------
#' Test for a single \code{FALSE} value.
#' 
#' @param x Test object.
#'      
#' @return \code{TRUE} if object is a single \code{FALSE} value;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isFALSE
isFALSE <- function(x) {
  return( identical(FALSE, x) )
}

# isFounderAllele --------------------------------------------------------------
#' Test if symbol is a valid founder allele.
#' 
#' @param x Test object.
#' 
#' @return \code{TRUE} if symbol is a valid founder allele;
#' \code{FALSE} otherwise.
#' 
#' @template section-geno-founder
#' 
#' @keywords internal
#' @rdname isFounderAllele
isFounderAllele <- function(x) {
  return( x %in% const$founder.allele.charset )
}

# isFounderGenotype ------------------------------------------------------------
#' Test if symbol is a valid founder genotype.
#' 
#' @param x Test object.
#' @param strict Return \code{TRUE} only for complete genotypes.
#' 
#' @return \code{TRUE} if symbol is a valid founder genotype;
#' \code{FALSE} otherwise.
#' 
#' @template section-geno-founder
#' 
#' @keywords internal
#' @rdname isFounderGenotype
isFounderGenotype <- function(x, strict=FALSE) {
  
  allele.list <- strsplit(x, '')
  
  strict.charset <- const$founder.allele.charset
  
  if (strict) {
    
    result <- unlist( lapply(allele.list, function(alleles)
      all( alleles %in% strict.charset ) ) )
    
  } else {
    
    extended.charset <- c(strict.charset, const$missing.value)
    result <- unlist( lapply(allele.list, function(alleles)
      all( alleles %in% extended.charset ) &&
        any( alleles %in% strict.charset ) ) )
  }
  
  return(result)
}

# isMarkerID -------------------------------------------------------------------
#' Test for marker IDs.
#'  
#' @param loc.ids Vector of locus IDs.
#'      
#' @return Logical vector indicating which elements are valid IDs that do
#' \emph{not} follow the pattern of a pseudomarker ID.
#' 
#' @template section-marker-ids
#' 
#' @export
#' @rdname isMarkerID
isMarkerID <- function(loc.ids) {
  return( isValidID(loc.ids) & ! grepl(const$pattern$pseudomarker.id, loc.ids) )
}

# isNegativeNumber -------------------------------------------------------------
#' Test for negative numbers.
#' 
#' @param n Test vector.
#' 
#' @return Logical vector indicating which elements are less than zero.
#' 
#' @keywords internal
#' @rdname isNegativeNumber
isNegativeNumber <- function(n) {
  return( is.numeric(n) & is.finite(n) & n < 0 )
}

# isNonNegativeNumber ----------------------------------------------------------
#' Test for non-negative numbers.
#' 
#' @param n Test vector.
#'      
#' @return Logical vector indicating which elements are greater than or equal to 
#' zero.
#' 
#' @keywords internal
#' @rdname isNonNegativeNumber
isNonNegativeNumber <- function(n) {
  return( is.numeric(n) & is.finite(n) & n >= 0 )
}

# isPositiveNumber -------------------------------------------------------------
#' Test for positive numbers.
#' 
#' @param n Test vector.
#'      
#' @return Logical vector indicating which elements are positive numbers.
#' 
#' @keywords internal
#' @rdname isPositiveNumber
isPositiveNumber <- function(n) {
  return( is.numeric(n) & is.finite(n) & n > 0 )
}

# isPositiveWholeNumber --------------------------------------------------------
#' Test for positive whole numbers.
#' 
#' @param n Test vector.
#' @param tol Numeric tolerance.
#'      
#' @return Logical vector indicating which elements are positive whole numbers.
#' 
#' @keywords internal
#' @rdname isPositiveWholeNumber
isPositiveWholeNumber <- function(n, tol=.Machine$double.eps^0.5) {
  stopifnot( tol > 0 )
  return( is.numeric(n) & is.finite(n) & n > tol & abs(n - round(n)) < abs(tol) )
}

# isProbability ----------------------------------------------------------------
#' Test for valid probabilities.
#' 
#' @param n Test vector.
#'      
#' @return Logical vector indicating which elements are valid probabilities.
#' 
#' @keywords internal
#' @rdname isProbability
isProbability <- function(n) {
  return( is.numeric(n) & is.finite(n) & n >= 0 & n <= 1 )
}

# isPseudomarkerID -------------------------------------------------------------
#' Test for pseudomarker IDs.
#'  
#' @param loc.ids Vector of locus IDs.
#'      
#' @return Logical vector indicating which elements are valid IDs that follow
#' the pattern of a pseudomarker ID.
#' 
#' @template section-pseudomarker-ids
#' 
#' @export
#' @rdname isPseudomarkerID
isPseudomarkerID <- function(loc.ids) {
  return( isValidID(loc.ids) & grepl(const$pattern$pseudomarker.id, loc.ids) )
}

# isRawAllele ------------------------------------------------------------------
#' Test if symbol is a valid raw SNP allele.
#' 
#' @param x Test object.
#' 
#' @return \code{TRUE} if symbol is a valid raw SNP allele;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isRawAllele
isRawAllele <- function(x) {
  return( x %in% const$raw.allele.charset )
}

# isRawGenotype ----------------------------------------------------------------
#' Test if symbol is a valid raw SNP genotype.
#' 
#' @param x Test object.
#' @param strict Return \code{TRUE} only for complete genotypes.
#' 
#' @return \code{TRUE} if symbol is a valid raw SNP genotype;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isRawGenotype
isRawGenotype <- function(x, strict=FALSE) {
  
  allele.list <- strsplit(x, '')
  
  strict.charset <- const$raw.allele.charset
  
  if (strict) {
    
    result <- unlist( lapply(allele.list, function(alleles)
      all( alleles %in% strict.charset ) ) )
    
  } else {
    
    extended.charset <- c(strict.charset, const$missing.value)
    result <- unlist( lapply(allele.list, function(alleles)
      all( alleles %in% extended.charset ) &&
        any( alleles %in% strict.charset ) ) )
  }
  
  return(result)
}

# isSingleChar -----------------------------------------------------------------
#' Test for a single character.
#' 
#' @param x Test object.
#'      
#' @return \code{TRUE} if the object is a single character;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isSingleChar
isSingleChar <- function(x) {
  return( length(x) == 1 && is.character(x) && ! is.na(x) && nchar(x) == 1 )
}

# isSingleFiniteNumber ---------------------------------------------------------
#' Test for a single finite number.
#' 
#' @param n Test object.
#' 
#' @return \code{TRUE} if the object is a single finite number;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isSingleFiniteNumber
isSingleFiniteNumber <- function(n) {
  return( length(n) == 1 && is.numeric(n) && is.finite(n) )
}

# isSingleNonNegativeNumber ----------------------------------------------------
#' Test for a single non-negative number.
#' 
#' @param n Test object.
#'      
#' @return \code{TRUE} if the object is a single non-negative number;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isSingleNonNegativeNumber
isSingleNonNegativeNumber <- function(n) {
  return( length(n) == 1 && is.numeric(n) && is.finite(n) && n >= 0 )
}

# isSingleNonNegativeWholeNumber -----------------------------------------------
#' Test for a single non-negative whole number.
#' 
#' @param n Test object.
#' @param tol Numeric tolerance.
#'      
#' @return \code{TRUE} if the object is a single non-negative whole number;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isSingleNonNegativeWholeNumber
isSingleNonNegativeWholeNumber <- function(n, tol=.Machine$double.eps^0.5) {
  return( length(n) == 1 && is.numeric(n) && is.finite(n) && 
            n >= 0 && abs(n - round(n)) < abs(tol) )
}

# isSinglePositiveNumber -------------------------------------------------------
#' Test for a single positive number.
#' 
#' @param n Test object.
#'      
#' @return \code{TRUE} if the object is a single positive number;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isSinglePositiveNumber
isSinglePositiveNumber <- function(n) {
  return( length(n) == 1 && is.numeric(n) && is.finite(n) && n > 0 )
}

# isSinglePositiveWholeNumber --------------------------------------------------
#' Test for a single positive whole number.
#' 
#' @param n Test object.
#' @param tol Numeric tolerance.
#'      
#' @return \code{TRUE} if the object is a single positive whole number;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isSinglePositiveWholeNumber
isSinglePositiveWholeNumber <- function(n, tol=.Machine$double.eps^0.5) {
  return( length(n) == 1 && is.numeric(n) && is.finite(n) && 
            n > tol && abs(n - round(n)) < abs(tol) )
}

# isSingleProbability ----------------------------------------------------------
#' Test for a single valid probability.
#' 
#' @param n Test object.
#'      
#' @return \code{TRUE} if the object is a single valid probability;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isSingleProbability
isSingleProbability <- function(n) {
  return( length(n) == 1 && is.numeric(n) && 
            is.finite(n) && n >= 0 && n <= 1 )
}

# isSingleString ---------------------------------------------------------------
#' Test for a single character string.
#' 
#' @param x Test object.
#'      
#' @return \code{TRUE} if the object is a single string;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isSingleString
isSingleString <- function(x) {
  return( length(x) == 1 && is.character(x) && ! is.na(x) )
}

# isSingleWholeNumber ----------------------------------------------------------
#' Test for a single whole number.
#' 
#' @param n Test object.
#' @param tol Numeric tolerance.
#'      
#' @return \code{TRUE} if the object is a single whole number;
#' \code{FALSE} otherwise.
#' 
#' @keywords internal
#' @rdname isSingleWholeNumber
isSingleWholeNumber <- function(n, tol=.Machine$double.eps^0.5) {
  return( length(n) == 1 && is.numeric(n) && is.finite(n) && 
            abs(n - round(n)) < abs(tol) )
}

# isValidAllele ----------------------------------------------------------------
#' Test if symbol is a valid allele.
#' 
#' @param x Test object.
#' 
#' @return \code{TRUE} if symbol is a valid enumerated or founder allele;
#' \code{FALSE} otherwise.
#' 
#' @template section-geno-enum
#' @template section-geno-founder
#' 
#' @keywords internal
#' @rdname isValidAllele
isValidAllele <- function(x) {
  return( isEnumAllele(x) | isFounderAllele(x) )
}

# isValidGenotype --------------------------------------------------------------
#' Test if symbol is a valid genotype.
#' 
#' @param x Test object.
#' @param strict Return \code{TRUE} only for complete genotypes.
#' 
#' @return \code{TRUE} if symbol is a valid enumerated or founder genotype;
#' \code{FALSE} otherwise.
#' 
#' @template section-geno-enum
#' @template section-geno-founder
#' 
#' @keywords internal
#' @rdname isValidGenotype
isValidGenotype <- function(x, strict=FALSE) {
  return( isEnumGenotype(x) | isFounderGenotype(x, strict=strict) )
}

# isValidID --------------------------------------------------------------------
#' Test identifier validity.
#'  
#' @param x Test vector.
#'      
#' @return Logical vector indicating which elements are valid identifiers.
#' 
#' @template section-item-ids
#' 
#' @export
#' @rdname isValidID
isValidID <- function(x) {
  return( is.character(x) & nzchar(x) & grepl(const$pattern$item.id, x) )
}

# isValidName ------------------------------------------------------------------
#' Test for syntactically valid names.
#'   
#' @param x Test vector.
#'      
#' @return Logical vector indicating which elements are syntactically valid R
#' names.
#' 
#' @template section-syntactic-names
#' 
#' @export
#' @rdname isValidName
isValidName <- function(x) {
  dotvar.pattern <- '^(?:[.]{3})|(?:[.]{2}[[:digit:]]+)$'
  return( is.character(x) & x == make.names(x) & ! grepl(dotvar.pattern, x) )
}

# isWholeNumber ----------------------------------------------------------------
#' Test for whole numbers.
#' 
#' @param n Test vector.
#' @param tol Numeric tolerance.
#'      
#' @return Logical vector indicating which elements are whole numbers.
#' 
#' @keywords internal
#' @rdname isWholeNumber
isWholeNumber <- function(n, tol=.Machine$double.eps^0.5) { 
  return( is.numeric(n) & is.finite(n) & abs(n - round(n)) < abs(tol) )
}

# loadChrTable -----------------------------------------------------------------
#' Load chromosome info.
#' 
#' This function loads standard chromosome info.
#' 
#' @return A \code{data.frame} object containing chromosome info.
#' 
#' @importFrom utils read.csv
#' @keywords internal
#' @rdname loadChrTable
loadChrTable <- function() {
  
  filepath <- requestPkgDataPath('extdata', 'genomes', 'chrtab.csv')
  
  column.classes <- c(seqids='character', seqnames='character', aliases='character', 
                      isCircular='logical', genome='character')
  
  chrtab <- utils::read.csv(filepath, quote='', stringsAsFactors=FALSE,
                            strip.white=TRUE, na.strings='', colClasses=column.classes)
  
  return(chrtab)
}

# loadMapping ------------------------------------------------------------------
#' Load a mapping from a line or file.
#' 
#' Given a specified line of text or filename, the input text is loaded
#' as a YAML mapping, and returned as a simple \code{mapping} object.
#' 
#' @param line Character string containing the line to be loaded.
#' @param file Character string containing the name of the file to be loaded.
#' 
#' @return A \code{mapping} of unique string keys to individual values.
#' 
#' @importFrom yaml yaml.load
#' @importFrom yaml yaml.load_file
#' @keywords internal
#' @rdname loadMapping
loadMapping <- function(line=NULL, file=NULL) {
  
  if ( ! xor( is.null(line), is.null(file) ) ) {
    stop("mapping must be loaded from either a line or a file, but not both")
  }
  
  if ( ! is.null(line) ) {
    
    stopifnot( isSingleString(line) )
    
    # Ensure line is enclosed by curly braces.
    first.char <- substr(line, 1, 1)
    last.char <- substr(line, nchar(line), nchar(line))
    if ( first.char != '{' || first.char != '}' ) {
      line <- paste0('{', line, '}', collapse='')
    }
    
    tryCatch({ # Load YAML object from line.
      x <- yaml::yaml.load(line)
    }, error=function(e) {
      stop("failed to load YAML from line - '", toString(line), "'")
    })
    
  } else if ( ! is.null(file) ) {
    
    stopifnot( isSingleString(file) )
    
    tryCatch({ # Load YAML object from file.
      x <- yaml::yaml.load_file(file)
    }, error=function(e) {
      stop("failed to load YAML from file - '", toString(file), "'")
    })
  }
  
  tryCatch({ # Convert YAML object to mapping.
    x <- mapping(x)
  }, error=function(e) {
    stop("failed to load mapping")
  })
  
  return(x)
}

# loadSeqTables ----------------------------------------------------------------
#' Load genome sequence info tables.
#' 
#' This function loads the sequence info tables of genomes for which package
#' data is available.
#' 
#' @return A \code{list} of \code{data.frame} objects, each element named for a
#' given genome and containing sequence info for that genome. Each such
#' \code{data.frame} object can be used to create, but is distinct from,
#' a \pkg{GenomeInfoDb} \code{Seqinfo} object.
#' 
#' @importFrom utils read.csv
#' @keywords internal
#' @rdname loadSeqTables
loadSeqTables <- function() {
  
  column.classes <- c(seqids='character', seqnames='character', 
                      seqlengths='integer', maplengths='numeric', isCircular='logical', 
                      genome='character')
  
  genome.root <- requestPkgDataPath('extdata', 'genomes')
  
  genomes <- list.dirs(genome.root, full.names=FALSE, recursive=FALSE)
  
  seqtab <- vector('list', length(genomes))
  names(seqtab) <- genomes
  
  for ( genome in genomes ) {
    
    filepath <- file.path(genome.root, genome, 'seqtab.csv')
    
    seqtab[[genome]] <- utils::read.csv(filepath, quote='',
                                        stringsAsFactors=FALSE, strip.white=TRUE, na.strings='', 
                                        colClasses=column.classes)
  }
  
  return(seqtab)
}

# loadVector -------------------------------------------------------------------
#' Load a vector from a line or file.
#' 
#' Given a specified line of text or filename, the input text is loaded
#' as a YAML list, and returned as an R \code{vector} object.
#' 
#' @param line Character string containing the line to be loaded.
#' @param file Character string containing the name of the file to be loaded.
#' 
#' @return Vector of loaded data.
#' 
#' @importFrom utils read.table
#' @importFrom yaml yaml.load
#' @importFrom yaml yaml.load_file
#' @keywords internal
#' @rdname loadVector
loadVector <- function(line=NULL, file=NULL, type=NULL) {
  
  if ( ! xor( is.null(line), is.null(file) ) ) {
    stop("vector must be loaded from either a line or a file, but not both")
  }
  
  if ( ! is.null(line) ) {
    
    stopifnot( isSingleString(line) )
    
    # Ensure line is enclosed by square brackets.
    first.char <- substr(line, 1, 1)
    last.char <- substr(line, nchar(line), nchar(line))
    if ( first.char != '[' || first.char != ']' ) {
      line <- paste0('[', line, ']', collapse='')
    }
    
    tryCatch({ # Load YAML list from line.
      x <- yaml::yaml.load(line)
      stopifnot( is.vector(x) && ! hasNames(x) )
    }, error=function(e) {
      stop("failed to load vector from line - '", toString(line), "'")
    })
    
  } else if ( ! is.null(file) ) {
    
    stopifnot( isSingleString(file) )
    
    tryCatch({
      
      # Read list file as a CSV file with one column.
      x <- utils::read.table(file, sep='\n', quote='', as.is=TRUE,
                             na.strings=NA, colClasses='character', strip.white=TRUE,
                             blank.lines.skip=FALSE, allowEscapes=TRUE,
                             stringsAsFactors=FALSE)
      stopifnot( ncol(x) == 1 )
      
      # Trim any blank rows from the bottom.
      x <- bstripBlankRows(x)
      
      # Convert to vector of YAML-loaded values.
      x <- lapply(getRowIndices(x), function(i) yaml::yaml.load(x[i, ]))
      stopifnot( all( lengths(x) == 1 ) )
      
      # Get data types of YAML-loaded elements.
      element.types <- unique( sapply(x, typeof) )
      
      if ( ! all( element.types %in%
                  c('character', 'logical', 'integer', 'double', 'numeric') ) ) {
        stop("unknown vector element types - '", toString(element.types), "'")
      }
      
      # If all elements are of the same type,
      # convert list to a vector of that type.
      if ( length(element.types) == 1 ) {
        
        if ( element.types == 'character' ) {
          x <- as.character(x)
        } else if ( element.types == 'logical' ) {
          x <- as.logical(x)
        } else if ( element.types == 'integer' ) {
          x <- as.integer(x)
        } else if ( element.types %in% c('double', 'numeric') ) {
          x <- as.numeric(x)
        }
      }
      
    }, error=function(e) {
      stop("failed to load vector from file - '", toString(file), "'")
    })
  }
  
  # If relevant, check type of vector against
  # that specified, and convert if appropriate.
  if ( ! is.null(type) ) {
    
    stopifnot( isSingleString(type) )
    
    err.type <- NULL
    if ( type == 'character' ) {
      
      if ( ! is.character(x) ) {
        x <- as.character(x)
      }
      
    } else if ( type == 'logical' ) {
      
      if ( ! is.logical(x) ) {
        err.type <- 'logical'
      }
      
    } else if ( type == 'integer' ) {
      
      if ( ! is.integer(x) ) {
        if ( is.numeric(x) ) {
          x <- as.integer(x)
        } else {
          err.type <- 'integer'
        }
      }
      
    } else if ( type %in% c('double', 'numeric') ) {
      
      if ( ! is.numeric(x) ) {
        err.type <- 'numeric'
      }
      
    } else if ( type == 'list' ) {
      
      if ( ! is.list(x) ) {
        x <- as.list(x)
      }
      
    } else {
      
      stop("unknown vector type - '", type, "'")
    }
    
    if ( ! is.null(err.type) ) {
      stop("vector is not of type - '", type, "'")
    }
  }
  
  return(x)
}

# makeAlleleSet ----------------------------------------------------------------
#' Make allele set from the given genotypes.
#' 
#' @param genotypes Sorted character vector of genotype symbols.
#' 
#' @return Character vector of allele symbols generated from the given genotypes.
#' 
#' @keywords internal
#' @rdname makeAlleleSet
makeAlleleSet <- function(genotypes) {
  validateGenotypeSet(genotypes)
  alleles <- sort( unique( unlist( strsplit(genotypes, '') ) ) )
  alleles <- alleles[ alleles != const$missing.value ]
  return(alleles)
}

# makeDefaultMarkerIDs ---------------------------------------------------------
#' Make default marker IDs for loci.
#' 
#' @description This function creates marker IDs from locus \code{data.frame}
#' \code{'loc'}, and validates that each locus has a map position that
#' is within range for the given reference sequence/chromosome. If any
#' locus position is found to be out-of-range, then it will stop with
#' an error message.
#' 
#' @param loc Locus \code{data.frame}, with columns \code{'chr'} and
#' \code{'pos'}, specifying physical map positions.
#' @param sep Separator between the two parts of a default marker ID.
#'      
#' @return Character vector of default marker IDs.
#' 
#' @template author-thomas-walsh
#' @template author-yue-hu
#' @template section-default-marker-ids
#' 
#' @export
#' @rdname makeDefaultMarkerIDs
makeDefaultMarkerIDs <- function(loc, sep=c(':', '-')) {
  
  stopifnot( is.data.frame(loc) )
  stopifnot( nrow(loc) > 0 )
  
  sep <- match.arg(sep)
  
  # Check map unit is physical, if present.
  map.unit <- getMapUnit(loc)
  if ( ! is.na(map.unit) ) {
    
    stopifnot( isPhysicalMapUnit(map.unit) )
    
    # Scale locus positions to base-pair units, if needed.
    basic.unit <- const$basic.map.unit[['pmap']]
    if ( map.unit != basic.unit ) {
      loc <- setMapUnit(loc, basic.unit)
    }
  }
  
  # Get locus sequences and positions.
  loc.seqs <- normSeq( pullLocusSeq(loc) )
  loc.pos <- pullLocusPos(loc)
  
  # Validate positions from genome sequence lengths.
  seqtab <- getSeqTable()
  loc.seq.lengths <- seqtab$seqlengths[ match(loc.seqs,seqtab$seqids) ]
  exrange <- loc.pos[ loc.pos < 1 | loc.pos > loc.seq.lengths ]
  if ( length(exrange) > 0 ) {
    stop("cannot make default marker IDs for positions '", toString(exrange), "'")
  }
  
  # Format locus positions.
  loc.pos <- sprintf('%07d', loc.pos)
  
  return( paste0('c', loc.seqs, sep, loc.pos) )
}

# makeDefaultQTLNames ----------------------------------------------------------
#' Make default QTL names for the specified loci.
#'  
#' @param loc Locus \code{data.frame}, with columns \code{'chr'} and
#' \code{'pos'}, specifying genetic map positions.
#' @param step Map step size.
#'      
#' @return Character vector of default QTL names.
#' 
#' @template section-qtl-names
#' 
#' @export
#' @rdname makeDefaultQTLNames
makeDefaultQTLNames <- function(loc, step=0) {
  
  stopifnot( is.data.frame(loc) )
  stopifnot( nrow(loc) > 0 )
  stopifnot( isSingleNonNegativeNumber(step) )
  
  # Check map unit is genetic, if present.
  map.unit <- getMapUnit(loc)
  if ( ! is.na(map.unit) ) {
    stopifnot( isGeneticMapUnit(map.unit) )
  }
  
  # Adjust number of digits in map position. 
  # NB: based on code in R/qtl::makeqtl
  digits <- 1
  if ( step > 0 ) {
    digits <- max( digits, -floor( log10(step) ) )
  }
  
  # Get locus sequences and positions.
  loc.seqs <- normSeq( pullLocusSeq(loc) )
  loc.pos <- pullLocusPos(loc)
  
  # Format locus positions.
  loc.pos <- sprintf(paste0('%.', digits, 'f'), loc.pos)
  
  return( paste(loc.seqs, loc.pos, sep='@') )
}

# makeGenotypeSet --------------------------------------------------------------
#' Make genotype set from the given alleles and cross type.
#' 
#' @param alleles Sorted character vector containing a set of allele symbols.
#' @param crosstype Cross type for which genotype set should be made.
#' @param strict Make genotype set containing only complete genotypes.
#' 
#' @return Character vector of genotype symbols generated from the
#' given alleles in a way determined by the specified cross type.
#' 
#' @keywords internal
#' @rdname makeGenotypeSet
makeGenotypeSet <- function(alleles, crosstype, strict=FALSE) {
  
  # TODO: implement for other cross types.
  # TODO: implement strict option. By default (strict=FALSE), partially
  # missing genotype symbols (e.g. 'A-') are included in the returned
  # genotype set. If strict is TRUE, the returned genotype set contains
  # only complete genotype symbols.
  
  validateAlleleSet(alleles)
  
  if ( crosstype == 'haploid' ) {
    
    # NB: strict option is irrelevant for a haploid cross type,
    # because partial genotypes cannot exist by definition.
    genotypes <- alleles
    
  } else if ( crosstype %in% const$supported.crosstypes ) {
    stop("cannot get genotype set for cross type - '", crosstype, "'")
  } else {
    stop("unsupported cross type - '", crosstype, "'")
  }
  
  return(genotypes)
}

# makeNumbers ------------------------------------------------------------------
#' Make numbers from numeric names.
#' 
#' @param x Vector of syntactically valid names that represent numbers.
#'      
#' @return Numeric vector of numbers corresponding to the input names.
#' 
#' @keywords internal
#' @rdname makeNumbers
makeNumbers <- function(x) {
  
  stopifnot( isValidName(x) )
  
  m <- regexec(const$pattern$numeric.name, x)
  matches <- regmatches(x, m)
  
  unconvertible <- x[ lengths(matches) == 0 ]
  if ( length(unconvertible) > 0 ) {
    stop("cannot make numbers from non-numeric names - '", 
         toString(unconvertible), "'")
  }
  
  is.negative <- sapply(matches, getElement, 2) == '.'
  numbers <- sapply(matches, function(regmatch) as.numeric(regmatch[3]))
  
  numbers[is.negative] <- -1 * numbers[is.negative]
  
  return(numbers)
}

# makePseudomarkerIDs ----------------------------------------------------------
#' Make pseudomarker IDs for the specified loci.
#'  
#' @param loc Locus \code{data.frame}, with columns \code{'chr'} and
#' \code{'pos'}, specifying genetic map positions.
#'      
#' @return Character vector of pseudomarker IDs.
#' 
#' @template section-pseudomarker-ids
#' 
#' @export
#' @rdname makePseudomarkerIDs
makePseudomarkerIDs <- function(loc) {
  
  stopifnot( is.data.frame(loc) )
  stopifnot( nrow(loc) > 0 )
  
  # Check map unit is genetic, if present.
  map.unit <- getMapUnit(loc)
  if ( ! is.na(map.unit) ) {
    stopifnot( isGeneticMapUnit(map.unit) )
  }
  
  # Get locus sequences and positions.
  loc.seqs <- normSeq( pullLocusSeq(loc) )
  loc.pos <- pullLocusPos(loc)
  
  return( paste0('c', loc.seqs, '.loc', loc.pos) )
}

# makeResultsOverview ----------------------------------------------------------
#' Make QTL analysis results overview.
#' 
#' @param phenotypes Vector of phenotype IDs.
#' @param analysis Name of analysis from which results were generated.
#' @param results Character vector summarising results of analysis. If the
#' results vector has names, results are matched to phenotypes by name.
#' 
#' @return A results overview \code{data.frame}, with two columns:
#' \code{'Phenotype'}, listing phenotypes for which the analysis was
#' done; and a second column summarising the results of that analysis.
#' 
#' @keywords internal
#' @rdname makeResultsOverview
makeResultsOverview <- function(phenotypes, analysis, results=NULL) {
  
  stopifnot( all( isValidID(phenotypes) ) )
  stopifnot( length(phenotypes) > 0 )
  stopif( anyDuplicated(phenotypes) )
  stopifnot( isSingleString(analysis) )
  stopifnot( analysis %in% names(const$supported.analyses) )
  
  phenotypes <- sort(phenotypes)
  
  if ( ! is.null(results) ) {
    
    stopifnot( is.character(results) )
    stopifnot( hasNames(results) )
    stopif( anyDuplicated( names(results) ) )
    stopifnot( setequal(names(results), phenotypes) )
    results <- unname(results[phenotypes])
    
  } else {
    
    results <- rep('', length(phenotypes))
  }
  
  columns <- structure(list(phenotypes, results),
                       names=c('Phenotype', analysis))
  
  overview <- as.data.frame(columns, stringsAsFactors=FALSE)
  
  return(overview)
}

# makeScanoneThresholdObject ---------------------------------------------------
#' Make \code{scanone} threshold object.
#' 
#' @param threshold LOD threshold value.
#' @param lodcolnames Character vector of LOD column names
#' to be used in the \code{scanone} threshold object.
#' @param alpha Significance level (alpha) associated with the specified LOD
#' threshold. (Incompatible with \code{fdr}.)
#' @param fdr False-discovery rate associated with the specified LOD
#' threshold. (Incompatible with \code{alpha}.)
#' 
#' @return If \code{alpha} is specified, a \code{summary.scanoneperm} object is
#' returned, with the given threshold at the specified \code{alpha} for each LOD
#' column. If \code{fdr} is specified, a \code{summary.scanonebins} object is
#' returned, with the given threshold at the specified \code{fdr} in the single
#' LOD column for all phenotypes. In any case, because no permutations are used
#' in generating the threshold, the returned object has attribute \code{'n.perm'}
#' containing the value \code{0}.
#' 
#' @keywords internal
#' @rdname makeScanoneThresholdObject
makeScanoneThresholdObject <- function(threshold, lodcolnames='lod', alpha=NULL,
                                       fdr=NULL) {
  
  stopifnot( isSingleNonNegativeNumber(threshold) )
  stopifnot( all( isValidName(lodcolnames) ) )
  stopifnot( length(lodcolnames) > 0 )
  
  if ( ! is.null(alpha) && ! is.null(fdr) ) {
    stop("scanone threshold cannot have both significance level (alpha) and FDR")
  } else if ( ! is.null(alpha) ) {
    
    stopifnot( isSingleProbability(alpha) )
    level <- paste0(as.character(100 * alpha), '%')
    threshold.class <- 'summary.scanoneperm'
    
  } else if ( ! is.null(fdr) ) {
    
    stopifnot( isSingleFiniteNumber(fdr) )
    stopifnot( fdr > 0 & fdr < 1 )
    level <- paste0(as.character(100 * fdr), '%')
    threshold.class <- 'summary.scanonebins'
    
  } else {
    stop("scanone threshold must have either significance level (alpha) or FDR")
  }
  
  if ( length(lodcolnames) == 1 && lodcolnames != 'lod' ) {
    warning("forcing single LOD column name to 'lod'")
    lodcolnames <- 'lod'
  }
  
  x <- matrix( rep(threshold, length(lodcolnames)), nrow=length(level),
               ncol=length(lodcolnames), dimnames=list(level, lodcolnames))
  class(x) <- threshold.class
  attr(x, 'n.perm') <- 0
  
  return(x)
}

# makeScantwoThresholdObject ---------------------------------------------------
#' Make \code{scantwo} threshold object.
#' 
#' @param thresholds Mapping or vector of \code{scantwo} LOD threshold values.
#' If a mapping or named vector, threshold values are accessed by name,
#' regardless of order. If an unnamed vector, threshold values are accessed
#' by position, and are assumed to be in the same order as those of a
#' \code{summary.scantwoperm} object.
#' @param phenotypes Character vector of phenotypes for which
#' the \code{scantwo} threshold object is being created.
#' @param alpha Significance level (alpha) associated
#' with the specified \code{scantwo} LOD thresholds.
#' 
#' @return A \code{summary.scantwoperm} object with the given \code{scantwo}
#' thresholds at the specified significance level (\code{alpha}) for each LOD
#' column. Because no permutations are used in generating the \code{scantwo}
#' thresholds object, the returned object has attribute \code{'n.perm'}
#' containing the value \code{0}.
#' 
#' Note that phenotype names are not included in the \code{summary.scantwoperm}
#' object. This is a feature, not a bug, and is intended to match the behaviour
#' of \pkg{R/qtl} in the case of a single significance level (\code{alpha}).
#' 
#' @keywords internal
#' @rdname makeScantwoThresholdObject
makeScantwoThresholdObject <- function(thresholds, phenotypes, alpha) {
  
  lodtypes <- const$scantwo.lodtypes$scantwoperm
  
  stopifnot( all( isValidID(phenotypes) ) )
  stopifnot( length(phenotypes) > 0 )
  stopifnot( isSingleProbability(alpha) )
  
  if ( ! is.mapping(thresholds) ) {
    if ( hasNames(thresholds) ) {
      thresholds <- as.mapping(thresholds)
    } else {
      thresholds <- mapping(values=thresholds, keys=lodtypes)
    }
  }
  
  invalid.names <- setdiff(names(thresholds), lodtypes)
  if ( length(invalid.names) > 0 ) {
    stop("invalid scantwo LOD threshold names - '", toString(invalid.names), "'")
  }
  
  missing.names <- setdiff(lodtypes, names(thresholds))
  if ( length(missing.names) > 0 ) {
    stop("missing scantwo LOD thresholds - '", toString(missing.names), "'")
  }
  
  non.numeric <- unlist( thresholds[ ! sapply(thresholds, is.numeric) ] )
  if ( length(non.numeric) > 0 ) {
    stop("non-numeric scantwo thresholds - '", toString(non.numeric), "'")
  }
  
  invalid.values <- unlist( thresholds[ ! sapply(thresholds, isNonNegativeNumber) ] )
  if ( length(invalid.values) > 0 ) {
    stop("invalid scantwo thresholds - '", toString(invalid.values), "'")
  }
  
  thresholds <- sapply(lodtypes, function(k) thresholds[[k]])
  
  level <- paste0(as.character(100 * alpha), '%')
  
  x <- lapply(lodtypes, function(k) matrix(
    rep(thresholds[k], length(phenotypes)), nrow=length(level),
    ncol=length(phenotypes), dimnames=list(level, NULL) ) )
  names(x) <- lodtypes
  
  class(x) <- c('summary.scantwoperm', 'list')
  
  attr(x, 'n.perm') <- 0
  
  return(x)
}

# otherattributes --------------------------------------------------------------
#' Get non-reserved object attributes. 
#' 
#' @details As with R base function \code{mostattributes}, this provides access
#' to non-reserved attributes.
#' 
#' @param x R object.
#' 
#' @return List of non-reserved object attributes.
#' 
#' @keywords internal
#' @rdname otherattributes
otherattributes <- function(x) {
  special.attributes <- getSpecialAttributeNames(x)
  return( attributes(x)[ ! names(attributes(x)) %in% special.attributes ] )
}

# `otherattributes<-` ----------------------------------------------------------
#' Set non-reserved object attributes. 
#' 
#' @details As with R base function \code{mostattributes}, this provides access
#' to non-reserved attributes.
#' 
#' @param x R object.
#' @param value List of attributes to set.
#' 
#' @return R object with non-reserved attributes set.
#' 
#' @keywords internal
#' @rdname otherattributes
`otherattributes<-` <- function(x, value) {
  
  stopifnot( hasNames(value) )
  
  special.attributes <- getSpecialAttributeNames(x)
  
  other.names <- names(value)[ ! names(value) %in% special.attributes ]
  for ( other.name in other.names ) {
    attr(x, other.name) <- value[[other.name]]
  }
  
  return(x)
}

# parseDefaultMarkerIDs --------------------------------------------------------
#' Parse default marker IDs.
#' 
#' This function parses an input vector of default marker IDs, and returns a
#' physical \code{mapframe} with the locus in each row derived from the
#' corresponding marker ID in the input vector. An error is raised if any
#' of the input values cannot be parsed as a default marker ID.
#'  
#' @param marker.ids Vector of default marker IDs.
#' 
#' @return Physical \code{mapframe} with loci
#' corresponding to the specified marker IDs.
#' 
#' @template section-default-marker-ids
#' 
#' @export
#' @rdname parseDefaultMarkerIDs
parseDefaultMarkerIDs <- function(marker.ids) {
  
  stopifnot( all( isDefaultMarkerID(marker.ids) ) )
  
  m <- regexec(const$pattern$default.marker.id, marker.ids)
  regmatch.list <- regmatches(marker.ids, m)
  
  marker.seqs <- sapply(regmatch.list, getElement, 2)
  marker.pos <- sapply(regmatch.list, function(x) as.numeric(x[3]))
  
  return( mapframe(chr=marker.seqs, pos=marker.pos,
                   row.names=marker.ids, map.unit='bp') )
}

# parseDefaultQTLNames ---------------------------------------------------------
#' Parse default QTL names.
#' 
#' This function parses an input vector of default QTL names, and returns a
#' genetic \code{mapframe} with the locus in each row derived from the
#' corresponding QTL name in the input vector. An error is raised if any
#' of the input values cannot be parsed as a default QTL name.
#' 
#' @param qtl.names Vector of default QTL names.
#' 
#' @return Genetic \code{mapframe} with loci
#' corresponding to the specified QTL names.
#' 
#' @template section-qtl-names
#' 
#' @export
#' @rdname parseDefaultQTLNames
parseDefaultQTLNames <- function(qtl.names) {
  
  stopifnot( all( isDefaultQTLName(qtl.names) ) )
  
  m <- regexec(const$pattern$default.qtl.name, qtl.names)
  regmatch.list <- regmatches(qtl.names, m)
  
  qtl.seqs <- sapply(regmatch.list, getElement, 2)
  qtl.pos <- sapply(regmatch.list, function(x) as.numeric(x[3]))
  
  return( gmapframe(chr=qtl.seqs, pos=qtl.pos, row.names=qtl.names) )
}

# parseFilenames ---------------------------------------------------------------
#' Parse filenames by pattern matching.
#' 
#' @param filenames Character vector of filenames to parse.
#' @param filename.pattern Filename pattern, which must be a valid Perl regex
#' with named capture groups. Neither the capture groups nor the regex itself
#' are required to match any given filename, but all capture groups must have
#' a unique name.
#' 
#' @return Character matrix in which each row contains the values of capture
#' groups for a given filename, and each column contains values of a given
#' capture group across files. Unmatched capture groups are represented by
#' \code{NA} values.
#' 
#' @keywords internal
#' @rdname parseFilenames
parseFilenames <- function(filenames, pattern) {
  
  stopifnot( is.character(filenames) )
  stopifnot( isSingleString(pattern) )
  
  parsed <- as.list( structure(rep(NA_character_,
                                   length(filenames)), names=filenames) )
  
  capture.names <- NULL
  
  # Parse each filename.
  for ( filename in filenames ) {
    
    # Apply regular expression to filename.
    pattern.match <- regexpr(pattern, filename, perl=TRUE)
    
    # Get capture-group names of pattern match.
    capture.names <- attr(pattern.match, 'capture.names')
    stopifnot( length(capture.names) > 0 )
    stopif( any( capture.names == '' ) )
    
    # Get first indices of capture groups.
    capture.first <- attr(pattern.match, 'capture.start')[1, ]
    capture.first[ capture.first %in% c(-1, 0) ] <- NA_integer_
    
    # Get lengths of capture groups.
    capture.length <- attr(pattern.match, 'capture.length')[1, ]
    capture.length[ capture.length %in% c(-1, 0) ] <- NA_integer_
    
    # Get last indices of capture groups.
    capture.last <- capture.first + capture.length - 1
    
    # Get capture-group strings.
    parsed[[filename]] <- substring(filename, capture.first, capture.last)
    names(parsed[[filename]]) <- attr(pattern.match, 'capture.names')
  }
  
  parsed <- matrix( unlist(parsed), nrow=length(parsed),
                    byrow=TRUE, dimnames=list(filenames, capture.names) )
  
  return(parsed)
}

# parsePseudomarkerIDs ---------------------------------------------------------
#' Parse pseudomarker IDs.
#' 
#' This function parses an input vector of pseudomarker IDs, and returns a
#' genetic \code{mapframe} with the locus in each row derived from the
#' corresponding pseudomarker ID in the input vector. An error is raised
#' if any of the input values cannot be parsed as a pseudomarker ID.
#' 
#' @param loc.ids Vector of pseudomarker IDs.
#' 
#' @return Genetic \code{mapframe} with loci
#' corresponding to the input pseudomarker IDs.
#' 
#' @template section-pseudomarker-ids
#' 
#' @export
#' @rdname parsePseudomarkerIDs
parsePseudomarkerIDs <- function(loc.ids) {
  
  stopifnot( all( isPseudomarkerID(loc.ids) ) )
  
  m <- regexec(const$pattern$pseudomarker.id, loc.ids)
  regmatch.list <- regmatches(loc.ids, m)
  
  loc.seqs <- sapply(regmatch.list, getElement, 2)
  loc.pos <- sapply(regmatch.list, function(x) as.numeric(x[3]))
  
  return( gmapframe(chr=loc.seqs, pos=loc.pos, row.names=loc.ids) )
}

# pull.alleles -----------------------------------------------------------------
#' Pull alleles from the given object.
#' 
#' @param x An \pkg{R/qtl} \code{cross} object,
#' or other object with alleles.
#' 
#' @return Vector of alleles in the given object.
#' 
#' @export
#' @family cross object functions
#' @rdname pull.alleles
pull.alleles <- function(x) {
  UseMethod('pull.alleles', x)
}

# pull.alleles.cross -----------------------------------------------------------
#' @export
#' @method pull.alleles cross
#' @rdname pull.alleles
pull.alleles.cross <- function(x) {
  return( attr(x, 'alleles') )
}

# pull.alleles.geno ------------------------------------------------------------
#' @export
#' @method pull.alleles geno
#' @rdname pull.alleles
pull.alleles.geno <- function(x) {
  return( attr(x, 'alleles') )
}

# pull.chr ---------------------------------------------------------------------
#' Pull chromosomes/sequences from the given object.
#' 
#' @param x An \pkg{R/qtl} \code{cross} object,
#' or other object with chromosomes/sequences.
#' 
#' @return Vector of chromosomes/sequences in the given object.
#' 
#' @export
#' @family cross object functions
#' @rdname pull.chr
pull.chr <- function(x) {
  UseMethod('pull.chr', x)
}

# pull.chr.cross ---------------------------------------------------------------
#' @export
#' @method pull.chr cross
#' @rdname pull.chr
pull.chr.cross <- function(x) {
  return( names(x$geno) )
}

# pull.chr.data.frame ----------------------------------------------------------
#' @export
#' @method pull.chr data.frame
#' @rdname pull.chr
pull.chr.data.frame <- function(x) {
  return( unique( pullLocusSeq(x) ) )
}

# pull.chr.geno ----------------------------------------------------------------
#' @export
#' @method pull.chr geno
#' @rdname pull.chr
pull.chr.geno <- function(x) {
  return( names(x) )
}

# pull.chr.list ----------------------------------------------------------------
#' @export
#' @method pull.chr list
#' @rdname pull.chr
pull.chr.list <- function(x) {
  return( unique( pullLocusSeq(x) ) )
}

# pull.chr.map -----------------------------------------------------------------
#' @export
#' @method pull.chr map
#' @rdname pull.chr
pull.chr.map <- function(x) {
  return( names(x) )
}

# pull.crosstype ---------------------------------------------------------------
#' Pull cross type from the given object.
#' 
#' @param x An \pkg{R/qtl} \code{cross} object,
#' or other object with cross type information.
#' 
#' @return Cross type of the given object.
#' 
#' @export
#' @family cross object functions
#' @rdname pull.crosstype
pull.crosstype <- function(x) {
  UseMethod('pull.crosstype', x)
}

# pull.crosstype.cross ---------------------------------------------------------
#' @export
#' @method pull.crosstype cross
#' @rdname pull.crosstype
pull.crosstype.cross <- function(x) {
  return( class(x)[1] )
}

# pull.crosstype.geno ----------------------------------------------------------
#' @export
#' @method pull.crosstype geno
#' @rdname pull.crosstype
pull.crosstype.geno <- function(x) {
  
  cross.info <- attr(x, 'info')
  
  if ( ! is.null(cross.info) ) {
    crosstype <- cross.info@crosstype
  } else {
    crosstype <- NULL
  }
  
  return(crosstype)
}

# pull.ind ---------------------------------------------------------------------
#' Pull individual sample IDs from the given object.
#' 
#' @param x An \pkg{R/qtl} \code{cross} object,
#' or other object with sample ID information.
#' 
#' @return Vector of individual sample IDs in the given object.
#' Returns \code{NULL} if no sample IDs are found.
#' 
#' @export
#' @family cross object functions
#' @rdname pull.ind
pull.ind <- function(x) {
  UseMethod('pull.ind', x)
}

# pull.ind.cross ---------------------------------------------------------------
#' @export
#' @method pull.ind cross
#' @rdname pull.ind
pull.ind.cross <- function(x) {
  
  cross.info <- attr(x, 'info')
  
  if ( ! is.null(cross.info) && hasSampleIDs(cross.info) ) {
    
    ind <- getSamples(cross.info)
    
  } else {
    
    id.col <- getIdColIndex(x)
    
    if ( ! is.null(id.col) ) {
      ind <- as.character(x$pheno[[id.col]])
    } else {
      ind <- NULL
    }
  }
  
  return(ind)
}

# pull.ind.geno ----------------------------------------------------------------
#' @export
#' @method pull.ind geno
#' @rdname pull.ind
pull.ind.geno <- function(x) {
  
  cross.info <- attr(x, 'info')
  
  if ( ! is.null(cross.info) && hasSampleIDs(cross.info) ) {
    ind <- getSamples(cross.info)
  } else {
    ind <- NULL
  }
  
  return(ind)
}

# pull.ind.pheno ---------------------------------------------------------------
#' @export
#' @method pull.ind pheno
#' @rdname pull.ind
pull.ind.pheno <- function(x) {
  
  cross.info <- attr(x, 'info')
  
  if ( ! is.null(cross.info) && hasSampleIDs(cross.info) ) {
    
    ind <- getSamples(cross.info)
    
  } else {
    
    id.col <- getIdColIndex(x)
    
    if ( ! is.null(id.col) ) {
      ind <- as.character(x[[id.col]])
    } else {
      ind <- NULL
    }
  }
  
  return(ind)
}

# removeColsNA -----------------------------------------------------------------
#' Remove columns containing only \code{NA} values.
#'
#' @param x A \code{data.frame} or \code{matrix}.
#'
#' @return Input object without those columns that contain only \code{NA} values.
#'
#' @keywords internal
#' @rdname removeColsNA
removeColsNA <- function(x) {
  stopifnot( is.data.frame(x) || is.matrix(x) )
  stopifnot( nrow(x) > 0 )
  mask <- ! apply( x, 2, function(column) allNA( as.vector(column) ) )
  x <- x[, mask, drop=FALSE]
  return(x)
}

# removeRowsNA -----------------------------------------------------------------
#' Remove rows containing only \code{NA} values.
#'
#' @param x A \code{data.frame} or \code{matrix}.
#'
#' @return Input object without those rows that contain only \code{NA} values.
#'
#' @keywords internal
#' @rdname removeRowsNA
removeRowsNA <- function(x) {
  stopifnot( is.data.frame(x) || is.matrix(x) )
  stopifnot( nrow(x) > 0 )
  mask <- ! apply( x, 1, function(row) allNA( as.vector(row) ) )
  x <- x[mask,, drop=FALSE]
  return(x)
}

# requestNodes -----------------------------------------------------------------
#' Request nodes for parallel execution.
#' 
#' @template param-n.cluster
#'     
#' @return Vector of node names.
#' 
#' @export
#' @rdname requestNodes
requestNodes <- function(n.cluster=1) {
  
  stopifnot( isSinglePositiveWholeNumber(n.cluster) )
  
  # Get the PBS node file, or NA value if not found.
  node.file <- Sys.getenv('PBS_NODEFILE', NA)
  
  # If PBS node file specified, create node list from node file..
  if ( ! is.na(node.file) && file.info(node.file)$size > 0 ) {
    
    # Read lines from node file.
    f <- file(node.file, 'r')
    node.lines <- readLines(f)
    close(f)
    
    # Create node list with N nodes, where N is the number of processes  
    # requested or the number of nodes in the file, whichever is smaller.
    node.list <- node.lines[1:n.cluster]
    node.list <- node.list[ ! is.na(node.list) ]
    
  } else { # ..otherwise create node list to run on local host. 
    
    # Get number of cores available on the local host.
    local.cores <- parallel::detectCores(logical=FALSE)
    
    # Cap the number of processes at the number of cores.
    n.cluster <- min(n.cluster, local.cores)
    
    # Set node list from number of processes.
    node.list <- rep('localhost', n.cluster)
  }
  
  return(node.list)
}

# requestPkgDataPath -----------------------------------------------------------
#' Request path of \pkg{shmootl} package data.
#' 
#' @param ... Path elements. These should be specified relative to the package 
#' data root. For the source package, this is the \pkg{shmootl} \code{'inst'} 
#' directory, while for the binary package, this is the same as the package 
#' root.
#' 
#' @return An absolute path to the requested \pkg{shmootl} data file/directory.
#' 
#' @keywords internal
#' @rdname requestPkgDataPath
requestPkgDataPath <- function(...) {
  
  datapath <- file.path(...)
  
  result <- system.file(datapath, package='shmootl')
  
  if ( result == '' ) {
    
    # NB: this assumes the current working directory is 'shmootl/R'.
    result <- normalizePath( file.path('..', 'inst', datapath), mustWork=TRUE )
    
    if ( ! file.exists(result) ) {
      stop("shmootl package data file/directory not found - '", 
           datapath, "'")
    }
  }
  
  return(result)
}

# resolveScantwoLodtypes -------------------------------------------------------
#' Resolve \code{scantwo} LOD types.
#' 
#' @param lodtypes Vector of \code{scantwo} LOD types. If this parameter is not
#' specified, or is \code{NULL}, a vector of all \code{scantwo} LOD types is
#' returned.
#' @param from One or more forms from which LOD types may be resolved.
#' @param to Form to which LOD types should be resolved.
#' 
#' @return Character vector with \code{scantwo} LOD types
#' resolved to the form specified by the \code{to} parameter.
#' 
#' @keywords internal
#' @rdname resolveScantwoLodtypes
resolveScantwoLodtypes <- function(lodtypes=NULL,
                                   from=c('scantwoperm', 'plot.scantwo', 'summary.scantwo', 'title'),
                                   to=c('scantwoperm', 'plot.scantwo', 'summary.scantwo', 'title') ) {
  
  from <- match.arg(from, several.ok=TRUE)
  to <- match.arg(to)
  
  if ( ! is.null(lodtypes) ) {
    
    stopifnot( is.character(lodtypes) )
    stopifnot( length(lodtypes) > 0 )
    
    # Get indices of rows containing a match for each specified LOD type.
    row.indices <- getRowIndices(const$scantwo.lodtypes)
    row.index.list <- lapply( lodtypes, function(lodtype)
      which( sapply( row.indices, function(row.index)
        lodtype %in% const$scantwo.lodtypes[row.index, from] ) ) )
    
    # Get number of rows matched by each specified LOD type.
    match.counts <- lengths(row.index.list)
    
    unresolved <- lodtypes[ match.counts == 0 ]
    if ( length(unresolved) > 0 ) {
      stop("cannot resolve '", to, "' form of LOD types - '",
           toString(unresolved), "'")
    }
    
    ambiguous <- lodtypes[ match.counts > 1 ]
    if ( length(ambiguous) > 0 ) {
      stop("cannot unambiguously resolve '", to, "' form of LOD types - '",
           toString(ambiguous), "'")
    }
    
    lodtypes <- const$scantwo.lodtypes[unlist(row.index.list), to]
    
  } else {
    
    lodtypes <- const$scantwo.lodtypes
  }
  
  return(lodtypes)
}

# rstripBlankCols --------------------------------------------------------------
#' Strip blank columns from right of \code{data.frame}.
#'
#' @param x A \code{data.frame} with columns of type \code{character}.
#'
#' @return Input object in which rightmost blank columns (i.e. columns
#' containing no non-whitespace characters) have been stripped.
#' 
#' @keywords internal
#' @rdname rstripBlankCols
rstripBlankCols <- function(x) {
  stopifnot( is.data.frame(x) )
  stopifnot( all( sapply(x, class) == 'character' ) )
  while( allWhite( as.character(x[, ncol(x)]) ) ) {
    x <- x[, -ncol(x), drop=FALSE]
  }
  return(x)
}

# setColumnFromRownames --------------------------------------------------------
#' Move object row names to the specified column.
#' 
#' @param x A \code{data.frame} or \code{matrix}.
#' @param col.name Name of the column to which row names are to be set.
#'     
#' @return Input object with row names inserted into the leftmost column, and 
#' original row names removed.
#' 
#' @keywords internal
#' @rdname setColumnFromRownames
setColumnFromRownames <- function(x, col.name='rownames') {
  
  stopifnot( is.data.frame(x) || is.matrix(x) )
  stopifnot( hasRownames(x) )
  stopif( col.name %in% colnames(x) )
  
  x <- insertColumn(x, 1, col.name=col.name, data=rownames(x))
  
  rownames(x) <- NULL
  
  return(x)
}

# setRownamesFromColumn --------------------------------------------------------
#' Move specified column to object row names.
#' 
#' @param x A \code{data.frame} or \code{matrix}.
#' @param col.name Name of the column from which row names are to be set.
#'     
#' @return Input object with the leftmost column removed, and with row names set
#' from its contents.
#' 
#' @keywords internal
#' @rdname setRownamesFromColumn
setRownamesFromColumn <- function(x, col.name='rownames') {
  
  stopifnot( is.data.frame(x) || is.matrix(x) )
  stopif( anyDuplicated( colnames(x) ) )
  stopifnot( col.name %in% colnames(x) )
  
  rownames(x) <- as.character( x[, col.name] )
  
  x <- deleteColumn(x, col.name=col.name)
  
  return(x)
}

# stopif -----------------------------------------------------------------------
#' Stop if the expression is \code{TRUE}.
#' 
#' @details This function is based on the R base function \code{stopifnot}, but  
#' is defined only for the simple case of a single expression. It should be used
#' to avoid confusing double negatives.
#'  
#' @param expression R expression to be tested.
#' 
#' @keywords internal
#' @rdname stopif
stopif <- function(expression) {
  
  if ( length(expression) > 0 ) {
    
    mc <- match.call()
    
    ch <- deparse(mc[[2]])
    
    if ( ! anyNA(expression) && all(expression) ) {
      
      if ( length(expression) > 1 ) {
        msg <- sprintf("%s are all TRUE", ch)
      } else {
        msg <- sprintf("%s is TRUE", ch)
      }
      
      stop(msg, call. = FALSE)
    }    
  }
  
  return( invisible() )
}

# stripWhite -------------------------------------------------------------------
#' Strip leading/trailing whitespace from elements of a character vector.
#' 
#' @param x Character vector whose elements are to be stripped.
#'     
#' @return Whitespace-stripped copy of the input character vector.
#' 
#' @keywords internal
#' @rdname stripWhite
stripWhite <- function(x) {
  stopifnot( is.character(x) )
  return( gsub( '^[[:space:]]+|[[:space:]]+$', '', x) )
}

# updateResultsOverview --------------------------------------------------------
#' Update QTL analysis results overview.
#' 
#' @param overview A \code{data.frame} containing
#' a QTL analysis results overview.
#' @param analysis Name of analysis from which results were generated.
#' @param results Named character vector summarising results of analysis.
#' Each result is matched to its corresponding phenotype by name.
#' 
#' @return A results overview \code{data.frame}, with two or more columns:
#' \code{'Phenotype'}, listing phenotypes for which analyses were done; and one
#' or more additional columns, each summarising the results of a QTL analysis.
#' 
#' @keywords internal
#' @rdname updateResultsOverview
updateResultsOverview <- function(overview, analysis, results=NULL) {
  
  supported.analyses <- names(const$supported.analyses)
  
  validateResultsOverview(overview)
  stopifnot( isSingleString(analysis) )
  stopifnot( analysis %in% supported.analyses )
  
  input.headings <- colnames(overview)
  input.analyses <- input.headings[-1]
  
  output.analyses <- supported.analyses[ supported.analyses %in%
                                           c(input.analyses, analysis) ]
  output.headings <- c('Phenotype', output.analyses)
  
  input.phenotypes <- overview[, 'Phenotype']
  
  if ( hasNames(results) ) {
    output.phenotypes <- sort( union( input.phenotypes, names(results) ) )
  } else {
    output.phenotypes <- sort(input.phenotypes)
  }
  
  columns <- as.list(overview)
  
  for ( output.analysis in output.analyses ) {
    
    output.results <- structure(rep('', length(output.phenotypes)),
                                names=output.phenotypes)
    
    if ( output.analysis %in% input.analyses ) {
      
      input.results <- structure(columns[[output.analysis]],
                                 names=input.phenotypes)
      
      for ( phenotype in names(input.results) ) {
        output.results[phenotype] <- input.results[phenotype]
      }
    }
    
    if ( output.analysis == analysis && ! is.null(results) ) {
      
      stopifnot( is.character(results) )
      
      if ( hasNames(results) ) {
        stopif( anyDuplicated( names(results) ) )
      } else {
        stopifnot( length(results) == length(input.phenotypes) )
        names(results) <- input.phenotypes
      }
      
      for ( phenotype in names(results) ) {
        output.results[phenotype] <- results[phenotype]
      }
    }
    
    columns[[output.analysis]] <- unname(output.results)
  }
  
  # Update phenotype column.
  columns[['Phenotype']] <- output.phenotypes
  
  # Order column list by output headings.
  columns <- columns[output.headings]
  
  # Create updated results overview.
  overview <- as.data.frame(columns, stringsAsFactors=FALSE)
  
  return(overview)
}

# validateAlleleSet ------------------------------------------------------------
#' Validate a set of alleles.
#' 
#' @param x Vector of alleles.
#' 
#' @return \code{TRUE} if allele set is valid; otherwise raises first error.
#' 
#' @keywords internal
#' @rdname validateAlleleSet
validateAlleleSet <- function(x) {
  
  stopifnot( is.character(x) )
  
  # Decompose allele symbols into different types.
  founder.symbols <- x[ isFounderAllele(x) ]
  enum.symbols <- x[ isEnumAllele(x) ]
  invalid.values <- x[ is.na(x) | ! isValidAllele(x) ]
  
  # Check for invalid values.
  if ( length(invalid.values) > 0 ) {
    if ( '' %in% invalid.values ) {
      stop("blank allele values")
    } else {
      stop("invalid allele symbols - '", toString(invalid.values), "'")
    }
  }
  
  # Check for duplicate allele symbols.
  dup.symbols <- unique( x[ duplicated(x) ] )
  if ( length(dup.symbols) > 0 ) {
    stop("duplicate allele symbols - '", toString(dup.symbols), "'")
  }
  
  # Check allele set order.
  if ( is.unsorted(x) ) {
    stop("allele set is not sorted - '", toString(x), "'")
  }
  
  # Check that symbols are either founder or enumerated alleles.
  if ( ! xor(length(founder.symbols) > 0, length(enum.symbols) > 0) ) {
    stop("alleles must be of enumerated or founder type, but not both")
  }
  
  # Verify that there are exactly two alleles.
  # TODO: handle more than two alleles.
  if ( length(x) != 2 ) {
    stop("unsupported number of alleles - '", length(x), "'")
  }
  
  return(TRUE)
}

# validateGenotypeSet ----------------------------------------------------------
#' Validate a set of genotypes.
#' 
#' @param x Vector of genotypes.
#' @param strict Restrict validity to complete genotypes only.
#' 
#' @return \code{TRUE} if genotype set is valid; otherwise raises first error.
#' 
#' @keywords internal
#' @rdname validateGenotypeSet
validateGenotypeSet <- function(x, strict=FALSE) {
  
  stopifnot( is.character(x) )
  
  # Decompose genotype symbols into different types.
  founder.symbols <- x[ isFounderGenotype(x, strict=strict) ]
  enum.symbols <- x[ isEnumGenotype(x) ]
  invalid.values <- x[ is.na(x) | ! isValidGenotype(x, strict=strict) ]
  
  # Check for invalid values.
  if ( length(invalid.values) > 0 ) {
    if ( '' %in% invalid.values ) {
      stop("blank genotype values")
    } else {
      stop("invalid genotype symbols - '", toString(invalid.values), "'")
    }
  }
  
  # Check that symbols are either founder or enumerated genotypes.
  if ( ! xor(length(founder.symbols) > 0, length(enum.symbols) > 0) ) {
    stop("genotypes must be of enumerated or founder type, but not both")
  }
  
  # Check for duplicate genotype symbols.
  dup.symbols <- unique( x[ duplicated(x) ] )
  if ( length(dup.symbols) > 0 ) {
    stop("duplicate genotype symbols - '", toString(dup.symbols), "'")
  }
  
  # Check genotype set order.
  if ( is.unsorted(x) ) {
    stop("genotype set is not sorted - '", toString(x), "'")
  }
  
  # Verify that genotypes are haploid.
  # TODO: handle other ploidies.
  if ( any( nchar(x) > 1 ) ) {
    stop("unsupported genotype ploidy")
  }
  
  # Verify that genotyes have consistent ploidy.
  if ( length( unique( nchar(x) ) ) != 1 ) {
    stop("inconsistent genotype ploidy")
  }
  
  # Verify that there are exactly two genotypes.
  # TODO: handle more than two genotypes.
  if ( length(x) != 2 ) {
    stop("unsupported number of genotypes - '", length(x), "'")
  }
  
  return(TRUE)
}

# validateResultsOverview ------------------------------------------------------
#' Validate a QTL analysis results overview.
#' 
#' @param overview Results overview \code{data.frame}.
#' 
#' @return \code{TRUE} if results overview is valid;
#' otherwise raises first error.
#' 
#' @keywords internal
#' @rdname validateResultsOverview
validateResultsOverview <- function(overview) {
  
  stopifnot( is.data.frame(overview) )
  stopifnot( all( sapply(overview, class) == 'character' ) )
  stopifnot(  ncol(overview) >= 2 )
  
  headings <- colnames(overview)
  if ( headings[1] != 'Phenotype' ) {
    stop("first column of results overview must be 'Phenotype' column")
  }
  
  analyses <- headings[-1]
  
  unknown <- analyses[ ! analyses %in% names(const$supported.analyses) ]
  if ( length(unknown) > 0 ) {
    stop("results overview contains unknown analyses - '", toString(unknown), "'")
  }
  
  phenotypes <- overview[, 'Phenotype']
  
  invalid <- phenotypes[ ! isValidID(phenotypes) ]
  if ( length(invalid) > 0 ) {
    stop("results overview phenotype IDs are invalid - '", toString(invalid), "'")
  }
  
  if ( is.unsorted(phenotypes, strictly=TRUE) ) { # NB: also ensures no duplicates
    stop("results overview phenotypes are out of order")
  }
  
  return(TRUE)
}

# End of util.R ################################################################

