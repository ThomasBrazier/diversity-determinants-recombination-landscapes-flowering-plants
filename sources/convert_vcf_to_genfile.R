# Start of convert_vcf_to_genfile.R

# convert_vcf_to_genfile
#' Convert VCF data to an \pkg{R/qtl} genotype file.
#'
#' Read SNP genotype data from one or more VCF \code{infiles}, and output
#' these to a \code{genfile} that may be accepted as input by \pkg{R/qtl}.
#'
#' @section What does this function do?:
#'
#' Given one or more VCF \code{infiles}, along with a set of identifiers for \code{samples} and
#' \code{founders} from a two-founder cross, this function first reads SNP allele and genotype
#' calls. Then, for the set of SNPs that are common to samples and founders, where the base
#' call of a sample allele matches that of a founder, the sample is assigned the corresponding
#' founder allele. The resulting sample genotype is made up of the combination of its founder
#' alleles. A SNP is retained only if each founder has a distinct homozygous genotype, and there
#' are at least two distinct genotypes among the samples for that SNP. Finally, the set of SNP
#' markers satisfying these conditions are output to \code{genfile}.
#'
#' @section What does this function \emph{not} do?:
#'
#' This function does not convert genotypes for a multi-founder cross.
#'
#' Variants are taken at face value, and there is no consideration of variant or genotype quality.
#' To exclude low-quality variants and genotypes, you may wish to perform a hard filter beforehand
#' using a suitable VCF or BCF toolkit.
#'
#' No attempt is made to collapse redundant markers. Depending on SNP density and linkage block
#' size, the output \code{genfile} may contain groups of linked SNPs that could be collapsed to
#' a single marker. If appropriate, this can be done with the data in the \code{genfile}.
#'
#' While care has been taken to ensure that the output \code{genfile} accurately reflects
#' the genotype data in the input VCF file(s), the genotype data is not validated with
#' respect to a specific cross type or experimental design. In general, you should review
#' the output file to ensure that the results are appropriate for your dataset.
#'
#' @section SNP marker IDs:
#' Each SNP marker is assigned an identifier of the form \code{'chr11:002160000'}, where
#' \code{'chr11'} is the chromosome identifier taken from the VCF \code{CHROM} field, and
#' \code{'002160000'} is a zero-padded number giving the position of the given SNP in the
#' reference genome, taken from the VCF \code{POS} field.
#'
#' The width to which the genomic position is padded is determined by the maximum sequence
#' length. Ideally reference sequence lengths are already specified in contig fields in the
#' header of the VCF. In that case the maximum sequence length is taken from the input VCF
#' file, and there is no need to specify a \code{max.seqlength} parameter.
#'
#' If the VCF header does not contain the necessary reference sequence length information,
#' then the \code{max.seqlength} parameter should be set to the appropriate value for the
#' given reference genome to ensure consistent formatting of SNP marker IDs across datasets.
#'
#' Sequence length information may be obtained by using one of the \code{ChromInfo} functions
#' provided by the \pkg{GenomeInfoDb} package. Alternatively, it may be obtained from a Picard
#' Tools sequence dictionary file using the function \code{\link{load_seq_dict}}, or directly
#' from a reference genome FASTA file using the \pkg{Biostrings} function \code{fasta.seqlengths}.
#' For more, see the \strong{Examples} section below.
#'
#' @param infiles Input VCF file paths.
#' @param genfile Output \pkg{R/qtl} genotype data file.
#' @param samples Cross sample IDs.
#' @param founders Founder sample IDs.
#' @param alleles Optional mapping of founder IDs to allele symbols (e.g. \code{utl::mapping(
#' c(DBVPG6044 = 'W', Y12 = 'S') )}). If this parameter is not specified, allele symbols are
#' taken from the letters of the alphabet (i.e. \code{'A'}, \code{'B'} etc.).
#' @param max.seqlength Optional parameter to indicate the maximum reference sequence length, which
#' is used to determine the zero-padded width of genomic positions in SNP marker IDs. Without this
#' information, SNP marker IDs may be formatted inconsistently in different datasets. For more
#' details, see the \strong{SNP marker IDs} section.
#' @param na.string String to replace \code{NA} values.
#'
#' @examples
#' \dontrun{
#'
#' # To get maximum sequence length from a sequence dictionary.
#' seqinfo <- utl::load_seq_dict('reference_genome.dict')
#' max.seqlength <- max(GenomeInfoDb::seqlengths(seqinfo))
#'
#' # To get maximum sequence length from reference sequence FASTA file.
#' max.seqlength <- max(Biostrings::fasta.seqlengths('reference_genome.fasta'))
#'
#' # Convert VCF to R/qtl genfile.
#' infiles <- c('samples.vcf', 'founders.vcf')
#' genfile <- 'geno.csv'
#' sample.ids <- utl::read_samples_from_vcf('samples.vcf')
#' founder.ids <- utl::read_samples_from_vcf('founders.vcf')
#' alleles <- utl::mapping(c(FOUNDER1='A', FOUNDER2='B'))
#' utl::convert_vcf_to_genfile(infiles, genfile, sample.ids, founder.ids,
#'                             alleles=alleles, max.seqlength=max.seqlength)
#'
#' }
#'
#' @seealso \href{https://bioconductor.org/packages/release/bioc/html/Biostrings.html}{Biostrings package}
#' @seealso \href{https://bioconductor.org/packages/release/bioc/html/GenomeInfoDb.html}{GenomeInfoDb package}
#' @seealso \href{https://broadinstitute.github.io/picard/}{Picard Tools documentation}
#' @seealso \href{http://www.rqtl.org}{R/qtl website}
#'
#' @export
#' @rdname convert_vcf_to_genfile
convert_vcf_to_genfile <- function(infiles, genfile, samples, founders, alleles=NULL,
                                   max.seqlength=NULL, na.string=c('-', 'NA')) {
  
  stopifnot( is_single_string(genfile) )
  na.string <- match.arg(na.string)
  
  clashing.ids <- intersect(samples, founders)
  if ( length(clashing.ids) > 0 ) {
    stop("sample/founder ID clash - '", toString(clashing.ids), "'")
  }
  
  sample.data <- read_snps_from_vcf(infiles, samples=samples, max.seqlength=max.seqlength,
                                    require.any=TRUE, require.polymorphic=TRUE)
  
  founder.data <- read_snps_from_vcf(infiles, samples=founders, max.seqlength=max.seqlength,
                                     require.all=TRUE, require.polymorphic=TRUE)
  
  geno.mat <- make_geno_matrix(sample.data, founder.data, alleles=alleles)
  snp.loc <- parse_snp_marker_ids(colnames(geno.mat))
  
  id.col <- c('id', '', rownames(geno.mat))
  geno.mat <- rbind(colnames(geno.mat), snp.loc$chr, geno.mat)
  geno.mat <- cbind(id.col, geno.mat)
  
  utils::write.table(geno.mat, file=genfile, na=na.string, sep=',',
                     quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  return( invisible() )
}

# End of convert_vcf_to_genfile.R

# Start of internal.R

# get_col_indices
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
#' @rdname get_col_indices
get_col_indices <- function(x, requested=NULL, strict=FALSE) {
  
  num.cols <- ncol(x)
  
  if ( ! is_single_nonnegative_whole_number(num.cols) ) {
    stop("cannot get column indices - object does not have columns")
  }
  
  available <- if ( num.cols > 0 ) { 1:num.cols } else { integer() }
  names(available) <- colnames(x)
  
  resolved <- get_indices(available, requested=requested, strict=strict)
  indices <- unname(available[resolved])
  
  return(indices)
}

# get_indices
#' Get indices of object elements.
#'
#' Get the indices of \code{x}, as constrained by the \code{requested} parameter. If using this
#' function without \code{requested} constraints, consider using the faster primitive R function
#' \code{seq_along} instead.
#'
#' @param x An object with elements that are accessible by an index.
#' @param requested A character vector of names, a logical vector of the same length as \code{x}, or
#' a numeric vector containing indices of \code{x}. If this parameter is not specified, all indices
#' are returned.
#' @param strict Option indicating that \code{requested}, if specified, must request indices that
#' are unique and in the same order as the corresponding elements of \code{x}.
#'
#' @return Integer vector containing indices of \code{x}.
#'
#' @keywords internal
#' @rdname get_indices
get_indices <- function(x, requested=NULL, strict=FALSE) {
  
  stopifnot( isTRUE(strict) || identical(FALSE, strict) )
  
  object.length <- length(x)
  
  if ( ! is_single_nonnegative_whole_number(object.length) ) {
    stop("cannot get object indices - object does not have length")
  }
  
  indices <- seq_along(x)
  
  if ( ! is.null(requested) ) {
    
    if ( is.numeric(requested) ) {
      
      nonintegers <- requested[ ! is_whole_number(requested) ]
      if ( length(nonintegers) > 0L ) {
        stop("requested indices are not integers - '", toString(nonintegers), "'")
      }
      
      exrange <- requested[ ! requested %in% indices ]
      if ( length(exrange) > 0L ) {
        stop("requested indices out of range - '", toString(exrange), "'")
      }
      
      indices <- indices[requested]
      
    } else if ( is.logical(requested) ) {
      
      if ( length(requested) != object.length ) {
        stop("cannot resolve indices by logical vector - length mismatch")
      }
      
      indices <- unname(which(requested))
      
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
      if ( length(unfound) > 0L ) {
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

# get_lod_col_index
#' Get LOD column index.
#'
#' @param x An \pkg{R/qtl} \code{scanone}, \code{scanoneperm}, or \code{summary.scanoneperm} object.
#' @param lodcolumn This parameter indicates for which LOD column an index should be returned. This
#' must be either a LOD column name or an index \emph{with respect to the set of LOD columns}. If no
#' LOD column is specified and one such column is found, the index of that column is returned by
#' default; otherwise a LOD column must be specified.
#'
#' @return LOD column index.
#'
#' @keywords internal
#' @rdname get_lod_col_index
get_lod_col_index <- function(x, lodcolumn=NULL) {
  
  lodcol.indices <- get_lod_col_indices(x, lodcolumns=lodcolumn)
  
  if ( length(lodcol.indices) > 1L ) {
    stop("object has multiple LOD columns - please choose one")
  } else if ( length(lodcol.indices) == 0L ) {
    stop("no LOD column found")
  }
  
  return(lodcol.indices[1L])
}

# get_lod_col_indices
#' Get LOD column indices.
#'
#' @param x An \pkg{R/qtl} \code{scanone}, \code{scanoneperm}, or \code{summary.scanoneperm} object.
#' @param lodcolumns This parameter indicates for which LOD columns indices should be returned. The
#' specified LOD columns must be a character vector of LOD column names, a logical vector with
#' length equal to the number of LOD columns, or a numeric vector of indices \emph{with respect
#' to the set of LOD columns}. If no LOD columns are specified, all LOD column indices are returned.
#' @param strict Option indicating that \code{lodcolumns}, if specified, must request LOD column
#' indices that are unique and in the same order as the LOD columns of \code{x}.
#'
#' @return Vector of LOD column indices.
#'
#' @keywords internal
#' @rdname get_lod_col_indices
get_lod_col_indices <- function(x, lodcolumns=NULL, strict=FALSE) {
  UseMethod('get_lod_col_indices', x)
}

# get_lod_col_indices.scanone
#' @method get_lod_col_indices scanone
#' @rdname get_lod_col_indices
get_lod_col_indices.scanone <- function(x, lodcolumns=NULL, strict=FALSE) {
  stopifnot( identical(colnames(x)[1:2], c('chr', 'pos')) )
  available <- seq_len(ncol(x))[-(1:2)]
  names(available) <- colnames(x)[available]
  resolved <- get_indices(available, requested=lodcolumns, strict=strict)
  indices <- unname(available[resolved])
  return(indices)
}

# get_lod_col_indices.scanoneperm
#' @method get_lod_col_indices scanoneperm
#' @rdname get_lod_col_indices
get_lod_col_indices.scanoneperm <- function(x, lodcolumns=NULL, strict=FALSE) {
  stopifnot( ncol(x) >= 1L )
  available <- seq_len(ncol(x))
  names(available) <- colnames(x)
  resolved <- get_indices(available, requested=lodcolumns, strict=strict)
  indices <- unname(available[resolved])
  return(indices)
}

# get_lod_col_indices.summary.scanoneperm
#' @export
#' @method get_lod_col_indices summary.scanoneperm
#' @rdname get_lod_col_indices
get_lod_col_indices.summary.scanoneperm <- function(x, lodcolumns=NULL, strict=FALSE) {
  stopifnot( ncol(x) >= 1L )
  available <- seq_len(ncol(x))
  names(available) <- colnames(x)
  resolved <- get_indices(available, requested=lodcolumns, strict=strict)
  indices <- unname(available[resolved])
  return(indices)
}

# get_run_index_list
#' Get index list of successive runs in a vector.
#'
#' @param x A vector.
#'
#' @return List of integer vectors, each containing indices for a run of repeated values in
#' \code{x}. Each list element takes its name from the corresponding repeated value. Returns
#' an empty list if \code{x} is of length zero.
#'
#' @keywords internal
#' @rdname get_run_index_list
get_run_index_list <- function(x, na.rm=FALSE) {
  
  stopifnot( is.vector(x) )
  stopifnot( isTRUE(na.rm) || identical(FALSE, na.rm) )
  
  if ( length(x) > 0L ) {
    
    # Get run-length encoding of vector.
    runs <- rle(x)
    
    # Set run names from RLE values.
    run.names <- runs$values
    
    # Get number of runs in RLE.
    num.runs <- unique(lengths(runs))
    
    # Get last index of each run.
    J <- cumsum(runs$lengths)
    
    # Get first index of each run.
    if ( num.runs > 1L ) {
      I <- c(1L, sapply(J[1:(length(J)-1)], function(j) j + 1))
    } else {
      I <- 1L
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

# infer_map_precision
#' Infer precision of map positions.
#'
#' @param x An \pkg{R/qtl} \code{scanone} or \code{map} object.
#' @param tol Tolerance for map position equality.
#'
#' @return Inferred precision of map in \code{x}. Map precision is taken from the map step size, if
#' the majority of the map steps in \code{x} are the same. Otherwise, map precision is set from the
#' smallest distance between any two adjacent loci.
#'
#' @keywords internal
#' @rdname infer_map_precision
infer_map_precision <- function(x, tol=.Machine$double.eps^0.5) {
  
  stopifnot( any( c('map', 'scanone') %in% class(x) ) )
  stopifnot( is_single_nonnegative_number(tol) )
  
  if ( 'scanone' %in% class(x) ) {
    map.steps <- unlist(lapply(unique(x$chr), function(map.seq)
      diff(x[x$chr == map.seq, 'pos'])))
  } else {  # 'map' %in% class(x)
    map.steps <- unlist(unname(lapply(x, function(map.pos) diff(as.numeric(map.pos)))))
  }
  
  # Get frequency table of map steps.
  step.freqs <- table(map.steps)
  
  # Get numeric step sizes.
  step.sizes <- as.numeric(names(step.freqs))
  
  # Get differences between step sizes.
  step.diffs <- diff(step.sizes)
  
  # Group step-size values that are very similar.
  size.groups <- vector('list')
  i <- 1L
  while ( i <= length(step.freqs) ) {
    
    j <- i
    
    while ( j < length(step.freqs) && step.diffs[j] < tol ) {
      j <- j + 1L
    }
    
    size.groups <- append(size.groups, list(i:j))
    
    i <- j + 1L
  }
  
  # Merge similar step values.
  merged.freqs <- integer(length=length(size.groups))
  for ( i in seq_along(size.groups) ) {
    g <- unlist(size.groups[i])
    merged.freqs[i] <- sum(step.freqs[g])
    names(merged.freqs)[i] <- names(sort(step.freqs[g], decreasing=TRUE))[1L]
  }
  
  # Sort steps by decreasing frequency.
  sorted.freqs <- sort(merged.freqs, decreasing=TRUE)
  
  # If the most frequent step is more frequent than all others combined, set as map precision..
  if ( length(sorted.freqs) == 1L || (length(sorted.freqs) > 1L &&
                                      sorted.freqs[1L] >= sum(sorted.freqs[2:length(sorted.freqs)])) ) {
    map.precision <- as.numeric(names(sorted.freqs)[1L])
    # ..otherwise set smallest step as map precision.
  } else {
    map.precision <- min(map.steps, na.rm=TRUE)
  }
  
  return(map.precision)
}

# is_single_char
#' Test for a single character.
#'
#' @param x Test object.
#'
#' @return \code{TRUE} if \code{x} is a single character; \code{FALSE} otherwise.
#'
#' @keywords internal
#' @rdname is_single_char
is_single_char <- function(x) {
  return( is.character(x) && length(x) == 1L && ! is.na(x) && nchar(x) == 1L )
}

# is_single_nonnegative_number
#' Test for a single non-negative number.
#'
#' @param n Test object.
#'
#' @return \code{TRUE} if \code{x} is a single non-negative number; \code{FALSE} otherwise.
#'
#' @keywords internal
#' @rdname is_single_nonnegative_number
is_single_nonnegative_number <- function(n) {
  return( is.numeric(n) && length(n) == 1L && is.finite(n) && n >= 0.0 )
}

# is_single_nonnegative_whole_number
#' Test for a single non-negative whole number.
#'
#' @param n Test object.
#' @param tol Numeric tolerance.
#'
#' @return \code{TRUE} if \code{x} is a single non-negative whole number; \code{FALSE} otherwise.
#'
#' @keywords internal
#' @rdname is_single_nonnegative_whole_number
is_single_nonnegative_whole_number <- function(n, tol=.Machine$double.eps^0.5) {
  return( is.numeric(n) && length(n) == 1L && is.finite(n) && n >= 0.0 &&
            abs(n - round(n)) < abs(tol) )
}

# is_single_positive_whole_number
#' Test for a single positive whole number.
#'
#' @param n Test object.
#' @param tol Numeric tolerance.
#'
#' @return \code{TRUE} if the object is a single positive whole number;
#' \code{FALSE} otherwise.
#'
#' @keywords internal
#' @rdname is_single_positive_whole_number
is_single_positive_whole_number <- function(n, tol=.Machine$double.eps^0.5) {
  return( length(n) == 1 && is.numeric(n) && is.finite(n) &&
            n > tol && abs(n - round(n)) < abs(tol) )
}

# is_single_probability
#' Test for a single valid probability.
#'
#' @param n Test object.
#'
#' @return \code{TRUE} if \code{x} is a single valid probability; \code{FALSE} otherwise.
#'
#' @keywords internal
#' @rdname is_single_probability
is_single_probability <- function(n) {
  return( is.numeric(n) && length(n) == 1L && is.finite(n) && n >= 0.0 && n <= 1.0 )
}

# is_single_string
#' Test for a single character string.
#'
#' @param x Test object.
#'
#' @return \code{TRUE} if \code{x} is a single string; \code{FALSE} otherwise.
#'
#' @keywords internal
#' @rdname is_single_string
is_single_string <- function(x) {
  return( is.character(x) && length(x) == 1L && ! is.na(x) )
}

# is_valid_id
#' Test identifier validity.
#'
#' An identifier is considered valid if it is a string
#' containing at least one non-whitespace character.
#'
#' @param x Test vector.
#'
#' @return Logical vector indicating which elements are valid identifiers.
#'
#' @keywords internal
#' @rdname is_valid_id
is_valid_id <- function(x) {
  return( is.character(x) & nzchar(gsub('[[:space:]]', '', x)) )
}

# is_whole_number
#' Test for whole numbers.
#'
#' @param n Test vector.
#' @param tol Numeric tolerance.
#'
#' @return Logical vector indicating which elements of \code{n} are whole numbers.
#'
#' @keywords internal
#' @rdname is_whole_number
is_whole_number <- function(n, tol=.Machine$double.eps^0.5) {
  return( is.numeric(n) & is.finite(n) & abs(n - round(n)) < abs(tol) )
}

# isBOOL
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

# isFALSE
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

# make_geno_matrix
#' Make genotype matrix from sample and founder genotype data.
#'
#' Given input genotype data for cross samples and their founder strains, this
#' function assigns a genotype symbol to each locus according to the inferred
#' combination of founder alleles.
#'
#' @param x Raw sample genotype \code{matrix}.
#' @param y Raw founder genotype \code{matrix}.
#' @param alleles Mapping of founder IDs to founder allele symbols.
#'
#' @return A founder genotype matrix, with genotypes encoded as integers and
#' their corresponding genotype symbols in the attribute \code{'genotypes'}.
#'
#' @keywords internal
#' @rdname make_geno_matrix
make_geno_matrix <- function(x, y, alleles=NULL) {
  
  founder.allele.charset <- c(LETTERS, letters)
  
  x.samples <- rownames(x)
  x.snps <- colnames(x)
  
  y.samples <- rownames(y)
  y.snps <- colnames(y)
  
  # Check upper limit of founder sample count.
  # NB: absolute upper limit is length of founder.allele.charset
  if ( length(y.samples) != 2 ) {
    stop("unsupported number of founders - '", length(y.samples), "'")
  }
  
  if ( ! is.null(alleles) ) {
    
    stopifnot( is_mapping(alleles) )
    
    founder.samples <- mapping_keys(alleles)  # NB: ensures unique sample IDs
    founder.symbols <- as.character( mapping_values(alleles) )
    
    unknown <- founder.samples[ ! founder.samples %in% y.samples ]
    if ( length(unknown) > 0 ) {
      stop("allele symbols specified for unknown founder samples - '",
           toString(unknown), "'")
    }
    
    unfound <- y.samples[ ! y.samples %in% founder.samples ]
    if ( length(unfound) > 0 ) {
      stop("allele symbols not specified for founder samples - '",
           toString(unfound), "'")
    }
    
    dup.alleles <- founder.symbols[ duplicated(founder.symbols) ]
    if ( length(dup.alleles) > 0 ) {
      stop("duplicate founder alleles - '", toString(dup.alleles), "'")
    }
    
    err.alleles <- founder.symbols[ ! founder.symbols %in% founder.allele.charset ]
    if ( length(err.alleles) > 0 ) {
      stop("invalid founder alleles - '", toString(err.alleles), "'")
    }
    
    names(founder.symbols) <- founder.samples
    allele.symbols <- sort( founder.symbols )
    
  } else {
    
    allele.symbols <- founder.allele.charset[ seq_along(y.samples) ]
    names(allele.symbols) <- y.samples
  }
  
  # Keep only loci common to samples and founders.
  common.snps <- intersect(x.snps, y.snps)
  x <- x[, common.snps, drop=FALSE]
  y <- y[, common.snps, drop=FALSE]
  
  # Init genotype matrix.
  geno.mat <- matrix(NA_character_, nrow=nrow(x), ncol=ncol(x),
                     dimnames=dimnames(x))
  
  for ( j in get_col_indices(geno.mat) ) {
    
    # Get sample symbols and genotypes for this
    # locus, skip if effectively monomorphic.
    x.symbols <- x[, j]
    uniq.x.symbols <- unique(x.symbols)
    uniq.x.geno <- uniq.x.symbols[ ! is.na(uniq.x.symbols) ]
    if ( length(uniq.x.geno) < 2 ) {
      next
    }
    
    # Get founder genotypes for this locus,
    # skip if any are missing.
    y.genotypes <- y[, j]
    if ( anyNA(y.genotypes) ) {
      next
    }
    
    # Get founder alleles, skip if monomorphic or heterozygous.
    y.alleles <- lapply(y.genotypes, function(gt)
      unique(strsplit(gt, '')[[1]]))
    if ( anyDuplicated(y.alleles) || any( lengths(y.alleles) > 1 ) ) {
      next
    }
    y.alleles <- unlist(y.alleles)
    
    # Assign locus genotypes from matching founder.
    x.calls <- lapply(x.symbols, function(gt) {
      x.alleles <- strsplit(gt, '')[[1]]
      x.indices <- match(x.alleles, y.alleles)
      if ( anyNA(x.indices) ) {
        x.genotype <- NA_character_
      } else {
        x.fdr.names <- names(y.alleles)[x.indices]
        # Sort genotype alleles, since 'BA'
        # is considered equivalent to 'AB'.
        x.fdr.geno <- sort(allele.symbols[x.fdr.names])
        x.genotype <- paste0(x.fdr.geno, collapse='')
      }
    })
    
    # Skip if monomorphic after matching to founder alleles.
    uniq.x.calls <- unique(x.calls[ ! is.na(x.calls) ])
    if ( length(uniq.x.calls) < 2 ) {
      next
    }
    
    geno.mat[, j] <- unlist(x.calls)
  }
  
  # Remove null loci.
  geno.mat <- remove_na_cols(geno.mat)
  
  if ( ncol(geno.mat) == 0 ) {
    stop("cannot make genotype matrix - no suitable loci found")
  }
  
  return(geno.mat)
}

# make_snp_marker_ids
#' Make SNP marker IDs for loci.
#'
#' @description This function creates SNP marker
#' IDs from locus \code{data.frame} \code{'loc'}.
#'
#' @param loc Locus \code{data.frame}, with columns \code{'chr'} and
#' \code{'pos'}, specifying physical map positions in base-pair units.
#' @param max.seqlength Maximum expected sequence length. This determines
#' the width to which SNP positions are zero-padded. By default, this is
#' the maximum position in the input \code{data.frame}.
#'
#' @return Character vector of SNP marker IDs.
#'
#' @author Thomas A. Walsh
#' @author Yue Hu
#'
#' @keywords internal
#' @rdname make_snp_marker_ids
make_snp_marker_ids <- function(loc, max.seqlength=NULL) {
  
  stopifnot( is.data.frame(loc) )
  stopifnot( nrow(loc) > 0 )
  
  loc.pos <- loc$pos
  
  max.pos <- max(loc.pos, na.rm=TRUE)
  stopifnot( is_single_positive_whole_number(max.pos) )
  max.pos.width <- floor(log10(max.pos)) + 1
  
  pad.width <- NULL
  if ( ! is.null(max.seqlength) ) {
    stopifnot( is_single_positive_whole_number(max.seqlength) )
    if ( max.pos <= max.seqlength ) {
      pad.width <- floor(log10(max.seqlength)) + 1
    } else {
      warning("input position (", format(max.pos, scientific=FALSE), "",
              ") exceeds specified maximum sequence length (",
              format(max.seqlength, scientific=FALSE), "), padding to ",
              max.pos.width, " digits")
    }
  }
  
  if ( is.null(pad.width) ) {
    pad.width <- max.pos.width
  }
  
  pad.fmt.str <- paste0('%0', pad.width, 'd')
  loc.pos <- sprintf(pad.fmt.str, loc.pos)
  
  return( paste0(loc$chr, ':', loc.pos) )
}

# parse_snp_marker_ids
#' Parse SNP marker IDs.
#'
#' This function parses an input vector of SNP marker IDs, and returns a \code{data.frame} with
#' the locus in each row derived from the corresponding marker ID in the input vector. An error
#' is raised if any of the input values cannot be parsed as a SNP marker ID.
#'
#' @param ids Vector of SNP marker IDs.
#'
#' @return A \code{data.frame} with loci
#' corresponding to the SNP marker IDs.
#'
#' @keywords internal
#' @rdname parse_snp_marker_ids
parse_snp_marker_ids <- function(ids) {
  
  snp_marker_id_patt = '^([^:]+):([[:digit:]]+)$'
  m <- regexec(snp_marker_id_patt, ids)
  regmatch.list <- regmatches(ids, m)
  
  valid.snp.marker.ids <- lengths(regmatch.list) > 0
  stopifnot( all(valid.snp.marker.ids) )
  
  snp.seqs <- sapply(regmatch.list, getElement, 2)
  snp.pos <- sapply(regmatch.list, function(x) as.numeric(x[3]))
  
  return( data.frame(chr=snp.seqs, pos=snp.pos, row.names=ids, stringsAsFactors=FALSE) )
}

# read_snps_from_vcf
#' Read SNP genotypes from VCF files.
#'
#' This function reads SNP genotype data from one or more VCF files, and returns
#' these as a \code{matrix} of SNP genotypes - effectively variant base calls.
#'
#' @param infiles Input VCF file paths.
#' @param samples Optional vector of samples for which SNP genotypes should be obtained.
#' If not specified, genotypes are returned for all available samples.
#' @param max.seqlength Optional parameter to indicate the maximum reference sequence length, which
#' is used to determine the zero-padded width of genomic positions in SNP marker IDs. Without this
#' information, SNP marker IDs may be formatted inconsistently in different datasets.
#' @param require.all Remove variants that are not completely genotyped
#' with respect to the given samples.
#' @param require.any Remove variants that do not have at least one genotype
#' call among the given samples.
#' @param require.polymorphic Remove variants that do not have at least two
#' different genotype calls among the given samples.
#'
#' @return A \code{matrix} object containing SNP genotype data.
#'
#' @keywords internal
#' @rdname read_snps_from_vcf
read_snps_from_vcf <- function(infiles, samples=NULL, max.seqlength=NULL,
                               require.all=FALSE, require.any=FALSE,
                               require.polymorphic=FALSE) {
  
  stopifnot( length(infiles) > 0 )
  stopifnot( isBOOL(require.all) )
  stopifnot( isBOOL(require.any) )
  stopifnot( isBOOL(require.polymorphic) )
  
  # Get samples in each input VCF file.
  sample.list <- lapply(infiles, read_samples_from_vcf)
  
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
  
  invalid.ids <- samples[ ! is_valid_id(samples) ]
  if ( length(invalid.ids) > 0 ) {
    stop("invalid sample IDs - '", toString(invalid.ids), "'")
  }
  
  # Keep only files relevant for the given samples.
  relevant <- lengths(sample.list) > 0
  sample.list <- sample.list[relevant]
  infiles <- infiles[relevant]
  
  # Init list of SNP data.
  snps <- vector('list', length(infiles))
  
  # Init list for mapping genotype strings to allele indices.
  geno.memo <- list()
  
  # Read SNP genotypes from each input VCF file.
  for ( i in seq_along(infiles) ) {
    
    file.samples <- sample.list[[i]]
    infile <- infiles[[i]]
    
    # Set VCF parameters.
    param <- VariantAnnotation::ScanVcfParam(fixed=c('ALT', 'QUAL'),
                                             geno='GT', samples=file.samples)
    
    # Read VCF.
    vcf <- VariantAnnotation::readVcf(infile, param=param)
    stopifnot( length(vcf) > 0 )
    
    # Get variant data from VCF.
    var.seqs <- as.character( GenomeInfoDb::seqnames(vcf) ) # reference sequence
    var.pos <- BiocGenerics::start( IRanges::ranges(vcf) )  # position
    var.ref <- as.character( VariantAnnotation::ref(vcf) )  # reference allele
    var.alt <- lapply( IRanges::CharacterList(              # alternative alleles
      VariantAnnotation::alt(vcf) ), as.character)
    var.geno <- VariantAnnotation::geno(vcf)                # genotype info
    geno.mat <- t( var.geno$GT )                            # genotype calls
    
    # Set indices of SNP variants.
    snp.indices <- which( var.ref %in% Biostrings::DNA_BASES &
                            sapply(var.alt, function(x) all(x %in% Biostrings::DNA_BASES)) )
    num.snps <- length(snp.indices)
    stopifnot( num.snps > 0 )
    
    # Filter genotype data to retain only SNP variants.
    geno.mat <- geno.mat[, snp.indices, drop=FALSE]
    
    # Combine REF and ALT alleles for each SNP.
    var.alleles <- lapply(snp.indices, function(j)
      c(var.ref[j], var.alt[[j]]))
    
    # Create dataframe with variant locus info.
    snp.loc <- data.frame(chr=var.seqs[snp.indices],
                          pos=var.pos[snp.indices])
    
    # If maximum reference sequence length not
    # specified, try to get that info from VCF.
    if ( is.null(max.seqlength) ) {
      seqinfo <- VariantAnnotation::seqinfo(vcf)
      obs.seqnames <- unique(var.seqs)
      known.seqnames <- GenomeInfoDb::seqnames(seqinfo)
      unknown.seqnames <- setdiff(obs.seqnames, known.seqnames)
      if ( ! anyNA(obs.seqnames) && length(unknown.seqnames) == 0 ) {
        known.seqlengths <- stats::setNames(GenomeInfoDb::seqlengths(seqinfo),
                                            known.seqnames)
        obs.seqlengths <- known.seqlengths[obs.seqnames]
        max.obs.seqlength <- max(obs.seqlengths, na.rm=TRUE)
        if ( is_single_positive_whole_number(max.obs.seqlength) ) {
          max.seqlength <- max.obs.seqlength
        }
      }
    }
    if ( is.null(max.seqlength) ) {
      warning("maximum reference sequence length is unknown - ",
              "marker IDs may be inconsistent between datasets")
    }
    
    # Get marker IDs for SNPs in this file.
    file.snps <- make_snp_marker_ids(snp.loc, max.seqlength=max.seqlength)
    
    # Check for multiple variants coinciding at same locus.
    # TODO: handle coinciding variants
    if ( anyDuplicated(file.snps) ) {
      stop("coinciding variants in file - '", infile, "'")
    }
    
    colnames(geno.mat) <- file.snps
    
    # Resolve raw genotypes for each variant as concatenated base calls.
    for ( j in 1:num.snps ) {
      
      # Get vector of alleles for this variant record.
      rec.alleles <- var.alleles[[j]]
      
      # Get genotype strings for each sample in this variant record.
      geno.data <- unname( geno.mat[, j] )
      
      # Resolve previously unseen genotype strings
      # to their corresponding allele indices.
      for ( unique.call in unique(geno.data) ) {
        if ( ! unique.call %in% names(geno.memo) ) {
          # Convert allele indices from 0-offset strings to
          # 1-offset integers, VCF missing '.' to NA value.
          indicesC <- unlist( strsplit(unique.call, '[/|]') )
          indices0 <- suppressWarnings( as.integer(indicesC) )
          indices1 <- indices0 + 1
          geno.memo[[unique.call]] <- indices1
        }
      }
      
      # Resolve genotype calls for this variant record.
      geno.calls <- lapply(geno.data, function(geno.str) {
        alel.idxs <- geno.memo[[geno.str]]
        alel.calls <- rec.alleles[alel.idxs]
        if ( anyNA(alel.calls) ) {
          geno.call <- NA_character_
        } else {
          geno.call <- paste0(alel.calls, collapse='')
        }
      })
      
      geno.mat[, j] <- unlist(geno.calls)
    }
    
    # Set file variant data, sorting columns
    # to facilitate subsequent merging.
    snps[[i]] <- geno.mat[, order(colnames(geno.mat)), drop=FALSE]
  }
  
  # If multiple relevant files, combine SNP genotypes,
  # taking the union of all variants in input files..
  if ( length(infiles) > 1 ) {
    
    # Prep combined variant data.
    combined.samples <- sort( unique( unlist( lapply(snps, rownames) ) ) )
    combined.snps <- sort( unique( unlist( lapply(snps, colnames) ) ) )
    combined.names <- list(combined.samples, combined.snps)
    
    # Set combined variant data from each input file.
    # NB: we previously checked for duplicate samples and coinciding
    # variants, so we can assume that there will be no conflicts.
    result <- matrix(NA_character_, nrow=length(combined.samples),
                     ncol=length(combined.snps), dimnames=combined.names)
    for ( i in seq_along(snps) ) {
      dest.idxs <- match(colnames(snps[[i]]), colnames(result))
      for ( sample.id in rownames(snps[[i]]) ) {
        result[sample.id, dest.idxs] <- snps[[i]][sample.id, ]
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
        ! anyNA(result[, i])
      )
    }
    
    if (require.any) { # Remove variants without genotypes.
      mask <- mask & sapply(1:num.snps, function(i)
        ! all( is.na(result[, i]) )
      )
    }
    
    if (require.polymorphic) { # Remove monomorphic variants.
      mask <- mask & sapply(1:num.snps, function(i) {
        geno.data <- result[, i]
        g.symbols <- unique(geno.data)
        genotypes <- g.symbols[ ! is.na(g.symbols) ]
        return( length(genotypes) > 1 )
      })
    }
    
    # Apply filter mask.
    result <- result[, mask, drop=FALSE]
  }
  
  return(result)
}

# remove_na_cols
#' Remove columns containing only \code{NA} values.
#'
#' @param x A \code{data.frame} or \code{matrix}.
#'
#' @return Input object without those columns that contain only \code{NA} values.
#'
#' @keywords internal
#' @rdname remove_na_cols
remove_na_cols <- function(x) {
  stopifnot( is.data.frame(x) || is.matrix(x) )
  stopifnot( nrow(x) > 0 )
  mask <- ! apply( x, 2, function(column) all( is.na(as.vector(column)) ) )
  x <- x[, mask, drop=FALSE]
  return(x)
}

# End of internal.R




# Start of read_samples_from_vcf.R

# read_samples_from_vcf
#' Read sample IDs from a VCF file.
#'
#' @param file A VCF file path.
#'
#' @return Vector of the sample IDs in the given VCF file.
#'
#' @export
#' @rdname read_samples_from_vcf
read_samples_from_vcf <- function(file) {
  stopifnot( is_single_string(file) )
  stopifnot( file.exists(file) )
  header <- VariantAnnotation::scanVcfHeader(file)
  return( VariantAnnotation::samples(header) )
}

# End of read_samples_from_vcf.R


