db_download_ncbi <- function(verbose = TRUE) {
  # set paths
  db_url <- 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip'
  db_path_file <- file.path(tdb_cache$cache_path_get(), 'taxdump.zip')
  db_path_dir <- file.path(tdb_cache$cache_path_get(), 'taxdump')
  ncbi_names_file <- file.path(db_path_dir, 'names.dmp')
  ncbi_nodes_file <- file.path(db_path_dir, 'nodes.dmp')
  final_file <- file.path(tdb_cache$cache_path_get(), 'NCBI.sql')
  
  if(file.exists(final_file)){
    mssg(verbose, "Database already exists, returning old file")
    return(final_file)
  }
  
  # make home dir if not already present
  tdb_cache$mkdir()
  # download data
  mssg(verbose, 'downloading...')
  curl::curl_download(db_url, db_path_file, quiet = TRUE)
  # unzip
  mssg(verbose, 'unzipping...')
  utils::unzip(db_path_file, files = c('names.dmp', 'nodes.dmp'), exdir = db_path_dir)
  
  # Taxonomy names file (names.dmp):
  #   tax_id       -- the id of node associated with this name
  #   name_txt     -- name itself
  #   unique name  -- the unique variant of this name if name not unique
  #   name class   -- (synonym, common name, ...)
  mssg(verbose, "loading 'names.dmp'...")
  ncbi_names <- readr::read_tsv(
    ncbi_names_file,
    col_names = c("tax_id", "name_txt", "unique_name", "name_class"),
    col_type = "i_c_c_c_",
    quote = ""
  )
  
  # nodes.dmp file consists of taxonomy nodes. The description for each node includes the following
  #   tax_id                                 -- node id in GenBank taxonomy database
  #   parent_tax_id                          -- parent node id in GenBank taxonomy database
  #   rank                                   -- rank of this node (superkingdom, kingdom, ...)
  #   embl_code                              -- locus-name prefix; not unique
  #   division_id                            -- see division.dmp file
  #   inherited_div_flag            (1 or 0) -- 1 if node inherits division from parent
  #   genetic_code_id                        -- see gencode.dmp file
  #   inherited_GC_flag             (1 or 0) -- 1 if node inherits genetic code from parent
  #   mitochondrial_genetic_code_id          -- see gencode.dmp file
  #   inherited_MGC_flag            (1 or 0) -- 1 if node inherits mitochondrial gencode from parent
  #   GenBank_hidden_flag           (1 or 0) -- 1 if name is suppressed in GenBank entry lineage
  #   hidden_subtree_root_flag      (1 or 0) -- 1 if this subtree has no sequence data yet
  #   comments                               -- free-text comments and citations
  mssg(verbose, "loading 'nodes.dmp'...")
  ncbi_nodes <- readr::read_tsv(
    ncbi_nodes_file,
    col_names=c(
      "tax_id",
      "parent_tax_id",
      "rank",
      "embl_code",
      "division_id",
      "inherited_div_flag",
      "genetic_code_id",
      "inherited_GC_flag",
      "mitochondrial_genetic_code_id",
      "inherited_MGC_flag",
      "GenBank_hidden_flag",
      "hidden_subtree_root_flag",
      "comments"
    ),
    col_types='i_i_c_c_i_i_i_i_i_i_i_i_c_',
    quote=""
  )
  
  mssg(verbose, 'building hierarchy table...')
  # will hold a table for every taxonomic level ascending from leaf to root the
  # length of this list will equal the depth of the taxonomic tree (e.g. 37
  # currently)
  hierarchs <- list()
  # set up the base table with columns 'tax_id', 'ancestor', and 'level', where
  # level is 1 for immediate parent
  hierarchs[[1]] <- ncbi_nodes[, c('tax_id', 'parent_tax_id')] %>%
    magrittr::set_names(c('tax_id', 'ancestor'))
  hierarchs[[1]]$level <- 1
  # make a child to parent map
  child2parent <- ncbi_nodes$parent_tax_id
  names(child2parent) <- ncbi_nodes$tax_id
  # Iteratively replace the ancestor column with the ancestor parent until all
  # lineages converge to root. Each iteration is stored in a new table with
  # level incremented.
  while(TRUE) {
    top <- tail(hierarchs, 1)[[1]]
    incomplete <- top$ancestor != 1L # 1 is the taxonomy root id
    top <- top[incomplete, ]
    if(nrow(top) == 0){
      break
    }
    hierarchs[[length(hierarchs)+1]] <- tibble::tibble(
      tax_id = top$tax_id,
      ancestor = child2parent[as.character(top$ancestor)],
      level = rep(length(hierarchs) + 1, nrow(top))
    )
  }
  # Bind all levels into one table.
  hierarchy <- do.call(rbind, hierarchs)
  hierarchy$level <- as.integer(hierarchy$level)
  
  
  mssg(verbose, 'building SQLite database...')
  
  db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname=final_file)
  
  # Create tables - I have to manually make the `names` table because I have to
  # set `COLLATE NOCASE` at the time of table creation.
  RSQLite::dbExecute(conn=db, "
    CREATE TABLE names (
      tax_id INTEGER,
      name_txt TEXT COLLATE NOCASE,
      unique_name TEXT,
      name_class TEXT
    )
    "
  )
  
  # Load tables
  RSQLite::dbWriteTable(
    conn   = db,
    name   = 'names',
    value  = as.data.frame(ncbi_names),
    append = TRUE # since I explicitly created the table above
  )
  RSQLite::dbWriteTable(
    conn  = db,
    name  = 'nodes',
    value = as.data.frame(ncbi_nodes),
  )
  RSQLite::dbWriteTable(
    conn  = db,
    name  = 'hierarchy',
    value = as.data.frame(hierarchy),
  )
  
  # Create indices on tax_id columns
  RSQLite::dbExecute(db,
                     'CREATE INDEX tax_id_index_names ON names (tax_id)'
  )
  RSQLite::dbExecute(db,
                     'CREATE INDEX name_txt_index_names ON names (name_txt COLLATE NOCASE)'
  )
  RSQLite::dbExecute(db,
                     'CREATE INDEX tax_id_index_nodes ON nodes (tax_id)'
  )
  RSQLite::dbExecute(db,
                     'CREATE INDEX tax_id_index_hierarchy ON hierarchy (tax_id)'
  )
  RSQLite::dbExecute(db,
                     'CREATE INDEX tax_id_ancestor_hierarchy ON hierarchy (ancestor)'
  )
  
  RSQLite::dbDisconnect(db)
  
  # cleanup
  mssg(verbose, 'cleaning up...')
  unlink(db_path_file)
  unlink(db_path_dir, recursive = TRUE)
  
  return(final_file)
}


mssg <- function(v, ...) if (v) message(...)

txdbc <- function(x) Filter(Negate(is.null), x)

cl <- function(x, y){
  if (is.null(y)) {
    ""
  } else {
    paste0(x, y)
  }
}

mkhome <- function(x) {
  dir.create(x, showWarnings = FALSE, recursive = FALSE)
}

# db_installed(x = 'psql')
db_installed <- function(x) {
  tmp <- Sys.which(x)
  if (any(tmp == "")) {
    nf <- paste0(names(tmp)[tmp == ""], collapse = ", ")
    stop(
      sprintf(
        "\n%s not found on your computer\nInstall the missing tool(s) and try again",
        nf))
  }
}

#' Retrieve the taxonomic hierarchies from local database
#'
#' This function is equivalent to the `taxize::classification()` function,
#' except that it uses a local database (so is much faster) and is currently
#' restricted to handling NCBI taxa. The output is identical to
#' `taxize::classification()`
#'
#' @export
#' @param x character) Vector of taxon keys for the given database
#' @param db character) The database to search
#' @param verbose (logical) Print verbose messages
#' @param ... Additional arguments passed to database specific classification
#' functions.
#' @return list of data.frames with the columns: name, rank, and id. This is
#' exactly equivalent to the output of `taxize::classification()`
#' @examples
#' \dontrun{
#' classification(c(3702, 9606))
#' }
classification.local <- function(x, db='ncbi', verbose=TRUE, ...){
  ap_dispatch(
    x       = x,
    db      = db,
    cmd     = 'classification',
    verbose = verbose,
    ...
  )
}

itis_classification <- function(src, x, ...){
  stop("The ITIS database is currently not supported")
}

tpl_classification <- function(src, x, ...){
  stop("The TPL database is currently not supported")
}

col_classification <- function(src, x, ...){
  stop("The COL database is currently not supported")
}

gbif_classification <- function(src, x, ...){
  stop("The GBIF database is currently not supported")
}

ncbi_classification <- function(src, x, ...){
  
  FUN <- function(src, x, ...){
    # Retrieve the hierarchy for each input taxon id
    cmd <- "SELECT tax_id, level, ancestor FROM hierarchy WHERE tax_id IN (%s)"
    query <- sprintf(cmd, sql_integer_list(x))
    tbl <- sql_collect(src, query)
    # If no IDs were found, return list of NA
    if(nrow(tbl) == 0){
      lineages <- as.list(as.logical(rep(NA, length(x))))
      names(lineages) <- x
      return(lineages)
    }
    
    # Add the query to the lineage as the lowest level
    rbind(tbl, tibble::tibble(
      tax_id   = x,
      ancestor = x,
      level    = rep(0L, length(x))
    )) %>%
      # NOTE: Remove the root node, for consistency with 'taxize'. The root
      # node really is important, though, because viruses are a thing.
      dplyr::filter(.data$ancestor != 1L) %>%
      # Add ranks (TODO: add taxid2rank function)
      merge({
        cmd <- "SELECT tax_id, rank FROM nodes WHERE tax_id IN (%s)"
        query <- sprintf(cmd, sql_integer_list(.$ancestor))
        sql_collect(src, query)
      }, by.x='ancestor', by.y='tax_id') %>%
      dplyr::mutate(
        # make taxon IDs character vectors (for consistency with taxize)
        ancestor = as.character(.data$ancestor),
        # add ancestor scientific name
        name = taxid2name(.data$ancestor)
      ) %>%
      split(f=.$tax_id) %>%
      lapply(function(d)
        dplyr::arrange(d, -.data$level) %>%
          # NOTE: Here I drop the 'level' column. I do this because it is not present
          # in the taxize::classification output. However, without the level column,
          # the ancestor order is encoded only in the row order of the data.frame,
          # which is not robost.
          dplyr::select(
            name = .data$name,
            rank = .data$rank,
            id   = .data$ancestor
          )
      )
  }
  
  ## TODO: probably the Right missing value is this:
  # missing = data.frame(
  #   name = character(),
  #   rank = character(),
  #   id   = character(),
  #   stringsAsFactors=FALSE
  # ),
  missing=NA
  
  ncbi_apply(src, x, FUN, missing=missing, ...)
}

### Internal utilities

# convert a vector to a comma separated string list suitable for SQL, e.g.
# c("Arabidopsis", "Peridermium sp. 'Ysgr-4'") -> "'Arabidopsis',
# 'Peridermium sp. ''Ysgr-4'''"
# Note that double quoting is the SQL convention for escaping quotes in strings
sql_character_list <- function(x){
  if(any(is.na(x))){
    stop("Cannot pass NA into SQL query")
  }
  as.character(x) %>%
    gsub(pattern="'", replacement="''") %>%
    sub(pattern="(.*)", replacement="'\\1'", perl=TRUE) %>%
    paste(collapse=", ")
}

sql_integer_list <- function(x){
  if(any(is.na(x))){
    stop("Cannot pass NA into SQL query")
  }
  x <- as.character(x)
  if(!all(grepl('^[0-9]+$', x, perl=TRUE))){
    stop("Found non-integer where integer required in SQL input")
  }
  paste(x, collapse=", ")
}

ap_vector_dispatch <- function(x, db, cmd, verbose=TRUE, empty=character(0), ...){
  if(is.null(x) || length(x) == 0){
    empty 
  } else {
    FUN <- paste0(db, "_", cmd)
    run_with_db(FUN=get(FUN), db=db, x=x, empty=empty, ...)
  }
}

ap_dispatch <- function(x, db, cmd, out_class=cmd, empty=list(), verbose=TRUE, ...){
  result <- if(is.null(x) || length(x) == 0){
    # For empty or NULL input, return empty list
    empty
  } else {
    FUN <- paste0(db, "_", cmd)
    run_with_db(FUN=get(FUN), db=db, x=x, ...)
  }
  
  attributes(result) <- list(names=names(result), class=out_class, db=db)
  
  if(verbose && all(is.na(result))){
    message("No results found. Consider checking the spelling or alternative classification")
  }
  
  result 
}

run_with_db <- function(FUN, db, ...){
  src <- if(db == 'ncbi'){
    src_ncbi(db_download_ncbi(verbose=FALSE))
  } else {
    stop("Database '", db, "' is not supported")
  }
  FUN(src, ...)
}

ncbi_apply <- function(src, x, FUN, missing=NA, die_if_ambiguous=TRUE, ...){
  # preserve original names (this is important when x is a name vector)
  namemap <- x
  names(namemap) <- x
  # find the entries that are named
  is_named <- !(grepl('^[0-9]+$', x, perl=TRUE) | is.na(x))
  # If x is not integrel, then we assume it is a name.
  if(any(is_named)){
    
    # This is not a pretty solution, since it makes an unnecessary call to the database
    if(die_if_ambiguous){
      # get a table mapping names to taxa
      d <- name2taxid(x[is_named], db='ncbi', out_type="summary")
      # find duplicated elements (ambiguous taxa)
      dups <- unique(d$name_txt[duplicated(d$name_txt)])
      # die if there are any
      if(length(dups) > 0){
        msg <- "The following taxa map to multiple taxonomy ids: "
        msg <- dplyr::group_by(d, .data$name_txt) %>%
          dplyr::filter(length(.data$name_txt) > 1) %>%
          dplyr::summarize(taxids = paste(.data$tax_id, collapse="|")) %>% {
            paste("    ", .$name_txt, " - ", .$taxids, collapse="\n")
          } %>%
          paste(msg, ., sep="\n")
        stop(msg)
      }
    }
    
    # get a table mapping names to taxa
    x[is_named] <- name2taxid(x[is_named], db='ncbi')
    names(namemap)[is_named] <- x[is_named]
  }
  # Remove any taxa that where not found, the missing values will be merged
  # into the final output later (with values of NA)
  x <- x[!is.na(x)]
  
  # If there isn't anything to consider, return a list of NA
  if(length(x) == 0){
    result <- as.list(as.logical(rep(NA, length(namemap))))
    result <- lapply(result, function(x) missing)
    names(result) <- namemap
    return(result)
  }
  
  
  # Run the given function on the clean x
  result <- FUN(src, x, ...)
  
  
  # Map the input names to them
  names(result) <- namemap[names(result)]
  # Add the missing values
  missing_names <- setdiff(namemap, names(result))
  if(length(missing_names) > 0){
    missing_values <- lapply(missing_names, function(x) missing)
    names(missing_values) <- missing_names
    result <- append(result, missing_values)
  }
  # Get result in the input order
  result <- result[as.character(namemap)]
  # Cleanup residual NULLs (if needed)
  result <- lapply(result, function(x) {
    if(is.null(x)){ NA } else { x }
  })
  result[is.na(names(result))] <- missing
  result
}

is_ambiguous <- function(scinames){
  # ambiguous terms (see taxize::ncbi_children.R)
  ambiguous_regex <- paste0(c("unclassified", "environmental", "uncultured", "unknown",
                              "unidentified", "candidate", "sp\\.", "s\\.l\\.", "sensu lato", "clone",
                              "miscellaneous", "candidatus", "affinis", "aff\\.", "incertae sedis",
                              "mixed", "samples", "libaries"), collapse="|")
  grepl(ambiguous_regex, scinames, perl=TRUE)
}


#' @export
#' @rdname src_taxizedb
src_ncbi <- function(path = db_path("ncbi"), ...) {
  stopifnot(file.exists(path))
  con <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = path, ...)
  dbplyr::src_dbi(con, auto_disconnect=TRUE)
}

#' Convert taxon IDs to scientific names
#'
#' @export
#' @param x (character) Vector of taxon keys for the given database
#' @param db (character) The database to search
#' @param verbose (logical) Print verbose messages
#' @param warn (logical) If `TRUE`, raise a warning if any taxon IDs can not
#' be found
#' @param ... Additional arguments passed to database specific classification
#' functions
#' @return character vector of scientific names
#' @examples
#' \dontrun{
#' taxid2name(c(3702, 9606))
#' }
taxid2name <- function(x, db='ncbi', verbose=TRUE, warn=TRUE, ...){
  result <- ap_vector_dispatch(
    x       = x,
    db      = db,
    cmd     = 'taxid2name',
    verbose = verbose,
    warn    = warn,
    empty   = character(0),
    ...
  )
  if(warn && any(is.na(result))){
    msg <- "No name found for %s of %s taxon IDs"
    msg <- sprintf(msg, sum(is.na(result)), length(result))
    if(verbose){
      msg <- paste0(msg, ". The followings are left unnamed: ", 
                    paste0(x[is.na(result)], collapse=', ') 
      )
    }
    warning(msg)
  }
  result 
}

itis_taxid2name <- function(src, x, ...){
  stop("The ITIS database is currently not supported")
}

tpl_taxid2name <- function(src, x, ...){
  stop("The TPL database is currently not supported")
}

col_taxid2name <- function(src, x, ...){
  stop("The COL database is currently not supported")
}

gbif_taxid2name <- function(src, x, ...){
  stop("The GBIF database is currently not supported")
}

ncbi_taxid2name <- function(src, x, ...){
  if(length(x) == 0){
    return(character(0))
  }
  query <- "SELECT tax_id, name_txt FROM names WHERE name_class == 'scientific name' AND tax_id IN (%s)"
  query <- sprintf(query, sql_integer_list(x))
  tbl <- sql_collect(src, query)
  as.character(tbl$name_txt[match(x, tbl$tax_id)])
}