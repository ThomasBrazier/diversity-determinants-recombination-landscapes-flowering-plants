# Analyses README

## Description


## Nomenclature

### File naming

All file names must be in snake_case.

Directories are named in lowercases.

Files are named in lowercases with uppercase for first letter of a word.

## Pipeline of analyses

### Structure and distribution of the recombination at large scale

Marey maps

### Genomic landscapes

'data-cleaned/genome/gene_positions' need to be computed in 'Genome_Architecture.R' first as it is the source for:
* Gene count in recombination maps ('gene_count()' in 'Gene_density.R')
* Estimated mean recombination rate per gene count, by 'estimmeanrecrat0e_genecount()' in 'Gene_density.R'
* Marey maps in gnee distance, by 'gene_distance.map()' in 'Genetic_Shuffling.R'

