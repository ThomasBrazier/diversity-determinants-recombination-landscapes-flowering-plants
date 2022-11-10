# Analyses README

This repository contains the code necessary to reproduce results in the paper "Diversity and determinants of recombination landscapes in flowering plants" published at PLOS Genetics.

https://doi.org/10.1371/journal.pgen.1010141

Please ask if you have any question on how to understand or re-use part of this code.


## Abstract

During meiosis, crossover rates are not randomly distributed along the chromosome and their location may have a strong impact on the functioning and evolution of the genome. To date, the broad diversity of recombination landscapes among plants has rarely been investigated and a formal comparative genomic approach is still needed to characterize and assess the determinants of recombination landscapes among species and chromosomes. We gathered genetic maps and genomes for 57 flowering plant species, corresponding to 665 chromosomes, for which we estimated large-scale recombination landscapes. We found that the number of crossover per chromosome spans a limited range (between one to five/six) whatever the genome size, and that there is no single relationship across species between genetic map length and chromosome size. Instead, we found a general relationship between the relative size of chromosomes and recombination rate, while the absolute length constrains the basal recombination rate for each species. At the chromosome level, we identified two main patterns (with a few exceptions) and we proposed a conceptual model explaining the broad-scale distribution of crossovers where both telomeres and centromeres play a role. These patterns correspond globally to the underlying gene distribution, which affects how efficiently genes are shuffled at meiosis. These results raised new questions not only on the evolution of recombination rates but also on their distribution along chromosomes.

## Description

The majority of scripts are in R and R markdown. The `reports/Results_final.Rmd`generate the output files with all results and many more tests and figures once all data has been generated. As it can be a heavy task to reproduce all results, feel free to ask intermediary files and precisions about any step.

The most important data are in `data-cleaned/marey_maps` which are the cleaned Marey maps used to estimate recombination rates at a 100 kb resolution using the Marey map approach described in the paper. These maps can be re-used to re-estimate recombination rates with your own methods. The recombination landscapes estimated are in `data-cleaned/recombination_maps`.

The original genetic and marey maps (raw) are in `data`.


