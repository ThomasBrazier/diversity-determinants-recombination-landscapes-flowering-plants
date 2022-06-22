#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GENOME ARCHITECTURE
Created on Fri Jan 31 15:11:31 2020

@author: tbrazier

A global pipeline gathering informations and statistics on the genome architecture of reference genomes

ARGUMENTS:
    Take a list of accession to process as input

VALUES:
    Return multiple data frames in txt files, with the genome features  of interest 
"""

import pandas
import os
import re
# import csv
from zipfile import ZipFile

# cwd = os.getcwd()
wd = "/Users/tbrazier/ownCloud/PhD/Analyses"
os.chdir(wd)
# Custom scripts
from sources.GenePosition import GenePosition

# The directory where to save output files
outputdir = "/output/genome_architecture/"
# Import the list of accession references


# List of genome features to assess; 0 if non required/1 if computation is required
genepos = 1

# Read the list of accessions to process
list_process = open("".join([cwd, "/data/Genome/ref_accessions_genomearchitecture.txt"]), 'r')

# For each accession, retrieve the path to the gff files
# Species, accession and chromosomes names are required
# All infos are given in the 'Genome_ressources.csv'
# The four subsequent vectors are linked, each position i targeting the same dataset 
species = []
accession = []
chromosome = []
gff = []
metadata = pandas.read_csv("".join([cwd, "/data/Genome/Genome_ressources.csv"]), sep = ";")
maps = pandas.read_csv("".join([cwd, "/data-cleaned/Marey_maps/AllMaps.txt"]), sep = ";")

list_process = list_process.readlines()
list_process = [line.rstrip('\n') for line in list_process]

for l in list_process:
    # print(l)
    # ACCESSION
    acc = l
    print(acc)
    # SPECIES
    sp = metadata.loc[metadata['accession'] == acc]
    sp = sp.iat[0,2]
    sp = re.sub("_", " ", sp)
    print(sp)
    # CHROMOSOMES

    # Retrieve list of chromosomes
    # gff files are save in .gff.gz with the pattern 'chromosome.[A-Za-z0-9]+.' in the directory /gff3
    # Get the list of gff files and retrieve all chromosome patterns
    
    # Extract only names in a list 'ch'
    # ch = 
    #
    # species.append(sp * len(ch))
    # accession.append(acc * len(ch))
    # chromosome.append(ch)

# Init results with empty text file
if genepos == 1:
    genepos_file = open("".join([outputdir, "gene_position.txt"]), "w")

# Then loop on all datasets and append results to files
for i in 0:len(species):
    s = species[i]
    a = accession[i]
    c = chromosome[i]
    # Unzip gff file
    zf = ZipFile('path_to_file/file.zip', 'r')
    zf.extractall('path_to_extract_folder')
     
    # Retrieve gene positions for reference genomes in the list
    if genepos == 1:
        results = GenePosition(path)
        # Add species, accession and chromosome columns to the results
        sp_col = [species] * results.count
        acc_col = [accession] * results.count
        chr_col = [chromosome] * results.count
        # Concatenate columns
        results = pandas.concat([sp_col, acc_col, chr_col, results], axis=1, sort=False)
        
        # Append data frame in a txt file
        genepos_file.append(results)
     
    
     
     # Remove unzipped file
     zf.close()



# END OF SCRIPT