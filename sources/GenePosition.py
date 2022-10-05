# -*- coding: utf-8 -*-
"""
Retrieve Gene positions in a sequence

Retrieve gene positions in a chromosome sequence coded in .gff
Return a data frame with gene id, start, end positions

Arguments:
    path = absolute path where to find genome sequence of the chromosome; must specify a gff file
    
Values:
    Return a data frame with three columns: gene id, start position, end position
"""


def GenePosition(path = ""):
    from BCBio import GFF
    # from Bio import SeqIO
    import re
    import pandas as pd
    
    #-------------------------------------------------------------
    # Import file and declare function + global variables
    in_file = "aegilops.gff3"
    # Open file
    in_handle = open(in_file)
    
    # Return the hierarchy in the GFF
    # pprint.pprint(examiner.parent_child_map(in_handle))
    
    # Return all features with counts
    # pprint.pprint(examiner.available_limits(in_handle))
    
    #-------------------------------------------------------------
    # returns a set of SeqRecord objects corresponding to the various IDs referenced in the file 
    # Limit to 'gene' features in gff_type
    limit_info = dict(
        gff_type = ["gene"])
    for rec in GFF.parse(in_handle, limit_info = limit_info):
        # print(rec.features)
        hits = rec.features
    
    # hits = GFF.parse(in_handle, limit_info = limit_info)
    # print(hits[0])
    
    # Init lists of results
    ids = []
    start = []
    end = []
    
    #-------------------------------------------------------------
    # Get a list of ids, start, end positions of all genes in the chromosome
    for i in range(0,len(hits)):
        # print(i)
        # print(hits[i])
        # Retrieve gene id
        record = hits[i]
        ids.append(re.sub("gene:", "", record.id))
        # Retrieve start & end positions
        # Given the +/- strand, start & end positions may be reversed
        # hence, start must always be the min of the two given positions
        # retrieve both positions in a list of two elements
        loc = str(record.location)
        match = re.findall("[0-9]+", loc)
        start.append(min(match))
        end.append(max(match))
    
        # print("-----------------------")
    
    # Return a data frame in a file
        data = {
        'geneID': ids,
        'start': start,
        'end': end
        }
    
        df = pd.DataFrame(data)
    
    #-------------------------------------------------------------
    # Close file
    in_handle.close()
    
    # return the data frame in terminal
    return df
