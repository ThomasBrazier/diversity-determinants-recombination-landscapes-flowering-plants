#!/bin/sh

# Retrieve centromere position in a gff file
# Arguments:
# Species
# Accession
# Chromosome

# GFF files
# Fields must be tab-separated. Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.'
# 
#     seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
#     source - name of the program that generated this feature, or the data source (database or project name)
#     feature - feature type name, e.g. Gene, Variation, Similarity
#     start - Start position of the feature, with sequence numbering starting at 1.
#     end - End position of the feature, with sequence numbering starting at 1.
#     score - A floating point value.
#     strand - defined as + (forward) or - (reverse).
#     frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
#     attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.



# Set your own directory
path=/Users/tbrazier/ownCloud/PhD/Analyses/data # Private path
cd "${path}"/Genome

# Metadata
# species="arabidopsis_thaliana"
# accession="GCA_902460315.1"
# chromosome=2
# Metadata
species=$1
accession=$2
chromosome=$3

#-----------------------------------------------------------------------------
# Find in annotation
#-----------------------------------------------------------------------------
dir=$(echo $species/$accession)
cd $dir/gff3
# Retrieve filename
filename=$(ls | grep $(echo "chromosome."$chromosome".gff3.gz"))
gunzip -k $filename
filename=$(sed 's/.gz//g' <<< $filename)
# hits="$(cat $filename | grep -a -e "[Cc]entromer" -e "[Kk]inetochore" | tee /dev/tty)"
# Multiple hits return multiple lines
if (cat $filename | grep -c -e "[Cc]entromer" -e "[Kk]inetochore")>0
    then (cat $filename | grep -a -e "[Cc]entromer" -e "[Kk]inetochore") | while read -r line ; do
    # echo "Processing ${line[0]}"
    # Retrieve values in hits
    link="NA"
    chr=$(echo $line | awk '{print $1}') # chromosome and start/end positions
    start=$(echo $line | awk '{print $4}') # chromosome and start/end positions
    end=$(echo $line | awk '{print $5}') # chromosome and start/end positions
    id=$(grep -o ';description=\(.*\)]' <<< $line | sed 's/;description=//' | sed 's/ /_/g')
    # Return a dataframe with hits formatted
    # Append to the local results in 
    echo "$link $chr $start $end $id"
done
fi

# (cat $filename | grep -a -e "[Cc]entromer" -e "[Kk]inetochore") | while read -r line ; do
#     # echo "Processing ${line[0]}"
#     # Retrieve values in hits
#     link="NA"
#     chr=$(echo $line | awk '{print $1}') # chromosome and start/end positions
#     start=$(echo $line | awk '{print $4}') # chromosome and start/end positions
#     end=$(echo $line | awk '{print $5}') # chromosome and start/end positions
#     id=$(grep -o ';description=\(.*\)]' <<< $line | sed 's/;description=//' | sed 's/ /_/g')
#     # Return a dataframe with hits formatted
#     # Append to the local results in 
#     echo "$link $chr $start $end $id"
# done
# hits=$(cat $filename | grep "[Kk]inetochore")
# cat $filename | grep "kinetochore"
# Retrieve values in hits
# link="NA"
# chr=$(echo $hits | awk '{print $1}') # chromosome and start/end positions
# start=$(echo $hits | awk '{print $4}') # chromosome and start/end positions
# end=$(echo $hits | awk '{print $5}') # chromosome and start/end positions
# id=$(grep -o ';description=\(.*\)gene_id=' <<< $hits | sed 's/;description=//' | sed 's/;gene_id=//' | sed 's/ /_/g')
# # Return a dataframe with hits formatted
# # Append to the local results in 
# echo "$link $chr $start $end $id"
# Save space and remove unzipped gff files
# gzip $filename
rm *.gff3