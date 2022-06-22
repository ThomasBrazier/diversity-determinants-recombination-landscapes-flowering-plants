# get chromosome names
gzcat *.fna.gz | grep '^>'
# Index which chromosome names to retain (trim scaffolds)
# Get sequence length
pip install pyfaidx
pip install BioPython
gunzip *.fna.gz
bgzip *.fna
faidx *.fna.gz -i chromsizes
# Trim chromosome names and sequence length vectors by index