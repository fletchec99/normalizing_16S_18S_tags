# Purpose:
This script merges 16S and 18S ASV tables generated through the qiime2 pipeline modified by Jesse McNichol by normalizing raw ASV counts against

1. The percent of reads passing DADA2, and
2. The sequencing platform's bias against longer 18S sequences, as determined empirically by the investigator.

# Inputs:
- 16S and 18S ASV tables with taxonomy ("02-PROKs/10-exports/all-16S-seqs.with-tax.tsv" and "02-EUKs/15-exports/all-18S-seqs.with-PR2-tax.tsv", respectively).
- User must specify the bias against 18S sequences. We have found that a bias of 2x correctly normalizes our mixed mock communities of both 16S and 18S sequences (Yeh et al., 2018). 2x is therefore the default bias against 18S sequences in this script. However, samples from an environment with a higher proportion of long 18S sequences may have a larger bias against 18S.
- User must also specify the prefix of the output file names. 

# Outputs:
- Two .tsv files will be generated, the normalized counts of 16S and 18S sequences and the counts converted to proportions (out of 16S+18S total sequences).

# Usage: 
`$ Rscript normalize_16S_18S.R --bias <number> --outputfile <name>`

# Dependent packages: 
- This script requires Rscript, base R, and optparse. 
- Dependent packages may be found within the conda .yml file in this repository.

A few conda environment basics for those who aren't familiar: 
- To create a conda environment from the .yml file: `$ conda env create --file R-env.yml`
- To activate the conda environment before running the script: `$ conda activate R-env`
- To deactivate the conda environment after running the script: `$ conda deactivate`
- For more information: https://conda.io/projects/conda/en/latest/user-guide/getting-started.html

# Version: 05.10.2022
# Author: Colette Fletcher-Hoppe, with help from Jesse McNichol. 
