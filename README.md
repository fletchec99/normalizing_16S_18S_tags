# Purpose:
This script merges 16S and 18S ASV tables generated through the [qiime2 pipeline](https://github.com/jcmcnch/eASV-pipeline-for-515Y-926R) modified by Jesse McNichol by normalizing raw ASV counts against:

1. The percent of reads passing DADA2's filters (accounting for variations likely due to quality differences between samples), and
2. The sequencing platform's bias against longer 18S sequences, as determined empirically by the investigator (accounting for instrument-specific biases against longer 18S sequences as described in [Yeh et al., 2021](https://sfamjournals.onlinelibrary.wiley.com/doi/10.1111/1462-2920.15553)).

# Inputs:
- 16S and 18S ASV count tables with taxonomy, as output by default from [the aforementioned pipeline](https://github.com/jcmcnch/eASV-pipeline-for-515Y-926R). *Do not* use tables that have been converted to relative abundance!
- Statistics for how many of the initial reads pass through the DADA2 pipeline (also automatically produced by the pipeline indicated above).
- A floating-point number representing the bias against 18S sequences. We have found that a bias of 2x correctly normalizes our mixed mock communities of both 16S and 18S sequences ([Yeh et al., 2021](https://sfamjournals.onlinelibrary.wiley.com/doi/10.1111/1462-2920.15553)). 2x is therefore the default bias against 18S sequences in this script. However, samples from an environment with a higher proportion of long 18S sequences may have a larger bias against 18S. This may also be taxa-dependent. For an example of how to do this calculation see [this repository](https://github.com/jcmcnch/P16N-S_2005-2006_CMAP/tree/main/18S-correction).
- User must also specify the prefix of the output file names.

# Outputs:
- Two .tsv files will be generated, the normalized counts of 16S and 18S sequences and the counts converted to proportions (out of 16S+18S total sequences).

# Usage: 

```
Rscript normalize_16S_18S.R --bias <float> \
	--inputproks <input prok count table, with taxonomy> \
	--prokstats <a tab-separated stats file for 16S> \
  	--inputeuks <input euk count table, with taxonomy> \
	--eukstats <a tab-separated stats file for 18S> \
	--outputfile <name>
```

# Dependent packages: 
- This script requires Rscript, base R, and optparse. 
- Dependent packages may be found within the conda .yml file in this repository.

A few conda environment basics for those who aren't familiar: 
- To create a conda environment from the .yml file: `conda env create --file R-env.yml`
- To activate the conda environment before running the script: `conda activate R-env`
- To deactivate the conda environment after running the script: `conda deactivate`
- For more information: https://conda.io/projects/conda/en/latest/user-guide/getting-started.html

Version: 05.10.2022
Author: Colette Fletcher-Hoppe, with help from Jesse McNichol. 
