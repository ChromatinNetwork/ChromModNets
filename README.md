# ChromModNets
Infer chromatin state-specific co-location networks from a set of ChIP-Seq or similar experiments. Detailed description of the method is available at Perner et al (http://nar.oxfordjournals.org/content/42/22/13689) and the implementation for chromatin segments is described at Juan et al (http://biorxiv.org/content/early/2015/08/05/008821). BAM files with the alignments are available at http://epistemnet.bioinfo.cnio.es/bam_files

## Usage:
Rscript main.R data/bamfilelist.txt
## Main.R
Starting point/usage-example for the analysis. Runs the counting, the Input normalization and the network inference on a list of matched Sample-Input bam files.

### Input
Several arguments are pre-set in main.R. Especially, make sure to have enough space available in outputdir and to set the correct number of available cores in ncores.

#### Command line arguments 
- path to list of bam-files (see data/bamfilelist.txt for an example)

#### Pre-set arguments in main.R
- outputdir: all intermediate and output files will be written here
- chromBins_filename: path to the file containing the annotation of the chromatin state ranges, last colum = chromatin state of bin (see data/all_chromstate_bins.txt for an example)
- mapInputSample_filename: filename of file in which each line connects hms/cms name, sample file name and corresponding input file name (see data/inputSampleMap.txt for an example)
- nrcores: the number of available cores for parallelization
- sdf: cutoff applied to the edge weights from Elastic net. A factor for multiplying the standard deviation of the features
- hms/cms: files containing the hm and cm names (see data/hmnames.txt and data/cmnames.txt for an example)

### Output
- counts-folder: intermediate bed-files containing summarized counts per chromatin segment per hms/cms sample and input (make sure there is enough space for this)
- all_data_norm.RData: the normalized count matrix
- chromnet.txt: the final chromatin state-specific networks
- supptable_enweights.xls: original edge-weights retrieved for each pair of hms/cms in each state in the Elastic Net
- supptable_pcor.xls: original edge-weights retrieved for each pair of hms/cms in each state in the SPCN

## Function.R
Contains all functions to run the network inference.
