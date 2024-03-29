# Chicane analysis

## results
Results for different runs. Chicane offers functionality around including zero 
count interactions in the model. As the default, only interactions that are 
detected at least once are included in the data. To also include zero counts, 
use the ‘include.zeros’ parameter: include.zeros = ‘cis’ includes all zero 
counts for bait/target combinations on the same chromosome; alternatively, 
include.zeros = ‘all’ includes all possible combinations. 

Only the default run using no zero interactions was considered. All results are
stored on the shared drive at:
/rds/general/project/neurogenomics-lab/live/Projects/radicl_seq/data/

* chicane.results.rdcl.txt - Chicane returned table from a run ont eh radicl-seq
data. This include no zero interactions with all other parameters kept at 
default.

* data_analy_inter_dist.RData - data produced when analysing the distance from 
the Transcriptional Start Site (TSS) of the noted radicl-seq interactions.

* ChromHMM_state_results.RData - data produced when using ChromHMM, Roadmap 
Adipose cell reference data to inspect the different regions the interactions 
were identified in

Input files for the Chicane run:

* frags_with_baits_cis.bed  - frag file for input to Chicane. For use with cis 
zero interactions.

* frags_with_baits_no_X_Y_M.bed - no zero interactions, RNA (baits) on X/Y 
chrom and Mitchondrial RNA removed. Whole genome split, including RNA.

* all_gene_location_no_X_Y_M.txt - RNA (baits) on X/Y chrom and Mitchondrial RNA
removed. All other RNA in the genome.

* 13.bam - Bam file for the RADICL-Seq experiment interactions.

* all_genes_location.txt - All RNA in the genome. Used as baits file for chicane

* frags_with_baits.bed - frag file for input to Chicane. For use with no zero 
interactions.


## scripts
R scripts used to run the analysis:

* **remove_overlapping_regions_rnas.R** - script and function to remove overlapping 
regions from RNA. Note this will remove from the end of the RNA pair in question so 
the transcriptional start site of the other RNA remains intact. Note this will also 
remove the RNA which are completely encompassed by another RNA. This function is 
useful for RADICL-Seq data when using CHICANE to analysis results: RNA are the baits 
file.

* **generate_frag_bed_file.R** - script to produce frags.bed, a file containing 
the whole genome segmented into random parts. This file is used as the fragment 
parameter in chicane.

* **generate_frag_bed_file_with_baits.R** - script to produce frags_with_baits.bed, 
a file containing the whole genome segmented into random parts excluding the RNA
baits and then adding on all RNA baits to the file. This file is used as the 
fragment parameter in chicane.

* **gen_frag_baits_wo_overlap.R** - script used to generate the input files which do not have
overlapping regions. Used for running chicane_wo_overlap.R.

* **running_chicane_wo_overlap.R** -script used to run chicane with input files which have
overlapping regrions removed from them.

* **running_chicane.R** - script used to complete chicane interaction calling for 
the three analysis types.

* **analysing_chicane_res.R** - early attempt at analysing the results from 
chicane interaction calling based on the qc files.

* **new.read.bed.R** - Update to the chicane read.bed() down to a bug in the 
author's code

* **ChromHMM_chrom_state_analysis.R** - Used ChromHMM, Roadmap Adipose cell 
reference data to inspect the different regions the interactions were identified
in

* **analy_inter_dist.R** - Analysing the distance from the Transcriptional Start
Site (TSS) of the noted radicl-seq interactions.

* **coding_non_coding_analysis** - Analysing trends for interactions based on 
the RNA type (coding, non-coding and subtypes)

## chicane_0.1.5.tar.gz
Zipped chicane package which was requested from the author to deal with a bug in the current CRAN release (11 June 2021)
