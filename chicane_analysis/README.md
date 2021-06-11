# Chicane analysis

## results
Results for different runs. Chicane offers functionality around including zero 
count interactions in the model. As the default, only interactions that are 
detected at least once are included in the data. To also include zero counts, 
use the ‘include.zeros’ parameter: include.zeros = ‘cis’ includes all zero 
counts for bait/target combinations on the same chromosome; alternatively, 
include.zeros = ‘all’ includes all possible combinations:

* **no_zero_int** - results for analysis using default approach where no zero 
count interactions were included.

* **cis_zero_int** - results for analysis where cis zero interactions were 
included.

* **all_zero_int** - results for analysis where both cis and trans zero 
interactions were included.

## scripts
R scripts used to run the analysis:

* **generate_frag_bed_file.R** - script to produce frags.bed, a file containing 
the whole genome segmented into random parts. This file is used as the fragment 
parameter in chicane.

* **generate_frag_bed_file_with_baits.R** - script to produce frags_with_baits.bed, 
a file containing the whole genome segmented into random parts excluding the RNA
baits and then adding on all RNA baits to the file. This file is used as the 
fragment parameter in chicane.

* **running_chicane.R** - script used to complete chicane interaction calling for 
the three analysis types.

* **analysing_chicane_res.R** - early attempt at analysing the results from 
chicane interaction calling based on the qc files.

## chicane_0.1.5.tar.gz
Zipped chicane package which was requested from the author to deal with a bug in the current CRAN release (11 June 2021)
