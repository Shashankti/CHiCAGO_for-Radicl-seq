#RADICL-Seq data info
#
#We have data from NiGa cells at day 0 and day 13 post-differentiation. NiGa 
#cells are a cell line isolated from subcutaneous abdominal white adipose tissue
#from a 16 year old male donor (https://academic.oup.com/jcem/article/98/3/E503/2536944).
#The D0 cells are nondifferentiated human adipose-derived stem cells (hADSCs). 
#The D13 cells are in vitro differentiated adipocytes.
#D13 data used
#
#The raw data contained the following columns.
#Col1: Chromosomal origin of the RNA tag
#Col2: Genomic location of the RNA tag start
#Col3: Genomic location of the RNA tag end
#Col4: Chromosomal origin of the DNA tag
#Col5: Genomic location of the DNA tag start
#Col6: Genomic location of the DNA tag end
#Col7: Cigar
#Col8: ignore
#Col9: RNA strand
#Col10: DNA strand (DNA ligation is not stranded so please ignore it).
#Col11: Map quality of RNA (STAR)
#Col12: Map quality of DNA (bwa)
#Col13: RNA ID

# removed the col7-12
#Also for Day13 there were 3 different files. I merged them into one
#The files were named like this: READme.rtf uniq.AGTTCC_D0_3.bedpe   
#uniq.CCGTCC_D13_2.bedpe  uniq.TGACCA_D13_3.bedpe uniq.AGTCAA_D0_1.bedpe  
#uniq.ATGTCA_D13_1.bedpe  uniq.CTTGTA_D0_2.bedpe. 
#So I assumed they represented data from same runs and same day, just split 
#keeping the size reasonabble.
#Q: is unique taken for PCR amplicifcation?? 

#Drawbacks:
#only RNA's which are bound considered
#only DNA which is bound considered
#RADICL-Seq is an all to all not a many to all like CHi-C


#chicane paper info
#CHi-C approaches generate data matrices representing unbalanced many-to-all 
#interaction profiles between baits and target fragments (also referred to as a 
#non-captured end), respectively, where baits are pre-defined but the target 
#fragments could be anywhere in the genome.
#
#Using CHiCANE, BAM files are processed into interaction tables by removing all 
#reads where neither end maps to a capture bait as well as those in which both 
#ends map to the same fragment.
#
#As the default, only interactions that are detected at least once are included
#in the data. To also include zero counts, use the ‘include.zeros’ parameter: 
#include.zeros = ‘cis’ includes all zero counts for bait/target combinations on 
#the same chromosome; alternatively, include.zeros = ‘all’ includes all possible
#combinations. The example below shows inclusion of cis zeros
#
#CHiCANE tries to split data into an optimal number of distance bins based on 
#the size of the dataset, such that the resulting datasets are large enough for
#the model to be fit. In order to override this default behavior, set the 
#desired number of distance bins with the ‘distance.bins’ parameter. For 
#example, the following call will result in the count model fitted separately in
#each of the 50 data bins (trans interactions are fitted separately from cis 
#interactions).

#setwd("~/radicl_seq_analysis/")
#install.packages("chicane") #doesn't work on linux ...
#install.packages("~/radicl_seq_analysis/chicane_0.1.5.tar.gz", 
#                  repos = NULL, type="source")
library(chicane)

#--------- vignette data
bam <- system.file('extdata', 'Bre80_2q35.bam', package = 'chicane');
baits <- system.file('extdata', '2q35.bed', package = 'chicane');
fragments <- system.file('extdata', 'GRCh38_HindIII_chr2.bed.gz', 
                         package = 'chicane'); # HindIII fragments on chromosome 2
bedtools.installed()
#inspecting files
bait.ids2 <- read.bed(baits)#[1] "chr2:217035649-217042355" "chr2:217430817-217437011"
fragment.ids2 <- read.bed(fragments)#[1] "chr2:0-15298"         "chr2:15298-16049"     "chr2:16049-18132" 
#--- look into fragment size
fragment.ids2dt <- data.table(id=fragment.ids2)
fragment.ids2dt[,chr:=gsub('(.*):(.*)', '\\1', id)]
fragment.ids2dt[,range:=sub(".*:", "", id)]
fragment.ids2dt[,start:=as.numeric(sub("\\-.*", "", range))] 
fragment.ids2dt[,end:=as.numeric(sub(".*-", "", range))]
fragment.ids2dt[,frag_length:=end-start]
#---

chicane.results <- chicane(
  bam = bam,
  baits = baits,
  fragments = fragments
);

print( chicane.results[ 1:10 ] );


#------------ chicane RADICL-Seq run
#outline
#Use all possible RNA as baits, not just those with an interaction
#Use all genome as fragments, not just those with an interaction
#don't create zero interactions*
#
# input files
# 13.bam - bam file from radicl-seq experiment containing noted RNA-DNA 
# interactions for the three replicates
# all_genes_location.txt - bp location of all genes in the genome of which a 
# subset were found to be interacting in 13.bam
# frags.bed - the whole genome randomly segmented into 16-45 bp segments. A 
# subset of these fragments will overlap with the other end fragments in the 
# 13.bam

bam_rdcl <- "~/radicl_seq_analysis/13.bam"
baits_rdcl <- "~/radicl_seq_analysis/all_genes_location.txt"
fragments_rdcl <- "~/radicl_seq_analysis/frags.bed"

#bait.ids_rdcl <- read.bed(baits_rdcl)
#fragment.ids_rdcl <- read.bed(fragments_rdcl)

chicane.results.rdcl <- chicane(
  bam = bam_rdcl,
  baits = baits_rdcl,
  fragments = fragments_rdcl,
  interim.data.dir= "~/radicl_seq_analysis/res/no_zero_int/qc/",#gives qc plots and info
  verbose=TRUE
)

fwrite(chicane.results.rdcl,
       "~/radicl_seq_analysis/res/no_zero_int/chicane.results.rdcl.txt")


#------------ chicane RADICL-Seq run cis
#outline
#Use all possible RNA as baits, not just those with an interaction
#Use all genome as fragments, not just those with an interaction
#use cis zero interactions
#
# input files
# 13.bam - bam file from radicl-seq experiment containing noted RNA-DNA 
# interactions for the three replicates
# all_genes_location.txt - bp location of all genes in the genome of which a 
# subset were found to be interacting in 13.bam
# frags.bed - the whole genome randomly segmented into 16-45 bp segments. A 
# subset of these fragments will overlap with the other end fragments in the 
# 13.bam
#try with baits missing X,Y,M chromosome
baits <- fread("~/radicl_seq_analysis/all_genes_location.txt")
baits <- baits[!V1 %in% c("chrX","chrY","chrM" ),]
baits <- fwrite(baits,"~/radicl_seq_analysis/all_genes_location_no_X_Y_M.txt",
                  col.names = FALSE)

bam_rdcl <- "~/radicl_seq_analysis/13.bam"
baits_rdcl <- "~/radicl_seq_analysis/all_genes_location.txt"
#use frags with baits added- no overlap, see generate_frag_bed_file_with_baits.R
fragments_rdcl <- "~/radicl_seq_analysis/frags_with_baits.bed"

#bait.ids_rdcl <- read.bed(baits_rdcl)
#fragment.ids_rdcl <- read.bed(fragments_rdcl)

#need to edit read.bed as issue with scientific notation - see new.read.bed.R
environment(new.read.bed) <- asNamespace('chicane')
assignInNamespace("read.bed", new.read.bed, ns = "chicane")

chicane.results.rdcl <- chicane(
  bam = bam_rdcl,
  baits = baits_rdcl,#"~/radicl_seq_analysis/all_genes_location_no_X_Y_M.txt"
  fragments = fragments_rdcl,
  interim.data.dir= "~/radicl_seq_analysis/res/cis_zero_int/qc/",#gives qc plots and info
  include.zeros = "cis",
  verbose=TRUE
)

fwrite(chicane.results.rdcl,
       "~/radicl_seq_analysis/res/cis_zero_int/chicane.results.rdcl.txt")
