library(chicane)

#Running chicane with the fragment file from the raw data and the overlap between the fragments and the baits file removed
# frag file-contains the otherEnds from the raw data file such that all overlapping otherEnds are reduced to remove the overlapping regions.
#the frag file is intersected with the baits file and overlapping regions are removed.
#the baits file contains overlapping regions
chicane.results <- chicane(
  bam = "Data/Processed/D13.bam",
  baits = "Data/Processed/graphs/chic_wo_baits.bed",
  fragments = "Data/Processed/graphs/frags_wo_ovrlap.bed",
  cores = 7
)
significant.results <- chicane.results[chicane.results$q.value<0.05]
fwrite(chicane.results,file = "Data/Processed/graphs/chicane.results.txt")




##Running chicane with fragment file from the raw data and the overlap within baits also removed.

chicane.results.wo.baits <- chicane(
  bam = "Data/Processed/D13.bam",
  baits = "Data/Processed/graphs/wo_overlap_baits.bed",
  fragments = "Data/Processed/graphs/frags_wo_bait_overlap.bed",
  cores = 7
)
fwrite(chicane.results.wo.baits,file = "Data/Processed/graphs/chicane.results.wo.baits.txt")

#Running chicane by usung the fragments file as a bins of 2kb for the enitre genome. The baits file contains overlapping region
#The overlap between the baits and the fragments is removed
chicane.results.binned <- chicane{
  bam = "Data/Processed/D13.bam",
  baits = "Data/Processed/graphs/chic_wo_baits.bed",
  fragments = "Data/Processed/graphs/binned.bed",
  cores=7
 }
 
#Runing chicace with binned fragment file and baits without overlaps,
chicane.results.bonned.wo <- chicane{
  bam = "Data/Processed/D13.bam",
  baits = "Data/Processed/graphs/wo_overlap_baits.bed",
  fragments = "Data/Processed/graphs/binned_wo_bait_overlap.bed",
  cores = 7
)
