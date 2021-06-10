chicane.results.rdcl <- 
  fread("~/radicl_seq_analysis/res/no_zero_int/chicane.results.rdcl.txt")


#how many significant results
chicane.results.rdcl[q.value<0.05,]
unique(chicane.results.rdcl$target.chr)

#----------- Analysising results
#following Example Workflow (T-47D breast cancer cell line)

# save full interaction calls table
write.table(
  chicane.results.rdcl, 
  file = '~/radicl_seq_analysis/res/no_zero_int/interaction_calls.txt', 
  row.names = FALSE, 
  quote = FALSE, 
  sep = '\t'
)

# filter for significant only and save these results separately
significant.results <- chicane.results.rdcl[q.value < 0.05]

write.table(
  significant.results, 
  file = '~/radicl_seq_analysis/res/no_zero_int/interaction_calls_significant.txt', 
  row.names = FALSE, 
  quote = FALSE, 
  sep = '\t'
)

# calculate proportions by interaction type
total <- nrow(significant.results)
trans.prop <- sum(is.na(significant.results$distance))/total
cis.prop <- sum(!is.na(significant.results$distance))/total
b2b.prop <- sum(significant.results$bait.to.bait)/total
int.data <- c('trans' = trans.prop, 'cis' = cis.prop, 'bait.to.bait' = b2b.prop)
print(int.data)

# get number of interactions by distance bins
binned.data <- NULL

distance.bins <- list(
  "0-10kbp" = c(0, 1e4),
  "10-100kbp" = c(1e4, 1e5),
  "100kbp-1Mbp" = c(1e5, 1e6),
  "1-10Mbp" = c(1e6, 1e7),
  "10-100Mbp" = c(1e7, 1e8),
  "100-1000Mbp" = c(1e8, 1e9)
)

for(dist.i in names(distance.bins)){
  bin.sum <- 
    length(which(significant.results$distance >= distance.bins[[dist.i]][1] & 
                    significant.results$distance < distance.bins[[dist.i]][2]))
  binned.data <- cbind(binned.data, bin.sum)
}

colnames(binned.data) <- names(distance.bins)
print(binned.data)

# make a browser compatible file from significant interactions
browser.file <- tempfile('significant_calls_standard.txt')
#function doesn't exist................
chicane:::create.standard.format(
  chicane.results = '~/radicl_seq_analysis/res/no_zero_int/interaction_calls_significant.txt', 
  file.name = browser.file
)


