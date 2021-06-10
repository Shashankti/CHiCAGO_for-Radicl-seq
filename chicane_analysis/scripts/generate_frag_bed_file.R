#rsplit is the rmap file which was previoulsy generated, basically the file 
#which was the columns4-6 from the bedpe file. I had to do this because I 
#couldnt find a way to generate random numbers as a sequence such that the 
#difference between consecutive terms is between 16-46. Please edit this if 
#you know a better way.
library(data.table)
# chrsize holds bp size for each chromosome
chrsize <- fread("~/radicl_seq_analysis/D13.sizes")
setnames(chrsize,c("chr","size"))
#exclude X and Y
chrsize <- chrsize[!chr %in% c("X","Y"),]
#make sure columns are numeric
chrsize[,chr:=as.numeric(chr)]
chrsize[,size:=as.numeric(size)]
#set key for qquick joins and reorders...
setkey(chrsize,chr)

#get rand numbers that add to chrom length
get_nums <- function(lower,upper,chrom_len,chr_num){
  nums <- vector()
  rng_strt <- vector()
  rng_end <- vector()
  sum <- 0
  i<-1
  #generate random numbers
  #worst case smallest sum is all lowers so only need to gen that many rand numbers
  rands <- sample(16:40,ceiling(chrom_len/lower),replace = TRUE)
  while(sum<chrom_len){
    #get start postion
    rng_strt[i] <- sum+1
    sum <- sum + rands[i]
    #get end postion
    rng_end[i] <- sum
    i <- i + 1
  }
  diff <- sum-chrom_len
  rng_end[i-1] <- rng_end[i-1]-diff
  return(data.table(chr=rep(paste0("chr",chr_num),length(rng_end)),
                    start=rng_strt,
                    end=rng_end))
}

#set seed so you can get same results again
set.seed(101)
rand_seg_chr <- 
  lapply(seq_along(chrsize$size),function(x) get_nums(16,40,chrsize$size[x],x))
#sense check to ensure end of last segment same length as chromosome
all.equal(unlist(lapply(rand_seg_chr,function(x) x[nrow(x),]$end)),
          chrsize$size)
#TRUE
#finally combine to make frag bed file
frags <- rbindlist(rand_seg_chr)
#write to bed file
fwrite(frags,file="~/radicl_seq_analysis/frags.bed",col.names = FALSE)
