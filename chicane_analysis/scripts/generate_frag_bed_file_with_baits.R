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
#baits holds all rna - needs to be added to frags file
baits <- fread("~/radicl_seq_analysis/all_genes_location.txt")

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
    #bait.chr==target.chr & bait.start<=target.end & bait.end>=target.start
    if(nrow(baits[V1==chr_num & V2<= rng_end[i] & V3>=rng_strt[i]])>0){
      rng_strt <- rng_strt[rng_strt!=rng_strt[i]]
      rng_end <- rng_end[rng_end!=rng_end[i]]
      i <- i - 1
    }
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
#combine to make frag bed file
frags <- rbindlist(rand_seg_chr)
#lastly, add baits to fragment file
frags <- rbindlist(frags,baits)
#write to bed file
fwrite(frags,file="~/radicl_seq_analysis/frags_with_baits.bed",
        col.names = FALSE)