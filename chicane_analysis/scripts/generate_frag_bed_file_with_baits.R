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
    baits_ovrlp <- baits[chr==chr_num & start<= rng_end[i] & end>=rng_strt[i]]
    if(nrow(baits_ovrlp)>0){
      rng_strt <- rng_strt[rng_strt!=rng_strt[i]]
      rng_end <- rng_end[rng_end!=rng_end[i]]
      i <- i - 1
      #we can skip the length of the gene
      sum <- max(baits_ovrlp$end)
    }
    i <- i + 1
  }
  print("chrom end")
  diff <- sum-chrom_len
  rng_end[i-1] <- rng_end[i-1]-diff
  rtrn <- data.table(chr=rep(chr_num,length(rng_end)),
             start=rng_strt,
             end=rng_end)
  return(rtrn)
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
frags <- rbindlist(list(frags,baits[,-4]))
#get correct order
frags[,ordr:=as.numeric(substr(chr,4,nchar(chr)))]
setorder(frags,ordr)
frags[,ordr:=NULL]
#write to bed file
fwrite(frags,file="~/radicl_seq_analysis/frags_with_baits.bed",
       col.names = FALSE)

                 
                 
#----------------- Alternative approach using multiple Cores for speed!
                 
#parallelise run for speed
library(parallel)
library(doParallel)
number_threads <- detectCores()-2
cl <- parallel::makeCluster(number_threads)
doParallel::registerDoParallel(cl)
set.seed(101)
frags <- 
  foreach::foreach(x = seq_along(chrsize$size)[!seq_along(chrsize$size) %in% c(1,2,3)],
                    .packages = c("data.table","stats")) %dopar% {
                      chr_num <- paste0("chr",x)
                      lower <- 16
                      upper <- 40
                      chrom_len <- chrsize$size[x]
                      nums <- vector()
                      rng_strt <- vector()
                      rng_end <- vector()
                      sum <- 0
                      i<-1
                      #generate random numbers
                      #worst case smallest sum is all lowers so only need to gen that many rand numbers
                      rands <- sample(lower:upper,ceiling(chrom_len/lower),replace = TRUE)
                      while(sum<chrom_len){
                        #get start postion
                        rng_strt[i] <- sum+1
                        sum <- sum + rands[i]
                        #get end postion
                        rng_end[i] <- sum
                        #bait.chr==target.chr & bait.start<=target.end & bait.end>=target.start
                        baits_ovrlp <- baits[chr==chr_num & start<= rng_end[i] & end>=rng_strt[i]]
                        if(nrow(baits_ovrlp)>0){
                          rng_strt <- rng_strt[rng_strt!=rng_strt[i]]
                          rng_end <- rng_end[rng_end!=rng_end[i]]
                          i <- i - 1
                          #we can skip the length of the gene
                          sum <- max(baits_ovrlp$end)
                        }
                        i <- i + 1
                      }
                      print("chrom end")
                      diff <- sum-chrom_len
                      rng_end[i-1] <- rng_end[i-1]-diff
                      rtrn <- data.table(chr=rep(chr_num,length(rng_end)),
                                         start=rng_strt,
                                         end=rng_end)
                      return(rtrn)
                  }

parallel::stopCluster(cl)

frags <-rbindlist(frags)
#lastly, add baits to fragment file
frags <- rbindlist(list(frags,baits[,-4]))
#get correct order
frags[,ordr:=as.numeric(substr(chr,4,nchar(chr)))]
setorder(frags,ordr)
frags[,ordr:=NULL]
fwrite(frags,file="~/radicl_seq_analysis/frags_with_baits.bed",
       col.names = FALSE)                 

                 
#----------------- Alternative approach BUT NEED BIGGER FRAGMENTS SO CHICANE WONT GIVE OUT - cis
                 
#parallelise run for speed
library(parallel)
library(doParallel)
number_threads <- detectCores()-2
cl <- parallel::makeCluster(number_threads)
doParallel::registerDoParallel(cl)
set.seed(101)
frags <- 
  foreach::foreach(x = seq_along(chrsize$size)[!seq_along(chrsize$size) %in% c(1,2,3)],
                    .packages = c("data.table","stats")) %dopar% {
                      chr_num <- paste0("chr",x)
                      lower <- 16
                      upper <- 40
                      chrom_len <- chrsize$size[x]
                      nums <- vector()
                      rng_strt <- vector()
                      rng_end <- vector()
                      sum <- 0
                      i<-1
                      #generate random numbers
                      #worst case smallest sum is all lowers so only need to gen that many rand numbers
                      rands <- sample(lower:upper,ceiling(chrom_len/lower),replace = TRUE)
                      while(sum<chrom_len){
                        #get start postion
                        rng_strt[i] <- sum+1
                        sum <- sum + rands[i]
                        #get end postion
                        rng_end[i] <- sum
                        #bait.chr==target.chr & bait.start<=target.end & bait.end>=target.start
                        baits_ovrlp <- baits[chr==chr_num & start<= rng_end[i] & end>=rng_strt[i]]
                        if(nrow(baits_ovrlp)>0){
                          rng_strt <- rng_strt[rng_strt!=rng_strt[i]]
                          rng_end <- rng_end[rng_end!=rng_end[i]]
                          i <- i - 1
                          #we can skip the length of the gene
                          sum <- max(baits_ovrlp$end)
                        }
                        i <- i + 1
                      }
                      print("chrom end")
                      diff <- sum-chrom_len
                      rng_end[i-1] <- rng_end[i-1]-diff
                      rtrn <- data.table(chr=rep(chr_num,length(rng_end)),
                                         start=rng_strt,
                                         end=rng_end)
                      return(rtrn)
                  }

parallel::stopCluster(cl)

frags <-rbindlist(frags)
#lastly, add baits to fragment file
frags <- rbindlist(list(frags,baits[,-4]))
#get correct order
frags[,ordr:=as.numeric(substr(chr,4,nchar(chr)))]
setorder(frags,ordr)
frags[,ordr:=NULL]
fwrite(frags,file="~/radicl_seq_analysis/frags_with_baits_cis.bed",
       col.names = FALSE)                   
