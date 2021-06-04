library(optparse)
library(doMC)
library(ggplot2)
library(genomation)
library(GenomicRanges)
library(data.table)
library(plyr)
library(dplyr)
library(stringr)
library(foreach)
library(doParallel)
#set number of cores
library(doParallel)
myCluster <- makeCluster(3, # number of cores to use
                         type = "FORK") # type of cluster
registerDoParallel(myCluster)


option_list = list(
  make_option(c("-b", "--bedpe"), type="character", default=NULL, 
              help="Path to bedpe file", metavar="character"),
  make_option(c("-r", "--rmap"), type="character", default=NULL, 
              help="Path to rmap file", metavar="character"),
  make_option(c("-p", "--baitmap"), type="character", default=NULL, 
              help="Path to baitmap file", metavar="character"),
  make_option(c("-s", "--size"), type="character", default=NULL, 
              help="Path to chromosome sizes file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.chinput", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$bedpe)){
  print_help(opt_parser)
  stop("Enter the files", call.=FALSE)
} else if (is.null(opt$baitmap)){
  print_help(opt_parser)
  stop("Enter the files", call.=FALSE)  
} else if (is.null(opt$rmap)) {
  print_help(opt_parser)
  stop("Enter the files", call.=FALSE)
} else if (is.null(opt$size)){
  print_help(opt_parser)
  stop("Enter the files", call.=FALSE)
}
#Read the bedpe input file
yyy = fread(opt$bedpe,header=FALSE,stringsAsFactors = FALSE,quote = "")
#colnames(yyy) <- c('rnachrom','rnachromStart','rnachromEnd','dnachrom','dnachromStart','dnachromEnd','cigar','ignore','rnaStrand','dnaStrand','rnQual','dnaQua','rnaID')
colnames(yyy) <- c('rnachrom','rnachromStart','rnachromEnd','dnachrom','dnachromStart','dnachromEnd','rnaID')

#read unprocessed Rmap file
rmap <- fread(opt$rmap,header=FALSE,stringsAsFactors = FALSE,quote = "")
colnames(rmap) <- c('chrom','start','end','Frag_id')
#read baitmap file
baitmap <- fread(opt$baitmap,header=FALSE,stringsAsFactors = FALSE,quote = "")
colnames(baitmap) <- c('chrom','start','end','Frag_id','bait')
baitmap <- baitmap[with(baitmap,order(chrom,start)),]




#Read size file
chrsize <- fread(opt$size,header=FALSE,stringsAsFactors = FALSE,quote = "")
chrsize <- chrsize[order(as.numeric(as.character(chrsize$V1)))]

##############
#correcting baitmap file
yelp <- as.data.frame(reduce(GRanges(baitmap$chrom,IRanges(baitmap$start,baitmap$end))))
#yelp$Frag_id <- baitmap$Frag_id[match(yelp$start,baitmap$start)]
yelp$bait <- baitmap$bait[match(yelp$Frag_id,baitmap$Frag_id)]

baitmap <- yelp[,c("seqnames","start","end","bait","width")]
colnames(baitmap)[1] <- "chrom"
colnames(baitmap)[6] <- "len"
btry <- baitmap[,c("chrom","start","end","len")]




#Correcting the rmap file
bres <- makeGRangesFromDataFrame(baitmap)
res <- makeGRangesFromDataFrame(rmap)
ov <- findOverlaps(res,bres)
grl <- extractList(bres,as(ov,"List"))
blip <- psetdiff(res,grl)
blip <- as.data.frame(blip)
blip$group = NULL
blip$group_name = NULL
rmap <- rbind(blip,as.data.frame(bres))
colnames(rmap)[1] <- "chrom"
rmap$Frag_id <- 1:nrow(rmap)
rmap$strand = NULL
rmap <-rmap[with(rmap,order(chrom,start)),]
rsplit <- split(rmap,rmap$chrom)
rmap$width = NULL
baitmap$Frag_id = rmap$Frag_id[match(baitmap$start,rmap$start)]
baitmap <- baitmap[,c(1,2,3,5,4)]
rsplit <- split(rmap,rmap$chrom)
#final processed rmap and baitmap files
write.table(baitmap,"Data/Processed/DesFile/D13.baitmap",sep = "\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
write.table(rmap,"Data/Processed/DesFile/D13.rmap",sep = "\t",col.names = FALSE,row.names = FALSE)







#generating the Chinput file
chkl <- split(yyy,yyy$dnachrom)
rm(yyy)
gc()
bin = NULL
# Bin the data 
for (i in 1:22){
  bin[[i]] <- transform(chkl[[i]], group = cut(dnachromStart,
                                               breaks=rsplit[[i]]$end))
}
rm(chkl)
gc(reset = TRUE)

for (i in 1:length(bin)) {
  bin[[i]] <-transform(bin[[i]], Freq = ave(seq(nrow(bin[[i]])),group, FUN = length))
  bin[[i]] <- bin[[i]][,-c(2:3)]
}

#Defining function to check integer(0)
is.integer0 <- function(x){
  is.integer(x) && length(x) == 0L
}

#getting otherendID column and distal length

k = NULL
value= NULL
otherID = NULL


value <- foreach( k =1:22 ) %:% 
  foreach ( i = 1:length(bin[[k]]$rnachrom), .combine = 'c'  ) %dopar% {
    ifelse(is.integer0(rsplit[[k]]$Frag_id[which(rsplit[[k]]$start <= bin[[k]]$dnachromStart[i] & rsplit[[k]]$end >= bin[[k]]$dnachromStart[i])]),NA,rsplit[[k]]$Frag_id[which(rsplit[[k]]$start <= bin[[k]]$dnachromStart[i] & rsplit[[k]]$end >= bin[[k]]$dnachromStart[i])])
  }
for(i in 1:22){
  otherID[[i]] <- value[[i]]
}
rm(value)
gc()

#######################
for (i in 1:22){
  bin[[i]]$otherID <- otherID[[i]]
}
#Defining otherEndlen,distSign
rmap$len <- rmap$end-rmap$start
all_chr_freq <- bind_rows(bin, .id = "column_label")
all_chr_freq <- all_chr_freq[with(all_chr_freq,order(as.numeric(column_label),as.numeric(group))),]

all_chr_freq$baitID <- baitmap$Frag_id[match(all_chr_freq$rnaID,baitmap$bait)]
all_chr_freq$otherLen <- rmap$len[match(all_chr_freq$otherID,rmap$Frag_id)]
all_chr_freq$baitStart <-(baitmap$start[match(all_chr_freq$baitID,baitmap$Frag_id)]+baitmap$end[match(all_chr_freq$baitID,baitmap$Frag_id)])/2 - (rmap$start[match(all_chr_freq$otherID,rmap$Frag_id)] + rmap$end[match(all_chr_freq$otherID,rmap$Frag_id)])/2
all_chr_freq$distLen <- ifelse(all_chr_freq$rnachrom != all_chr_freq$dnachrom, NA, all_chr_freq$baitStart)
chinput <- data.frame(all_chr_freq$baitID,all_chr_freq$otherID,all_chr_freq$Freq,all_chr_freq$otherLen,all_chr_freq$distLen)
chinput <- chinput[with(chinput,order(as.numeric(chinput$all_chr_freq.baitID))),]
chinput <- chinput[!is.na(chinput$all_chr_freq.otherID),]
chinput <- chinput[!is.na(chinput$all_chr_freq.baitID),]
chinput[sapply(chinput, is.infinite )] <- NA
chinput$all_chr_freq.distLen <- as.numeric(chinput$all_chr_freq.distLen)
#Define the output file
write.table(chinput,opt$out,sep = "\t",col.names = FALSE,row.names = FALSE)




