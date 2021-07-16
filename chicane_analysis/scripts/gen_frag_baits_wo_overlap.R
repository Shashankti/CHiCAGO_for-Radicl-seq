#make fragment and baits file

library(data.table)
library(GenomicRanges)
options(scipen = 999)


#read initial file
#Reading the transript file
yyy = fread("Data/Raw/Day13_done.bedpe",header=FALSE,stringsAsFactors = FALSE,quote = "")
#colnames(yyy) <- c('rnachrom','rnachromStart','rnachromEnd','dnachrom','dnachromStart','dnachromEnd','cigar','ignore','rnaStrand','dnaStrand','rnQual','dnaQua','rnaID')
colnames(yyy) <- c('rnachrom','rnachromStart','rnachromEnd','dnachrom','dnachromStart','dnachromEnd','rnaID')
unqID = unique(yyy$rnaID)
IDs = fread("Data/Raw/all_genes_location.txt",header=FALSE,stringsAsFactors = FALSE,quote = "")
#Read chrsize file
chrsize <- fread("Data/Raw/D13.sizes",header=FALSE,stringsAsFactors = FALSE,quote = "")

#### Generating trial baitmap file by intersecting the bait names with the raw data
data = IDs[IDs$V4 %in% unqID]

data = data[!data$V1 %in% c("chrX","chrY"),]

data <- data[with(data,order(as.numeric(V1,V2))),]


#Naming columns in baitmap
colnames(data) <- c('chr','start','end','baitAnnotation')
data$start <- as.numeric(data$start)
data$end <- as.numeric(data$end)
# Removing overlap from baits file
reduced_data <- as.data.frame(reduce(GRanges(data$schr,IRanges(data$start,data$end))))



#Generating frags file with the raw data.
frag <- yyy[,4:6]
frag <- frag[!duplicated(frag)]
colnames(frag) <- c("chrom","start","end")
#reducing the frags to remove overlaps
yelp <- as.data.frame(reduce(GRanges(frag$chrom,IRanges(frag$start,frag$end))))
rmap <- yelp[,1:3]
rmap$seqnames <- paste0("chr",rmap$seqnames)
#convert baits to genomic ranges
bres <- makeGRangesFromDataFrame(data)
res <- makeGRangesFromDataFrame(rmap)
#finding overlaps and creating the set difference
ov <- findOverlaps(res,bres)
grl <- extractList(bres,as(ov,"List"))
blip <- psetdiff(res,grl)
blip <- as.data.frame(blip)
blip$group = NULL
blip$group_name = NULL
bres_df <- as.data.frame(bres)
rmap <- rbind(blip,bres_df)
rmap <- rmap[order(rmap$seqnames,rmap$start),]
colnames(rmap)[1] <- "chrom"
rmap$strand = NULL
rmap$width = NULL
rmap <-rmap[with(rmap,order(chrom,start)),] 


fwrite(rmap,file="Data/Processed/graphs/frags_wo_ovrlap.bed",col.names = FALSE,row.names = FALSE,sep = "\t")
baits <- data[,1:3]
fwrite(baits,file = "Data/Processed/graphs/chic_wo_baits.bed",sep="\t",col.names = FALSE,row.names = FALSE)

#Using reduced baits

reduced_bait <- makeGRangesFromDataFrame(reduced_data)
#finding overlaps and creating the set difference
ov <- findOverlaps(res,reduced_bait)
grl <- extractList(reduced_bait,as(ov,"List"))
blip <- psetdiff(res,grl)
blip <- as.data.frame(blip)
blip$group = NULL
blip$group_name = NULL
bres_df <- as.data.frame(reduced_bait)
rmap <- rbind(blip,bres_df)
rmap <- rmap[order(rmap$seqnames,rmap$start),]
colnames(rmap)[1] <- "chrom"
rmap$strand = NULL
rmap$width = NULL
rmap <-rmap[with(rmap,order(chrom,start)),] 
fwrite(rmap,file = "Data/Processed/graphs/frags_wo_bait_overlap.bed",col.names = FALSE,row.names = FALSE,sep = "\t")
baits <- reduced_data[,1:3]
fwrite(baits,file = "Data/Processed/graphs/wo_overlap_baits.bed",sep="\t",col.names = FALSE,row.names = FALSE)

#Using 2kb regions without removing bait overlap
binned_frag <- fread(file)
binned_res <- makeGRangesFromDataFrame(binned_frag)
#finding overlaps and creating the set difference
ov <- findOverlaps(binned_res,bres)
grl <- extractList(bres,as(ov,"List"))
blip <- psetdiff(binned_res,grl)
blip <- as.data.frame(blip)
blip$group = NULL
blip$group_name = NULL
bres_df <- as.data.frame(bres)
rmap <- rbind(blip,bres_df)
rmap <- rmap[order(rmap$seqnames,rmap$start),]
colnames(rmap)[1] <- "chrom"
rmap$strand = NULL
rmap$width = NULL
rmap <-rmap[with(rmap,order(chrom,start)),] 
fwrite(rmap,file = "Data/Processed/graphs/frags_binned_with_bait_overlap.bed",col.names = FALSE,row.names = FALSE,sep = "\t")

#Without bait overlap

#finding overlaps and creating the set difference
ov <- findOverlaps(binned_res,reduced_bait)
grl <- extractList(reduced_bait,as(ov,"List"))
blip <- psetdiff(binned_res,grl)
blip <- as.data.frame(blip)
blip$group = NULL
blip$group_name = NULL
bres_df <- as.data.frame(reduced_bait)
rmap <- rbind(blip,bres_df)
rmap <- rmap[order(rmap$seqnames,rmap$start),]
colnames(rmap)[1] <- "chrom"
rmap$strand = NULL
rmap$width = NULL
rmap <-rmap[with(rmap,order(chrom,start)),] 
fwrite(rmap,file = "Data/Processed/graphs/frags_binned_wo_bait_overlap.bed",col.names = FALSE,row.names = FALSE,sep = "\t")



