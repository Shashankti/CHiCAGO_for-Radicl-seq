library(optparse)
library(doMC)
library(ggplot2)

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

args<-commandArgs(TRUE)

#Reading the transript file
yyy = fread(args[1],header=FALSE,stringsAsFactors = FALSE,quote = "")
#colnames(yyy) <- c('rnachrom','rnachromStart','rnachromEnd','dnachrom','dnachromStart','dnachromEnd','cigar','ignore','rnaStrand','dnaStrand','rnQual','dnaQua','rnaID')
colnames(yyy) <- c('rnachrom','rnachromStart','rnachromEnd','dnachrom','dnachromStart','dnachromEnd','rnaID')
unqID = unique(yyy$rnaID)
IDs = fread(args[2],header=FALSE,stringsAsFactors = FALSE,quote = "")
#Read chrsize file
chrsize <- fread(args[3],header=FALSE,stringsAsFactors = FALSE,quote = "")

#### Generating trial baitmap file


data = filter(IDs,V4 %in% unqID)
data = data %>% mutate_all(~gsub("chr","",.))

data[,sapply(V1,is.numeric)]
data = data[!data$V1 %in% c("X","Y"),]

data <- data[with(data,order(as.numeric(V1,V2))),]


#Naming columns in baitmap
colnames(data) <- c('chr','start','end','baitAnnotation')
######
data$chr <- as.numeric(data$chr)
data$start <- as.numeric(data$start)
data$end <- as.numeric(data$end)


###
#Create rmap

bait_try <- data.frame(data[,1:3])
set.seed(2323)
windowSize = 20000
map = NULL
some = NULL
l = NULL
chkl <- split(yyy,yyy$dnachrom)
rm(yyy)
gc()


lit1 = NULL
lit2 = NULL
lit3 = NULL
for(i in 1:22){
  lit1[[i]] = sample.int(chrsize$V2[i],length(rsplit[[i]]$chrom))
  lit2[[i]] <- append(lit1[[i]],chrsize$V2[i])
  lit1[[i]] <- append(lit1[[i]],0)
  lit1[[i]] <- lit1[[i]][order(lit1[[i]])]
  lit2[[i]] <- lit2[[i]][order(lit2[[i]])]
}

for (i in 1:22) {
  lit3[[i]] <- rep(i,length(lit1[[i]]))
  rrmap[[i]] = data.frame("chr" = lit3[[i]],"start" = lit1[[i]], "end" = lit2[[i]])
}
df <- bind_rows(rrmap, .id = "column_label")
df$column_label <- NULL
colnames(df)[1] <- 'chr'




write.table(data,args[4],sep = "\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
write.table(df,args[5],sep = "\t",col.names = FALSE,row.names = FALSE)
