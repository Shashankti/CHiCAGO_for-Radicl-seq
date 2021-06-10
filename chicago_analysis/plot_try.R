#import libraries
library(doMC)
library(ggplot2)
registerDoMC(cores = 4)
library(data.table)
library(dplyr)
library(plyr)
#Read file
yyy = fread("/home/master/chicago-chinput-file/Day13_done.bedpe",header=FALSE,stringsAsFactors = FALSE,quote = "")
#colnames(yyy) <- c('rnachrom','rnachromStart','rnachromEnd','dnachrom','dnachromStart','dnachromEnd','cigar','ignore','rnaStrand','dnaStrand','rnQual','dnaQua','rnaID')
colnames(yyy) <- c('rnachrom','rnachromStart','rnachromEnd','dnachrom','dnachromStart','dnachromEnd','rnaID')

#Make length colums
yyy$len <- abs((yyy$dnachromStart+yyy$dnachromEnd)/2-(yyy$rnachromStart+yyy$rnachromEnd)/2)
mean(yyy$len)

#Remove trans inetractions
yyy$match <- NULL
yyy[yyy$rnachrom == yyy$dnachrom,"match"] <-"1"
yyy[yyy$rnachrom != yyy$dnachrom,"match"] <-"0"
vc <- '1'
schr <- yyy[is.element(yyy$match,vc),]
#########



ggplot(yyy,aes(x = len))+
  geom_histogram(aes(y = stat(count)/sum(count)))+
  scale_y_continuous(labels = scales::percent)+
  xlab("Distance of interaction")+
  ylab("Probability")+
  ggtitle("Probability of detecting an interaction based on distance from the RNA's genomic location")

#Sort the file
kls <- split(schr,schr$rnaID)
rm(yyy)
rm(schr)
gc()
kls$GNB1

######################### Get the plots
###################
setwd("~/Project/Plots/D13/")

plot_list = list()
for (i in 1:length(kls)){
  if (length(kls[[i]]$rnachrom)>100){
    p = ggplot(kls[[i]],aes(x = len))+
      geom_histogram(aes(y = stat(count)/sum(count)))+
      scale_y_continuous(labels = scales::percent)+
      ggtitle(kls[[i]]$rnaID[1])
    plot_list[[i]] = p
  } else {
    p = ggplot(kls[[i]],aes(x = len))+
      geom_freqpoly(aes(y = stat(count)/sum(count)))+
      scale_y_continuous(labels = scales::percent)+
      ggtitle(kls[[i]]$rnaID[1])
    plot_list[[i]] = p
  }
}
for (i in 1:length(kls)){
  file_name = paste("distance_plot_", kls[[i]]$rnaID[1], ".tiff", sep="")
  tiff(file_name)
  print(plot_list[[i]])
  dev.off()
}
pdf("plots.pdf")
for (i in 1:length(kls)){
  print(plot_list[[i]])
}
dev.off


#######
chr1 <- yyy[yyy$dnachrom == "chr1",]
chr1
ichr1 <- chr1[with(chr1,order(dnachromStart)),]
ichr1
size_chr1 <- 248956422
binchr1 <- transform(ichr1, group = cut(dnachromStart,
                                        breaks=seq(from = 0, to = 250000000, by = 20000 )))
binchr1
klbin <- split(binchr1,binchr1$group)
klbin[[12400]]

dt <- data.table(binchr1)
chr1_freq <- dt[,list(Freq = .N), by = list(rnaID,group)]
as.matrix(chr1_freq)
##############
p = ggplot(kls$PIWIL2,aes(x = len))+
  geom_histogram(aes(y = stat(count)/sum(count)))+
  scale_y_continuous(labels = scales::percent)+
  ggtitle(kls$PIWIL2$rnaID[1])+
  xlab("Distance of interaction")+
  ylab("Probability")
plot(p)
  