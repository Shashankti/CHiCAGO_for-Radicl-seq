#Analysis of ibed file
#Read ibed
library(optparse)
library(gridExtra)
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
myCluster <- makeCluster(3, # number of cores to use
                         type = "FORK") # type of cluster
registerDoParallel(myCluster)
ibed <- fread("Data/Processed/InputFile/vignetteOutput.ibed",header=TRUE,stringsAsFactors = FALSE,quote = "")
ibed$dist <- ifelse(ibed$bait_chr != ibed$otherEnd_chr,NA,abs((ibed$bait_start+ibed$bait_end)/2-(ibed$otherEnd_start+ibed$otherEnd_end)/2))
#Read Example Dataset
ex <- fread("Data/Processed/InputFile/exampleOutput.ibed",header=TRUE,stringsAsFactors = FALSE,quote = "")
ex$dist <- ifelse(ex$bait_chr != ex$otherEnd_chr,NA,abs((ex$bait_start+ex$bait_end)/2-(ex$otherEnd_start+ex$otherEnd_end)/2))
#Read size file
chrsize <- fread("/home/master/chicago-chinput-file/D13.sizes",header=FALSE,stringsAsFactors = FALSE,quote = "")
chrsize <- chrsize[order(as.numeric(as.character(chrsize$V1)))]

#Function to determine what percentage of genome is covered in significant interactions
Percentage_significant <- function(data,threshold){
  #Significant Interactions
  signif <- data[data$score > threshold] 
  #Find what percentage of genome is covered in significant interactions
  nif1 <- signif[,c(5,6,7,10)]
  type <- makeGRangesFromDataFrame(nif1)
  type1 <- data.frame(unique(type))
  #Percent of significant intercations is given by
  per <- 100 *sum(type1$width)/sum(chrsize$V2[1:22])
  print(per)
}
#Find how far are they generally from transcription start sites.
#Make plot for one bait.
jist <- split(signif,signif$bait_name)
#Considering only non NA interactions
nona <- ibed[ibed$score > 5]
nona <- nona[!is.na(nona$dist)]
up <-  split(nona,nona$bait_name)

p = ggplot(up[[20]],aes(x = dist,y=N_reads,color=score))+
  geom_point()
plot(p)
#List of plots


plot_list = list()
for (i in 1:length(up)){
  p = ggplot(up[[i]],aes(x=dist,y=score,color=N_reads))+
    ggtitle(up[[i]]$bait_name[1])+
    geom_point()
  plot_list[[i]]=p
}
for (i in 1:length(up)){
  file_name = paste("distance_plot_", up[[i]]$bait_name[1], ".tiff", sep="")
  tiff(file_name)
  print(plot_list[[i]])
  dev.off()
}
pdf("plots.pdf")
for (i in 1:length(unq)){
  print(sig_list[[i]])
}
dev.off
##

ggplot(cd@x, aes(x = distSign, y = score)) +
  geom_point(aes(color = ifelse(N>70, 'red', 'black'))) +
  scale_colour_manual(labels = c("<70", ">70"), values=c('black', 'red')) + 
  ylim(0,10000)+
  labs(color = "Values")



pdf("Plots/plots.pdf", onefile = TRUE)
for (i in seq(length(sig_list))) {
  do.call("grid.arrange", sig_list[[i]])  
}
dev.off
#Plot which dont follow distance relation
sig_list = list()
#Sort data from nona
like <- nona[nona$dist>2e+07 & nona$score>65]
unq <- unique(like$bait_name)
for (i in 1:length(unq)){
  p = ggplot(nona[nona$bait_name==unq[[i]]],aes(x=dist,y=score,color=N_reads))+
    ggtitle(nona[nona$bait_name==unq[[i]]]$bait_name[1])+
    geom_point()
 sig_list[[i]]=p
}



#Plots from example data








#Plot number of interactions based on the distance and significant interactions
q = ggplot(ibed[ibed$bait_name=="ABI2"],aes(x=dist,color=score))+
  geom_point(aes(y=N_reads))
plot(q)

