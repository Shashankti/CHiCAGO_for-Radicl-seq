BiocManager::install("genomation")

#read in chicane results
res <- 
    fread("/rds/general/project/neurogenomics-lab/live/Projects/radicl_seq/data/chicane.results.rdcl.txt")

#E025 Adipose Derived Mesenchymal Stem Cell Cultured Cells - hg38
#from: https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/
#This bed file only contains one metadata column so genomation::readBed will fail
#Approach below uses code from in this function
chrHMM_url <- 
    "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/E025_15_coreMarks_hg38lift_segments.bed.gz"
file <- genomation:::compressedAndUrl2temp(chrHMM_url)
#read in as dataframe
df<-readr::read_delim(file,skip=0,col_names=FALSE, delim="\t")
names(df) <- c("chr","start","end","ChromHMM") 
#make GRanges obj
g = GenomicRanges::makeGRangesFromDataFrame(
    df, 
    keep.extra.columns=FALSE, 
    starts.in.df.are.0based=FALSE,
    ignore.strand=is.null(NULL))
col.names1=list(chr=1,start=2,end=3,strand=NULL)
black.names=c("seqnames", "ranges", "strand", "seqlevels", "seqlengths",
              "isCircular", "start", "end", "width", "element")
my.mcols = df[,-unlist(col.names1),drop=FALSE]
#add meta data column
GenomicRanges::mcols(g) = my.mcols[, !colnames(my.mcols) %in% black.names]
#split GRange object by the different names
chrHMM_list <- GenomicRanges::split(g, g$ChromHMM, drop = TRUE)
#get state names for IDs
sort(unique(df$ChromHMM))
#annotations found at:
# https://www.nature.com/articles/nature14309
#1_TssA	active transcription start site (TSS)  -proximal promoter states 
#2_TssAFlnk Flanking active transcription start site (TSS). -proximal promoter states 
#3_TxFlnk a transcribed state at the 5′ and 3′ end of genes -showing both promoter and enhancer signatures
#4_Tx Strong transcription
#5_TxWk Weak transcription
#6_EnhG Genic enhancer -genic enhancers occur more frequently in gene bodies & in exons, while enhancers don't
#7_Enh enhancers
#8_ZNF/Rpts zinc finger protein genes and repeats. -ZNFs are involved in the regulation of several cellular processes (v abundant)
#9_Het heterochromatin. - tightly packed DNA, inactive states 
#10_TssBiv bivalent/poised transcription start site. -both activation and repression (low expression), poised to express
#11_BivFlnk flanking bivalent transcription start site/enhancer
#12_EnhBiv bivalent enhancer
#13_ReprPC repressed Polycomb. -gene silencing
#14_ReprPCWk weak repressed Polycomb
#15_Quies quiescent/low - inactivity/dormitary
rm(g,df)

#make results into GRanges
#need to specify the DNA segments as the start and end columns
setnames(res,c("target.chr","target.start","target.end"),
            c("chr","start","end"))
#first NA's in q-values, 50k of them
res <- res[!is.na(q.value),]
#make bins split res by q.value bins - go with 0.05 separate as rel high num and so you get q.val<0.05
#split rest by .2
res[,q_bin:=cut(q.value,breaks=c(0,0.05,0.2,0.4,0.6,0.8,1))]
#check counts of bins to ensure enough in each
table(res$q_bin) #fine..
#make GRanges
gres = GenomicRanges::makeGRangesFromDataFrame(
    res, 
    keep.extra.columns=TRUE, 
    starts.in.df.are.0based=FALSE,
    ignore.strand=TRUE)
#now split based on q-value bins
gres_list <- GenomicRanges::split(gres, gres$q_bin, drop = TRUE)
rm(gres)
#now get overlap with ChromHMM annotation data - target is radicl-seq data
annotation <- 
    genomation::annotateWithFeatures(gres_list, chrHMM_list)
#interested in percentage of target elements overlapping with features:
#make a heatmap
#don't produce heatmap, take matrix, standardise and replot
annot_mtrx <- genomation::heatTargetAnnotation(annotation,plot=FALSE)
#change column order
col.order <- paste0("E",1:15)
annot_mtrx <- annot_mtrx[,col.order]
#plot raw first
#convert to df
annot_df<-reshape2::melt(annot_mtrx)
names(annot_df) <- c("q-value","State","Overlap")
#Add better explaining state names
library(dplyr)
annot_df <- annot_df%>% 
    mutate(State_name = as.factor(case_when(
        .$State == "E1" ~ "1.Active.TSS",
        .$State == "E2" ~ "2.Flanking.Active.TSS",
        .$State == "E3" ~ "3.Flanking.Transcribed.State",
        .$State == "E4" ~ "4.Strong.Transcription",
        .$State == "E5" ~ "5.Weak.Transcription",
        .$State == "E6" ~ "6.Genic.Enhancer",
        .$State == "E7" ~ "7.Enhancer",
        .$State == "E8" ~ "8.ZNC.Finger.Prot.&.Repeats",
        .$State == "E9" ~ "9.Heterochromatin",
        .$State == "E10" ~ "10.Bivalent/Poised.TSS",
        .$State == "E11" ~ "11.Flanking.Bivalent.TSS",
        .$State == "E12" ~ "12.Bivalent.Enhancer",
        .$State == "E13" ~ "13.Repressed.Polycomb",
        .$State == "E14" ~ "14.Weak.Repressed.Polycomb",
        .$State == "E15" ~ "15.Quiescent")))
library(ggplot2)
library(cowplot)
library(viridis)
# Reorder factor levels for plot
annot_df$State_name <- 
    factor(annot_df$State_name,
            levels = levels(annot_df$State_name)[c(1,8:15,2:7)])
ggplot(annot_df,aes(x=State_name,y=`q-value`,fill=Overlap))+
    geom_tile()+theme_cowplot()+scale_fill_viridis_c()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("ChromHMM State")+
    ylab("Adj. P-value Bin")
#very little notable difference in each state for diff q-values
#try standardising across states
scl_annot_mtrx <- apply(annot_mtrx,2,function(x) scale(x)[,1])
#plot
scl_annot_df<-reshape2::melt(scl_annot_mtrx)
names(scl_annot_df) <- c("q-value","State","Overlap")
#Add better explaining state names
library(dplyr)
scl_annot_df <- scl_annot_df%>% 
    mutate(State_name = as.factor(case_when(
        .$State == "E1" ~ "1.Active.TSS",
        .$State == "E2" ~ "2.Flanking.Active.TSS",
        .$State == "E3" ~ "3.Flanking.Transcribed.State",
        .$State == "E4" ~ "4.Strong.Transcription",
        .$State == "E5" ~ "5.Weak.Transcription",
        .$State == "E6" ~ "6.Genic.Enhancer",
        .$State == "E7" ~ "7.Enhancer",
        .$State == "E8" ~ "8.ZNC.Finger.Prot.&.Repeats",
        .$State == "E9" ~ "9.Heterochromatin",
        .$State == "E10" ~ "10.Bivalent/Poised.TSS",
        .$State == "E11" ~ "11.Flanking.Bivalent.TSS",
        .$State == "E12" ~ "12.Bivalent.Enhancer",
        .$State == "E13" ~ "13.Repressed.Polycomb",
        .$State == "E14" ~ "14.Weak.Repressed.Polycomb",
        .$State == "E15" ~ "15.Quiescent")))
# Reorder factor levels for plot
scl_annot_df$State_name <- 
    factor(scl_annot_df$State_name,
           levels = levels(scl_annot_df$State_name)[c(1,8:15,2:7)])
ggplot(scl_annot_df,aes(x=State_name,y=`q-value`,fill=Overlap))+
    geom_tile()+theme_cowplot()+scale_fill_viridis_c()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("ChromHMM State")+
    ylab("Adj. P-value Bin")
save(annotation,annot_df,scl_annot_df,file=
     "/rds/general/project/neurogenomics-lab/live/Projects/radicl_seq/data/ChromHMM_state_results.RData")
saveRDS(annotation,"/rds/general/project/neurogenomics-lab/live/Projects/radicl_seq/data/ChromHMM_state_annotation.rds")
