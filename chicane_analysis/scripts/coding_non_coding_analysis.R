res <- 
    fread("/rds/general/project/neurogenomics-lab/live/Projects/radicl_seq/data/chicane.results.rdcl.txt")

res <- res[!is.na(p.value)]
library(ggplot2)
library(cowplot)
library(viridis)
#make bins split res by q.value bins - go with 0.05 separate as rel high num and so you get q.val<0.05
#split rest by .2
res[,q_bin:=cut(q.value,breaks=c(0,0.05,0.2,0.4,0.6,0.8,1))]

#add cis, trans indicator
res[is.na(distance),interaction:="trans"]
res[!is.na(distance),interaction:="cis"]

#add in gene type i.e. coding, non-coding
#get gene names added on to results
baits <- 
    fread("/rds/general/project/neurogenomics-lab/live/Projects/radicl_seq/data/all_genes_location.txt",
            header = FALSE,stringsAsFactors = FALSE)
colnames(baits) <- c("chr","start","end","bait_name")    
res[,bait_name := baits$bait_name[match(res$bait.start,baits$start)]]

genes <- unique(res$bait_name)
library(biomaRt)
mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "useast.ensembl.org")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
    mart = mart,
    attributes = c(
        "hgnc_symbol",
        "entrezgene_id",
        "ensembl_gene_id",
        "gene_biotype"),
    filter = "hgnc_symbol",
    values = genes,
    uniqueRows=TRUE,
    useCache = FALSE)

annotLookup <- as.data.table(annotLookup)

#gene_biotype
setnames(annotLookup,
            c("bait_name","entrezgene_id","ensembl_gene_id","gene_biotype"))
setkey(res,bait_name)
setkey(annotLookup,bait_name)
res[annotLookup,gene_biotype:= i.gene_biotype]
#clear memory
rm(genes,annotLookup,mart,baits)
#10,004,838 with 1,716,157 without
res_biotyp <- res[!is.na(gene_biotype),]

#info here: http://www.ensembl.org/info/genome/genebuild/biotypes.html
#
#Pseudogene: A gene that has homology to known protein-coding genes but contain 
#a frameshift and/or stop codon(s) which disrupts the ORF. Thought to have 
#arisen through duplication followed by loss of function
#
#Protein coding: Gene/transcipt that contains an open reading frame (ORF).
#
#TR gene: T cell receptor gene that undergoes somatic recombination, annotated 
#in collaboration with IMGT http://www.imgt.org/
#
#TEC (To be Experimentally Confirmed): Regions with EST clusters that have polyA
#features that could indicate the presence of protein coding genes. These 
#require experimental validation, either by 5' RACE or RT-PCR to extend the 
#transcripts, or by confirming expression of the putatively-encoded peptide with
#specific antibodies.
#
#ncRNA: A non-coding gene + Long non-coding RNA (lncRNA) A non-coding 
#gene/transcript >200bp in length
#
#Mt_RNA: mitochondrial RNA
#
#Ribozyme: ribozymal RNA
#
#DEF BELOW FROM WIKAPEDIA**
#scRNA: A small conditional RNA (scRNA) is a small RNA molecule or complex 
#(typically less than approximately 100 nt) engineered to interact and change 
#conformation conditionally in response to cognate molecular inputs so as to 
#perform signal transduction in vitro, in situ, or in vivo.
#
#IG gene: Immunoglobulin gene that undergoes somatic recombination, annotated in
#collaboration with IMGT http://www.imgt.org/.
#
res_biotyp[grepl( "protein_coding", gene_biotype),gene_biotype_hl:=gene_biotype]
res_biotyp[grepl( "pseudogene", gene_biotype),gene_biotype_hl:="pseudogene"]
res_biotyp[grepl( "TR_", gene_biotype),gene_biotype_hl:="TR_gene"]
res_biotyp[grepl( "TEC", gene_biotype),gene_biotype_hl:=gene_biotype]
res_biotyp[gene_biotype %in% c("miRNA","miscRNA","piRNA","rRNA","siRNA","snRNA",
                                "snoRNA","tRNA","vaultRNA","lncRNA","misc_RNA",
                                "scaRNA"),
           gene_biotype_hl:="ncRNA"]
res_biotyp[grepl( "Mt_", gene_biotype),gene_biotype_hl:="Mt_RNA"]
res_biotyp[gene_biotype=="ribozyme",gene_biotype_hl:=gene_biotype]
res_biotyp[gene_biotype=="scRNA",gene_biotype_hl:=gene_biotype]
res_biotyp[grepl( "IG_", gene_biotype),gene_biotype_hl:="IG_gene"]

#plot sig with higher and lower level RNA types
ggplot(res_biotyp[q.value<0.05], aes(x=gene_biotype_hl,fill=gene_biotype_hl)) +
    geom_bar(aes(y = (..count..)/sum(..count..)))+
    facet_wrap(~interaction)+
    theme_cowplot()+
    scale_fill_viridis(discrete=T)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position="none")+
    ylab("")+
    xlab("Proportion of significant cis/trans interactions")

ggplot(res_biotyp[q.value<0.05], aes(x=gene_biotype,fill=gene_biotype)) +
    geom_bar(aes(y = (..count..)/sum(..count..)))+
    facet_wrap(~interaction)+
    theme_cowplot()+
    scale_fill_viridis(discrete=T)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position="none")+
    xlab("")+
    ylab("Proportion of significant cis/trans interactions")

#finally for cis what is the difference in distance of interaction across
#diff RNA types
#Remove IG TR, only 2,3 interactions respectively
ggplot(res_biotyp[q.value<0.05 & interaction=="cis" & 
                      (!gene_biotype_hl %in% c("IG_gene","TR_gene"))], 
       aes(x=gene_biotype_hl,y=log(distance,base=10),fill=gene_biotype_hl)) +
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.1)+
    theme_cowplot()+
    scale_fill_viridis(discrete=T,alpha=0.5)+
    theme(legend.position="none")+
    ylab("Log 10 Distance from Interaction to RNA")+
    xlab("Interaction Type")


#snoRNA,misc_RNA removed too low cases - 2,3
#removed pseudogene RNAs too to compare coding and non-coding directly
ggplot(res_biotyp[q.value<0.05 & interaction=="cis" & 
                      (!gene_biotype_hl %in% c("IG_gene","TR_gene",
                                                "pseudogene")) &
                      (!gene_biotype %in% 
                           c("snoRNA","misc_RNA"))], 
       aes(x=gene_biotype,y=log(distance,base=10),fill=gene_biotype)) +
    geom_violin(trim=FALSE)+
    geom_boxplot(width=0.1)+
    theme_cowplot()+
    scale_fill_viridis(discrete=T,alpha=0.5)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position="none")+
    ylab("Log 10 Distance from Interaction to RNA")+
    xlab("Interaction Type")

#Next I'd like to annotate the DNA these RNA Bind to (similar to ChromHMM but
#giving......
