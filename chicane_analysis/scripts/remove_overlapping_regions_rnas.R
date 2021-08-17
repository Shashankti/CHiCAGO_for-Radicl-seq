#' Function to remove overlapping regions from RNA. Note this will remove from
#' the end of the RNA pair in question so the transcriptional start site of the
#' other RNA remains intact. Note this will also remove the RNA which are 
#' completely encompassed by another RNA. This function is useful for RADICL-Seq
#' data when using CHICANE to analysis results: RNA are the baits file.
#' 
#' @param rna_path file path to RNA txt file containing chr,start,end of RNA
#' @return Path to updated rna file
#' @importFrom data.table fread fwrite setnames setkey foverlaps setorder
remove_overlapping_regions_rnas <- function(rna_path){
    #read in RNA file
    baits <- data.table::fread(baits_rdcl)
    data.table::setnames(baits,c("chr","start","end","gene"))
    #create a unique ID field as gene name here isn't unique
    baits[,gene_id:=paste0(chr,":",start,"-",end)]
    data.table::setkey(baits,chr,start,end)
    
    #First - identify RNA which are entirely encompased by other RNA
    #these can not be dealt with in the model and should be removed
    within_ovrlps <- data.table::foverlaps(baits,baits,type="within",nomatch=0L)
    #remove overlaps with same region
    within_ovrlps <- within_ovrlps[gene_id!=i.gene_id,]
    # 20464 genes are encompassed by another gene and need to be removed
    encompassed_genes <- unique(within_ovrlps[,i.gene_id])
    baits <- baits[!(gene_id %in% encompassed_genes)]
    data.table::setkey(baits,chr,start,end)
    #now check for overlapping regions
    ovrlps <- data.table::foverlaps(baits,baits,nomatch=0L)
    #remove overlaps with same region
    ovrlps <- ovrlps[gene_id!=i.gene_id,]
    #remove overlapping part of end of RNA pairing
    # both a-b and b-a present, make sure only a-b kept where a.start<b.start
    ovrlps <- ovrlps[start<i.start]
    # a can be overlapping two other rna's take the one with the earlier start site
    # this will give the shortest distance of a
    data.table::setorder(ovrlps,chr,start,end,i.start)
    ovrlps <- ovrlps[!duplicated(ovrlps$gene_id),]
    #7.3k gene pairs have overlap to be removed
    ovrlps[,new_end:=i.start-1]
    data.table::setkey(ovrlps,chr,start,end)
    #Now update baits file
    baits[ovrlps,new_end:=i.new_end]
    #Now update end if needed
    baits[!is.na(new_end),end:=new_end]
    #clean up file
    baits[,new_end:=NULL]
    baits[,gene_id:=NULL]
    #save to path with updated name
    #remove extension
    ext<- sub('.*\\.', '', rna_path)
    rna_path_no_ovrlps <-
        sub(pattern = "(.*)\\..*$", replacement = "\\1", rna_path)
    rna_path_no_ovrlps <-paste0(rna_path_no_ovrlps,"_no_ovrlps.",ext)
    data.table::fwrite(baits,file=rna_path_no_ovrlps)
    return(rna_path_no_ovrlps)
}

baits_rdcl <- "~/Documents/ICL/radicl_seq_analysis/all_genes_location.txt"
remove_overlapping_regions_rnas(baits_rdcl)
