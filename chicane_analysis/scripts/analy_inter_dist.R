res <- 
    fread("/rds/general/project/neurogenomics-lab/live/Projects/radicl_seq/data/chicane.results.rdcl.txt")

res <- res[!is.na(p.value)]
library(ggplot2)
library(cowplot)
library(viridis)
#make bins split res by q.value bins - go with 0.05 separate as rel high num and so you get q.val<0.05
#split rest by .2
res[,q_bin:=cut(q.value,breaks=c(0,0.05,0.2,0.4,0.6,0.8,1))]

#need to split out separate row for each count in count column to get density plot
#remove trans for space....
res_stretch_cis <- res[!is.na(distance),.(count2=1:count),names(res)]

ggplot(res_stretch_cis, aes(log(distance,base=10))) +
    geom_density(alpha = 0.1)+theme_cowplot()+
    scale_fill_viridis(name="Adj. P-value Bin",discrete=T)+
    scale_colour_viridis(name="Adj. P-value Bin",discrete=T)+
    geom_vline(xintercept=log(60000,base=10), colour="grey") +
    geom_vline(xintercept=log(300000,base=10), colour="grey") +
    annotate(x=log(60000,base=10),y=2,label="60Kb",vjust=2,geom="label")+
    annotate(x=log(300000,base=10),y=2,label="300Kb",vjust=2,geom="label")+
    xlab("Log 10 Distance from Interaction to RNA")

#split by groups
ggplot(res_stretch_cis, aes(log(distance,base=10), 
                                fill = q_bin, colour = q_bin)) +
    geom_density(alpha = 0.1)+theme_cowplot()+
    scale_fill_viridis(name="Adj. P-value Bin",discrete=T)+
    scale_colour_viridis(name="Adj. P-value Bin",discrete=T)+
    geom_vline(xintercept=log(60000,base=10), colour="grey") +
    geom_vline(xintercept=log(300000,base=10), colour="grey") +
    annotate(x=log(60000,base=10),y=2,label="60Kb",vjust=2,geom="label")+
    annotate(x=log(300000,base=10),y=2,label="300Kb",vjust=2,geom="label")+
    xlab("Log 10 Distance from Interaction to RNA")


#try standardising by transcript length?
res_stretch_cis[,bait.length:=bait.end-bait.start]
res_stretch_cis[,stnd_dist:=distance/bait.length]

ggplot(res_stretch_cis, aes(log(stnd_dist,base=10))) +
    geom_density(alpha = 0.1)+theme_cowplot()+
    scale_fill_viridis(name="Adj. P-value Bin",discrete=T)+
    scale_colour_viridis(name="Adj. P-value Bin",discrete=T)+
    geom_vline(xintercept=log(1,base=10), colour="grey") +
    geom_vline(xintercept=log(120,base=10), colour="grey") +
    annotate(x=log(3,base=10),y=.8,label="RNA Length",vjust=2,geom="label")+
    annotate(x=log(400,base=10),y=.8,label="120x RNA Length",vjust=2,
             geom="label")+
    xlab("Log 10 Distance from Interaction to RNA of RNA Length")


#split by groups
ggplot(res_stretch_cis, aes(log(stnd_dist,base=10), 
                            fill = q_bin, colour = q_bin)) +
    geom_density(alpha = 0.1)+theme_cowplot()+
    scale_fill_viridis(name="Adj. P-value Bin",discrete=T)+
    scale_colour_viridis(name="Adj. P-value Bin",discrete=T)+
    geom_vline(xintercept=log(1,base=10), colour="grey") +
    geom_vline(xintercept=log(120,base=10), colour="grey") +
    annotate(x=log(3,base=10),y=1.8,label="RNA Length",vjust=2,geom="label")+
    annotate(x=log(400,base=10),y=1.8,label="120x RNA Length",vjust=2,
                geom="label")+
    xlab("Log 10 Distance from Interaction to RNA of RNA Length")
#second peak is gone!!


#is there a relationship between RNA length and q-value?
cor_test <- cor.test(res_stretch_cis$bait.length,
                     res_stretch_cis$q.value,method = "spearman",exact = F)
#less sig than count for either p/q value so doesn't hold as much a relationship

#is there a relationship between distance and q-value?
cor_test <- cor.test(res_stretch_cis$distance,
                     res_stretch_cis$p.value,method = "spearman",exact = F)
#rho .46 - as distance inc so too does p-value
save(res_stretch_cis,
     file="/rds/general/project/neurogenomics-lab/live/Projects/radicl_seq/data/data_analy_inter_dist.RData")
