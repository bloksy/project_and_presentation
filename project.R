### Packages ####
library(dplyr)
library(tidyverse)
library(stringr)
#source("https://bioconductor.org/biocLite.R")
#biocLite()
require(DESeq2) #analysis itself
library(data.table)
library(ggplot2)
library(gridExtra)


##### Data #####
setwd("D:/GoogleDrive/Labor/DR/Andmed/DESeq2")


mirna_hor=read.table("mirna_hor.txt", sep="\t", header=T)
mirna_ver<-read.table("mirna_ver.txt",header=T)

# ### #mirna-mirna correlation analysis - mirmir ####
# 
# if(dir.exists("mir_mir")){
#   dir.create("mir_mir")
# }
# 
# 
# setwd("D:/GoogleDrive/Labor/DR/Andmed/DESeq2/mir_mir")
# 
# 
# colNR=2
# for(i in row.names(mirna_hor)){
#   for(j in row.names(mirna_ver)){
#     colNR=colNR+1
#     b<-colnames(mirna_hor[colNR])
#     testTable=data.frame(mirna_hor[,colNR],mirna_hor$type)
#     a<-str_replace(j,"-",".")
#     colnames(testTable)=c("mirna","type")
#     dds=DESeqDataSetFromMatrix(countData=mirna_ver, colData=testTable, design=~type+mirna)
#     dds=DESeq(dds,test="LRT",reduced=~1)
#     ddsReplace <- replaceOutliers(dds)
#     ddsReplace <- DESeq(ddsReplace)
#     res=results(ddsReplace)
#     write.table(res, file=paste0("mirna_",b,"_test_DESeq2.txt"))
#     print(b)
#   }
# }
# 
# 
# ###tulbad järjest#######
# 
setwd("D:/GoogleDrive/Labor/DR/Andmed/DESeq2/mir_mir")
# 
# failid=list.files()
# koik_andmed=read.table(fail, header=T) ###algse matrixi suuruse esitamine
# for (fail in failid) {
#   id=unlist(strsplit(fail, "[_]"))[2]
#   print(id)
#   andmed=read.table(fail, header=T)
#   colnames(andmed)[colnames(andmed)=="padj"]<-id
#   koik_andmed=cbind(koik_andmed, andmed[6])
# }
# 
# write.table(koik_andmed, "mir_mir_koos.tsv", sep="\t", row.names=T, quote=F)
# 
mir_mir <- read.table("mir_mir_koos.tsv", head=T)
mir_mir <- mir_mir[!(names(mir_mir) %in% c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))]
mir_mir <- setDT(mir_mir, keep.rownames = T)
mir_mir[1:10,1:7]


### parema p-valuega mirna-mirna####
colSums(mir_mir < 0.00005)

### plot top 10 miR.30d mirnat
head(mir_mir$miR.003)

mir_mir<-arrange(mir_mir, miR.003)

mir_mir[1:10,1]



####plot mirna-mirna#####

top_mirna <- mirna_hor[mir_mir[1,1]]


mir_mir[1:10,1]
top10_2 <- mirna_hor[,mir_mir[2,1]]
top10_3 <- mirna_hor[,mir_mir[3,1]]
top10_4 <- mirna_hor[,mir_mir[4,1]]
top10_5 <- mirna_hor[,mir_mir[5,1]]
top10_6 <- mirna_hor[,mir_mir[6,1]]
top10_7 <- mirna_hor[,mir_mir[7,1]]
top10_8 <- mirna_hor[,mir_mir[8,1]]
top10_9 <- mirna_hor[,mir_mir[9,1]]
top10_10 <- mirna_hor[,mir_mir[10,1]]


p1 <- mirna_hor %>% ggplot(aes(top_mirna, top10_2 	))+geom_point()+geom_smooth(method="lm")
p2 <- mirna_hor %>% ggplot(aes(top_mirna, top10_3 	))+geom_point()+geom_smooth(method="lm")
p3 <- mirna_hor %>% ggplot(aes(top_mirna, top10_4 	))+geom_point()+geom_smooth(method="lm")
p4 <- mirna_hor %>% ggplot(aes(top_mirna, top10_5 	))+geom_point()+geom_smooth(method="lm")
p5 <- mirna_hor %>% ggplot(aes(top_mirna, top10_6 	))+geom_point()+geom_smooth(method="lm")
p6 <- mirna_hor %>% ggplot(aes(top_mirna, top10_7 	))+geom_point()+geom_smooth(method="lm")
p7 <- mirna_hor %>% ggplot(aes(top_mirna, top10_8 	))+geom_point()+geom_smooth(method="lm")
p8 <- mirna_hor %>% ggplot(aes(top_mirna, top10_9 	))+geom_point()+geom_smooth(method="lm")
p9 <- mirna_hor %>% ggplot(aes(top_mirna, top10_10 	))+geom_point()+geom_smooth(method="lm")


grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9,ncol=3)

###mirna-gene correlation analysis ####

####DATA####
# rm(list = setdiff(ls(), lsf.str()))

setwd("D:/GoogleDrive/Labor/DR/Andmed/DESeq2")
mirna_hor=read.table("mirna_hor.txt", sep="\t", header=T)
gene_ver<-read.table("gene_ver.txt", header=T)



# if(dir.exists("mir_gene")){
#   dir.create("mir_gene")
# }

setwd("D:/GoogleDrive/Labor/DR/Andmed/DESeq2/mir_gene")

# 
# 
# 
# colNR=2
# 
# for(i in row.names(mirna_hor)){
#   for(j in row.names(gene_ver)){
#     colNR=colNR+1
#     b<-colnames(mirna_hor[colNR])
#     testTable=data.frame(mirna_hor[,colNR],mirna_hor$type)
#     a<-str_replace(j,"-",".")
#     colnames(testTable)=c("gene","type")
#     dds=DESeqDataSetFromMatrix(countData=gene_ver, colData=testTable, design=~type+gene)
#     dds=DESeq(dds,test="LRT",reduced=~1)
#     ddsReplace <- replaceOutliers(dds)
#     ddsReplace <- DESeq(ddsReplace)
#     res=results(ddsReplace)
#     write.table(res, file=paste0("mirna_",b,"_test_DESeq2.txt"))
#     print(b)
#   }
# }


####tulbad järjest#######
# 
# failid=list.files()
# koik_andmed=read.table(fail, header=T) ###algse matrixi suuruse esitamine
# for (fail in failid) {
#   id=unlist(strsplit(fail, "[_]"))[2]
#   print(id)
#   andmed=read.table(fail, header=T)
#   colnames(andmed)[colnames(andmed)=="padj"]<-id
#   koik_andmed=cbind(koik_andmed, andmed[6])
# }
# 
# write.table(koik_andmed, "mir_gene_koos.tsv", sep="\t", row.names=T, quote=F)

mir_gene <- read.table("mir_gene_koos.tsv", head=T, sep="\t")
mir_gene=mir_gene[!(names(mir_gene) %in% c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))]
mir_gene <- setDT(mir_gene, keep.rownames = T)
rn<-mir_gene$rn

mir_gene[1:10,1:7]

#
#
#### parema p-valuega mirna-gene####
colSums(mir_gene < 0.00005)

### plot####
mir_gene<-arrange(mir_gene, miR.003)

mir_gene[1:10,1]


gene_hor = t(gene_ver)


top10_2 <- gene_hor[,mir_gene[1,1]]
top10_3 <- gene_hor[,mir_gene[2,1]]
top10_4 <- gene_hor[,mir_gene[3,1]]
top10_5 <- gene_hor[,mir_gene[4,1]]
top10_6 <- gene_hor[,mir_gene[5,1]]
top10_7 <- gene_hor[,mir_gene[6,1]]
top10_8 <- gene_hor[,mir_gene[7,1]]
top10_9 <- gene_hor[,mir_gene[8,1]]
top10_10 <- gene_hor[,mir_gene[9,1]]


p1 <- mirna_hor %>% ggplot(aes(top_mirna, top10_2 	))+geom_point()+geom_smooth(method="lm")
p2 <- mirna_hor %>% ggplot(aes(top_mirna, top10_3 	))+geom_point()+geom_smooth(method="lm")
p3 <- mirna_hor %>% ggplot(aes(top_mirna, top10_4 	))+geom_point()+geom_smooth(method="lm")
p4 <- mirna_hor %>% ggplot(aes(top_mirna, top10_5 	))+geom_point()+geom_smooth(method="lm")
p5 <- mirna_hor %>% ggplot(aes(top_mirna, top10_6 	))+geom_point()+geom_smooth(method="lm")
p6 <- mirna_hor %>% ggplot(aes(top_mirna, top10_7 	))+geom_point()+geom_smooth(method="lm")
p7 <- mirna_hor %>% ggplot(aes(top_mirna, top10_8 	))+geom_point()+geom_smooth(method="lm")
p8 <- mirna_hor %>% ggplot(aes(top_mirna, top10_9 	))+geom_point()+geom_smooth(method="lm")
p9 <- mirna_hor %>% ggplot(aes(top_mirna, top10_10 	))+geom_point()+geom_smooth(method="lm")


grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9,ncol=3)


#####neg-correlation#####

if(dir.exists("mir_gene_pearson")){
  dir.create("mir_gene_pearson")
}



setwd("D:/GoogleDrive/Labor/DR/Andmed/DESeq2/mir_gene_pearson")



# 
# mirnaNR=0
# 
# 
# for(i in row.names(mirna_ver)){
#   mirnaNR=mirnaNR+1
#   print(mirnaNR)
#   sink(paste0("geen_",i,"_test_pearson.txt"),type="output")
#   print(i)
#   for(j in row.names(gene_ver)){
#     suhe.lm=cor(as.numeric(mirna_ver[i,]), as.numeric(gene_ver[j,],method="pearson"))
#     print(suhe.lm)
#   }
#   sink()
# }
# 
# 
# 
# ####tulbad järjest#######
# 
# failid=list.files()
# koik_andmed=read.table(fail, header=T) ###algse matrixi suuruse esitamine
# for (fail in failid) {
#   id=unlist(strsplit(fail, "[_]"))[2]
#   print(id)
#   andmed=read.table(fail, header=T)
#   koik_andmed=cbind(koik_andmed, andmed[2])
# }
# 
# write.table(koik_andmed, "mir_gene_pearson_koos.tsv", sep="\t", row.names=T, quote=F)

mir_gene_pearson <- read.table("mir_gene_pearson_koos.tsv", head=T, sep="\t")
mir_gene_pearson=mir_gene_pearson[!(names(mir_gene_pearson) %in% c("X.1.", "miR.001.1"))]
mir_gene_pearson <- setDT(mir_gene_pearson, keep.rownames = T)
mir_gene_pearson$rn <- rn

head(mir_gene_pearson[1:10,1:7])

plot(gene_hor[,"Geen0003"], mirna_hor$miR.003)
suhe.lm=cor(as.numeric(mirna_hor$miR.003), as.numeric(gene_hor[,"Geen0003"],method="pearson"))


#### parema p-valuega mirna-gene####
colSums(mir_gene_pearson < 0)

### plot####
mir_gene_pearson<-arrange(mir_gene_pearson, miR.003)

mir_gene_pearson[1:10,1:7]

top_mirna <- mirna_hor[mir_mir[1,1]]

top10_2 <- gene_hor[,mir_gene_pearson[1,1]]
top10_3 <- gene_hor[,mir_gene_pearson[2,1]]
top10_4 <- gene_hor[,mir_gene_pearson[3,1]]
top10_5 <- gene_hor[,mir_gene_pearson[4,1]]
top10_6 <- gene_hor[,mir_gene_pearson[5,1]]
top10_7 <- gene_hor[,mir_gene_pearson[6,1]]
top10_8 <- gene_hor[,mir_gene_pearson[7,1]]
top10_9 <- gene_hor[,mir_gene_pearson[8,1]]
top10_10 <- gene_hor[,mir_gene_pearson[9,1]]


p1 <- mirna_hor %>% ggplot(aes(top_mirna, top10_2 	))+geom_point()+geom_smooth(method="lm")
p2 <- mirna_hor %>% ggplot(aes(top_mirna, top10_3 	))+geom_point()+geom_smooth(method="lm")
p3 <- mirna_hor %>% ggplot(aes(top_mirna, top10_4 	))+geom_point()+geom_smooth(method="lm")
p4 <- mirna_hor %>% ggplot(aes(top_mirna, top10_5 	))+geom_point()+geom_smooth(method="lm")
p5 <- mirna_hor %>% ggplot(aes(top_mirna, top10_6 	))+geom_point()+geom_smooth(method="lm")
p6 <- mirna_hor %>% ggplot(aes(top_mirna, top10_7 	))+geom_point()+geom_smooth(method="lm")
p7 <- mirna_hor %>% ggplot(aes(top_mirna, top10_8 	))+geom_point()+geom_smooth(method="lm")
p8 <- mirna_hor %>% ggplot(aes(top_mirna, top10_9 	))+geom_point()+geom_smooth(method="lm")
p9 <- mirna_hor %>% ggplot(aes(top_mirna, top10_10 	))+geom_point()+geom_smooth(method="lm")


grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9,ncol=3)

