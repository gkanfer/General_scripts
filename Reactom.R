# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("clusterProfiler")
# 
# install_github("GuangchuangYu/bitr")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# bitr(geneList$V2)
# BiocManager::install("ReactomePA")
# #devtools::install_github("r-lib/rlang", build_vignettes = TRUE)
# library(ReactomePA)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("org.Hs.eg.db")

library(BiocManager)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)

load("~/Desktop/Gil_LabWork/AIMS/AIMS_090518/Tfeb_NGS_analysis/second_analysis_H1_H3_H4/data_072919.RData")
df.100000$Gene<-as.character(df.100000$Gene)
gene.ind<-which(df.100000$PV < 0.1)
setwd("~/Desktop/Gil_LabWork/AIMS/AIMS_090518/Tfeb_NGS_analysis/file_list/list_analysis_3/")
Genes<-df.100000[gene.ind,1]
#saveRDS(Genes,"genelist.txt")
geneList_vec<-readRDS("genelist.txt")

data(geneList)
####converting using david
library("dplyr")
geneList<-read.table(file = "From_david.txt"       ,fill = T,stringsAsFactors = F)
df<-geneList %>% filter(V3=="Homo")
row.names(df)<-df$V1
df.filter<-df[,c(1,2)]
df.filter<-as.data.frame(df.filter[,-1],stringsAsFactors=F)
row.names(df.filter)<-df$V1

#df.filter<-t(df.filter)

# # # geneList<-read.csv(file = "uniprot_table.csv" , stringsAsFactors = F)
# # #
# # # str(geneList)
# # # nrow(geneList)
# # #
# # # human<-rep("HUMAN",1081)
# # # geneList<-
# # # geneList_vec<-geneList$Entry
# # # EID<-bitr(geneList_vec,"UNIPROT","ENTREZID",org.Hs.eg.db)
# # # str(EID)
# # write.table(geneList_vec,"uniprot.txt",row.names = F,quote = F)
# geneList<-read.table(file = "conv_entreasID_1.txt",fill = T,stringsAsFactors = F)
# str(geneList)
# EID.vec<-EID$ENTREZID

EID.vec<-as.numeric(EID.vec)
dee <- names(df.filter[,1])[abs(df.filter) > 1.5]
head(de)
de<-df.filter$`df.filter[, -1]`

x <- enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
head(as.data.frame(x))

barplot(x, showCategory=5)
dotplot(x, showCategory=5)
emapplot(x)
cnetplot(x, categorySize="pvalue", circular=F,showCategory = 5)
cnetplot.enrichResult(x, showCategory = 5, foldChange = NULL,layout = "kk", colorEdge = FALSE, circular = FALSE, node_label = TRUE)
