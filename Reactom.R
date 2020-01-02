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


# load("~/Desktop/Gil_LabWork/AIMS/AIMS_090518/Tfeb_NGS_analysis/second_analysis_H1_H3_H4/data_072919.RData")
# df.100000$Gene<-as.character(df.100000$Gene)
# gene.ind<-which(df.100000$PV < 0.1)
setwd("~/Desktop/Gil_LabWork/AIMS/AIMS_090518/Tfeb_NGS_analysis/122919_bow1_UMI_tryhard_redo/reactom/")
setwd("~/Desktop/Gil_LabWork/AIMS/AIMS_090518/Tfeb_NGS_analysis/122919_bow1_UMI_tryhard_redo/roast/")
#rm(roast)
dir()
roast.UMI<-rbind(readRDS("roast_h1.rds"),
                 readRDS("roast_h2.rds"),
                 readRDS("roast_h3.rds"),
                 readRDS("roast_h4.rds"),
                 readRDS("roast_h5.rds"),
                 readRDS("roast_h6.rds"),
                 readRDS("roast_h7.rds"))

roast.sel<-roast.UMI %>% filter(FC > 0) %>%  filter(FDR < 0.2)

getwd()
setwd("/Users/kanferg/Desktop/Gil_LabWork/AIMS/AIMS_090518/Tfeb_NGS_analysis/122919_bow1_UMI_tryhard_redo/reactom/")
Genes<-roast.sel$Gene
saveRDS(Genes,"genelist.rds")
geneList_vec<-readRDS("genelist.rds")
write.table(geneList_vec,"genelist.txt",quote = F,row.names = F,col.names = F)
data(geneList)
####converting using david
library("dplyr")
dir()
geneList<-read.table(file = "From_david.txt",fill = T,stringsAsFactors = F)
df<-geneList %>% filter(V3=="Homo")
row.names(df)<-df$V1
df.filter<-df[,c(1,2)]
df.filter<-as.data.frame(df.filter[,-1],stringsAsFactors=F)
row.names(df.filter)<-df$V1
str(df.filter)
#df.filter<-t(df.filter)

# # # geneList<-read.csv(file = "uniprot_table.csv" , stringsAsFactors = F)
# # #
# # # str(geneList)
# # # nrow(geneList)
# # #
# # # human<-rep("HUMAN",1081)
# # # geneList<-
# geneList_vec<-as.data.frame(geneList_vec,stringsAsFactor=F)
# str(geneList_vec)
# geneList_vec$geneList_vec<-as.character(geneList_vec$geneList_vec)
# 
# EID<-bitr(geneList_vec$geneList_vec,"UNIPROT","ENTREZID",org.Hs.eg.db)
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

barplot(x, showCategory=20)
dotplot(x, showCategory=20)
emapplot(x)
cnetplot(x, categorySize="pvalue", circular=F,showCategory = 5)
cnetplot.enrichResult(x, showCategory = 5, foldChange = NULL,layout = "kk", colorEdge = FALSE, circular = FALSE, node_label = TRUE)
