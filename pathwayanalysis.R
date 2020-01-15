# source("http://bioconductor.org/biocLite.R")
# biocLite("pathview")
library(pathview)
data(gse16873.d)
data(demo.paths)
require(org.Hs.eg.db)
library(dplyr)
library(edgeR)
#------------------------------------------------------------------------------
# detect gene symble
#------------------------------------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("mygene")
library(mygene)
setwd("/Users/kanferg/Desktop/Gil_LabWork/AIMS/AIMS_090518/Tfeb_NGS_analysis/011020_summary/Export/reactom/")
geneList_vec<-readRDS("genelist.rds")
length(geneList_vec)
getGenes(geneList_vec, fields='symbol')
library("org.Hs.eg.db")
x<-columns(org.Hs.eg.db)

for (i in 1:length(x)){
  geneSymbols <- try(mapIds(org.Hs.eg.db, keys=geneList_vec, column="SYMBOL", keytype=x[i], multiVals="first"))
  if (inherits(geneSymbols, "try-error")) {
    next()
  } else {
    print(x[i])
    print(geneSymbols)}
}

geneSymbols <- mapIds(org.Hs.eg.db, keys=geneList_vec, column="ENTREZID", keytype="ALIAS", multiVals="first") 
inds <- which(!is.na(geneSymbols))
found_genes <- geneSymbols[inds]  

#------------------------------------------------------------------------------
# creat a matrix contains the enraseID in the rows and FC and FDR in the columns
#------------------------------------------------------------------------------
#loading roast export
setwd("~/Desktop/Gil_LabWork/AIMS/AIMS_090518/Tfeb_NGS_analysis/011020_summary/export/roast/")
rm(roast)
dir()
roast.export<-rbind(readRDS("roast_h1.rds"),
                    readRDS("roast_h2.rds"),
                    readRDS("roast_h3.rds"),
                    readRDS("roast_h4.rds"),
                    readRDS("roast_h5.rds"),
                    readRDS("roast_h6.rds"),
                    readRDS("roast_h7.rds"))


roast.sel<-roast.export %>% filter(FDR < 0.3, FC>1)
str(roast.sel)
# creat genesymble matrix

geneSymbols <- mapIds(org.Hs.eg.db, keys=roast.sel$Gene, column="ENTREZID", keytype="ALIAS", multiVals="first") 
inds <- which(!is.na(geneSymbols))
found_genes <- geneSymbols[inds]  



de<-NULL
for (i in 1:nrow(roast.sel)){
    #browser()
    ind.sel<-which(roast.sel$Gene %in% attr(geneSymbols,"names")[i])
    if (length(ind.sel)==1){
      gene<-attr(geneSymbols,"names")[i]
      entrezID<-geneSymbols[[i]]
      FC<-as.numeric(roast.sel[ind.sel,10])
      FDR<-as.numeric(roast.sel[ind.sel,7])
      de<-rbind(de,c(gene,entrezID,FC,FDR))
      rm(gene,entrezID,FC,FDR)
    }else
      { next() }
}

de<-as.data.frame(de,stringsAsFactor=F)
str(de)
de[,1]<-as.character(de[,1])
de[,2]<-as.character(de[,2])
de[,3]<-as.numeric(as.character(de[,3]))
de[,4]<-as.numeric(as.character(de[,4]))
mat.de<-matrix(NA,nrow = nrow(de),ncol = 2)
for (i in 1:2){
  mat.de[,i]<-de[,3:4][,1]
}

mat.de<-as.matrix(de[,c(3:4)])
rownames(mat.de)<-de[,2]
str(mat.de)
mat.de[,1]<-as.numeric(mat.de[,1])
mat.de[,2]<-as.numeric(mat.de[,2])
keg <- kegga(de$V2,species="Hs")
topKEGG(keg, n=c(10))
setwd("/Users/kanferg/Desktop/Gil_LabWork/AIMS/AIMS_090518/Tfeb_NGS_analysis/011020_summary/Export/roast/kegg/")
#dir.create("kegg")
pv.out <- pathview(gene.data = mat.de[,1], pathway.id = "01521",species = "hsa")
# saveRDS(de,"DE_table.rds")
# saveRDS(mat.de,"DE_matrix.rds")


#------------------------------------------------------------------------------
# creat a matrix contains the enraseID in the rows and FC and FDR in the columns
#------------------------------------------------------------------------------
#loading roast Import
setwd("~/Desktop/Gil_LabWork/AIMS/AIMS_090518/Tfeb_NGS_analysis/011020_summary/Import/roast/")
rm(roast)
dir()
roast.import<-rbind(readRDS("roast_h1.rds"),
                    readRDS("roast_h2.rds"),
                    readRDS("roast_h3.rds"),
                    readRDS("roast_h4.rds"),
                    readRDS("roast_h5.rds"),
                    readRDS("roast_h6.rds"),
                    readRDS("roast_h7.rds"))
setwd("~/Desktop/Gil_LabWork/AIMS/AIMS_090518/Tfeb_NGS_analysis/011020_summary/Import/roast/kegg/")
#dir.create("kegg")

roast.sel<-roast.import %>% filter(FDR < 0.3, FC>1)
nrow(roast.sel)

geneSymbols <- mapIds(org.Hs.eg.db, keys=roast.sel$Gene, column="ENTREZID", keytype="ALIAS", multiVals="first") 
inds <- which(!is.na(geneSymbols))
found_genes <- geneSymbols[inds] 


de<-NULL
for (i in 1:nrow(roast.sel)){
  #browser()
  ind.sel<-which(roast.sel$Gene %in% attr(geneSymbols,"names")[i])
  if (length(ind.sel)==1){
    gene<-attr(geneSymbols,"names")[i]
    entrezID<-geneSymbols[[i]]
    FC<-roast.sel[ind.sel,10]
    FDR<-roast.sel[ind.sel,7]
    de<-rbind(de,c(gene,entrezID,FC,FDR))
    rm(gene,entrezID,FC,FDR)
  }else
  { next() }
}

de<-as.data.frame(de,stringsAsFactor=F)
str(de)
de$V1<-as.character(de$V1)
de$V2<-as.character(de$V2)
de$V3<-as.numeric(as.character(de$V3))
de$V4<-as.numeric(as.character(de$V4))

keg <- kegga(de$V2, species="Hs")
topKEGG(keg, n=15)
setwd("/Users/kanferg/Desktop/Gil_LabWork/AIMS/AIMS_090518/Tfeb_NGS_analysis/011020_summary/Import/roast/kegg/")
#dir.create("kegg")
pv.out <- pathview(gene.data = mat.de[,1], pathway.id = "05230",species = "hsa")




 

