setwd("~/Desktop/temp")
source("https://bioconductor.org/biocLite.R")
#biocLite("ShortRead")
library(ShortRead)
library(ggseqlogo)
library(ggplot2)

setwd("/Users/kanferg/Desktop/Gil_LabWork/AIMS/deep_seq/Sequance_alighment_100219/")
fastqPath <- list.files("/Users/kanferg/Desktop/Gil_LabWork/AIMS/deep_seq/Sequance_alighment_100219/", pattern = "092719_207bp.fastq"  , full = TRUE)[1]
reads <- readFastq(fastqPath)
fqFile <- FastqFile(fastqPath)
fqFile
reads <- readFastq(fqFile)
seq<-sread(reads)[1:1340]
x<-1:length(seq)
for (i in 1:length(seq)) {
  x[i]<-paste(unlist(seq[i]))
  print(x[i])
}
# data(ggseqlogo_sample)
# str(df_h7)
# test<-df_h7[1:100,]

ggplot() + geom_logo( x ) + theme_logo()
ggplot() + geom_logo(x) + xlim(100,250)
ggplot() + xlim(180,210) + geom_logo(x)


constring<-consensusString(x) #get the consensue of char list of sequances
segment<-gsub("\\?","N",constring)



m <- consensusMatrix(x)[1:5,] # nucl. x position counts
df.m<-cbind(row.names(m),m)
df.concen<-vector()
for (i in 2:ncol(df.m)){
  x<-as.numeric(df.m[,i])
  max.x<-max(x)
  ind<-which(x==max.x)
  if (ind==1){temp<-"A"}
  if (ind==2){temp<-"C"}
  if (ind==3){temp<-"G"}
  if (ind==4){temp<-"N"}
  if (ind==5){temp<-"T"}
  df.concen<-c(df.concen,temp)
  sequance<-paste(df.concen,collapse = "")
}


polymorphic <- which(colSums(m != 0) > 1)

m[, polymorphic]
df.m<-as.data.frame(m,stringsASFactor=F)

string<-vector()
nucleotide<-row.names(m)
for (i in 1:ncol(df.m)){
  df.temp<-df.m[,i]
  ind<-which(df.temp>400)
  if ((length(ind))==0){
    string[i]<-nucleotide[4]
  } else {
    string[i]<-nucleotide[ind]
  } 
  
}
paste0(string,collapse = "")
