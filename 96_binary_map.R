rm(list=ls())
library(dplyr)
library(tidyr)


#setwd("/Volumes/SH_colocolisation/200319_DKO_TKO/20200319_123402_655")
setwd("/Volumes/SH_colocolisation/200319_DKO_TKO/20200319_120619_533/")


txt<-dir()[grep(".txt$",dir())]

df.plate1<-NULL
for (i in 1:length(txt)){
  df<-try(read.table(txt[i],stringsAsFactors = F))
  if (inherits(df, "try-error")) {next()}
  else {
    #str(df)
    print(i)
    #print(df)
    #browser()
    colnames(df)<-df[1,]
    df<-df[c(2:nrow(df)),]
    row.names(df)<-c(1:nrow(df))
    df[,c(1:2)]<-apply(df[,c(1:2)],2,as.numeric)
    #quntile.25<-quantile(df$GFP_intensity)["25%"]
    # nrow(df)
    # #df<-df[which(df$GFP_intensity > quntile.25),]
    # #df<-df[which(df$m2 > 0),]
    # df<-arrange(df,desc(m1))
    # df<-df[c(1:10),]
    #str(df)
    ncell<-nrow(df)
    p<-length(which(df$Foci_binery==1))
    df.plate1<-rbind(df.plate1,c(txt[i],ncell,p))
    rm(df,ncell,p)    
  }
  
}
str(df.plate1)
df.plate1<-as.data.frame(df.plate1,stringsAsFactors = F)

df.plate1[,c(2:ncol(df.plate1))]<-apply(df.plate1[,c(2:ncol(df.plate1))],2,as.numeric)
colnames(df.plate1)<-c("Well","cell_number","postives")
df.plate1$Well<-gsub("_Point.*","",df.plate1$Well)
df.plate1$Well<-gsub("Well","",df.plate1$Well)

############################################################################################
#' compere the two samples df.palte2 / df.plate1
############################################################################################  
mat.m1<-matrix(NA, nrow = 6, ncol = 10)
row<-c("B","C","D","E","F","G")
column<-c("02","03","04","05","06","07","08","09","10","11")

rownames(mat.m1)<-row
colnames(mat.m1)<-column

mat.postive_presantage<-mat.m1
mat.sd_Postive_pracantage<-mat.m1
mat.cell_number<-mat.m1

str(df.plate1)

uniquewell<-unique(df.plate1$Well)
for (m in uniquewell){
  txt<-m
  c<-gsub("\\D","",txt,perl = T)
  r<-gsub("\\d","",txt,perl = T)
  ind.r<-grep(r,rownames(mat.m1)) 
  ind.c<-grep(c,colnames(mat.m1)) 
  temp<-df.plate1[df.plate1$Well==m,]
  if (nrow(temp) < 1){ 
    mat.m2[ind.r,ind.c]<-0
    mat.m1[ind.r,ind.c]<-0
  } else {
    precn<-mean(temp$postives/temp$cell_number,na.rm = F)
    STD<-sd(temp$postives/temp$cell_number,na.rm = F)
    cn<-sum(temp$cell_number)
    mat.postive_presantage[ind.r,ind.c]<-precn
    mat.sd_Postive_pracantage[ind.r,ind.c]<-STD
    mat.cell_number[ind.r,ind.c]<-cn
    rm(precn,STD)
  }
  
}
setwd("/Volumes/data/temp/SHS/output/")
write.csv(mat.postive_presantage,"Plate1_032820_postive_precentage.csv",quote = F,row.names = T,col.names = T)
write.csv(mat.sd_Postive_pracantage,"Plate1_032820_postive_precentage_SD.csv",quote = F,row.names = T,col.names = T)
write.csv(mat.cell_number,"Plate1_032820_cellnumber.csv",quote = F,row.names = T,col.names = T)




