setwd("/Users/kanferg/Desktop/Gil_LabWork/AIMS/deep_seq/rafferance_table/reffrancetable/weissman/")
#creat index rafrance seq
library("data.table")
#install.packages("seqRFLP")
library("seqRFLP")
#Extract Kinome (H1)
load("~/Desktop/Gil_LabWork/AIMS/deep_seq/rafferance_table/reffrancetable/weissman/data.RData")


library(dplyr)
# bob %>% mutate_if(is.factor, as.character) -> bob


df$seq<-df[,4]

df_h7<-df[Sublibrary=="h7"]
df_h7<-df_h7  %>% mutate_if(is.factor, as.character)


for (i in 1:length(df_h7[,1])){
  df_h7[i,12]<-paste(seqright,df_h7[i,4],seqleft,sep="")
  print(i)
}


gene_name<-c(df_h7[,1])
seq<-c(c(df_h7[,12]))

dfff_h7 <- data.frame(gene_name,seq)
dfff_h7$gene_name<- paste(">",dfff_h7$gene_name,sep="")
dfff_h7$seq<-as.character(dfff_h7$seq)

vec.h7<-c(dfff_h7[1,1],dfff_h7[1,2])
for (i in 1:nrow(dfff_h7)){
  vec.temp<-c(dfff_h7[i,1],dfff_h7[i,2])
  vec.h7<-c(vec.h7,vec.temp)
  print(i)
}



fileConn<-file("h7.fasta")
writeLines(vec.h7, fileConn)
close(fileConn)

file.show("h7.fasta")


#---------------------------------------------------------------
#-

####great

##secondtry
# setwd("/Users/kanferg/Desktop/Hawk")
# a<-read.csv("Mouse_GeCKOv2_Library_A_09Mar2015.csv")
# a$id<-paste(a$UID,a$gene_id,sep="_")         
# d<-a[,c(4,3)] # this reorganises the table to be unique_id and sequence
# write.table(d,"d.tab",sep="\t",quote=F,row.names=F) 

#try the samething
df_h2_new<-df_h2[,c(1,4)]
write.table(d,"df_h2.fast",sep="\t",quote=F,row.names=F) 
saveRDS(df_h2, "Library_table_h2.rds")

#copy all the sequances
string<-paste(df_h1[1,4],sep = "")
for (i in 2:100){
sstring<-paste(df_h1[i,4],sep = "")
string<-paste(string,sstring,sep = "")
    
}
write.table(string,"string.txt")


#####bed file compresion
setwd("~/Desktop/Gil_LabWork/AIMS/deep_seq/Analysis/032718")
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("rtracklayer", "AnnotationHub", "Rsamtools"))
library(rtracklayer)
library(AnnotationHub)
library(Rsamtools)

test_new<-read.table("new.bed")
str(test)

library(ggplot2)
gnew<-ggplot(test_new,(aes(x=V1)))+
  geom_histogram(binwidth=10)

?hist
hist(test[,2])  
setwd("~/Desktop/Gil_LabWork/AIMS/deep_seq/Analysis/032718/index10and14")
test<-read.table("new14.bed")
str(test)
write.csv(test,"kinome.csv")
g<-ggplot(test,(aes(x=V1)))+
  geom_histogram(binwidth=1)

gnew14<-ggplot(test,(aes(x=V1))) + 
        geom_density(aes(color = "tomato3", fill="tomato3", alpha=0.3))+
        theme(axis.title.x =element_text(size = 14))+
        theme(axis.title.y =element_text(size = 14))+
        theme(axis.text = element_text(size=12))+
        theme(axis.line = element_line(colour = "black",size=1), panel.border = element_blank(),panel.background = element_blank())+
        labs(y= "Frequency")+
        labs(x="reads")+
        theme(panel.border = element_blank(),panel.background = element_blank())
        theme(panel.grid.minor = element_line(colour = "white"))
  
        
        ggplot(test,(aes(x=V1, Y=test[,2]))) + 
          geom_bar(aes(color = "tomato3", fill="tomato3", alpha=0.3))+
          theme(axis.title.x =element_text(size = 14))+
          theme(axis.title.y =element_text(size = 14))+
          theme(axis.text = element_text(size=12))+
          theme(axis.line = element_line(colour = "black",size=1), panel.border = element_blank(),panel.background = element_blank())+
          labs(y= "Frequency")+
          labs(x="reads")+
          theme(panel.border = element_blank(),panel.background = element_blank())+
        theme(panel.grid.minor = element_line(colour = "white"))        

x<-rep(0,752)        
y<-rep("zero",752)
name<-names(test)
zerotable<-data.frame(x,y)        
colnames(zerotable)<-name
  
  
  
test_zero<-rbind(test,zerotable)        

ggplot(test_zero,(aes(x=V1)))+
  geom_density(aes(color = "tomato3", fill="tomato3", alpha=0.3))+
  theme(axis.title.x =element_text(size = 14))+
  theme(axis.title.y =element_text(size = 14))+
  theme(axis.text = element_text(size=12))+
  theme(axis.line = element_line(colour = "black",size=1), panel.border = element_blank(),panel.background = element_blank())+
  labs(y= "Frequency")+
  labs(x="reads")+
  theme(panel.border = element_blank(),panel.background = element_blank())
theme(panel.grid.minor = element_line(colour = "white"))



gnew<-ggplot(test_new,(aes(x=V1)))+
  geom_density(aes(color = "tomato3", fill="tomato3", alpha=0.3))+
  theme(axis.title.x =element_text(size = 14))+
  theme(axis.title.y =element_text(size = 14))+
  theme(axis.text = element_text(size=12))+
  theme(axis.line = element_line(colour = "black",size=1), panel.border = element_blank(),panel.background = element_blank())+
  labs(y= "Frequency")+
  labs(x="reads")+
  theme(panel.border = element_blank(),panel.background = element_blank())
theme(panel.grid.minor = element_line(colour = "white"))

        

        
        
        
        
# 
# 
# g<-ggplot(TSS, aes(x=name, y=TSS$`Number of cells`, fill=TSS$class))+
#   geom_bar(stat = "identity",position=position_dodge(),colour="black" )+
#   geom_text(aes(label=TSS$`Number of cells`), vjust=1.6, color="white", position = position_dodge(0.9), size=3.5)+
#   scale_fill_manual(values=c("palegreen4","tomato3"))+
#   # scale_fill_manual(values=c("palegreen4", "tomato3"))+
#   theme(axis.title.x =element_text(size = 14))+
#   theme(axis.title.y =element_text(size = 14))+
#   theme(axis.text = element_text(size=12))+
#   labs(y= expression(atop("Cells number", paste("(Cytosol Parkin, log)"))))+
#   labs(x="Precentage of sgPink1 expressing cells ")+
#   # scale_y_continuous( trans='log2')+
#   # theme(axis.line = element_line(colour = "black",size=1), panel.border = element_blank(),panel.background = element_blank())+
#   theme(panel.border = element_blank(),panel.background = element_blank())+
#   # coord_capped_cart(bottom="none",left=brackets_vertical())+
#   annotate(x=0, xend=0, y=8, yend=10000,  lwd=0.75, geom="segment")+
#   theme(panel.grid.minor = element_line(colour = "white"))
# ###change the factors in order to change the order
# TSS$name<-factor(TSS$name, levels=c("0.2%","0.5%", "2%", "10%"))
# TSS$class<-factor(TSS$class, levels=c("Cyto", "Cyto Pred"))
# 
# g<-ggplot(TSS, aes(x=name, y=TSS$`Number of cells`, fill=TSS$class))+
#   geom_bar(stat = "identity",position=position_dodge(),colour="black" )+
#   geom_text(aes(label=TSS$`Number of cells`), vjust=1.6, color="white", position = position_dodge(0.9), size=3.5)+
#   scale_fill_manual(values=c("palegreen4","tomato3"))+
#   theme(axis.title.x =element_text(size = 14))+
#   theme(axis.title.y =element_text(size = 14))+
#   theme(axis.text = element_text(size=12))+
#   labs(y= expression(atop("Cells number", paste("(Cytosol Parkin, log)"))))+
#   labs(x="Precentage of sgPink1 expressing cells ")+
#   coord_cartesian(ylim=c(0,1500))+
#   #scale_y_continuous( trans='log2')+
#   # theme(axis.line = element_line(colour = "black",size=1), panel.border = element_blank(),panel.background = element_blank())+
#   theme(panel.border = element_blank(),panel.background = element_blank())+
#   # coord_capped_cart(bottom="none",left=brackets_vertical())+
#   annotate(x=0, xend=0, y=8, yend=1500,  lwd=0.75, geom="segment")+
#   theme(panel.grid.minor = element_line(colour = "white"))
# 
# 
# g<-ggplot(temp, aes( y=temp[1,]))+
#   geom_point() +
#   theme_bw(base_size = 14) 
# geom_smooth(method = "auto", level=0)+
#   labs(y= expression(atop("Number of cells", paste("(% of Cytosol Parkin)"))))+labs(x="Time(min)") +
#   ggtitle("Phenotypic Window")+ 
#   geom_rect(aes(xmin = 110, xmax = 230, ymin = -Inf, ymax = Inf),fill = "gray95", alpha = 0.03)+
#   geom_point() + 
#   theme_bw(base_size = 12) +
#   geom_smooth(method = "auto", level=0)+
#   theme(plot.title = element_text(size = 12, face = "bold"))+
#   theme(plot.title = element_text(hjust = 0.3))+
#   theme(axis.title.x =element_text(size = 14))+
#   theme(axis.title.y =element_text(size = 14))+
#   theme(axis.text = element_text(size=12))+
#   scale_y_continuous(breaks = round(seq(min(dt$postive), max(dt$postive), by = 0.005),1),labels = scales::percent)+
#   scale_x_continuous(breaks = round(seq(min(dt$min), max(dt$min), by = 50),220))+
#   theme(axis.line = element_line(colour = "black",size=1),panel.border = element_blank(),panel.background = element_blank())+
#   theme(panel.grid.major = element_line(colour = "white"))+
#   theme(panel.grid.minor = element_line(colour = "white"))+
#   geom_vline(xintercept = c(110,230), linetype="dotted", size=  0.5, color = "red")  ###add vertical lines
# ggsave("test1.svg", plot = p, width =4 ,height =3.5)



