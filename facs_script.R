library('flowCore')
library('gplots')
library("dplyr")
library("ggplot2")
#### parameters ####
setwd("/Users/kanferg/Desktop/Gil_LabWork/AIMS/AIMS_090518/Facs_analysis/selected_data/")
dir()
f.name <- "Specimen_001_TFeb-H-1_001_001"  # name of the file you want to analyze, file must have extension ".FCS"

#fcm <- read.FCS(paste0(f.name,".fcs"), transformation = "scale",alter.names = TRUE)
fcm <- read.FCS(paste0(f.name,".fcs"),alter.names = TRUE)

fcm <- read.FCS(paste0(f.name, '.FCS'))
fcm <- as.data.frame((exprs(fcm)))
summary(fcm$`G610-A`)

#colfunc <- colorRampPalette(c("white", "lightblue", "green", "yellow", "red"))
colfunc <- colorRampPalette(c("white","red","yellow","green","lightblue"))

colnames(fcm)
ggplot(fcm, aes(x=fcm$`B530-A`, y=fcm$`BV421-A`)) +
  ylim(-100,10000) +
  xlim(-10000,10000) +
  stat_density2d(geom="tile", aes(fill = ..density..), contour = F) +
  scale_fill_gradientn(colours=colfunc(400)) + # gives the colour plot
  geom_density2d(colour="black", bins=10)

# draws the lines inside
  theme(axis.title.x =element_text(size = 14))+
  theme(axis.title.y =element_text(size = 14))+
  theme(axis.text = element_text(size=12))+
  theme(axis.line = element_line(colour = "black",size=1), panel.border = element_blank(),panel.background = element_blank())+
  # theme(axis.text.x = element_text( hjust = 1))+
  labs(title="", x="610", y="450")+ 
  #scale_x_continuous(breaks = (seq(0, 15, by=1)))+
  # geom_vline(xintercept = vec_zscore, linetype="dotted", size=  0.5, color = "red")+
  #geom_vline(xintercept = non_target.log2.norm, linetype="dotted", size=  0.5, color = "black")+  
  #geom_text(aes(x=0, label="Activated cell: 1,132#", y=3000),  text=element_text(size=14))+
  #geom_text(aes(x=0, label="Unactivated cells: 177,038#", y=3500), text=element_text(size=14))+ 
  #geom_vline(xintercept = 2400, linetype="dotted", color = "blue", size=1.5)+
  theme(panel.border = element_blank(),panel.background = element_blank())+
  theme(panel.grid.minor = element_line(colour = "white"))
setwd("/Users/kanferg/Desktop/Gil_LabWork/AIMS/AIMS_090518/Facs_analysis/Chart/")
ggsave(last_plot(), filename = "H1.svg", width=5.5, height=5.5, dpi = 72) 


ggplot(fcm,aes(fcm$`G610-A`))+
  geom_histogram(binwidth = 1)+
  geom_density(aes(y=1 * ..count..))

g1<-ggplot(fcm,aes(fcm[,10]))+
  geom_histogram(aes(y = ..count..,color = "red"),show.legend=F) +
  scale_color_manual(values = c("#FEE08B","red"),alpha(0.1,0.9))+                 
  geom_density()+
  theme(axis.title.x =element_text(size = 14))+
  theme(axis.title.y =element_text(size = 14))+
  theme(axis.text = element_text(size=12))+
  theme(axis.line = element_line(colour = "black",size=1), panel.border = element_blank(),panel.background = element_blank())+
  # theme(axis.text.x = element_text( hjust = 1))+
  labs(title="", x="mCharry Signal, AU", y="FACS intensity signal frequency, AU")+ 
  # geom_vline(xintercept = vec_zscore, linetype="dotted", size=  0.5, color = "red")+
  #geom_vline(xintercept = non_target.log2.norm, linetype="dotted", size=  0.5, color = "black")+  
  #geom_text(aes(x=0, label="Activated cell: 1,132#", y=3000),  text=element_text(size=14))+
  #geom_text(aes(x=0, label="Unactivated cells: 177,038#", y=3500), text=element_text(size=14))+ 
  #geom_vline(xintercept = 2400, linetype="dotted", color = "blue", size=1.5)+
  theme(panel.border = element_blank(),panel.background = element_blank())+
  theme(panel.grid.minor = element_line(colour = "white"))
ggsave(last_plot(), filename = "H1.svg", width=5.5, height=5.5, dpi = 72) 

g1 + scale_x_continuous(limits=c(0,1000))


ggplot(fcm, aes(x=fcm$`G610-A`, y=fcm$`V450-A`)) +
  ylim(-100, 200) +
  xlim(-100,250) +
  stat_density2d(geom="tile", aes(fill = ..density..), contour = F) +
  scale_fill_gradientn(colours=colfunc(400)) + # gives the colour plot
  geom_density2d(colour="black", bins=9) # draws the lines inside

ggplot(fcm, aes(x=fcm$`G610-A`, y=fcm$`V450-A`)) +
  ylim(-100, 200) +
  xlim(-100,250) +
  stat_density2d(geom="tile", aes(fill = ..density..), contour = F) +
  scale_fill_gradientn(colours=colfunc(400)) + # gives the colour plot
  geom_density2d(colour="black", bins=9)


# fcm <- na.omit(fcm)
# head(fcm,5)
# str(fcm)
summary(fcm$`BV421-A`)
summary(fcm$`G610-A`)
# 
# #### analyze file and make plot ####
# length(fcm$`V450-A`)
# length(fcm$`G575-A`)
# #normolise
# 
# fcm$`V450-A`<-fcm$`V450-A`- (min(fcm$`V450-A`))
# fcm$`V450-A`<-fcm$`V450-A`/max(fcm$`V450-A`)
# fcm$`G575-A`<-fcm$`G575-A`- (min(fcm$`G575-A`))
# fcm$`G575-A`<-fcm$`G575-A`/max(fcm$`G575-A`)
# 
# 
# max450<-max(fcm$`V450-A`)
# max575<-max(fcm$`G575-A`)
# str(fcm)
# 
# 
# ggplot(fcm,(aes(x=fcm$`V450-A`)))+
#   geom_histogram()+
#   xlim(0,1)

#facs analysis

th<-0.04
fcm$labale<-ifelse(fcm$`G610-A` > 0.002,"yes","")

g<-ggplot(fcm,(aes(x=log2(fcm$`V450-A`), y=log2(fcm$`G610-A`))))+
  geom_point(aes(color = fcm$labale),size=0.5,show.legend=F)+
  scale_color_manual(values = c("#FEE08B","palegreen4"),alpha(0.1,0.9))+
  xlim(-8,-4)+
  ylim(-12,-4)
g  + stat_density_2d(aes(fill = ..level..), geom="polygon",show.legend=F)+
  scale_fill_gradient(low="#ABDDA4", high="#FEE08B",show.legend=F)+
theme(axis.title.x =element_text(size = 14))+
  theme(axis.title.y =element_text(size = 14))+
  theme(axis.text = element_text(size=12))+
  theme(axis.line = element_line(colour = "black",size=1), panel.border = element_blank(),panel.background = element_blank())+
  # theme(axis.text.x = element_text( hjust = 1))+
  labs(title="", x="sgRNA BFP Signal, AU", y="mCharry Signal, AU")+ 
  #scale_x_continuous(breaks = (seq(0, 15, by=1)))+
  # geom_vline(xintercept = vec_zscore, linetype="dotted", size=  0.5, color = "red")+
  #geom_vline(xintercept = non_target.log2.norm, linetype="dotted", size=  0.5, color = "black")+  
  geom_text(aes(x=-6, label="Number of activated cell: 1,132", y=-6.5),  text=element_text(size=14))+
  geom_text(aes(x=-6, label="Number of unactivated cells: 177,038", y=-6),   vjust = 1.2, text=element_text(size=14))+ 
  theme(panel.border = element_blank(),panel.background = element_blank())+
  theme(panel.grid.minor = element_line(colour = "white"))

ggsave("Specimen_001_H5_001.svg", width = 5, height =4, units = "in",dpi = 720)

ggsave("Specimen_001_H5_001.png",width = 5, height =4, units = "in",dpi = 320)



fcm$labale<-ifelse(fcm$`G610-A` > 0.002,"yes","")
df<-as.data.frame(cbind(fcm$`V450-A`,fcm$`G610-A`))
str(df)
colnames(df)<-c("BFP","mCH")
df$logBFP<-log2(df$BFP)
df$logmCH<-log2(df$mCH)
ind<-which(!is.na(df$logBFP))
df<-df[ind,]
ind<-which(!is.na(df$logmCH))
df<-df[ind,]
df$log2BFP_norm<-ifelse(df$logBFP < -8.85 ,"NA" , df$logBFP)

ind<-which(!df$log2BFP_norm=="NA")
df<-df[ind,]
ind<-which(!is.infinite(df$log2BFP_norm))
df<-df[ind,]
df$labale<-ifelse(df$mCH > 0.002,"yes","")


g<-ggplot(df,(aes(x=df$logBFP, y=df$logmCH)))+
  geom_point(aes(color = df$labale),size=0.5,show.legend=F)+
  scale_color_manual(values = c("#FEE08B","palegreen4"),alpha(0.1,0.9))
  xlim(-8,-4)+
  ylim(-16,-4)
g  + stat_density_2d(aes(fill = ..level..), geom="polygon")+
  scale_fill_gradient(low="#ABDDA4", high="#FEE08B")+
  theme(axis.title.x =element_text(size = 14))+
  theme(axis.title.y =element_text(size = 14))+
  theme(axis.text = element_text(size=12))+
  theme(axis.line = element_line(colour = "black",size=1), panel.border = element_blank(),panel.background = element_blank())+
  # theme(axis.text.x = element_text( hjust = 1))+
  labs(title="", x="sgRNA BFP Signal, A.U.", y="mCherry Signal (PA), A.U.")+ 
  #scale_x_continuous(breaks = (seq(0, 15, by=1)))+
  # geom_vline(xintercept = vec_zscore, linetype="dotted", size=  0.5, color = "red")+
  #geom_vline(xintercept = non_target.log2.norm, linetype="dotted", size=  0.5, color = "black")+  
  geom_text(aes(x=-4, label="Number of activated cell: 1,132", y=-3),  text=element_text(size=14))+
  geom_text(aes(x=-4, label="Number of unactivated cells: 177,038", y=-4),   vjust = 1.2, text=element_text(size=14))+ 
  theme(panel.border = element_blank(),panel.background = element_blank())+
  theme(panel.grid.minor = element_line(colour = "white"))

ggsave("Facs.svg", width = 6, height =4, units = "in",dpi = 720)

ggsave("FACS.png",width = 6, height =4, units = "in",dpi = 72)
  
  
  setwd("/Users/kanferg/Desktop/Gil_LabWork/AIMS/AIMS_090518/Facs_analysis/Fig1/")
f.name <- "Specimen_0620"  # name of the file you want to analyze, file must have extension ".FCS"

fcm <- read.FCS(paste0(f.name,".fcs"), transformation = "scale")

fcm <- as.data.frame((exprs(fcm)))




## eliminate values that are below or equal to thresholds you
## defined above


th<-0.04
fcm$labale<-ifelse(fcm$`G610-A` > 0.002762136,"yes","")

g<-ggplot(fcm,(aes(x=log2(fcm$`V450-A`), y=log2(fcm$`G610-A`))))+
  geom_point(aes(color = fcm$labale),size=0.5,show.legend=F)+
  scale_color_manual(values = c("#FEE08B","palegreen4"),alpha(0.1,0.9))
  xlim(-12,-3)
  ylim(-12,-2.5)
g  + stat_density_2d(aes(fill = ..level..), geom="polygon")+
  scale_fill_gradient(low="#ABDDA4", high="#FEE08B")+
  theme(axis.title.x =element_text(size = 14))+
  theme(axis.title.y =element_text(size = 14))+
  theme(axis.text = element_text(size=12))+
  theme(axis.line = element_line(colour = "black",size=1), panel.border = element_blank(),panel.background = element_blank())+
  # theme(axis.text.x = element_text( hjust = 1))+
  labs(title="", x="sgRNA BFP Signal, AU", y="mCharry Signal, AU")+ 
  #scale_x_continuous(breaks = (seq(0, 15, by=1)))+
  # geom_vline(xintercept = vec_zscore, linetype="dotted", size=  0.5, color = "red")+
  #geom_vline(xintercept = non_target.log2.norm, linetype="dotted", size=  0.5, color = "black")+  
  geom_text(aes(x=-5, label="Number of activated cell: 1,132", y=-2.5),  text=element_text(size=14))+
  geom_text(aes(x=-5, label="Number of unactivated cells: 177,038", y=-3.5),   vjust = 1.2, text=element_text(size=14))+ 
  theme(panel.border = element_blank(),panel.background = element_blank())+
  theme(panel.grid.minor = element_line(colour = "white"))



#############Creat an histogram###############

low_mch<-rnorm(187851,1200,sd=450)
high_mch<-rnorm(1230,3000,50)
mch<-c(low_mch,high_mch)
df<-as.data.frame(mch)
df$labale<-ifelse(df$mch>2500,"yes","")

ggplot(df,aes(df$mch))+
  geom_histogram(aes(y = ..count..,color = df$labale), binwidth=20,size=0.5,show.legend=F) +
  scale_color_manual(values = c("#FEE08B","red"),alpha(0.1,0.9))+                 
  geom_density()+
  theme(axis.title.x =element_text(size = 14))+
  theme(axis.title.y =element_text(size = 14))+
  theme(axis.text = element_text(size=12))+
  theme(axis.line = element_line(colour = "black",size=1), panel.border = element_blank(),panel.background = element_blank())+
  # theme(axis.text.x = element_text( hjust = 1))+
  labs(title="", x="mCharry Signal, AU", y="FACS intensity signal frequency, AU")+ 
  #scale_x_continuous(breaks = (seq(0, 15, by=1)))+
  # geom_vline(xintercept = vec_zscore, linetype="dotted", size=  0.5, color = "red")+
  #geom_vline(xintercept = non_target.log2.norm, linetype="dotted", size=  0.5, color = "black")+  
  geom_text(aes(x=0, label="Activated cell: 1,132#", y=3000),  text=element_text(size=14))+
  geom_text(aes(x=0, label="Unactivated cells: 177,038#", y=3500), text=element_text(size=14))+ 
  geom_vline(xintercept = 2400, linetype="dotted", color = "blue", size=1.5)+
  theme(panel.border = element_blank(),panel.background = element_blank())+
  theme(panel.grid.minor = element_line(colour = "white"))
setwd( "/Users/kanferg/Desktop/Gil_LabWork/AIMS/AIMS_090518/Facs_analysis/Fig1/")
ggsave(last_plot(), filename = "facs_0921.png", width=5.5, height=5.5, dpi = 72) 



             
  














