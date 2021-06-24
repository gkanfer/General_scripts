#install.packages('VennDiagram')
library(VennDiagram)



topUMI<-read.table("topTirGene_UMIanalysis.txt")
topUndetermind<-read.table("topTirGene_undertermaind.txt")



venn.diagram(
  x = list(as.vector(topUMI) , as.vector(topUndetermind)),
  category.names = c("UMI","Undetrmiand"),
  filename = "vennn_diagramm.png",
  # output = TRUE ,
  # imagetype="png" ,
  # height = 480 ,
  # width = 480 ,
  # resolution = 300,
  # compression = "lzw",
  # lwd = 2,
  # lty = 'blank',
  fill = c("palegreen4","lightblue"),
  cex = 1,
  fontface = "bold"
  # fontfamily = "sans",
  # cat.cex = 0.6,
  # cat.fontface = "bold",
  # cat.default.pos = "outer",
  # cat.pos = c(-27, 27, 135),
  # cat.dist = c(0.055, 0.055, 0.085),
  # cat.fontfamily = "sans",
  # rotation = 1
)
#where on the plot the over lapping genes are?
#first find the overlapp


intersect.genes<-Reduce(intersect, list(topUMI,topUndetermind))
