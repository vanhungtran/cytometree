


logTrans <- function(data)
{
  min_data <- min(data)
  if(min_data > 0)
  {
    return(log(data))
  }
  log(data + abs(min_data) + 1)
}





# set work drive 


#setwd("E:/vhtran/mydf/Bcell")



# Import the required packages

# if (!requireNamespace("BiocManager"))
#   install.packages("BiocManager")
# BiocManager::install()
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager") 
# a
# BiocManager::install("flowCore", version = "3.8")
# a
# BiocManager::install("openCyto", version = "3.8") 
# a
# BiocManager::install("ggcyto", version = "3.8")
# a
library(flowCore) 
library(ggcyto)
library(tidyr)
library(cytometree)
library(flowWorkspace)


rm(list=ls())



# file name for working


# file.name = system.file("ACC2001_Tube_001", 
#              "ACC2001_Tube_002",
#              "ACC2001_Tube_003",package="flowCore")
#########################################################################
# file names with problem , unimporttable 

#file.name = "ACC2093_Tube_001.fcs"
# ACC2123_Tube_008.fcs  /nData may be truncated!
# ACC2139_Tube_009.fcs 
# ACC2255_Tube_005.fcs'
#####################################################################"

#ACC2695_Tube_002.fcs

# read FCS file 
#file.name = "ID01_02_104_ACC2540_Tube_002.fcs" # tube 3

file.name = "ID01_03_000_ACC2432_Tube_002.fcs"




FCS <- read.FCS("D:/Bcell/ID01_03_000_ACC2432_Tube_002.fcs")
T2 <- FCS@exprs
head(T2)


colnames(T2)[which(colnames(T2)=="FSC-A")] <- "FSC-A"
colnames(T2)[which(colnames(T2)=="FITC-A")] <- "IgD"
colnames(T2)[which(colnames(T2)=="PE-A")] <- "CD24"
colnames(T2)[which(colnames(T2)=="PE-Texas Red-A")] <- "CD19"
colnames(T2)[which(colnames(T2)=="PE-Cy5-A")] <- "CD20"
colnames(T2)[which(colnames(T2)=="PE-Cy7-A")] <- "CD38"
colnames(T2)[which(colnames(T2)=="APC-A")] <- "CD138"
colnames(T2)[which(colnames(T2)=="APC-Cy7-A")] <- "CD27"
colnames(T2)[which(colnames(T2)=="Pacific Orange-A")] <- "CD45"



# P <- ncol(T2)
# par(mfrow=c(3,3))
# for(i in 1:P)
# {
#   print(sum(!T2[,i])/dim(T2)[1] )
#   plot(density(T2[,i]))
# }


T2 = T2[,1:10]

# Tree <- CytomeTree(T2, minleaf = nrow(T2)*0.01)
# save(Tree, file = "CytoFTree.Rdata")
# Annot <- Annotation(Tree,plot=FALSE, K2markers = colnames(T2))
# save(Annot, file = "CytoFAnnot.Rdata")
# 

M = T2
library(mclust)
CytEMRes <- CytEM(M, 1:1967160, 1, level, 0.1)
CytEMRes
summary(CytEMRes$mark_not_dis)
length(M)


