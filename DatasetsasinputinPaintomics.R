###########################################################
###########   DatasetsasinputinPaintomics.R     ###########
###########################################################
# Author: Leandro Balzano-Nogueira
# Genetics Institute, University of Florida (Gainesville)

# This script is to create the tables to include in paintomics for pathways analyzes

###########################################################
homedir<- "/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/" # Home directory where all your results are going to be contained
setwd(homedir)
getwd()
###########################################################
# Functions:
source("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/ScriptsForNComm/Tools/TEDDYtools2.R") # Location of TEDDYtools

PaintomicsSubsetter<- function (X=LM_globalTargets, Y=LM_globalEXPMAT, Z= metabolomicstargets, ZZ=globalmetabolomicsMatrix,
                                ageg=2, FAAb=0, timepoints=3, timeframe=c(12,9,6,3,0) ) {
  if (ageg==0 & FAAb==0) {
    print("Both matrices are exactly the same, No Subset was created")
    targetsGE = X
    datasetGE = Y[,as.character(targetsGE$sample_mask_id)]
    targetsMET = Z
    datasetMET = ZZ[,as.character(targetsMET$Sample.Mask.Id)]
  }
  if (ageg!=0 & FAAb==0) {
    targetsGE = filter(X,agegroup==ageg)
    datasetGE = Y[,as.character(targetsGE$sample_mask_id)]
    targetsMET = filter(Z,AgeGroup==ageg)
    datasetMET = ZZ[,as.character(targetsMET$Sample.Mask.Id)]
  }
  if (ageg==0 & FAAb!=0) {
    targetsGE = filter(X,FirstAAb==FAAb)
    datasetGE = Y[,as.character(targetsGE$sample_mask_id)]
    targetsMET = filter(Z,FirstAAbCC==FAAb)
    datasetMET = ZZ[,as.character(targetsMET$Sample.Mask.Id)]
  }
  if (ageg!=0 & FAAb!=0) {
    print("Warning: Subsetting by two variables could reduce too much the number of Invividuals")
    targetsGE = filter(X,agegroup==ageg,FirstAAb==FAAb)
    datasetGE = Y[,as.character(targetsGE$sample_mask_id)]
    targetsMET = filter(Z,AgeGroup==ageg,FirstAAbCC==FAAb)
    datasetMET = ZZ[,as.character(targetsMET$Sample.Mask.Id)]
    
  }
  # Get Patients with at least 3 time points  (N=136)
  #GeneExpression
  patientTable = table(targetsGE$mask_id)
  patients = names(patientTable)[patientTable>=timepoints]
  #length(patients)
  Xsmall = filter(X,mask_id %in% patients)
  
  #Metabolomics
  patientTableMET = table(targetsMET$Mask.Id)
  patients = names(patientTableMET)[patientTableMET>=timepoints]
  #length(patients)
  Zsmall = filter(Z,Mask.Id %in% patients)
  
  
  # Paintomics Genes tables
  V2_pltDF = c()
  V2_pltDFMET = c()
  for(timepoint in timeframe){
    #GeneExpression
    testTargets = filter(Xsmall,time==timepoint)
    meansDFv2 = getGroupMean(groupcolumn = "outcome",stargets = testTargets,sexprmat = Y)
    LogRv2 = as.numeric(meansDFv2[,3]) - as.numeric(meansDFv2[,2])
    names(LogRv2) = rownames(meansDFv2)
    V2_pltDF = cbind(V2_pltDF, LogRv2)
    
    #Metabolomics
    testTargetsMET = filter(Zsmall,Sconv==timepoint*-1)
    meansDFv2MET = getGroupMeanMET(groupcolumn = "Outcome",stargets = testTargetsMET,sexprmat = ZZ)
    #head(meansDFv2MET)
    LogRv2MET = as.numeric(meansDFv2MET[,3]) - as.numeric(meansDFv2MET[,2])
    names(LogRv2MET) = rownames(meansDFv2MET)
    V2_pltDFMET = cbind(V2_pltDFMET, LogRv2MET)
    
  }
  columnnames=NULL
  for (i in 1:length(timeframe)) {
    names<-paste("Time",timeframe[i], sep="_-")
    columnnames<-c(columnnames,names)
  }
  colnames(V2_pltDF) = columnnames
  colnames(V2_pltDFMET) = columnnames
  result<-list(targets=targets,datasetGE=datasetGE,datasetMET=datasetMET, paintomicsGETable= V2_pltDF,paintomicsMETTable= V2_pltDFMET)
  return(result)
  
}

###########################################################
# Data:
GE_Processed<-read.csv ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/GeneExpression/GE_Processed.csv",header = TRUE)
GE_Processed[1:10,1:10]

wholedescriptivevars<-read.csv ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/LM_globalTargets.csv",header = TRUE, row.names = 1)
wholedescriptivevars
AllmetabolitesConverted <- read.csv("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/AllmetabolitesConverted.txt", sep= "\t")

Allmetabolomicstogether <- read.delim("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/Allmetabolomicstogether.txt")

LM_globalTargets<-read.csv ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/LM_globalTargets.csv",header = TRUE, row.names = 1)
LM_globalTargets

###########################################################
# Libraries:
require("abind")
library(gplots)
library(data.table)
library("gdata")
require("VennDiagram")

###########################################################
#Retrieve Expression Data
GE_Processed[1:10,1:10]
preDat<-GE_Processed[-c(2:4),]
dim(preDat)
preDat[1:10,1:10]
rownames(preDat)<-preDat[,1];preDat<-preDat[,-1]
colnames(preDat)<-preDat[1,];preDat<-preDat[-1,]
preDat[1:10,1:10]

AgeFeatures = list("agegroup"=unique(LM_globalTargets$agegroup))
FAABFeatures = list("FirstAAb"=c(1,2,5))

###########################################################
# SUBSET HERE
dim(LM_globalTargets)# 742 x 20
LM_globalTargets[1:3,1:10]
LM_globalEXPMAT[1:3,1:10]
colnames(Allmetabolomicstogether[1:40])
Allmetabolomicstogether[1:10,1:10]
dim(Allmetabolomicstogether)
Allmetabolomicstogether[1:5,50:60]
rownames(Allmetabolomicstogether)<-Allmetabolomicstogether[,1]
metabolomicstargets<-Allmetabolomicstogether[,1:35]
dim(metabolomicstargets)
metabolomicstargets[1:10,1:10]
globalmetabolomicsMatrix<-t(Allmetabolomicstogether[,-(1:35)])
dim(globalmetabolomicsMatrix)
globalmetabolomicsMatrix[1:3,1:10]
globalmetabolomicsMatrix<-log2(globalmetabolomicsMatrix)
globalmetabolomicsMatrix[1:3,1:10]
namesrows<-gsub("._neglip", "!_neglip",rownames(globalmetabolomicsMatrix))
namesrows2<-gsub("._poslip", "!_poslip",namesrows)
rownames(globalmetabolomicsMatrix)<-namesrows2

###########################################################
# For AgeGroup:
test<-PaintomicsSubsetter(X=LM_globalTargets, Y=LM_globalEXPMAT,  Z= metabolomicstargets, ZZ=globalmetabolomicsMatrix,
                          ageg=0, FAAb=0,timepoints=3, timeframe=c(12,9,6,3,0) )

test$paintomicsGETable
test$paintomicsMETTable
test$targets
test$datasetGE
###########################################################
# Merging with AllmetabolitesConverted
test$paintomicsMETTable
head(AllmetabolitesConverted)
rownames(AllmetabolitesConverted)<-AllmetabolitesConverted[,1];AllmetabolitesConverted<-AllmetabolitesConverted[-1]
log2ratMetabolsm12tom0originalnames<-merge(AllmetabolitesConverted,test$paintomicsMETTable, by=0)
dim(AllmetabolitesConverted);dim(test$paintomicsMETTable);dim(log2ratMetabolsm12tom0originalnames);
head(log2ratMetabolsm12tom0originalnames)
#write.table(log2ratMetabolsm12tom0originalnames, "/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/Paintomics/log2ratMetabolsm12tom0originalnames.txt", sep="\t")

###########################################################
# Rename with FunCatsamename all dataset
#colnames(log2ratMetabolsm12tom0originalnames)
funcatsn<-log2ratMetabolsm12tom0originalnames[,c(5,7:11)]
head(funcatsn)

#  Make Unique based on the names
UniqueSamemetnamesm12tom0from1321<-aggregate(.~FunCatsameName,funcatsn,mean)
dim(UniqueSamemetnamesm12tom0from1321)
head(UniqueSamemetnamesm12tom0from1321)
rownames(UniqueSamemetnamesm12tom0from1321)<-UniqueSamemetnamesm12tom0from1321[,1]

UniqueSamemetnamesm12tom0from1321<-UniqueSamemetnamesm12tom0from1321[,-1]
#write.table(UniqueSamemetnamesm12tom0from1321, "/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/Paintomics/UniqueSamemetnamesm12tom0from1321.txt", sep="\t")
#a<-grep("rol", UniqueSamemetnamesm12tom0from1321[,1], value=TRUE)
#fila<-UniqueSamemetnamesm12tom0from1321[is.element(UniqueSamemetnamesm12tom0from1321[,1],a),]
#fila

###########################################################
#####       So for Paintomics we have to create      ######
###########################################################
# 1) PaintomicsGlobal
# GeneExpression:
# test$paintomicsGETable
# 1) PaintGlobalGenes<- test$paintomicsGETable
#write.table(PaintGlobalGenes, "/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/Paintomics/Stratification/PaintGlobalGenes.txt", sep="\t")

# 2) PaintDefinitiveGeneListNewData.txt # This is the list of 862 genes selected by NPLSDA

# Metabolomics:
# UniqueSamemetnamesm12tom0from1321
# 3) PaintGlobalMetabs<-UniqueSamemetnamesm12tom0from1321
#write.table(PaintGlobalMetabs, "/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/Paintomics/Stratification/PaintGlobalMetabs.txt", sep="\t")
# 4) PaintlistUniqueMetabolites.txt  # This is the list of Unique names of the selected 245 metabolites

# Use the four datasets as input in Paintomics