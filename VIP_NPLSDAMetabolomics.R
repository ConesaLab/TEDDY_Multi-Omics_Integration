###########################################################
############     VIP_NPLSDAMetabolomics.R     #############
###########################################################
# Author: Leandro Balzano-Nogueira
# Genetics Institute, University of Florida (Gainesville)

# This script is to create the Metabolomics  dataset to calculate the NPLSDA and VIP selection

###########################################################
homedir<- "/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/Metabolomics/" # Home directory where all your results are going to be contained
setwd(homedir)
getwd()
###########################################################
# Functions:
# "/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/ScriptsForNComm/Tools/" # Location of TEDDYtools
source ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/ScriptsForNComm/Tools/TEDDYtools2.R") # These are the functions to reformat the data
source ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/ScriptsForNComm/Tools/NPLSDAfunctionsApr11.R") # These are the functions created to perform the NPLSDA

###########################################################
# Data:
# Metabolomics RAW
GCTOFRaw<-read.csv ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/Metabolomics/Raw/GCTOFRaw.csv",header = TRUE)
GCTOFRaw[1:10,1:10]

NegLipRaw<-read.csv ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/Metabolomics/Raw/NegLipRaw.csv",header = TRUE)
NegLipRaw[1:10,1:10]

PosLipRaw<-read.csv ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/Metabolomics/Raw/PosLipRaw.csv",header = TRUE)
PosLipRaw[1:10,1:10]

# Metabolomics Processed
GCTOFProcessed<-read.csv ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/Metabolomics/Processed/GCTOF.csv",header = TRUE)
GCTOFProcessed[1:10,1:10]
 
NegLipProcessed<-read.csv ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/Metabolomics/Processed/NegLip.csv",header = TRUE)
NegLipProcessed[1:10,1:10]
 
PosLipProcessed<-read.csv ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/Metabolomics/Processed/PosLip.csv",header = TRUE)
PosLipProcessed[1:10,1:10]


# Response Variable
CohortData<-read.csv ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/CohortData.csv",header = TRUE)
CohortData[1:10,]

# List of Cases with at least 3 out of 5 time points with data
patients3tps<-data.frame(V1=CohortData$Individual.Id[CohortData$Model.or.Validation=="Model"])
patients3tps

###########################################################
# Libraries:
library("gdata")
require("VennDiagram")

###########################################################
# GCTOFRaw
GCTOFv2<-t(GCTOFRaw[-1,])
GCTOFv2[1:5,1:5]
colnames(GCTOFv2)<- GCTOFv2[1,]; GCTOFv2<-GCTOFv2[-1,]

GCTOFv2 <- transform(GCTOFv2, Individual.Id =as.numeric(as.character(Individual.Id)),
                           Age.in.Months = as.numeric(as.character(Age.in.Months)), 
                           Time.to.IA = as.numeric(as.character(Time.to.IA)))

GCTOFv2[1:10,1:10]

# Subsetting the data to the 136 pairs with values reported in 3 0ut 0f 5 time points
GCTOFv2136<- GCTOFv2[is.element(GCTOFv2$Individual.Id, patients3tps[,1]),]
dim(GCTOFv2136)
length(unique(GCTOFv2136$Individual.Id))  # 136 PERFECT


# NegLipRaw
NegLipRawv2<-t(NegLipRaw)
NegLipRawv2[1:5,1:5]
colnames(NegLipRawv2)<- NegLipRawv2[1,]; NegLipRawv2<-NegLipRawv2[-1,]

NegLipRawv2 <- transform(NegLipRawv2, Individual.Id =as.numeric(as.character(Individual.Id)),
                     Age.in.Months = as.numeric(as.character(Age.in.Months)), 
                     Time.to.IA = as.numeric(as.character(Time.to.IA)))

NegLipRawv2[1:10,1:10]

# Subsetting the data to the 136 pairs with values reported in 3 0ut 0f 5 time points
NegLipRawv2136<- NegLipRawv2[is.element(NegLipRawv2$Individual.Id, patients3tps[,1]),]
dim(NegLipRawv2136)
length(unique(NegLipRawv2136$Individual.Id))  # 136 PERFECT


# PosLipRaw
PosLipRawv2<-t(PosLipRaw)
PosLipRawv2[1:5,1:5]
colnames(PosLipRawv2)<- PosLipRawv2[1,]; PosLipRawv2<-PosLipRawv2[-1,]

PosLipRawv2 <- transform(PosLipRawv2, Individual.Id =as.numeric(as.character(Individual.Id)),
                         Age.in.Months = as.numeric(as.character(Age.in.Months)), 
                         Time.to.IA = as.numeric(as.character(Time.to.IA)))

PosLipRawv2[1:10,1:10]

# Subsetting the data to the 136 pairs with values reported in 3 0ut 0f 5 time points
PosLipRawv2136<- PosLipRawv2[is.element(PosLipRawv2$Individual.Id, patients3tps[,1]),]
dim(PosLipRawv2136)
length(unique(PosLipRawv2136$Individual.Id))  # 136 PERFECT

###########################################################
# Data in 3D structure
# GCTOF
dim(GCTOFv2136)   # 570 x 367
colnames(GCTOFv2136[1:23])
GCTOFX<-GCTOFv2136[,c(1,4:367)]       # Just ID and the variables


GCTOFX12<-GCTOFX[GCTOFv2136$Time.to.IA == "-12",]; dim(GCTOFX12);dim(GCTOFX)
GCTOFX9<-GCTOFX[GCTOFv2136$Time.to.IA == "-9",];dim(GCTOFX9)
GCTOFX6<-GCTOFX[GCTOFv2136$Time.to.IA == "-6",];dim(GCTOFX6)
GCTOFX3<-GCTOFX[GCTOFv2136$Time.to.IA == "-3",];dim(GCTOFX3)
GCTOFX0<-GCTOFX[GCTOFv2136$Time.to.IA == "0",];dim(GCTOFX0)

dim(GCTOFX12);dim(GCTOFX9);dim(GCTOFX6);dim(GCTOFX3);dim(GCTOFX0)


# NegLip
dim(NegLipRawv2136)   # 570 x 446
colnames(NegLipRawv2136[1:23])
NegLipRawX<-NegLipRawv2136[,c(1,4:446)]       # Just ID and the variables


NegLipX12<-NegLipRawX[NegLipRawv2136$Time.to.IA == "-12",]; dim(NegLipX12);dim(NegLipRawX)
NegLipX9<-NegLipRawX[NegLipRawv2136$Time.to.IA == "-9",];dim(NegLipX9)
NegLipX6<-NegLipRawX[NegLipRawv2136$Time.to.IA == "-6",];dim(NegLipX6)
NegLipX3<-NegLipRawX[NegLipRawv2136$Time.to.IA == "-3",];dim(NegLipX3)
NegLipX0<-NegLipRawX[NegLipRawv2136$Time.to.IA == "0",];dim(NegLipX0)

dim(NegLipX12);dim(NegLipX9);dim(NegLipX6);dim(NegLipX3);dim(NegLipX0)


# PosLip
dim(PosLipRawv2136)   # 570 x 446
colnames(PosLipRawv2136[1:23])
PosLipRawX<-PosLipRawv2136[,c(1,4:517)]       # Just ID and the variables


PosLipX12<-NegLipRawX[PosLipRawv2136$Time.to.IA == "-12",]; dim(PosLipX12);dim(PosLipRawX)
PosLipX9<-NegLipRawX[PosLipRawv2136$Time.to.IA == "-9",];dim(PosLipX9)
PosLipX6<-NegLipRawX[PosLipRawv2136$Time.to.IA == "-6",];dim(PosLipX6)
PosLipX3<-NegLipRawX[PosLipRawv2136$Time.to.IA == "-3",];dim(PosLipX3)
PosLipX0<-NegLipRawX[PosLipRawv2136$Time.to.IA == "0",];dim(PosLipX0)

dim(PosLipX12);dim(PosLipX9);dim(PosLipX6);dim(PosLipX3);dim(PosLipX0)

###########################################################
# Merging with all cases
# GCTOF
patients3tps2<-patients3tps
colnames(patients3tps2)<- "Individual.Id"
GCTOFX12total<-merge(patients3tps2,GCTOFX12, by="Individual.Id", all.x = T);dim(GCTOFX12total);dim(patients3tps2);dim(GCTOFX12)
GCTOFX12total[1:5,1:5]
rownames(GCTOFX12total)<- GCTOFX12total[,1]
GCTOFX12total<- as.matrix(GCTOFX12total[,c(-1)])
dim(GCTOFX12total)

GCTOFX9total<-merge(patients3tps2,GCTOFX9, by="Individual.Id", all.x = T);dim(GCTOFX9total);dim(patients3tps2);dim(GCTOFX9)
GCTOFX9total[1:5,1:5]
rownames(GCTOFX9total)<- GCTOFX9total[,1]
GCTOFX9total<- as.matrix(GCTOFX9total[,c(-1)])
dim(GCTOFX9total)

GCTOFX6total<-merge(patients3tps2,GCTOFX6, by="Individual.Id", all.x = T);dim(GCTOFX6total);dim(patients3tps2);dim(GCTOFX6total)
GCTOFX6total[1:5,1:5]
rownames(GCTOFX6total)<- GCTOFX6total[,1]
GCTOFX6total<- as.matrix(GCTOFX6total[,c(-1)])
dim(GCTOFX6total)

GCTOFX3total<-merge(patients3tps2,GCTOFX3, by="Individual.Id", all.x = T);dim(GCTOFX3total);dim(patients3tps2);dim(GCTOFX3)
GCTOFX3total[1:5,1:5]
rownames(GCTOFX3total)<- GCTOFX3total[,1]
GCTOFX3total<- as.matrix(GCTOFX3total[,c(-1)])
dim(GCTOFX3total)

GCTOFX0total<-merge(patients3tps2,GCTOFX0, by="Individual.Id", all.x = T);dim(GCTOFX0total);dim(patients3tps2);dim(GCTOFX0)
GCTOFX0total[1:5,1:5]
rownames(GCTOFX0total)<- GCTOFX0total[,1]
GCTOFX0total<- as.matrix(GCTOFX0total[,c(-1)])
dim(GCTOFX0total)
# Dimensions are 136*364*5


#NegLip
patients3tps2<-patients3tps
colnames(patients3tps2)<- "Individual.Id"
NegLipX12total<-merge(patients3tps2,NegLipX12, by="Individual.Id", all.x = T);dim(NegLipX12total);dim(patients3tps2);dim(NegLipX12)
NegLipX12total[1:5,1:5]
rownames(NegLipX12total)<- NegLipX12total[,1]
NegLipX12total<- as.matrix(NegLipX12total[,c(-1)])
dim(NegLipX12total)

NegLipX9total<-merge(patients3tps2,NegLipX9, by="Individual.Id", all.x = T);dim(NegLipX9total);dim(patients3tps2);dim(NegLipX9)
NegLipX9total[1:5,1:5]
rownames(NegLipX9total)<- NegLipX9total[,1]
NegLipX9total<- as.matrix(NegLipX9total[,c(-1)])
dim(NegLipX9total)

NegLipX6total<-merge(patients3tps2,NegLipX6, by="Individual.Id", all.x = T);dim(NegLipX6total);dim(patients3tps2);dim(NegLipX6)
NegLipX6total[1:5,1:5]
rownames(NegLipX6total)<- NegLipX6total[,1]
NegLipX6total<- as.matrix(NegLipX6total[,c(-1)])
dim(NegLipX6total)

NegLipX3total<-merge(patients3tps2,NegLipX3, by="Individual.Id", all.x = T);dim(NegLipX3total);dim(patients3tps2);dim(NegLipX3)
NegLipX3total[1:5,1:5]
rownames(NegLipX3total)<- NegLipX3total[,1]
NegLipX3total<- as.matrix(NegLipX3total[,c(-1)])
dim(NegLipX3total)

NegLipX0total<-merge(patients3tps2,NegLipX0, by="Individual.Id", all.x = T);dim(NegLipX0total);dim(patients3tps2);dim(NegLipX0)
NegLipX0total[1:5,1:5]
rownames(NegLipX0total)<- NegLipX0total[,1]
NegLipX0total<- as.matrix(NegLipX0total[,c(-1)])
dim(NegLipX0total)
# Dimensions are 136*443*5


#PosLip
patients3tps2<-patients3tps
colnames(patients3tps2)<- "Individual.Id"
PosLipX12total<-merge(patients3tps2,PosLipX12, by="Individual.Id", all.x = T);dim(PosLipX12total);dim(patients3tps2);dim(PosLipX12)
PosLipX12total[1:5,1:5]
rownames(PosLipX12total)<- PosLipX12total[,1]
PosLipX12total<- as.matrix(PosLipX12total[,c(-1)])
dim(PosLipX12total)

PosLipX9total<-merge(patients3tps2,PosLipX9, by="Individual.Id", all.x = T);dim(PosLipX9total);dim(patients3tps2);dim(PosLipX9)
PosLipX9total[1:5,1:5]
rownames(PosLipX9total)<- PosLipX9total[,1]
PosLipX9total<- as.matrix(PosLipX9total[,c(-1)])
dim(PosLipX9total)

PosLipX6total<-merge(patients3tps2,PosLipX6, by="Individual.Id", all.x = T);dim(PosLipX6total);dim(patients3tps2);dim(PosLipX6)
PosLipX6total[1:5,1:5]
rownames(PosLipX6total)<- PosLipX6total[,1]
PosLipX6total<- as.matrix(PosLipX6total[,c(-1)])
dim(PosLipX6total)

PosLipX3total<-merge(patients3tps2,PosLipX3, by="Individual.Id", all.x = T);dim(PosLipX3total);dim(patients3tps2);dim(PosLipX3)
PosLipX3total[1:5,1:5]
rownames(PosLipX3total)<- PosLipX3total[,1]
PosLipX3total<- as.matrix(PosLipX3total[,c(-1)])
dim(PosLipX3total)

PosLipX0total<-merge(patients3tps2,PosLipX0, by="Individual.Id", all.x = T);dim(PosLipX0total);dim(patients3tps2);dim(PosLipX0)
PosLipX0total[1:5,1:5]
rownames(PosLipX0total)<- PosLipX0total[,1]
PosLipX0total<- as.matrix(PosLipX0total[,c(-1)])
dim(PosLipX0total)
# Dimensions are 136*443*5

###########################################################
# Array structure GCTOF
dim(GCTOFX12total)
arrayGCTOFX <- array(data = NA, dim = c(136,364,5),dimnames = list(NULL, NULL, c("-12","-9","-6", "-3", "0")))

arrayGCTOFX
arrayGCTOFX[,,1] <- GCTOFX12total
arrayGCTOFX[,,2] <- GCTOFX9total
arrayGCTOFX[,,3] <- GCTOFX6total
arrayGCTOFX[,,4] <- GCTOFX3total
arrayGCTOFX[,,5] <- GCTOFX0total

rownames(arrayGCTOFX)<-rownames(GCTOFX12total)
colnames(arrayGCTOFX)<-colnames (GCTOFX12total)
arrayGCTOFX[1:50,1:5,1]
dim(arrayGCTOFX)   # 136 * 364 * 5


# Array structure NegLip
dim(NegLipX12total)
arrayNegLipX <- array(data = NA, dim = c(136,443,5),dimnames = list(NULL, NULL, c("-12","-9","-6", "-3", "0")))

arrayNegLipX
arrayNegLipX[,,1] <- NegLipX12total
arrayNegLipX[,,2] <- NegLipX9total
arrayNegLipX[,,3] <- NegLipX6total
arrayNegLipX[,,4] <- NegLipX3total
arrayNegLipX[,,5] <- NegLipX0total

rownames(arrayNegLipX)<-rownames(NegLipX12total)
colnames(arrayNegLipX)<-colnames (NegLipX12total)
arrayNegLipX[1:50,1:5,1]
dim(arrayNegLipX)   # 136 * 443 * 5


# Array structure PosLip
dim(PosLipX12total)
arrayPosLipX <- array(data = NA, dim = c(136,443,5),dimnames = list(NULL, NULL, c("-12","-9","-6", "-3", "0")))

arrayPosLipX
arrayPosLipX[,,1] <- PosLipX12total
arrayPosLipX[,,2] <- PosLipX9total
arrayPosLipX[,,3] <- PosLipX6total
arrayPosLipX[,,4] <- PosLipX3total
arrayPosLipX[,,5] <- PosLipX0total

rownames(arrayPosLipX)<-rownames(PosLipX12total)
colnames(arrayPosLipX)<-colnames (PosLipX12total)
arrayPosLipX[1:50,1:5,1]
dim(arrayPosLipX)   # 136 * 443 * 5

###########################################################
# Imputation (For convenience it must be done in an HPC)
# 1) Determining the best fitted model

modelGCTOFX136<-bestfittedmodel (X=arrayGCTOFX,centering=0) # 0= No centering; 1= centering by Individuals; 2= centering by Variables;3= centering by Time
# The best model was 4,2,4 

modelNegLipX136<-bestfittedmodel (X=arrayNegLipX,centering=0) # 0= No centering; 1= centering by Individuals; 2= centering by Variables;3= centering by Time
# The best model was 4,2,4 

modelPosLipXnuevo<-bestfittedmodel (X=arrayPosLipX,centering=0) # 0= No centering; 1= centering by Individuals; 2= centering by Variables;3= centering by Time
# The best model was 4,4,4 


# 2) Imputing the best fitted model data
# GCTOF
fullarrayGCTOF136<-Imputemethod(X=arrayGCTOFX, fac=c(4, 2, 4), conver = 1e-07, max.iter = 1000)

# NegLip
fullarrayNegLip136<-Imputemethod(X=arrayNegLipX, fac=c(4, 2, 4), conver = 1e-07, max.iter = 1000)

# PosLip
fullarrayPosLip136<-Imputemethod(X=arrayPosLipX, fac=c(4, 4, 4), conver = 1e-07, max.iter = 1000)

###########################################################
#   GCTOF:
NPLSDA136GCTOF<-NPLSDAmod (XN=fullarrayGCTOF136, YN=outcomedummyarray136,factors = 2,COMP= c(4,2,4), conver = 1e-16, max.iteration = 10000,
                           centering=0)

ploteoNPLSDAGCTOFtotal136indvs<- plotNPLSDAmod (X=NPLSDA136GCTOF, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                                cutoff = 20, factors=2, penalty=1) 


###########################################################
#   NegLip:
NPLSDA136NegLip<-NPLSDAmod (XN=fullarrayNegLip136, YN=outcomedummyarray136,factors = 2,COMP= c(4,2,4), conver = 1e-16, max.iteration = 10000,
                            centering=0)

ploteoNPLSDANegLiptotal136indvs<- plotNPLSDAmod (X=NPLSDA136NegLip, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                                 cutoff = 20, factors=2, penalty=1) 

###########################################################
#   PosLip:
NPLSDA136PosLip<-NPLSDAmod (XN=fullarrayPosLip136, YN=outcomedummyarray136,factors = 2,COMP= c(4,2,4), conver = 1e-16, max.iteration = 10000,
                            centering=0)

ploteoNPLSDAPosLiptotal136indvs<- plotNPLSDAmod (X=NPLSDA136PosLip, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                                 cutoff = 20, factors=2, penalty=1) 

###########################################################
############    Variable Selection GCTOF     ##############
###########################################################
# Variable selection by VIP3Dmodel2
summary(NPLSDA136GCTOF)
NPLSDA136GCTOF$VIP3Dmodel2

### Desde aca
vipsoutcomemet<-data.frame(NPLSDA136GCTOF$VIP3Dmodel2)

thrs<-95
vipsoutcomemetsubset95<-NULL
temp<-NULL
vipsoutcomemetsubset<-list()
for (i in 1:length(vipsoutcomemet)){
  a<-quantile(vipsoutcomemet[,i], probs = c(50, 95,99, thrs)/100)
  print(i)
  print(a)
  temp<-vipsoutcomemet[vipsoutcomemet[,i]>=a[4],]
  namesubset<-paste('t',i,sep='')
  vips<-list(temp)
  vipsoutcomemetsubset[[namesubset]]<-vips
}
vipsoutcomemetsubset
rm(vaya)
#rownames(vipsoutcomemetsubset$t1[[1]])
vaya<-c (as.vector(rownames(vipsoutcomemetsubset$t1[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t2[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t3[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t4[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t5[[1]]))) 
duplicated(vaya)
PosLipselVars<-vaya[!duplicated(vaya)]
length(PosLipselVars)
######################################################################################################
# Retain just these variables in GCTOF 
dim(fullarrayGCTOF136)
PosLipselVars
#selectedGCTOFbydifVIPs<-list(VIP3Dmodel2.95p.85v=PosLipselVars)

#summary(selectedGCTOFbydifVIPs)
#selectedgenesbydifVIPs<-c(selectedgenesbydifVIPs,tulo="pelotas")

FullarrayGenesVIPSelVars<-fullarrayGCTOF136[,is.element(colnames(fullarrayGCTOF136),
                                                        PosLipselVars),]
dim(FullarrayGenesVIPSelVars)

#FullarrayGenesVIPSelVarsVIP3D.2<-FullarrayGenesVIPSelVars
#save(FullarrayGenesVIPSelVarsVIP3D.2,file="FullarrayGenesVIPSelVarsVIP3D.2.RData")

### Now NPLSDA and Graph
NPLSDAarraygenesVIPselected<-NPLSDAmod(XN=FullarrayGenesVIPSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

#NPLSDAarraygenesVIPselected$NPLSDAvariatesperMode3
#NPLSDAarraygenesVIPselected$NPLSDAvariates$X
#dim(NPLSDAarraygenesVIPselected$NPLSDAvariatesperMode3)
# Plotting
#pdf("VIP3Dmodel2.99p.1013g.pdf")
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 3), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=2, penalty=2) 
#dev.off()

ploteoNPLSDAVIPselectedgenes

######################################################################################################################
# Variable selection by VIP3Dmodel1

summary(NPLSDA136GCTOF)
NPLSDA136GCTOF$VIP3Dmodel1
dim(NPLSDA136GCTOF$VIP3Dmodel1)



### Desde aca
vipsoutcomemet<-data.frame(NPLSDA136GCTOF$VIP3Dmodel1)

thrs<-95
vipsoutcomemetsubset95<-NULL
temp<-NULL
vipsoutcomemetsubset<-list()
for (i in 1:2){
  a<-quantile(vipsoutcomemet[,i], probs = c(50, 95,99, thrs)/100)
  print(i)
  print(a)
  temp<-vipsoutcomemet[vipsoutcomemet[,i]>=a[4],]
  namesubset<-paste('t',i,sep='')
  vips<-list(temp)
  vipsoutcomemetsubset[[namesubset]]<-vips
}
vipsoutcomemetsubset
rm(vaya)
#rownames(vipsoutcomemetsubset$t1[[1]])
vaya<-c (as.vector(rownames(vipsoutcomemetsubset$t1[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t2[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t3[[1]]))) 
duplicated(vaya)
PosLipselVars<-vaya[!duplicated(vaya)]
length(PosLipselVars)  # 2050
######################################################################################################
# Retain just these variables in Gene Expression 2D
dim(fullarrayGCTOF136)
PosLipselVars

FullarrayGenesVIPSelVars<-fullarrayGCTOF136[,is.element(colnames(fullarrayGCTOF136),
                                                        PosLipselVars),]
dim(FullarrayGenesVIPSelVars)
#FullarrayGenesVIPSelVarsVIP1.comp1<-FullarrayGenesVIPSelVars
#FullarrayGenesVIPSelVarsVIP1.comps12<-FullarrayGenesVIPSelVars
#save(FullarrayGenesVIPSelVarsVIP1.comp1,file="FullarrayGenesVIPSelVarsVIP1.comp1.RData")
#save(FullarrayGenesVIPSelVarsVIP1.comps12,file="FullarrayGenesVIPSelVarsVIP1.comps12.RData")


### Now NPLSDA and Graph
NPLSDAarraygenesVIPselected<-NPLSDAmod(XN=FullarrayGenesVIPSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=3, centering=0) 

# Plotting
#pdf("VIP3Dmodel1.99p.1comp.213g.pdf")
#pdf("VIP3Dmodel1.99p.12comp.336g.pdf")
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=3, penalty=2) 
#dev.off()


#summary(selectedgenesbydifVIPs)
#selectedgenesbydifVIPs<-append(selectedgenesbydifVIPs,list(VIP3Dmodel1.99p.1comp.213g=PosLipselVars))
#selectedgenesbydifVIPs<-append(selectedgenesbydifVIPs,list(VIP3Dmodel1.99p.12comp.336g=PosLipselVars))


######################################################################################################
# Gene Expression USING VIP2D

summary(NPLSDA136GCTOF)
NPLSDA136GCTOF$VIP2D


### Desde aca
vipsoutcome2D<-data.frame(NPLSDA136GCTOF$VIP2D)

colp1<-paste(rownames(NPLSDA136GCTOF$VIP3Dmodel1),"-12", sep="!_____")
colp2<-paste(rownames(NPLSDA136GCTOF$VIP3Dmodel1),"-9", sep="!_____")
colp3<-paste(rownames(NPLSDA136GCTOF$VIP3Dmodel1),"-6", sep="!_____")
colp4<-paste(rownames(NPLSDA136GCTOF$VIP3Dmodel1),"-3", sep="!_____")
colp5<-paste(rownames(NPLSDA136GCTOF$VIP3Dmodel1),"0", sep="!_____")


rownames(vipsoutcome2D)<-c(colp1,colp2,colp3,colp4,colp5)

#apply(vipsoutcomemet, 2, function(x) is.numeric(x))
vipsoutcomemet<-vipsoutcome2D
thrs<-95
vipsoutcomemetsubset95<-NULL
temp<-NULL
vipsoutcomemetsubset<-list()
for (i in 1){
  a<-quantile(vipsoutcomemet[,i], probs = c(50, 95,99, thrs)/100)
  print(i)
  print(a)
  temp<-vipsoutcomemet[vipsoutcomemet[,i]>=a[4],]
  namesubset<-paste('t',i,sep='')
  vips<-list(temp)
  vipsoutcomemetsubset[[namesubset]]<-vips
}
vipsoutcomemetsubset
rm(vaya)
#rownames(vipsoutcomemetsubset$t1[[1]])
vaya<-c (as.vector(rownames(vipsoutcomemetsubset$t1[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t2[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t3[[1]]))) 
duplicated(vaya)
PosLipselVars<-vaya[!duplicated(vaya)]
length(PosLipselVars)
string<-strsplit (PosLipselVars,"!_____")
string<-listTOdata.frame(string)
genelist<- string[,1]
duplicated (genelist)
genes<-genelist[!duplicated(genelist)]
length(genes)

######################################################################################################
# Retain just these variables in Gene Expression 2D
dim(fullarrayGCTOF136)
genes
#as.vector(genes)
FullarrayGenesVIPSelVars<-fullarrayGCTOF136[,is.element(colnames(fullarrayGCTOF136),
                                                        genes),]
dim(FullarrayGenesVIPSelVars)
#FullarrayGenesVIPSelVarsVIP2Dcomp2<-FullarrayGenesVIPSelVars
#save(FullarrayGenesVIPSelVarsVIP2Dcomp2,file="FullarrayGenesVIPSelVarsVIP2Dcomp2.RData")

### Now NPLSDA and Graph
NPLSDAarraygenesVIPselected<-NPLSDAmod(XN=FullarrayGenesVIPSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf ("VIP2D.99p.comp2.1022g.pdf")
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=2, penalty=2) 
#dev.off()

#summary(selectedGCTOFbydifVIPs)
#selectedGCTOFbydifVIPs<-append(selectedGCTOFbydifVIPs,list(VIP2D.95p.comp12.82v= as.vector(genes)))
#save(selectedGCTOFbydifVIPs,file="selectedGCTOFbydifVIPs.RData")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/selectedGCTOFbydifVIPs.RData")
###########################################################################################################
# Venn Diagram of the four genelists #

summary(selectedGCTOFbydifVIPs)

# If this is a bit confusing you can also write a function and then use it in lapply() 
removeEMPTYstrings <- function(x) {
  
  newVectorWOstrings <- x[x != ""]
  return(newVectorWOstrings)
  
}
geneLS2 <- lapply(selectedGCTOFbydifVIPs, removeEMPTYstrings)

summary(geneLS2)



lapply(geneLS2, tail) # Both methods return the same results

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:
require("VennDiagram")

VENN.LIST <- geneLS2
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkmagenta", "darkblue"),#,"darkseagreen4","goldenrod2"), 
                          alpha=c(0.5,0.5), cex = 2, cat.fontface=4, 
                          category.names=c("VIP3D.2", "VIP2D.comp12"), main="NPLSDAGCTOFLists",
                          cex.main=2
)

# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
grid.draw(venn.plot)

#pdf("VenndiagramsselectedGCTOF.pdf")
grid.draw(venn.plot)
#dev.off()
####################################################################################################################
# Test the Union and the intersection between the two metabolites sets

theunion<- union(selectedGCTOFbydifVIPs[[1]],selectedGCTOFbydifVIPs[[2]])
length(theunion) #91
theintersection<- intersect(selectedGCTOFbydifVIPs[[1]],selectedGCTOFbydifVIPs[[2]])
length(theintersection) # 76

####
FullarrayUNIONGCTOFSelVars<-fullarrayGCTOF136[,is.element(colnames(fullarrayGCTOF136),
                                                          theunion),]
dim(FullarrayUNIONGCTOFSelVars)

### Now NPLSDA and Graph
NPLSDAUnionGCTOF<-NPLSDAmod(XN=FullarrayUNIONGCTOFSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf ("UNIONGCTOF.pdf")
ploteoNPLSDAUnionGCTOF<- plotNPLSDAmod (X=NPLSDAUnionGCTOF, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                        cutoff = 20, factors=2, penalty=2) 
#dev.off()

################################################################################################################
FullarrayINTERSECTGCTOFSelVars<-fullarrayGCTOF136[,is.element(colnames(fullarrayGCTOF136),
                                                              theintersection),]
dim(FullarrayINTERSECTGCTOFSelVars)

### Now NPLSDA and Graph
NPLSDAIntersectionGCTOF<-NPLSDAmod(XN=FullarrayINTERSECTGCTOFSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf ("IntersectionGCTOF.pdf")
ploteoNPLSDAIntersectionGCTOF<- plotNPLSDAmod (X=NPLSDAIntersectionGCTOF, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                               cutoff = 20, factors=2, penalty=2) 
#dev.off()


#############################################################################################################
###################################    Variable Selection NegLip     ########################################
#############################################################################################################
# Variable selection by VIP3Dmodel2
summary(NPLSDA136NegLip)
NPLSDA136NegLip$VIP3Dmodel2

### Desde aca
vipsoutcomemet<-data.frame(NPLSDA136NegLip$VIP3Dmodel2)

thrs<-95
vipsoutcomemetsubset95<-NULL
temp<-NULL
vipsoutcomemetsubset<-list()
for (i in 1:length(vipsoutcomemet)){
  a<-quantile(vipsoutcomemet[,i], probs = c(50, 95,99, thrs)/100)
  print(i)
  print(a)
  temp<-vipsoutcomemet[vipsoutcomemet[,i]>=a[4],]
  namesubset<-paste('t',i,sep='')
  vips<-list(temp)
  vipsoutcomemetsubset[[namesubset]]<-vips
}
vipsoutcomemetsubset
rm(vaya)
#rownames(vipsoutcomemetsubset$t1[[1]])
vaya<-c (as.vector(rownames(vipsoutcomemetsubset$t1[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t2[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t3[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t4[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t5[[1]]))) 
duplicated(vaya)
PosLipselVars<-vaya[!duplicated(vaya)]
length(PosLipselVars)
######################################################################################################
# Retain just these variables in GCTOF 
dim(fullarrayNegLip136)
PosLipselVars
#selectedNegLipbydifVIPs<-list(VIP3Dmodel2.95p.98v=PosLipselVars)

#summary(selectedGCTOFbydifVIPs)
#selectedgenesbydifVIPs<-c(selectedgenesbydifVIPs,tulo="pelotas")

FullarrayGenesVIPSelVars<-fullarrayNegLip136[,is.element(colnames(fullarrayNegLip136),
                                                         PosLipselVars),]
dim(FullarrayGenesVIPSelVars)

#FullarrayGenesVIPSelVarsVIP3D.2<-FullarrayGenesVIPSelVars
#save(FullarrayGenesVIPSelVarsVIP3D.2,file="FullarrayGenesVIPSelVarsVIP3D.2.RData")

### Now NPLSDA and Graph
NPLSDAarraygenesVIPselected<-NPLSDAmod(XN=FullarrayGenesVIPSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

#NPLSDAarraygenesVIPselected$NPLSDAvariatesperMode3
#NPLSDAarraygenesVIPselected$NPLSDAvariates$X
#dim(NPLSDAarraygenesVIPselected$NPLSDAvariatesperMode3)
# Plotting
#pdf("VIP3Dmodel2.99p.1013g.pdf")
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=2, penalty=2) 
#dev.off()

ploteoNPLSDAVIPselectedgenes

######################################################################################################################
# Variable selection by VIP3Dmodel1

summary(NPLSDA136NegLip)
NPLSDA136NegLip$VIP3Dmodel1
dim(NPLSDA136NegLip$VIP3Dmodel1)



### Desde aca
vipsoutcomemet<-data.frame(NPLSDA136NegLip$VIP3Dmodel1)

thrs<-95
vipsoutcomemetsubset95<-NULL
temp<-NULL
vipsoutcomemetsubset<-list()
for (i in 1:2){
  a<-quantile(vipsoutcomemet[,i], probs = c(50, 95,99, thrs)/100)
  print(i)
  print(a)
  temp<-vipsoutcomemet[vipsoutcomemet[,i]>=a[4],]
  namesubset<-paste('t',i,sep='')
  vips<-list(temp)
  vipsoutcomemetsubset[[namesubset]]<-vips
}
vipsoutcomemetsubset
rm(vaya)
#rownames(vipsoutcomemetsubset$t1[[1]])
vaya<-c (as.vector(rownames(vipsoutcomemetsubset$t1[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t2[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t3[[1]]))) 
duplicated(vaya)
PosLipselVars<-vaya[!duplicated(vaya)]
length(PosLipselVars)  # 2050
######################################################################################################
# Retain just these variables in Gene Expression 2D
dim(fullarrayNegLip136)
PosLipselVars

FullarrayGenesVIPSelVars<-fullarrayNegLip136[,is.element(colnames(fullarrayNegLip136),
                                                         PosLipselVars),]
dim(FullarrayGenesVIPSelVars)
#FullarrayGenesVIPSelVarsVIP1.comp1<-FullarrayGenesVIPSelVars
#FullarrayGenesVIPSelVarsVIP1.comps12<-FullarrayGenesVIPSelVars
#save(FullarrayGenesVIPSelVarsVIP1.comp1,file="FullarrayGenesVIPSelVarsVIP1.comp1.RData")
#save(FullarrayGenesVIPSelVarsVIP1.comps12,file="FullarrayGenesVIPSelVarsVIP1.comps12.RData")


### Now NPLSDA and Graph
NPLSDAarraygenesVIPselected<-NPLSDAmod(XN=FullarrayGenesVIPSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf("VIP3Dmodel1.99p.1comp.213g.pdf")
#pdf("VIP3Dmodel1.99p.12comp.336g.pdf")
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=2, penalty=2) 
#dev.off()


summary(selectedgenesbydifVIPs)
#selectedgenesbydifVIPs<-append(selectedgenesbydifVIPs,list(VIP3Dmodel1.99p.1comp.213g=PosLipselVars))
selectedgenesbydifVIPs<-append(selectedgenesbydifVIPs,list(VIP3Dmodel1.99p.12comp.336g=PosLipselVars))


######################################################################################################
# Gene Expression USING VIP2D

summary(NPLSDA136NegLip)
NPLSDA136NegLip$VIP2D


### Desde aca
vipsoutcome2D<-data.frame(NPLSDA136NegLip$VIP2D)

colp1<-paste(rownames(NPLSDA136NegLip$VIP3Dmodel1),"-12", sep="!_____")
colp2<-paste(rownames(NPLSDA136NegLip$VIP3Dmodel1),"-9", sep="!_____")
colp3<-paste(rownames(NPLSDA136NegLip$VIP3Dmodel1),"-6", sep="!_____")
colp4<-paste(rownames(NPLSDA136NegLip$VIP3Dmodel1),"-3", sep="!_____")
colp5<-paste(rownames(NPLSDA136NegLip$VIP3Dmodel1),"0", sep="!_____")


rownames(vipsoutcome2D)<-c(colp1,colp2,colp3,colp4,colp5)

#apply(vipsoutcomemet, 2, function(x) is.numeric(x))
vipsoutcomemet<-vipsoutcome2D
thrs<-95
vipsoutcomemetsubset95<-NULL
temp<-NULL
vipsoutcomemetsubset<-list()
for (i in 1:2){
  a<-quantile(vipsoutcomemet[,i], probs = c(50, 95,99, thrs)/100)
  print(i)
  print(a)
  temp<-vipsoutcomemet[vipsoutcomemet[,i]>=a[4],]
  namesubset<-paste('t',i,sep='')
  vips<-list(temp)
  vipsoutcomemetsubset[[namesubset]]<-vips
}
vipsoutcomemetsubset
rm(vaya)
#rownames(vipsoutcomemetsubset$t1[[1]])
vaya<-c (as.vector(rownames(vipsoutcomemetsubset$t1[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t2[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t3[[1]]))) 
duplicated(vaya)
PosLipselVars<-vaya[!duplicated(vaya)]
length(PosLipselVars)
string<-strsplit (PosLipselVars,"!_____")
string<-listTOdata.frame(string)
genelist<- string[,1]
duplicated (genelist)
genes<-genelist[!duplicated(genelist)]
length(genes)

######################################################################################################
# Retain just these variables in Gene Expression 2D
dim(fullarrayNegLip136)
genes
#as.vector(genes)
FullarrayGenesVIPSelVars<-fullarrayNegLip136[,is.element(colnames(fullarrayNegLip136),
                                                         genes),]
dim(FullarrayGenesVIPSelVars)
#FullarrayGenesVIPSelVarsVIP2Dcomp2<-FullarrayGenesVIPSelVars
#save(FullarrayGenesVIPSelVarsVIP2Dcomp2,file="FullarrayGenesVIPSelVarsVIP2Dcomp2.RData")

### Now NPLSDA and Graph
NPLSDAarraygenesVIPselected<-NPLSDAmod(XN=FullarrayGenesVIPSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf ("VIP2D.99p.comp2.1022g.pdf")
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=2, penalty=2) 
#dev.off()

#summary(selectedNegLipbydifVIPs)
#selectedNegLipbydifVIPs<-append(selectedNegLipbydifVIPs,list(VIP2D.95p.comp12.153v= as.vector(genes)))
#save(selectedNegLipbydifVIPs,file="selectedNegLipbydifVIPs.RData")
###########################################################################################################
# We selected VIP3D/2 95% and VIP2D 95% Comps 1&2
load("/media/data/leobalzano/ScriptsForTEDDY/Data/selectedNegLipbydifVIPs.RData")

# Venn Diagram of the four genelists #
summary(selectedNegLipbydifVIPs)

# If this is a bit confusing you can also write a function and then use it in lapply() 
removeEMPTYstrings <- function(x) {
  
  newVectorWOstrings <- x[x != ""]
  return(newVectorWOstrings)
  
}
geneLS2 <- lapply(selectedNegLipbydifVIPs, removeEMPTYstrings)

summary(geneLS2)



lapply(geneLS2, tail) # Both methods return the same results

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:
require("VennDiagram")

VENN.LIST <- geneLS2
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkmagenta", "darkblue"),#,"darkseagreen4","goldenrod2"), 
                          alpha=c(0.5,0.5), cex = 2, cat.fontface=4, 
                          category.names=c("VIP3D.2", "VIP2D.comp12"), main="NPLSDANegLipLists",
                          cex.main=2
)

# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
plot.new()
grid.draw(venn.plot)

#pdf("VenndiagramsselectedNegLip.pdf")
#grid.draw(venn.plot)
#dev.off()
####################################################################################################################
# Test the Union and the intersection between the two metabolites sets

theunionNegLip<- union(selectedNegLipbydifVIPs[[1]],selectedNegLipbydifVIPs[[2]])
length(theunionNegLip) #160
theintersectionNegLip<- intersect(selectedNegLipbydifVIPs[[1]],selectedNegLipbydifVIPs[[2]])
length(theintersectionNegLip) # 91

####
FullarrayUNIONNegLipSelVars<-fullarrayNegLip136[,is.element(colnames(fullarrayNegLip136),
                                                            theunionNegLip),]
dim(FullarrayUNIONNegLipSelVars)

### Now NPLSDA and Graph
NPLSDAUnionNegLip<-NPLSDAmod(XN=FullarrayUNIONNegLipSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf ("UNIONNegLip.pdf")
ploteoNPLSDAUnionNegLip<- plotNPLSDAmod (X=NPLSDAUnionNegLip, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                         cutoff = 20, factors=2, penalty=2) 
#dev.off()

################################################################################################################
FullarrayINTERSECTNegLipSelVars<-fullarrayNegLip136[,is.element(colnames(fullarrayNegLip136),
                                                                theintersectionNegLip),]
dim(FullarrayINTERSECTNegLipSelVars)

### Now NPLSDA and Graph
NPLSDAIntersectionNegLip<-NPLSDAmod(XN=FullarrayINTERSECTNegLipSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf ("IntersectionNegLip.pdf")
ploteoNPLSDAIntersectionNegLip<- plotNPLSDAmod (X=NPLSDAIntersectionNegLip, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                                cutoff = 20, factors=2, penalty=2) 
#dev.off()
#summary(selectedNegLipbydifVIPs)
#save(selectedNegLipbydifVIPs, file="selectedNegLipbydifVIPs.RData")
################################################################################################################

#############################################################################################################
###################################    Variable Selection PosLip     ########################################
#############################################################################################################
# Variable selection by VIP3Dmodel2
summary(NPLSDA136PosLip)
NPLSDA136PosLip$VIP3Dmodel2

### Desde aca
vipsoutcomemet<-data.frame(NPLSDA136PosLip$VIP3Dmodel2)

thrs<-99
vipsoutcomemetsubset95<-NULL
temp<-NULL
vipsoutcomemetsubset<-list()
for (i in 1:length(vipsoutcomemet)){
  a<-quantile(vipsoutcomemet[,i], probs = c(50, 95,99, thrs)/100)
  print(i)
  print(a)
  temp<-vipsoutcomemet[vipsoutcomemet[,i]>=a[4],]
  namesubset<-paste('t',i,sep='')
  vips<-list(temp)
  vipsoutcomemetsubset[[namesubset]]<-vips
}
vipsoutcomemetsubset
rm(vaya)
#rownames(vipsoutcomemetsubset$t1[[1]])
vaya<-c (as.vector(rownames(vipsoutcomemetsubset$t1[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t2[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t3[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t4[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t5[[1]]))) 
duplicated(vaya)
PosLipselVars<-vaya[!duplicated(vaya)]
length(PosLipselVars)
######################################################################################################
# Retain just these variables in GCTOF 
dim(fullarrayPosLip136)
PosLipselVars
#selectedPosLipbydifVIPs<-list(VIP3Dmodel2.99p.28v=PosLipselVars)

#summary(selectedGCTOFbydifVIPs)
#selectedgenesbydifVIPs<-c(selectedgenesbydifVIPs,tulo="pelotas")

FullarrayGenesVIPSelVars<-fullarrayPosLip136[,is.element(colnames(fullarrayPosLip136),
                                                         PosLipselVars),]
dim(FullarrayGenesVIPSelVars)

#FullarrayGenesVIPSelVarsVIP3D.2<-FullarrayGenesVIPSelVars
#save(FullarrayGenesVIPSelVarsVIP3D.2,file="FullarrayGenesVIPSelVarsVIP3D.2.RData")

### Now NPLSDA and Graph
NPLSDAarraygenesVIPselected<-NPLSDAmod(XN=FullarrayGenesVIPSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

#NPLSDAarraygenesVIPselected$NPLSDAvariatesperMode3
#NPLSDAarraygenesVIPselected$NPLSDAvariates$X
#dim(NPLSDAarraygenesVIPselected$NPLSDAvariatesperMode3)
# Plotting
#pdf("VIP3Dmodel2.99p.1013g.pdf")
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=2, penalty=2) 
#dev.off()

ploteoNPLSDAVIPselectedgenes

######################################################################################################################
# Variable selection by VIP3Dmodel1

summary(NPLSDA136PosLip)
NPLSDA136PosLip$VIP3Dmodel1
dim(NPLSDA136PosLip$VIP3Dmodel1)



### Desde aca
vipsoutcomemet<-data.frame(NPLSDA136PosLip$VIP3Dmodel1)

thrs<-95
vipsoutcomemetsubset95<-NULL
temp<-NULL
vipsoutcomemetsubset<-list()
for (i in 1:2){
  a<-quantile(vipsoutcomemet[,i], probs = c(50, 95,99, thrs)/100)
  print(i)
  print(a)
  temp<-vipsoutcomemet[vipsoutcomemet[,i]>=a[4],]
  namesubset<-paste('t',i,sep='')
  vips<-list(temp)
  vipsoutcomemetsubset[[namesubset]]<-vips
}
vipsoutcomemetsubset
rm(vaya)
#rownames(vipsoutcomemetsubset$t1[[1]])
vaya<-c (as.vector(rownames(vipsoutcomemetsubset$t1[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t2[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t3[[1]]))) 
duplicated(vaya)
PosLipselVars<-vaya[!duplicated(vaya)]
length(PosLipselVars)  # 2050
######################################################################################################
# Retain just these variables in Gene Expression 2D
dim(fullarrayPosLip136)
PosLipselVars

FullarrayGenesVIPSelVars<-fullarrayPosLip136[,is.element(colnames(fullarrayPosLip136),
                                                         PosLipselVars),]
dim(FullarrayGenesVIPSelVars)
#FullarrayGenesVIPSelVarsVIP1.comp1<-FullarrayGenesVIPSelVars
#FullarrayGenesVIPSelVarsVIP1.comps12<-FullarrayGenesVIPSelVars
#save(FullarrayGenesVIPSelVarsVIP1.comp1,file="FullarrayGenesVIPSelVarsVIP1.comp1.RData")
#save(FullarrayGenesVIPSelVarsVIP1.comps12,file="FullarrayGenesVIPSelVarsVIP1.comps12.RData")


### Now NPLSDA and Graph
NPLSDAarraygenesVIPselected<-NPLSDAmod(XN=FullarrayGenesVIPSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf("VIP3Dmodel1.99p.1comp.213g.pdf")
#pdf("VIP3Dmodel1.99p.12comp.336g.pdf")
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=2, penalty=2) 
#dev.off()


#summary(selectedPosLipbydifVIPs)
#selectedPosLipbydifVIPs<-append(selectedPosLipbydifVIPs,list(VIP3Dmodel1.95p.comp12.46v=PosLipselVars))

#save(selectedPosLipbydifVIPs, file="selectedPosLipbydifVIPs.RData")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/selectedPosLipbydifVIPs.RData")
######################################################################################################
# Gene Expression USING VIP2D

summary(NPLSDA136PosLip)
NPLSDA136PosLip$VIP2D


### Desde aca
vipsoutcome2D<-data.frame(NPLSDA136PosLip$VIP2D)

colp1<-paste(rownames(NPLSDA136PosLip$VIP3Dmodel1),"-12", sep="!_____")
colp2<-paste(rownames(NPLSDA136PosLip$VIP3Dmodel1),"-9", sep="!_____")
colp3<-paste(rownames(NPLSDA136PosLip$VIP3Dmodel1),"-6", sep="!_____")
colp4<-paste(rownames(NPLSDA136PosLip$VIP3Dmodel1),"-3", sep="!_____")
colp5<-paste(rownames(NPLSDA136PosLip$VIP3Dmodel1),"0", sep="!_____")


rownames(vipsoutcome2D)<-c(colp1,colp2,colp3,colp4,colp5)

#apply(vipsoutcomemet, 2, function(x) is.numeric(x))
vipsoutcomemet<-vipsoutcome2D
thrs<-95
vipsoutcomemetsubset95<-NULL
temp<-NULL
vipsoutcomemetsubset<-list()
for (i in 1:2){
  a<-quantile(vipsoutcomemet[,i], probs = c(50, 95,99, thrs)/100)
  print(i)
  print(a)
  temp<-vipsoutcomemet[vipsoutcomemet[,i]>=a[4],]
  namesubset<-paste('t',i,sep='')
  vips<-list(temp)
  vipsoutcomemetsubset[[namesubset]]<-vips
}
vipsoutcomemetsubset
rm(vaya)
#rownames(vipsoutcomemetsubset$t1[[1]])
vaya<-c (as.vector(rownames(vipsoutcomemetsubset$t1[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t2[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t3[[1]]))) 
duplicated(vaya)
PosLipselVars<-vaya[!duplicated(vaya)]
length(PosLipselVars)
string<-strsplit (PosLipselVars,"!_____")
string<-listTOdata.frame(string)
genelist<- string[,1]
duplicated (genelist)
genes<-genelist[!duplicated(genelist)]
length(genes)

######################################################################################################
# Retain just these variables in Gene Expression 2D
dim(fullarrayPosLip136)
genes
#as.vector(genes)
FullarrayGenesVIPSelVars<-fullarrayPosLip136[,is.element(colnames(fullarrayPosLip136),
                                                         genes),]
dim(FullarrayGenesVIPSelVars)
#FullarrayGenesVIPSelVarsVIP2Dcomp2<-FullarrayGenesVIPSelVars
#save(FullarrayGenesVIPSelVarsVIP2Dcomp2,file="FullarrayGenesVIPSelVarsVIP2Dcomp2.RData")

### Now NPLSDA and Graph
NPLSDAarraygenesVIPselected<-NPLSDAmod(XN=FullarrayGenesVIPSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf ("VIP2D.99p.comp2.1022g.pdf")
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=2, penalty=2) 
#dev.off()

#summary(selectedNegLipbydifVIPs)
#selectedNegLipbydifVIPs<-append(selectedNegLipbydifVIPs,list(VIP2D.95p.comp12.153v= as.vector(genes)))
#save(selectedNegLipbydifVIPs,file="selectedNegLipbydifVIPs.RData")
###########################################################################################################
# We selected for PosLip: VIP3D/2 99% and VIP3D/1 95%

# Venn Diagram of the four genelists #

summary(selectedPosLipbydifVIPs)

# If this is a bit confusing you can also write a function and then use it in lapply() 
removeEMPTYstrings <- function(x) {
  
  newVectorWOstrings <- x[x != ""]
  return(newVectorWOstrings)
  
}
geneLS2 <- lapply(selectedPosLipbydifVIPs, removeEMPTYstrings)

summary(geneLS2)
lapply(geneLS2, tail) # Both methods return the same results

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:


VENN.LIST <- geneLS2
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkmagenta", "darkblue"),#,"darkseagreen4","goldenrod2"), 
                          alpha=c(0.5,0.5), cex = 2, cat.fontface=4, 
                          category.names=c("VIP3D.2", "VIP2D.comp12"), main="NPLSDAPosLipLists",
                          cex.main=2
)

# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
grid.draw(venn.plot)

#pdf("VenndiagramsselectedPosLip.pdf")
#grid.draw(venn.plot)
#dev.off()
####################################################################################################################
# Test the Union and the intersection between the two metabolites sets

theunionPosLip<- union(selectedPosLipbydifVIPs[[1]],selectedPosLipbydifVIPs[[2]])
length(theunionPosLip) #63
theintersectionPosLip<- intersect(selectedPosLipbydifVIPs[[1]],selectedPosLipbydifVIPs[[2]])
length(theintersectionPosLip) # 11

####
FullarrayUNIONPosLipSelVars<-fullarrayPosLip136[,is.element(colnames(fullarrayPosLip136),
                                                            theunionPosLip),]
dim(FullarrayUNIONPosLipSelVars)

### Now NPLSDA and Graph
NPLSDAUnionPosLip<-NPLSDAmod(XN=FullarrayUNIONPosLipSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf ("UNIONPosLip.pdf")
ploteoNPLSDAUnionPosLip<- plotNPLSDAmod (X=NPLSDAUnionPosLip, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                         cutoff = 20, factors=2, penalty=2) 
#dev.off()

################################################################################################################
FullarrayINTERSECTPosLipSelVars<-fullarrayPosLip136[,is.element(colnames(fullarrayPosLip136),
                                                                theintersectionPosLip),]
dim(FullarrayINTERSECTPosLipSelVars)

### Now NPLSDA and Graph
NPLSDAIntersectionPosLip<-NPLSDAmod(XN=FullarrayINTERSECTPosLipSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf ("IntersectionPosLip.pdf")
ploteoNPLSDAIntersectionPosLip<- plotNPLSDAmod (X=NPLSDAIntersectionPosLip, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                                cutoff = 20, factors=2, penalty=2) 
#dev.off()





#####################################################
######     Merge of all Metabolomics Data     #######
#####################################################
#Merge winners

###########################################################
###########################################################
#################     Checkpoint 5     ####################
###########################################################
###########################################################
# From this point on you can continue just by importing FullarrayUNIONGCTOFSelVars, 
# FullarrayINTERSECTNegLipSelVars and FullarrayUNIONPosLipSelVars
###########################################################
dim(FullarrayUNIONGCTOFSelVars)          # 91
#save(FullarrayUNIONGCTOFSelVars, file= "FullarrayUNIONGCTOFSelVars.RData")
dim(FullarrayINTERSECTNegLipSelVars)     # 91
sort(table(colnames(FullarrayINTERSECTNegLipSelVars)))
#save(FullarrayINTERSECTNegLipSelVars, file= "FullarrayINTERSECTNegLipSelVars.RData")
dim(FullarrayUNIONPosLipSelVars)         # 63
sort(table(colnames(FullarrayUNIONPosLipSelVars)))
#save(FullarrayUNIONPosLipSelVars, file="FullarrayUNIONPosLipSelVars.RData")

###########################################################
# There are repeated names in PosLip and NegLip, so we are adding a little pieace of name
FullarrayINTERSECTNegLipSelVars
colp1<-paste(colnames(FullarrayINTERSECTNegLipSelVars),"NegLip", sep="!_")
colnames(FullarrayINTERSECTNegLipSelVars)<-colp1

FullarrayUNIONPosLipSelVars
colp2<-paste(colnames(FullarrayUNIONPosLipSelVars),"PosLip", sep="!_")
colnames(FullarrayUNIONPosLipSelVars)<-colp2


posneglip<-abind(FullarrayUNIONPosLipSelVars,FullarrayINTERSECTNegLipSelVars,along=2);dim(FullarrayUNIONPosLipSelVars);dim(FullarrayINTERSECTNegLipSelVars);dim(posneglip)
Allmetabolomics136<-abind(posneglip,FullarrayUNIONGCTOFSelVars,along=2);dim(posneglip);dim(FullarrayUNIONGCTOFSelVars);dim(Allmetabolomics136)
#
#save(Allmetabolomics136,file="Allmetabolomics136.RData")
DefinitiveMetabolomicsListFeb2<-as.data.frame(colnames(Allmetabolomics136))
#write.table(DefinitiveMetabolomicsListFeb2, "/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/DefinitiveMetabolomicsListFeb2.txt", sep="\t")
#DefinitiveMetabolomicsListFeb2 <- read.csv("~/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/DefinitiveMetabolomicsListFeb2.txt", sep="")

###########################################################
# NPLSDA and Plotting
dim(outcomedummyarray136)
dim(Allmetabolomics136)
#colnames(Allmetabolomics136)
#sort(table(colnames(Allmetabolomics136)))

NPLSDAallmetabolomicsSelVars<-NPLSDAmod(XN=Allmetabolomics136, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 


#save(NPLSDAallmetabolomicsSelVars,file="NPLSDAallmetabolomicsSelVars.RData")
# Plotting
#pdf ("NPLSDAallmetabolomicsSelVars.pdf")
ploteo<- plotNPLSDAmod (X=NPLSDAallmetabolomicsSelVars, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                        cutoff = 20, factors=2, penalty=2)
#dev.off()

###########################################################
# Repeat the same procedure for processed data