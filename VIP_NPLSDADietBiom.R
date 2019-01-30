###########################################################
##############     VIP_NPLSDADietBiom.R     ###############
###########################################################
# Author: Leandro Balzano-Nogueira
# Genetics Institute, University of Florida (Gainesville)

# This script is to create the Dietary Biomarkers  dataset to calculate the NPLSDA and VIP selection

###########################################################
homedir<- "/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/DietaryBiomarkers/" # Home directory where all your results are going to be contained
setwd(homedir)
getwd()
###########################################################
# Functions:
"/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/ScriptsForNComm/Tools/" # Location of TEDDYtools
source ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/ScriptsForNComm/Tools/TEDDYtools2.R") # These are the functions to reformat the data
source ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/ScriptsForNComm/Tools/NPLSDAfunctionsApr11.R") # These are the functions created to perform the NPLSDA

###########################################################
# Data:
# Dietary Biomarkers
DBRaw<-read.csv ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/DietaryBiomarkers/DBRaw.csv",header = FALSE)
DBRaw[1:10,1:10]

DBprocessed<-read.csv ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/DietaryBiomarkers/DBprocessed.csv",header = FALSE)
DBprocessed[1:10,1:10]

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
# Gene Expression:
DietBiomRaw<-t(DBRaw)
colnames(DietBiomRaw)<- DietBiomRaw[1,]; DietBiomRaw<-DietBiomRaw[-1,]

DietBiomRaw <- transform(DietBiomRaw, Individual.Id =as.numeric(as.character(Individual.Id)),
                           Age.in.Months = as.numeric(as.character(Age.in.Months)), 
                           Time.to.IA = as.numeric(as.character(Time.to.IA)))

DietBiomRaw[1:10,1:10]

# Subsetting the data to the 136 pairs with values reported in 3 0ut 0f 5 time points
DietBiomRaw136<- DietBiomRaw[is.element(DietBiomRaw$Individual.Id, patients3tps[,1]),]
dim(DietBiomRaw136)
length(unique(DietBiomRaw136$Individual.Id))

###########################################################
# GE data in 3D structure
dim(DietBiomRaw136)   # 276 x 46
colnames(DietBiomRaw136[1:23])
DBX<-DietBiomRaw136[,c(1,4:46)]       # Just ID and the variables


DBX12<-DBX[DietBiomRaw136$Time.to.IA == "-12",]; dim(DBX12);dim(DBX)
DBX9<-DBX[DietBiomRaw136$Time.to.IA == "-9",];dim(DBX9)
DBX6<-DBX[DietBiomRaw136$Time.to.IA == "-6",];dim(DBX6)
DBX3<-DBX[DietBiomRaw136$Time.to.IA == "-3",];dim(DBX3)
DBX0<-DBX[DietBiomRaw136$Time.to.IA == "0",];dim(DBX0)


###########################################################
# Merging with all cases
patients3tps2<-patients3tps
colnames(patients3tps2)<- "Individual.Id"
DBX12total<-merge(patients3tps2,DBX12, by="Individual.Id", all.x = T);dim(DBX12total);dim(patients3tps2);dim(DBX12)
DBX12total[1:5,1:5]
rownames(DBX12total)<- DBX12total[,1]
DBX12total<- as.matrix(DBX12total[,c(-1)])
dim(DBX12total)

DBX9total<-merge(patients3tps2,DBX9, by="Individual.Id", all.x = T);dim(DBX9total);dim(patients3tps2);dim(DBX9)
DBX9total[1:5,1:5]
rownames(DBX9total)<- DBX9total[,1]
DBX9total<- as.matrix(DBX9total[,c(-1)])
dim(DBX9total)

DBX6total<-merge(patients3tps2,DBX6, by="Individual.Id", all.x = T);dim(DBX6total);dim(patients3tps2);dim(DBX6total)
DBX6total[1:5,1:5]
rownames(DBX6total)<- DBX6total[,1]
DBX6total<- as.matrix(DBX6total[,c(-1)])
dim(DBX6total)

DBX3total<-merge(patients3tps2,DBX3, by="Individual.Id", all.x = T);dim(DBX3total);dim(patients3tps2);dim(DBX3)
DBX3total[1:5,1:5]
rownames(DBX3total)<- DBX3total[,1]
DBX3total<- as.matrix(DBX3total[,c(-1)])
dim(DBX3total)

DBX0total<-merge(patients3tps2,DBX0, by="Individual.Id", all.x = T);dim(DBX0total);dim(patients3tps2);dim(DBX0)
DBX0total[1:5,1:5]
rownames(DBX0total)<- DBX0total[,1]
DBX0total<- as.matrix(DBX0total[,c(-1)])
dim(DBX0total)
# Dimensions are 136*43*5

###########################################################
dim(DBX0total)
arrayDBXMarch <- array(data = NA, dim = c(136,43,5),dimnames = list(NULL, NULL, c("-12","-9","-6", "-3", "0")))

arrayDBXMarch
arrayDBXMarch[,,1] <- DBX12total
arrayDBXMarch[,,2] <- DBX9total
arrayDBXMarch[,,3] <- DBX6total
arrayDBXMarch[,,4] <- DBX3total
arrayDBXMarch[,,5] <- DBX0total

rownames(arrayDBXMarch)<-rownames(DBX12total)
colnames(arrayDBXMarch)<-colnames (DBX12total)
arrayDBXMarch[1:50,1:5,1]
arrayDBXMarch136<-arrayDBXMarch
dim(arrayDBXMarch136)   # 136 * 43 * 5

# This is how the data array for DB expression was created 

###########################################################
# Imputation (For convenience it must be done in an HPC)
# 1) Determining the best fitted model

modelDBXnuevo<-bestfittedmodel (X=arrayDBXMarch136,centering=0) # 0= No centering; 1= centering by Individuals; 2= centering by Variables;3= centering by Time
# The best model was 4,4,3 

# 2) Imputing the best fitted model data
FullarrayGEMARCH136<-Imputemethod(X=arrayDBXMarch136,fac=c(4, 4, 3), conver = 1e-07, max.iter = 1000)

summary(FullarrayGEMARCH136)
dim(FullarrayGEMARCH136)


# 3) NPLSDA (For convenience it must be done in an HPC)
NPLSDAFullarrayGEMARCH136<-NPLSDAmod(XN=FullarrayGEMARCH136, YN=outcomedummyarray136, outcome.Y=NULL, factors=3, centering=0) 

# 4) Plotting

ploteoNPLSDAFullarrayGEMARCH136<- plotNPLSDAmod (X=NPLSDAFullarrayGEMARCH136, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                                 cutoff = 20, factors=2, penalty=1) 

###########################################################
# Repeat the same procedure for processed data
# The process allow us to select Vitamin C, Vitamin D and Alpha-tocopherol as variables for the final model
