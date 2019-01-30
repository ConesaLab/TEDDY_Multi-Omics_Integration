###########################################################
##########     PartialCorrelationNetworks.R     ###########
###########################################################
# Author: Leandro Balzano-Nogueira
# Genetics Institute, University of Florida (Gainesville)

# This script is to perform the Partial Correlation networks reported
# We Perform ParCor to the cases only, then subset the correlation matrix and plot the network
# Also, this analysis, we are selecting the central core of correlation by a Shrinkage procedure...
# ... based on the lowest root mean square error of prediction (RMSEP) value from...
# ... Ridge, Lasso or ElasticNet methods to guarantee the best fitted model.


###########################################################
homedir<- "/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/" # Home directory where all your results are going to be contained
setwd(homedir)
getwd()
###########################################################
# Data: 
wholedescriptivevars<-read.csv ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/LM_globalTargets.csv",header = TRUE, row.names = 1)
wholedescriptivevars

load("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/DeltaMatrices_FullArray.Rdata")


load("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/FullarrayDB.RData")
AllmetabolitesConverted <- read.csv("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/AllmetabolitesConverted.txt", sep= "\t")
definitivemetabolitesSelectedConvertedmarch9 <- read.csv("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/definitivemetabolitesSelectedConvertedmarch9.txt", sep="\t")
load("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/FullarrayGEMETABDBBestModel.RData")

# List of Cases with at least 3 out of 5 time points with data
patients3tps<-data.frame(V1=CohortData$Individual.Id[CohortData$Model.or.Validation=="Model"])
patients3tps

CohortData<-read.csv ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/CohortData.csv",header = TRUE)
CohortData[1:10,]

###########################################################
# Libraries:
require("ropls")
library(reshape2)
library("glmnet", lib.loc="/usr/local/lib/R/site-library")
library("qgraph")
library("parcor")
require(reshape)

###########################################################
CohortData2<-CohortData[,c(2,8)]
rownames(CohortData2)<-CohortData2[,1]
namecolumn<-"Outcome"
CohortData2[,namecolumn]<-NA

for (i in 1:nrow(CohortData2)){
  if(CohortData2$Case.or.Control[i] == "Case" ) {
    print("Es Case")
    CohortData2$Outcome[i] <- 1
  }
  if(CohortData2$Case.or.Control[i] == "Control" ) {
    print("Es Control")
    CohortData2$Outcome[i] <- 0
  }
}
CohortData2
dim(CohortData2)
CohortData2[1:10,]
outcomedummyarray306<-as.data.frame(CohortData2[,-c(1,2)])
rownames(outcomedummyarray306)<-rownames(CohortData2)
colnames(outcomedummyarray306)<-"Outcome"
outcomedummyarray306

outcomedummyarray136<-outcomedummyarray306[rownames(outcomedummyarray306) %in% patients3tps$V1,,drop=FALSE]

###########################################################
class(DeltaMatrices_FullArray)
summary(DeltaMatrices_FullArray)
dim(DeltaMatrices_FullArray$`Period 12-9`)   # 136 x 1107
DeltaMatrices_FullArray$`Period 12-9`[1:10,1:10]
###########################################################
# First we have to add to each table the deltas from the three DB variables
dim(FullarrayDB)

FullarrayDB136<- FullarrayDB[rownames(FullarrayDB) %in% rownames(outcomedummyarray136),,]
dim(FullarrayDB136)
colnames(FullarrayDB136) # Retain just vitD (15),VitC (43) and alphaTocopherol(16)
PieceFullarrayDB136<-FullarrayDB136[,c(15,16,43),]
dim(PieceFullarrayDB136)
PieceFullarrayDB136[1:10,,1]
###########################################################
# Changing the colnames of the metabolomics variables
colnames(DeltaMatrices_FullArray$`Period 12-9`)
###  Time period -12 -9
tulo<- DeltaMatrices_FullArray$`Period 12-9`[,863:1107]
dim(tulo)
colnames(tulo)
colnames(tulo)<- definitivemetabolitesSelectedConvertedmarch9[ match( colnames( tulo ) ,definitivemetabolitesSelectedConvertedmarch9[ , "ID_Var" ]),
                                                               "FunCat" ]

DeltaMatrices_FullArraym12m9<-merge(DeltaMatrices_FullArray$`Period 12-9`[,1:862],tulo, by=0)
dim(DeltaMatrices_FullArray$`Period 12-9`); dim (tulo);dim(DeltaMatrices_FullArraym12m9)
rownames(DeltaMatrices_FullArraym12m9)<-DeltaMatrices_FullArraym12m9[,1]
DeltaMatrices_FullArraym12m9 <-DeltaMatrices_FullArraym12m9[,-1]
dim(DeltaMatrices_FullArraym12m9)
DeltaMatrices_FullArraym12m9[1:10,900:910]

###  Time period -9 -6
tulo<- DeltaMatrices_FullArray$`Period 9-6`[,863:1107]
dim(tulo)
colnames(tulo)
colnames(tulo)<- definitivemetabolitesSelectedConvertedmarch9[ match( colnames( tulo ) ,definitivemetabolitesSelectedConvertedmarch9[ , "ID_Var" ]),
                                                               "FunCat" ]

DeltaMatrices_FullArraym9m6<-merge(DeltaMatrices_FullArray$`Period 9-6`[,1:862],tulo, by=0)
dim(DeltaMatrices_FullArray$`Period 9-6`); dim (tulo);dim(DeltaMatrices_FullArraym9m6)
rownames(DeltaMatrices_FullArraym9m6)<-DeltaMatrices_FullArraym9m6[,1]
DeltaMatrices_FullArraym9m6 <-DeltaMatrices_FullArraym9m6[,-1]
dim(DeltaMatrices_FullArraym9m6)
DeltaMatrices_FullArraym9m6[1:10,1100:1107]

###  Time period -6 -3
tulo<- DeltaMatrices_FullArray$`Period 6-3`[,863:1107]
dim(tulo)
colnames(tulo)
colnames(tulo)<- definitivemetabolitesSelectedConvertedmarch9[ match( colnames( tulo ) ,definitivemetabolitesSelectedConvertedmarch9[ , "ID_Var" ]),
                                                               "FunCat" ]

DeltaMatrices_FullArraym6m3<-merge(DeltaMatrices_FullArray$`Period 6-3`[,1:862],tulo, by=0)
dim(DeltaMatrices_FullArray$`Period 6-3`); dim (tulo);dim(DeltaMatrices_FullArraym6m3)
rownames(DeltaMatrices_FullArraym6m3)<-DeltaMatrices_FullArraym6m3[,1]
DeltaMatrices_FullArraym6m3 <-DeltaMatrices_FullArraym6m3[,-1]
dim(DeltaMatrices_FullArraym6m3)
DeltaMatrices_FullArraym6m3[1:10,1100:1107]

###  Time period -3 to 0
tulo<- DeltaMatrices_FullArray$`Period 3-0`[,863:1107]
dim(tulo)
colnames(tulo)
colnames(tulo)<- definitivemetabolitesSelectedConvertedmarch9[ match( colnames( tulo ) ,definitivemetabolitesSelectedConvertedmarch9[ , "ID_Var" ]),
                                                               "FunCat" ]

DeltaMatrices_FullArraym3m0<-merge(DeltaMatrices_FullArray$`Period 3-0`[,1:862],tulo, by=0)
dim(DeltaMatrices_FullArray$`Period 3-0`); dim (tulo);dim(DeltaMatrices_FullArraym3m0)
rownames(DeltaMatrices_FullArraym3m0)<-DeltaMatrices_FullArraym3m0[,1]
DeltaMatrices_FullArraym3m0 <-DeltaMatrices_FullArraym3m0[,-1]
dim(DeltaMatrices_FullArraym3m0)
DeltaMatrices_FullArraym3m0[1:10,1100:1107]

###########################################################
# Substraction of tf-to per Dietary biomarker variable
tablem12m9<-PieceFullarrayDB136[,,2] - PieceFullarrayDB136[,,1]
PieceFullarrayDB136[1:10,,1]
PieceFullarrayDB136[1:10,,2]
tablem12m9 [1:10,] # perfect
dim(tablem12m9)
anyNA(tablem12m9)
tablem9m6<-PieceFullarrayDB136[,,3] - PieceFullarrayDB136[,,2]
tablem6m3<-PieceFullarrayDB136[,,4] - PieceFullarrayDB136[,,3]
tablem3m0<-PieceFullarrayDB136[,,5] - PieceFullarrayDB136[,,4]
dim(tablem12m9);dim(tablem9m6);dim(tablem6m3);dim(tablem3m0)

###########################################################
# Merging this tables to the array
### Time period m9 -m12
Deltaallvarsm12m9<- merge (DeltaMatrices_FullArraym12m9,tablem12m9, by= 0)
dim(DeltaMatrices_FullArraym12m9);dim(tablem12m9); dim(Deltaallvarsm12m9)
Deltaallvarsm12m9[1:10,1:10]
rownames(Deltaallvarsm12m9)<-Deltaallvarsm12m9[,1]
Deltaallvarsm12m9 <-Deltaallvarsm12m9[,-1]
Deltaallvarsm12m9[1:10,1:10]
dim(Deltaallvarsm12m9)
colnames(Deltaallvarsm12m9)
Deltaallvarsm12m9[1:10,1100:1110]

### Time period m9 -m6
Deltaallvarsm9m6<- merge (DeltaMatrices_FullArraym9m6,tablem9m6, by= 0)
dim(DeltaMatrices_FullArraym9m6);dim(tablem9m6); dim(Deltaallvarsm9m6)
Deltaallvarsm9m6[1:10,1:10]
rownames(Deltaallvarsm9m6)<-Deltaallvarsm9m6[,1]
Deltaallvarsm9m6 <-Deltaallvarsm9m6[,-1]
Deltaallvarsm9m6[1:10,1:10]
dim(Deltaallvarsm9m6)
colnames(Deltaallvarsm9m6)
### Time period m6 -m3
Deltaallvarsm6m3<- merge (DeltaMatrices_FullArraym6m3,tablem6m3, by= 0)
dim(DeltaMatrices_FullArraym6m3);dim(tablem6m3); dim(Deltaallvarsm6m3)
Deltaallvarsm6m3[1:10,1:10]
rownames(Deltaallvarsm6m3)<-Deltaallvarsm6m3[,1]
Deltaallvarsm6m3 <-Deltaallvarsm6m3[,-1]
Deltaallvarsm6m3[1:10,1:10]
dim(Deltaallvarsm6m3)
colnames(Deltaallvarsm6m3)
### Time period m3 - m0
Deltaallvarsm3m0<- merge (DeltaMatrices_FullArraym3m0,tablem3m0, by= 0)
dim(DeltaMatrices_FullArraym3m0);dim(tablem3m0); dim(Deltaallvarsm3m0)
Deltaallvarsm3m0[1:10,1:10]
rownames(Deltaallvarsm3m0)<-Deltaallvarsm3m0[,1]
Deltaallvarsm3m0 <-Deltaallvarsm3m0[,-1]
Deltaallvarsm3m0[1:10,1:10]
dim(Deltaallvarsm3m0)
colnames(Deltaallvarsm3m0)

###########################################################
###########################################################
#################     Shrinkage method     ################
###########################################################
###########################################################
# This have to be done with the CASES ONLY
dim(AgeGroupdummycasesfrom136)

Deltaallvarsm12m9cases<-Deltaallvarsm12m9 [is.element(rownames (Deltaallvarsm12m9),rownames(AgeGroupdummycasesfrom136)),];dim(Deltaallvarsm12m9cases)

Deltaallvarsm9m6cases<-Deltaallvarsm9m6 [is.element(rownames (Deltaallvarsm9m6),rownames(AgeGroupdummycasesfrom136)),];dim(Deltaallvarsm9m6cases)
Deltaallvarsm6m3cases<-Deltaallvarsm6m3 [is.element(rownames (Deltaallvarsm6m3),rownames(AgeGroupdummycasesfrom136)),];dim(Deltaallvarsm6m3cases)
Deltaallvarsm3m0cases<-Deltaallvarsm3m0 [is.element(rownames (Deltaallvarsm3m0),rownames(AgeGroupdummycasesfrom136)),];dim(Deltaallvarsm3m0cases)
###########################################################
# Creating the FAAb vector
colnames(wholedescriptivevars)
preFAAbcases<-wholedescriptivevars[wholedescriptivevars$Mask.Id %in% rownames(Deltaallvarsm12m9cases),]
dim(preFAAbcases)
FAAbcases<-unique(preFAAbcases[,c(1,4)])
dim(FAAbcases)
dim(Deltaallvarsm12m9cases)

FAAbcasesGADAIAA<-FAAbcases[FAAbcases$FirstAAb != 4 & FAAbcases$FirstAAb != 5 & FAAbcases$FirstAAb != 7,]
dim(FAAbcasesGADAIAA)
table(FAAbcasesGADAIAA$FirstAAb)

###########################################################

dim(FullarrayGEMETABDBBestModel)
FullarrayGEMETABDBBestModel [1:10,1:10,1]
FullarrayGEMETABDBBestModelCASESGADAIAA<-FullarrayGEMETABDBBestModel[rownames(FullarrayGEMETABDBBestModel) %in% FAAbcasesGADAIAA$Mask.Id,,]
dim(FullarrayGEMETABDBBestModelCASESGADAIAA)
FAAbcases$FirstAAb

###########################################################
# For Loop to detect the best alpha value per time point


Lambda.minm12<-Lambda.minm9<-Lambda.minm6<-Lambda.minm3<-Lambda.minm0<- NULL
alphaValuesm12<-alphaValuesm9<-alphaValuesm6<-alphaValuesm3<-alphaValuesm0<-NULL
alpha_lambdam12<-alpha_lambdam9<-alpha_lambdam6<-alpha_lambdam3<-alpha_lambdam0<-NULL
Lambda.minm12tita<-Lambda.minm9tita<-Lambda.minm6tita<-Lambda.minm3tita<-Lambda.minm0tita<-NULL
alphaWinnerm12tita<-alphaWinnerm9tita<-alphaWinnerm6tita<-alphaWinnerm3tita<-alphaWinnerm0tita<-NULL
alphaWinnerm12TOTA<-alphaWinnerm9TOTA<-alphaWinnerm6TOTA<-alphaWinnerm3TOTA<-alphaWinnerm0TOTA<-NULL
alphaWinnerm12<-alphaWinnerm9<-alphaWinnerm6<-alphaWinnerm3<-alphaWinnerm0<-NULL

#i=1
for (i in 1:11) {
  print (paste ("iteration to obtain best alphas for Feature selection by Elastic Net: ",i, sep=" "))
  for(alpha_val in seq(0,1, by=0.2)){
    ############
    # TIME -12 #
    ############
    #for(alpha_val in 0.2){
    #alpha_val=0.4
    #alpha_val=seq(0,1, by=0.2)
    cfit<-cv.glmnet(as.matrix(FullarrayGEMETABDBBestModelCASESGADAIAA [,,1]),FAAbcasesGADAIAA$FirstAAb, 
                    standardize=TRUE, family="binomial", 
                    alpha = alpha_val,
                    nfolds = 5,
                    type.measure = "auc") # cross validation with  
    #plot(cfit)
    #print(paste("alpham12 =",alpha_val, ", lambda.minm12 =", cfit$lambda.min))
    Lambda.minm12tita <-c(Lambda.minm12tita,cfit$lambda.min)
    alphaValuesm12<- c(alphaValuesm12,alpha_val)
    alpha_lambdam12ROW<-cbind(alphaValuesm12,Lambda.minm12tita)
    alpha_lambdam12<-rbind(alpha_lambdam12,alpha_lambdam12ROW)
    alpha_lambdam12<-data.frame(alpha_lambdam12)
    #class(alpha_lambdam12)
    alpha_lambdam12<-alpha_lambdam12[order (alpha_lambdam12$Lambda.minm12tita),]
    #dim(alpha_lambdam12)
    Lambda.minm12 <- alpha_lambdam12$alphaValuesm12 [which(alpha_lambdam12$Lambda.minm12tita==min(alpha_lambdam12$Lambda.minm12tita))]
    #print(paste("Lambda.minm12 =",Lambda.minm12))
    alpha_lambdam12ROW<-alphaValuesm12<-Lambda.minm12tita<-NULL
    ###########
    # TIME -9 #
    ###########
    cfit<-cv.glmnet(as.matrix(FullarrayGEMETABDBBestModelCASESGADAIAA [,,2]),FAAbcasesGADAIAA$FirstAAb, 
                    standardize=TRUE, family="binomial", 
                    alpha = alpha_val,
                    nfolds = 5,
                    type.measure = "auc") # cross validation with  
    #plot(cfit)
    #print(paste("alpham12 =",alpha_val, ", lambda.minm12 =", cfit$lambda.min))
    Lambda.minm9tita <-c(Lambda.minm9tita,cfit$lambda.min)
    alphaValuesm9<- c(alphaValuesm9,alpha_val)
    alpha_lambdam9ROW<-cbind(alphaValuesm9,Lambda.minm9tita)
    alpha_lambdam9<-rbind(alpha_lambdam9,alpha_lambdam9ROW)
    alpha_lambdam9<-data.frame(alpha_lambdam9)
    #class(alpha_lambdam9)
    alpha_lambdam9<-alpha_lambdam9[order (alpha_lambdam9$Lambda.minm9tita),]
    #dim(alpha_lambdam9)
    Lambda.minm9 <- alpha_lambdam9$alphaValuesm9 [which(alpha_lambdam9$Lambda.minm9tita==min(alpha_lambdam9$Lambda.minm9tita))]
    #print(paste("Lambda.minm9 =",Lambda.minm9))
    alpha_lambdam9ROW<-alphaValuesm9<-Lambda.minm9tita<-NULL
    ###########
    # TIME -6 #
    ###########
    cfit<-cv.glmnet(as.matrix(FullarrayGEMETABDBBestModelCASESGADAIAA [,,3]),FAAbcasesGADAIAA$FirstAAb, 
                    standardize=TRUE, family="binomial", 
                    alpha = alpha_val,
                    nfolds = 5,
                    type.measure = "auc") # cross validation with  
    #plot(cfit)
    #print(paste("alpham12 =",alpha_val, ", lambda.minm12 =", cfit$lambda.min))
    Lambda.minm6tita <-c(Lambda.minm6tita,cfit$lambda.min)
    alphaValuesm6<- c(alphaValuesm6,alpha_val)
    alpha_lambdam6ROW<-cbind(alphaValuesm6,Lambda.minm6tita)
    alpha_lambdam6<-rbind(alpha_lambdam6,alpha_lambdam6ROW)
    alpha_lambdam6<-data.frame(alpha_lambdam6)
    #class(alpha_lambdam9)
    alpha_lambdam6<-alpha_lambdam6[order (alpha_lambdam6$Lambda.minm6tita),]
    #dim(alpha_lambdam6)
    Lambda.minm6 <- alpha_lambdam6$alphaValuesm6 [which(alpha_lambdam6$Lambda.minm6tita==min(alpha_lambdam6$Lambda.minm6tita))]
    #print(paste("Lambda.minm6 =",Lambda.minm6))
    alpha_lambdam6ROW<-alphaValuesm6<-Lambda.minm6tita<-NULL
    
    ###########
    # TIME -3 #
    ###########
    cfit<-cv.glmnet(as.matrix(FullarrayGEMETABDBBestModelCASESGADAIAA [,,4]),FAAbcasesGADAIAA$FirstAAb, 
                    standardize=TRUE, family="binomial", 
                    alpha = alpha_val,
                    nfolds = 5,
                    type.measure = "auc") # cross validation with  
    #plot(cfit)
    #print(paste("alpham12 =",alpha_val, ", lambda.minm12 =", cfit$lambda.min))
    Lambda.minm3tita <-c(Lambda.minm3tita,cfit$lambda.min)
    alphaValuesm3<- c(alphaValuesm3,alpha_val)
    alpha_lambdam3ROW<-cbind(alphaValuesm3,Lambda.minm3tita)
    alpha_lambdam3<-rbind(alpha_lambdam3,alpha_lambdam3ROW)
    alpha_lambdam3<-data.frame(alpha_lambdam3)
    #class(alpha_lambdam9)
    alpha_lambdam3<-alpha_lambdam3[order (alpha_lambdam3$Lambda.minm3tita),]
    #dim(alpha_lambdam3)
    Lambda.minm3 <- alpha_lambdam3$alphaValuesm3 [which(alpha_lambdam3$Lambda.minm3tita==min(alpha_lambdam3$Lambda.minm3tita))]
    #print(paste("Lambda.minm3 =",Lambda.minm3))
    alpha_lambdam3ROW<-alphaValuesm3<-Lambda.minm3tita<-NULL
    
    ##########
    # TIME 0 #
    ##########
    cfit<-cv.glmnet(as.matrix(FullarrayGEMETABDBBestModelCASESGADAIAA [,,5]),FAAbcasesGADAIAA$FirstAAb, 
                    standardize=TRUE, family="binomial", 
                    alpha = alpha_val,
                    type.measure = "auc") # cross validation with  
    #plot(cfit)
    #print(paste("alpham12 =",alpha_val, ", lambda.minm12 =", cfit$lambda.min))
    Lambda.minm0tita <-c(Lambda.minm0tita,cfit$lambda.min)
    alphaValuesm0<- c(alphaValuesm0,alpha_val)
    alpha_lambdam0ROW<-cbind(alphaValuesm0,Lambda.minm0tita)
    alpha_lambdam0<-rbind(alpha_lambdam0,alpha_lambdam0ROW)
    alpha_lambdam0<-data.frame(alpha_lambdam0)
    #class(alpha_lambdam9)
    alpha_lambdam0<-alpha_lambdam0[order (alpha_lambdam0$Lambda.minm0tita),]
    #dim(alpha_lambdam0)
    Lambda.minm0 <- alpha_lambdam0$alphaValuesm0 [which(alpha_lambdam0$Lambda.minm0tita==min(alpha_lambdam0$Lambda.minm0tita))]
    #print(paste("Lambda.minm6 =",Lambda.minm0))
    alpha_lambdam0ROW<-alphaValuesm0<-Lambda.minm0tita<-NULL
  }
  ############
  # TIME -12 #
  ############
  print(alpha_lambdam12)
  alphaWinnerm12tita<-alpha_lambdam12[1,1]
  print(paste("alphaWinnerm12tita=",alphaWinnerm12tita))
  alphaWinnerm12TOTA<-c(alphaWinnerm12TOTA,alphaWinnerm12tita)
  
  print(paste("alphaWinnerm12TOTA=",alphaWinnerm12TOTA))
  alphaWinnerm12tita<-NULL
  alpha_lambdam12<-NULL
  print(alpha_lambdam12)
  ###########
  # TIME -9 #
  ###########
  print(alpha_lambdam9)
  alphaWinnerm9tita<-alpha_lambdam9[1,1]
  print(paste("alphaWinnerm9tita=",alphaWinnerm9tita))
  alphaWinnerm9TOTA<-c(alphaWinnerm9TOTA,alphaWinnerm9tita)
  
  print(paste("alphaWinnerm9TOTA=",alphaWinnerm9TOTA))
  alphaWinnerm9tita<-NULL
  alpha_lambdam9<-NULL
  print(alpha_lambdam9)
  ###########
  # TIME -6 #
  ###########
  print(alpha_lambdam6)
  alphaWinnerm6tita<-alpha_lambdam6[1,1]
  print(paste("alphaWinnerm6tita=",alphaWinnerm6tita))
  alphaWinnerm6TOTA<-c(alphaWinnerm6TOTA,alphaWinnerm6tita)
  
  print(paste("alphaWinnerm6TOTA=",alphaWinnerm6TOTA))
  alphaWinnerm6tita<-NULL
  alpha_lambdam6<-NULL
  print(alpha_lambdam6)
  ###########
  # TIME -3 #
  ###########
  print(alpha_lambdam3)
  alphaWinnerm3tita<-alpha_lambdam3[1,1]
  print(paste("alphaWinnerm3tita=",alphaWinnerm3tita))
  alphaWinnerm3TOTA<-c(alphaWinnerm3TOTA,alphaWinnerm3tita)
  
  print(paste("alphaWinnerm3TOTA=",alphaWinnerm3TOTA))
  alphaWinnerm3tita<-NULL
  alpha_lambdam3<-NULL
  print(alpha_lambdam3)
  ##########
  # TIME 0 #
  ##########
  print(alpha_lambdam0)
  alphaWinnerm0tita<-alpha_lambdam0[1,1]
  print(paste("alphaWinnerm0tita=",alphaWinnerm0tita))
  alphaWinnerm0TOTA<-c(alphaWinnerm0TOTA,alphaWinnerm0tita)
  
  print(paste("alphaWinnerm0TOTA=",alphaWinnerm0TOTA))
  alphaWinnerm0tita<-NULL
  alpha_lambdam0<-NULL
  print(alpha_lambdam0)
}
############
# TIME -12 #
############
tabelitam12<-data.frame(table(alphaWinnerm12TOTA))
tabelitam12<-tabelitam12[order(-tabelitam12$Freq),]
alphaWinnerm12 <- as.numeric(as.character(tabelitam12[1,1]))
print(paste("alphaWinnerm12",alphaWinnerm12))  

###########
# TIME -9 #
###########
tabelitam9<-data.frame(table(alphaWinnerm9TOTA))
tabelitam9<-tabelitam9[order(-tabelitam9$Freq),]
alphaWinnerm9 <- as.numeric(as.character(tabelitam9[1,1]))
print(paste("alphaWinnerm9",alphaWinnerm9))  

###########
# TIME -6 #
###########
tabelitam6<-data.frame(table(alphaWinnerm6TOTA))
tabelitam6<-tabelitam6[order(-tabelitam6$Freq),]
alphaWinnerm6 <- as.numeric(as.character(tabelitam6[1,1]))
print(paste("alphaWinnerm6",alphaWinnerm6))  

###########
# TIME -3 #
###########
tabelitam3<-data.frame(table(alphaWinnerm3TOTA))
tabelitam3<-tabelitam3[order(-tabelitam3$Freq),]
alphaWinnerm3 <- as.numeric(as.character(tabelitam3[1,1]))
print(paste("alphaWinnerm3",alphaWinnerm3))  

##########
# TIME 0 #
##########
tabelitam0<-data.frame(table(alphaWinnerm0TOTA))
tabelitam0<-tabelitam0[order(-tabelitam0$Freq),]
alphaWinnerm0 <- as.numeric(as.character(tabelitam0[1,1]))
print(paste("alphaWinnerm0",alphaWinnerm0))  



alphaWinnerm12
alphaWinnerm9
alphaWinnerm6
alphaWinnerm3
alphaWinnerm0

###########################################################
#####     STEP 2: Iteration with the best alphas     ######
###########################################################
SelVarsm12<-SelVarsm9<-SelVarsm6<-SelVarsm3<-SelVarsm0<- NULL

for (i in 1:1000) {
  print (paste ("iteration ",i, sep=" "))
  # TIME POINT -12 #
  cfit<-cv.glmnet(as.matrix(FullarrayGEMETABDBBestModelCASESGADAIAA [,,1]),FAAbcasesGADAIAA$FirstAAb, 
                  standardize=TRUE, family="binomial", 
                  alpha = alphaWinnerm12,
                  nfolds = 5,
                  #alpha = 0.6,
                  type.measure = "auc") # cross validation with  
  
  #plot(cfit)
  #print(paste("alpha =",alpha_val, ", lambda.min =", cfit$lambda.min))
  
  CoefficientsVars<-as.matrix(coef(cfit, s = cfit$lambda.min))
  SelVarsPositions<-which(CoefficientsVars != 0)
  SelVarsm12small<-data.frame(CoefficientsVars=CoefficientsVars[which(CoefficientsVars != 0),])
  SelVarsm12smalllista<- rownames(SelVarsm12small)[-1]
  SelVarsm12<-c(SelVarsm12,SelVarsm12smalllista)
  SelVarsm12<-unique(SelVarsm12)
  # TIME POINT -9 #
  cfit<-cv.glmnet(as.matrix(FullarrayGEMETABDBBestModelCASESGADAIAA [,,2]),FAAbcasesGADAIAA$FirstAAb, 
                  standardize=TRUE, 
                  alpha = alphaWinnerm9,
                  nfolds = 5,
                  family="binomial", type.measure = "auc") # cross validation with  
  CoefficientsVars<-as.matrix(coef(cfit, s = cfit$lambda.min))
  SelVarsPositions<-which(CoefficientsVars != 0)
  SelVarsm9small<-data.frame(CoefficientsVars=CoefficientsVars[which(CoefficientsVars != 0),])
  SelVarsm9smalllista<- rownames(SelVarsm9small)[-1]
  SelVarsm9<-c(SelVarsm9,SelVarsm9smalllista)
  SelVarsm9<-unique(SelVarsm9)
  # TIME POINT -6 #
  cfit<-cv.glmnet(as.matrix(FullarrayGEMETABDBBestModelCASESGADAIAA [,,3]),FAAbcasesGADAIAA$FirstAAb,
                  standardize=TRUE, 
                  alpha = alphaWinnerm6,
                  nfolds = 5,
                  family="binomial", type.measure = "auc") # cross validation with  
  CoefficientsVars<-as.matrix(coef(cfit, s = cfit$lambda.min))
  SelVarsPositions<-which(CoefficientsVars != 0)
  SelVarsm6small<-data.frame(CoefficientsVars=CoefficientsVars[which(CoefficientsVars != 0),])
  SelVarsm6smalllista<- rownames(SelVarsm6small)[-1]
  SelVarsm6<-c(SelVarsm6,SelVarsm6smalllista)
  SelVarsm6<-unique(SelVarsm6)
  # TIME POINT -3 #
  cfit<-cv.glmnet(as.matrix(FullarrayGEMETABDBBestModelCASESGADAIAA [,,4]),FAAbcasesGADAIAA$FirstAAb,
                  standardize=TRUE, 
                  alpha = alphaWinnerm3,
                  nfolds = 5,
                  family="binomial", type.measure = "auc") # cross validation with  
  CoefficientsVars<-as.matrix(coef(cfit, s = cfit$lambda.min))
  SelVarsPositions<-which(CoefficientsVars != 0)
  SelVarsm3small<-data.frame(CoefficientsVars=CoefficientsVars[which(CoefficientsVars != 0),])
  SelVarsm3smalllista<- rownames(SelVarsm3small)[-1]
  SelVarsm3<-c(SelVarsm3,SelVarsm3smalllista)
  SelVarsm3<-unique(SelVarsm3)
  # TIME POINT 0 #
  cfit<-cv.glmnet(as.matrix(FullarrayGEMETABDBBestModelCASESGADAIAA [,,5]),FAAbcasesGADAIAA$FirstAAb,
                  standardize=TRUE, 
                  alpha = alphaWinnerm0,
                  nfolds = 5,
                  family="binomial", type.measure = "auc") # cross validation with  
  CoefficientsVars<-as.matrix(coef(cfit, s = cfit$lambda.min))
  SelVarsPositions<-which(CoefficientsVars != 0)
  SelVarsm0small<-data.frame(CoefficientsVars=CoefficientsVars[which(CoefficientsVars != 0),])
  SelVarsm0smalllista<- rownames(SelVarsm0small)[-1]
  SelVarsm0<-c(SelVarsm0,SelVarsm0smalllista)
  SelVarsm0<-unique(SelVarsm0)
}

ENSelVarswithDBGADAIAA<-unique(c(SelVarsm12,SelVarsm9,SelVarsm6,SelVarsm3,SelVarsm0))
ENSelVarswithDBGADAIAA
#setwd("/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/NetworksAnalyzes")
#write.table(ENSelVarswithDBGADAIAA,file="/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/NetworksAnalyzes/ENSelVarswithDBGADAIAA.txt")

###########################################################
###########################################################
####################     NETWORKS     #####################
###########################################################
###########################################################
# 1.- Perform ParCor to the cases only, then subset the correlation matrix and plot the network
# Time Period -12 to -9

#################################################
#############     CORRELATIONS     ##############
#################################################
# Compute correlations:
##############################
### Time Period -12 to -9 ####
##############################
Correlation9to12cases <- cor(Deltaallvarsm12m9cases, use = "pairwise")  

class(Correlation9to12cases)
isSymmetric(Correlation9to12cases)
PCorMat9to12cases <- cor2pcor(Correlation9to12cases)
class(PCorMat9to12cases)
isSymmetric(PCorMat9to12cases)
PCorMat9to12cases<-round(PCorMat9to12cases,digits = 6)
isSymmetric(PCorMat9to12cases)
colnames(PCorMat9to12cases)<-colnames(Correlation9to12cases)
rownames(PCorMat9to12cases)<-rownames(Correlation9to12cases)
dim(PCorMat9to12cases)
min(PCorMat9to12cases)
max(PCorMat9to12cases)
########################################################
# Genes= 862, Metabolites= 245, DB= 3
colgroup<-c(rep("dodgerblue2",862),rep("gold2",245), rep("deeppink1",3))
length(colgroup)==length(colnames(PCorMat9to12cases))
namegroup<-c(rep("Genes",862),rep("Metabolites",245), rep("deeppink1",3))
length(namegroup)==length(colnames(PCorMat9to12cases))

################################################
#  Subsetting the ParCor matrix to the mentioned actors
length(ENSelVarswithDBGADAIAA) # 315

PCorMat9to12cases2<-PCorMat9to12cases[rownames(PCorMat9to12cases) %in% ENSelVarswithDBGADAIAA$x,]
dim(PCorMat9to12cases2)
PCorMat9to12casesSizeEN<-PCorMat9to12cases2[,colnames(PCorMat9to12cases2) %in% ENSelVarswithDBGADAIAA$x]
dim(PCorMat9to12casesSizeEN)
PCorMat9to12casesSizeEN[1:10,1:10]
################################################
# Melting the matrix
upperTriangle = upper.tri(PCorMat9to12casesSizeEN, diag=F) #turn into a upper triangle
correlations.upperTriangle = PCorMat9to12casesSizeEN #take a copy of the original cor-mat
correlations.upperTriangle[!upperTriangle]<-NA #set everything not in upper triangle o NA
PCorMat9to12caseswithDB_melted<-na.omit(melt(correlations.upperTriangle, value.name ="Correlation")) 
dim(PCorMat9to12caseswithDB_melted)
head(PCorMat9to12caseswithDB_melted)
colnames(PCorMat9to12caseswithDB_melted)<- c("Var1", "Var2","Correlations")
PCorMat9to12cases_meltedGADAIAA<-PCorMat9to12caseswithDB_melted
#write.table(PCorMat9to12cases_meltedGADAIAA,file="/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/NetworksAnalyzes/FAAb/PCorMat9to12cases_meltedGADAIAA.txt", sep="\t")

##########################################################################################
# Three vectors to assign the category GE (1),MEtab (2) or DB (3)
Deltaallvarsm12m9cases[1:10,1:10]
dim(Deltaallvarsm12m9cases)
vectitodeGE<-colnames(Deltaallvarsm12m9cases)[1:862]
vectitodeMetab<-colnames(Deltaallvarsm12m9cases)[863:1107]
vectitodeDB<-colnames(Deltaallvarsm12m9cases)[1108:1110]
#########################################################################
# Creation of the Node Table

NodeTable9to12caseswithDB<-PCorMat9to12caseswithDB_melted[,-3]

NodeTable9to12caseswithDB2<- data.frame(Var=c(as.vector(NodeTable9to12caseswithDB[,1]),as.vector(NodeTable9to12caseswithDB[,2])))
length(as.vector(NodeTable9to12caseswithDB[,1])); length(as.vector(NodeTable9to12caseswithDB[,2]));dim(NodeTable9to12caseswithDB2)
head(NodeTable9to12caseswithDB2)
#########################################################################
# Creation of clasifier of what type of dataset, the feature is
NodeTable9to12caseswithDB2$Omic<-""
head(NodeTable9to12caseswithDB2)
NodeTable9to12caseswithDB2[1:50,]

NodeTable9to12caseswithDB2$Omic [NodeTable9to12caseswithDB2$Var %in% vectitodeGE] <- 1
NodeTable9to12caseswithDB2$Omic [NodeTable9to12caseswithDB2$Var %in% vectitodeMetab] <- 2
NodeTable9to12caseswithDB2$Omic [NodeTable9to12caseswithDB2$Var %in% vectitodeDB] <- 3
table(NodeTable9to12caseswithDB2$Omic)
#sum(table(NodeTable9to12caseswithDB2$Omic))
NodeTable9to12caseswithDB2[850:1065,]
dim(NodeTable9to12caseswithDB2)

##########################################################################################
# Calculate the weight of each variable

thrs<-0.7
table(abs(PCorMat9to12caseswithDB_melted$Correlations)>=thrs) # 428 true
dim(PCorMat9to12caseswithDB_melted)
PCorMat9to12caseswithDB_melted_CO0.7<-PCorMat9to12caseswithDB_melted[abs(PCorMat9to12caseswithDB_melted$Correlations)>=thrs,]
dim(PCorMat9to12caseswithDB_melted_CO0.7)
head(PCorMat9to12caseswithDB_melted_CO0.7)
PCorMat9to12caseswithDB_melted_CO0.72<- data.frame(Var=c(as.vector(PCorMat9to12caseswithDB_melted_CO0.7[,1]),as.vector(PCorMat9to12caseswithDB_melted_CO0.7[,2])))
length(as.vector(PCorMat9to12caseswithDB_melted_CO0.7[,2]));dim(PCorMat9to12caseswithDB_melted_CO0.72)
dim(unique(PCorMat9to12caseswithDB_melted_CO0.72))

tulo<-as.data.frame(sort(table(PCorMat9to12caseswithDB_melted_CO0.72$Var),decreasing = TRUE))

NodeTable9to12caseswithDB2$Weight2<-tulo$Freq [match(NodeTable9to12caseswithDB2$Var,tulo$Var1)]
NodeTable9to12caseswithDB2[1:50,]

NodeTable9to12caseswithDB2[NodeTable9to12caseswithDB2$Var=="KCNG4",]
dim(NodeTable9to12caseswithDB2[NodeTable9to12caseswithDB2$Var=="KCNG4",])
##########################################################################################
# Adjusting of Weights to remove NA's
NodeTable9to12caseswithDB2$Weight2[is.na(NodeTable9to12caseswithDB2$Weight2)]<-0
NodeTable9to12caseswithDB2$Weight<-NodeTable9to12caseswithDB2$Weight2 +1

##########################################################################################
# Creating a abbreviated name to plot
NodeTable9to12caseswithDB2$Nameito<- abbreviate(NodeTable9to12caseswithDB2$Var, 7, method = "both")
NodeTable9to12caseswithDB2[NodeTable9to12caseswithDB2$Omic==2,]
##########################################################################################
# Creating Expression Level column to know if the feature is up or downregulated
Deltaallvarsm12m9cases[1:10,1:10]
dim(Deltaallvarsm12m9cases)
vectitodebehaviormedio<-data.frame(ExpressionLevel=colMeans(Deltaallvarsm12m9cases))
head(vectitodebehaviormedio)
dim(vectitodebehaviormedio)


NodeTable9to12caseswithDB2$ExpressionLevel<-vectitodebehaviormedio$ExpressionLevel[match(NodeTable9to12caseswithDB2$Var,rownames(vectitodebehaviormedio))]

##########################################################################################
# Definitive table #
NodeTable9to12caseswithDB<-unique (NodeTable9to12caseswithDB2)
length(unique (NodeTable9to12caseswithDB2$Var)) # 406
dim(NodeTable9to12caseswithDB)
#write.table(NodeTable9to12caseswithDB,file="/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/NetworksAnalyzes/withDB/NodeTable9to12caseswithDB", sep="\t")
NodeTable9to12caseswithDB <- read.delim("~/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/NetworksAnalyzes/withDB/NodeTable9to12caseswithDB", sep="\t")
NodeTable9to12caseswithDB
dim(AgeGroupdummycasesfrom136)
Cases12to9<- FullarrayGEMETABDBBestModel[,,2][rownames(FullarrayGEMETABDBBestModel) %in% rownames(AgeGroupdummycasesfrom136),]
dim(Cases12to9)
Cases12to9SMALL<-Cases12to9[,colnames(Cases12to9) %in% NodeTable9to12caseswithDB$Var]
dim(Cases12to9SMALL)
vectitodeExpression<-data.frame(colMeans(Cases12to9SMALL))
rownames(NodeTable9to12caseswithDB)<-NodeTable9to12caseswithDB$Var

table12to9<-merge(NodeTable9to12caseswithDB,vectitodeExpression, by=0)
colnames(table12to9)[9]<-"ExpressionTF"
#write.table(table12to9,file="/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/NetworksAnalyzes/withDB/table12to9.txt", sep="\t")

##########################################################################################
# To Cytoscape
# 1) give a tab at the beginning of the first row (the colnames row)
# 2) Erase all "" in both datasets
# 3) Import PCorMat9to12caseswithDB_melted_CO0.7 as Network specifying...
#     3.1) In Network Collection: Create a new network collection
#     3.2) Then hit select none and then...
#     3.3) Assign Var 1 as the Source Node
#     3.4) Assign Var2 as the target Node
#     3.5) Assign Correlations as the Edge Attribute
# 4) To the node table add NodeTable9to12caseswithDB as follows:
#     4.1) File -> Import -> Table -> File -> NodeTable9to12caseswithDB
#     4.2) In Where to Import Table Data: Select To selected networks only
#     4.3) In Network List: Select PCorMat9to12caseswithDB_melted_CO0.7.txt
#     4.4) In Import Data as: Select Node Table Columns
#     4.5) Then hit select none and then...
#     4.6) Assign Var as the Key
#     4.7) Assign Omic, Weight and Nameito as Attribute
# 5) Go to Tools-> NetworkAnalyzer -> Network Analysis -> Analyze network->

##############################
### Time Period -9 to -6 #####
##############################

Correlation9to6cases <- cor(Deltaallvarsm9m6cases, use = "pairwise")  

class(Correlation9to6cases)
isSymmetric(Correlation9to6cases)
PCorMat9to6cases <- cor2pcor(Correlation9to6cases)
class(PCorMat9to6cases)
isSymmetric(PCorMat9to6cases)
PCorMat9to6cases<-round(PCorMat9to6cases,digits = 6)
isSymmetric(PCorMat9to6cases)
colnames(PCorMat9to6cases)<-colnames(Correlation9to6cases)
rownames(PCorMat9to6cases)<-rownames(Correlation9to6cases)
dim(PCorMat9to6cases)
min(PCorMat9to6cases)
max(PCorMat9to6cases)
########################################################
# Genes= 862, Metabolites= 245, DB= 3
colgroup<-c(rep("dodgerblue2",862),rep("gold2",245), rep("deeppink1",3))
length(colgroup)==length(colnames(PCorMat9to12cases))
namegroup<-c(rep("Genes",862),rep("Metabolites",245), rep("deeppink1",3))
length(namegroup)==length(colnames(PCorMat9to12cases))

################################################
#  Subsetting the ParCor matrix to the mentioned actors
length(ENSelVarswithDB) # 406

PCorMat9to6cases2<-PCorMat9to6cases[rownames(PCorMat9to6cases) %in% ENSelVarswithDB,]
dim(PCorMat9to6cases2)
PCorMat9to6casesSizeEN<-PCorMat9to6cases2[,colnames(PCorMat9to6cases2) %in% ENSelVarswithDB]
dim(PCorMat9to6casesSizeEN)
PCorMat9to6casesSizeEN[1:10,1:10]
################################################
# Melting the matrix
upperTriangle = upper.tri(PCorMat9to6casesSizeEN, diag=F) #turn into a upper triangle
correlations.upperTriangle = PCorMat9to6casesSizeEN #take a copy of the original cor-mat
correlations.upperTriangle[!upperTriangle]<-NA #set everything not in upper triangle o NA
PCorMat9to6caseswithDB_meltedapproach1<-na.omit(melt(correlations.upperTriangle, value.name ="Correlation")) 
dim(PCorMat9to6caseswithDB_meltedapproach1)
head(PCorMat9to6caseswithDB_meltedapproach1)
colnames(PCorMat9to6caseswithDB_meltedapproach1)<- c("Var1", "Var2","Correlations")
#write.table(PCorMat9to6caseswithDB_meltedapproach1,file="/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/NetworksAnalyzes/withDB/PCorMat9to6caseswithDB_meltedapproach1.txt", sep="\t")
##########################################################################################
# Three vectors to assign the category GE (1),MEtab (2) or DB (3)
Deltaallvarsm12m9cases[1:10,1:10]
dim(Deltaallvarsm12m9cases)
vectitodeGE<-colnames(Deltaallvarsm12m9cases)[1:862]
vectitodeMetab<-colnames(Deltaallvarsm12m9cases)[863:1107]
vectitodeDB<-colnames(Deltaallvarsm12m9cases)[1108:1110]
#########################################################################
# Creation of the Node Table

NodeTable9to6caseswithDB<-PCorMat9to6caseswithDB_meltedapproach1[,-3]

NodeTable9to6caseswithDB2<- data.frame(Var=c(as.vector(NodeTable9to6caseswithDB[,1]),as.vector(NodeTable9to6caseswithDB[,2])))
length(as.vector(NodeTable9to6caseswithDB[,1])); length(as.vector(NodeTable9to6caseswithDB[,2]));dim(NodeTable9to6caseswithDB2)
head(NodeTable9to6caseswithDB2)
#########################################################################
# Creation of clasifier of what type of dataset, the feature is
NodeTable9to6caseswithDB2$Omic<-""
head(NodeTable9to6caseswithDB2)
NodeTable9to6caseswithDB2[1:50,]

NodeTable9to6caseswithDB2$Omic [NodeTable9to6caseswithDB2$Var %in% vectitodeGE] <- 1
NodeTable9to6caseswithDB2$Omic [NodeTable9to6caseswithDB2$Var %in% vectitodeMetab] <- 2
NodeTable9to6caseswithDB2$Omic [NodeTable9to6caseswithDB2$Var %in% vectitodeDB] <- 3
table(NodeTable9to6caseswithDB2$Omic)
#sum(table(NodeTable9to12caseswithDB2$Omic))
NodeTable9to6caseswithDB2[850:1065,]
dim(NodeTable9to6caseswithDB2)

##########################################################################################
# Calculate the weight of each variable

thrs<-0.7
table(abs(PCorMat9to6caseswithDB_meltedapproach1$Correlations)>=thrs) # 428 true
dim(PCorMat9to6caseswithDB_meltedapproach1)
PCorMat9to6caseswithDB_melted_CO0.7<-PCorMat9to6caseswithDB_meltedapproach1[abs(PCorMat9to6caseswithDB_meltedapproach1$Correlations)>=thrs,]
dim(PCorMat9to6caseswithDB_melted_CO0.7)
head(PCorMat9to6caseswithDB_melted_CO0.7)
PCorMat9to6caseswithDB_melted_CO0.72<- data.frame(Var=c(as.vector(PCorMat9to6caseswithDB_melted_CO0.7[,1]),as.vector(PCorMat9to6caseswithDB_melted_CO0.7[,2])))
length(as.vector(PCorMat9to6caseswithDB_melted_CO0.7[,2]));dim(PCorMat9to6caseswithDB_melted_CO0.72)
dim(unique(PCorMat9to6caseswithDB_melted_CO0.72))

tulo<-as.data.frame(sort(table(PCorMat9to6caseswithDB_melted_CO0.72$Var),decreasing = TRUE))

NodeTable9to6caseswithDB2$Weight2<-tulo$Freq [match(NodeTable9to6caseswithDB2$Var,tulo$Var1)]
NodeTable9to6caseswithDB2[1:50,]

NodeTable9to6caseswithDB2[NodeTable9to6caseswithDB2$Var=="KCNG4",]
dim(NodeTable9to6caseswithDB2[NodeTable9to6caseswithDB2$Var=="KCNG4",])
##########################################################################################
# Adjusting of Weights to remove NA's
NodeTable9to6caseswithDB2$Weight2[is.na(NodeTable9to6caseswithDB2$Weight2)]<-0
NodeTable9to6caseswithDB2$Weight<-NodeTable9to6caseswithDB2$Weight2 +1

##########################################################################################
# Creating a abbreviated name to plot
NodeTable9to6caseswithDB2$Nameito<- abbreviate(NodeTable9to6caseswithDB2$Var, 7, method = "both")
NodeTable9to6caseswithDB2[NodeTable9to6caseswithDB2$Omic==2,]
##########################################################################################
# Creating Expression Level column to know if the feature is up or downregulated
Deltaallvarsm9m6cases[1:10,1:10]
dim(Deltaallvarsm9m6cases)
vectitodebehaviormedio<-data.frame(ExpressionLevel=colMeans(Deltaallvarsm9m6cases))
head(vectitodebehaviormedio)
dim(vectitodebehaviormedio)


NodeTable9to6caseswithDB2$ExpressionLevel<-vectitodebehaviormedio$ExpressionLevel[match(NodeTable9to6caseswithDB2$Var,rownames(vectitodebehaviormedio))]

##########################################################################################
# Definitive table #
load("/media/data/leobalzano/ScriptsForTEDDY/Data/FullarrayGEMETABDBBestModel.RData")
AgeGroupdummycasesfrom136 <- read.delim("/media/data/leobalzano/ScriptsForTEDDY/Data/AgeGroupdummycasesfrom136.txt", row.names=1)
NodeTable9to6caseswithDB<-unique (NodeTable9to6caseswithDB2)
length(unique (NodeTable9to6caseswithDB2$Var)) # 386
dim(NodeTable9to6caseswithDB)
#write.table(NodeTable9to6caseswithDB,file="/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/NetworksAnalyzes/withDB/NodeTable9to6caseswithDB.txt", sep="\t")
NodeTable9to6caseswithDB<-read.delim(file="/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/NetworksAnalyzes/withDB/NodeTable9to6caseswithDB.txt", sep="\t")
NodeTable9to6caseswithDB
dim(AgeGroupdummycasesfrom136)
Cases9to6<- FullarrayGEMETABDBBestModel[,,3][rownames(FullarrayGEMETABDBBestModel) %in% rownames(AgeGroupdummycasesfrom136),]
dim(Cases9to6)
Cases9to6SMALL<-Cases9to6[,colnames(Cases9to6) %in% NodeTable9to6caseswithDB$Var]
dim(Cases9to6SMALL)
vectitodeExpression<-data.frame(colMeans(Cases9to6SMALL))
rownames(NodeTable9to6caseswithDB)<-NodeTable9to6caseswithDB$Var

table9to6<-merge(NodeTable9to6caseswithDB,vectitodeExpression, by=0)
colnames(table9to6)[9]<-"ExpressionTF"
#write.table(table9to6,file="/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/NetworksAnalyzes/withDB/table9to6.txt", sep="\t")
##########################################################################################
# To Cytoscape
# 1) give a tab at the beginning of the first row (the colnames row)
# 2) Erase all "" in both datasets
# 3) Import PCorMat9to12caseswithDB_melted_CO0.7 as Network specifying...
#     3.1) In Network Collection: Create a new network collection
#     3.2) Then hit select none and then...
#     3.3) Assign Var 1 as the Source Node
#     3.4) Assign Var2 as the target Node
#     3.5) Assign Correlations as the Edge Attribute
# 4) To the node table add NodeTable9to12caseswithDB as follows:
#     4.1) File -> Import -> Table -> File -> NodeTable9to12caseswithDB
#     4.2) In Where to Import Table Data: Select To selected networks only
#     4.3) In Network List: Select PCorMat9to12caseswithDB_melted_CO0.7.txt
#     4.4) In Import Data as: Select Node Table Columns
#     4.5) Then hit select none and then...
#     4.6) Assign Var as the Key
#     4.7) Assign Omic, Weight and Nameito as Attribute
# 5) Go to Tools-> NetworkAnalyzer -> Network Analysis -> Analyze network->

##############################
### Time Period -6 to -3 #####
##############################

Correlation6to3cases <- cor(Deltaallvarsm6m3cases, use = "pairwise")  

class(Correlation6to3cases)
isSymmetric(Correlation6to3cases)
PCorMat6to3cases <- cor2pcor(Correlation6to3cases)
class(PCorMat6to3cases)
isSymmetric(PCorMat6to3cases)
PCorMat6to3cases<-round(PCorMat6to3cases,digits = 6)
isSymmetric(PCorMat6to3cases)
colnames(PCorMat6to3cases)<-colnames(Correlation6to3cases)
rownames(PCorMat6to3cases)<-rownames(Correlation6to3cases)
dim(PCorMat6to3cases)
min(PCorMat6to3cases)
max(PCorMat6to3cases)
PCorMat6to3cases[1:10,1:10]

########################################################
# Genes= 862, Metabolites= 245, DB= 3
colgroup<-c(rep("dodgerblue2",862),rep("gold2",245), rep("deeppink1",3))
length(colgroup)==length(colnames(PCorMat9to12cases))
namegroup<-c(rep("Genes",862),rep("Metabolites",245), rep("deeppink1",3))
length(namegroup)==length(colnames(PCorMat9to12cases))

################################################
#  Subsetting the ParCor matrix to the mentioned actors
length(ENSelVarswithDB) # 386

PCorMat6to3cases2<-PCorMat6to3cases[rownames(PCorMat6to3cases) %in% ENSelVarswithDB,]
dim(PCorMat6to3cases2)
PCorMat6to3casesSizeEN<-PCorMat6to3cases2[,colnames(PCorMat6to3cases2) %in% ENSelVarswithDB]
dim(PCorMat6to3casesSizeEN)
PCorMat6to3casesSizeEN[1:10,1:10]
################################################
# Melting the matrix
upperTriangle = upper.tri(PCorMat6to3casesSizeEN, diag=F) #turn into a upper triangle
correlations.upperTriangle = PCorMat6to3casesSizeEN #take a copy of the original cor-mat
correlations.upperTriangle[!upperTriangle]<-NA #set everything not in upper triangle o NA
PCorMat6to3caseswithDB_meltedapproach1<-na.omit(melt(correlations.upperTriangle, value.name ="Correlation")) 
dim(PCorMat6to3caseswithDB_meltedapproach1)
head(PCorMat6to3caseswithDB_meltedapproach1)
colnames(PCorMat6to3caseswithDB_meltedapproach1)<- c("Var1", "Var2","Correlations")
#write.table(PCorMat6to3caseswithDB_meltedapproach1,file="/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/NetworksAnalyzes/withDB/PCorMat6to3caseswithDB_meltedapproach1.txt", sep="\t")
##########################################################################################
# Three vectors to assign the category GE (1),MEtab (2) or DB (3)
Deltaallvarsm12m9cases[1:10,1:10]
dim(Deltaallvarsm12m9cases)
vectitodeGE<-colnames(Deltaallvarsm12m9cases)[1:862]
vectitodeMetab<-colnames(Deltaallvarsm12m9cases)[863:1107]
vectitodeDB<-colnames(Deltaallvarsm12m9cases)[1108:1110]
#########################################################################
# Creation of the Node Table
NodeTable6to3caseswithDB<-PCorMat6to3caseswithDB_meltedapproach1[,-3]

NodeTable6to3caseswithDB2<- data.frame(Var=c(as.vector(NodeTable6to3caseswithDB[,1]),as.vector(NodeTable6to3caseswithDB[,2])))
length(as.vector(NodeTable6to3caseswithDB[,1])); length(as.vector(NodeTable6to3caseswithDB[,2]));dim(NodeTable6to3caseswithDB2)
head(NodeTable6to3caseswithDB2)
#########################################################################
# Creation of clasifier of what type of dataset, the feature is
NodeTable6to3caseswithDB2$Omic<-""
head(NodeTable6to3caseswithDB2)
NodeTable6to3caseswithDB2[1:50,]

NodeTable6to3caseswithDB2$Omic [NodeTable6to3caseswithDB2$Var %in% vectitodeGE] <- 1
NodeTable6to3caseswithDB2$Omic [NodeTable6to3caseswithDB2$Var %in% vectitodeMetab] <- 2
NodeTable6to3caseswithDB2$Omic [NodeTable6to3caseswithDB2$Var %in% vectitodeDB] <- 3
table(NodeTable6to3caseswithDB2$Omic)
#sum(table(NodeTable6to3caseswithDB2$Omic))
NodeTable6to3caseswithDB2[850:1065,]
dim(NodeTable6to3caseswithDB2)

##########################################################################################
# Calculate the weight of each variable

thrs<-0.7
table(abs(PCorMat6to3caseswithDB_meltedapproach1$Correlations)>=thrs) # 428 true
dim(PCorMat6to3caseswithDB_meltedapproach1)
PCorMat6to3caseswithDB_melted_CO0.7<-PCorMat6to3caseswithDB_meltedapproach1[abs(PCorMat6to3caseswithDB_meltedapproach1$Correlations)>=thrs,]
dim(PCorMat6to3caseswithDB_melted_CO0.7)
head(PCorMat6to3caseswithDB_melted_CO0.7)
PCorMat6to3caseswithDB_melted_CO0.72<- data.frame(Var=c(as.vector(PCorMat6to3caseswithDB_melted_CO0.7[,1]),as.vector(PCorMat6to3caseswithDB_melted_CO0.7[,2])))
length(as.vector(PCorMat6to3caseswithDB_melted_CO0.7[,2]));dim(PCorMat6to3caseswithDB_melted_CO0.72)
dim(unique(PCorMat6to3caseswithDB_melted_CO0.72))

tulo<-as.data.frame(sort(table(PCorMat6to3caseswithDB_melted_CO0.72$Var),decreasing = TRUE))

NodeTable6to3caseswithDB2$Weight2<-tulo$Freq [match(NodeTable6to3caseswithDB2$Var,tulo$Var1)]
NodeTable6to3caseswithDB2[1:50,]

NodeTable6to3caseswithDB2[NodeTable6to3caseswithDB2$Var=="KCNG4",]
dim(NodeTable6to3caseswithDB2[NodeTable6to3caseswithDB2$Var=="KCNG4",])
##########################################################################################
# Adjusting of Weights to remove NA's
NodeTable6to3caseswithDB2$Weight2[is.na(NodeTable6to3caseswithDB2$Weight2)]<-0
NodeTable6to3caseswithDB2$Weight<-NodeTable6to3caseswithDB2$Weight2 +1

##########################################################################################
# Creating a abbreviated name to plot
NodeTable6to3caseswithDB2$Nameito<- abbreviate(NodeTable6to3caseswithDB2$Var, 7, method = "both")
NodeTable6to3caseswithDB2[NodeTable6to3caseswithDB2$Omic==2,]
##########################################################################################
# Creating Expression Level column to know if the feature is up or downregulated
Deltaallvarsm6m3cases[1:10,1:10]
dim(Deltaallvarsm6m3cases)
vectitodebehaviormedio<-data.frame(ExpressionLevel=colMeans(Deltaallvarsm6m3cases))
head(vectitodebehaviormedio)
dim(vectitodebehaviormedio)


NodeTable6to3caseswithDB2$ExpressionLevel<-vectitodebehaviormedio$ExpressionLevel[match(NodeTable6to3caseswithDB2$Var,rownames(vectitodebehaviormedio))]

##########################################################################################
# Definitive table #
NodeTable6to3caseswithDB<-unique (NodeTable6to3caseswithDB2)
length(unique (NodeTable6to3caseswithDB2$Var)) # 386
dim(NodeTable6to3caseswithDB)
#write.table(NodeTable6to3caseswithDB,file="/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/NetworksAnalyzes/withDB/NodeTable6to3caseswithDB.txt", sep="\t")
NodeTable6to3caseswithDB <- read.delim("~/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/NetworksAnalyzes/withDB/NodeTable6to3caseswithDB.txt", sep="\t")
NodeTable6to3caseswithDB
dim(AgeGroupdummycasesfrom136)
Cases6to3<- FullarrayGEMETABDBBestModel[,,4][rownames(FullarrayGEMETABDBBestModel) %in% rownames(AgeGroupdummycasesfrom136),]
dim(Cases6to3)
Cases6to3SMALL<-Cases6to3[,colnames(Cases6to3) %in% NodeTable6to3caseswithDB$Var]
dim(Cases6to3SMALL)
vectitodeExpression<-data.frame(colMeans(Cases6to3SMALL))
rownames(NodeTable6to3caseswithDB)<-NodeTable6to3caseswithDB$Var

table6to3<-merge(NodeTable6to3caseswithDB,vectitodeExpression, by=0)
colnames(table6to3)[9]<-"ExpressionTF"
#write.table(table6to3,file="/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/NetworksAnalyzes/withDB/table6to3.txt", sep="\t")

##########################################################################################
# To Cytoscape
# 1) give a tab at the beginning of the first row (the colnames row)
# 2) Erase all "" in both datasets
# 3) Import PCorMat9to12caseswithDB_melted_CO0.7 as Network specifying...
#     3.1) In Network Collection: Create a new network collection
#     3.2) Then hit select none and then...
#     3.3) Assign Var 1 as the Source Node
#     3.4) Assign Var2 as the target Node
#     3.5) Assign Correlations as the Edge Attribute
# 4) To the node table add NodeTable9to12caseswithDB as follows:
#     4.1) File -> Import -> Table -> File -> NodeTable9to12caseswithDB
#     4.2) In Where to Import Table Data: Select To selected networks only
#     4.3) In Network List: Select PCorMat9to12caseswithDB_melted_CO0.7.txt
#     4.4) In Import Data as: Select Node Table Columns
#     4.5) Then hit select none and then...
#     4.6) Assign Var as the Key
#     4.7) Assign Omic, Weight and Nameito as Attribute
# 5) Go to Tools-> NetworkAnalyzer -> Network Analysis -> Analyze network->

##############################
### Time Period -3 to 0 ######
##############################

Correlation3to0cases <- cor(Deltaallvarsm3m0cases, use = "pairwise")  

class(Correlation3to0cases)
isSymmetric(Correlation3to0cases)
PCorMat3to0cases <- cor2pcor(Correlation3to0cases)
class(PCorMat3to0cases)
isSymmetric(PCorMat3to0cases)
PCorMat3to0cases<-round(PCorMat3to0cases,digits = 6)
isSymmetric(PCorMat3to0cases)
colnames(PCorMat3to0cases)<-colnames(Correlation6to3cases)
rownames(PCorMat3to0cases)<-rownames(Correlation6to3cases)
dim(PCorMat3to0cases)
min(PCorMat3to0cases)
max(PCorMat3to0cases)
########################################################
# Genes= 862, Metabolites= 245, DB= 3
colgroup<-c(rep("dodgerblue2",862),rep("gold2",245), rep("deeppink1",3))
length(colgroup)==length(colnames(PCorMat9to12cases))
namegroup<-c(rep("Genes",862),rep("Metabolites",245), rep("deeppink1",3))
length(namegroup)==length(colnames(PCorMat9to12cases))

################################################
#  Subsetting the ParCor matrix to the mentioned actors
length(ENSelVarswithDB) # 386

PCorMat3to0cases2<-PCorMat3to0cases[rownames(PCorMat3to0cases) %in% ENSelVarswithDB,]
dim(PCorMat3to0cases2)
PCorMat3to0casesSizeEN<-PCorMat3to0cases2[,colnames(PCorMat3to0cases2) %in% ENSelVarswithDB]
dim(PCorMat3to0casesSizeEN)
PCorMat3to0casesSizeEN[1:10,1:10]
################################################
# Melting the matrix
upperTriangle = upper.tri(PCorMat3to0casesSizeEN, diag=F) #turn into a upper triangle
correlations.upperTriangle = PCorMat3to0casesSizeEN #take a copy of the original cor-mat
correlations.upperTriangle[!upperTriangle]<-NA #set everything not in upper triangle o NA
PCorMat3to0caseswithDB_meltedapproach1<-na.omit(melt(correlations.upperTriangle, value.name ="Correlation")) 
dim(PCorMat3to0caseswithDB_meltedapproach1)
head(PCorMat3to0caseswithDB_meltedapproach1)
colnames(PCorMat3to0caseswithDB_meltedapproach1)<- c("Var1", "Var2","Correlations")
#write.table(PCorMat3to0caseswithDB_meltedapproach1,file="/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/NetworksAnalyzes/withDB/PCorMat3to0caseswithDB_meltedapproach1.txt", sep="\t")
##########################################################################################
# Three vectors to assign the category GE (1),MEtab (2) or DB (3)
Deltaallvarsm12m9cases[1:10,1:10]
dim(Deltaallvarsm12m9cases)
vectitodeGE<-colnames(Deltaallvarsm12m9cases)[1:862]
vectitodeMetab<-colnames(Deltaallvarsm12m9cases)[863:1107]
vectitodeDB<-colnames(Deltaallvarsm12m9cases)[1108:1110]
#########################################################################
# Creation of the Node Table
NodeTable3to0caseswithDB<-PCorMat3to0caseswithDB_meltedapproach1[,-3]

NodeTable3to0caseswithDB2<- data.frame(Var=c(as.vector(NodeTable3to0caseswithDB[,1]),as.vector(NodeTable3to0caseswithDB[,2])))
length(as.vector(NodeTable3to0caseswithDB[,1])); length(as.vector(NodeTable3to0caseswithDB[,2]));dim(NodeTable3to0caseswithDB2)
head(NodeTable3to0caseswithDB2)
#########################################################################
# Creation of clasifier of what type of dataset, the feature is
NodeTable3to0caseswithDB2$Omic<-""
head(NodeTable3to0caseswithDB2)
NodeTable3to0caseswithDB2[1:50,]

NodeTable3to0caseswithDB2$Omic [NodeTable3to0caseswithDB2$Var %in% vectitodeGE] <- 1
NodeTable3to0caseswithDB2$Omic [NodeTable3to0caseswithDB2$Var %in% vectitodeMetab] <- 2
NodeTable3to0caseswithDB2$Omic [NodeTable3to0caseswithDB2$Var %in% vectitodeDB] <- 3
table(NodeTable3to0caseswithDB2$Omic)
#sum(table(NodeTable3to0caseswithDB2$Omic))
NodeTable3to0caseswithDB2[850:1065,]
dim(NodeTable3to0caseswithDB2)

##########################################################################################
# Calculate the weight of each variable

thrs<-0.7
table(abs(PCorMat3to0caseswithDB_meltedapproach1$Correlations)>=thrs) # 428 true
dim(PCorMat3to0caseswithDB_meltedapproach1)
PCorMat3to0caseswithDB_melted_CO0.7<-PCorMat3to0caseswithDB_meltedapproach1[abs(PCorMat3to0caseswithDB_meltedapproach1$Correlations)>=thrs,]
dim(PCorMat3to0caseswithDB_melted_CO0.7)
head(PCorMat3to0caseswithDB_melted_CO0.7)
PCorMat3to0caseswithDB_melted_CO0.72<- data.frame(Var=c(as.vector(PCorMat3to0caseswithDB_melted_CO0.7[,1]),as.vector(PCorMat3to0caseswithDB_melted_CO0.7[,2])))
length(as.vector(PCorMat3to0caseswithDB_melted_CO0.7[,2]));dim(PCorMat3to0caseswithDB_melted_CO0.72)
dim(unique(PCorMat3to0caseswithDB_melted_CO0.72))

tulo<-as.data.frame(sort(table(PCorMat3to0caseswithDB_melted_CO0.72$Var),decreasing = TRUE))

NodeTable3to0caseswithDB2$Weight2<-tulo$Freq [match(NodeTable3to0caseswithDB2$Var,tulo$Var1)]
NodeTable3to0caseswithDB2[1:50,]

NodeTable3to0caseswithDB2[NodeTable3to0caseswithDB2$Var=="KCNG4",]
dim(NodeTable3to0caseswithDB2[NodeTable3to0caseswithDB2$Var=="KCNG4",])
##########################################################################################
# Adjusting of Weights to remove NA's
NodeTable3to0caseswithDB2$Weight2[is.na(NodeTable3to0caseswithDB2$Weight2)]<-0
NodeTable3to0caseswithDB2$Weight<-NodeTable3to0caseswithDB2$Weight2 +1

##########################################################################################
# Creating a abbreviated name to plot
NodeTable3to0caseswithDB2$Nameito<- abbreviate(NodeTable3to0caseswithDB2$Var, 7, method = "both")
NodeTable3to0caseswithDB2[NodeTable3to0caseswithDB2$Omic==2,]
##########################################################################################
# Creating Expression Level column to know if the feature is up or downregulated
Deltaallvarsm3m0cases[1:10,1:10]
dim(Deltaallvarsm3m0cases)
vectitodebehaviormedio<-data.frame(ExpressionLevel=colMeans(Deltaallvarsm3m0cases))
head(vectitodebehaviormedio)
dim(vectitodebehaviormedio)


NodeTable3to0caseswithDB2$ExpressionLevel<-vectitodebehaviormedio$ExpressionLevel[match(NodeTable3to0caseswithDB2$Var,rownames(vectitodebehaviormedio))]

##########################################################################################
# Definitive table #
NodeTable3to0caseswithDB<-unique (NodeTable3to0caseswithDB2)
length(unique (NodeTable3to0caseswithDB2$Var)) # 386
dim(NodeTable3to0caseswithDB)
#write.table(NodeTable3to0caseswithDB,file="/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/NetworksAnalyzes/withDB/NodeTable3to0caseswithDB.txt", sep="\t")
NodeTable3to0caseswithDB <- read.delim("~/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/NetworksAnalyzes/withDB/NodeTable3to0caseswithDB.txt", sep="\t")
NodeTable3to0caseswithDB
dim(AgeGroupdummycasesfrom136)
Cases3to0<- FullarrayGEMETABDBBestModel[,,5][rownames(FullarrayGEMETABDBBestModel) %in% rownames(AgeGroupdummycasesfrom136),]
dim(Cases3to0)
Cases3to0SMALL<-Cases3to0[,colnames(Cases3to0) %in% NodeTable3to0caseswithDB$Var]
dim(Cases3to0SMALL)
vectitodeExpression<-data.frame(colMeans(Cases3to0SMALL))
rownames(NodeTable3to0caseswithDB)<-NodeTable3to0caseswithDB$Var

table3to0<-merge(NodeTable3to0caseswithDB,vectitodeExpression, by=0)
colnames(table3to0)[9]<-"ExpressionTF"
#write.table(table3to0,file="/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/NetworksAnalyzes/withDB/table3to0.txt", sep="\t")



##########################################################################################
# To Cytoscape
# 1) give a tab at the beginning of the first row (the colnames row)
# 2) Erase all "" in both datasets
# 3) Import PCorMat9to12caseswithDB_melted_CO0.7 as Network specifying...
#     3.1) In Network Collection: Create a new network collection
#     3.2) Then hit select none and then...
#     3.3) Assign Var 1 as the Source Node
#     3.4) Assign Var2 as the target Node
#     3.5) Assign Correlations as the Edge Attribute
# 4) To the node table add NodeTable9to12caseswithDB as follows:
#     4.1) File -> Import -> Table -> File -> NodeTable9to12caseswithDB
#     4.2) In Where to Import Table Data: Select To selected networks only
#     4.3) In Network List: Select PCorMat9to12caseswithDB_melted_CO0.7.txt
#     4.4) In Import Data as: Select Node Table Columns
#     4.5) Then hit select none and then...
#     4.6) Assign Var as the Key
#     4.7) Assign Omic, Weight, Nameito and ExpressionLevel as Attribute
# 5) Go to Select, select from -0.7 to 0.7 and hide nodes and edges
# 6) Go to Layout -> yFiles Layout -> Organic

##################################################################################################
# Plotting proportions of the different classes of features

genesNetworks<-c(139,91,96,84)
metabsNetworks<-c(47,64,65,60)
DBNetworks<-c(3,2,3,2)

tabelinha<-rbind(genesNetworks,metabsNetworks,DBNetworks)
colnames(tabelinha)<-c("-12/-9","-9/-6","-6/-3","-3/0")


tabelinha
tabelinha2<-tabelinha
ttabelinha<-t(tabelinha)
melted<-melt(tabelinha)

tabelinha2<-rbind(tabelinha2,colSums(tabelinha2))


proportions<-prop.table(tabelinha)*100
meltedproportions<-melt(proportions)



library(ggplot2)
pdf("NumberofNodesperTime-periods.pdf")
ggplot(melted, aes(x=factor(melted$X2, level= c("-12/-9","-9/-6","-6/-3","-3/0")), y=melted$value, group=melted$X1)) + 
  labs(title="Number of Nodes per Time-periods",x="Time periods", y = "Number of Nodes",legend_title="Feature Class") +
  geom_line(aes(color=melted$X1)) + geom_point(aes(color=melted$X1))#+theme_classic()
dev.off()

#
ggplot(meltedproportions, aes(x=factor(meltedproportions$X2, level= c("-12/-9","-9/-6","-6/-3","-3/0")), y=meltedproportions$value, 
                              group=meltedproportions$X1)) + 
  labs(title="Number of Nodes per Time-periods",x="Time periods", y = "Proportions of nodes") +
  geom_line(aes(color=meltedproportions$X1)) + geom_point(aes(color=meltedproportions$X1))#+theme_classic()