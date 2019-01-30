###########################################################
#################   Permutationtest.R     #################
###########################################################
# Author: Leandro Balzano-Nogueira
# Genetics Institute, University of Florida (Gainesville)

# This script is to perform the permutation tests applied to determine if the features
# selected by NPLSDA gathers information concerning to the disease.

###########################################################
homedir<- "/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/" # Home directory where all your results are going to be contained
setwd(homedir)
getwd()
###########################################################
# Functions:
source ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/ScriptsForNComm/Tools/NPLSDAfunctionsApr11.R") # These are the functions created to perform the NPLSDA

###########################################################
# Data: 
# Gene Expression
load("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/GeneExpression/Ful1larrayGEMARCH306.RData")
load("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/GeneExpression/FullarrayGEMARCH136.RData")
load("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/GeneExpression/FullarrayGE862x136.RData")

# Metabolomics
load("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/Metabolomics/FullarrayMetabolomics136.RData")  # 136 x 1321 x 5
load("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/Metabolomics/FullarrayMetabolomics306x245.RData")  # 306 x 245 x 5
load("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/Metabolomics/Allmetabolomics136.RData")  # 136 x 245 x 5

# Response Variable
CohortData<-read.csv ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/CohortData.csv",header = TRUE)
CohortData[1:10,]

AllPatients<-data.frame(V1=CohortData$Individual.Id)

# List of Cases with at least 3 out of 5 time points with data
patients3tps<-data.frame(V1=CohortData$Individual.Id[CohortData$Model.or.Validation=="Model"])
patients3tps

###########################################################
# Libraries:
require("ropls")
library(reshape2)

###########################################################
###############     Gene Expression     ###################
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

# List of permutation tables and p values per timepoint 
# Summary table of R2s and Q2s
Full306<-c(0.5664,0.244)
Full136<-c(0.7644,0.472)
GE862x306<-c(0.5216,0.356)
GE862x136<-c(0.9553,0.761)

R2andQ2GE<-rbind(Full306,Full136,GE862x306,GE862x136)
colnames (R2andQ2GE)<-c("R2Y", "Q2")
###########################################################
# Random selection of 862 genes from the total for 136

permuts2<-list()
pvalsprom2<-list()
for (s in 1:2){
  muestra<-sample(colnames(FullarrayGEMARCH136), size=862, replace = FALSE, prob = NULL)
  muestraarray<-FullarrayGEMARCH136[,is.element(colnames(FullarrayGEMARCH136),muestra),]
  times<-c("-12","-9","-6","-3","0")
  pvals<-NULL
  listpvals<-list()
  
  perms<-NULL
  listperms<-list()
  
  for (i in 1:dim(muestraarray)[3]){
    TEDDYGEplsda <- opls(muestraarray[,,i], outcomedummyarray136,predI = 1, permI = 10,printL=FALSE,
                         plotL=FALSE) # For PLSDA
    pvalstita<-TEDDYGEplsda$summaryDF
    namepvals<-paste("time:",times[i],sep = "")
    pvals<-list(pvalstita)
    listpvals[namepvals]<-pvals
    
    permstita<-TEDDYGEplsda$suppLs$permMN[,c(2,3,7)]
    namesperms<-paste("time:",times[i],sep = "")
    perms<-list(permstita)
    listperms[namesperms]<-perms
  }
  
  permuts<-(listperms$`time:-12` + listperms$`time:-9`+ listperms$`time:-6` + listperms$`time:-3` + listperms$`time:0`)/5
  pvalsprom<-(listpvals$`time:-12` +listpvals$`time:-9`+listpvals$`time:-6`+listpvals$`time:-3`+listpvals$`time:0`)/5
  
  ##### permuts2
  permutssmall<-permuts
  namespermuts2<-paste("Iter",s,sep = "")
  per<-list(permutssmall)
  permuts2[namespermuts2]<-per
  
  ##### pvalsprom2
  pvalsSmall<-pvalsprom
  namespvalsprom2<-paste("Iter",s,sep = "")
  perpval<-list(pvalsSmall)
  pvalsprom2[namespermuts2]<-perpval
  
  
}

permuts3<-rbind(permuts2$`Iter1`[1:nrow(permuts2$`Iter1`),],permuts2$`Iter2`[1:nrow(permuts2$`Iter1`),],
                permuts2$`Iter3`[1:nrow(permuts2$`Iter1`),],permuts2$`Iter4`[1:nrow(permuts2$`Iter1`),],
                permuts2$`Iter5`[1:nrow(permuts2$`Iter1`),],permuts2$`Iter6`[1:nrow(permuts2$`Iter1`),],
                permuts2$`Iter7`[1:nrow(permuts2$`Iter1`),],permuts2$`Iter8`[1:nrow(permuts2$`Iter1`),],
                permuts2$`Iter9`[1:nrow(permuts2$`Iter1`),],permuts2$`Iter10`[1:nrow(permuts2$`Iter1`),])


pvalsprom3<-(pvalsprom2$`Iter1`+pvalsprom2$`Iter2`+pvalsprom2$`Iter3`+pvalsprom2$`Iter4`+pvalsprom2$`Iter5`+
               pvalsprom2$`Iter6`+pvalsprom2$`Iter7`+pvalsprom2$`Iter8`+pvalsprom2$`Iter9`+pvalsprom2$`Iter10`)/10


permutsx<-cbind(permuts3[,3],permuts3[,3])
permutsy<-permuts3[,c(1:2)]

###########################################################
# For 862x136 model
times<-c("-12","-9","-6","-3","0")
pvals<-NULL
listpvals<-list()

perms<-NULL
listperms<-list()
for (i in 1:dim(FullarrayGE862x136)[3]){
  #for (i in 1:2){
  TEDDYGEplsda <- opls(FullarrayGE862x136[,,i], outcomedummyarray136,predI = 1, permI = 100,printL=FALSE,
                       plotL=FALSE) # For PLSDA
  pvalstita<-TEDDYGEplsda$summaryDF
  namepvals<-paste("time:",times[i],sep = "")
  pvals<-list(pvalstita)
  listpvals[namepvals]<-pvals
  
  permstita<-TEDDYGEplsda$suppLs$permMN[,c(2,3,7)]
  namesperms<-paste("time:",times[i],sep = "")
  perms<-list(permstita)
  listperms[namesperms]<-perms
  
}
modelpermuts<-(listperms$`time:-12` + listperms$`time:-9`+ listperms$`time:-6` + listperms$`time:-3` + listperms$`time:0`)/5
modelpvalsprom<-(listpvals$`time:-12` +listpvals$`time:-9`+listpvals$`time:-6`+listpvals$`time:-3`+listpvals$`time:0`)/5

modelpermutsx<-cbind(modelpermuts[,3],modelpermuts[,3])
modelpermutsy<-modelpermuts[,c(1:2)]
modelpermutsy[1,]<-R2andQ2GE[4,]
modelpermutsy

###########################################################
# Calculating pvalues
# permuts3 for random
head(permuts3)
muR2random<-mean(permuts3[,1]) 
sdR2random<-sd(permuts3[,1])
muQ2random<-mean(permuts3[,2]) 
sdQ2random<-sd(permuts3[,2])
ntotalrandom<-dim(permuts3)[1]
randomNopermuted<-permuts3[permuts3[,3]==1,]
ntotalrandomNopermuted<-dim(randomNopermuted)[1]
zbarR2random<-mean(randomNopermuted[,1])
zbarQ2random<-mean(randomNopermuted[,2])

###########################################################
# modelpermuts for model

muR2model<-mean(modelpermuts[,1]) 
sdR2model<-sd(modelpermuts[,1])
muQ2model<-mean(modelpermuts[,2]) 
sdQ2model<-sd(modelpermuts[,2])
ntotalmodel<-dim(modelpermuts)[1]

zbarR2model<-modelpermuts[1,1]
zbarQ2model<-modelpermuts[1,2]
###########################################################
#### Calculations of p-values
# Verticals
tR2randomselection<- (zbarR2model - zbarR2random)/ (sdR2random/sqrt(ntotalrandomNopermuted)) # These are the vertical
tQ2randomselection<- (zbarQ2model - zbarQ2random)/ (sdQ2random/sqrt(ntotalrandomNopermuted))  # These are the vertical

pvalR2randomselection<-2* pt(-abs(tR2randomselection), df = ntotalrandomNopermuted-1)
pvalR2randomselection<-formatC(pvalR2randomselection, format = "e", digits = 2)
pvalQ2randomselection<-2* pt(-abs(tQ2randomselection), df = ntotalrandomNopermuted-1)
pvalQ2randomselection<-formatC(pvalQ2randomselection, format = "e", digits = 2)

###########################################################
# Horizontals
tR2modelpermuted<- (zbarR2model - muR2model)/ (sdR2model/sqrt(ntotalmodel)) # These are the horizontal
tQ2modelpermuted<- (zbarQ2model - muQ2model)/ (sdQ2model/sqrt(ntotalmodel))  # These are the horizontal
pvalR2modelpermuted<-2* pt(-abs(tR2modelpermuted), df = ntotalmodel-1)
pvalR2modelpermuted<-formatC(pvalR2modelpermuted, format = "e", digits = 2)
pvalQ2modelpermuted<-2* pt(-abs(tQ2modelpermuted), df = ntotalmodel-1)
pvalQ2modelpermuted<-formatC(pvalQ2modelpermuted, format = "e", digits = 2)

###############################
par(mfrow =c(1,1))
#pdf("permutationtest136x862.pdf")
plot(x=modelpermutsx,y=modelpermutsy, col="transparent", 
     
     main = paste("Permutation Plot 136x862x5\n pR2Random= ", pvalR2randomselection,", pQ2Random= ",pvalQ2randomselection,
                  " pR2Model= ", pvalR2modelpermuted,", pQ2Model= ",pvalQ2modelpermuted), 
     xlab = "Similarity(permuted y/y)" ,ylab ="R2/Q2" )
points(x=permutsx[,1],y=permutsy[,1], bg="navajowhite", col="black",pch=21,cex=1.5)
points(x=permutsx[,2],y=permutsy[,2], bg="tan1", col="black",pch=24,cex=1.5)
lineaR2<-lm(permutsy[,1]~permutsx[,1])
lineaQ2<-lm(permutsy[,2]~permutsx[,2])
# Slopes: 
YR2<-permutsy[1,1]
bR2<-lineaR2$coefficients[1]
mR2=YR2-bR2
abline(a=lineaR2$coefficients[1],b=mR2, col="navajowhite", lty=2)
#
YQ2<-permutsy[1,2]
bQ2<-lineaQ2$coefficients[1]
mQ2=YQ2-bQ2
abline(a=lineaQ2$coefficients[1],b=mQ2,col="tan1", lty=1)

###############################
points(x=modelpermutsx[,1],y=modelpermutsy[,1], bg="slategray2", col="black",pch=21,cex=1.5)
points(x=modelpermutsx[,2],y=modelpermutsy[,2], bg="steelblue3", col="black",pch=24,cex=1.5)
lineaR2<-lm(modelpermutsy[,1]~modelpermutsx[,1])
lineaQ2<-lm(modelpermutsy[,2]~modelpermutsx[,2])
# Slopes: 
modelYR2<-modelpermutsy[1,1]
bR2<-lineaR2$coefficients[1]
mR2=modelYR2-bR2
abline(a=lineaR2$coefficients[1],b=mR2, col="slategray2", lty=2)
#
YQ2<-modelpermutsy[1,2]
bQ2<-lineaQ2$coefficients[1]
mQ2=YQ2-bQ2
abline(a=lineaQ2$coefficients[1],b=mQ2,col="steelblue3", lty=1)
legend("bottomright",c("R2Random","Q2Random","R2Model","Q2Model"), col=c("navajowhite","tan1","slategray2","steelblue3"), 
       pch = c(16,17,16,17))
#dev.off()

###########################################################
##################     Metabolomics     ###################
###########################################################
# List of permutation tables and p values per timepoint 
# Summary table of R2s and Q2s
Full306<-c(0.4736,0.254)
Full136<-c(0.803,0.471)
Metab245x306<-c(0.598,0.264)
Metab245x136<-c(0.8151,0.305)

R2andQ2GE<-rbind(Full306,Full136,Metab245x306,Metab245x136)
colnames (R2andQ2GE)<-c("R2Y", "Q2")
###########################################################
# Random selection of 245 metabolites from the total for 136 individuals
# FullarrayMetabolomics136

permuts2<-list()
pvalsprom2<-list()
for (s in 1:10){
  muestra<-sample(colnames(FullarrayMetabolomics136), size=245, replace = FALSE, prob = NULL)
  muestraarray<-FullarrayMetabolomics136[,is.element(colnames(FullarrayMetabolomics136),muestra),]
  times<-c("-12","-9","-6","-3","0")
  pvals<-NULL
  listpvals<-list()
  
  perms<-NULL
  listperms<-list()
  
  for (i in 1:dim(muestraarray)[3]){
    #for (i in 1:2){
    TEDDYGEplsda <- opls(muestraarray[,,i], outcomedummyarray136,predI = 1, permI = 100,printL=FALSE,
                         plotL=FALSE) # For PLSDA
    pvalstita<-TEDDYGEplsda$summaryDF
    namepvals<-paste("time:",times[i],sep = "")
    pvals<-list(pvalstita)
    listpvals[namepvals]<-pvals
    
    permstita<-TEDDYGEplsda$suppLs$permMN[,c(2,3,7)]
    namesperms<-paste("time:",times[i],sep = "")
    perms<-list(permstita)
    listperms[namesperms]<-perms
  }
  
  permuts<-(listperms$`time:-12` + listperms$`time:-9`+ listperms$`time:-6` + listperms$`time:-3` + listperms$`time:0`)/5
  pvalsprom<-(listpvals$`time:-12` +listpvals$`time:-9`+listpvals$`time:-6`+listpvals$`time:-3`+listpvals$`time:0`)/5
  
  ##### permuts2
  permutssmall<-permuts
  namespermuts2<-paste("Iter",s,sep = "")
  per<-list(permutssmall)
  permuts2[namespermuts2]<-per
  
  ##### pvalsprom2
  pvalsSmall<-pvalsprom
  namespvalsprom2<-paste("Iter",s,sep = "")
  perpval<-list(pvalsSmall)
  pvalsprom2[namespermuts2]<-perpval
  
  
}

permuts3<-rbind(permuts2$`Iter1`[1:nrow(permuts2$`Iter1`),],permuts2$`Iter2`[1:nrow(permuts2$`Iter1`),],
                permuts2$`Iter3`[1:nrow(permuts2$`Iter1`),],permuts2$`Iter4`[1:nrow(permuts2$`Iter1`),],
                permuts2$`Iter5`[1:nrow(permuts2$`Iter1`),],permuts2$`Iter6`[1:nrow(permuts2$`Iter1`),],
                permuts2$`Iter7`[1:nrow(permuts2$`Iter1`),],permuts2$`Iter8`[1:nrow(permuts2$`Iter1`),],
                permuts2$`Iter9`[1:nrow(permuts2$`Iter1`),],permuts2$`Iter10`[1:nrow(permuts2$`Iter1`),])


pvalsprom3<-(pvalsprom2$`Iter1`+pvalsprom2$`Iter2`+pvalsprom2$`Iter3`+pvalsprom2$`Iter4`+pvalsprom2$`Iter5`+
               pvalsprom2$`Iter6`+pvalsprom2$`Iter7`+pvalsprom2$`Iter8`+pvalsprom2$`Iter9`+pvalsprom2$`Iter10`)/10

permutsx<-cbind(permuts3[,3],permuts3[,3])
permutsy<-permuts3[,c(1:2)]

###########################################################
# For 245x136 model
times<-c("-12","-9","-6","-3","0")
pvals<-NULL
listpvals<-list()

perms<-NULL
listperms<-list()
for (i in 1:dim(Allmetabolomics136)[3]){
  #for (i in 1:2){
  TEDDYGEplsda <- opls(Allmetabolomics136[,,i], outcomedummyarray136,predI = 1, permI = 100,printL=FALSE,
                       plotL=FALSE) # For PLSDA
  pvalstita<-TEDDYGEplsda$summaryDF
  namepvals<-paste("time:",times[i],sep = "")
  pvals<-list(pvalstita)
  listpvals[namepvals]<-pvals
  
  permstita<-TEDDYGEplsda$suppLs$permMN[,c(2,3,7)]
  namesperms<-paste("time:",times[i],sep = "")
  perms<-list(permstita)
  listperms[namesperms]<-perms
  
}
modelpermuts<-(listperms$`time:-12` + listperms$`time:-9`+ listperms$`time:-6` + listperms$`time:-3` + listperms$`time:0`)/5
modelpvalsprom<-(listpvals$`time:-12` +listpvals$`time:-9`+listpvals$`time:-6`+listpvals$`time:-3`+listpvals$`time:0`)/5

###########################################################
# Calculating pvalues

# permuts3 for random
head(permuts3)
muR2random<-mean(permuts3[,1]) 
sdR2random<-sd(permuts3[,1])
muQ2random<-mean(permuts3[,2]) 
sdQ2random<-sd(permuts3[,2])
ntotalrandom<-dim(permuts3)[1]
randomNopermuted<-permuts3[permuts3[,3]==1,]
ntotalrandomNopermuted<-dim(randomNopermuted)[1]
zbarR2random<-mean(randomNopermuted[,1])
zbarQ2random<-mean(randomNopermuted[,2])


# modelpermuts for model

muR2model<-mean(modelpermuts[,1]) 
sdR2model<-sd(modelpermuts[,1])
muQ2model<-mean(modelpermuts[,2]) 
sdQ2model<-sd(modelpermuts[,2])
ntotalmodel<-dim(modelpermuts)[1]
zbarR2model<-modelpermuts[1,1]
zbarQ2model<-modelpermuts[1,2]
###########################################################
#### Calculations of p-values
# Verticals
tR2randomselection<- (zbarR2model - zbarR2random)/ (sdR2random/sqrt(ntotalrandomNopermuted)) # These are the vertical
tQ2randomselection<- (zbarQ2model - zbarQ2random)/ (sdQ2random/sqrt(ntotalrandomNopermuted))  # These are the vertical

pvalR2randomselection<-2* pt(-abs(tR2randomselection), df = ntotalrandomNopermuted-1)
pvalR2randomselection<-formatC(pvalR2randomselection, format = "e", digits = 2)
pvalQ2randomselection<-2* pt(-abs(tQ2randomselection), df = ntotalrandomNopermuted-1)
pvalQ2randomselection<-formatC(pvalQ2randomselection, format = "e", digits = 2)

###########################################################
# Horizontals
tR2modelpermuted<- (zbarR2model - muR2model)/ (sdR2model/sqrt(ntotalmodel)) # These are the horizontal
tQ2modelpermuted<- (zbarQ2model - muQ2model)/ (sdQ2model/sqrt(ntotalmodel))  # These are the horizontal

pvalR2modelpermuted<-2* pt(-abs(tR2modelpermuted), df = ntotalmodel-1)
pvalR2modelpermuted<-formatC(pvalR2modelpermuted, format = "e", digits = 2)
pvalQ2modelpermuted<-2* pt(-abs(tQ2modelpermuted), df = ntotalmodel-1)
pvalQ2modelpermuted<-formatC(pvalQ2modelpermuted, format = "e", digits = 2)

modelpermutsx<-cbind(modelpermuts[,3],modelpermuts[,3])
modelpermutsy<-modelpermuts[,c(1:2)]
modelpermutsy[1,]<-R2andQ2GE[4,]
###############################

par(mfrow =c(1,1))
#pdf("permutationtestmetabolomics136x245.pdf")
plot(x=modelpermutsx,y=modelpermutsy, col="transparent", 
     #main = paste("Metabolomics 136x245x5"),
     main = paste("Permutation Plot Metabolomics 136x245x5\n pR2Random= ", pvalR2randomselection,", pQ2Random= ",pvalQ2randomselection,
                  " pR2Model= ", pvalR2modelpermuted,", pQ2Model= ",pvalQ2modelpermuted
     ), 
     xlab = "Similarity(permuted y/y)" ,ylab ="R2/Q2" )
points(x=permutsx[,1],y=permutsy[,1], bg="navajowhite", col="black",pch=21,cex=1.5)
points(x=permutsx[,2],y=permutsy[,2], bg="tan1", col="black",pch=24,cex=1.5)
lineaR2<-lm(permutsy[,1]~permutsx[,1])
lineaQ2<-lm(permutsy[,2]~permutsx[,2])
# Slopes: 
YR2<-permutsy[1,1]
bR2<-lineaR2$coefficients[1]
mR2=YR2-bR2
abline(a=lineaR2$coefficients[1],b=mR2, col="navajowhite", lty=2)
#
YQ2<-permutsy[1,2]
bQ2<-lineaQ2$coefficients[1]
mQ2=YQ2-bQ2
abline(a=lineaQ2$coefficients[1],b=mQ2,col="tan1", lty=1)
###############################
points(x=modelpermutsx[,1],y=modelpermutsy[,1], bg="slategray2", col="black",pch=21,cex=1.5)
points(x=modelpermutsx[,2],y=modelpermutsy[,2], bg="steelblue3", col="black",pch=24,cex=1.5)
lineaR2<-lm(modelpermutsy[,1]~modelpermutsx[,1])
lineaQ2<-lm(modelpermutsy[,2]~modelpermutsx[,2])
# Slopes: 
modelYR2<-modelpermutsy[1,1]
bR2<-lineaR2$coefficients[1]
mR2=modelYR2-bR2
abline(a=lineaR2$coefficients[1],b=mR2, col="slategray2", lty=2)
#
YQ2<-modelpermutsy[1,2]
bQ2<-lineaQ2$coefficients[1]
mQ2=YQ2-bQ2
abline(a=lineaQ2$coefficients[1],b=mQ2,col="steelblue3", lty=1)
legend("bottomright",c("R2Random","Q2Random","R2Model","Q2Model"), col=c("navajowhite","tan1","slategray2","steelblue3"), 
       pch = c(16,17,16,17))

#dev.off()