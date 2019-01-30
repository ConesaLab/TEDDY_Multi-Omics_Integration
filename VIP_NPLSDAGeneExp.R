###########################################################
###############     VIP_NPLSDAGeneExp.R     ###############
###########################################################
# Author: Leandro Balzano-Nogueira
# Genetics Institute, University of Florida (Gainesville)

# This script is to create the GeneExpression dataset to calculate the NPLSDA and VIP selection

###########################################################
homedir<- "/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/GeneExpression/" # Home directory where all your results are going to be contained
setwd(homedir)
getwd()
###########################################################
# Functions:
"/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/ScriptsForNComm/Tools/" # Location of TEDDYtools
source ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/ScriptsForNComm/Tools/TEDDYtools2.R") # These are the functions to reformat the data
source ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/ScriptsForNComm/Tools/NPLSDAfunctionsApr11.R") # These are the functions created to perform the NPLSDA

###########################################################
# Data:
# Gene Expression
GE_Raw<-read.csv ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/GeneExpression/GE_Raw.csv",header = TRUE)
GE_Raw[1:10,1:10]

GE_Processed<-read.csv ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/GeneExpression/GE_Processed.csv",header = FALSE)
GE_Processed[1:10,1:10]

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
Genexpression<-t(GE_Raw[-1,])
colnames(Genexpression)<- Genexpression[1,]; Genexpression<-Genexpression[-1,]

Genexpression <- transform(Genexpression, Individual.Id =as.numeric(as.character(Individual.Id)),
                           Age.in.Months = as.numeric(as.character(Age.in.Months)), 
                           Time.to.IA = as.numeric(as.character(Time.to.IA)))

Genexpression[1:10,1:10]

# Subsetting the data to the 136 pairs with values reported in 3 0ut 0f 5 time points
GeneExpression136<- Genexpression[is.element(Genexpression$Individual.Id, patients3tps[,1]),]
dim(GeneExpression136)
length(unique(GeneExpression136$Individual.Id))  # 136 PERFECT

###########################################################
# GE data in 3D structure
dim(GeneExpression136)   # 476 x 21288
colnames(GeneExpression136[1:23])
GENEX<-GeneExpression136[,c(1,4:21288)]       # Just ID and the variables


GENEX12<-GENEX[GeneExpression136$Time.to.IA == "-12",]; dim(GENEX12);dim(GENEX)
GENEX9<-GENEX[GeneExpression136$Time.to.IA == "-9",];dim(GENEX9)
GENEX6<-GENEX[GeneExpression136$Time.to.IA == "-6",];dim(GENEX6)
GENEX3<-GENEX[GeneExpression136$Time.to.IA == "-3",];dim(GENEX3)
GENEX0<-GENEX[GeneExpression136$Time.to.IA == "0",];dim(GENEX0)

dim(GENEX12);dim(GENEX9);dim(GENEX6);dim(GENEX3);dim(GENEX0)

###########################################################
# Merging with all cases
patients3tps2<-patients3tps
colnames(patients3tps2)<- "Individual.Id"
GENEX12total<-merge(patients3tps2,GENEX12, by="Individual.Id", all.x = T);dim(GENEX12total);dim(patients3tps2);dim(GENEX12)
GENEX12total[1:5,1:5]
rownames(GENEX12total)<- GENEX12total[,1]
GENEX12total<- as.matrix(GENEX12total[,c(-1)])
dim(GENEX12total)

GENEX9total<-merge(patients3tps2,GENEX9, by="Individual.Id", all.x = T);dim(GENEX9total);dim(patients3tps2);dim(GENEX9)
GENEX9total[1:5,1:5]
rownames(GENEX9total)<- GENEX9total[,1]
GENEX9total<- as.matrix(GENEX9total[,c(-1)])
dim(GENEX9total)

GENEX6total<-merge(patients3tps2,GENEX6, by="Individual.Id", all.x = T);dim(GENEX6total);dim(patients3tps2);dim(GENEX6total)
GENEX6total[1:5,1:5]
rownames(GENEX6total)<- GENEX6total[,1]
GENEX6total<- as.matrix(GENEX6total[,c(-1)])
dim(GENEX6total)

GENEX3total<-merge(patients3tps2,GENEX3, by="Individual.Id", all.x = T);dim(GENEX3total);dim(patients3tps2);dim(GENEX3)
GENEX3total[1:5,1:5]
rownames(GENEX3total)<- GENEX3total[,1]
GENEX3total<- as.matrix(GENEX3total[,c(-1)])
dim(GENEX3total)

GENEX0total<-merge(patients3tps2,GENEX0, by="Individual.Id", all.x = T);dim(GENEX0total);dim(patients3tps2);dim(GENEX0)
GENEX0total[1:5,1:5]
rownames(GENEX0total)<- GENEX0total[,1]
GENEX0total<- as.matrix(GENEX0total[,c(-1)])
dim(GENEX0total)
# Dimensions are 136*21285*5

###########################################################
arrayGENEXMarch <- array(data = NA, dim = c(136,21285,5),dimnames = list(NULL, NULL, c("-12","-9","-6", "-3", "0")))

arrayGENEXMarch
arrayGENEXMarch[,,1] <- GENEX12total
arrayGENEXMarch[,,2] <- GENEX9total
arrayGENEXMarch[,,3] <- GENEX6total
arrayGENEXMarch[,,4] <- GENEX3total
arrayGENEXMarch[,,5] <- GENEX0total

rownames(arrayGENEXMarch)<-rownames(GENEX12total)
colnames(arrayGENEXMarch)<-colnames (GENEX12total)
arrayGENEXMarch[1:50,1:5,1]
arrayGENEXMarch136<-arrayGENEXMarch
dim(arrayGENEXMarch136)   # 136 * 21285 * 5

# This is how the data array for Gene expression was created 

###########################################################
# Imputation (For convenience it must be done in an HPC)
# 1) Determining the best fitted model

modelGENEXnuevo<-bestfittedmodel (X=arrayGENEXMarch136,centering=0) # 0= No centering; 1= centering by Individuals; 2= centering by Variables;3= centering by Time
# The best model was 4,4,3 

# 2) Imputing the best fitted model data
FullarrayGEMARCH136<-Imputemethod(X=arrayGENEXMarch136,fac=c(4, 4, 3), conver = 1e-07, max.iter = 1000)

summary(FullarrayGEMARCH136)
dim(FullarrayGEMARCH136)


# 3) NPLSDA (For convenience it must be done in an HPC)
NPLSDAFullarrayGEMARCH136<-NPLSDAmod(XN=FullarrayGEMARCH136, YN=outcomedummyarray136, outcome.Y=NULL, factors=3, centering=0) 

# 4) Plotting

ploteoNPLSDAFullarrayGEMARCH136<- plotNPLSDAmod (X=NPLSDAFullarrayGEMARCH136, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                                 cutoff = 20, factors=2, penalty=1) 

###########################################################
#################    Variable Selection     ###############
###########################################################
# Variable selection by VIP3Dmodel2

summary(NPLSDAFullarrayGEMARCH136)
NPLSDAFullarrayGEMARCH136$VIP3Dmodel2


### From here
vipsoutcomemet<-data.frame(NPLSDAFullarrayGEMARCH136$VIP3Dmodel2)
#apply(vipsoutcomemet, 2, function(x) is.numeric(x))

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
###########################################################
# Retain just these variables in Gene Expression 
dim(FullarrayGEMARCH136)
PosLipselVars
#selectedgenesbydifVIPs<-list(VIP3Dmodel2.99p.1013g=PosLipselVars)
#summary(selectedgenesbydifVIPs)
#selectedgenesbydifVIPs<-c(selectedgenesbydifVIPs,tulo="pelotas")

FullarrayGenesVIPSelVars<-FullarrayGEMARCH136[,is.element(colnames(FullarrayGEMARCH136),
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
#pdf("VIP3Dmodel2genesmarch.99p.1027g.pdf")
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=2, penalty=2) 
#dev.off()

ploteoNPLSDAVIPselectedgenes
#This is the first selection

#selectedgenesbydifVIPs<-list(VIP3Dmodel2genesmarch.99p.1027g=PosLipselVars)
#summary(selectedgenesbydifVIPs)
###############selectedgenesbydifVIPs<-c(selectedgenesbydifVIPs,tulo="pelotas")

###########################################################
# Variable selection by VIP3Dmodel1
summary(NPLSDAFullarrayGEMARCH136)
NPLSDAFullarrayGEMARCH136$VIP3Dmodel1
dim(NPLSDAFullarrayGEMARCH136$VIP3Dmodel1)


### Desde aca
vipsoutcomemet<-data.frame(NPLSDAFullarrayGEMARCH136$VIP3Dmodel1)
#apply(vipsoutcomemet, 2, function(x) is.numeric(x))

thrs<-99
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
length(PosLipselVars)  # 296
###########################################################
# Retain just these variables in Gene Expression 2D
dim(FullarrayGEMARCH136)
PosLipselVars



FullarrayGenesVIPSelVars<-FullarrayGEMARCH136[,is.element(colnames(FullarrayGEMARCH136),
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
#pdf("VIP3Dmodel1.99p.1comp.213g.pdf")
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=2, penalty=2) 
#dev.off()


summary(selectedgenesbydifVIPs)
#selectedgenesbydifVIPs<-append(selectedgenesbydifVIPs,list(VIP3Dmodel1.99p.1comp.213g=PosLipselVars))

###########################################################
# Gene Expression USING VIP2D

summary(NPLSDAFullarrayGEMARCH136)
NPLSDAFullarrayGEMARCH136$VIP2D



### Desde aca
vipsoutcome2D<-data.frame(NPLSDAFullarrayGEMARCH136$VIP2D)

colp1<-paste(rownames(NPLSDAFullarrayGEMARCH136$VIP3Dmodel1),"-12", sep="_")
colp2<-paste(rownames(NPLSDAFullarrayGEMARCH136$VIP3Dmodel1),"-9", sep="_")
colp3<-paste(rownames(NPLSDAFullarrayGEMARCH136$VIP3Dmodel1),"-6", sep="_")
colp4<-paste(rownames(NPLSDAFullarrayGEMARCH136$VIP3Dmodel1),"-3", sep="_")
colp5<-paste(rownames(NPLSDAFullarrayGEMARCH136$VIP3Dmodel1),"0", sep="_")


rownames(vipsoutcome2D)<-c(colp1,colp2,colp3,colp4,colp5)

#apply(vipsoutcomemet, 2, function(x) is.numeric(x))
vipsoutcomemet<-vipsoutcome2D
thrs<-99
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
string<-strsplit (PosLipselVars,"_")
string<-listTOdata.frame(string)
genelist<- string[,1]
duplicated (genelist)
genes<-genelist[!duplicated(genelist)]
length(genes)

###########################################################
# Retain just these variables in Gene Expression 2D
dim(FullarrayGEMARCH136)
genes
#as.vector(genes)
FullarrayGenesVIPSelVars<-FullarrayGEMARCH136[,is.element(colnames(FullarrayGEMARCH136),
                                                          genes),]
dim(FullarrayGenesVIPSelVars)
#FullarrayGenesVIPSelVarsVIP2Dcomp2<-FullarrayGenesVIPSelVars
#save(FullarrayGenesVIPSelVarsVIP2Dcomp2,file="FullarrayGenesVIPSelVarsVIP2Dcomp2.RData")

### Now NPLSDA and Graph
NPLSDAarraygenesVIPselected<-NPLSDAmod(XN=FullarrayGenesVIPSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf ("VIP2D.99p.comp2.1030g.pdf")
#pdf ("VIP2D.99p.comp12.1179g.pdf")
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=2, penalty=2) 
#dev.off()
###########################################################
summary(selectedgenesbydifVIPs)
#selectedgenesbydifVIPs<-append(selectedgenesbydifVIPs,list(VIP2D.99p.comp2.1030g= as.vector(genes)))
#selectedgenesbydifVIPs<-append(selectedgenesbydifVIPs,list(VIP2D.99p.comp12.1179g= as.vector(genes)))
selectedgenesbydifVIPsMarch<-selectedgenesbydifVIPs
#save(selectedgenesbydifVIPsMarch,file="selectedgenesbydifVIPsMarch.RData")

###########################################################
# Venn Diagram of the four genelists #
summary(selectedgenesbydifVIPsMarch)

removeEMPTYstrings <- function(x) {
  
  newVectorWOstrings <- x[x != ""]
  return(newVectorWOstrings)
  
}
geneLS2 <- lapply(selectedgenesbydifVIPsMarch, removeEMPTYstrings)

tail(geneLS2[[3]])

summary(geneLS2)

lapply(geneLS2, tail) # Both methods return the same results

###########################################################
# plot Venn diagram
VENN.LIST <- geneLS2
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkmagenta", "darkblue","darkseagreen4","goldenrod2"), 
                          alpha=c(0.5,0.5,0.5,0.5), cex = 2, cat.fontface=4, 
                          category.names=c("VIP3D.2","VIP3D.1.comp1","VIP2D.comp2","VIP2D.comp12"), main="NPLSDAGeneLists",
                          cex.main=2
)

# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
grid.draw(venn.plot)
# As before, VIP2D. comp2 shares all its genes with VIP2D.comp12, so i can remove it from the analysis
###########################################################

summary(geneLS2)
geneLS3<-geneLS2[-3]
summary(geneLS3)

VENN.LIST2 <- geneLS3
venn.plot2 <- venn.diagram(VENN.LIST2 , NULL, fill=c("darkmagenta", "darkseagreen4","goldenrod2"), 
                           alpha=c(0.5,0.5,0.5), cex = 2, cat.fontface=4, 
                           category.names=c("VIP3D.2", "VIP3D.1.comp1","VIP2D.comp12"), main="NPLSDAGeneLists",
                           cex.main=2
)

grid.draw(venn.plot2)

############################################################
#     Selection of the genes present in the signatures     #
############################################################
# We are talking about the three signatures, because the other is completely overlapped so it would be the VENN.LIST2
summary(geneLS3)

VENN.LIST2 <- geneLS3
venn.plot2 <- venn.diagram(VENN.LIST2 , NULL, fill=c("darkmagenta", "darkseagreen4","goldenrod2"), 
                           alpha=c(0.5,0.5,0.5), cex = 2, cat.fontface=4, 
                           category.names=c("VIP3D.2", "VIP3D.1.comp12","VIP2D.comp2"), main="NPLSDAGeneLists",
                           cex.main=2
)

grid.draw(venn.plot2)

# To get the list of gene present in each Venn compartment we can use the gplots package
a <- venn(VENN.LIST2, show.plot=FALSE)
# You can inspect the contents of this object with the str() function
str(a)

# By inspecting the structure of the a object created, 
# you notice two attributes: 1) dimnames 2) intersections
# We can store the intersections in a new object named inters
inters <- attr(a,"intersections")

# We can summarize the contents of each venn compartment, as follows:
# in 1) ConditionA only, 2) ConditionB only, 3) ConditionA & ConditionB
lapply(inters, head) 

# 1) So for the union is the unique of the three balls
str(geneLS3)
UnionGeneList<-c(geneLS3$VIP3Dmodel2genesmarch.99p.1027g,geneLS3$VIP3Dmodel1.99p.1comp.213g,geneLS3$VIP2D.99p.comp12.1179g)
length(UnionGeneList)
length(unique(UnionGeneList))
####
# 2) Is the 129 genes of the intersection
IntersectionGenelist<-c(inters$`VIP3Dmodel2genesmarch.99p.1027g:VIP3Dmodel1.99p.1comp.213g:VIP2D.99p.comp12.1179g`)
length(IntersectionGenelist)
####
# 3) 2outof3 group, So it would be the one with 709,129,8,and 16 genes
twooutofthreeGeneList<-c(inters$`VIP3Dmodel1.99p.1comp.213g:VIP2D.99p.comp12.1179g`,inters$`VIP3Dmodel2genesmarch.99p.1027g:VIP2D.99p.comp12.1179g`,
                         inters$`VIP3Dmodel2genesmarch.99p.1027g:VIP3Dmodel1.99p.1comp.213g`,inters$`VIP3Dmodel2genesmarch.99p.1027g:VIP3Dmodel1.99p.1comp.213g:VIP2D.99p.comp12.1179g`
)
length(twooutofthreeGeneList)   # 862 genes
###########################################################
#     Testing DefinitiveGeneList with NPLSDA with 136 individuals     #
dim(FullarrayGEMARCH136)
dim(outcomedummyarray136)

dim(Outcomedummyarray)
# Union
arrayuniongenelist136<-FullarrayGEMARCH136[,is.element(colnames(FullarrayGEMARCH136),
                                                       UnionGeneList),]
dim(arrayuniongenelist136)


# Intersection
arrayintersectiongenelist136<-FullarrayGEMARCH136[,is.element(colnames(FullarrayGEMARCH136),
                                                              IntersectionGenelist),]
dim(arrayintersectiongenelist136)

#2 outofthree
arraytwooutofthreeGeneList136<-FullarrayGEMARCH136[,is.element(colnames(FullarrayGEMARCH136),
                                                               twooutofthreeGeneList),]
dim(arraytwooutofthreeGeneList136)
FullarrayGE862x136<-arraytwooutofthreeGeneList136
#save(FullarrayGE862x136, file = "FullarrayGE862x136.RData")

###########################################################
#   NPLSDA's
# Union:
NPLSDAarrayuniongenelist136<-NPLSDAmod(XN=arrayuniongenelist136,YN=outcomedummyarray136,factors = 2,COMP= c(2,2,2), conver = 1e-16, max.iteration = 10000,
                                       centering=0)

ploteoNPLSDAarrayuniongenelist136<- plotNPLSDAmod (X=NPLSDAarrayuniongenelist136, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                                   cutoff = 20, factors=2, penalty=2) 


# Intersection:
NPLSDAarrayintersectiongenelist136<-NPLSDAmod(XN=arrayintersectiongenelist136,YN=outcomedummyarray136,factors = 2,COMP= c(2,2,2), conver = 1e-16, max.iteration = 10000,
                                              centering=0)

ploteoNPLSDAarrayintersectiongenelist136<- plotNPLSDAmod (X=NPLSDAarrayintersectiongenelist136, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                                          cutoff = 20, factors=2, penalty=2) 

# Two out of three:
NPLSDAarraytwooutofthreeGeneList136<-NPLSDAmod(XN=arraytwooutofthreeGeneList136,YN=outcomedummyarray136,factors = 2,COMP= c(2,2,2), conver = 1e-16, max.iteration = 10000,
                                               centering=0)

ploteoNPLSDAarraytwooutofthreeGeneList136<- plotNPLSDAmod (X=NPLSDAarraytwooutofthreeGeneList136, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                                           cutoff = 20, factors=2, penalty=2) 

# The model with the higher explained variance resulted to be the one with 862 genes. is the one reported in the paper and in the milestones

# Repeat the same procedure for processed data
