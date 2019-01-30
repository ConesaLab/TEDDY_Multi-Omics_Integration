###########################################################
#######   MetabolomicsHeatmapandSummaryTable.R     ########
###########################################################
# Author: Leandro Balzano-Nogueira
# Genetics Institute, University of Florida (Gainesville)

# This script is to create the heatmap of the metabolites present in metabolomics dataset

###########################################################
homedir<- "/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/" # Home directory where all your results are going to be contained
setwd(homedir)
getwd()
###########################################################
# Data:
# Metabolomics
AllmetabolitesConverted <- read.csv("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/AllmetabolitesConverted.txt", sep= "\t")
definitivemetabolitesSelectedConvertedmarch9 <- read.csv("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/definitivemetabolitesSelectedConvertedmarch9.txt", sep="\t")

load("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/FullarrayGEMETABenXFuncat.RData")

# List of Cases with at least 3 out of 5 time points with data
patients3tps<-data.frame(V1=CohortData$Individual.Id[CohortData$Model.or.Validation=="Model"])
patients3tps
###########################################################
# Libraries:
require("abind")
library(gplots)

###########################################################
Allmetabs<-AllmetabolitesConverted
the245<-definitivemetabolitesSelectedConvertedmarch9
names(the245)
length (unique(the245$Macromolecule)) # 15
sort(unique(the245$Macromolecule))
table(the245$Macromolecule)
vectorchico<-data.frame(table(the245$Macromolecule))
length(unique(Allmetabs$Macromolecule)) # 22
vectorall<-data.frame(table(Allmetabs$Macromolecule))
###########################################################
# Creating the Summary table

categorizedmetabolites<-merge(vectorall,vectorchico, by="Var1", all=TRUE)
categorizedmetabolites[is.na(categorizedmetabolites)]<-0
colnames(categorizedmetabolites)<-c("Category", "Total_Number", "Selection")
categorizedmetabolites

##################################################################
# Calculating The heatmap for selected metabolites per timepoint #
##################################################################
# This part is to calculate the heatmap per time point of the selected metabolites
#1)  We have to subset the original array to the selected features
#2) Change the colnames so they have the same name and by so create a mean per case and a mean per control
#3) With that, perform the heatmap
##################################################################
dim(FullarrayGEMETABenXFuncat)
colnames(FullarrayGEMETABenXFuncat)[863:1107]
metabs<-FullarrayGEMETABenXFuncat[,863:1107,]
dim(metabs)

colnames(metabs)<- definitivemetabolitesSelectedConvertedmarch9[ match( colnames( metabs ) ,
                                                                        definitivemetabolitesSelectedConvertedmarch9[ , "FunCat" ]),"Macromolecule" ]

###########################################################
# Step1: Subsetting all uniques

tabla<-table(colnames(metabs))
nombresunicos<-names(tabla[tabla==1])
dim(metabs)

testunicos<-metabs[,match(nombresunicos,colnames(metabs)),]
dim(testunicos)

###########################################################
# Step2: Grab all non-unique and calculate the mean
# declare the column names
nombresrepetidos<-names(tabla[tabla!=1])
dim(metabs)

metabsm12<-data.frame(metabs[,,1]);colnames(metabsm12)<-colnames(metabs)
metabsm9<-data.frame(metabs[,,2]);colnames(metabsm9)<-colnames(metabs)
metabsm6<-data.frame(metabs[,,3]);colnames(metabsm6)<-colnames(metabs)
metabsm3<-data.frame(metabs[,,4]);colnames(metabsm3)<-colnames(metabs)
metabsm0<-data.frame(metabs[,,5]);colnames(metabsm0)<-colnames(metabs)

metabsrepetidosm12<-sapply(nombresrepetidos, function(x) rowMeans(metabsm12 [, grep(x, colnames(metabsm12))] )  )
metabsrepetidosm9<-sapply(nombresrepetidos, function(x) rowMeans(metabsm9 [, grep(x, colnames(metabsm9))] )  )
metabsrepetidosm6<-sapply(nombresrepetidos, function(x) rowMeans(metabsm6 [, grep(x, colnames(metabsm6))] )  )
metabsrepetidosm3<-sapply(nombresrepetidos, function(x) rowMeans(metabsm3 [, grep(x, colnames(metabsm3))] )  )
metabsrepetidosm0<-sapply(nombresrepetidos, function(x) rowMeans(metabsm0 [, grep(x, colnames(metabsm0))] )  )

###########################################################
# Array for the repeated
dim(metabsrepetidosm12)[1]

arrayrepetidos<- array(data=NA, dim = c(dim(metabsrepetidosm12)[1],dim(metabsrepetidosm12)[2],5), dimnames=list(NULL,NULL, c("-12","-9","-6", "-3", "0")))

arrayrepetidos
arrayrepetidos[,,1] <- metabsrepetidosm12
arrayrepetidos[,,2] <- metabsrepetidosm9
arrayrepetidos[,,3] <- metabsrepetidosm6
arrayrepetidos[,,4] <- metabsrepetidosm3
arrayrepetidos[,,5] <- metabsrepetidosm0

rownames(arrayrepetidos)<-rownames(metabsrepetidosm12)
colnames(arrayrepetidos)<-colnames (metabsrepetidosm12)
arrayrepetidos[1:5,1,1]
dim(arrayrepetidos)   # 136 * 11 * 5

###########################################################
# Merging arrays
Metabolitesmeanforheatmap<- abind (testunicos,arrayrepetidos, along=2)
dim(Metabolitesmeanforheatmap)

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
outcomedummyarray136<-as.data.frame(CohortData2[,-c(1,2)])
rownames(outcomedummyarray136)<-rownames(CohortData2)
outcomedummyarray136

Metabolitesmeanforheatmapwithoutcome<- abind (outcomedummyarray136,Metabolitesmeanforheatmap, along=2)
dim(Metabolitesmeanforheatmapwithoutcome)
Metabolitesmeanforheatmapwithoutcome[1:20,,1]

###########################################################
# Decomposing the array to a data.frame
# This must be done just for cases
Metabolitesmeanforheatmapwithoutcomem12<-data.frame(Metabolitesmeanforheatmapwithoutcome[,,1])
Metabolitesmeanforheatmapwithoutcomem12

# Time -12
Casesm12<- Metabolitesmeanforheatmapwithoutcomem12[Metabolitesmeanforheatmapwithoutcomem12$Outcome==1,]

MeanCasesm12<-data.frame(colMeans(Casesm12)); 
MeanCasesm12<-data.frame(rownames(MeanCasesm12),MeanCasesm12)
MeanCasesm12<-data.frame(MeanCasesm12[-1,])
colnames(MeanCasesm12)<-c("Metabolite","time_12")
rownames(MeanCasesm12)<-MeanCasesm12[,1]
MeanCasesm12

# Time -9
Metabolitesmeanforheatmapwithoutcomem9<-data.frame(Metabolitesmeanforheatmapwithoutcome[,,2])
Casesm9<- Metabolitesmeanforheatmapwithoutcomem9[Metabolitesmeanforheatmapwithoutcomem9$Outcome==1,]

MeanCasesm9<-data.frame(colMeans(Casesm9))
MeanCasesm9<-data.frame(rownames(MeanCasesm9),MeanCasesm9)
MeanCasesm9<-MeanCasesm9[-1,]
colnames(MeanCasesm9)<-c("Metabolite","time_9")
rownames(MeanCasesm9)<-MeanCasesm9[,1]
MeanCasesm9

# Time -6
Metabolitesmeanforheatmapwithoutcomem6<-data.frame(Metabolitesmeanforheatmapwithoutcome[,,3])
Casesm6<- Metabolitesmeanforheatmapwithoutcomem6[Metabolitesmeanforheatmapwithoutcomem6$Outcome==1,]

MeanCasesm6<-data.frame(colMeans(Casesm6))
MeanCasesm6<-data.frame(rownames(MeanCasesm6),MeanCasesm6)
MeanCasesm6<-MeanCasesm6[-1,]
colnames(MeanCasesm6)<-c("Metabolite","time_6")
rownames(MeanCasesm6)<-MeanCasesm6[,1]
MeanCasesm6

# Time -3
Metabolitesmeanforheatmapwithoutcomem3<-data.frame(Metabolitesmeanforheatmapwithoutcome[,,4])
Casesm3<- Metabolitesmeanforheatmapwithoutcomem3[Metabolitesmeanforheatmapwithoutcomem3$Outcome==1,]

MeanCasesm3<-data.frame(colMeans(Casesm3))
MeanCasesm3<-data.frame(rownames(MeanCasesm3),MeanCasesm3)
MeanCasesm3<-MeanCasesm3[-1,]
colnames(MeanCasesm3)<-c("Metabolite","time_3")
rownames(MeanCasesm3)<-MeanCasesm3[,1]
MeanCasesm3

# Time 0
Metabolitesmeanforheatmapwithoutcomem0<-data.frame(Metabolitesmeanforheatmapwithoutcome[,,5])
Casesm0<- Metabolitesmeanforheatmapwithoutcomem0[Metabolitesmeanforheatmapwithoutcomem0$Outcome==1,]

MeanCasesm0<-data.frame(colMeans(Casesm0))
MeanCasesm0<-data.frame(rownames(MeanCasesm0),MeanCasesm0)
MeanCasesm0<-MeanCasesm0[-1,]
colnames(MeanCasesm0)<-c("Metabolite","time_0")
rownames(MeanCasesm0)<-MeanCasesm0[,1]
MeanCasesm0
###########################################################
# MERGE
Casestotalm12m9<-merge (MeanCasesm12,MeanCasesm9, by="Metabolite")
Casestotalm12m9
Casestotalm12m9m6<-merge (Casestotalm12m9,MeanCasesm6, by="Metabolite")
Casestotalm12m9m6
Casestotalm12m9m6m3<-merge (Casestotalm12m9m6,MeanCasesm3, by="Metabolite")
Casestotalm12m9m6m3

Casestotal<-merge (Casestotalm12m9m6m3,MeanCasesm0, by="Metabolite")
Casestotal
rownames(Casestotal)<-Casestotal[,1]
Casestotal<-Casestotal[,-1]
Casestotal<-as.matrix(Casestotal)
Casestotal

###########################################################
###### Heatmap
range(Casestotal)
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 149)

col_breaks = c(seq(-10000,-101,length=50),  # for red
               seq(-100,100,length=50),           # for white
               seq(101,3000,length=50))        # for blue 

#pdf("MetabolomicsHeatmap.pdf",height = 12,width = 10)
heatmap.2 (Casestotal,tracecol="transparent", margins=c(10,15),
           col=my_palette,
           key=F,
           #sepcolor="white",colsep=1:ncol(Casestotal),rowsep=1:nrow(Casestotal),
           breaks=col_breaks,
           #RowSideColors=ifelse(sorteddummy$sorteddummy==0,"blue","red"),
           Rowv = TRUE,
           Colv = FALSE
)
#dev.off()
