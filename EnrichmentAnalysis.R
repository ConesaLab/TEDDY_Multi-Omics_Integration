###########################################################
###############   EnrichmentAnalysis.R     ################
###########################################################
# Authors: Ricardo Ramirez Flores, Leandro Balzano-Nogueira
# Genetics Institute, University of Florida (Gainesville)

# This script performs a time-specific directional gene set enrichment analysis 
# using PIANO with KEGG pathways and MSigDB's Hallmarks

# PIANO calculates the directional Enrichment Analysis, using the consensus of the results of 
# different enrichment methods by aggregating
# their FDR corrected p values through the median

# This script describes step by step what was done to obtain the heatmap of the paper

###########################################################
homedir<- "/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/" # Home directory where all your results are going to be contained
setwd(homedir)
getwd()
###########################################################
# Functions:
source ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/ScriptsForNComm/Tools/loadTEDDYtools.R") # Location of TEDDYtools
####
plot_PIANO_consensus = function(GlobalResultsFolder, GSA_ALL, med_Pval = 0.3){
  # Input
  # -GSA_ALL object, a median pval threshold and a folder where to generate the plots
  #
  # Output
  # -Consensus P value complete matrices
  # -Folders with heatmaps
  
  GSA_ALL_consensus = GSA_ALL
  
  for(tag in names(GSA_ALL)){
    print(tag)
    
    tagFolder = paste(GlobalResultsFolder,tag,"/",sep="")
    
    bash_line = paste("mkdir",tagFolder,sep=" ")
    system(bash_line)
    
    TAG_master = GSA_ALL[[tag]]
    
    for(sGLS in names(TAG_master)){
      print(sGLS)
      file_name = paste(tagFolder,sGLS,".pdf",sep="")
      
      gsaList = TAG_master[[sGLS]]
      pdf(file = file_name,width = 10,height = 10)
      
      GSA_ALL_consensus[[tag]][[sGLS]] = consensusHeatmap(gsaList,cutoff=30,method="mean",adjusted=T,ncharLabel=50,cellnote = "medianPvalue",colorkey = F)
      
      pMat = na.omit(GSA_ALL_consensus[[tag]][[sGLS]][["pMat"]])
      
      pMat = pMat[rowSums(pMat<=0.3)>=1,]
      
      heatmap.2(pMat,Rowv = T, Colv = F, distfun = dist, dendrogram = "row",
                trace = "none",density.info = "none",
                cellnote = round(pMat,2),notecol = "black",
                margins=c(7,12),notecex = .7,cexRow = .6,cexCol = .6)
      
      dev.off()
    }
  }
  
  return(GSA_ALL_consensus)
  
}

####

#
# In order to generate the whole dynamic of the gene sets we need to be sure that we have results for all of them in all tests.
# this function allows you to keep all the sets that appear in every test for a given subgroup
#

getConsistentSets = function(GSA_ALL_consensus, GSET_reference){
  for(AG in names(GSA_ALL_consensus)){
    
    AG_genes = GSET_reference
    
    for(sGLS in names(GSA_ALL_consensus[[AG]])){
      
      pMat = na.omit(GSA_ALL_consensus[[AG]][[sGLS]][["pMat"]])
      
      AG_genes = intersect(AG_genes, rownames(pMat))
      
    }
    
    GSA_ALL_consensus[[AG]][["features"]] = AG_genes
  }
  
  return(GSA_ALL_consensus)
}

####

plotGSET_dynamics = function(GSA_ALL_consensus, med_Pval = 0.3, output_pdf){
  
  pdf(output_pdf,height = 11,width = 9)
  
  for(AG in names(GSA_ALL_consensus)[names(GSA_ALL_consensus)!="FirstAAb_4"]){
    print(AG)
    AG_consensus = GSA_ALL_consensus[[AG]]
    TimeLabels = c("time_0","time_3","time_6","time_9","time_12")
    Dir_Mat = c()
    Pval_Mat = c()
    
    for(TL in TimeLabels){
      print(TL)
      Tix = grep(TL,names(AG_consensus))
      if(identical(Tix,integer(0))){
        next
      }else{
        pMat = na.omit(AG_consensus[[names(AG_consensus)[Tix]]][["pMat"]])
      }
      
      # For each gset in a pMAT get the direction of change
      # pMat_Directionality =data.frame(apply(pMat,1,function(x){
      #   return(names(which(x<0.2))[1])
      #   
      # }),stringsAsFactors = F)
      
      pMat_Directionality =t(data.frame(apply(pMat,1,function(x){
        return(cbind(x[which.min(x)][1],names(x[which.min(x)][1])))
      }),stringsAsFactors = F,check.names = F))
      
      pMat_Directionality = (as.matrix(pMat_Directionality))[AG_consensus$features,,drop=F]
      
      PvalDF = pMat_Directionality[,1,drop=F]
      DirDF = pMat_Directionality[,2,drop=F]
      
      colnames(PvalDF) = colnames(DirDF) = TL
      
      Dir_Mat = cbind(DirDF,Dir_Mat)
      Pval_Mat = cbind(PvalDF,Pval_Mat)
      
    }
    class(Pval_Mat) <- "numeric"
    
    #First identify classes without a significant pval = 0.2 and cluster them...
    GOOD_SETS = !(rowSums(Pval_Mat <= med_Pval)==0)
    Pval_Mat = Pval_Mat[GOOD_SETS,]
    Dir_Mat = Dir_Mat[GOOD_SETS,]
    Dir_Mat = Dir_Mat[do.call(order, as.data.frame(Dir_Mat)),]
    Pval_Mat = Pval_Mat[rownames(Dir_Mat),]
    
    # MELT MATRICES
    M_Pval_Mat = melt(Pval_Mat)
    M_Dir_Mat = melt(Dir_Mat)
    
    #Create plotting DF
    pianoDF = cbind(M_Pval_Mat,M_Dir_Mat[,3])
    colnames(pianoDF) = c("GSET","TIME","PVAL","DIRECTION")
    #Modify PVAL so the higher the better
    pianoDF$PVALcomp = 1 - pianoDF$PVAL
    
    #MOdify levels, so the colors are consistent troughout the heatmaps
    pianoDF$DIRECTION = as.character(pianoDF$DIRECTION)
    #### This is for Version 1 
    #direction_classes = c("Mixed-directional (up)","Distinct-directional (up)","Distinct-directional (dn)","Mixed-directional (dn)","Non-directional")
    #directioncolors = c("steelblue", "seagreen3", "firebrick1","darkorchid1","bisque3")
    #directioncolors = c("brown1", "red", "blue1","deepskyblue1","white")
    #This is for Version2
    direction_classes = c("Distinct-directional (up)","Mixed-directional (up)","Non-directional","Mixed-directional (dn)","Distinct-directional (dn)")
    directioncolors = c("red","brown1", "white","deepskyblue1","blue1")
    LEGEND_IX = which(direction_classes%in%unique(pianoDF$DIRECTION))
    
    direction_classes = direction_classes[LEGEND_IX]
    pianoDF$DIRECTION = factor(pianoDF$DIRECTION, levels = direction_classes)
    
    
    
    #Assign the same colors every time
    directioncolors = directioncolors[LEGEND_IX]
    colorends = 1:(length(direction_classes)*2)
    colorends[(1:(length(direction_classes)) * 2) - 1] = "white"
    colorends[(1:(length(direction_classes)) * 2) ] = directioncolors
    
    # Reescaling values so we can plot different categories
    
    Nclasses = sort((unique(100 * (as.numeric(pianoDF$DIRECTION)-1 ))))
    
    pianoDF$PVALresc = pianoDF$PVALcomp + (100 * (as.numeric(pianoDF$DIRECTION)-1 ))
    
    scalerange <- range(pianoDF$PVALcomp)
    gradientends <- scalerange + rep(Nclasses, each=2)
    
    
    pianoDF2 <- pianoDF[order(pianoDF$TIME,pianoDF$PVALresc),]
    pianoDF2
    pianoDF2$GSET <- factor(pianoDF2$GSET, levels = rev(unique(as.character(pianoDF2$GSET))))
    
    
    p = ggplot(pianoDF, aes(TIME, GSET)) + 
      geom_tile(aes(fill = PVALresc), colour = "white") + 
      scale_fill_gradientn(colours = colorends, values = rescale(gradientends)) + 
      scale_x_discrete("", expand = c(0, 0)) + 
      scale_y_discrete("", expand = c(0, 0)) + 
      geom_text(aes(label = round(PVAL, 2))) +
      theme_grey(base_size = 9) + 
      theme(legend.position = "none",
            axis.ticks = element_blank(), 
            axis.text.x = element_text(angle = 330, hjust = 0))  + ggtitle(AG)
    
    print(p)
    
    pv3 = ggplot(pianoDF2, aes(TIME, GSET)) + 
      geom_tile(aes(fill = PVALresc), colour = "white") + 
      scale_fill_gradientn(colours = colorends, values = rescale(gradientends)) + 
      scale_x_discrete("", expand = c(0, 0)) + 
      scale_y_discrete("", expand = c(0, 0)) + 
      geom_text(aes(label = round(PVAL, 2))) +
      theme_grey(base_size = 9) + 
      theme(legend.position = "none",
            axis.ticks = element_blank(), 
            axis.text.x = element_text(angle = 330, hjust = 0))  + ggtitle(AG)
    
    print(pv3)
    
  }
  write.table(pianoDF, file="pianoDFvariante.txt")
  dev.off()
  
  
}

Whicharethepathways<- function (data=pianoDFsmallSubset, threshold=thrs, Do_you_want_to_know="yes") {
  for (thr in 1:length(threshold)) {
    #print(thr)}
    significants<-data[data$PVAL<=threshold[thr],];print(paste(length(unique(significants$GSET)), " significant pathways out of", length(unique(data$GSET)), "in at least 1 timepoint with a threshold of", threshold[thr]))
    
    if (Do_you_want_to_know == "yes") {
      print(unique(as.vector(significants$GSET)) )
    } else print ("None pathways showed because you don't wanna know who they are dude!")
  }
}

###########################################################
# Data:

# Gene Expression
GE_Processed<-read.csv ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/GeneExpression/GE_Processed.csv",header = TRUE)
GE_Processed[1:10,1:5]

# Misc:
# Response Variable
CohortData<-read.csv ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/CohortData.csv",header = TRUE)
CohortData[1:10,]

LM_globalTargets<-read.csv(file = "/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/LM_globalTargets.csv", row.names = 1)
colnames(LM_globalTargets)<-c("sample_mask_id", "mask_id", "outcome","time","FirstAAb")

# List of Cases with at least 3 out of 5 time points with data
patients3tps<-data.frame(V1=CohortData$Individual.Id[CohortData$Model.or.Validation=="Model"])
patients3tps

###########################################################
# Libraries:
library("ggplot2")

###########################################################
# First step: Calculating KEGG Gene Set Enrichment Analysis through Piano
# This part requires time and is computationally demanding.

###########################################################
####################     KEGG     #########################
###########################################################
# Get Expression Data
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
# Time Specific Differential Expression Analysis using Limma to calculate logFC, AveExpr, t, P-value, adj P-value and B
# This is to compare the results between limma and Piano

ALL_limma_results = getSignatures(cMAT = preDat,targetMAT = LM_globalTargets,Features = NULL,time = c(0,3,6,9,12),g.pval = 0.05,s.pval=0.05,model_type = "limma")
Age_signatures_limma_results = getSignatures(cMAT = LM_globalEXPMAT,targetMAT = LM_globalTargets,Features = AgeFeatures,time = c(0,3,6,9,12),g.pval = 0.05,s.pval=0.05,model_type = "limma")
FAAB_signatures_limma_results = getSignatures(cMAT = LM_globalEXPMAT,targetMAT = LM_globalTargets,Features = FAABFeatures,time = c(0,3,6,9,12),g.pval = 0.05,s.pval=0.05,model_type = "limma")

###########################################################
# Get specific results from each test (discard global results from Age and FAAB analysis)

ALL_signatures = ALL_limma_results$RESULTS
ALL_signatures = list("ALL" = ALL_signatures)
head(ALL_signatures$ALL$g_time_0)
dim(ALL_signatures$ALL$g_time_0)


Age_signatures = Age_signatures_limma_results$RESULTS
subgroups = names(Age_signatures)
subgroups = subgroups[grep("^g_",subgroups)*-1]
Age_signatures = Age_signatures[subgroups]

FAAB_signatures = FAAB_signatures_limma_results$RESULTS
subgroups = names(FAAB_signatures)
subgroups = subgroups[grep("^g_",subgroups)*-1]
FAAB_signatures = FAAB_signatures[subgroups]

###########################################################
# Classify signatures by agegroup, FAAB - - GENERATE LIST OF GENE LEVEL STATISTICS

GLS_compilation = c(Age_signatures,FAAB_signatures)
GLS_ALL = list()
GLS_labs = names(GLS_compilation)
tag_vec = strsplit2(GLS_labs,split = ".time_")[,1]

for(tag in unique(tag_vec)){
  sel_GLS = GLS_labs[tag_vec==tag]
  GLS_ALL[[tag]] = GLS_compilation[sel_GLS]
}

GLS_ALL = c(ALL_signatures,GLS_ALL)
summary(GLS_ALL)


###########################################################

#KEGG_REST_sets
GSC = c()
for(GS in names(KEGG_REST_sets)){
  GS_list = cbind(KEGG_REST_sets[[GS]],GS)
  GSC = rbind(GSC,GS_list)
}
GSC = data.frame(GSC,stringsAsFactors = F)
GSC = loadGSC(GSC) # This is the Piano script to load a gene set collection
###########################################################
# RUN MULTI GSEA (This is the piano script to generate GSA from 4 different manners, mean, median, sum and maxmean)

GSA_ALL = GLS_ALL

for(tag in names(GLS_ALL)){
  print(tag)
  TAG_master = GLS_ALL[[tag]]
  for(sGLS in names(TAG_master)){
    print(sGLS)
    GLS = TAG_master[[sGLS]]
    GSA_ALL[[tag]][[sGLS]] = runMultiGSA(GSC = GSC,GLS = GLS)
  }
}

KEGG_GSA_ALL = GSA_ALL

dim(GLS_ALL$ALL$g_time_0)

# This take a lot of time, please go to next check point to advance faster
#save(KEGG_GSA_ALL,file= "/media/data/leobalzano/ScriptsForTEDDY/Data/Piano/KEGG_GSA_ALL.ro")

###########################################################
#     KEGG Analysis:
# This is the Gene set enrichment analysis, without including a genelist
# KEGG_GSA_ALL is the results from Piano through runMultiGSA strategy calculating by 4 different manners, mean, median, sum and maxmean)
# KEGG_REST_sets Are the list of genes inside KEGG pathways
summary(KEGG_GSA_ALL)
summary(KEGG_GSA_ALL$ALL)
summary(KEGG_GSA_ALL$ALL$g_time_0$median)
KEGG_GSA_ALL$ALL$g_time_0$mean$nGenesUp


summary(KEGG_REST_sets)
str(KEGG_REST_sets)
summary(KEGG_REST_sets$`Glycolysis_/_Gluconeogenesis`)
KEGG_REST_sets$`Glycolysis_/_Gluconeogenesis`

###########################################################
#KEGG_REST_sets
GSC = c()
for(GS in names(KEGG_REST_sets)){
  GS_list = cbind(KEGG_REST_sets[[GS]],GS)
  GSC = rbind(GSC,GS_list)
}
GSC = data.frame(GSC,stringsAsFactors = F)

#KEGG_GSA_consensus = plot_PIANO_consensus(GlobalResultsFolder = "/media/data/leobalzano/ScriptsForTEDDY/Data/Piano/Results/",GSA_ALL = KEGG_GSA_ALL,med_Pval = 1)
# Generate Consensus analysis
#KEGG_GSA_consensus = getConsistentSets(GSA_ALL_consensus = KEGG_GSA_consensus, GSET_reference = unique(GSC[,2]))
#summary(KEGG_GSA_consensus)
#summary(KEGG_GSA_consensus$ALL)
#summary(KEGG_GSA_consensus$ALL$g_time_12)
#plotGSET_dynamics(GSA_ALL_consensus = KEGG_GSA_consensus,med_Pval = 1,output_pdf = "/media/data/leobalzano/ScriptsForTEDDY/Data/Piano/Results/TESTGSets_dynamics.pdf")


##########################################################
# Enriched processes in the gene list selected by NPLSDA #
##########################################################
# This is to calculate the results just including the list of genes versus Kegg
KEGGbig<- TEDDY_geneSets$KEGG

##########################################################
# We performed the hypergeometric test subsetting the list of genes to the ones inside TEDDY genes 
# Loop to know which are the genes in the genes list, inside the KEGG pathways
GenesinsidePaths<-list()
GenesinsidePathstita<-NULL

for (i in 1:length(KEGGbig)) {
  GenesinsidePathstita<-NPLSDA_genes[NPLSDA_genes %in% KEGGbig[[i]]]
  NAMEGenesinsidePathstita<- names(KEGGbig[i])
  GenesinsidePathstitaTEMP<-list(GenesinsidePathstita)
  GenesinsidePaths[NAMEGenesinsidePathstita]<-GenesinsidePathstitaTEMP
}
GenesinsidePaths
summary(GenesinsidePaths)

DF1<-data.frame(summary(GenesinsidePaths))
DF2<-DF1[DF1$Var2=="Length",]
DF2$Freq<-as.numeric(as.character((DF2$Freq)))
str(DF2)
DF3<-DF2[DF2$Freq > 0,]

rownames(DF3)<-DF3[,1];DF3<-DF3[,-1]
rownames(DF3)

GenesinsidePathsnoEmpties<-GenesinsidePaths[match (rownames(DF3), names(GenesinsidePaths))]
GenesinsidePathsnoEmpties #These are the paths with at least one gene of the selection 

length(GenesinsidePaths)
length(GenesinsidePathsnoEmpties)  # 160 pathways out of 220 has at least one gene 
sum(DF2$Freq) 
vec<-data.frame(unlist(GenesinsidePathsnoEmpties))
vecgenesNPLSDA<-vec$unlist.GenesinsidePathsnoEmpties.
length(unique(vecgenesNPLSDA)) # Just 184 genes out of 862 are in the 160 paths out of the 220 that KEGG has


##########################################################
# Here we retained all the genes present in all these 160 paths

KEGGsmall<- KEGGbig[match (rownames(DF3),names(KEGGbig))]
length(KEGGbig) # 220
length(KEGGsmall)  # 160
vecallgenesofPathstep1<-data.frame(unlist(KEGGsmall))
vecallgenesofPath<-vecallgenesofPathstep1$unlist.KEGGsmall.
length(unique(vecallgenesofPath)) # Just 5245 genes are present in the 160 pathways

###########################################################
######     Hypergeometric test of all subsetted     #######
###########################################################
# Now the analysis is performed with the reduced list of genes related to the pathways with at least one gene inside a pathway

###########################################################
# Subset the data to the mentioned size 5245 genes
vecallgenesofPath2<-as.vector(unique(vecallgenesofPath))  # 5245
length(vecallgenesofPath2)

SubsetEXPMAT<-LM_globalEXPMAT[rownames(LM_globalEXPMAT) %in% vecallgenesofPath2,]
dim(SubsetEXPMAT)  

###########################################################
# Now subset by the 136 individuals
smallSubsetSampleMaskIDs<-wholedescriptivevars[ wholedescriptivevars$Mask.Id %in% rownames(outcomedummyarray136),]
dim(wholedescriptivevars)
dim(smallSubsetSampleMaskIDs)
length(unique(smallSubsetSampleMaskIDs$sample_mask_id))

###########################################################
SmallerSubsetEXPMAT<-SubsetEXPMAT[,colnames(SubsetEXPMAT) %in% smallSubsetSampleMaskIDs$sample_mask_id]
dim(SubsetEXPMAT)
dim(SmallerSubsetEXPMAT)
###########################################################
dim(LM_globalTargets)


LM_globalTargetsSmallerSubset<-LM_globalTargets[LM_globalTargets$sample_mask_id %in% smallSubsetSampleMaskIDs$sample_mask_id,]
dim(LM_globalTargets)
dim(LM_globalTargetsSmallerSubset)

###################################################
###############     GSEA test     #################
###################################################

ALL_limma_results = getSignatures(cMAT = SmallerSubsetEXPMAT,targetMAT = LM_globalTargetsSmallerSubset,Features = NULL,time = c(0,3,6,9,12),g.pval = 0.05,s.pval=0.05,model_type = "limma")

ALL_signatures = ALL_limma_results$RESULTS
ALL_signatures = list("ALL" = ALL_signatures)
names(ALL_signatures)
head(ALL_signatures$ALL$g_time_0)



###############################################
#######     Restructuring the data     ########
###############################################

GLS_ALL = list()
GLS_ALL = c(ALL_signatures,GLS_ALL)
summary(GLS_ALL$ALL)
dim(GLS_ALL$ALL$g_time_0)

###############################################
########     GENERATE GENE SETS     ###########
###############################################
# The gene sets would be KEGGsmall
length(KEGGsmall)
KEGGsmall



GSC = c()
for(GS in names(KEGGsmall)){
  GS_list = cbind(KEGGsmall[[GS]],GS)
  GSC = rbind(GSC,GS_list)
}
GSC = data.frame(GSC,stringsAsFactors = F)
GSC_KEGGsmall = loadGSC(GSC)
str(GSC_KEGGsmall)

#save(GSC_KEGGsmall,file= "/media/data/leobalzano/ScriptsForTEDDY/Data/Piano/Results/GSC_KEGGsmall.RData")

###############################################
# RUN MULTI GSEA with Piano
GLS_ALL
GSA_ALL = ALL_signatures


for(tag in names(GLS_ALL)){
  print(tag)
  TAG_master = GLS_ALL[[tag]]
  for(sGLS in names(TAG_master)){
    print(sGLS)
    GLS = TAG_master[[sGLS]]
    GSA_ALL[[tag]][[sGLS]] = runMultiGSA(GSC = GSC_KEGGsmall,GLS = GLS)
  }
}

KEGGGSAALLsubset = GSA_ALL
# This again takes an important amount of time
#save(KEGGGSAALLsubset,file= "/Users/leobalzano/Desktop/TEDDY/TeddyToolsV2/Piano/KEGG_Results/SUBSET/KEGGGSAALLsubset.RData")

###########################################################
# Create Heatmaps of dynamics with the subsetted data
summary(KEGGGSAALLsubset)
summary(KEGGGSAALLsubset$ALL)
summary(KEGGGSAALLsubset$ALL$g_time_0$mean)
head(KEGGGSAALLsubset$ALL$g_time_0)


# Plotting separately per time point
KEGG_GSA_consensus = plot_PIANO_consensus(GlobalResultsFolder = "/media/data/leobalzano/ScriptsForTEDDY/Data/Piano/Results/KEGG_GSA_consensus",GSA_ALL = KEGGGSAALLsubset,med_Pval = 1)

###########################################################
# Generate Consensus analysis
KEGG_REST_sets<-GSC_KEGGsmall$gsc

GSC = c()
for(GS in names(KEGG_REST_sets)){
  GS_list = cbind(KEGG_REST_sets[[GS]],GS)
  GSC = rbind(GSC,GS_list)
}
GSC = data.frame(GSC,stringsAsFactors = F)

summary (GSC)
unique(GSC[,2])

KEGG_GSA_consensus2 = getConsistentSets(GSA_ALL_consensus = KEGG_GSA_consensus, GSET_reference = unique(GSC[,2]))

#summary(KEGG_GSA_consensus2$ALL$g_time_0)

plotGSET_dynamics(GSA_ALL_consensus = KEGG_GSA_consensus2,med_Pval = 0.8,
                  output_pdf = "/media/data/leobalzano/ScriptsForTEDDY/Data/Piano/Results/GSets_dynamics_p0.8.pdf")


###########################################################
####     Obtaining the plot of the desired pathways     ###
###########################################################
# This is performed to create the heatmap with the pathways that we considered are most important
# for the Type 1 Diabetes islet of autoimmunity analysis.

##################################################################
# Hypergeometrictest with less pathways     ##
##################################################################
##################################################################
# We decided to remove pathways that are not considered related to the disease previous to the Hypergeometric test.
# This means first we are going to remove the pathways
# Second we calculate the test
###############################################
########     GENERATE GENE SETS     ###########
###############################################
# The gene sets would be from KEGGsmall
# We have to create even smaller pathway list

length(KEGGsmall)
KEGGsmall2<-data.frame(summary(KEGGsmall))
KEGGsmall2
KEGGsmall3<-KEGGsmall2[KEGGsmall2$Var2=="Length",]
dim(KEGGsmall3)
vectordepathways<-as.vector (KEGGsmall3$Var1)
selection<-c(
  "RNA degradation", 
  #"Leishmaniasis" , 
  #"Toxoplasmosis" , 
  #"African trypanosomiasis", 
  ##"Olfactory transduction" , 
  "Type I diabetes mellitus", 
  "Endocytosis" , 
  "Phagosome" , 
  "Inositol phosphate metabolism", 
  "Cell cycle" , 
  #"Osteoclast differentiation" , 
  "Fructose and mannose metabolism" , 
  "Drug metabolism" , 
  #"Long-term potentiation" , 
  "Cell adhesion molecules (CAMs)" , 
  #"Oocyte meiosis" , 
  "Lysine degradation" , 
  "Other types of O-glycan biosynthesis" , 
  "Tight junction" , 
  "RNA transport" , 
  #"Pathways in cancer" , 
  "Steroid hormone biosynthesis" , 
  "Cytokine-cytokine receptor interaction" , 
  "Wnt signaling pathway" , 
  #"Neuroactive ligand-receptor interaction" , 
  #"Fc gamma R-mediated phagocytosis" , 
  #"Huntington's disease" , 
  #"Neurotrophin signaling pathway" , 
  "Natural killer cell mediated cytotoxicity" , 
  "Chemokine signaling pathway" , 
  #"Cardiac muscle contraction" , 
  "Starch and sucrose metabolism" , 
  "Insulin signaling pathway" , 
  "Jak-STAT signaling pathway" , 
  "Pentose phosphate pathway" , 
  "Adipocytokine signaling pathway" , 
  #"Gastric acid secretion" , 
  "B cell receptor signaling pathway" ,
  "ECM-receptor interaction" ,
  #"Glioma" ,
  "Glycerophospholipid metabolism" ,
  "Focal adhesion" ,
  #"Prostate cancer" ,
  "ErbB signaling pathway" ,
  #"Hepatitis C" ,
  "mRNA surveillance pathway" ,
  "Retinol metabolism" ,
  "Calcium signaling pathway" ,
  "Endometrial cancer" ,
  "Arachidonic acid metabolism" ,
  "Phosphatidylinositol signaling system" ,
  "Hematopoietic cell lineage" ,
  "p53 signaling pathway" ,
  #"Porphyrin and chlorophyll metabolism" ,
  "Protein processing in endoplasmic reticulum" ,
  #"Chagas disease (American trypanosomiasis)" ,
  #"Alzheimer's disease" ,
  "Glutathione metabolism" ,
  "Regulation of actin cytoskeleton" ,
  #"Vasopressin-regulated water reabsorption" ,
  "Antigen processing and presentation" ,
  "Fatty acid elongation" ,
  "Metabolism of xenobiotics by cytochrome P450" ,
  "Ribosome" ,
  "Leukocyte transendothelial migration" ,
  "Oxidative phosphorylation" ,
  #"Melanogenesis" ,
  #"Tryptophan metabolism" ,
  "Aldosterone-regulated sodium reabsorption" ,
  #"Malaria" ,
  #"Colorectal cancer" ,
  #"Parkinson's disease" ,
  "mTOR signaling pathway" ,
  "Ubiquitin mediated proteolysis" ,
  #"Taste transduction" ,
  "MAPK signaling pathway" ,
  "Ascorbate and aldarate metabolism" ,
  ##"SNARE interactions in vesicular transport" ,
  "Biosynthesis of unsaturated fatty acids" ,
  "Fc epsilon RI signaling pathway" ,
  "Purine metabolism" ,
  "Complement and coagulation cascades" ,
  #"Viral myocarditis" ,
  "Intestinal immune network for IgA production" ,
  #"Progesterone-mediated oocyte maturation" ,
  #"Staphylococcus aureus infection" ,
  "Fatty acid degradation" ,
  "Proteasome" ,
  "Butanoate metabolism" ,
  "Notch signaling pathway" ,
  #"Salivary secretion" ,
  "Pancreatic secretion" ,
  "Citrate cycle (TCA cycle)" ,
  #"Shigellosis" ,
  "Glyoxylate and dicarboxylate metabolism" ,
  "Toll-like receptor signaling pathway" ,
  "Gap junction" ,
  #"Acute myeloid leukemia" ,
  "Base excision repair" ,
  #"Melanoma" ,
  "PPAR signaling pathway" ,
  "Lysosome" ,
  "Systemic lupus erythematosus" ,
  "Type II diabetes mellitus" ,
  #"Bile secretion" ,
  "Glycine, serine and threonine metabolism" ,
  "Asthma" ,
  "Pyrimidine metabolism" ,
  "One carbon pool by folate" ,
  "Aminoacyl-tRNA biosynthesis" ,
  "Cysteine and methionine metabolism" ,
  "Galactose metabolism" ,
  "Glycolysis / Gluconeogenesis" ,
  #"Basal cell carcinoma" ,
  "Selenocompound metabolism" ,
  "Folate biosynthesis" ,
  "Spliceosome" ,
  "NOD-like receptor signaling pathway" ,
  "DNA replication" ,
  "Protein digestion and absorption" ,
  "Apoptosis" ,
  #"Arrhythmogenic right ventricular cardiomyopathy (ARVC)" ,
  "Arginine and proline metabolism" ,
  #"Terpenoid backbone biosynthesis" ,
  ##"Hedgehog signaling pathway" ,
  "Fat digestion and absorption" ,
  "Pathogenic Escherichia coli infection" ,
  #"Axon guidance" ,
  #"Caffeine metabolism" ,
  "TGF-beta signaling pathway" ,
  "Glycerolipid metabolism" ,
  #"Vibrio cholerae infection" ,
  "Ether lipid metabolism" ,
  "Tyrosine metabolism" ,
  "Epithelial cell signaling in Helicobacter pylori infection",
  "Glycosphingolipid biosynthesis" ,
  "Peroxisome" ,
  "N-Glycan biosynthesis" ,
  "Carbohydrate digestion and absorption" ,
  "Alanine, aspartate and glutamate metabolism" ,
  ##"Bacterial invasion of epithelial cells" ,
  "beta-Alanine metabolism" ,
  "Primary immunodeficiency" ,
  #"Collecting duct acid secretion" ,
  "RNA polymerase" ,
  "Pyruvate metabolism" ,
  #"Non-homologous end-joining" ,
  "Cyanoamino acid metabolism" ,
  #"Phototransduction" ,
  "Amino sugar and nucleotide sugar metabolism" ,
  "Histidine metabolism" ,
  "Taurine and hypotaurine metabolism" ,
  "Phenylalanine metabolism" ,
  "Pantothenate and CoA biosynthesis" ,
  "Proximal tubule bicarbonate reclamation" ,
  #"Primary bile acid biosynthesis" ,
  #"Butirosin and neomycin biosynthesis" ,
  "Vitamin B6 metabolism" ,
  #"Biotin metabolism" ,
  "D-Arginine and D-ornithine metabolism" 
)


#
KEGGsmaller<-KEGGsmall [match (selection,names(KEGGsmall))]
length(KEGGsmall)
length(KEGGsmaller)   # 110 pathways


#############################################################
KEGG_REST_sets<- KEGGsmaller


GSC = c()
for(GS in names(KEGG_REST_sets)){
  GS_list = cbind(KEGG_REST_sets[[GS]],GS)
  GSC = rbind(GSC,GS_list)
}
GSC = data.frame(GSC,stringsAsFactors = F)
GSC = loadGSC(GSC)
GSCKEGGsmaller<-GSC
#save(GSCKEGGsmaller, file= "/Users/leobalzano/Desktop/TEDDY/TeddyToolsV2/Piano/KEGG_Results/SUBSET/GSCKEGGsmaller.RData")


##################################################################
#     Taking the total N of the genes

vecallgenesofPathstep1<-data.frame(unlist(KEGGsmaller))
vecallgenesofPath<-vecallgenesofPathstep1$unlist.KEGGsmaller.
length(unique(vecallgenesofPath)) # This is our new total N. 
##################################################################
# Now the analysis is performed with the reduced list of genes related to the pathways with at least 
# One gene inside the 862
##################################################################
# Subset the data to the mentioned genes size
vecallgenesofPath2<-as.vector(unique(vecallgenesofPath))  

SubsetEXPMAT<-LM_globalEXPMAT[rownames(LM_globalEXPMAT) %in% vecallgenesofPath2,]
dim(SubsetEXPMAT)  
##################################################################
colnames(SubsetEXPMAT)
SmallerSubsetEXPMAT<-SubsetEXPMAT[,colnames(SubsetEXPMAT) %in% smallSubsetSampleMaskIDs$sample_mask_id]
dim(SubsetEXPMAT)
dim(SmallerSubsetEXPMAT)
#############
dim(LM_globalTargets)

LM_globalTargetsSmallerSubset<-LM_globalTargets[LM_globalTargets$sample_mask_id %in% smallSubsetSampleMaskIDs$sample_mask_id,]
dim(LM_globalTargets)
dim(LM_globalTargetsSmallerSubset)


######################################################
#################     GSEA test     ##################
######################################################

ALL_limma_results = getSignatures(cMAT = SmallerSubsetEXPMAT,targetMAT = LM_globalTargetsSmallerSubset,Features = NULL,time = c(0,3,6,9,12),g.pval = 0.05,s.pval=0.05,model_type = "limma")

ALL_signatures = ALL_limma_results$RESULTS
ALL_signatures = list("ALL" = ALL_signatures)

###############################################
#######     Restructuring the data     ########
###############################################

GLS_ALL = list()
GLS_ALL = c(ALL_signatures,GLS_ALL)
summary(GLS_ALL$ALL)
dim(GLS_ALL$ALL$g_time_12)


##################################################
#############     RUN MULTI GSEA     #############
##################################################

for(tag in names(GLS_ALL)){
  print(tag)
  TAG_master = GLS_ALL[[tag]]
  for(sGLS in names(TAG_master)){
    print(sGLS)
    GLS = TAG_master[[sGLS]]
    GSA_ALL[[tag]][[sGLS]] = runMultiGSA(GSC = GSC,GLS = GLS)
  }
}



KEGGGSAALLSmallersubsetbyANA = GSA_ALL


###########################################################
###########################################################
##################     Check point 3     ##################
###########################################################
###########################################################

###########################################################
####     Obtaining the plot of the desired pathways     ###
###########################################################
# Create Heatmaps of dynamics
summary(KEGGGSAALLSmallersubsetbyANA)
summary(KEGGGSAALLSmallersubsetbyANA$ALL)
summary(KEGGGSAALLSmallersubsetbyANA$ALL$g_time_0$mean)

# Plotting separately the results per time point
KEGG_GSA_consensus = plot_PIANO_consensus(GlobalResultsFolder = "/media/data/leobalzano/ScriptsForTEDDY/Data/Piano/Results/Smaller/",GSA_ALL = KEGGGSAALLSmallersubsetbyANA,med_Pval = 1)

# Generate Consensus analysis

KEGG_REST_sets<-GSCKEGGsmaller$gsc


GSC = c()
for(GS in names(KEGG_REST_sets)){
  GS_list = cbind(KEGG_REST_sets[[GS]],GS)
  GSC = rbind(GSC,GS_list)
}
GSC = data.frame(GSC,stringsAsFactors = F)

summary (GSC)
unique(GSC[,2])

KEGG_GSA_consensus2 = getConsistentSets(GSA_ALL_consensus = KEGG_GSA_consensus, GSET_reference = unique(GSC[,2]))

summary(KEGG_GSA_consensus2$ALL$g_time_0)

plotGSET_dynamics(GSA_ALL_consensus = KEGG_GSA_consensus2,med_Pval = 1,output_pdf = "/media/data/leobalzano/ScriptsForTEDDY/Data/Piano/Results/Smaller/GSets_dynamics_p1.pdf")


###########################################################
# This is performed to create the heatmap with the pathways that we considered are most important
# for the Type 1 Diabetes islet of autoimmunity analysis.

# Creating Piano Data Frame 

dim(KEGG_GSA_consensus2$ALL$g_time_0$pMat) # 85 Pathways
dim(KEGG_GSA_consensus2$ALL$g_time_3$pMat) # 94 Pathways
dim(KEGG_GSA_consensus2$ALL$g_time_6$ pMat) # 109 Pathways
dim(KEGG_GSA_consensus2$ALL$g_time_9$pMat) # 86 Pathways
dim(KEGG_GSA_consensus2$ALL$g_time_12$pMat) # 87 Pathways

length(unique(c(rownames(KEGG_GSA_consensus2$ALL$g_time_12$pMat),
                rownames(KEGG_GSA_consensus2$ALL$g_time_9$pMat),
                rownames(KEGG_GSA_consensus2$ALL$g_time_6$pMat),
                rownames(KEGG_GSA_consensus2$ALL$g_time_3$pMat),
                rownames(KEGG_GSA_consensus2$ALL$g_time_0$pMat))))


AG<-names(KEGG_GSA_consensus)



GSA_ALL_consensus = KEGG_GSA_consensus2
med_Pval = 0.8
AG_consensus = GSA_ALL_consensus[[AG]]
TimeLabels = c("time_0","time_3","time_6","time_9","time_12")
dim(GSA_ALL_consensus$ALL$g_time_0$pMat)[1]
allpaths<-data.frame(rownames(KEGG_GSA_consensus2$ALL$g_time_6$ pMat))
rownames(allpaths)<-allpaths$rownames.KEGG_GSA_consensus2.ALL.g_time_6.pMat.
n <- max(dim(GSA_ALL_consensus$ALL$g_time_0$pMat)[1],
         dim(GSA_ALL_consensus$ALL$g_time_3$pMat)[1],
         dim(GSA_ALL_consensus$ALL$g_time_6$pMat)[1],
         dim(GSA_ALL_consensus$ALL$g_time_9$pMat)[1],
         dim(GSA_ALL_consensus$ALL$g_time_12$pMat)[1])

#####
Dir_Mat = NULL
Pval_Mat = NULL

for(TL in TimeLabels){
  print(TL)
  Tix = grep(TL,names(AG_consensus))
  if(identical(Tix,integer(0))){
    next
  }else{
    pMat = na.omit(AG_consensus[[names(AG_consensus)[Tix]]][["pMat"]])
    pMat2<-merge (allpaths,pMat,by=0, all.x=TRUE) 
    rownames(pMat2)<-pMat2[,1]; pMat2<-pMat2[,-c(1,2)]
  }
  
  pMat_Directionality =t(data.frame(apply(pMat,1,function(x){
    return(cbind(x[which.min(x)][1],names(x[which.min(x)][1])))
  }),stringsAsFactors = F,check.names = F))
  pMat_Directionality2<-merge(allpaths,pMat_Directionality, by=0, all.x=TRUE)
  rownames(pMat_Directionality2)<-pMat_Directionality2[,1]; pMat_Directionality2<-pMat_Directionality2[,-c(1,2)]
  dim(pMat_Directionality2)
  
  PvalDF = as.matrix(pMat_Directionality2[,1,drop=F])
  DirDF = as.matrix(pMat_Directionality2[,2,drop=F])
  
  colnames(PvalDF) = colnames(DirDF) = TL
  
  Dir_Mat = cbind(DirDF,Dir_Mat)
  Pval_Mat = cbind(PvalDF,Pval_Mat)
  
}
class(Pval_Mat) <- "numeric"

#First identify classes without a significant pval = 0.2 and cluster them...
Pval_Mat
Pval_Mat[is.na(Pval_Mat)] <- 1

GOOD_SETS = !(rowSums(Pval_Mat <= med_Pval)==0)
Pval_Mat = Pval_Mat[GOOD_SETS,]
Dir_Mat = Dir_Mat[GOOD_SETS,]
Dir_Mat = Dir_Mat[do.call(order, as.data.frame(Dir_Mat)),]
Pval_Mat = Pval_Mat[rownames(Dir_Mat),]

# MELT MATRICES
M_Pval_Mat = melt(Pval_Mat)
Dir_Mat[is.na(Dir_Mat)]<-"Non-directional"
M_Dir_Mat = melt(Dir_Mat)


#Create plotting DF
pianoDF = cbind(M_Pval_Mat,M_Dir_Mat[,3])
colnames(pianoDF) = c("GSET","TIME","PVAL","DIRECTION")
#Modify PVAL so the higher the better
pianoDF$PVALcomp = 1 - pianoDF$PVAL

#MOdify levels, so the colors are consistent troughout the heatmaps
pianoDF$DIRECTION = as.character(pianoDF$DIRECTION)
unique(pianoDF$GSET)

pianoDFsmallerSubset<-pianoDF
pianoDFsmallerSubset

direction_classes = c("Distinct-directional (up)","Mixed-directional (up)","Non-directional","Mixed-directional (dn)","Distinct-directional (dn)")
directioncolors = c("red","brown1", "white","deepskyblue1","blue1")
LEGEND_IX = which(direction_classes%in%unique(pianoDFsmallerSubset$DIRECTION))

direction_classes = direction_classes[LEGEND_IX]
pianoDFsmallerSubset$DIRECTION = factor(pianoDFsmallerSubset$DIRECTION, levels = direction_classes)


directioncolors = directioncolors[LEGEND_IX]
colorends = 1:(length(direction_classes)*2)
colorends[(1:(length(direction_classes)) * 2) - 1] = "white"
colorends[(1:(length(direction_classes)) * 2) ] = directioncolors

# Reescaling values so we can plot different categories

Nclasses = sort((unique(100 * (as.numeric(pianoDFsmallerSubset$DIRECTION)-1 ))))

pianoDFsmallerSubset$PVALresc = pianoDFsmallerSubset$PVALcomp + (100 * (as.numeric(pianoDFsmallerSubset$DIRECTION)-1 ))

scalerange <- range(pianoDFsmallerSubset$PVALcomp)
gradientends <- scalerange + rep(Nclasses, each=2)

#save(pianoDFsmallerSubset, file="/Users/leobalzano/Desktop/TEDDY/TeddyToolsV2/Piano/KEGG_Results/SUBSET/SmallerSubset/pianoDFsmallerSubset.RData")
pianoDFsmallerSubset
AG<-names(KEGG_GSA_consensus)

#p = ggplot(pianoDFsmallerSubset, aes(TIME, GSET)) + 
#  geom_tile(aes(fill = PVALresc), colour = "white") + 
#  scale_fill_gradientn(colours = colorends, values = rescale(gradientends)) + 
#  scale_x_discrete("", expand = c(0, 0)) + 
#  scale_y_discrete("", expand = c(0, 0)) + 
#  geom_text(aes(label = round(PVAL, 2))) +
#  theme_grey(base_size = 9) + 
#  theme(legend.position = "none",
#        axis.ticks = element_blank(), 
#        axis.text.x = element_text(angle = 330, hjust = 0))  + ggtitle(AG)

#print(p)

summary(pianoDFsmallerSubset)
pianoDFsmallerSubset
vectitodepathwaysSmallerSubsetbyANA<-as.vector(unique(pianoDFsmallerSubset$GSET))
#write.table(vectitodepathwaysSmallerSubsetbyANA, file = "/Users/leobalzano/Desktop/TEDDY/TeddyToolsV2/Piano/KEGG_Results/SUBSET/SmallerSubset/vectitodepathwaysSmallerSubsetbyANA.txt",sep="\t" )
#vectitodepathwaysSmallerSubsetbyANA <- read.delim("~/Desktop/TEDDY/TeddyToolsV2/Piano/KEGG_Results/SUBSET/SmallerSubset/vectitodepathwaysSmallerSubsetbyANA.txt", sep="\t")
vectitodepathwaysSmallerSubsetbyANA


selection<-c(
  #"Ubiquitin mediated proteolysis", #sig:-3D
  #"RNA degradation", #sig:-6D
  #"Spliceosome", #sig:-9D,-6D
  #"Apoptosis", #sig:-3D,0D
  "B cell receptor signaling pathway", #sig:-6U,0D
  #"Pantothenate and CoA biosynthesis", #sig:-3D
  #"mRNA surveillance pathway", 
  "Citrate cycle (TCA cycle)", #sig: 0U
  #"PPAR signaling pathway",
  
  #####"Cell cycle", #sig:0U
  #"Antigen processing and presentation", #sig:-3D
  #"One carbon pool by folate", #sig:0U
  "Cytokine-cytokine receptor interaction", #sig:-6U,-3D,0D
  "Natural killer cell mediated cytotoxicity",#sig:-3D
  "Pancreatic secretion",#sig:-9U,-6U
  "Gap junction",#sig:-6U
  "mTOR signaling pathway",#sig:0D
  
  #####"MAPK signaling pathway",#sig:-3D,0D
  #"Fat digestion and absorption",
  #"Adipocytokine signaling pathway",
  #"p53 signaling pathway",
  
  #####"NOD-like receptor signaling pathway",#sig:-3D
  #####"Endocytosis",#sig:-6U,0D
  #####"Toll-like receptor signaling pathway",#sig:-6U,-3D,0D
  #####"Hematopoietic cell lineage",#sig:0D
  "Jak-STAT signaling pathway",#sig:-3D,0D
  #"Folate biosynthesis",
  "Phagosome",#sig:-6U,0D
  #"Peroxisome",
  "Base excision repair", #sig:-9D,-3U,0U
  
  #####"DNA replication",#sig:-3U,0U
  "Proteasome",#sig:-3U,0U
  
  #####"Pyrimidine metabolism",#sig:-3U
  "Intestinal immune network for IgA production",#sig:0D
  #"Cysteine and methionine metabolism",
  #"Fatty acid degradation",
  #"Biosynthesis of unsaturated fatty acids",
  
  #####"RNA polymerase",#sig:-6U,-3U
  #####"Ribosome", #sig:-12U,-9U,-6D,-3U,0U
  #"Regulation of actin cytoskeleton",
  "Glycolysis / Gluconeogenesis",#sig:-3U,0U
  "Oxidative phosphorylation",#sig:-12U,-9U,-6U,-3U,0U
  
  #####"Proximal tubule bicarbonate reclamation",#sig:0U
  #"Primary immunodeficiency",
  #"Arginine and proline metabolism",
  #"Glycerolipid metabolism",
  "Glutathione metabolism", #sig:0U
  #"Steroid hormone biosynthesis",
  #"Tight junction",
  #"Glyoxylate and dicarboxylate metabolism",
  #"Carbohydrate digestion and absorption",
  
  #####"Taurine and hypotaurine metabolism",#sig:-6U,0D
  "Epithelial cell signaling in Helicobacter pylori infection",#sig:-6U
  
  #####"Fructose and mannose metabolism", #sig:-6U,-3U
  #"Tyrosine metabolism",
  #"Metabolism of xenobiotics by cytochrome P450",
  #"Butanoate metabolism",
  #"Protein processing in endoplasmic reticulum",
  "Type I diabetes mellitus",#sig:-3D,0D
  #"Glycosphingolipid biosynthesis",
  #"Lysine degradation",
  "Pentose phosphate pathway",#sig:-6U,-3U
  #"Pathogenic Escherichia coli infection",
  "Asthma",#sig:0D
  #"beta-Alanine metabolism",
  
  #####"ECM-receptor interaction",#sig:0D
  #####"Amino sugar and nucleotide sugar metabolism",#sig:0U
  #####"Aminoacyl-tRNA biosynthesis",#sig:-9D,-6D,-3U
  #"Lysosome",
  #"Type II diabetes mellitus"#,
  
  #####"Phenylalanine metabolism", #sig:0D
  #"Histidine metabolism",
  #"Cell adhesion molecules (CAMs)",
  #"Selenocompound metabolism",
  #"Glycine, serine and threonine metabolism",
  #"Starch and sucrose metabolism",
  #"Fc epsilon RI signaling pathway",#sig:-6U,0D
  "Insulin signaling pathway", #sig:-6U
  #"N-Glycan biosynthesis",
  #"Vitamin B6 metabolism",
  #"Other types of O-glycan biosynthesis",
  #"Calcium signaling pathway",#sig:0D
  
  #####"Drug metabolism",#sig:-3U
  #"Retinol metabolism",
  #"Ascorbate and aldarate metabolism",#sig:0D
  #"Purine metabolism",
  "Pyruvate metabolism",#sig:-3U
  
  #####"RNA transport",#sig:-12D
  #"Endometrial cancer",
  
  "ErbB signaling pathway",#sig:0D
  #####"Inositol phosphate metabolism",#sig:-6U,-3D
  #####"Phosphatidylinositol signaling system",#sig:-6U,-3D
  #"Complement and coagulation cascades",
  "Chemokine signaling pathway",#sig:-12U,-9U,-6D
  "Systemic lupus erythematosus",#sig:-12U
  
  #####"Aldosterone-regulated sodium reabsorption",#sig:-6U
  "Leukocyte transendothelial migration",#sig:-9U,-6U
  "TGF-beta signaling pathway",#sig:-9U,-3D
  #"Fatty acid elongation",
  
  #####"Focal adhesion",#sig:0D
  #"Glycerophospholipid metabolism",
  
  "Wnt signaling pathway",#sig:0D
  #####"Protein digestion and absorption",#sig:0U
  #"Galactose metabolism",
  #"Ether lipid metabolism",
  #"Arachidonic acid metabolism",
  #"Alanine, aspartate and glutamate metabolism",
  
  "Notch signaling pathway"#, #sig:-3D,0D
  #"Cyanoamino acid metabolism" 
)

#####################################################################################
# After discussion we determine that it is better to instead of having 5 color frames,
# One for each timepoint, it is better to have just one
AG<-names(KEGG_GSA_consensus2)
direction_classes = c("Distinct-directional (up)","Mixed-directional (up)","Non-directional","Mixed-directional (dn)","Distinct-directional (dn)")
directioncolors = c("red","red", "white","blue","blue")
LEGEND_IX = which(direction_classes%in%unique(pianoDFsmallerSubset$DIRECTION))

direction_classes = direction_classes[LEGEND_IX]
pianoDFsmallerSubset$DIRECTION = factor(pianoDFsmallerSubset$DIRECTION, levels = direction_classes)


directioncolors = directioncolors[LEGEND_IX]
colorends = 1:(length(direction_classes)*2)
colorends[(1:(length(direction_classes)) * 2) - 1] = "white"
colorends[(1:(length(direction_classes)) * 2) ] = directioncolors

# Reescaling values so we can plot different categories
Nclasses = sort((unique(100 * (as.numeric(pianoDFsmallerSubset$DIRECTION)-1 ))))

pianoDFsmallerSubset$PVALresc = pianoDFsmallerSubset$PVALcomp + (100 * (as.numeric(pianoDFsmallerSubset$DIRECTION)-1 ))

scalerange <- range(pianoDFsmallerSubset$PVALcomp)
gradientends <- scalerange + rep(Nclasses, each=2)

pianoDFALLSmallersubsetselection<-pianoDFsmallerSubset[pianoDFsmallerSubset$GSET %in% selection,]


############################################
# Selection sorted by 
direction_classes = c("Distinct-directional (up)","Mixed-directional (up)","Non-directional","Mixed-directional (dn)","Distinct-directional (dn)")
pianoDFALLSmallersubsetselection$directionnumbers[pianoDFALLSmallersubsetselection$DIRECTION=="Distinct-directional (up)"]<-1
pianoDFALLSmallersubsetselection$directionnumbers[pianoDFALLSmallersubsetselection$DIRECTION=="Mixed-directional (up)"]<-1
pianoDFALLSmallersubsetselection$directionnumbers[pianoDFALLSmallersubsetselection$DIRECTION=="Non-directional"]<-2
pianoDFALLSmallersubsetselection$directionnumbers[pianoDFALLSmallersubsetselection$DIRECTION=="Distinct-directional (dn)"]<-3
pianoDFALLSmallersubsetselection$directionnumbers[pianoDFALLSmallersubsetselection$DIRECTION=="Mixed-directional (dn)"]<-3

# This sorting i liked it 
pianoDFALLselectionSorted <- pianoDFALLSmallersubsetselection[order(pianoDFALLSmallersubsetselection$directionnumbers, pianoDFALLSmallersubsetselection$PVAL),]


#pianoDFALLselectionSorted
pianoDFALLselectionSorted$GSET <- factor(pianoDFALLselectionSorted$GSET, levels = rev(unique(as.character(pianoDFALLselectionSorted$GSET))))
pselectionSorted = ggplot(pianoDFALLselectionSorted, aes(TIME, GSET)) + 
  geom_tile(aes(fill = PVALresc), colour = "white") + 
  scale_fill_gradientn(colours = colorends, values = rescale(gradientends)) + 
  scale_x_discrete("", expand = c(0, 0)) + 
  scale_y_discrete("", expand = c(0, 0)) + 
  #geom_text(aes(label = round(PVAL, 2))) +
  geom_text(aes(label = ifelse(PVAL<0.2,"*","")), size=9) + # To put asterisk in significant
  theme_grey(base_size = 18) +
  #theme_grey(base_size = 5) +
  theme(legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text.x = element_text(angle = 330, hjust = 0))  + ggtitle(AG)

#pdf (file="pianoDFSMALLERSUBSETByANASorted171017.pdf",height = 11,width = 9)
print(pselectionSorted)
#dev.off()
# This is the heatmap reported in the paper

###########################################################