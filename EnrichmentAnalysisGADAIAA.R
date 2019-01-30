###########################################################
###########   EnrichmentAnalysisGADAIAA.R     #############
###########################################################
# Authors: Ricardo Ramirez Flores, Leandro Balzano-Nogueira
# Genetics Institute, University of Florida (Gainesville)

# This script performs a time-specific directional gene set enrichment analysis of IAA-first vs GADA-first individuals 
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
      
      GSA_ALL_consensus[[tag]][[sGLS]] = consensusHeatmap(gsaList,cutoff=30,method="median",adjusted=T,ncharLabel=50,cellnote = "medianPvalue",colorkey = F)
      
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
###########################################################
Whicharethepathways<- function (data, threshold=thrs, Do_you_want_to_know="yes") {
  for (thr in 1:length(threshold)) {
    #print(thr)}
    significants<-data[data$PVAL<=threshold[thr],];print(paste(length(unique(significants$GSET)), " significant pathways out of", length(unique(data$GSET)), "in at least 1 timepoint with a threshold of", threshold[thr]))
    
    if (Do_you_want_to_know == "yes") {
      print(unique(as.vector(significants$GSET)) )
    } else print ("None pathways showed because you don't wanna know who they are dude!")
  }
}

################################################################
Whicharethepathways2<- function (data, threshold=thrs, Do_you_want_to_know="yes") {
  for (thr in 1:length(threshold)) {
    #print(thr)
    significants<-data[data$PVAL<=threshold[thr],]
    Sigs<-print(paste(length(unique(significants$GSET)), " significant pathways out of", length(unique(data$GSET)), "in at least 1 timepoint with a threshold of", threshold[thr]))
    if (Do_you_want_to_know == "yes") {
      Paths<-unique(as.vector(significants$GSET))
      #print(unique(as.vector(significants$GSET)) )
    } else print ("None pathways showed because you don't wanna know who they are dude!")
  }
  result<-list(Sigs=Sigs,
               Pathways=Paths)
  return (result)
}
############################################################
# Plotting all significant based on a threshold


Unsupervisedselection<- function (data, threshold=0.01, Heatmaptitle="Insert title here") {
  Whicharethepathways (data=data, threshold = threshold, Do_you_want_to_know="yes")
  selection2<-data[data$PVAL<=threshold,]
  selection<-selection2$GSET
  #selection[[2]]
  pianoDFALLSmallersubsetselectionUnsupervised<-data[data$GSET %in% selection,]
  #length(unique(pianoDFALLSmallersubsetselectionUnsupervised$GSET)) # PERFECT
  ###########################
  # Disordered Heatmap
  pselection = ggplot(pianoDFALLSmallersubsetselectionUnsupervised, aes(TIME, GSET)) + 
    geom_tile(aes(fill = PVALresc), colour = "transparent") + 
    scale_fill_gradientn(colours = colorends, values = rescale(gradientends)) + 
    scale_x_discrete("", expand = c(0, 0)) + 
    scale_y_discrete("", expand = c(0, 0)) + 
    #geom_text(aes(label = round(PVAL, 2))) +   #Turned off for noPval figure
    geom_text(aes(label = ifelse(PVAL<threshold,"*","")), size=8) + # To put asterisk in significant
    theme_grey(base_size = 18) + 
    theme(legend.position = "none",
          axis.ticks = element_blank(), 
          axis.text.x = element_text(angle = 330, hjust = 0))  + ggtitle(Heatmaptitle)
  
  #pdf (file="pianoDFALLselectionLE0afterdiscussionNoPvalwithasterisks.pdf",height = 11,width = 9)
  print(pselection)
  #dev.off()
  ############################################
  # Selection sorted #
  ############################################
  
  direction_classes = c("Distinct-directional (up)","Mixed-directional (up)","Non-directional","Mixed-directional (dn)","Distinct-directional (dn)")
  pianoDFALLSmallersubsetselectionUnsupervised$directionnumbers[pianoDFALLSmallersubsetselectionUnsupervised$DIRECTION=="Distinct-directional (up)"]<-1
  pianoDFALLSmallersubsetselectionUnsupervised$directionnumbers[pianoDFALLSmallersubsetselectionUnsupervised$DIRECTION=="Mixed-directional (up)"]<-1
  pianoDFALLSmallersubsetselectionUnsupervised$directionnumbers[pianoDFALLSmallersubsetselectionUnsupervised$DIRECTION=="Non-directional"]<-2
  pianoDFALLSmallersubsetselectionUnsupervised$directionnumbers[pianoDFALLSmallersubsetselectionUnsupervised$DIRECTION=="Distinct-directional (dn)"]<-3
  pianoDFALLSmallersubsetselectionUnsupervised$directionnumbers[pianoDFALLSmallersubsetselectionUnsupervised$DIRECTION=="Mixed-directional (dn)"]<-3
  
  # This sorting i liked it 
  pianoDFALLselectionUnsupervisedSorted <- pianoDFALLSmallersubsetselectionUnsupervised[order(pianoDFALLSmallersubsetselectionUnsupervised$directionnumbers, pianoDFALLSmallersubsetselectionUnsupervised$PVAL),]
  
  ############################################
  #pianoDFALLselectionUnsupervisedSorted
  pianoDFALLselectionUnsupervisedSorted$GSET <- factor(pianoDFALLselectionUnsupervisedSorted$GSET, levels = rev(unique(as.character(pianoDFALLselectionUnsupervisedSorted$GSET))))
  ############################################
  pselectionUnsupervisedSorted = ggplot(pianoDFALLselectionUnsupervisedSorted, aes(TIME, GSET)) + 
    geom_tile(aes(fill = PVALresc), colour = "white") + 
    scale_fill_gradientn(colours = colorends, values = rescale(gradientends)) + 
    scale_x_discrete("", expand = c(0, 0)) + 
    scale_y_discrete("", expand = c(0, 0)) + 
    #geom_text(aes(label = round(PVAL, 2))) + #Turned off for noPval figure
    geom_text(aes(label = ifelse(PVAL<threshold,"*","")), size=8) + # To put asterisk in significant
    theme_grey(base_size = 18) +
    #theme_grey(base_size = 5) +
    theme(legend.position = "none",
          axis.ticks = element_blank(), 
          axis.text.x = element_text(angle = 330, hjust = 0))  + ggtitle(Heatmaptitle)
  
  #pdf (file="pianoDFSMALLERSUBSETSorted290917.pdf",height = 11,width = 9)
  print(pselectionUnsupervisedSorted)
  #dev.off()
}

###########################################################
# Data:
CohortData<-read.csv ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/CohortData2.csv",header = TRUE)
CohortData[1:10,]

LM_globalTargets<-read.csv ("/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/LM_globalTargets.csv",header = TRUE, row.names = 1)
LM_globalTargets

# List of Cases with at least 3 out of 5 time points with data
patients3tps<-data.frame(V1=CohortData$Individual.Id[CohortData$Model.or.Validation=="Model"])
patients3tps
###########################################################
# Libraries:
library("ggplot2")

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
# Get specific results from each test (discard global results from Age and FAAB analysis)

ALL_signatures = ALL_limma_results$RESULTS
ALL_signatures = list("ALL" = ALL_signatures)
summary(ALL_signatures)

Age_signatures = Age_signatures_limma_results$RESULTS
subgroups = names(Age_signatures)
subgroups = subgroups[grep("^g_",subgroups)*-1]
Age_signatures = Age_signatures[subgroups]
summary(Age_signatures)

FAAB_signatures = FAAB_signatures_limma_results$RESULTS
subgroups = names(FAAB_signatures)
subgroups = subgroups[grep("^g_",subgroups)*-1]
FAAB_signatures = FAAB_signatures[subgroups]
summary(FAAB_signatures$FirstAAb_1.time_0)
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
summary(GLS_ALL$FirstAAb_1)

###########################################################
#KEGG_REST_sets
GSC = c()
for(GS in names(KEGG_REST_sets)){
  GS_list = cbind(KEGG_REST_sets[[GS]],GS)
  GSC = rbind(GSC,GS_list)
}
GSC = data.frame(GSC,stringsAsFactors = F)
GSC = loadGSC(GSC)

KEGG_GSA_consensus = plot_PIANO_consensus(GlobalResultsFolder = "/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/",GSA_ALL = KEGG_GSA_ALL,med_Pval = 1)

# Generate Consensus analysis
KEGG_GSA_consensus = getConsistentSets(GSA_ALL_consensus = KEGG_GSA_consensus, GSET_reference = unique(GSC[,2]))

summary(KEGG_GSA_consensus)
summary(KEGG_GSA_consensus$ALL)
summary(KEGG_GSA_consensus$ALL$g_time_12)

summary(KEGG_GSA_consensus)
summary(KEGG_GSA_consensus$FirstAAb_1)
summary(KEGG_GSA_consensus$FirstAAb_1$FirstAAb_1.time_12)


plotGSET_dynamics(GSA_ALL_consensus = KEGG_GSA_consensus,med_Pval = 1,output_pdf = "/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/TESTGSets_dynamics.pdf")

##########################################################
# RUN MULTI GSEA

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

KEGG_GSA_ALLEFaab = GSA_ALL

# This take a lot of time, please go to next check point to advance faster
#save(KEGG_GSA_ALLEFaab,file= "/media/data/leobalzano/EnrichmentFaab/KEGG_GSA_ALLEFaab.ro")

###########################################################
#     KEGG Analysis:
# This is the Gene set enrichment analysis, without including a genelist

summary(KEGG_GSA_ALLEFaab)
summary(KEGG_GSA_ALLEFaab$ALL)
summary(KEGG_GSA_ALLEFaab$ALL$g_time_0$mean)
KEGG_GSA_ALLEFaab$ALL$g_time_0$mean$geneStatType

summary(KEGG_REST_sets)
str(KEGG_REST_sets)
summary(KEGG_REST_sets$`Glycolysis_/_Gluconeogenesis`)
KEGG_REST_sets$`Glycolysis_/_Gluconeogenesis`

###########################################################



##########################################################
# Enriched processes in the gene list selected by NPLSDA #
##########################################################
# Genes selected by VIP-NPLSDA
HyperResult = GSE_analysis(geneList = NPLSDA_genes,Annotation_DB = TEDDY_geneSets$KEGG)
str(TEDDY_geneSets$KEGG)
dim(HyperResult) # 220
head(HyperResult) # This is the result of the hypergeometric test


KEGGbig<- TEDDY_geneSets$KEGG

##########################################################
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
length(GenesinsidePathsnoEmpties)  # 160 pathways de 220 tienen al menos un gen de los 
sum(DF2$Freq) 
vec<-data.frame(unlist(GenesinsidePathsnoEmpties))
vecgenesNPLSDA<-vec$unlist.GenesinsidePathsnoEmpties.
length(unique(vecgenesNPLSDA)) # Just 184 genes out of 862 are in the 160 paths out of the 220 that KEGG has

##########################################################
# Grab all the genes present in all these 160 paths

KEGGsmall<- KEGGbig[match (rownames(DF3),names(KEGGbig))]
length(KEGGbig) # 220
length(KEGGsmall)  # 160
vecallgenesofPathstep1<-data.frame(unlist(KEGGsmall))
vecallgenesofPath<-vecallgenesofPathstep1$unlist.KEGGsmall.
length(unique(vecallgenesofPath)) # Just 5245 genes are present in the 160 pathways

###########################################################
######     Hypergeometric test of all subsetted     #######
###########################################################
# Now the analysis is performed with the reduced list of genes related to the pathways with at least 
# One gene inside the 862
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
colnames(SubsetEXPMAT)
SmallerSubsetEXPMAT<-SubsetEXPMAT[,colnames(SubsetEXPMAT) %in% smallSubsetSampleMaskIDs$sample_mask_id]
dim(SubsetEXPMAT)
dim(SmallerSubsetEXPMAT)
###########################################################
dim(LM_globalTargets)
LM_globalTargets[1:10,1:10]

LM_globalTargetsSmallerSubset<-LM_globalTargets[LM_globalTargets$sample_mask_id %in% smallSubsetSampleMaskIDs$sample_mask_id,]
dim(LM_globalTargets)
dim(LM_globalTargetsSmallerSubset)

###################################################
###############     GSEA test     #################
###################################################

ALL_limma_results = getSignatures(cMAT = SmallerSubsetEXPMAT,targetMAT = LM_globalTargetsSmallerSubset,Features = NULL,time = c(0,3,6,9,12),g.pval = 0.05,s.pval=0.05,model_type = "limma")

Age_signatures_limma_results = getSignatures(cMAT = LM_globalEXPMAT,targetMAT = LM_globalTargetsSmallerSubset,Features = AgeFeatures,time = c(0,3,6,9,12),g.pval = 0.05,s.pval=0.05,model_type = "limma")

FAAB_signatures_limma_results = getSignatures(cMAT = SmallerSubsetEXPMAT,targetMAT = LM_globalTargetsSmallerSubset,Features = FAABFeatures,time = c(0,3,6,9,12),g.pval = 0.05,s.pval=0.05,model_type = "limma")

# Get specific results from each test (discard global results from Age and FAAB analysis)

#ALL_signatures = ALL_limma_results$RESULTS
#ALL_signatures = list("ALL" = ALL_signatures)

#Age_signatures = Age_signatures_limma_results$RESULTS
#subgroups = names(Age_signatures)
#subgroups = subgroups[grep("^g_",subgroups)*-1]
#Age_signatures = Age_signatures[subgroups]

FAAB_signatures = FAAB_signatures_limma_results$RESULTS
subgroups = names(FAAB_signatures)
#subgroups = subgroups[grep("^g_",subgroups)*-1]
#FAAB_signatures = FAAB_signatures[subgroups]

# Classify signatures by agegroup, FAAB - - GENERATE LIST OF GENE LEVEL STATISTICS
###############################################
###     Restructuring the data for FAAB    ####
###############################################

#GLS_ALL = list()
#GLS_ALL = c(ALL_signatures,GLS_ALL)
#summary(GLS_ALL)

GLS_ALL = list()
GLS_ALL = c(FAAB_signatures,GLS_ALL)
summary(GLS_ALL)
dim(GLS_ALL$FirstAAb_1.time_0)

###############################################
########     GENERATE GENE SETS     ###########
###############################################
# The gene sets would be KEGGsmall
length(KEGGsmall)
KEGGsmall
KEGG_REST_sets<- KEGGsmall


GSC = c()
for(GS in names(KEGG_REST_sets)){
  GS_list = cbind(KEGG_REST_sets[[GS]],GS)
  GSC = rbind(GSC,GS_list)
}
GSC = data.frame(GSC,stringsAsFactors = F)
GSC_Faab = loadGSC(GSC)

#save(GSC_Faab,file= "/media/data/leobalzano/EnrichmentFaab/GSC_Faab.RData")
summary(GSC_Faab)

# RUN MULTI GSEA
GLS_ALL
#### Bypass for this to work
class (GLS_ALL)
GLS_ALL2<-GLS_ALL
names(GLS_ALL2)

MyGLS_ALL2<-list(FirstAAb_1.time_0=GLS_ALL2$FirstAAb_1.time_0,
                 FirstAAb_1.time_3=GLS_ALL2$FirstAAb_1.time_3,
                 FirstAAb_1.time_6=GLS_ALL2$FirstAAb_1.time_6,
                 FirstAAb_1.time_9=GLS_ALL2$FirstAAb_1.time_9,
                 FirstAAb_1.time_12=GLS_ALL2$FirstAAb_1.time_12,
                 FirstAAb_2.time_0=GLS_ALL2$FirstAAb_2.time_0,
                 FirstAAb_2.time_3=GLS_ALL2$FirstAAb_2.time_3,
                 FirstAAb_2.time_6=GLS_ALL2$FirstAAb_2.time_6,
                 FirstAAb_2.time_9=GLS_ALL2$FirstAAb_2.time_9,
                 FirstAAb_2.time_12=GLS_ALL2$FirstAAb_2.time_12
)
names(MyGLS_ALL2)

#MyGLS_ALL2<-list(FirstAAb_1.time_0=GLS_ALL2$FirstAAb_1.time_0)
####

GSA_ALL = MyGLS_ALL2

for(tag in names(MyGLS_ALL2)){
  print(tag)
  TAG_master = MyGLS_ALL2[tag]
  for(sGLS in names(TAG_master)){
    print(sGLS)
    GLS = TAG_master[[sGLS]]
    GSA_ALL[tag][[sGLS]] = runMultiGSA(GSC = GSC_Faab,GLS = GLS)
  }
}

KEGGGSAFAAbsubset = GSA_ALL
#save(KEGGGSAFAAbsubset,file= "/media/data/leobalzano/EnrichmentFaab/KEGGGSAFAAbsubset.RData")

###########################################################
# Create Heatmaps of dynamics
summary(KEGGGSAFAAbsubset)


FirstAAb_1<-list(FirstAAb_1.time_0=KEGGGSAFAAbsubset$FirstAAb_1.time_0,
                 FirstAAb_1.time_3=KEGGGSAFAAbsubset$FirstAAb_1.time_3,
                 FirstAAb_1.time_6=KEGGGSAFAAbsubset$FirstAAb_1.time_6,
                 FirstAAb_1.time_9=KEGGGSAFAAbsubset$FirstAAb_1.time_9,
                 FirstAAb_1.time_12=KEGGGSAFAAbsubset$FirstAAb_1.time_12)
summary(FirstAAb_1)

FirstAAb_2<-list(FirstAAb_2.time_0=KEGGGSAFAAbsubset$FirstAAb_2.time_0,
                 FirstAAb_2.time_3=KEGGGSAFAAbsubset$FirstAAb_2.time_3,
                 FirstAAb_2.time_6=KEGGGSAFAAbsubset$FirstAAb_2.time_6,
                 FirstAAb_2.time_9=KEGGGSAFAAbsubset$FirstAAb_2.time_9,
                 FirstAAb_2.time_12=KEGGGSAFAAbsubset$FirstAAb_2.time_12)
summary(FirstAAb_2)



MyKEGGGSAFAAbsubset<-list(FirstAAb_1=FirstAAb_1,FirstAAb_2=FirstAAb_2)

summary(MyKEGGGSAFAAbsubset)




# Step 1 
KEGG_GSA_consensus = plot_PIANO_consensus(GlobalResultsFolder = "/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/",GSA_ALL = MyKEGGGSAFAAbsubset,med_Pval = 1)

# Generate Consensus analysis
KEGG_REST_sets<-GSC_Faab$gsc


GSC = c()
for(GS in names(KEGG_REST_sets)){
  GS_list = cbind(KEGG_REST_sets[[GS]],GS)
  GSC = rbind(GSC,GS_list)
}
GSC = data.frame(GSC,stringsAsFactors = F)

summary (GSC)
unique(GSC[,2])

KEGG_GSA_consensus2 = getConsistentSets(GSA_ALL_consensus = KEGG_GSA_consensus, GSET_reference = unique(GSC[,2]))

summary(KEGG_GSA_consensus2$FirstAAb_1$FirstAAb_1.time_12)

plotGSET_dynamics(GSA_ALL_consensus = KEGG_GSA_consensus2,med_Pval = 0.1,output_pdf = "/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/GSets_dynamics_p01.pdf")

###########################################################
####     Obtaining the plot of the desired pathways     ###
###########################################################
# This is performed to create the heatmap with the pathways that we considered are most important
# for the Type 1 Diabetes islet of autoimmunity analysis.
#########################
#    GADAFirst     #
#########################
# Creating Piano Data Frame hacking script so we can play with it at will for GADA first
dim(KEGG_GSA_consensus2$FirstAAb_1$FirstAAb_1.time_12$pMat) # 154 Pathways
dim(KEGG_GSA_consensus2$FirstAAb_1$FirstAAb_1.time_9$pMat) # 88 Pathways
dim(KEGG_GSA_consensus2$FirstAAb_1$FirstAAb_1.time_6$pMat) # 156 Pathways
dim(KEGG_GSA_consensus2$FirstAAb_1$FirstAAb_1.time_3$pMat) # 156 Pathways
dim(KEGG_GSA_consensus2$FirstAAb_1$FirstAAb_1.time_0$pMat) # 156 Pathways




length(unique(c(rownames(KEGG_GSA_consensus2$FirstAAb_1$FirstAAb_1.time_12$pMat),
                rownames(KEGG_GSA_consensus2$FirstAAb_1$FirstAAb_1.time_9$pMat),
                rownames(KEGG_GSA_consensus2$FirstAAb_1$FirstAAb_1.time_6$pMat),
                rownames(KEGG_GSA_consensus2$FirstAAb_1$FirstAAb_1.time_3$pMat),
                rownames(KEGG_GSA_consensus2$FirstAAb_1$FirstAAb_1.time_0$pMat))))



AG="FirstAAb_1"
GSA_ALL_consensus = KEGG_GSA_consensus2
med_Pval = 0.8
AG_consensus = GSA_ALL_consensus[[AG]]
summary(AG_consensus)

TimeLabels = c("time_0","time_3","time_6","time_9","time_12")
#dim(GSA_ALL_consensus$ALL$g_time_0$pMat)[1]
allpaths<-data.frame(unique(c(rownames(KEGG_GSA_consensus2$FirstAAb_1$FirstAAb_1.time_12$pMat),
                              rownames(KEGG_GSA_consensus2$FirstAAb_1$FirstAAb_1.time_9$pMat),
                              rownames(KEGG_GSA_consensus2$FirstAAb_1$FirstAAb_1.time_6$pMat),
                              rownames(KEGG_GSA_consensus2$FirstAAb_1$FirstAAb_1.time_3$pMat),
                              rownames(KEGG_GSA_consensus2$FirstAAb_1$FirstAAb_1.time_0$pMat))))
rownames(allpaths)<-allpaths$unique.c.rownames.KEGG_GSA_consensus2.FirstAAb_1.FirstAAb_1.time_12.pMat...
n <- max(dim(GSA_ALL_consensus$FirstAAb_1$FirstAAb_1.time_0$pMat)[1],
         dim(GSA_ALL_consensus$FirstAAb_1$FirstAAb_1.time_3$pMat)[1],
         dim(GSA_ALL_consensus$FirstAAb_1$FirstAAb_1.time_6$pMat)[1],
         dim(GSA_ALL_consensus$FirstAAb_1$FirstAAb_1.time_9$pMat)[1],
         dim(GSA_ALL_consensus$FirstAAb_1$FirstAAb_1.time_12$pMat)[1])

#####
Dir_Mat = NULL
Pval_Mat = NULL

for(TL in TimeLabels){
  print(TL)
  Tix = grep(TL,names(AG_consensus))
  #print(Tix)
  if(identical(Tix,integer(0))){
    next
  }else{
    pMat = na.omit(AG_consensus[[names(AG_consensus)[Tix]]][["pMat"]])
    print(pMat)
    pMat2<-merge (allpaths,pMat,by=0, all.x=TRUE) 
    rownames(pMat2)<-pMat2[,1]; pMat2<-pMat2[,-c(1,2)]
    print(pMat2)
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
selection<-unique(pianoDF$GSET)

pianoDFGADAfirst<-pianoDF
pianoDFGADAfirst

#####################################################################################
# After discussion we determine that it is better to instead of having 5 color frames,
# One for each timepoint, it is better to have just one
AG<-"FirstAAb_1"
direction_classes = c("Distinct-directional (up)","Mixed-directional (up)","Non-directional","Mixed-directional (dn)","Distinct-directional (dn)")
directioncolors = c("red","red", "white","blue","blue")
LEGEND_IX = which(direction_classes%in%unique(pianoDFGADAfirst$DIRECTION))

direction_classes = direction_classes[LEGEND_IX]
pianoDFGADAfirst$DIRECTION = factor(pianoDFGADAfirst$DIRECTION, levels = direction_classes)


directioncolors = directioncolors[LEGEND_IX]
colorends = 1:(length(direction_classes)*2)
colorends[(1:(length(direction_classes)) * 2) - 1] = "white"
colorends[(1:(length(direction_classes)) * 2) ] = directioncolors

# Reescaling values so we can plot different categories
Nclasses = sort((unique(100 * (as.numeric(pianoDFGADAfirst$DIRECTION)-1 ))))

pianoDFGADAfirst$PVALresc = pianoDFGADAfirst$PVALcomp + (100 * (as.numeric(pianoDFGADAfirst$DIRECTION)-1 ))

scalerange <- range(pianoDFGADAfirst$PVALcomp)
gradientends <- scalerange + rep(Nclasses, each=2)

pianoDFGADAfirstSelection<-pianoDFGADAfirst[pianoDFGADAfirst$GSET %in% selection,]


# Heatmap #
pselection = ggplot(pianoDFGADAfirstSelection, aes(TIME, GSET)) + 
  geom_tile(aes(fill = PVALresc), colour = "transparent") + 
  #geom_tile(aes(fill = PVALresc), colour = ifelse(pianoDFALLselection$PVAL<0.3,"black","transparent"), size=1) +  
  scale_fill_gradientn(colours = colorends, values = rescale(gradientends)) + 
  scale_x_discrete("", expand = c(0, 0)) + 
  scale_y_discrete("", expand = c(0, 0)) + 
  #geom_text(aes(label = round(PVAL, 2))) +   #Turned off for noPval figure
  geom_text(aes(label = ifelse(PVAL<0.3,"*","")), size=10) + # To put asterisk in significant
  #theme_grey(base_size = 18) + 
  theme(legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text.x = element_text(angle = 330, hjust = 0))  + ggtitle(AG)

#pdf (file="pianoDFALLselectionLE0afterdiscussionNoPvalwithasterisks.pdf",height = 11,width = 9)
print(pselection)
#dev.off()
############################################
# Selection sorted by 
direction_classes = c("Distinct-directional (up)","Mixed-directional (up)","Non-directional","Mixed-directional (dn)","Distinct-directional (dn)")
pianoDFGADAfirstSelection$directionnumbers[pianoDFGADAfirstSelection$DIRECTION=="Distinct-directional (up)"]<-1
pianoDFGADAfirstSelection$directionnumbers[pianoDFGADAfirstSelection$DIRECTION=="Mixed-directional (up)"]<-1
pianoDFGADAfirstSelection$directionnumbers[pianoDFGADAfirstSelection$DIRECTION=="Non-directional"]<-2
pianoDFGADAfirstSelection$directionnumbers[pianoDFGADAfirstSelection$DIRECTION=="Distinct-directional (dn)"]<-3
pianoDFGADAfirstSelection$directionnumbers[pianoDFGADAfirstSelection$DIRECTION=="Mixed-directional (dn)"]<-3

# This sorting i liked it 
pianoDFGADAfirstSelectionSorted <- pianoDFGADAfirstSelection[order(pianoDFGADAfirstSelection$directionnumbers, pianoDFGADAfirstSelection$PVAL),]


#pianoDFALLselectionSorted
pianoDFGADAfirstSelectionSorted$GSET <- factor(pianoDFGADAfirstSelectionSorted$GSET, levels = rev(unique(as.character(pianoDFGADAfirstSelectionSorted$GSET))))
pselectionSorted = ggplot(pianoDFGADAfirstSelectionSorted, aes(TIME, GSET)) + 
  geom_tile(aes(fill = PVALresc), colour = "white") + 
  scale_fill_gradientn(colours = colorends, values = rescale(gradientends)) + 
  scale_x_discrete("", expand = c(0, 0)) + 
  scale_y_discrete("", expand = c(0, 0)) + 
  #geom_text(aes(label = round(PVAL, 2))) +
  geom_text(aes(label = ifelse(PVAL<0.3,"*","")), size=10) + # To put asterisk in significant
  theme_grey(base_size = 18) +
  #theme_grey(base_size = 5) +
  theme(legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text.x = element_text(angle = 330, hjust = 0))  + ggtitle(AG)

#pdf (file="pianoDFSMALLSUBSETSorted210917.pdf",height = 11,width = 9)
print(pselectionSorted)
#dev.off()
###########################################################
#########################
####    IAAFirst     ####
#########################
# Creating Piano Data Frame hacking script so we can play with it at will for IAA first
dim(KEGG_GSA_consensus2$FirstAAb_2$FirstAAb_2.time_12$pMat) # 104 Pathways
dim(KEGG_GSA_consensus2$FirstAAb_2$FirstAAb_2.time_9$pMat) # 306 Pathways
dim(KEGG_GSA_consensus2$FirstAAb_2$FirstAAb_2.time_6$pMat) # 306 Pathways
dim(KEGG_GSA_consensus2$FirstAAb_2$FirstAAb_2.time_3$pMat) # 110 Pathways
dim(KEGG_GSA_consensus2$FirstAAb_2$FirstAAb_2.time_0$pMat) # 109 Pathways

length(unique(c(rownames(KEGG_GSA_consensus2$FirstAAb_2$FirstAAb_2.time_12$pMat),
                rownames(KEGG_GSA_consensus2$FirstAAb_2$FirstAAb_2.time_9$pMat),
                rownames(KEGG_GSA_consensus2$FirstAAb_2$FirstAAb_2.time_6$pMat),
                rownames(KEGG_GSA_consensus2$FirstAAb_2$FirstAAb_2.time_3$pMat),
                rownames(KEGG_GSA_consensus2$FirstAAb_2$FirstAAb_2.time_0$pMat))))



AG="FirstAAb_2"
GSA_ALL_consensus = KEGG_GSA_consensus2
med_Pval = 0.8
AG_consensus = GSA_ALL_consensus[[AG]]
summary(AG_consensus)

TimeLabels = c("time_0","time_3","time_6","time_9","time_12")
#dim(GSA_ALL_consensus$ALL$g_time_0$pMat)[1]
allpaths<-data.frame(unique(c(rownames(KEGG_GSA_consensus2$FirstAAb_2$FirstAAb_2.time_12$pMat),
                              rownames(KEGG_GSA_consensus2$FirstAAb_2$FirstAAb_2.time_9$pMat),
                              rownames(KEGG_GSA_consensus2$FirstAAb_2$FirstAAb_2.time_6$pMat),
                              rownames(KEGG_GSA_consensus2$FirstAAb_2$FirstAAb_2.time_3$pMat),
                              rownames(KEGG_GSA_consensus2$FirstAAb_2$FirstAAb_2.time_0$pMat))))
rownames(allpaths)<-allpaths$unique.c.rownames.KEGG_GSA_consensus2.FirstAAb_2.FirstAAb_2.time_12.pMat...
n <- max(dim(GSA_ALL_consensus$FirstAAb_2$FirstAAb_2.time_0$pMat)[1],
         dim(GSA_ALL_consensus$FirstAAb_2$FirstAAb_2.time_3$pMat)[1],
         dim(GSA_ALL_consensus$FirstAAb_2$FirstAAb_2.time_6$pMat)[1],
         dim(GSA_ALL_consensus$FirstAAb_2$FirstAAb_2.time_9$pMat)[1],
         dim(GSA_ALL_consensus$FirstAAb_2$FirstAAb_2.time_12$pMat)[1])

#####
Dir_Mat = NULL
Pval_Mat = NULL

for(TL in TimeLabels){
  print(TL)
  Tix = grep(TL,names(AG_consensus))
  #print(Tix)
  if(identical(Tix,integer(0))){
    next
  }else{
    pMat = na.omit(AG_consensus[[names(AG_consensus)[Tix]]][["pMat"]])
    print(pMat)
    pMat2<-merge (allpaths,pMat,by=0, all.x=TRUE) 
    rownames(pMat2)<-pMat2[,1]; pMat2<-pMat2[,-c(1,2)]
    print(pMat2)
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
selection<-unique(pianoDF$GSET)

pianoDFIAAfirst<-pianoDF
pianoDFIAAfirst

#####################################################################################
# After discussion we determine that it is better to instead of having 5 color frames,
# One for each timepoint, it is better to have just one
AG<-"FirstAAb_2"
direction_classes = c("Distinct-directional (up)","Mixed-directional (up)","Non-directional","Mixed-directional (dn)","Distinct-directional (dn)")
directioncolors = c("red","red", "white","blue","blue")
LEGEND_IX = which(direction_classes%in%unique(pianoDFIAAfirst$DIRECTION))

direction_classes = direction_classes[LEGEND_IX]
pianoDFIAAfirst$DIRECTION = factor(pianoDFIAAfirst$DIRECTION, levels = direction_classes)


directioncolors = directioncolors[LEGEND_IX]
colorends = 1:(length(direction_classes)*2)
colorends[(1:(length(direction_classes)) * 2) - 1] = "white"
colorends[(1:(length(direction_classes)) * 2) ] = directioncolors

# Reescaling values so we can plot different categories
Nclasses = sort((unique(100 * (as.numeric(pianoDFIAAfirst$DIRECTION)-1 ))))

pianoDFIAAfirst$PVALresc = pianoDFIAAfirst$PVALcomp + (100 * (as.numeric(pianoDFIAAfirst$DIRECTION)-1 ))

scalerange <- range(pianoDFIAAfirst$PVALcomp)
gradientends <- scalerange + rep(Nclasses, each=2)

pianoDFIAAfirstSelection<-pianoDFIAAfirst[pianoDFIAAfirst$GSET %in% selection,]


# Heatmap #
pselection = ggplot(pianoDFIAAfirstSelection, aes(TIME, GSET)) + 
  geom_tile(aes(fill = PVALresc), colour = "transparent") + 
  #geom_tile(aes(fill = PVALresc), colour = ifelse(pianoDFALLselection$PVAL<0.3,"black","transparent"), size=1) +  
  scale_fill_gradientn(colours = colorends, values = rescale(gradientends)) + 
  scale_x_discrete("", expand = c(0, 0)) + 
  scale_y_discrete("", expand = c(0, 0)) + 
  #geom_text(aes(label = round(PVAL, 2))) +   #Turned off for noPval figure
  geom_text(aes(label = ifelse(PVAL<0.3,"*","")), size=10) + # To put asterisk in significant
  #theme_grey(base_size = 18) + 
  theme(legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text.x = element_text(angle = 330, hjust = 0))  + ggtitle(AG)

#pdf (file="pianoDFALLselectionLE0afterdiscussionNoPvalwithasterisks.pdf",height = 11,width = 9)
print(pselection)
#dev.off()
############################################
# Selection sorted by 
direction_classes = c("Distinct-directional (up)","Mixed-directional (up)","Non-directional","Mixed-directional (dn)","Distinct-directional (dn)")
pianoDFIAAfirstSelection$directionnumbers[pianoDFIAAfirstSelection$DIRECTION=="Distinct-directional (up)"]<-1
pianoDFIAAfirstSelection$directionnumbers[pianoDFIAAfirstSelection$DIRECTION=="Mixed-directional (up)"]<-1
pianoDFIAAfirstSelection$directionnumbers[pianoDFIAAfirstSelection$DIRECTION=="Non-directional"]<-2
pianoDFIAAfirstSelection$directionnumbers[pianoDFIAAfirstSelection$DIRECTION=="Distinct-directional (dn)"]<-3
pianoDFIAAfirstSelection$directionnumbers[pianoDFIAAfirstSelection$DIRECTION=="Mixed-directional (dn)"]<-3

# This sorting i liked it 
pianoDFIAAfirstSelectionSorted <- pianoDFIAAfirstSelection[order(pianoDFIAAfirstSelection$directionnumbers, pianoDFIAAfirstSelection$PVAL),]


#pianoDFALLselectionSorted
pianoDFIAAfirstSelectionSorted$GSET <- factor(pianoDFIAAfirstSelectionSorted$GSET, levels = rev(unique(as.character(pianoDFIAAfirstSelectionSorted$GSET))))
pselectionSorted = ggplot(pianoDFIAAfirstSelectionSorted, aes(TIME, GSET)) + 
  geom_tile(aes(fill = PVALresc), colour = "white") + 
  scale_fill_gradientn(colours = colorends, values = rescale(gradientends)) + 
  scale_x_discrete("", expand = c(0, 0)) + 
  scale_y_discrete("", expand = c(0, 0)) + 
  #geom_text(aes(label = round(PVAL, 2))) +
  geom_text(aes(label = ifelse(PVAL<0.3,"*","")), size=10) + # To put asterisk in significant
  theme_grey(base_size = 18) +
  #theme_grey(base_size = 5) +
  theme(legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text.x = element_text(angle = 330, hjust = 0))  + ggtitle(AG)

#pdf (file="pianoDFSMALLSUBSETSorted210917.pdf",height = 11,width = 9)
print(pselectionSorted)
#dev.off()

###########################################################
# Test to know how many significant pathways there are depending on the threshold
pianoDFIAAfirst
summary(pianoDFIAAfirst)
head(pianoDFIAAfirst)

####

Whicharethepathways (data=pianoDFIAAfirst, threshold = 0.05, Do_you_want_to_know="yes")

thrs<-c(0.01,0.05,0.1)
Whicharethepathways (data=pianoDFIAAfirst, threshold = thrs, Do_you_want_to_know="yes")




Unsupervisedselection (data=pianoDFGADAfirst, threshold=0.01,Heatmaptitle="GADA_First")
Unsupervisedselection (data=pianoDFIAAfirst, threshold=0.01,Heatmaptitle="IAA_First")
################################################################
# Now we are going to plot a set of figures

pdf (file="/home/leobalzano/Dropbox (UFL)/TEDDY/Paper1/NatureCommFormat/NCommV1/NCommV1ToShare/SupplementaryData/testing.pdf",height = 11,width = 9)
#### GADA 0.01
plot.new()
tulo<-Whicharethepathways2 (data=pianoDFGADAfirst, threshold = 0.01, Do_you_want_to_know="yes")
tulo2<-as.data.frame(tulo$Pathways)
mtext(tulo$Sigs)
for (i in 1:nrow(tulo2)) {
  mtext(paste(as.vector(tulo2[i, ]), collapse = "       "), line = -(i +1), adj = 0)
} 
Unsupervisedselection (data=pianoDFGADAfirst, threshold=0.01,Heatmaptitle="GADA_First alpha=0.01")
plot.new()

#### IAA 0.01
plot.new()
tulo<-Whicharethepathways2 (data=pianoDFIAAfirst, threshold = 0.01, Do_you_want_to_know="yes")
tulo2<-as.data.frame(tulo$Pathways)
mtext(tulo$Sigs)
for (i in 1:nrow(tulo2)) {
  mtext(paste(as.vector(tulo2[i, ]), collapse = "       "), line = -(i +1), adj = 0)
} 
Unsupervisedselection (data=pianoDFIAAfirst, threshold=0.01,Heatmaptitle="IAA_First alpha=0.01")
plot.new()

#### GADA 0.05
plot.new()
tulo<-Whicharethepathways2 (data=pianoDFGADAfirst, threshold = 0.05, Do_you_want_to_know="yes")
tulo2<-as.data.frame(tulo$Pathways)
mtext(tulo$Sigs)
for (i in 1:nrow(tulo2)) {
  mtext(paste(as.vector(tulo2[i, ]), collapse = "       "), line = -(i +1), adj = 0)
} 
Unsupervisedselection (data=pianoDFGADAfirst, threshold=0.05,Heatmaptitle="GADA_First alpha=0.05")
plot.new()
#### IAA 0.05
plot.new()
tulo<-Whicharethepathways2 (data=pianoDFIAAfirst, threshold = 0.05, Do_you_want_to_know="yes")
tulo2<-as.data.frame(tulo$Pathways)
mtext(tulo$Sigs)
for (i in 1:nrow(tulo2)) {
  mtext(paste(as.vector(tulo2[i, ]), collapse = "       "), line = -(i +1), adj = 0)
} 
Unsupervisedselection (data=pianoDFIAAfirst, threshold=0.05,Heatmaptitle="IAA_First alpha=0.05")
plot.new()


#### GADA 0.1
plot.new()
tulo<-Whicharethepathways2 (data=pianoDFGADAfirst, threshold = 0.1, Do_you_want_to_know="yes")
tulo2<-as.data.frame(tulo$Pathways)
mtext(tulo$Sigs)
for (i in 1:nrow(tulo2)) {
  mtext(paste(as.vector(tulo2[i, ]), collapse = "       "), line = -(i +1), adj = 0)
} 
Unsupervisedselection (data=pianoDFGADAfirst, threshold=0.1,Heatmaptitle="GADA_First alpha=0.1")
plot.new()
#### IAA 0.1
plot.new()
tulo<-Whicharethepathways2 (data=pianoDFIAAfirst, threshold = 0.1, Do_you_want_to_know="yes")
tulo2<-as.data.frame(tulo$Pathways)
mtext(tulo$Sigs)
for (i in 1:nrow(tulo2)) {
  mtext(paste(as.vector(tulo2[i, ]), collapse = "       "), line = -(i +1), adj = 0)
} 
Unsupervisedselection (data=pianoDFIAAfirst, threshold=0.1,Heatmaptitle="IAA_First alpha=0.1")
plot.new()

####GADA 0.2
plot.new()
tulo<-Whicharethepathways2 (data=pianoDFGADAfirst, threshold = 0.2, Do_you_want_to_know="yes")
tulo2<-as.data.frame(tulo$Pathways)
mtext(tulo$Sigs)
for (i in 1:nrow(tulo2)) {
  mtext(paste(as.vector(tulo2[i, ]), collapse = "       "), line = -(i +1), adj = 0)
} 
Unsupervisedselection (data=pianoDFGADAfirst, threshold=0.2,Heatmaptitle="GADA_First alpha=0.2")
plot.new()

#### IAA 0.2
plot.new()
tulo<-Whicharethepathways2 (data=pianoDFIAAfirst, threshold = 0.2, Do_you_want_to_know="yes")
tulo2<-as.data.frame(tulo$Pathways)
mtext(tulo$Sigs)
for (i in 1:nrow(tulo2)) {
  mtext(paste(as.vector(tulo2[i, ]), collapse = "       "), line = -(i +1), adj = 0)
} 
Unsupervisedselection (data=pianoDFIAAfirst, threshold=0.2,Heatmaptitle="IAA_First alpha=0.2")
plot.new()
####


dev.off()