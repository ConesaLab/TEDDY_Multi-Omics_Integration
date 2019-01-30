###########################################################
###############     loadTEDDYtools.R     ##################
###########################################################
# Authors: Ricardo Ramirez Flores
# Genetics Institute, University of Florida (Gainesville)

# TEDDY tools: Requirements
# Loads all the functions, objects and collection of functions of TEDDYtools

###########################################################
setwd("/media/data/leobalzano/ScriptsForTEDDY")

# The location of TEDDYtoolsV2
#setwd("/home/rramirez/Dropbox (UFL)/myTEDDY_V2/TEDDYtoolsV2/")

# Global Variables : RAW DATA

# Gene Expression

load("Data/Piano/GlobalData/GEXtargets.ro") #Target Matrix
load("Data/Piano/GlobalData/rawTEDDY_GEX.ro") #Count Matrix (Original Table)
load("Data/Piano/GlobalData/rawTEDDY_GEX_inversevst.ro") #New Matrix

# Metabolomics

load("Data/Piano/GlobalData/MetabolicTargets.ro") #Target Matrix
load("Data/Piano/GlobalData/Metabolic_Counts_Raw.ro") #Complete Matrix
load("Data/Piano/GlobalData/GCTOFRAW_counts.ro") #GCTOF
load("Data/Piano/GlobalData/negLipRAW_counts.ro") #NegLip
load("Data/Piano/GlobalData/posLipRAW_counts.ro") #PosLip

# GeneSets

load("Data/Piano/GlobalData/TEDDY_geneSets.ro")


# External Libraries
library(dplyr)
library(data.table)
library(illuminaHumanv4.db)
library(GO.db)
library(splines)
library(limma)
library(piano)
library(matrixStats)
library(maSigPro)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(gridExtra)
library(cluster)
library(gtools)
library(stringdist)
library(fpc)
library(mclust)
library(scales)

# Load TEDDYtools Suite

print("Loading Processing and filtering tools")
source("Data/Piano/Suite/Processing_Filtering_TT.R")
print("Loading Count Manipulation tools")
source("Data/Piano/Suite/CountManipulation_TT.R")
print("Loading GSEA tools")
source("Data/Piano/Suite/GSEA_TT.R")
print("Loading Linear Modelling tools")
source("Data/Piano/Suite/LinearModels_DEA_TT.R")
print("Loading Visual Manager tools")
source("Data/Piano/Suite/VisualManager_TT.R")
print("Loading Transcriptional Signature Comparison tools")
source("Data/Piano/Suite/TranscriptionalSignComp_TT.R")
print("Loading Clustering tools")
source("Data/Piano/Suite/Clustering_TT.R")

