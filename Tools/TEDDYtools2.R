###################################################
###############     TEDDY TOOLS     ###############
###################################################
# Author: Ricardo Ramirez Flores, Leandro Balzano-Nogueira
# Genetics Institute, University of Florida (Gainesville)

# Set of R functions for TEDDY analysis 

###################################################
setwd("/media/data/leobalzano/ScriptsForTEDDY")
getwd()

###################################################
# Set of libraries and R objects needed to run linear model analyses
###################################################

library(dplyr)
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




load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/targetFile.ro") #Load sample annotations (targets)
load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/rawTEDDYexprdata.ro") #Load RAW expression data ()
load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/GLSAAblimma.ro") #Load Gene Level Statistics of limma models in all samples and in antibody specific samples
#load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/MetabolicTargets.ro") #Load Metabolomic annotations
#load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/MetabolicCounts.ro") #Load Metabolomic counts
load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/GLSAAbcomplimma.ro") #Load Gene Level Statistics of limma models in all samples and in antibody specific samples, antibody comparison
load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/global_summary.ro") #Load Gene Level Statistics of limma models in all samples and in antibody specific samples, antibody comparison
load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/DYN_summaryres.ro") #Load summary lists for all time-specific expression analyses using dynamics instead of time points
load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/summary_res.ro") #Load Gene Level Statistics of limma models in all samples and in antibody specific samples, antibody comparison
load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/summary_res_woAG.ro") #Load summary lists for all time-specific expression analyses using time-points without age-group distinction
load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/GSA_All_ResultsList.ro") #Load a list containing Piano analysis with 7 methods for all the 90 models #Origin: FunctionalCharacSignR
load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/GSA_All_ResultsList_consensus.ro") #Load a list containing the results of the consensus Piano analysis done to the 90 models
load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/GSA_All_ResultsList_consensus_significant.ro") #Load a list containing the results of the consensus Piano analysis done to the 90 models grouped by TAG (20 major subgroups)
load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/AllTEDDY_signatures.ro") #Load all signatures #To be modified in case samples are taken
load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/MatrixLMSignatures_ES.ro") #Load EnrichmentScore Matrix, with leading edges...
#Metabolic Data
load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/negLip_counts.ro")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/negLip_annotations.ro")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/posLip_counts.ro")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/posLip_annotations.ro")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/GCTOF_counts.ro")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/GCTOF_annotations.ro")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/Metabolic_Counts_Raw.ro")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/TEDDY_toolsLeo/Metabolic_targets.ro")


# TARGET MANIPULATION: Tools to select patient-specific and time-specific data // Tools to add new data to original target tables
#
#


loadTargets = function(SELdisease,time_point="All",complete_case=NULL,targets){
  #
  # Function that provides a target table of the cases you are interested
  # call:
  # loadTargets(SELdisease = c("IA","T1D"), time_point = c("All", vector with specific times -integers-), complete_case = c("strict","flexible",NULL))
  # 
  # complete_case: It allows you to get only the samples with a corresponding control in a given time 
  # if NULL: this feature is deactivated
  # if flexible: it keeps all the patients that have at leats 1 case-control for the defined time points 
  # if strict: it keeps the patients that have 1 case-control in all specified time points
  #
  # Example:
  # test = loadTargets(SELdisease = "IA",time_point = "All",complete_case = "strict")
  #
  
  #Select disease
  selTargets = filter(targets,disease==SELdisease)
  
  #Select time point
  if(time_point=="All"){
    time_point= c(12,9,6,3,0)
    selTargets =  filter(selTargets,time %in% time_point)
  }else{
    selTargets =  filter(selTargets,time %in% time_point)
  }
  
  #Select data with complete information
  if(is.null(complete_case)){
    return(selTargets)
  }else{
    smids = c()
    for(case_control in unique(selTargets$CASE_IND)){
      strictFlag = TRUE
      cases = filter(selTargets,CASE_IND==case_control,outcome==1)
      controls = filter(selTargets,CASE_IND==case_control,outcome==0)
      ccsmids = c()
      
      for(t in time_point){
        
        caseFlag = t %in% cases$time
        controlFlag = t %in% controls$time
        time_flag = caseFlag & controlFlag
        
        if(time_flag){
          casetmp = filter(cases,time==t)
          controltmp = filter(controls,time==t)
          ccsmids = c(ccsmids,casetmp$sample_mask_id,controltmp$sample_mask_id)
        }else{
          if(complete_case == "strict"){
            strictFlag = FALSE
            break
          }
        }
      }
      
      if(strictFlag){
        smids = c(smids,ccsmids)
      }
      
    }
    
    selTargets =  filter(selTargets,sample_mask_id %in% smids)
    selTargets = arrange(selTargets,CASE_IND,time,outcome)
    if(sum(duplicated(selTargets$sample_mask_id))>0){
      bids = selTargets$sample_mask_id[duplicated(selTargets$sample_mask_id)]
      bcind = unique((filter(selTargets,sample_mask_id %in% bids))$CASE_IND)
      selTargets = filter(selTargets,!(CASE_IND %in% bcind))
    }
    
    return(selTargets)
  }
}

addFeatures = function(newDFdata,TARGETS,selCol){
  #
  # Generates a new target table with a new column of data
  #
  # inputs:
  # newDFdata = New -data frame- to merge
  # TARGETS = Original target table (labels as in Teddy website)
  # selCol = name of the column used for merging
  #
  # output:
  # New target file
  #
  newTARGETS = left_join(TARGETS,newDFdata,by =selCol)
  return(newTARGETS)
}

#
#
# Expression Data Processing:
# Tools for count normalisation and gene expression processing
#
#

normGEX = function(expMAT,targetTable){
  #
  # Normalises gene expression by substracting for each time specific case-control case, the mean signal
  # Note: Only works to normalise the original 5 time points (12,9,6,3) 
  # 
  # inputs:
  # -expMAT: Any count matrix organised by case-control
  # -targetTable = Target table that corresponds to the expression matrix (time must be a column in this table as)
  #
  # outputs:
  # -Normalized expression matrix
  #
  for(case_control in unique(targetTable$CASE_IND)){
    subsetTargets = filter(targetTable,CASE_IND==case_control)
    for(t in c(12,9,6,3,0)){
      if(t %in% subsetTargets$time){
        sids_t =  as.character((filter(targetTable,CASE_IND==case_control,time==t))$sample_mask_id)
        subsetMat = expMAT[,sids_t]
        subsetMat = subsetMat - rowMeans(subsetMat)
        expMAT[,sids_t] = subsetMat
      }
    }
  }
  return(expMAT)
}

normGEX_v2 = function(expMAT,targetTable,time = c(12,9,6,3,0)){
  #
  # Normalises gene expression by substracting for each time specific case-control case, the mean signal
  # Note: Only works to normalise the original 5 time points (12,9,6,3) 
  # 
  # inputs:
  # -expMAT: Any count matrix organised by case-control
  # -targetTable = Target table that corresponds to the expression matrix (time must be a column in this table as)
  #
  # outputs:
  # -Normalized expression matrix
  #
  for(case_control in unique(targetTable$CASE_IND)){
    subsetTargets = filter(targetTable,CASE_IND==case_control)
    for(t in time){
      if(t %in% subsetTargets$time){
        
        case_id = as.character((filter(targetTable,CASE_IND==case_control,time==t,outcome==1))$sample_mask_id)
        control_id = as.character((filter(targetTable,CASE_IND==case_control,time==t,outcome==0))$sample_mask_id)
        
        if(length(control_id)>1){
          control_mean = rowMeans(expMAT[,control_id],na.rm = TRUE)
          mean_mat = cbind(expMAT[,case_id],control_mean)
          
          #sids_t =  as.character((filter(targetTable,CASE_IND==case_control,time==t))$sample_mask_id)
          subsetMat = expMAT[,c(case_id,control_id)]
          subsetMat = subsetMat - rowMeans(mean_mat,na.rm = TRUE)
          expMAT[,c(case_id,control_id)] = subsetMat
        } else if(length(control_id)==1){
          
          sids_t =  c(case_id,control_id)
          subsetMat = expMAT[,sids_t]
          subsetMat = subsetMat - rowMeans(subsetMat,na.rm = TRUE)
          expMAT[,sids_t] = subsetMat
          
        }
        
      }
    }
  }
  return(expMAT)
}


gnormGEX = function(expMAT,targetTable){
  #
  # Normalizes gene expression by substracting for each outcome group, the mean signal
  # call:
  # normGEX(expMAT=expressionMatrix,targetTable=target table generated by loadTargets or a table with CASE_IND,time columns)
  # 
  # inputs:
  # -expMAT: Any count matrix organised by case-control
  # -targetTable = Target table that corresponds to the expression matrix (time must be a column in this table as)
  #
  # outputs:
  # -Normalized expression matrix
  #
  for(oc in unique(targetTable$outcome)){
    subsetTargets = filter(targetTable,outcome==oc)
    for(t in c(12,9,6,3,0)){
      if(t %in% subsetTargets$time){
        sids_t =  as.character((filter(subsetTargets,time==t))$sample_mask_id)
        subsetMat = expMAT[,sids_t]
        subsetMat = subsetMat - rowMeans(subsetMat)
        expMAT[,sids_t] = subsetMat
      }
    }
  }
  return(expMAT)
}

annotateGEX = function(expMAT,annotation_set){
  #
  # Annotates a gene expression matrix and takes the mean of repeated symbols
  # inputs:
  # -expMAT = expression matrix
  # -annotation_set = your microarray database (AnnotationDBi)
  #
  # output:
  # -Annotated expression matrix
  #
  probe_names = rownames(expMAT)
  probe_symbol = na.omit(data.frame(mapIds(annotation_set, keys=probe_names, keytype = "PROBEID", column="SYMBOL",multiVals = "first"),stringsAsFactors = F))
  #FIlter exprmat with known probes
  expMAT = expMAT[rownames(probe_symbol),]
  rownames(expMAT) = probe_symbol[,1]
  expMAT = avereps(expMAT,ID=rownames(expMAT))
  return(expMAT)
}

getExprmat = function(RAWEXPRMAT, STARGETS, ciNORM = TRUE, gNORM = FALSE){
  #
  # Gets an expression matrix for a given target file (specific tables)
  # input:
  # RAWEXPRMAT: Non annotated expression matrix
  # STARGETS: target file of selected samples
  # ciNORM: case-control normalization (TRUE-FALSE)
  # gNORM: group(outcome) normalization (TRUE-FALSE)
  # 
  # output:
  # -FIltered and processed expression matrix
  #
  
  # Filter expression matrix
  expmat = RAWEXPRMAT[,as.character(STARGETS$sample_mask_id)]
  
  # Normalization (substracting mean case-control signal to values)
  if(ciNORM){
    expmat = normGEX(expMAT=expmat,targetTable=STARGETS) 
  }
  # Normalization by group (Sick-Healthy)
  if(gNORM){
    expmat = gnormGEX(expMAT=expmat,targetTable=STARGETS) 
  }
  # Annotation of matrices
  expmat = annotateGEX(expMAT = expmat,annotation_set = illuminaHumanv4.db)
  
  return(expmat)
}

annotateGEX_v2 = function(expMAT, annotation_set, annotation_column){
  #
  # Annotates a gene expression matrix and takes the mean of repeated symbols
  # inputs:
  # -expMAT = expression matrix
  # -annotation_set = your microarray database (AnnotationDBi)
  #
  # output:
  # -Annotated expression matrix
  #
  probe_names = rownames(expMAT)
  probe_symbol = na.omit(data.frame(mapIds(annotation_set, keys=probe_names, keytype = "PROBEID", column = annotation_column,multiVals = "first"),stringsAsFactors = F))
  #FIlter exprmat with known probes
  expMAT = expMAT[rownames(probe_symbol),]
  rownames(expMAT) = probe_symbol[,1]
  expMAT = avereps(expMAT,ID=rownames(expMAT))
  return(expMAT)
}

getExprmat_v2 = function(RAWEXPRMAT, STARGETS, ciNORM = TRUE, gNORM = FALSE, annotation_column = "SYMBOL"){
  #
  # Gets an expression matrix for a given target file (specific tables)
  # input:
  # RAWEXPRMAT: Non annotated expression matrix
  # STARGETS: target file of selected samples
  # ciNORM: case-control normalization (TRUE-FALSE)
  # gNORM: group(outcome) normalization (TRUE-FALSE)
  # 
  # output:
  # -FIltered and processed expression matrix
  #
  
  # Filter expression matrix
  expmat = RAWEXPRMAT[,as.character(STARGETS$sample_mask_id)]
  
  # Normalization (substracting mean case-control signal to values)
  if(ciNORM){
    expmat = normGEX(expMAT=expmat,targetTable=STARGETS) 
  }
  # Normalization by group (Sick-Healthy)
  if(gNORM){
    expmat = gnormGEX(expMAT=expmat,targetTable=STARGETS) 
  }
  # Annotation of matrices
  expmat = annotateGEX_v2(expMAT = expmat,annotation_set = illuminaHumanv4.db, annotation_column = annotation_column)
  
  return(expmat)
}

FilterTargets = function(targetTable,query_class,query_values,OR = FALSE){
  #
  # A filter function that gives you the rows from targetTable that have all the features defined in query_class/query_values (OR=FALSE)
  # or at least one feature (OR=TRUE)
  #
  
  # Make individuals
  if(length(query_class) == length(query_values)){
    EvaluationMat = c()
    for(i in 1:length(query_class)){
      EvalVector = as.character(targetTable[,query_class[i]]) == query_values[i]
      EvaluationMat = cbind(EvaluationMat,EvalVector)
    }
  }else{
    print("Length of query class and query values are not the same")
    return()
  }
  
  if(OR==FALSE){
    
    EvalVector = rowSums(EvaluationMat)
    EvalVector = EvalVector == length(query_values)
    return(targetTable[EvalVector,])
    
  }else{
    EvalVector = rowSums(EvaluationMat)
    EvalVector = EvalVector > 0
    return(targetTable[EvalVector,])
  }
}

###########################################
# Generalization of TEDDYtools processing #
###########################################

loadTargets_v2 = function(SELdisease,time_point="All",complete_case=NULL,targets){
  #
  # Function that provides a target table of the cases you are interested
  # call:
  # loadTargets(SELdisease = c("IA","T1D"), time_point = c("All", vector with specific times -integers-), complete_case = c("strict","flexible",NULL))
  # 
  # complete_case: It allows you to get only the samples with a corresponding control in a given time 
  # if NULL: this feature is deactivated
  # if flexible: it keeps all the patients that have at leats 1 case-control for the defined time points 
  # if strict: it keeps the patients that have 1 case-control in all specified time points
  #
  # Example:
  # test = loadTargets(SELdisease = "IA",time_point = "All",complete_case = "strict")
  #
  
  #Select disease
  selTargets = filter(targets,disease==SELdisease)
  
  #Select time point
  if(time_point[1]=="All"){
    time_point= unique(targets$time)
    selTargets =  filter(selTargets,time %in% time_point)
  }else{
    selTargets =  filter(selTargets,time %in% time_point)
  }
  
  #Select data with complete information
  if(is.null(complete_case)){
    return(selTargets)
  }else{
    smids = c()
    for(case_control in unique(selTargets$CASE_IND)){
      strictFlag = TRUE
      cases = filter(selTargets,CASE_IND==case_control,outcome==1)
      controls = filter(selTargets,CASE_IND==case_control,outcome==0)
      ccsmids = c()
      
      for(t in time_point){
        
        caseFlag = t %in% cases$time
        controlFlag = t %in% controls$time
        time_flag = caseFlag & controlFlag
        
        if(time_flag){
          casetmp = filter(cases,time==t)
          controltmp = filter(controls,time==t)
          ccsmids = c(ccsmids,casetmp$sample_mask_id,controltmp$sample_mask_id)
        }else{
          if(complete_case == "strict"){
            strictFlag = FALSE
            break
          }
        }
      }
      
      if(strictFlag){
        smids = c(smids,ccsmids)
      }
      
    }
    
    selTargets =  filter(selTargets,sample_mask_id %in% smids)
    selTargets = arrange(selTargets,CASE_IND,time,outcome)
    if(sum(duplicated(selTargets$sample_mask_id))>0){
      bids = selTargets$sample_mask_id[duplicated(selTargets$sample_mask_id)]
      bcind = unique((filter(selTargets,sample_mask_id %in% bids))$CASE_IND)
      selTargets = filter(selTargets,!(CASE_IND %in% bcind))
    }
    
    return(selTargets)
  }
}

###########################################

getExprmat_v3 = function(RAWEXPRMAT, STARGETS, ciNORM = TRUE, gNORM = FALSE, time,annotation_column = "SYMBOL"){
  #
  # Gets an expression matrix for a given target file (specific tables)
  # input:
  # RAWEXPRMAT: Non annotated expression matrix
  # STARGETS: target file of selected samples
  # ciNORM: case-control normalization (TRUE-FALSE)
  # gNORM: group(outcome) normalization (TRUE-FALSE)
  # 
  # output:
  # -FIltered and processed expression matrix
  #
  
  # Filter expression matrix
  expmat = RAWEXPRMAT[,as.character(STARGETS$sample_mask_id)]
  
  # Normalization (substracting mean case-control signal to values)
  if(ciNORM){
    expmat = normGEX_v2(expMAT=expmat,targetTable=STARGETS,time = time) 
  }
  # Normalization by group (Sick-Healthy)
  if(gNORM){
    expmat = gnormGEX_v2(expMAT=expmat,targetTable=STARGETS,timepoints = time) 
  }
  # Annotation of matrices
  expmat = annotateGEX_v2(expMAT = expmat,annotation_set = illuminaHumanv4.db, annotation_column = annotation_column)
  
  return(expmat)
}

normGEX_v2 = function(expMAT,targetTable,time = c(12,9,6,3,0)){
  #
  # Normalises gene expression by substracting for each time specific case-control case, the mean signal
  # X WITHIN
  # 
  # inputs:
  # -expMAT: Any count matrix organised by case-control
  # -targetTable = Target table that corresponds to the expression matrix (time must be a column in this table as)
  #
  # outputs:
  # -Normalized expression matrix
  #
  for(case_control in unique(targetTable$CASE_IND)){
    subsetTargets = filter(targetTable,CASE_IND==case_control)
    for(t in time){
      if(t %in% subsetTargets$time){
        
        case_id = as.character((filter(targetTable,CASE_IND==case_control,time==t,outcome==1))$sample_mask_id)
        control_id = as.character((filter(targetTable,CASE_IND==case_control,time==t,outcome==0))$sample_mask_id)
        
        # Verify how many controls you have, just to decide which version of the normalization is going to be used
        
        if(length(control_id)>1){
          # Get the means of the controls, and then the mean of means
          control_mean = rowMeans(expMAT[,control_id],na.rm = TRUE)
          mean_mat = cbind(expMAT[,case_id],control_mean)
          
          #sids_t =  as.character((filter(targetTable,CASE_IND==case_control,time==t))$sample_mask_id)
          subsetMat = expMAT[,c(case_id,control_id)]
          subsetMat = subsetMat - rowMeans(mean_mat,na.rm = TRUE)
          expMAT[,c(case_id,control_id)] = subsetMat
        } else if(length(control_id)==1){
          
          sids_t =  c(case_id,control_id)
          subsetMat = expMAT[,sids_t]
          subsetMat = subsetMat - rowMeans(subsetMat,na.rm = TRUE)
          expMAT[,sids_t] = subsetMat
          
        }
        
      }
    }
  }
  return(expMAT)
}

gnormGEX_v2 = function(expMAT,targetTable,timepoints){
  #
  # Normalizes gene expression by substracting for each outcome group, the mean signal
  # call:
  # normGEX(expMAT=expressionMatrix,targetTable=target table generated by loadTargets or a table with CASE_IND,time columns)
  # 
  # inputs:
  # -expMAT: Any count matrix organised by case-control
  # -targetTable = Target table that corresponds to the expression matrix (time must be a column in this table as)
  #
  # outputs:
  # -Normalized expression matrix
  #
  for(oc in unique(targetTable$outcome)){
    subsetTargets = filter(targetTable,outcome==oc)
    for(t in timepoints){
      if(t %in% subsetTargets$time){
        sids_t =  as.character((filter(subsetTargets,time==t))$sample_mask_id)
        subsetMat = expMAT[,sids_t]
        subsetMat = subsetMat - rowMeans(subsetMat)
        expMAT[,sids_t] = subsetMat
      }
    }
  }
  return(expMAT)
}


############################################################################################################
# Performs gene-specific linear models to identify differences in dynamics between sick and healthy patients
# call:
# annotateGEX(expMAT = expression matrix,annotation_set = your microarray database)
#
# Example:
# test = annotateGEX(expMAT = expmat,annotation_set = illuminaHumanv4.db)
#

genedynamics = function(EXPRMAT, GENES, STARGETS, TRSH = 0.05, SIGNIFICANCE = FALSE,CASE_CONTROL=FALSE){
  #Load libraries
  library(gridExtra)
  library(ggplot2)
  
  EXPRMAT = t(EXPRMAT)
  GENES = GENES[GENES%in%colnames(EXPRMAT)]
  
  if(!identical(GENES, character(0))){
    #Create data.frame with selected genes
    genedframe = cbind(STARGETS,EXPRMAT[as.character(STARGETS$sample_mask_id),GENES])
    colnames(genedframe)[(ncol(genedframe)-(length(GENES)-1)):ncol(genedframe)] = GENES
    genedframe$mask_id = as.factor(genedframe$mask_id)
    genedframe$outcome = as.factor(ifelse(genedframe$outcome==0,"Healthy","Sick"))
    genedframe$CASE_IND = as.factor(as.character(genedframe$CASE_IND))
    #For each gene at each time adjust a linear model and save results
    lm_results = c()
    for(GENE in GENES){
      for(t in unique(genedframe$time)){
        t_GEX = genedframe[genedframe$time==t,]
        t_lm = summary(lm(t_GEX[,GENE]~t_GEX[,"outcome"]), data = t_GEX)
        coef = coefficients(t_lm)[2,1]
        pval = coefficients(t_lm)[2,4]
        lm_results =  rbind(lm_results,c(t,GENE,coef,pval))
      }
    }
    #Manipulate pval results
    colnames(lm_results) = c("TIME","GENE","COEF","PVAL")
    lm_results = data.frame(lm_results,stringsAsFactors = F)
    lm_results[,2] = as.factor(lm_results[,2]) 
    lm_results[,3] = as.numeric(lm_results[,3])
    lm_results[,4] = p.adjust(p=as.numeric(lm_results[,4]),method = "fdr")
    
    for(GENE in GENES){
      print(GENE)
      significant = FALSE
      gene_results = lm_results[lm_results$GENE==GENE,]
      
      if(SIGNIFICANCE==TRUE){
        if(sum(gene_results$PVAL < TRSH) > 0){
          significant = TRUE
        }
      }else{
        significant = TRUE
      }
      
      gene_results$PVAL = -log(gene_results$PVAL)
      
      if(significant){
        
        vp = ggplot(data=gene_results, aes(x=COEF, y=PVAL, label=TIME),guide=FALSE) + geom_point(colour="white", fill="lightblue", shape=21) +
          geom_text(size=5,colour="black") + scale_size_continuous(range = c(8, 18)) + scale_x_continuous(name="Coefficient") +
          scale_y_continuous(name="-log(pval)") + geom_hline(yintercept=-log(TRSH), linetype="dashed") + theme(plot.margin = unit(c(1,3.3,1,1), "cm"),axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12))
        
        if(CASE_CONTROL){
          GEX_plot <- ggplot(genedframe, aes(x = genedframe[,"time"] * -1, y = genedframe[,GENE], color = CASE_IND, group = mask_id, linetype = outcome)) + 
            labs(y = GENE) + geom_point() + geom_line() + scale_x_continuous(name="Time",breaks=c(-12,-9,-6,-3,0)) + ylab("Expression") + theme(plot.margin = unit(c(1,0.6,1,0.6), "cm"), axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12))
        }else{
          GEX_plot <- ggplot(genedframe, aes(x = genedframe[,"time"] * -1, y = genedframe[,GENE], color = outcome, group = mask_id)) + 
            labs(y = GENE) + geom_point() + geom_line() + scale_x_continuous(name="Time",breaks=c(-12,-9,-6,-3,0)) + ylab("Expression") + theme(plot.margin = unit(c(1,0.6,1,0.6), "cm"),axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12))
        }
        
        fgenedframe = genedframe
        fgenedframe$time = factor(genedframe$time,levels=c(12,9,6,3,0))
        gp = ggplot(fgenedframe, aes(x = time, y = fgenedframe[,GENE])) + geom_boxplot(aes(fill = outcome), alpha = 0.5) + 
          stat_summary(fun.y=mean, geom="line", aes(group=outcome, col=outcome),size=1.5)  + ylab("Expression") +
          theme(axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12), plot.margin = unit(c(1,0.6,1,0.6), "cm"))
        
        
        grid.arrange(GEX_plot,gp,vp, nrow=3,top = GENE)
        
      }
    }
    
    return(lm_results)
  }
  print("We failed to find genes in expression matrix") 
  return(NULL)
}

############################################################################################################
# Performs gene-specific linear models to identify differences in dynamics between sick and healthy patients
# call:
# annotateGEX(expMAT = expression matrix,annotation_set = your microarray database)
#
# Example:
# test = annotateGEX(expMAT = expmat,annotation_set = illuminaHumanv4.db)
#

plotgenetrajectory = function(EXPRMAT, GENES, STARGETS, CASE_CONTROL=FALSE){
  
  #Load libraries
  library(gridExtra)
  library(ggplot2)
  
  EXPRMAT = t(EXPRMAT)
  GENES = GENES[GENES%in%colnames(EXPRMAT)]
  
  if(!identical(GENES, character(0))){
    
    #Create data.frame with selected genes
    genedframe = cbind(STARGETS,EXPRMAT[as.character(STARGETS$sample_mask_id),GENES])
    colnames(genedframe)[(ncol(genedframe)-(length(GENES)-1)):ncol(genedframe)] = GENES
    genedframe$mask_id = as.factor(genedframe$mask_id)
    genedframe$outcome = as.factor(ifelse(genedframe$outcome==0,"Healthy","Sick"))
    genedframe$CASE_IND = as.factor(as.character(genedframe$CASE_IND))
    
    for(GENE in GENES){
      
      if(CASE_CONTROL){
        GEX_plot <- ggplot(genedframe, aes(x = genedframe[,"time"] * -1, y = genedframe[,GENE], color = CASE_IND, group = mask_id, linetype = outcome)) + 
          labs(y = GENE) + geom_point() + geom_line() + theme(plot.margin = unit(c(1,1,1,1), "cm"))
        print(GEX_plot)
      }else{
        GEX_plot <- ggplot(genedframe, aes(x = genedframe[,"time"] * -1, y = genedframe[,GENE], color = outcome, group = mask_id)) + 
          labs(y = GENE) + geom_point() + geom_line() + theme(plot.margin = unit(c(1,1,1,1), "cm"))
        print(GEX_plot)
      }
      
    }
    
  }
  else{
    print("We failed to find genes in expression matrix") 
    return(NULL)
  }
}

#######################################

plotgenetrajectorybp = function(EXPRMAT,STARGETS,GENES){
  
  EXPRMAT = t(EXPRMAT)
  GENES = GENES[GENES%in%colnames(EXPRMAT)]
  
  if(!identical(GENES, character(0))){
    
    #Create data.frame with selected genes
    genedframe = cbind(STARGETS,EXPRMAT[as.character(STARGETS$sample_mask_id),GENES])
    colnames(genedframe)[(ncol(genedframe)-(length(GENES)-1)):ncol(genedframe)] = GENES
    genedframe$mask_id = as.factor(genedframe$mask_id)
    genedframe$outcome = as.factor(ifelse(genedframe$outcome==0,"Healthy","Sick"))
    genedframe$CASE_IND = as.factor(as.character(genedframe$CASE_IND))
    genedframe$time = factor(genedframe$time,levels=c(12,9,6,3,0))
    
    for(GENE in GENES){
      
      gp = ggplot(genedframe, aes(x = time, y = genedframe[,GENE])) + geom_boxplot(aes(fill = outcome), alpha = 0.5) + 
        stat_summary(fun.y=median, geom="line", aes(group=outcome, col=outcome),size=1.5)  + labs(title = GENE) + ylab("Expression") +
        theme(axis.title.x = element_text(size = 14,face="bold"), axis.title.y = element_text(size = 14,face="bold"),axis.text.x  = element_text(size=14),axis.text.y  = element_text(size=14))
      print(gp)
    }
    
  }
  else{
    print("We failed to find genes in expression matrix") 
    return(NULL)
  }
}

########################

timelm = function(EXPRMAT, GENES, STARGETS, YRESP = "outcome"){
  #
  # This function adjusts a time-specific linear model to identify significant changes in expression/presence 
  # of a given feature between various responses (usually outcome)
  #
  # input:
  # -EXPRMAT: Any count matrix with features in rows and sample_mask_ids in columns
  # -STARGETS: Targets that correspond to the Expression Mat
  # -GENES: a group of genes of interest
  # -YRESP: Response variable
  #
  # output:
  # -Results from the linear model, with corrected pvalue
  #
  genedframe = cbind(STARGETS,EXPRMAT[GENES,as.character(STARGETS$sample_mask_id)])
  colnames(genedframe)[(ncol(genedframe)-(length(GENES)-1)):ncol(genedframe)] = GENES
  genedframe$mask_id = as.factor(genedframe$mask_id)
  genedframe$outcome = as.factor(ifelse(genedframe$outcome==0,"Healthy","Sick"))
  genedframe$CASE_IND = as.factor(as.character(genedframe$CASE_IND))
  #For each gene at each time adjust a linear model and save results
  lm_results = c()
  for(GENE in GENES){
    for(t in unique(genedframe$time)){
      t_GEX = genedframe[genedframe$time==t,]
      t_lm = summary(lm(t_GEX[,GENE]~t_GEX[,YRESP]), data = t_GEX)
      coef = coefficients(t_lm)[2,1]
      pval = coefficients(t_lm)[2,4]
      lm_results =  rbind(lm_results,c(t,GENE,coef,pval))
    }
  }
  #Manipulate pval results
  colnames(lm_results) = c("TIME","GENE","COEF","PVAL")
  lm_results = data.frame(lm_results,stringsAsFactors = F)
  lm_results[,2] = as.factor(lm_results[,2]) 
  lm_results[,3] = as.numeric(lm_results[,3])
  lm_results[,4] = p.adjust(p=as.numeric(lm_results[,4]),method = "fdr")
  
  return(lm_results)
}

plottimelm = function(gene_results,TRSH,GENE){
  #
  # Plots a volcano plot of pvals and coefficients
  #
  # input:
  # gene_results= timelm results of an specific gene
  # output:
  # plot object
  
  gene_results$PVAL = -log(gene_results$PVAL)
  
  vp = ggplot(data=gene_results, aes(x=COEF, y=PVAL, label=TIME),guide=FALSE) + geom_point(colour="white", fill="lightblue", shape=21) +
    geom_text(size=5,colour="black") + scale_size_continuous(range = c(8, 18)) + scale_x_continuous(name="Coefficient") +
    scale_y_continuous(name="-log(pval)") + geom_hline(yintercept=-log(TRSH), linetype="dashed") + ggtitle(paste(GENE,"LINEAR MODEL")) + 
    theme(plot.margin = unit(c(1,3.3,1,1), "cm"))
  
  return(vp)
  
}

plotHMps = function(genes,expmat,targetdf,time.point){
  #
  # Plots time-specific heatmaps of a count matrix
  #
  # inputs:
  # -genes: a list of genes to be clustered and plotted
  # -expmat: a count matrix that contains the subset of data wanted to be plotted
  # -
  color = ifelse(targetdf$outcome == 0,"lightcoral","lightblue")
  targetdf = cbind(targetdf,color)
  
  # creates a own color palette from red to green
  my_palette <- colorRampPalette(c("darkred", "lightgrey", "darkblue"))(n = 299)
  
  # Divide the matrix in time points
  if(is.null(time.point)){
    for(t in sort(unique(targetdf$time))){
      ttargets = filter(targetdf,time==t)
      tmat = expmat[genes,as.character(ttargets$sample_mask_id)]
      
      heatmap.2(tmat,
                main=paste("Time ",as.character(t*-1),sep=""),
                Rowv = TRUE,
                Colv = TRUE,
                distfun = dist,
                hclustfun = hclust,
                trace = "none",
                sepwidth=c(0.0005,0.0005),
                dendrogram = "column",
                col=my_palette,       # use on color palette defined earlier
                ColSideColors = as.character(ttargets$color)
      )
      par(xpd=TRUE)
      legend("bottomleft",      # location of the legend on the heatmap plot
             legend = c("control","case"), # category labels
             col = c("lightcoral","lightblue"),  # color key
             lty= 1,             # line style
             lwd = 10,            # line width
             inset=c(-0.1,0)
      )
      
    }
  }else{
    for(t in time.point){
      #t = time.point
      ttargets = filter(targetdf,time==t)
      tmat = expmat[genes,as.character(ttargets$sample_mask_id)]
      
      heatmap.2(tmat,
                main=paste("Time ",as.character(t*-1),sep=""),
                Rowv = TRUE,
                Colv = TRUE,
                distfun = dist,
                hclustfun = hclust,
                trace = "none",
                dendrogram = "column",
                col=my_palette,       # use on color palette defined earlier
                ColSideColors = as.character(ttargets$color)
      )
      par(xpd=TRUE)
      legend("bottomleft",      # location of the legend on the heatmap plot
             legend = c("Control","Case"), # category labels
             col = c("lightcoral","lightblue"),  # color key
             lty= 1,             # line style
             lwd = 10,           # line width
             inset=c(-0.1,0)
      )
      
    }
  }
}

#
# LIMMA and gene expression pipelines
# (MasigPro processing)
#

runTEDDYlimma = function(signatureTargets, sigEXPMAT){
  
  geneList = list()
  
  Sclass = paste(ifelse(signatureTargets$outcome==0,"Healthy","Sick"),signatureTargets$time,sep=".")
  signatureTargets = cbind(signatureTargets,Sclass)
  
  f = factor(signatureTargets$Sclass,levels = sort(unique(signatureTargets$Sclass)))
  design = model.matrix(~0+f)
  colnames(design) =  levels(f)
  
  fit <- lmFit(sigEXPMAT, design)
  cont.IA = makeContrasts(
    Dif03 = (Sick.3 - Sick.0) - (Healthy.3 - Healthy.0),
    Dif36 = (Sick.6 - Sick.3) - (Healthy.6 - Healthy.3),
    Dif69 = (Sick.9 - Sick.6) - (Healthy.9 - Healthy.6),
    Dif129 = (Sick.12 - Sick.9) - (Healthy.12 - Healthy.9),
    levels = design
  )
  fit2 = contrasts.fit(fit,cont.IA)
  fit2 = eBayes(fit2)
  
  lmres = topTableF(fit2,adjust.method = "BH",number = Inf)
  geneList[["global"]] = lmres
  
  for(i in 1:4){
    pr = paste("Period", as.character(i * -1))
    lmres = topTable(fit2,coef = i,adjust.method = "BH",number = Inf)
    geneList[[pr]] = lmres
  }
  
  return(geneList)
}

runTimeCourse_test = function(countMAT, annotationDF, compCOL){
  # It only supports 1-to-1 comparisons
  # A - B is always compared
  
  geneList = list()
  
  Gs = unique(annotationDF[,compCOL])
  ClassA = Gs[1]
  ClassB = Gs[2]
  
  Sclass = paste(ifelse(annotationDF[,compCOL]==ClassA,"ClassA","ClassB"),annotationDF$time,sep=".")
  annotationDF = cbind(annotationDF,Sclass)
  
  f = factor(annotationDF$Sclass,levels = sort(unique(annotationDF$Sclass)))
  design = model.matrix(~0+f)
  colnames(design) =  levels(f)
  
  countMAT = countMAT[,as.character(annotationDF$sample_mask_id)]
  
  fit <- lmFit(countMAT, design)
  cont.IA = makeContrasts(
    Dif03 = (ClassA.3 - ClassA.0) - (ClassB.3 - ClassB.0),
    Dif36 = (ClassA.6 - ClassA.3) - (ClassB.6 - ClassB.3),
    Dif69 = (ClassA.9 - ClassA.6) - (ClassB.9 - ClassB.6),
    Dif129 = (ClassA.12 - ClassA.9) - (ClassB.12 - ClassB.9),
    levels = design
  )
  fit2 = contrasts.fit(fit,cont.IA)
  fit2 = eBayes(fit2)
  
  lmres = topTableF(fit2,adjust.method = "BH",number = Inf)
  
  geneList[["global"]] = lmres
  
  for(i in 1:4){
    pr = paste("Period", as.character(i * -1))
    lmres = topTable(fit2,coef = i,adjust.method = "BH",number = Inf)
    geneList[[pr]] = lmres
  }
  
  geneList[["ClassA"]] = ClassA
  geneList[["ClassB"]] = ClassB
  
  return(geneList)
}

runTEDDYlimma_basic = function(stargets, sEXPMAT){
  #
  # Functions fit a simple Sick vs Healthy contrast
  #
  
  geneList = list()
  
  Sclass = ifelse(stargets$outcome==0,"Healthy","Sick")
  stargets = cbind(stargets,Sclass)
  
  f = factor(stargets$Sclass,levels = sort(unique(stargets$Sclass)))
  design = model.matrix(~0+f)
  colnames(design) =  levels(f)
  
  fit <- lmFit(sEXPMAT, design)
  cont.IA = makeContrasts(
    Dif = Sick - Healthy,
    levels = design
  )
  fit2 = contrasts.fit(fit,cont.IA)
  fit2 = eBayes(fit2)
  
  lmres = topTable(fit2,coef = 1,adjust.method = "BH",number = Inf)
  
  return(lmres)
}


ExhaustiveLimma = function(stargets,sEXPMAT){
  #
  # Functions to obtain time specific gene expression changes
  #
  
  resultList = list()
  
  for(i in unique(stargets$time)){
    
    tstargets = filter(stargets,time==i)
    if(nrow(tstargets)>0){
      
      # At least 2 samples per outcome
      
      ncase = sum(tstargets$outcome==1)
      ncontrol = sum(tstargets$outcome==0)
      
      if((ncase>1) & (ncontrol>1)){
        
        tsEXPMAT = sEXPMAT[,as.character(tstargets$sample_mask_id)]
        resultList[[paste("Time",as.character(i),sep="_")]] = runTEDDYlimma_basic(stargets = tstargets,sEXPMAT = tsEXPMAT)
        
      }
    }
  }
  
  return(resultList)
  
}

#####################
#
#
#

filterLimmaGLSv1 = function(GLSlist,TRSH){
  
  geneList = list()
  
  for(contName in names(GLSlist)){
    GLS = GLSlist[[contName]]
    Tgenes = rownames(GLS[GLS$adj.P.Val<=TRSH,])
    geneList[[contName]] = Tgenes
  }
  
  return(geneList)
}

filterLimmaGLSv3_list = function(GLSlist,TRSH){
  
  geneList = list()
  
  for(contName in names(GLSlist)){
    GLS = GLSlist[[contName]]
    Tgenes = filterLimmaGLSv3(GLS = GLS,pTRSH = TRSH,FC_col = "logFC")
    geneList[[contName]] = Tgenes
  }
  
  return(geneList)
}


filterLimmaGLS = function(GLS,pTRSH=0.05,fcTRSH="adj",global=FALSE){
  
  if(global==TRUE){
    Tgenes = rownames(GLS[GLS$adj.P.Val<=pTRSH,])
  }else{
    
    if(fcTRSH == "adj"){
      lfc = GLS$logFC
      fcTRSHup = mean(lfc) + (2 * sd(lfc))
      fcTRSHdown = mean(lfc) + (2 * sd(lfc))
    }else{
      fcTRSHup = fcTRSH
      fcTRSHdown = fcTRSH * -1
    }
    
    plog = GLS$adj.P.Val<=pTRSH
    fclog = GLS$logFC<=fcTRSHdown | GLS$logFC>=fcTRSHup
    Tgenes = rownames(GLS[plog & fclog,])
  }
  
  return(Tgenes) 
}

filterLimmaGLSv3 = function(GLS,pTRSH=0.05,FC_col = "logFC"){
  
  filt_GLS = GLS[GLS$adj.P.Val<=pTRSH,]
  up_genes = rownames(filt_GLS[filt_GLS[,FC_col]>0,])
  down_genes = rownames(filt_GLS[filt_GLS[,FC_col]<0,])
  
  return(list("up"=up_genes,"down"=down_genes))
  
}

#Filter general linear model
filterGLM = function(GLMres,THRS){
  
  SIG = GLMres[GLMres$adj_P_val<=THRS,]
  UP = SIG$Hallmark[SIG$Coefficient>0]
  DOWN = SIG$Hallmark[SIG$Coefficient<0]
  
  return(list("up"=UP,"down"=DOWN))
}

writeSummary = function(Rdf){
  #
  # Writes a summary table from a DF which contains a list of features in the first column and a class of sample in the second
  #
  Rdf = na.omit(Rdf)
  smat = matrix(0,nrow = length(unique(Rdf[,1])), ncol = length(unique(Rdf[,2])))
  rownames(smat) = (unique(Rdf[,1]))
  colnames(smat) = (unique(Rdf[,2]))
  
  for(i in 1:nrow(Rdf)){
    smat[Rdf[i,1],Rdf[i,2]] = 1
  }
  
  smat = smat[order(rownames(smat)),]
  smat = smat[,order(colnames(smat))]
  
  return(smat)
  
}

writeDEX_list = function(DEX_Rlist,TRSH=0.05,fcol = "adj.P.Val",f=""){
  #
  # Writes the results from a set of different time-specific experiments
  # input:
  # DEX_Rlist = a list of objects generated from Exhaustive limma or any list of lists of gene level statistics
  # TRSH = pval threshold
  # fcol = column to evaluate TRSH
  #
  fname= paste("mkdir",f,sep=" ")
  system(fname)
  gene_list = c()
  kegg_list = c()
  go_list = c()
  
  for(aab_ag in names(DEX_Rlist)){
    print(aab_ag)
    sDEXlist = DEX_Rlist[[aab_ag]]
    fname = paste(f,"/",aab_ag,sep="")
    system(paste("mkdir",fname,sep=" "))
    
    if(length(sDEXlist)>0){
      
      ESfname= paste(fname,"ES", sep = "/")
      DEXfname= paste(fname,"DEX",sep = "/")
      system(paste("mkdir",ESfname,sep=" "))
      system(paste("mkdir",DEXfname,sep=" "))
      
      for(t in names(sDEXlist)){
        print(t)
        
        tag = paste(aab_ag,t,sep="_")
        filename = paste(DEXfname,"/",tag,".txt",sep="")
        goESfilename = paste(ESfname,"/",tag,"_goES",".txt",sep="")
        keggESfilename = paste(ESfname,"/",tag,"_keggES",".txt",sep="")
        
        tsDEXlist = sDEXlist[[t]]
        
        tsDEXlist = tsDEXlist[tsDEXlist[,fcol] <= TRSH,]
        if(nrow(tsDEXlist)>0){
          write.table(tsDEXlist,file=filename,quote = F,row.names = T,col.names = T,sep = "\t")
          genes = cbind(rownames(tsDEXlist),tag)
          gene_list = rbind(gene_list,genes)
          
          probe_symbol = na.omit(data.frame(mapIds(illuminaHumanv4.db, keys=genes, keytype = "SYMBOL", column="ENTREZID",multiVals = "first"),stringsAsFactors = F))
          keggES_R = kegga(probe_symbol[,1])
          keggES_R = keggES_R[keggES_R$P.DE<=TRSH,]
          if(nrow(keggES_R)>0){
            write.table(keggES_R,file = keggESfilename,quote = F,sep = "\t",row.names = F,col.names = F)
            terms = cbind(keggES_R[,1],tag)
            kegg_list = rbind(kegg_list,terms)
          }
          goES_R = goana(probe_symbol[,1])
          goES_R = goES_R[goES_R$P.DE<=TRSH,]
          if(nrow(goES_R)>0){
            write.table(goES_R,file = goESfilename,quote = F,sep = "\t",row.names = F,col.names = F)
            terms = cbind(goES_R[,1],tag)
            go_list = rbind(go_list,terms)
          }
        }
      }
    }
  }
  
  sumname = paste(f,"/","gene_summary",".txt",sep="")
  goESsumname = paste(f,"/","goES_summary",".txt",sep="")
  keggESsumname = paste(f,"/","keggES_summary",".txt",sep="")
  
  gene_summary = writeSummary(Rdf = gene_list)
  write.table(gene_summary,file = sumname,quote = F,sep = "\t",row.names = T,col.names = T)
  
  goES_summary = writeSummary(Rdf = go_list)
  write.table(goES_summary,file = goESsumname,quote = F,sep = "\t",row.names = T,col.names = T)
  
  keggES_summary = writeSummary(Rdf = kegg_list)
  write.table(keggES_summary,file = keggESsumname,quote = F,sep = "\t",row.names = T,col.names = T)
  
  return(list("genes"=gene_list,"kegg"=kegg_list,"go"=go_list))
  
}

getConservedSig = function(summobj, gdivision){
  ConservedSig = list()
  for(g in gdivision){
    selcols = grep(g,colnames(summobj))
    selmat = summobj[,selcols]
    cDF = data.frame(sort(rowSums(selmat),decreasing = T))
    colnames(cDF) = c("counts")
    cDF = cDF[cDF$counts>0,,drop=F]
    N = ncol(selmat)
    cDF = cbind(cDF,N)
    ConservedSig[[g]] = cDF
  }
  return(ConservedSig)
}

#
# Enrichment
#

runGESlimma = function(signatureTargets,sigEXPMAT,EnrichData,TRSH){
  #First, transform symbols into ENTREZ
  probe_names = rownames(sigEXPMAT)
  probe_symbol = na.omit(data.frame(mapIds(illuminaHumanv4.db, keys=probe_names, keytype = "SYMBOL", column="ENTREZID",multiVals = "first"),stringsAsFactors = F))
  #FIlter exprmat with known probes
  sigEXPMAT = sigEXPMAT[rownames(probe_symbol),]
  rownames(sigEXPMAT) = probe_symbol[,1]
  sigEXPMAT = avereps(sigEXPMAT,ID=rownames(sigEXPMAT))
  
  geneList = list()
  
  Sclass = paste(ifelse(signatureTargets$outcome==0,"Healthy","Sick"),signatureTargets$time,sep=".")
  signatureTargets = cbind(signatureTargets,Sclass)
  
  f = factor(signatureTargets$Sclass,levels = sort(unique(signatureTargets$Sclass)))
  design = model.matrix(~0+f)
  colnames(design) =  levels(f)
  
  fit <- lmFit(sigEXPMAT, design)
  cont.IA = makeContrasts(
    Dif03 = (Sick.3 - Sick.0) - (Healthy.3 - Healthy.0),
    Dif36 = (Sick.6 - Sick.3) - (Healthy.6 - Healthy.3),
    Dif69 = (Sick.9 - Sick.6) - (Healthy.9 - Healthy.6),
    Dif129 = (Sick.12 - Sick.9) - (Healthy.12 - Healthy.9),
    levels = design
  )
  fit2 = contrasts.fit(fit,cont.IA)
  fit2 = eBayes(fit2)
  
  lmres = topTableF(fit2,adjust.method = "BH",number = Inf)
  Tgenes = rownames(lmres[lmres$adj.P.Val<=TRSH,])
  
  if(EnrichData=="GO"){
    GSEAres = goana(Tgenes,species = "Hs")
    geneList[["global"]] = GSEAres[GSEAres$P.DE<0.05,]
    for(i in 1:4){
      pr = paste("Period", as.character(i * -1))
      GSEAres = goana(fit2, coef = i,FDR = 0.05)
      geneList[[pr]] = filter(GSEAres,P.Up<0.05 | P.Down<0.05)
    }
  }else if(EnrichData=="KEGG"){
    GSEAres = kegga(Tgenes,species = "Hs")
    geneList[["global"]] = GSEAres[GSEAres$P.DE<0.05,]
    for(i in 1:4){
      pr = paste("Period", as.character(i * -1))
      GSEAres = kegga(fit2, coef = i,FDR = 0.05)
      geneList[[pr]] = filter(GSEAres,P.Up<0.05 | P.Down<0.05)
    }
  }
  
  return(geneList)
}


#####################
#
#

MasigPro_targets = function(TARGETS){
  
  # First, for a given target table, drop NAs from FirstAAb
  TARGETS = TARGETS[!is.na(TARGETS$FirstAAb),]
  
  AAblist = TARGETS$FirstAAb
  AAblist[TARGETS$outcome==0] = 0
  
  # Determine classes
  Sclass = as.numeric(factor(paste(ifelse(TARGETS$outcome==0,"Healthy","Sick"),TARGETS$time,sep=".")))
  
  # Create column names of AAb 
  AAbclass = paste(as.character(AAblist),"AAb",sep="")
  AAbcols = sort(unique(AAbclass))
  AAbmatrix = c()
  for(AAb in AAbcols){
    AAbmatrix = cbind(AAbmatrix,ifelse(AAbclass==AAb,1,0))
  }
  AAbcols[AAbcols=="0AAb"] = "Control"
  
  # Creating target file
  mspTargets = cbind(TARGETS$time,Sclass,AAbmatrix)
  colnames(mspTargets) = c("Time","Replicate",AAbcols)
  rownames(mspTargets) = TARGETS$sample_mask_id 
  
  return(mspTargets)
}

masigPro_edesign = function(target_mat,factor_column,factor_name,time_column,sample_column){
  
  group_cols = sort(unique(target_mat[,factor_column]))
  group_names = paste(factor_name,group_cols,sep="_")
  
  group_mat = matrix(0,nrow = NROW(target_mat),ncol = length(group_cols))
  edesign = cbind(target_mat[,c(sample_column,time_column)],paste(target_mat[,c(time_column)],target_mat[,c(factor_column)],sep="_"),group_mat)
  
  # Format edesign
  rownames(edesign) = edesign[,1]
  edesign = edesign[,-1]
  colnames(edesign) = c("time","replicates",group_names)
  edesign[,"replicates"] = as.numeric(factor(edesign[,"replicates"]))
  
  for(i in 1:length(group_cols)){
    
    group_val = group_cols[i]
    group_val_name = group_names[i]
    selected_samples = as.character((target_mat[,sample_column])[target_mat[,factor_column] == group_val])
    edesign[selected_samples,group_val_name] = 1
    
  }
  
  return(edesign)
}

#####################
# DE genes searching tool
# Note, this only works in 4 classifications: ALL, AAb.1,AAb.2,AAb.5
#####################

getDEgenes = function(AAb, t_period, pTRSH, fcTRSH = "adj"){
  if(AAb %in% names(AAbGLSResults)){
    AAb_GLS = AAbGLSResults[[AAb]]
    if(t_period %in% names(AAb_GLS)){
      GLS = AAb_GLS[[t_period]]
      if(t_period == "global"){
        T_genes = filterLimmaGLS(GLS = GLS,pTRSH=pTRSH,fcTRSH=fcTRSH,global=TRUE)
      }else{
        T_genes = filterLimmaGLS(GLS = GLS,pTRSH=pTRSH,fcTRSH=fcTRSH,global=FALSE)
      }
      return(T_genes)
    }
  }
}

####
# Obtain differentially expressed genes between AAb
####
# Nomenclature = AAbxvsAAby.pz
# where x<y E c(1,2,5)
# and
# z E c(1,2,3,4) which represents comparisons time specific
#
# or
# type global

getDEgenesAAbcomp = function(comp, pTRSH, fcTRSH = "adj"){
  if(comp %in% names(GLS_AAbcomparison)){
    GLS = GLS_AAbcomparison[[comp]]
    if(comp == "global"){
      T_genes = filterLimmaGLS(GLS = GLS,pTRSH=pTRSH,fcTRSH=fcTRSH,global=TRUE)
    }else{
      T_genes = filterLimmaGLS(GLS = GLS,pTRSH=pTRSH,fcTRSH=fcTRSH,global=FALSE)
    }
    return(T_genes)
  }
}


#
# GSEA Enrichment tools
# Originally created by the Broad Institute
#

GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {  
  #
  # Computes the weighted GSEA score of gene.set in gene.list. 
  # The weighted score type is the exponent of the correlation 
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
  # necessary to input the correlation vector with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
  #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
  #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1) 
  #   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain") 
  #   RES: Numerical vector containing the running enrichment score for all locations in the gene list 
  #   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list 
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  
  tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
  no.tag.indicator <- 1 - tag.indicator 
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 
  if (weighted.score.type == 0) {
    correl.vector <- rep(1, N)
  }
  alpha <- weighted.score.type
  correl.vector <- abs(correl.vector**alpha)
  sum.correl.tag  <- sum(correl.vector[tag.indicator == 1])
  norm.tag    <- 1.0/sum.correl.tag
  norm.no.tag <- 1.0/Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)      
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (max.ES > - min.ES) {
    #      ES <- max.ES
    ES <- signif(max.ES, digits = 5)
    arg.ES <- which.max(RES)
  } else {
    #      ES <- min.ES
    ES <- signif(min.ES, digits=5)
    arg.ES <- which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))    
}


OLD.GSEA.EnrichmentScore <- function(gene.list, gene.set) {  
  #
  # Computes the original GSEA score from Mootha et al 2003 of gene.set in gene.list 
  #
  # Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1) 
  #   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain") 
  #   RES: Numerical vector containing the running enrichment score for all locations in the gene list 
  #   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list 
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  
  
  tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
  no.tag.indicator <- 1 - tag.indicator 
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 
  
  norm.tag    <- sqrt((N - Nh)/Nh)
  norm.no.tag <- sqrt(Nh/(N - Nh))
  
  RES <- cumsum(tag.indicator * norm.tag - no.tag.indicator * norm.no.tag)      
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (max.ES > - min.ES) {
    ES <- signif(max.ES, digits=5)
    arg.ES <- which.max(RES)
  } else {
    ES <- signif(min.ES, digits=5)
    arg.ES <- which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))    
}

GSEA.EnrichmentScore2 <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {  
  #
  # Computes the weighted GSEA score of gene.set in gene.list. It is the same calculation as in 
  # GSEA.EnrichmentScore but faster (x8) without producing the RES, arg.RES and tag.indicator outputs.
  # This call is intended to be used to asses the enrichment of random permutations rather than the 
  # observed one.
  # The weighted score type is the exponent of the correlation 
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
  # necessary to input the correlation vector with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
  #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
  #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1) 
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 
  
  loc.vector <- vector(length=N, mode="numeric")
  peak.res.vector <- vector(length=Nh, mode="numeric")
  valley.res.vector <- vector(length=Nh, mode="numeric")
  tag.correl.vector <- vector(length=Nh, mode="numeric")
  tag.diff.vector <- vector(length=Nh, mode="numeric")
  tag.loc.vector <- vector(length=Nh, mode="numeric")
  
  loc.vector[gene.list] <- seq(1, N)
  tag.loc.vector <- loc.vector[gene.set]
  
  tag.loc.vector <- sort(tag.loc.vector, decreasing = F)
  
  if (weighted.score.type == 0) {
    tag.correl.vector <- rep(1, Nh)
  } else if (weighted.score.type == 1) {
    tag.correl.vector <- correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } else if (weighted.score.type == 2) {
    tag.correl.vector <- correl.vector[tag.loc.vector]*correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } else {
    tag.correl.vector <- correl.vector[tag.loc.vector]**weighted.score.type
    tag.correl.vector <- abs(tag.correl.vector)
  }
  
  norm.tag <- 1.0/sum(tag.correl.vector)
  tag.correl.vector <- tag.correl.vector * norm.tag
  norm.no.tag <- 1.0/Nm
  tag.diff.vector[1] <- (tag.loc.vector[1] - 1) 
  tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
  tag.diff.vector <- tag.diff.vector * norm.no.tag
  peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
  valley.res.vector <- peak.res.vector - tag.correl.vector
  max.ES <- max(peak.res.vector)
  min.ES <- min(valley.res.vector)
  ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)
  
  return(list(ES = ES))
  
}

GSEA_ES = function(eMat, gene.set){
  # GSEA.EnrichmentScore Applied to a matrix
  # Input: 
  # eMat = matrix of expression profiles
  # gene.set =  vector of genes to  be tested
  # Output:
  # -vector with enrichment scores per sample 
  #
  ES_list = c()
  for(i in 1:ncol(eMat)){
    gene.list = names(sort(eMat[,i],decreasing = T))
    ES_list = c(ES_list,GSEA.EnrichmentScore(gene.list = gene.list,gene.set = gene.set,weighted.score.type = 0)$ES)
  }
  names(ES_list) = colnames(eMat)
  return(ES_list)
}

GSEA_ES_v2 = function(eMat, eTargets,gene.set){
  # GSEA.EnrichmentScore Applied to a matrix
  # Input: 
  # eMat = matrix of expression profiles
  # gene.set =  vector of genes to  be tested
  # eTargets = targets of the sample analysed
  # Output:
  # -vector with enrichment scores per sample 
  #
  ES_list = c()
  O = GSEA.GeneRanking(A=eMat, class.labels=eTargets$outcome, gene.labels=rownames(eMat), nperm=2, permutation.type = 0, sigma.correction = "GeneCluster", fraction=1.0, replace=F, reverse.sign= F)
  obs.correl.matrix = O$obs.s2n.matrix
  obs.s2n = apply(obs.correl.matrix, 1, median)  # using median to assign enrichment scores
  obs.s2n = sort(obs.s2n, decreasing=T)   
  
  for(i in 1:ncol(eMat)){
    gene.list = names(sort(eMat[,i],decreasing = T))
    signal.noise.scores = obs.s2n[gene.list]
    ES_list = c(ES_list,GSEA.EnrichmentScore(gene.list = gene.list,gene.set = gene.set,weighted.score.type = 1,correl.vector = signal.noise.scores)$ES)
  }
  names(ES_list) = colnames(eMat)
  return(ES_list)
}

GSEA_ES_v3 = function(eMat, eTargets, gene.set, weighted = TRUE){
  # For a given set of samples (expression matrix) calculates the enrichment score of a gene.set
  # If you decide to weight the list, Signal2Noise is calculated. (Broad implementation: difference of means scaled by the standard deviation)
  #
  ES_list = c()
  LeadingEdge_list = list()
  
  if(weighted == FALSE){
    
    for(i in 1:ncol(eMat)){
      sampleName = colnames(eMat)[i]
      gene.list = names(sort(eMat[,i],decreasing = T))
      #ADD NULL robustness
      int.gene.set = intersect(gene.list,gene.set)
      if(length(int.gene.set>0)){
        ES_results = GSEA.EnrichmentScore(gene.list = gene.list,gene.set = int.gene.set,weighted.score.type = 0)
        ES = ES_results$ES
        ES_list = c(ES_list,ES)
        peak = ES_results$arg.ES
        if(ES>0){
          coordinates = as.logical(ES_results$indicator[1:peak])
          leadingedge = gene.list[1:peak]
          leadingedge =  leadingedge[coordinates]
          LeadingEdge_list[[sampleName]] = leadingedge
        }else{
          coordinates = as.logical(ES_results$indicator[peak:length(gene.list)])
          leadingedge = gene.list[peak:length(gene.list)]
          leadingedge =  leadingedge[coordinates]
          LeadingEdge_list[[sampleName]] = leadingedge
        }
      } else{
        ES_list = c(ES_list,NA)
        LeadingEdge_list[[sampleName]] = NA
      }
      
    }
    names(ES_list) = colnames(eMat)
    return(list(ES_vector = ES_list, LeadingEdge = LeadingEdge_list))
    
    
  }
  
  if(weighted == TRUE){
    
    O = GSEA.GeneRanking(A=eMat, class.labels=eTargets$outcome, gene.labels=rownames(eMat), nperm=2, permutation.type = 0, sigma.correction = "GeneCluster", fraction=1.0, replace=F, reverse.sign= F)
    obs.correl.matrix = O$obs.s2n.matrix
    obs.s2n = apply(obs.correl.matrix, 1, median)  # using median to assign enrichment scores
    obs.s2n = sort(obs.s2n, decreasing=T)   
    
    for(i in 1:ncol(eMat)){
      sampleName = colnames(eMat)[i]
      gene.list = names(sort(eMat[,i],decreasing = T))
      signal.noise.scores = obs.s2n[gene.list]
      int.gene.set = intersect(gene.list,gene.set)
      ES_results = GSEA.EnrichmentScore(gene.list = gene.list,gene.set = int.gene.set,weighted.score.type = 1,correl.vector = signal.noise.scores)
      ES = ES_results$ES
      ES_list = c(ES_list,ES)
      peak = ES_results$arg.ES
      if(ES>0){
        coordinates = as.logical(ES_results$indicator[1:peak])
        leadingedge = gene.list[1:peak]
        leadingedge =  leadingedge[coordinates]
        LeadingEdge_list[[sampleName]] = leadingedge
      }else{
        coordinates = as.logical(ES_results$indicator[peak:length(gene.list)])
        leadingedge = gene.list[peak:length(gene.list)]
        leadingedge =  leadingedge[coordinates]
        LeadingEdge_list[[sampleName]] = leadingedge
      }
    }
    names(ES_list) = colnames(eMat)
    return(list(ES_vector = ES_list, LeadingEdge = LeadingEdge_list))
    
  }
  
}

plotES = function(ES_list,sigTargets,TRSH=0.05,main_title){
  # This function compares the Enrichment Scores of patients and controls, with a T-test
  # Input:
  # -ES_list: Results from GSEA_ES
  # -sigTargets =  Annotation of the matrix analised in GSEA_ES
  # -TRSH =  Significance threshold of the T-test
  # -main_title = plot Title
  # Output:
  # -Plot of ES dynamics
  # -Boxplot of groups across time
  plotdf = data.frame(cbind(ES_list,sigTargets$time,ifelse(sigTargets$outcome==0,"Healthy","Sick")),stringsAsFactors = F)
  colnames(plotdf) = c("ES","TIME","OUTCOME")
  plotdf$ES = as.numeric(plotdf$ES) 
  plotdf$TIME = as.numeric(plotdf$TIME) * -1
  
  new_plotdf = c()
  
  for(t in unique(plotdf$TIME)){
    
    tsp_df = filter(plotdf,TIME==t)
    tsp_df$OUTCOME = as.factor(tsp_df$OUTCOME)
    ttestRes = t.test(tsp_df$ES~tsp_df$OUTCOME)
    if(ttestRes$p.value<=TRSH){
      pval = "significant"
      tsp_df = cbind(tsp_df,pval)
    }else{
      pval = "not-significant"
      tsp_df = cbind(tsp_df,pval)
    }
    
    new_plotdf = rbind(new_plotdf,tsp_df)
  }
  
  new_plotdf$pval = factor(new_plotdf$pval,levels = c("significant","not-significant"))
  new_plotdf$TIME = as.numeric(new_plotdf$TIME)
  
  ESplot = ggplot(new_plotdf, aes(x=TIME, y=ES, color=OUTCOME, shape = pval)) + geom_point(size=2.5) + scale_x_continuous(name="Time",breaks=c(-12,-9,-6,-3,0)) + ggtitle(main_title) + theme(axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12))
  print(ESplot)
  
  new_plotdf$TIME = as.factor(new_plotdf$TIME)
  gp = ggplot(new_plotdf, aes(x = TIME, y = ES)) + geom_boxplot(aes(fill = OUTCOME,linetype = pval), alpha = 0.5) + ylab("Enrichment Score") + theme(axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12), plot.margin = unit(c(1,0.6,1,0.6), "cm")) + ggtitle(main_title)
  print(gp)
  
}

getFunctionMatrix = function(Emat, Signature_Dict){
  Fmat = c()
  for(geneSet_names in names(Signature_Dict)){
    print(geneSet_names)
    geneSet = Signature_Dict[[geneSet_names]]
    ES_vector = GSEA_ES(eMat = Emat,gene.set = geneSet)
    Fmat = rbind(Fmat, ES_vector)
  }
  rownames(Fmat) = names(Signature_Dict)
  return(Fmat)
}

getFunctionMatrix_v2 = function(eMat, eTargets, Signature_Dict, weighted=TRUE){
  Fmat = c()
  Ledge = list()
  for(geneSet_names in names(Signature_Dict)){
    print(geneSet_names)
    geneSet = Signature_Dict[[geneSet_names]]
    ES_results = GSEA_ES_v3(eMat = eMat,gene.set = geneSet,eTargets = eTargets,weighted = weighted)
    ES_vector = ES_results$ES_vector
    Fmat = rbind(Fmat, ES_vector)
    Ledge[[geneSet_names]] = ES_results$LeadingEdge
  }
  rownames(Fmat) = names(Signature_Dict)
  return(list(ESmat = Fmat, LeadingEdge = Ledge))
}

######################################
#
# TOOLS for paintomics
#

getGroupMedian = function(groupcolumn, stargets, sexprmat){
  #Inputs:
  #groupcolumn: attribute to be used for the separation
  #stargets: set of samples to be reduced
  #sexprmat: expression matrix related to selected annotations
  meanEX = c()
  colLabs = c()
  for(f in sort(unique(stargets[,groupcolumn]))){
    colLabs = c(colLabs,f)
    fsamples = sexprmat[,as.character(stargets$sample_mask_id[stargets[,groupcolumn]==f])]
    meanEX = cbind(meanEX, rowMedians(fsamples))
  }
  colLabs = as.character(colLabs)
  rownames(meanEX) = rownames(sexprmat)
  colnames(meanEX) = c(paste(groupcolumn,colLabs,sep = "_"))
  return(meanEX)
}

getGroupMean = function(groupcolumn, stargets, sexprmat){
  #Inputs:
  #groupcolumn: attribute to be used for the separation
  #stargets: set of samples to be reduced
  #sexprmat: expression matrix related to selected annotations
  meanEX = c()
  colLabs = c()
  for(f in unique(stargets[,groupcolumn])){
    colLabs = c(colLabs,f)
    fsamples = sexprmat[,as.character(stargets$sample_mask_id[stargets[,groupcolumn]==f])]
    meanEX = cbind(meanEX, rowMeans(fsamples))
  }
  meanEX = cbind(rownames(meanEX),meanEX)
  colnames(meanEX) = c("feature_name",paste(groupcolumn,colLabs,sep = "_"))
  return(meanEX)
}

getGroupMeanMET = function(groupcolumn, stargets, sexprmat){
  #Inputs:
  #groupcolumn: attribute to be used for the separation
  #stargets: set of samples to be reduced
  #sexprmat: expression matrix related to selected annotations
  meanEX = c()
  colLabs = c()
  for(f in unique(stargets[,groupcolumn])){
    colLabs = c(colLabs,f)
    fsamples = sexprmat[,as.character(stargets$Sample.Mask.Id[stargets[,groupcolumn]==f])]
    #fsamples = sexprmat[,as.character(stargets$Sample.Mask.Id[stargets[,groupcolumn]==0])]
    #fsamples["id_fa__16_0__d3_3_22_258_25._neglip",]
    #table(is.infinite(fsamples))
    fsamples[is.infinite(fsamples)]<-0
    meanEX = cbind(meanEX, rowMeans(fsamples, na.rm = TRUE))
    head(meanEX)
  }
  meanEX = cbind(rownames(meanEX),meanEX)
  colnames(meanEX) = c("feature_name",paste(groupcolumn,colLabs,sep = "_"))
  return(meanEX)
}


##############################################
#
# Tools for signature comparison
#
##############################################

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

# Functions to plot intersection
similix = function(vectorA, vectorB, n = 14302, significancePercentage = .50){
  I = length(intersect(vectorA,vectorB))
  similix = I/n
  
  referenceV = min(c(length(vectorA)),length(vectorB))
  
  if(I >= round(referenceV*significancePercentage)){
    sig = 1
  }else{
    sig = 0
  }
  
  return(c(similix,sig))
}

jaccardIx = function(set1,set2){
  intersect_val = length(intersect(set1,set2))
  union_val = length(union(set1,set2))
  return(intersect_val/union_val)
}


getIntersect_matrix = function(SignatureDictionary, significancePercentage = .50){
  # Define matrices that are used for plotting (Directed, non-directed and their corresponding significance matrix)
  # NOTE: All NaN are converted to 0 and the similarity indexes are not rezised
  # Inputs: 
  # -SIgnatureDictionary: List with signatures to be compared
  # -significancePercentage: how much an intersection should cover from the smallest signature in a given comparison, to be considered significant
  # Outputs:
  # -list containing 4 matrices, 
  
  sigNames = names(SignatureDictionary)
  N = length(unique(unlist(SignatureDictionary)))
  
  sigmat = matrix(0,nrow = length(sigNames),ncol = length(sigNames))
  rownames(sigmat) = sigNames
  colnames(sigmat) = sigNames
  
  sigmat_sig = matrix(0,nrow = length(sigNames),ncol = length(sigNames))
  rownames(sigmat_sig) = sigNames
  colnames(sigmat_sig) = sigNames
  
  dir_sigmat = matrix(0,nrow = length(sigNames),ncol = length(sigNames))
  rownames(dir_sigmat) = sigNames
  colnames(dir_sigmat) = sigNames
  
  dir_sigmat_sig = matrix(0,nrow = length(sigNames),ncol = length(sigNames))
  rownames(dir_sigmat_sig) = sigNames
  colnames(dir_sigmat_sig) = sigNames
  
  # Fill the matrices
  
  for(i in 1:(length(sigNames)-1)){
    
    #Fill identities
    
    i_sig = SignatureDictionary[[sigNames[i]]]
    
    i_sig_up = i_sig$up
    i_sig_down = i_sig$down
    i_sig = unique(unlist(i_sig))
    
    simil_score = similix(i_sig, i_sig, n = N, significancePercentage = significancePercentage)
    
    sigmat[sigNames[i],sigNames[i]] = simil_score[1]
    sigmat_sig[sigNames[i],sigNames[i]] = simil_score[2]
    
    dir_sigmat[sigNames[i],sigNames[i]] = simil_score[1]
    dir_sigmat_sig[sigNames[i],sigNames[i]] = simil_score[2]
    
    for(j in (i+1):length(sigNames)){
      
      #Fill comparisons
      
      j_sig = SignatureDictionary[[sigNames[j]]]
      
      j_sig_up = j_sig$up
      j_sig_down = j_sig$down
      j_sig = unlist(j_sig)
      
      jix = similix(i_sig,j_sig, n = N, significancePercentage = significancePercentage)
      up_jix = similix(i_sig_up,j_sig_up, n = N, significancePercentage = significancePercentage)
      down_jix = similix(i_sig_down,j_sig_down, n = N, significancePercentage = significancePercentage)
      
      #Fill undirected
      sigmat[sigNames[i],sigNames[j]] = jix[1]
      sigmat[sigNames[j],sigNames[i]] = jix[1]
      sigmat_sig[sigNames[i],sigNames[j]] = jix[2]
      sigmat_sig[sigNames[j],sigNames[i]] = jix[2]
      
      #Fill directed
      dir_sigmat[sigNames[i],sigNames[j]] = up_jix[1]
      dir_sigmat_sig[sigNames[i],sigNames[j]] = up_jix[2]
      
      dir_sigmat[sigNames[j],sigNames[i]] = down_jix[1]
      dir_sigmat_sig[sigNames[j],sigNames[i]] = down_jix[2]
    }
  }
  
  i_sig = SignatureDictionary[[sigNames[j]]]
  
  i_sig_up = i_sig$up
  i_sig_down = i_sig$down
  i_sig = unique(unlist(i_sig))
  
  simil_score = similix(i_sig, i_sig, n = N, significancePercentage = significancePercentage)
  
  sigmat[sigNames[j],sigNames[j]] = simil_score[1]
  sigmat_sig[sigNames[j],sigNames[j]] = simil_score[2]
  
  dir_sigmat[sigNames[j],sigNames[j]] = simil_score[1]
  dir_sigmat_sig[sigNames[j],sigNames[j]] = simil_score[2]
  
  # Delete all NAN
  
  sigmat[is.nan(sigmat)] = 0
  sigmat_sig[is.nan(sigmat_sig)] = 0
  dir_sigmat[is.nan(dir_sigmat)] = 0
  dir_sigmat_sig[is.nan(dir_sigmat_sig)] = 0
  
  return(list("undirected" = sigmat,"undirected_significance" = sigmat_sig,"directed" = dir_sigmat,"directed_significance" = dir_sigmat_sig))
}


plotIntersection = function(Signature_mats){
  # Normalises similarity matrix and creates plots of the matrices, directed and not directed 
  # INPUT: 
  # - Signature_mats: getIntersect_matrix output
  #
  directed = Signature_mats$directed
  directed_sig = Signature_mats$directed_significance
  
  undirected = Signature_mats$undirected
  undirected_sig = Signature_mats$undirected_significance
  
  # First Plot directed mats
  
  up_regulated = get_upper_tri(directed)
  regulation="upregulated"
  up_df = cbind(melt(up_regulated,na.rm=TRUE),regulation)
  
  up_regulated_sig = get_upper_tri(directed_sig)
  up_df_sig = melt(up_regulated_sig,na.rm=TRUE)
  
  down_regulated = get_lower_tri(directed)
  regulation="downregulated"
  down_df = cbind(melt(down_regulated,na.rm=TRUE),regulation)
  
  down_regulated_sig = get_lower_tri(directed_sig)
  down_df_sig = melt(down_regulated_sig,na.rm=TRUE)
  
  significance = c(up_df_sig$value,down_df_sig$value)
  mat_df = rbind(up_df,down_df)
  mat_df = cbind(mat_df, significance)
  
  #Manipulate data
  colnames(mat_df) = c("SigA","SigB","Similarity_Ix","Regulation","Significance")
  mat_df$Regulation = as.character(mat_df$Regulation)
  mat_df$Regulation[mat_df$SigA == mat_df$SigB] = "identity"
  mat_df = unique(mat_df)
  
  mat_df$Similarity_Ix[mat_df$Similarity_Ix==0]=NA
  mat_df$Similarity_Ix[mat_df$Similarity_Ix==Inf]=NA
  mat_df$Similarity_Ix[is.nan(mat_df$Similarity_Ix)]=NA
  mat_df$Similarity_Ix = mat_df$Similarity_Ix/max(mat_df$Similarity_Ix,na.rm = T)
  
  mat_df$Regulation = as.factor(mat_df$Regulation)
  mat_df$Significance = factor(mat_df$Significance,levels=c(1,0))
  
  p = ggplot(data = mat_df, aes(x=SigA, y=SigB, fill=Regulation)) + geom_tile(colour="black") + geom_point(aes(size = Similarity_Ix,colour=Significance)) + scale_colour_manual(values=c("gold","black")) +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p)
  
  # Plot undirected mats
  
  up_regulated = get_upper_tri(undirected)
  up_df = melt(up_regulated,na.rm=TRUE)
  
  up_regulated_sig = get_upper_tri(undirected_sig)
  up_df_sig = melt(up_regulated_sig,na.rm=TRUE)
  significance = up_df_sig$value
  
  mat_df = cbind(up_df, significance)
  
  #Manipulate data
  colnames(mat_df) = c("SigA","SigB","Similarity_Ix","Significance")
  
  mat_df$Similarity_Ix[mat_df$Similarity_Ix==0]=NA
  mat_df$Similarity_Ix[mat_df$Similarity_Ix==Inf]=NA
  mat_df$Similarity_Ix[is.nan(mat_df$Similarity_Ix)]=NA
  mat_df$Similarity_Ix = mat_df$Similarity_Ix/max(mat_df$Similarity_Ix,na.rm = T)
  
  mat_df$Significance = factor(mat_df$Significance,levels=c(1,0))
  
  p = ggplot(data = mat_df, aes(x=SigA, y=SigB)) + geom_tile(colour="black",fill="lightgrey") + geom_point(aes(size = Similarity_Ix,colour=Significance)) + scale_colour_manual(values=c("red", "black")) +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(p)
  
}

#
# Plot similarity of selected tags
# 

SignatureIntersect = function(SignatureDictionary, collapse = FALSE, tags = c("Time"), significancePercentage = .50){
  # For a given subset of signatures, create directed and not directed similarity matrices and plot them.
  # You can create tags to create subgroups of signatures and collapse them, or simply just to reduce the set of signatures in a signature dictionary
  #
  # Inputs: 
  # -SignatureDictionary: A list of signatures
  # -collapse: TRUE if you like to create and evaluate TAG signatures
  # -tags: strings or regular expressions to subset the signatures. It is relative to names(SignatureDictionary)
  # -significancePercentage: how much an intersection should cover from the smallest signature in a given comparison, to be considered significant
  #
  # Outputs:
  # -list with similarity matrices and their respective significance matrices
  # -plot of the matrices
  #
  #
  #
  selectedSignatures = c()
  
  if(length(tags)<2){
    collapse == FALSE
    print("Collapse has been set to FALSE: Only 1 tag was provided")
  }
  
  if(collapse == TRUE){
    collapseDict = list()
  }
  
  for(tag in tags){
    selectedNames = names(SignatureDictionary)[grep(tag,names(SignatureDictionary))]
    selectedSignatures = c(selectedSignatures, selectedNames)
    
    if(collapse == TRUE){
      
      collapseDict[[tag]] = list("up"=c(),"down"=c())
      
      up_genes = c()
      down_genes = c()
      
      for(n in selectedNames){
        up_genes = c(up_genes,SignatureDictionary[[n]]$up)
        down_genes = c(up_genes,SignatureDictionary[[n]]$down)
      }
      
      up_genes = unique(up_genes)
      down_genes = unique(down_genes)
      
      #Not allow to have genes in both signatures
      #up_genes = setdiff(up_genes, down_genes)
      #down_genes = setdiff(down_genes, up_genes)
      
      collapseDict[[tag]]$up = up_genes
      collapseDict[[tag]]$down = down_genes
      
    }
  }
  
  if(collapse == TRUE){
    SignatureDictionary = collapseDict
  }else{
    SignatureDictionary = SignatureDictionary[selectedSignatures]
  }
  
  #Signatures matrix
  
  Signature_mats = getIntersect_matrix(SignatureDictionary = SignatureDictionary, significancePercentage = significancePercentage)
  
  #Plot Signatures
  plotIntersection(Signature_mats)
  
  return(Signature_mats)
  
}

################################################
#
#  Functions to functionally annotate GLS
#
################################################

runMultiGSA = function(GSC, GLS){
  pval = data.frame(GLS$adj.P.Val)
  rownames(pval) = rownames(GLS)
  Fchange = data.frame(GLS$logFC)
  rownames(Fchange) = rownames(GLS)
  tval = data.frame(GLS$t)
  rownames(tval) = rownames(GLS)
  
  gsaRes1 = runGSA(tval,gsc=GSC,geneSetStat = "mean",gsSizeLim=c(5,300))
  gsaRes2 = runGSA(tval,gsc=GSC,geneSetStat = "median",gsSizeLim=c(5,300))
  gsaRes3 = runGSA(tval,gsc=GSC,geneSetStat = "sum",gsSizeLim=c(5,300))
  gsaRes4 = runGSA(tval,gsc=GSC,geneSetStat = "maxmean",gsSizeLim=c(5,300))
  gsaRes5 = runGSA(pval,Fchange,gsc=GSC,geneSetStat = "fisher",gsSizeLim=c(5,300))
  gsaRes6 = runGSA(pval,Fchange,gsc=GSC,geneSetStat = "stouffer",gsSizeLim=c(5,300))
  gsaRes7 = runGSA(pval,Fchange,gsc=GSC,geneSetStat = "tailStrength",gsSizeLim=c(5,300))
  
  resList <- list(gsaRes1,gsaRes2,gsaRes3,gsaRes4,gsaRes5,gsaRes6,gsaRes7)
  names(resList) <- c("mean","median","sum","maxmean","fisher","stouffer","tailStrength")
  
  return(resList)
}

######################################
#
# Funtions to analyse Hallmarks
#
#####################################

fitGLM = function(cMAT,cTargets){
  # This function is used to fit a linear model at each time point for a given count matrix (it corrects with FDR)
  # Inputs:
  # cMAT=count matrix
  # cTargets=annotation table
  # Outputs:
  # - Table with results
  # - Matrix with the data use for the analysis
  
  ann_filtMAT = c()
  lmResults = c()
  rownames(cTargets) = cTargets$sample_mask_id
  
  for(t in c(0,3,6,9,12)){
    
    samples = as.character(filter(cTargets,time==t)$sample_mask_id)
    
    if(length(samples)>=4){ #Enough samples for model ?
      outcome = cTargets[samples,c("outcome","time")]
      filtMAT = data.frame(cbind(outcome,t(cMAT[,samples])))
      filtMAT$outcome = factor(filtMAT$outcome)
      
      if(sum(table(filtMAT$outcome) >= 2) == 2){ #Enough samples for model?
        
        for(bf in colnames(filtMAT)[3:(ncol(filtMAT))]){
          t_lm = summary(lm(filtMAT[,bf]~filtMAT[,"outcome"]), data = filtMAT)
          coef = coefficients(t_lm)[2,1]
          pval = coefficients(t_lm)[2,4]
          lmResults =  rbind(lmResults,c(t,bf,coef,pval))
        }
        
        ann_filtMAT = rbind(ann_filtMAT, filtMAT)
        
      }
      
    }
    
  }
  
  if(!is.null(lmResults)){
    
    lmResults = data.frame(lmResults,stringsAsFactors = F)
    colnames(lmResults) = c("Time","Hallmark","Coefficient","P_val")
    lmResults$Time = factor(lmResults$Time)
    lmResults$Coefficient = as.numeric(lmResults$Coefficient)
    lmResults$P_val = as.numeric(lmResults$P_val)
    adj_P_val = p.adjust(lmResults$P_val,method = "fdr")
    lmResults = cbind(lmResults,adj_P_val)
    
    return(list("lmResults" = lmResults,"lmData" = ann_filtMAT))
    
  }
  
}

fitGLM_v2 = function(cMAT,cTargets){
  # This function is used to fit a linear model at each time point for a given count matrix (it corrects with FDR)
  # Inputs:
  # cMAT=count matrix
  # cTargets=annotation table
  # Outputs:
  # - Table with results
  # - Matrix with the data use for the analysis
  
  lmResults = c()
  samples =  as.character(cTargets$sample_mask_id)
  rownames(cTargets) = samples
  outcome = cTargets[samples,c("outcome","time")]
  filtMAT = data.frame(cbind(outcome,t(cMAT[,samples])))
  filtMAT$outcome = factor(filtMAT$outcome)
  
  for(bf in colnames(filtMAT)[3:(ncol(filtMAT))]){
    cases_values = filtMAT[filtMAT$outcome == 1, bf]
    cases_values = sum(!is.na(cases_values))>1
    
    control_values = filtMAT[filtMAT$outcome == 0, bf]
    control_values = sum(!is.na(control_values))>1
    
    if(cases_values & control_values){
      
      naFREE_FMat = na.omit(filtMAT[,c(bf,"outcome")])
      
      t_lm = summary(lm(naFREE_FMat[,bf]~naFREE_FMat[,"outcome"]), data = naFREE_FMat)
      coef = coefficients(t_lm)[2,1]
      pval = coefficients(t_lm)[2,4]
      lmResults =  rbind(lmResults,c(bf,coef,pval))
      
    } 
  }
  
  lmResults = data.frame(lmResults,stringsAsFactors = F)
  colnames(lmResults) = c("Hallmark","Coefficient","P_val")
  lmResults$Coefficient = as.numeric(lmResults$Coefficient)
  lmResults$P_val = as.numeric(lmResults$P_val)
  adj_P_val = p.adjust(lmResults$P_val,method = "fdr")
  lmResults = cbind(lmResults,adj_P_val)
  
  return(lmResults)
  
}



plotLM_results = function(FitLM, FitLMdata){
  # Plot individual results from fitGLM: Results and data
  
  #Create Dictionary of significance
  
  for(sign_id in unique(FitLM$Hallmark)){
    
    significance_list = FitLM$Time[FitLM$Hallmark==sign_id]
    plot_df = FitLMdata[,c(1,2,grep(sign_id,colnames(FitLMdata)))]
    
    colnames(plot_df) = c("outcome","time","ES")
    plot_df$outcome = factor(plot_df$outcome)
    plot_df$time = factor(plot_df$time,levels=c(12,9,6,3,0))
    gp = ggplot(plot_df, aes(x = time, y = ES)) + geom_boxplot(aes(fill = outcome), alpha = 0.5) + ylab("Enrichment Score") + theme(axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12), plot.margin = unit(c(1,0.6,1,0.6), "cm")) + ggtitle(sign_id)
    
    gp = gp + annotate("text", x = factor(significance_list), y = max(plot_df$ES) + 0.03, label = "*", size=12, colour="blue")
    
    print(gp)
    
  }
}


#######################################

GSA_comparison = function(GSA_A,GSA_B,comp_type = "union"){
  #
  # Gets the union or intersection of the results observed in GSA analysis (Piano and LM)
  #
  #
  #
  comp_dic = list()
  
  if(comp_type=="union"){
    
    HMRKS = union(GSA_A$Hallmark,GSA_B$Hallmark)
    
  }else if(comp_type=="intersection"){
    
    HMRKS = intersect(GSA_A$Hallmark,GSA_B$Hallmark)
    
  }
  
  for(h in HMRKS){
    
    hA = filter(GSA_A,Hallmark==h)
    if(!is.null(hA)){
      th_A = hA$Time 
    }else{
      th_A = c()
    }
    
    hB = filter(GSA_B,Hallmark==h)
    if(!is.null(hB)){
      th_B = hB$Time 
    }else{
      th_B = c()
    }
    
    t_vector = union(th_A,th_B)
    
    comp_dic[[h]] = list("A"=as.numeric(as.character(th_A)),"B"=as.numeric(as.character(th_B)))
    
  }
  
  return(comp_dic)
}


# Boxplots for comparison between piano and lm, gSets comes from GSA_comparison
plotGSETboxplot = function(gSets,countDF){
  
  for(gs in names(gSets)){
    plot_df = countDF[,c(1,2,grep(gs,colnames(countDF)))]
    colnames(plot_df) = c("outcome","time","ES")
    plot_df$outcome = factor(plot_df$outcome)
    plot_df$time = factor(plot_df$time,levels=sort(plot_df$time))
    gp = ggplot(plot_df, aes(x = time, y = ES)) + geom_boxplot(aes(fill = outcome), alpha = 0.5) + ylab("Enrichment Score") + theme(axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12), plot.margin = unit(c(1,0.6,1,0.6), "cm")) + ggtitle(gs)
    
    
    if(!identical(gSets[[gs]]$A,numeric(0))){
      
      gp = gp + annotate("text", x = factor(gSets[[gs]]$A), y = max(plot_df$ES) + 0.03, label = "*", size=12, colour="blue")
      
    }
    
    if(!identical(gSets[[gs]]$B,numeric(0))){
      
      gp = gp + annotate("text", x = factor(gSets[[gs]]$B), y = max(plot_df$ES) + 0.06, label = "*", size=12, colour="black")
      
    }
    
    print(gp)
    
    
  }
  
}




######################################################################################################################

GSEA.GeneRanking <- function(A, class.labels, gene.labels, nperm, permutation.type = 0, sigma.correction = "GeneCluster", fraction=1.0, replace=F, reverse.sign= F) { 
  
  # This function ranks the genes according to the signal to noise ratio for the actual phenotype and also random permutations and bootstrap  
  # subsamples of both the observed and random phenotypes. It uses matrix operations to implement the signal to noise calculation 
  # in stages and achieves fast execution speed. It supports two types of permutations: random (unbalanced) and balanced. 
  # It also supports subsampling and bootstrap by using masking and multiple-count variables.  When "fraction" is set to 1 (default)
  # the there is no subsampling or boostrapping and the matrix of observed signal to noise ratios will have the same value for 
  # all permutations. This is wasteful but allows to support all the multiple options with the same code. Notice that the second 
  # matrix for the null distribution will still have the values for the random permutations 
  # (null distribution). This mode (fraction = 1.0) is the defaults, the recommended one and the one used in the examples.
  # It is also the one that has be tested more thoroughly. The resampling and boostrapping options are intersting to obtain 
  # smooth estimates of the observed distribution but its is left for the expert user who may want to perform some sanity 
  # checks before trusting the code.
  #
  # Inputs:
  #   A: Matrix of gene expression values (rows are genes, columns are samples) 
  #   class.labels: Phenotype of class disticntion of interest. A vector of binary labels having first the 1's and then the 0's 
  #   gene.labels: gene labels. Vector of probe ids or accession numbers for the rows of the expression matrix 
  #   nperm: Number of random permutations/bootstraps to perform 
  #   permutation.type: Permutation type: 0 = unbalanced, 1 = balanced. For experts only (default: 0) 
  #   sigma.correction: Correction to the signal to noise ratio (Default = GeneCluster, a choice to support the way it was handled in a previous package) 
  #   fraction: Subsampling fraction. Set to 1.0 (no resampling). For experts only (default: 1.0) 
  #   replace: Resampling mode (replacement or not replacement). For experts only (default: F) 
  #   reverse.sign: Reverse direction of gene list (default = F)
  #
  # Outputs:
  #   s2n.matrix: Matrix with random permuted or bootstraps signal to noise ratios (rows are genes, columns are permutations or bootstrap subsamplings
  #   obs.s2n.matrix: Matrix with observed signal to noise ratios (rows are genes, columns are boostraps subsamplings. If fraction is set to 1.0 then all the columns have the same values
  #   order.matrix: Matrix with the orderings that will sort the columns of the obs.s2n.matrix in decreasing s2n order
  #   obs.order.matrix: Matrix with the orderings that will sort the columns of the s2n.matrix in decreasing s2n order
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  A <- A + 0.00000001
  
  N <- length(A[,1])
  Ns <- length(A[1,])
  
  subset.mask <- matrix(0, nrow=Ns, ncol=nperm)
  reshuffled.class.labels1 <- matrix(0, nrow=Ns, ncol=nperm)
  reshuffled.class.labels2 <- matrix(0, nrow=Ns, ncol=nperm)
  class.labels1 <- matrix(0, nrow=Ns, ncol=nperm)
  class.labels2 <- matrix(0, nrow=Ns, ncol=nperm)
  
  order.matrix <- matrix(0, nrow = N, ncol = nperm)
  obs.order.matrix <- matrix(0, nrow = N, ncol = nperm)
  s2n.matrix <- matrix(0, nrow = N, ncol = nperm)
  obs.s2n.matrix <- matrix(0, nrow = N, ncol = nperm)
  
  obs.gene.labels <- vector(length = N, mode="character")
  obs.gene.descs <- vector(length = N, mode="character")
  obs.gene.symbols <- vector(length = N, mode="character")
  
  M1 <- matrix(0, nrow = N, ncol = nperm)
  M2 <- matrix(0, nrow = N, ncol = nperm)
  S1 <- matrix(0, nrow = N, ncol = nperm)
  S2 <- matrix(0, nrow = N, ncol = nperm)
  
  gc()
  
  C <- split(class.labels, class.labels)
  class1.size <- length(C[[1]])
  class2.size <- length(C[[2]])
  class1.index <- seq(1, class1.size, 1)
  class2.index <- seq(class1.size + 1, class1.size + class2.size, 1)
  
  for (r in 1:nperm) {
    class1.subset <- sample(class1.index, size = ceiling(class1.size*fraction), replace = replace)
    class2.subset <- sample(class2.index, size = ceiling(class2.size*fraction), replace = replace)
    class1.subset.size <- length(class1.subset)
    class2.subset.size <- length(class2.subset)
    subset.class1 <- rep(0, class1.size)
    for (i in 1:class1.size) {
      if (is.element(class1.index[i], class1.subset)) {
        subset.class1[i] <- 1
      }
    }
    subset.class2 <- rep(0, class2.size)
    for (i in 1:class2.size) {
      if (is.element(class2.index[i], class2.subset)) {
        subset.class2[i] <- 1
      }
    }
    subset.mask[, r] <- as.numeric(c(subset.class1, subset.class2))
    fraction.class1 <- class1.size/Ns
    fraction.class2 <- class2.size/Ns
    
    if (permutation.type == 0) { # random (unbalanced) permutation
      full.subset <- c(class1.subset, class2.subset)
      label1.subset <- sample(full.subset, size = Ns * fraction.class1)
      reshuffled.class.labels1[, r] <- rep(0, Ns)
      reshuffled.class.labels2[, r] <- rep(0, Ns)
      class.labels1[, r] <- rep(0, Ns)
      class.labels2[, r] <- rep(0, Ns)
      for (i in 1:Ns) {
        m1 <- sum(!is.na(match(label1.subset, i)))
        m2 <- sum(!is.na(match(full.subset, i)))
        reshuffled.class.labels1[i, r] <- m1
        reshuffled.class.labels2[i, r] <- m2 - m1
        if (i <= class1.size) {
          class.labels1[i, r] <- m2
          class.labels2[i, r] <- 0
        } else {
          class.labels1[i, r] <- 0
          class.labels2[i, r] <- m2
        }
      }
      
    } else if (permutation.type == 1) { # proportional (balanced) permutation
      
      class1.label1.subset <- sample(class1.subset, size = ceiling(class1.subset.size*fraction.class1))
      class2.label1.subset <- sample(class2.subset, size = floor(class2.subset.size*fraction.class1))
      reshuffled.class.labels1[, r] <- rep(0, Ns)
      reshuffled.class.labels2[, r] <- rep(0, Ns)
      class.labels1[, r] <- rep(0, Ns)
      class.labels2[, r] <- rep(0, Ns)
      for (i in 1:Ns) {
        if (i <= class1.size) {
          m1 <- sum(!is.na(match(class1.label1.subset, i)))
          m2 <- sum(!is.na(match(class1.subset, i)))
          reshuffled.class.labels1[i, r] <- m1
          reshuffled.class.labels2[i, r] <- m2 - m1
          class.labels1[i, r] <- m2
          class.labels2[i, r] <- 0
        } else {
          m1 <- sum(!is.na(match(class2.label1.subset, i)))
          m2 <- sum(!is.na(match(class2.subset, i)))
          reshuffled.class.labels1[i, r] <- m1
          reshuffled.class.labels2[i, r] <- m2 - m1
          class.labels1[i, r] <- 0
          class.labels2[i, r] <- m2
        }
      }
    }
  }
  
  # compute S2N for the random permutation matrix
  
  P <- reshuffled.class.labels1 * subset.mask
  n1 <- sum(P[,1])         
  M1 <- A %*% P
  M1 <- M1/n1      
  gc()
  A2 <- A*A        
  S1 <- A2 %*% P   
  S1 <- S1/n1 - M1*M1    
  S1 <- sqrt(abs((n1/(n1-1)) * S1))   
  gc()
  P <- reshuffled.class.labels2 * subset.mask
  n2 <- sum(P[,1])           
  M2 <- A %*% P           
  M2 <- M2/n2          
  gc()
  A2 <- A*A           
  S2 <- A2 %*% P      
  S2 <- S2/n2 - M2*M2 
  S2 <- sqrt(abs((n2/(n2-1)) * S2))
  rm(P)
  rm(A2)
  gc()
  
  if (sigma.correction == "GeneCluster") {  # small sigma "fix" as used in GeneCluster
    S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
    S2 <- ifelse(S2 == 0, 0.2, S2)
    S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
    S1 <- ifelse(S1 == 0, 0.2, S1)
    gc()
  }
  
  M1 <- M1 - M2
  rm(M2)
  gc()
  S1 <- S1 + S2
  rm(S2)
  gc()
  
  s2n.matrix <- M1/S1
  
  if (reverse.sign == T) {
    s2n.matrix <- - s2n.matrix
  }
  gc()
  
  for (r in 1:nperm) {
    order.matrix[, r] <- order(s2n.matrix[, r], decreasing=T)            
  }
  
  # compute S2N for the "observed" permutation matrix
  
  P <- class.labels1 * subset.mask
  n1 <- sum(P[,1])         
  M1 <- A %*% P
  M1 <- M1/n1      
  gc()
  A2 <- A*A        
  S1 <- A2 %*% P   
  S1 <- S1/n1 - M1*M1    
  S1 <- sqrt(abs((n1/(n1-1)) * S1))   
  gc()
  P <- class.labels2 * subset.mask
  n2 <- sum(P[,1])           
  M2 <- A %*% P           
  M2 <- M2/n2          
  gc()
  A2 <- A*A           
  S2 <- A2 %*% P      
  S2 <- S2/n2 - M2*M2 
  S2 <- sqrt(abs((n2/(n2-1)) * S2))
  rm(P)
  rm(A2)
  gc()
  
  if (sigma.correction == "GeneCluster") {  # small sigma "fix" as used in GeneCluster
    S2 <- ifelse(0.2*abs(M2) < S2, S2, 0.2*abs(M2))
    S2 <- ifelse(S2 == 0, 0.2, S2)
    S1 <- ifelse(0.2*abs(M1) < S1, S1, 0.2*abs(M1))
    S1 <- ifelse(S1 == 0, 0.2, S1)
    gc()
  } 
  
  M1 <- M1 - M2
  rm(M2)
  gc()
  S1 <- S1 + S2
  rm(S2)
  gc()
  
  obs.s2n.matrix <- M1/S1
  gc()
  
  if (reverse.sign == T) {
    obs.s2n.matrix <- - obs.s2n.matrix
  }
  
  for (r in 1:nperm) {
    obs.order.matrix[,r] <- order(obs.s2n.matrix[,r], decreasing=T)            
  }
  
  return(list(s2n.matrix = s2n.matrix, 
              obs.s2n.matrix = obs.s2n.matrix, 
              order.matrix = order.matrix,
              obs.order.matrix = obs.order.matrix))
}

#################################################
#
# Functions to cluster patients
#
#

#This function creates a matrix that merges time point measurements in a single vector
meltCountMatrix = function(targetDF,countMat,selFeat,meltType = "both"){
  
  meltMat = c()
  
  if(meltType=="case"){
    targetDF = filter(targetDF,outcome==1)
  }else if(meltType=="control"){
    targetDF = filter(targetDF,outcome==0)
  }
  
  patients = unique(targetDF$mask_id)
  
  for(patient in patients){
    
    patient_row = c()
    patient_samples = filter(targetDF,mask_id == patient)
    
    for(t in unique(sort(targetDF$time))){
      
      if(t %in% patient_samples$time){
        
        sample_id = patient_samples[patient_samples$time==t,"sample_mask_id"]
        expr_vec = countMat[selFeat,as.character(sample_id)]
        names(expr_vec) = paste(selFeat,"_t",as.character(t),sep="")
        patient_row = c(patient_row, expr_vec)
        
      } else{
        
        expr_vec = rep(NA, length(selFeat))
        names(expr_vec) = paste(selFeat,"_t",as.character(t),sep="")
        patient_row = c(patient_row, expr_vec)
      }
      
    }
    
    meltMat = rbind(meltMat,patient_row)
    
  }
  
  rownames(meltMat) = patients
  return(meltMat)
}

# Get results from hclust
getClusters = function(hclust_res,k){
  clusters = cutree(hclust_res,k=k)
  kList = list()
  for(i in 1:4){
    K = names(clusters)[clusters==i]
    kList[[i]] = K
  }
  return(kList)
}


compareKdynamics = function(meltMat,annotation_table,FeatSet){
  #meltMat: obtained from meltCountMatrix
  #annotation_table: it contains mask_ids as rownames, columns should be: mask_id, outcome, cluster number 
  for(f in FeatSet){
    
    plotDF = c()
    FeatMat = meltMat[,grep(paste("^",f,sep=""),colnames(meltMat))]
    Time_vec = strsplit2(colnames(FeatMat),"_t")[,2]
    
    for(C in 1:ncol(FeatMat)){
      timeDF = cbind(cbind(annotation_table,FeatMat[rownames(annotation_table),C]),Time_vec[C])
      plotDF = rbind(plotDF,timeDF)
    }
    
    colnames(plotDF) = c("patient","outcome","cluster","EnrichmentS","Time")
    plotDF = data.frame(plotDF,stringsAsFactors = F)
    plotDF$patient = factor(plotDF$patient)
    plotDF$outcome = factor(plotDF$outcome)
    plotDF$cluster = factor(plotDF$cluster)
    plotDF$Time = as.numeric(as.character(plotDF$Time)) * -1
    
    GEX_plot <- ggplot(plotDF, aes(x = Time, y = EnrichmentS, color = outcome, group = patient)) + 
      labs(y = "Enrichment Score") + geom_point() + geom_line() + scale_x_continuous(name="Time",breaks=c(-12,-9,-6,-3,0)) + ylab("Expression") + theme(plot.margin = unit(c(1,0.6,1,0.6), "cm"),axis.text.x  = element_text(size=12),axis.text.y  = element_text(size=12)) + ggtitle(f) + facet_grid(cluster ~ .)
    
    print(GEX_plot)
  }
  
}

getclusterDictionary = function(dir.name,class.name){
  # This function takes a folder with cluster results and creates a dictionary of clusters divided by partition
  # Remember to use the results from "clusteringPatients.R"
  
  f.names = list.files(dir.name)
  f.names = f.names[grep(".txt",f.names)]
  c.dic = list()
  
  for(f.name in f.names){
    
    partition = strsplit(f.name,"_")[[1]]
    partition = gsub(".txt","",partition[length(partition)])
    
    kdata = read.delim(file=paste(dir.name,f.name,sep="/"),header = T,stringsAsFactors = F)
    
    for(k in unique(kdata$K)){
      
      klabel = paste(class.name,paste("p",partition,sep=""),paste("k",(as.character(k)),sep=""),sep = "_")
      c.dic[[klabel]] = kdata$mask_id[kdata$K==k]
      
    }
    
  }
  
  return(c.dic)
  
}

###########################################################
# Methods to identify radical expression values per patient
###########################################################

t.test.p.value = function(mu,dist.vec){
  t.test.res = t.test(dist.vec,mu=mu)
  return(t.test.res$p.value)
}

t.test.statistic = function(mu,dist.vec){
  t.test.res = t.test(dist.vec,mu=mu)
  return(t.test.res$statistic)
}

radical.test = function(Sval,cutoff.vals){
  #cutoff.vals = quantile(dist.vec,c(0.05,0.95))
  if(Sval <= cutoff.vals[1]){
    return(-1)
  }else if(Sval >= cutoff.vals[2]){
    return(1)
  }else{
    return(0)
  }
}

###############################################################3
# Method that generates exhaustively feature signatures 
###############################################################

getSignatures = function(cMAT,targetMAT,Features=NULL,time,g.pval,s.pval = NULL,model_type = "limma"){
  # For limma and GLM: It creates block of patients by Features + their combinations
  # Inputs:
  # - Count_matrix (cMAT) = count matrix containing the samples you wish to subset and model
  # - Annotation_file (targetMAT) = annotation of the count matrix
  # - Features = dictionary that contains the name of the features and the values to be used in the combinations
  # - Time = times to be used in the analysis
  # - globalTreshold = pval to be used
  #
  # Outputs:
  # - Signatures and direction of the change
  # - Results from the linear model
  # - Data used for each model and its annotation
  #
  
  if(is.null(Features)){ #If features is NULL, it creates time-specific models for the selected data
    # Objects that store the information of each test (Data,annotation,results,signatures)
    signatureDictionary = list()
    resultsDictionary = list()
    dataDictionary = list()
    annotDictionary = list()
    tag = "g"
    
    for(tm in time){ #Subselection divided by time
      tTargets = FilterTargets(targetTable = targetMAT,query_class= c("time"), query_values = c(tm),OR = FALSE) #Function that makes the filtering
      Ncases = sum(tTargets$outcome==1)
      Ncontrols = sum(tTargets$outcome==0)
      
      if(Ncases>2 & Ncontrols>2){
        tdata = cMAT[,as.character(tTargets$sample_mask_id)]
        dataDictionary[[paste(tag,"time",as.character(tm),sep="_")]] = tdata
        annotDictionary[[paste(tag,"time",as.character(tm),sep="_")]] = tTargets
      }
    }
    
    
    if(model_type == "limma"){
      
      for(dselection in names(dataDictionary)){
        print(paste("Fitting model for...",dselection,"using",model_type))
        resultsDictionary[[dselection]] = runTEDDYlimma_basic(stargets = annotDictionary[[dselection]], sEXPMAT = dataDictionary[[dselection]])
        signatureDictionary[[dselection]] = filterLimmaGLSv3(GLS = resultsDictionary[[dselection]], pTRSH = g.pval,FC_col = "logFC")
      }
      
    } else if(model_type == "glm"){
      
      for(dselection in names(dataDictionary)){
        print(paste("Fitting model for...",dselection,"using",model_type))
        resultsDictionary[[dselection]] = fitGLM_v2(cTargets = annotDictionary[[dselection]], cMAT = dataDictionary[[dselection]])
        signatureDictionary[[dselection]] = filterGLM(GLMres = resultsDictionary[[dselection]],THRS = g.pval)
      }
      
    }
    
    return(list("DATA"=dataDictionary,"ANNOTATIONS"=annotDictionary,"RESULTS"=resultsDictionary,"SIGNATURES"=signatureDictionary))
    
  } else{ #If youy give a list of features for combinatorial analysis
    
    #1st, make a selection of all data and annotation
    
    signatureDictionary = list()
    resultsDictionary = list()
    dataDictionary = list()
    annotDictionary = list()
    tag = "g"
    
    for(tm in time){ #Subselection divided by time
      tTargets = FilterTargets(targetTable = targetMAT,query_class= c("time"), query_values = c(tm),OR = FALSE) #Function that makes the filtering
      Ncases = sum(tTargets$outcome==1)
      Ncontrols = sum(tTargets$outcome==0)
      
      if(Ncases>2 & Ncontrols>2){
        tdata = cMAT[,as.character(tTargets$sample_mask_id)]
        dataDictionary[[paste(tag,"time",as.character(tm),sep="_")]] = tdata
        annotDictionary[[paste(tag,"time",as.character(tm),sep="_")]] = tTargets
      }
    }
    
    #2nd, generate all the combinations
    
    Feat_names = names(Features)
    for(i in 1:length(Feat_names)){
      #Generate combinations of feauture classes, containing 1 to n number of features
      feature_class = combinations((length(Feat_names)),i,Feat_names)
      
      for(ir in 1:nrow(feature_class)){ #For each combination of features, obtain a list and get all the combinations of their values
        
        useful_vars = Features[feature_class[ir,]]
        useful_vars = c(useful_vars,list("time"=time))
        combination_df = expand.grid(useful_vars,stringsAsFactors = F)
        filtFeatures = colnames(combination_df)
        
        for(iq in 1:nrow(combination_df)){ #for every combination, filter patients
          
          query_values = as.character(combination_df[iq,])
          tTargets = FilterTargets(targetTable = targetMAT,query_class= filtFeatures, query_values = query_values,OR = FALSE) #Function that makes the filtering
          
          Ncases = sum(tTargets$outcome==1)
          Ncontrols = sum(tTargets$outcome==0)
          
          if(Ncases>2 & Ncontrols>2){
            
            tag = paste(paste(filtFeatures,query_values,sep="_"),collapse = ".")
            tdata = cMAT[,as.character(tTargets$sample_mask_id)]
            dataDictionary[[tag]] = tdata
            annotDictionary[[tag]] = tTargets
          }
        }
      }
    }
    
    # Apply_methods
    
    if(model_type == "limma"){
      
      for(dselection in names(dataDictionary)){
        print(paste("Fitting model for...",dselection,"using",model_type))
        
        if(identical(grep("^g",dselection),integer(0))){
          resultsDictionary[[dselection]] = runTEDDYlimma_basic(stargets = annotDictionary[[dselection]], sEXPMAT = dataDictionary[[dselection]])
          signatureDictionary[[dselection]] = filterLimmaGLSv3(GLS = resultsDictionary[[dselection]], pTRSH = s.pval,FC_col = "logFC")
        } else{
          resultsDictionary[[dselection]] = runTEDDYlimma_basic(stargets = annotDictionary[[dselection]], sEXPMAT = dataDictionary[[dselection]])
          signatureDictionary[[dselection]] = filterLimmaGLSv3(GLS = resultsDictionary[[dselection]], pTRSH = g.pval,FC_col = "logFC")
        }
      }
      
    } else if(model_type == "glm"){
      
      for(dselection in names(dataDictionary)){
        print(paste("Fitting model for...",dselection,"using",model_type))
        if(identical(grep("^g",dselection),integer(0))){
          resultsDictionary[[dselection]] = fitGLM_v2(cTargets = annotDictionary[[dselection]], cMAT = dataDictionary[[dselection]])
          signatureDictionary[[dselection]] = filterGLM(GLMres = resultsDictionary[[dselection]],THRS = s.pval)
        } else{
          resultsDictionary[[dselection]] = fitGLM_v2(cTargets = annotDictionary[[dselection]], cMAT = dataDictionary[[dselection]])
          signatureDictionary[[dselection]] = filterGLM(GLMres = resultsDictionary[[dselection]],THRS = g.pval)
        }
        
      }
      
    }
    
    return(list("DATA"=dataDictionary,"ANNOTATIONS"=annotDictionary,"RESULTS"=resultsDictionary,"SIGNATURES"=signatureDictionary))
    
  }
  
}

unlistSignatures = function(Signature_Dictionary){
  
  All_signatures_dictionary_unlist = list()
  for(sign in names(Signature_Dictionary)){
    
    upT = paste(sign,"up",sep="_")
    downT = paste(sign,"down",sep="_")
    Flevel = Signature_Dictionary[[sign]]
    
    if(!is.null(Flevel)){
      if(!identical(Flevel[["up"]],character(0))){
        All_signatures_dictionary_unlist[[upT]] = Flevel[["up"]]
      }
      if(!identical(Flevel[["down"]],character(0))){
        All_signatures_dictionary_unlist[[downT]] = Flevel[["down"]]
      }
    }
  }
  
  return(All_signatures_dictionary_unlist)
  
}

plot_MDS_Tv = function(meltMat,targetDF,main){
  
  #dist_pat = as.matrix(daisy(meltMat,metric = "gower"))
  dist_pat = as.matrix(dist(meltMat))
  dist_pat_temp = dist_pat
  test_coef = as.numeric(names(table(rowSums(is.na(dist_pat)))))
  
  maxSamples = c()
  for(coef in test_coef){
    dist_pat = dist_pat_temp
    dist_pat = dist_pat[rowSums(is.na(dist_pat))<coef,rowSums(is.na(dist_pat))<coef]
    dist_pat = dist_pat[rowSums(is.na(dist_pat))<1,rowSums(is.na(dist_pat))<1]
    maxSamples = c(maxSamples,dim(dist_pat)[1])
  }
  
  coef = test_coef[which(maxSamples == max(maxSamples))]
  dist_pat = dist_pat_temp
  dist_pat = dist_pat[rowSums(is.na(dist_pat))<coef,rowSums(is.na(dist_pat))<coef]
  dist_pat = dist_pat[rowSums(is.na(dist_pat))<1,rowSums(is.na(dist_pat))<1]
  
  mds = cmdscale(as.dist(dist_pat))
  
  pat_annotation = as.matrix(unique(targetDF[,c("mask_id","outcome")]))
  rownames(pat_annotation) = pat_annotation[,1] 
  pat_annotation = pat_annotation[,-1]
  pat_annotation = ifelse(pat_annotation==1,"case","control")
  pat_color = ifelse(pat_annotation=="case","red","blue")
  
  
  plot(mds, type = 'p',pch=21,main=main,col = pat_color[rownames(mds)],bg = pat_color[rownames(mds)])
  #text(mds[, 1], mds[, 2], pat_annotation[rownames(mds)],col = pat_color[rownames(mds)])
  
}

### Alphabet Alignment ###

equalSets = function(setA,setB){
  # Comparison of sets of equal length
  
  if(length(setA)==length(setB)){
    EVAL_val = length(intersect(setA,setB)) == length(setA)
    return(EVAL_val)
  }else{
    print("SETS don't have the same length")
    return()
  }
  
}

getAlphabet = function(pat_counts,CountValue = "max"){
  if(CountValue == "max"){
    max_Letter = names(Alphabet_dict)[sapply(Alphabet_dict,equalSets,setB = names(pat_counts[which(pat_counts == max(pat_counts,na.rm = T))]))]
    return(max_Letter)
  }else if(CountValue == "min"){
    min_Letter = names(Alphabet_dict)[sapply(Alphabet_dict,equalSets,setB = names(pat_counts[which(pat_counts == min(pat_counts,na.rm = T))]))]
    return(min_Letter)
  }
}

plotMDS_stringV = function(str_matrix,targetDF,main,color = "outcome",shape = "FirstAAb"){
  
  StringVector = apply(str_matrix,1,paste,collapse="")
  dist_pat = as.matrix(stringdistmatrix(StringVector))
  colnames(dist_pat) = rownames(str_matrix)
  rownames(dist_pat) = rownames(str_matrix)
  
  mds = cmdscale(as.dist(dist_pat))
  colnames(mds) = c("X","Y")
  pat_annotation = na.omit(as.matrix(unique(targetDF[,c("mask_id","outcome","FirstAAb","agegroup","gender")])))
  rownames(pat_annotation) = pat_annotation[,1] 
  pat_annotation = pat_annotation[,-1]
  
  
  plotDF = data.frame(cbind(mds,StringVector,pat_annotation[rownames(mds),]))
  plotDF$X = as.numeric(as.character(plotDF$X))
  plotDF$Y = as.numeric(as.character(plotDF$Y))
  
  sp = ggplot(plotDF, aes(x = X, y = Y,label = StringVector)) + geom_point(aes(shape = plotDF[,shape],color = plotDF[,color],size=5)) + geom_text(hjust = 0, nudge_x = 0.05) + labs(title=main)
  print(sp)
}

# Functions used for statistics annotations of words

getLetterStats = function(combination_name,typeof="max"){
  str_matrix = group_results[[combination_name]][[typeof]]
  targetDF = sample_dictionary[[combination_name]]
  case_mat = str_matrix[unique(as.character(targetDF$mask_id[targetDF$outcome==1])),]
  case_recurrence = lapply(lapply(apply(case_mat,2,table),sort),nameStats)
  
  names(case_recurrence) = unlist(Alphabet_dict[names(case_recurrence)])
  
  control_mat = str_matrix[unique(as.character(targetDF$mask_id[targetDF$outcome==0])),]
  control_recurrence = lapply(lapply(apply(control_mat,2,table),sort),nameStats)
  names(control_recurrence) = unlist(Alphabet_dict[names(control_recurrence)])
  
  return(list("case_recurrence"=case_recurrence,"control_recurrence"=control_recurrence))
}

nameStats = function(l){
  names(l) = unlist(Alphabet_dict[names(l)])
  return(l)
}

#################################################

read_DAVID = function(file_path){
  
  DavidAnnotation = read.delim(file_path,blank.lines.skip = TRUE,comment.char = "#",header = F,stringsAsFactors = F)
  colnames(DavidAnnotation) = DavidAnnotation[1,]
  DavidAnnotation = DavidAnnotation[!apply(DavidAnnotation,1,function(x){
    return(sum(x==DavidAnnotation[1,])==ncol(DavidAnnotation))
  }),]
  
  return(DavidAnnotation)
}

#################
binning_targets = function(target_file){
  new_time = round(target_file[,"agemonth"] - target_file[, "case_endptage_month"])
  new_time_binning = round(new_time/3)
  return(new_time_binning)
}

################

##################
#
# MASIGPRO FUNCTIONS
#
##################

runMASIGPRO_pop = function(targetTable, factor_column = "outcome",factor_name = "OUTCOME",time_column = "time",sample_column = "sample_mask_id", 
                           ExpressionCounts, Dgree, FDR, R2, K, min.obs = 4,
                           postAnalysis = TRUE, testTargets,testCounts, 
                           generatePlots=TRUE,foldername,fileID){
  
  # Generates edesign matrix: targetTable, factor_column = "outcome",factor_name = "OUTCOME",time_column = "time",sample_column = "sample_mask_id"
  # Runs MASIGPRO: ExpressionCounts, Dgree, FDR, R2
  # Generates the number of clusters you need
  # Runs ENrichment analysis of defined clusters
  # Runs GLMs to show significant differences in expression
  #
  
  # Create eDESIGN matrix
  edesign = masigPro_edesign(target_mat = targetTable,factor_column = "outcome",factor_name = "OUTCOME",time_column = "time",sample_column = "sample_mask_id")
  ExpressionCounts = ExpressionCounts[,as.character(rownames(edesign))]
  
  # Run MSpro
  design <- make.design.matrix(edesign, degree = Dgree)
  fit <- p.vector(ExpressionCounts, design, Q = FDR, MT.adjust = "BH", min.obs = min.obs)
  tstep = T.fit(fit, step.method = "backward", alfa = FDR)
  sigs <- get.siggenes(tstep, rsq = R2, vars = "groups")
  
  clusters = see.genes(sigs$sig.genes[[2]], show.fit = T, dis =design$dis,cluster.method="hclust" ,cluster.data = 1, k = K)
  
  
  SignatureDictionary = list()
  for(i in unique(clusters$cut)){
    tagN = paste("K",as.character(i),sep="_")
    SignatureDictionary[[tagN]] = names(clusters$cut)[clusters$cut==i]
  }
  
  # Automatically generates list of genes
  
  folderPATH = paste(getwd(),foldername,fileID,sep="/")
  
  lapply(names(SignatureDictionary),function(x){
    #AnnotationTable_LD = cbind(sort(DiabetesStages_unlist$Late_Down),data.frame(mapIds(illuminaHumanv4.db, keys=sort(DiabetesStages_unlist$Late_Down), keytype = "SYMBOL", column = "GENENAME"),stringsAsFactors = F))
    write.table(file=paste(folderPATH,"_",x,"_ENTREZ",".txt",sep=""),
                (na.omit(data.frame(mapIds(illuminaHumanv4.db, keys=SignatureDictionary[[x]], keytype = "SYMBOL", column = "ENTREZID",multiVals = "first"),stringsAsFactors = F)))[,1],
                quote = F,sep = "\t",row.names = F,col.names = F)
    write.table(file=paste(folderPATH,"_",x,".txt",sep=""), SignatureDictionary[[x]],quote = F,sep = "\t",row.names = F,col.names = F)
  })
  
  
  #################################################
  # After obtaining gene sets with similar dynamics
  # I calculate the Enrichment Scores for these signatures
  # in all time points!
  
  if(postAnalysis){
    
    # Create Function Matrix
    
    LM_globalFMAT = getFunctionMatrix_v2(eMat = testCounts,eTargets = testTargets,Signature_Dict = SignatureDictionary,weighted = TRUE)
    LM_globalFMAT = LM_globalFMAT$ESmat
    
    # Filter times with low number of samples
    
    completeD = names(table(testTargets$time))[table(testTargets$time) >=8]
    LM_globalTargets = filter(testTargets, time %in% completeD)
    LM_globalFMAT = LM_globalFMAT[,as.character(LM_globalTargets$sample_mask_id)]
    
    # Get results from T-test
    masigClusters = getSignatures(cMAT = LM_globalFMAT,targetMAT = LM_globalTargets,Features = NULL,time = sort(unique(LM_globalTargets$time)),g.pval = 0.05,s.pval=0.05,model_type = "glm")
    # Create matrix of results
    resultData = masigClusters$RESULTS
    
    if(generatePlots){
      pdf(file = paste(folderPATH,"Dynamics.pdf",sep = "_"),height = 10,width = 12)
      for(k in rownames(LM_globalFMAT)){
        plotES(ES_list = LM_globalFMAT[k,],sigTargets = LM_globalTargets,TRSH = 0.05,main_title = k)
      }
      dev.off()
      
      pdf(file = paste(folderPATH,"LMresults.pdf",sep = "_"),height = 8,width = 8)
      SummaryTables = summaryTestTable(resultData)
      dev.off()
    }
    
    return(list("ClusterList" = SignatureDictionary, "ES_Matrix" = LM_globalFMAT,"GLM_RESULTS" = resultData))
    
  } else{
    
    return(SignatureDictionary)
    
  }
  
}

#########################################

summaryTestTable = function(resultData){
  #Function that plots the coefficients and pvalues of a GLM Result list object
  PVAL_mat = matrix(1,nrow = nrow((resultData[[1]])), ncol = length(resultData))
  COEF_mat = matrix(1,nrow = nrow((resultData[[1]])), ncol = length(resultData))
  
  rownames(PVAL_mat) = rownames(COEF_mat) = resultData[[1]]$Hallmark
  colnames(PVAL_mat) = colnames(COEF_mat) = names(resultData)
  
  for(Tres in names(resultData)){
    
    COEF_mat[,Tres] = resultData[[Tres]]$Coefficient
    PVAL_mat[,Tres] = resultData[[Tres]]$adj_P_val
    
  }
  
  COEF_mat = round(COEF_mat,2)
  PVAL_mat = round(PVAL_mat,2)
  
  my_palette <- colorRampPalette(c("lightcoral", "white", "lightgreen"))(n = 20)
  
  PVALplot = heatmap.2(PVAL_mat,Rowv = FALSE,Colv = FALSE,dendrogram = 'none',cellnote = PVAL_mat,trace = "none",notecol = "black",col=my_palette,lwid=  c(0.05,.5),lhei= c(0.05,.30),margins(30,15),key = F,srtCol=35)
  COEFplot = heatmap.2(COEF_mat,Rowv = FALSE,Colv = FALSE,dendrogram = 'none',cellnote = COEF_mat,trace = "none",notecol = "black",col=my_palette,lwid=  c(0.05,.5),lhei= c(0.05,.30),margins(30,15),key = F,srtCol=35)
  
  print(PVALplot)
  print(COEFplot)
  
  return(list("pval" = PVAL_mat,"coef" = COEF_mat))
}

# Modificantions of MASIGPRO

make.design.matrix = function (edesign, degree = 2, time.col = 1, repl.col = 2, group.cols = c(3:ncol(edesign))) 
{
  #TO DO: CAST THE MATRIX INTO A NUMERIC MATRIX
  control.label <- colnames(edesign)[group.cols][1]
  if (dim(as.matrix(edesign))[2] > 3) {
    dummy.cols <- group.cols[2:length(group.cols)]
    treatm.label <- paste(colnames(edesign)[dummy.cols], 
                          "vs", control.label, sep = "")
    groups.label <- c(control.label, treatm.label)
    matrix.dummy <- as.matrix(edesign[, dummy.cols])
    dummy <- NULL
    j = 0
    origen <- min(edesign[, time.col])
    origen <- edesign[edesign[, 1] == origen,,drop=F]
    for (i in 1:length(dummy.cols)) {
      share <- apply(origen[, c(3, dummy.cols[i]),drop=F], 1, sum)
      if (!is.element(TRUE, share > 1)) {
        j = j + 1
        dummy <- cbind(dummy, matrix.dummy[, i])
        colnames(dummy)[j] <- treatm.label[i]
      }
    }
    time <- as.matrix(edesign[, time.col])
    colnames(time) <- colnames(edesign)[time.col]
    dis <- cbind(dummy, time)
    rownames(dis) <- rownames(edesign)
    groups.vector <- c(colnames(dummy), control.label)
    colnames.dis <- colnames(dis)
    dis <- cbind(dis, dis[, ncol(dis)] * matrix.dummy)
    colnames(dis) <- c(colnames.dis, paste(colnames(edesign)[time.col], 
                                           "x", colnames(edesign)[dummy.cols], sep = ""))
    groups.vector <- c(groups.vector, treatm.label)
    if (degree >= 2) {
      for (i in 2:degree) {
        colnames.dis <- colnames(dis)
        dis <- cbind(dis, edesign[, time.col]^i, edesign[, 
                                                         time.col]^i * edesign[, dummy.cols])
        colnames(dis) <- c(colnames.dis, paste(colnames(edesign)[time.col], 
                                               i, sep = ""), paste(colnames(edesign)[time.col], 
                                                                   "", i, "x", colnames(edesign)[dummy.cols], 
                                                                   sep = ""))
        groups.vector <- c(groups.vector, groups.label)
      }
    }
  }
  else {
    dis <- as.matrix(edesign[, time.col])
    colnames(dis) <- colnames(edesign)[time.col]
    rownames(dis) <- rownames(edesign)
    if (degree > 1) {
      for (i in 2:degree) {
        colnames.dis <- colnames(dis)
        dis <- cbind(dis, edesign[, time.col]^i)
        colnames(dis) <- c(colnames.dis, paste(colnames(edesign)[time.col], 
                                               i, sep = ""))
      }
    }
    groups.vector <- rep(colnames(edesign)[group.cols], 
                         degree)
  }
  output <- list(dis, groups.vector, edesign)
  names(output) <- c("dis", "groups.vector", "edesign")
  output
}

##############################

see.genesv2 = function (showPlots=FALSE,data, edesign = data$edesign, time.col = 1, repl.col = 2, 
                        group.cols = c(3:ncol(edesign)), names.groups = colnames(edesign)[3:ncol(edesign)], 
                        cluster.data = 1, groups.vector = data$groups.vector, k = 9, 
                        m = 1.45, cluster.method = "hclust", distance = "cor", agglo.method = "ward.D", 
                        show.fit = FALSE, dis = NULL, step.method = "backward", 
                        min.obs = 3, alfa = 0.05, nvar.correction = FALSE, show.lines = TRUE, 
                        iter.max = 500, summary.mode = "median", color.mode = "rainbow", 
                        cexlab = 1, legend = TRUE, newX11 = TRUE, ylim = NULL, main = NULL, 
                        ...) 
{
  time = edesign[, time.col]
  repvect = edesign[, repl.col]
  groups = edesign[, group.cols]
  narrays <- length(time)
  if (!is.null(dim(data))) {
    dat <- as.data.frame(data)
    clusterdata <- data
  }
  else {
    clusterdata <- data[[cluster.data]]
    dat <- as.data.frame(data$sig.profiles)
  }
  clusterdata <- clusterdata
  if (nrow(dat) > 1) {
    dat <- as.data.frame(dat[, (ncol(dat) - length(time) + 
                                  1):ncol(dat)])
    count.na <- function(x) length(x[is.na(x)])
    NAs <- apply(as.matrix(dat), 1, count.na)
    count.noNa <- function(x) (length(x) - length(x[is.na(x)]))
    dat <- dat[which(apply(as.matrix(dat), 1, count.noNa) >= 
                       2), ]
  }
  else {
    NAs <- 1
  }
  kdata <- NULL
  out <- TRUE
  if (nrow(dat) > 1) {
    if (cluster.data != 1 || cluster.data != "sig.profiles") {
      if (any(is.na(clusterdata))) 
        clusterdata[is.na(clusterdata)] <- 0
    }
    else if (is.na(all(dist(clusterdata) > 0)) || (cluster.method == 
                                                   "kmeans" & any(is.na(clusterdata))) || (distance == 
                                                                                           "cor" & any(sd(t(clusterdata), na.rm = TRUE) == 
                                                                                                       0))) {
      if (!is.null(kdata)) {
        clusterdata <- kdata
      }
      else {
        clusterdata <- NULL
      }
    }
    clusterdata <- clusterdata
    if (!is.null(clusterdata)) {
      k <- min(k, nrow(dat), na.rm = TRUE)
      if (cluster.method == "hclust") {
        if (distance == "cor") {
          dcorrel <- matrix(rep(1, nrow(clusterdata)^2), 
                            nrow(clusterdata), nrow(clusterdata)) - 
            cor(t(clusterdata), use = "pairwise.complete.obs")
          clust <- hclust(as.dist(dcorrel), method = agglo.method)
          c.algo.used = paste(cluster.method, "cor", 
                              agglo.method, sep = "_")
        }
        else {
          clust <- hclust(dist(clusterdata, method = distance), 
                          method = agglo.method)
          c.algo.used = paste(cluster.method, distance, 
                              agglo.method, sep = "_")
        }
        cut <- cutree(clust, k = k)
      }
      else if (cluster.method == "kmeans") {
        cut <- kmeans(clusterdata, k, iter.max)$cluster
        c.algo.used = paste("kmeans", k, iter.max, sep = "_")
      }
      else if (cluster.method == "mfuzz") {
        n <- dim(clusterdata)[2]
        clusterdata[is.na(clusterdata)] <- 0
        temp <- tempfile()
        write.table(clusterdata, temp, quote = FALSE, 
                    sep = "\t", row.names = TRUE, col.names = TRUE)
        signif <- readExpressionSet(temp)
        cl <- mfuzz(signif, c = k, m = m)
        clus <- acore(signif, cl = cl, min.acore = (1/k))
        for (i in 1:k) {
          clus[[i]] <- transform(clus[[i]], cluster = i)
        }
        cut0 <- clus[[1]][, c(1, 3)]
        for (i in 2:k) {
          cut0 <- rbind(cut0, clus[[i]][, c(1, 3)])
        }
        cut <- transform(clusterdata, name = "")
        cut <- transform(cut, cluster = 0)
        cut <- cut[, c(n + 1, n + 2)]
        cut[, 1] <- rownames(cut)
        for (i in 1:dim(clusterdata)[1]) {
          cut[i, 2] <- cut0[cut[i, 1], 2]
        }
        cut <- cut[, 2]
        c.algo.used = paste("mfuzz", k, m, sep = "_")
      }
      else stop("Invalid cluster algorithm")
      
      groups <- as.matrix(groups)
      colnames(groups) <- names.groups
      
      if(showPlots){
        if (newX11)
          X11()
        if (k <= 4)
          par(mfrow = c(2, 2))
        else if (k <= 6)
          par(mfrow = c(3, 2))
        else if (k > 6)
          par(mfrow = c(3, 3))
        for (i in 1:(k)) {
          PlotProfiles(data = dat[cut == i, ], repvect = repvect,
                       main = i, ylim = ylim, color.mode = color.mode,
                       cond = rownames(edesign), ...)
        }
        if (newX11)
          X11()
        if (k <= 4) {
          par(mfrow = c(2, 2))
          cexlab = 0.6
        }
        else if (k <= 6) {
          par(mfrow = c(3, 2))
          cexlab = 0.6
        }
        else if (k > 6) {
          par(mfrow = c(3, 3))
          cexlab = 0.35
        }
        for (j in 1:(k)) {
          PlotGroups(data = dat[cut == j, ], show.fit = show.fit,
                     dis = dis, step.method = step.method, min.obs = min.obs,
                     alfa = alfa, nvar.correction = nvar.correction,
                     show.lines = show.lines, time = time, groups = groups,
                     repvect = repvect, summary.mode = summary.mode,
                     xlab = "time", main = paste("Cluster", j,
                                                 sep = " "), ylim = ylim, cexlab = cexlab,
                     legend = legend, groups.vector = groups.vector,
                     ...)
        }
      }
      
    }
    else {
      print("warning: impossible to compute hierarchical clustering")
      c.algo.used <- NULL
      cut <- 1
    }
  }
  else if (nrow(dat) == 1) {
    if(showPlots){
      if (newX11)
        X11()
      PlotProfiles(data = dat, repvect = repvect, main = NULL,
                   ylim = ylim, color.mode = color.mode, cond = rownames(edesign),
                   ...)
      if (newX11)
        X11()
      PlotGroups(data = dat, show.fit = show.fit, dis = dis,
                 step.method = step.method, min.obs = min.obs, alfa = alfa,
                 nvar.correction = nvar.correction, show.lines = show.lines,
                 time = time, groups = groups, repvect = repvect,
                 summary.mode = summary.mode, xlab = "time", main = main,
                 ylim = ylim, cexlab = cexlab, legend = legend, groups.vector = groups.vector,
                 ...)
    }
    c.algo.used <- NULL
    cut <- 1
  }
  else {
    print("warning: NULL data. No visualization possible")
    c.algo.used <- NULL
    cut <- NULL
  }
  OUTPUT <- list(cut, c.algo.used, groups)
  names(OUTPUT) <- c("cut", "cluster.algorithm.used", "groups")
  OUTPUT
}

################################
# Clustering Tools
################################

myHclust = function(testCounts, DIST = TRUE){
  #Uses 1-correlation, as the metric for HC
  
  if(DIST == TRUE){
    dcorrel = dist(testCounts)
  } else {
    # Uses correlation as distance
    dcorrel = as.dist(matrix(rep(1, nrow(testCounts)^2), nrow(testCounts), nrow(testCounts)) - cor(t(testCounts), use = "pairwise.complete.obs"))
  }
  
  
  # Get best number of clusters
  asw = numeric(8)
  for (k in 2:8){
    asw[[k]] = pam(dcorrel, k)$silinfo$avg.width
  }
  k.best = which.max(asw)
  Gclass = pam(x = dcorrel,k = k.best, cluster.only=TRUE)
  #Gclust = hclust(dcorrel)
  #Gclass = cutree(Gclust,k = k.best) 
  # Get lists of genes
  Glist = list()
  Nclust = unique(Gclass)
  for(i in Nclust){
    Glist[[paste("Clust",as.character(i),sep="_")]] = names(Gclass)[Gclass==i]
  }
  return(Glist)
  
}

testSignatureDynamics = function(testCounts,testTargets,SignatureDictionary,folderPATH){
  
  LM_globalFMAT = getFunctionMatrix_v2(eMat = testCounts,eTargets = testTargets,Signature_Dict = SignatureDictionary,weighted = TRUE)
  LM_globalFMAT = LM_globalFMAT$ESmat
  
  # Filter times with low number of samples
  
  completeD = names(table(testTargets$time))[table(testTargets$time) >=8]
  LM_globalTargets = filter(testTargets, time %in% completeD)
  LM_globalFMAT = LM_globalFMAT[,as.character(LM_globalTargets$sample_mask_id)]
  
  # Get results from T-test
  masigClusters = getSignatures(cMAT = LM_globalFMAT,targetMAT = LM_globalTargets,Features = NULL,time = sort(unique(LM_globalTargets$time)),g.pval = 0.05,s.pval=0.05,model_type = "glm")
  # Create matrix of results
  resultData = masigClusters$RESULTS
  
  pdf(file = paste(folderPATH,"Dynamics.pdf",sep = "_"),height = 10,width = 12)
  for(k in rownames(LM_globalFMAT)){
    plotES(ES_list = LM_globalFMAT[k,],sigTargets = LM_globalTargets,TRSH = 0.05,main_title = k)
  }
  dev.off()
  
  pdf(file = paste(folderPATH,"LMresults.pdf",sep = "_"),height = 8,width = 8)
  SummaryTables = summaryTestTable(resultData)
  dev.off()
  
  return(list("ESmat"=LM_globalFMAT,"GLM_Res"=resultData))
  
}

plot_CompareDynamics = function(MedianResps, main){
  
  plotDF = c()
  for(TG in names(MedianResps)){
    plotDF = rbind(plotDF, cbind(melt(MedianResps[[TG]]),TestGroup = TG))
  }
  colnames(plotDF) = c("ClusterName","Time","MedianES","TestGroup")
  plotDF = cbind(plotDF,Tclass = paste(plotDF$ClusterName,plotDF$TestGroup,sep = "_"))
  
  plotDF$Time =  as.numeric(gsub("time_","",plotDF$Time))
  plotDF$ClusterName = factor(plotDF$ClusterName)
  plotDF$MedianES = as.numeric(plotDF$MedianES)
  plotDF$TestGroup = factor(plotDF$TestGroup)
  plotDF$Tclass = factor(plotDF$Tclass)
  GEX_plot <- ggplot(plotDF, aes(x = plotDF[,"Time"] * -1, y = plotDF[,"MedianES"], group = Tclass, color = TestGroup)) + 
    scale_x_continuous(name="Time",breaks=sort(unique(plotDF[,"Time"] * -1))) + ylab("Median Expression") + labs(title = main) +
    theme(axis.title.x = element_text(size = 14,face="bold"), axis.title.y = element_text(size = 14,face="bold"),axis.text.x  = element_text(size=14),axis.text.y  = element_text(size=14)) +
    geom_point() + geom_line()
  print(GEX_plot)
  #return(GEX_plot)
}

##########################################

testDynamics_Exhaustive = function(TestCounts_list, TestAnnotations_list, globalResultFolder, clusterList){
  
  MedianValues = list()
  
  TestNames = names(TestCounts_list)
  
  Tset_Results = list()
  
  for(Tset in TestNames){
    
    print(paste("Testing",Tset))
    
    countT = TestCounts_list[[Tset]]
    annT = TestAnnotations_list[[Tset]]
    
    
    folderPATH = paste(globalResultFolder,Tset,"/",sep = "")
    system(paste("mkdir ",folderPATH))
    
    #1st part, Make an enrichment matrix for each testset
    #Compare if there's is a difference between cases and controls, using Enrichment Scores
    print(paste("Results at", folderPATH))
    dyn_RES = testSignatureDynamics(testCounts = countT,testTargets = annT,SignatureDictionary = clusterList, folderPATH = folderPATH)
    Tset_Results[[Tset]] = dyn_RES
    
    # Then, I calculated the median of the enrichment score of all the cases from the test set at all times
    testcases = filter(annT, outcome==1, sample_mask_id %in% colnames(dyn_RES$ESmat))
    medianResponse = getGroupMedian(groupcolumn = "time",stargets = testcases,sexprmat = dyn_RES$ESmat[,as.character(testcases$sample_mask_id)])
    
    MedianValues[[Tset]] = medianResponse
    
  }
  
  #Plot all medians together
  pdf(paste(globalResultFolder,"MedianComparison.pdf",sep=""),width = 12, height = 10)
  for(K in names(clusterList)){
    KmedianResps = list()
    for(Tset in TestNames){
      KmedianResps[[Tset]] = MedianValues[[Tset]][K,,drop=F]
    }
    plot_CompareDynamics(MedianResps = KmedianResps,main = K)
  }
  dev.off()
  
  #Plot all combinations of groups
  
  AgeCombinations = combn(TestNames,2)
  
  apply(AgeCombinations,2,function(x){
    Agecomb = paste(x,collapse = "_")
    G1 = x[1]
    G2 = x[2]
    filename = paste(globalResultFolder,Agecomb,".pdf",sep="")
    MedianResp = MedianValues[x]
    
    pdf(filename,width = 12, height = 10)
    
    for(K in names(clusterList)){
      MR = list()
      MR[[G1]] = MedianResp[[G1]][K,,drop=F]
      MR[[G2]] = MedianResp[[G2]][K,,drop=F]
      plot_CompareDynamics(MedianResps = MR,main = K)
    }
    
    dev.off()
  })
  
  return(Tset_Results)
}

write_clusters = function(globalResultFolder, clusterList){
  folderPATH = paste(globalResultFolder,"ClusterList","/",sep = "")
  system(paste("mkdir ",folderPATH))
  for(k in names(clusterList)){
    Elist = (na.omit(data.frame(mapIds(illuminaHumanv4.db, keys= clusterList[[k]], keytype = "SYMBOL", column = "ENTREZID",multiVals = "first"),
                                stringsAsFactors = F)))[,1]
    write.table(clusterList[[k]], file = paste(folderPATH,k,".txt",sep=""), sep = "\t",row.names = F,col.names = F,quote = F)
    write.table(Elist, file = paste(folderPATH,k,"_ENTREZ",".txt",sep=""), sep = "\t",row.names = F,col.names = F,quote = F)
  }
}

GetDangerMatrix = function(DangerZone, testDynamics_Results, clusterList){
  
  DangerMat = matrix(0,nrow = length(names(testDynamics_Results)), ncol = length(names(clusterList)))
  colnames(DangerMat) = names(clusterList)
  rownames(DangerMat) = names(testDynamics_Results)
  for(Tset in names(testDynamics_Results)){
    
    GLM_Res = testDynamics_Results[[Tset]]$GLM_Res
    DangerPoints = names(GLM_Res)[names(GLM_Res) %in% DangerZone]
    GLM_Res = GLM_Res[DangerPoints]
    
    for(DP in DangerPoints){
      spcGLM_Resp = GLM_Res[[DP]]
      FClust = spcGLM_Resp$Hallmark
      Vals = ifelse(spcGLM_Resp$adj_P_val<0.05,1,0)
      DangerMat[Tset,FClust] = DangerMat[Tset,FClust] + Vals
    }
  }
  
  return(DangerMat)
  
}



