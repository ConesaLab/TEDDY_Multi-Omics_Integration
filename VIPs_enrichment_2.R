##### PIANO ANALYSIS FOR VIPS                                    ###
#### Script used in revision for Genome Biology 30 october 2020
###################################################################

library(piano)
library(RColorBrewer)
library("d3heatmap")
library("tools")
library("metaseqR")
library("ggplot2")

# VIP data per time-point
VIPs <- read.delim("Data/VIP3D2.txt", as.is = TRUE, row.names = 1) ; head (VIPs)
VIPgs <- read.delim("Data/VIP3D1.txt", as.is = TRUE, row.names = 1) ; head (VIPgs)
VIPgs_add <- cbind(VIPgs[,2], VIPgs[,2]+ VIPgs[,3],VIPgs[,2]+VIPgs[,3],+VIPgs[,4])
VIPs_list <- list()  # I make ONE list with VIPs per time point (12 to 0), and combined for 3 components
 for (i in 1:5) {
   vip <- VIPs[,(i+1)]
   names(vip) <- VIPs[,1] 
   vip <- sort(vip, decreasing = TRUE)
   VIPs_list[[i]] <- vip
 }

for (i in 1:3) {
  vip <- VIPgs_add[,(i)]
  names(vip) <- VIPgs[,1] 
  vip <- sort(vip, decreasing = TRUE)
  VIPs_list[[i+5]] <- vip
}
#############
# Functions #
#############
# Calculates Enrichment with PIANO, non directional, with mean and median and return lowest pvalue
min.pval <- function (VIPs_list,annotations, pval) {
 med <- lapply(VIPs_list, function (x) {runGSA(x, gsc=annotations, nPerm=500, geneSetStat = "median")})
 sum  <- lapply(VIPs_list, function (x) {runGSA(x, gsc=annotations, nPerm=500, geneSetStat = "sum")})
 pvals <- list()
 for (i in 1:length(VIPs_list)) {
   a<- names(med[[i]]$gsc)
   b <- as.numeric(med[[i]]$pNonDirectional[,1])
   c <- as.numeric(sum[[i]]$pNonDirectional[,1])
   p <- apply(cbind(b,c), 1, min)
   d <- data.frame (term = a,med = b, sum = c, min.pval= p)
   d <- d[d$min.pval < pval,]
   pvals[[i]] <- d
  }
 pvals 
}

# Combines pvalues of differnet time points with the Fisher Method
pvals.table <- function (pvals, from , to , filter = 0.05, times = 1) {
  # pvalues in different times are put in the same data frame
  all.p <-merge(pvals[[from ]], pvals[[(from+1)]], by.x = 1, by.y = 1, all = TRUE)
  colnames(all.p)[2:3] <- names(pvals)[1:2]
  for (i in (from+2):to) {
    all.p <- merge(all.p, pvals[[i]], by.x = 1, by.y = 1, all = TRUE)
    colnames(all.p)[2:(i+1)] <- names(pvals)[1:i]
  }
  rownames(all.p) <- all.p[,1] ; all.p <- as.matrix( all.p[,-1])
  # we calculate now a combined pvalue following the fisher method. All possible combinatios of time points
  fisher.combinations <- list()
  for (i in 2:5) {
    fisher.combinations <- c(fisher.combinations, as.list(as.data.frame(combn(c(1:5), i, FUN = NULL))))
  } # these are all possible combinations are indexes
  yp <- sqrt(length(fisher.combinations))
  fisher.p <- fisher.padjust<- matrix(0, nrow = nrow(all.p), ncol = 1)
  rownames(fisher.p) <- rownames(fisher.padjust) <- rownames(all.p)
  for (i in fisher.combinations) { # apply fisher method to all combinations
    fm <- fisher.method(all.p[,i], p.corr = "BH")
    comby <- paste(colnames(all.p[,i]), collapse = "_")
    colnames(fm) <- paste(colnames(fm), comby, sep = "_")
    fisher.p <- cbind(fisher.p , fm[,3]); colnames(fisher.p)[ncol(fisher.p )] <- comby # p.value
    fisher.padjust <- cbind(fisher.padjust, fm[,4]/yp); colnames(fisher.padjust)[ncol(fisher.p )] <- comby #adjust p.value
  }
 
  # Now calculating the lowest combined pvalues of all possible time combinations
  combined.p <- as.matrix(fisher.p[,-1])
  combined.p <- cbind(combined.p, combined.p.min = rowMins(combined.p, na.rm = TRUE))
  combined.padjust <- as.matrix(fisher.padjust[,-1])
  combined.padjust <- cbind(combined.padjust, combined.padjust.min = rowMins(combined.padjust, na.rm = TRUE))
  single_gsa <- cbind(all.p, lowest.p = rowMins(all.p, na.rm = TRUE))
  single_gsa <- single_gsa[rownames(combined.p ),] # having same order of rows
  # Identify the time combination with the lowest combined pvalue
  CC <- apply(combined.p, 1,function (x) which(x[1:(length(x)-1)] == x[length(x)])) ; CC <- sapply(CC,function(x) x[1])
  times.lowest.p <- colnames(combined.p)[unlist(CC)]
  LL <- apply(combined.padjust, 1,function (x) which(x[1:(length(x)-1)] == x[length(x)])) ; LL <- sapply(LL,function(x) x[1])
  times.lowest.padj <- colnames(combined.padjust)[unlist(LL)]
  MM <- apply(single_gsa, 1,function (x) which(x[1:(length(x)-1)] == x[length(x)])) ; MM <- sapply(MM,function(x) x[1])
  time.lowest.single_gsa <- colnames(single_gsa)[unlist(MM)]
 
  # find the min pvalue for single time points
  # gather results
  all.results <- data.frame(Single.p = all.p, 
                            time.lowest.single_gsa  = time.lowest.single_gsa,
                            FisherMethod.p = rowMins(combined.p, na.rm = TRUE),
                            times.lowest.p = times.lowest.p,
                            FisherMethod.padjust = rowMins(combined.padjust, na.rm = TRUE), 
                            times.lowest.padj = times.lowest.padj)
  rownames(all.results) <- rownames(combined.p)
  
  all.results <- all.results[order(all.results$FisherMethod.padjust, decreasing = FALSE),] # order
  
  #final result
  all.results
}

# Returns the expression values for all the genes in one particular pathway
check.pathway <-function (matrix= data.mean, pathway= "mTOR", annotation = TEDDY_geneSets) {
  paths <- lapply(TEDDY_geneSets, function (x) x[grep(names(x), ignore.case = TRUE, pattern = pathway)])
  paths <<-  paths[unlist(lapply (paths, function (x) length(x) > 0))]
  path.expression <- list()
  par(mfrow=c(4,4))
  h = 0
  for ( i in 1:length(paths)) {
    for (j in 1:length(paths[[i]])) {
      h = h+1
      genes.path <- paths[[i]][[j]]
      #print(genes.path)
      path.expression[[h]] <- data.mean[genes.path,]
      print (names(paths[[i]])[j])
      heatmap.2(as.matrix(path.expression[[h]]), col = rev(brewer.pal(9,"RdBu")),
                Colv=FALSE, key = FALSE, Rowv = TRUE, dendrogram="none", density.info="none", trace="none",
                main = names(paths[[i]])[j])
    }
  }
}

# format data for ggplot heatmap
format2gplot <- function (x) {
  x <- x[nrow(x):1,]
  term = rep(rownames(x),5)
  pval = c(x[,1],x[,2],x[,3],x[,4],x[,5])
  direction = sign(pval)
  time = rep(colnames(x), each = nrow(x))
  data.gplot = data.frame(term = factor(term, levels = rownames(x)), 
                          direction = direction,
                          pval = pval,
                          time = time)
  data.gplot$time <- factor(data.gplot$time,levels=colnames(x))
  return (data.gplot)
}

# Converts string to sentence case
Caps <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

# This function applies both min.pval and pvals.table to VIPs_list with the given gene set.
wrapper_time_course_pathway_enrichment <- function (sets, exclude, VIPs_list, name, fdr = 0.05,from = 1, to = 5) {
  sets.sel <- sets[-c(exclude)]
  GSC = c()
  for(GS in names(sets.sel)){
    GS_list = cbind(sets.sel[[GS]],GS)
    GSC = rbind(GSC,GS_list)
  }
  GSCv = data.frame(GSC,stringsAsFactors = F)
  GSC = loadGSC(GSCv) # This is the Piano object of the gene set collection
  min.pval <- min.pval(VIPs_list, GSC, 1)
  min.pval2<-lapply(min.pval, function (x) {x[,c(1,4)]})
  names(min.pval2) <- c("M-12", "M-9", "M-6", "M-3", "M-0", "First", "Second", "Third")
  save(min.pval2, file = paste("Piano", name, "selected_pathways.RData", sep = "_"))
  print("Piano results saved!")
  print("Calculating Fisher Method")
  min.pval.table <- pvals.table(min.pval2, from = from, to = to) # obtain time combined pvalues for 5 time points
  print("Fisher Method done!")
  final <- min.pval.table[min.pval.table$FisherMethod.padjust < fdr,]  # select by padjust < -.05
  result <- list(final, GSCv) ; names(result) <- paste(c("enriched_results", "GSCv"), name, sep = "_")
  result
}

####################################
# Piano analysis for KEGG database #
####################################

load("Data/Piano/KEGG_REST_Annotation.ro") # This was downloaded by KEGG database
exclude.kegg <- c(19,37,43:46,50:52,62:64,73,75,80,83,88,89,94:96,99,101,108,113:115,131,133,136,138:140, 
                  143:144,150:156,158,162,164:165,167,
                  174,182,197:210,216:222,224,231,238,240:242,246,249:317)  # removing non-applicable pathways
KEGG_result <- wrapper_time_course_pathway_enrichment(KEGG_REST_sets, exclude.kegg, VIPs_list, name = "KEGG", fdr = 0.05)
write.table(KEGG_result[[1]], file = "KEGG_result.txt", row.names = T, col.names = T, sep = "\t", quote = F )

##################################################
#####  Piano analysis for all Genesets           #
load("Data/Piano/GlobalData/TEDDY_geneSets.ro")  #
##################################################

####### Biological process
BP <- TEDDY_geneSets[["MSIGDB_GO_BIOLPROC"]]
excluding_words <- c("Cancer", "Disease", "Alzheimer", "Neuro", "Parkinson", "Malaria", "Cangas", "Bacteria", "Fung", "Alzheimer", "HIV")
exclude.BP <- unlist(sapply(excluding_words, function (x) grep(x = names(BP), pattern = x, ignore.case = TRUE)))
BP_result <- wrapper_time_course_pathway_enrichment(BP, exclude = exclude.BP , VIPs_list, name = "BP", fdr = 0.05)
write.table(BP_result[[1]], file = "BP_result.txt", row.names = T, col.names = T, sep = "\t", quote = F)

####### Molecular function
MF <- TEDDY_geneSets[["MSIGDB_GO_MOLFUNC"]]
exclude.MF <- unlist(sapply(excluding_words, function (x) grep(x = names(MF), pattern = x, ignore.case = TRUE)))
MF_result <- wrapper_time_course_pathway_enrichment(MF, exclude = exclude.MF , VIPs_list, name = "MF")
write.table(MF_result[[1]], file = "MF_result.txt", row.names = T, col.names = T, sep = "\t", quote = F)

####### Biocarta
CART <- TEDDY_geneSets[["PID_BIOCARTA"]]
exclude.CART <- c(1,5,20,101,137,142,216,220,226,228,231,248)
CART_result <- wrapper_time_course_pathway_enrichment(CART, exclude = exclude.CART , VIPs_list, name = "CART", fdr = 0.1)
write.table(CART_result[[1]], file = "CART_result.txt", row.names = T, col.names = T, sep = "\t", quote = F)

#######  Reactome
REACTOME <- TEDDY_geneSets[["REACTOME"]]
excluding_words <- c("Cancer", "Disease", "Alzheimer", "Neuro", "Parkinson", "Malaria", "Cangas", "Bacteria", "Fung", "Alzheimer", "HIV")
exclude.REACTOME <- unlist(sapply(excluding_words, function (x) grep(x = names(REACTOME), pattern = x, ignore.case = TRUE)))
REACTOME_result <- wrapper_time_course_pathway_enrichment(REACTOME, exclude = exclude.REACTOME , VIPs_list, name = "REACTOME", fdr = 0.05)
write.table(REACTOME_result[[1]], file = "REACTOME_result.txt", row.names = T, col.names = T, sep = "\t", quote = F)

#######  HMARKS
HMARKS <- TEDDY_geneSets[["MSIGDB_HMARKS"]]
exclude.HMARKS <- c(37,39,42)
HMARKS_result <- wrapper_time_course_pathway_enrichment(HMARKS, exclude = exclude.HMARKS , VIPs_list, name = "HMARKS")
write.table(HMARKS_result[[1]], file = "HMARKS_result.txt", row.names = T, col.names = T, sep = "\t", quote = F)

#######  Recover gene expression data
# We identify if for each significant pathway the gene set tends to be up or down regualted

# Gene expression data
data <- read.csv ('GeneMeans.csv', as.is = T, row.names = 1) # for each gene, mean, median and VIP of cases and controls
data.mean <- data[,c(1,4,7,10,13)]
data.med <- data[,c(2,5,8,11,14)]
data.vip <- data[,c(3,6,9,12,15)]
data.mean.vip <- data.mean * data.vip # genes are weighted by thier VIP score
data.med.vip <- data.med * data.vip # genes are weighted by thier VIP score

## Final gene set
GSC.all <- rbind (KEGG_result$GSCv_KEGG,
                  BP_result$GSCv_BP,
                  MF_result$GSCv_MF,
                  REACTOME_result$GSCv_REACTOME,
                  CART_result$GSCv_CART,
                  HMARKS_result$GSCv_HMARKS) # all annotations

final_gene_set <- rbind(KEGG_result[[1]][,1:5], 
                        BP_result[[1]][,1:5], 
                        MF_result[[1]][,1:5], 
                        REACTOME_result[[1]][,1:5],
                        CART_result[[1]][,1:5],
                        HMARKS_result[[1]][,1:5]
                        ) # all significant gene sets
final_gene_set.pos <- 1-final_gene_set
final_path_vip.mean  <- final_path_vip.med <-  final_gene_set.pos
for (i in 1: nrow(final_path_vip.mean)) {   # this calculates the mean of the gene set
  path <- rownames(final_gene_set)[i]
  genes <- GSC.all[GSC.all[,2] == path,1]
  path.data1 <-data.mean.vip[genes,]
  path.data2 <-data.med.vip[genes,]
  final_path_vip.mean[i,]  <- apply(path.data1, 2, mean, na.rm = TRUE)
  final_path_vip.med[i,]  <- apply(path.data2, 2, median, na.rm = TRUE)
}

final_path.pos_mean <- final_gene_set.pos *sign(final_path_vip.mean) # we combine the sign of the geneset expression and its pvalue
final_path.pos_med <- final_gene_set.pos *sign(final_path_vip.med)

# How to look for a specific pathway
check.pathway(matrix= data.mean, pathway= "Spliceosome", annotation = TEDDY_geneSets)


###########################
####### Creating Heatmaps #
###########################
save(final_path.pos_mean, file = "heatmap_data_mean_revision.Rdata")
save(final_path.pos_med, file = "heatmap_data_median_revision.Rdata")
write.table(final_path.pos_mean, file = "heatmap_data_mean_revision.txt", row.names = T, col.names = T, sep = "\t", quote = F )
write.table(final_path.pos_med, file = "heatmap_data_mean_revision.txt", row.names = T, col.names = T, sep = "\t", quote = F )
final_path_clean <- read.delim ("heatmap_data_Fisher.txt", row.names = 2, header = T, as.is = T ) # in Excel we edited names and add a pathway type
head(final_path_clean)
type.path <- final_path_clean$type
final_path_clean.data <- final_path_clean[,c(3:7)]

##### Metabolic Pathways
## Using d3heatmap

data.heat = as.matrix(final_path_clean.data [type.path == "m", ])
data.heat.ordered <- data.heat[c(4,14,11,7,10,6,12,13,15,18,3,1,5,19,2,9,8,16,17),]
d3heatmap(data.heat.ordered, colors = rev(brewer.pal(9,"RdBu")),Colv = FALSE,
          dendrogram = 'none', Rowv = FALSE, cexRow = 0.9, cexCol = 0.9, 
          yaxis_width=300)


# Using ggplot
data.heat.ordered.gplot <- format2gplot(data.heat.ordered)
max(sapply(as.character(data.heat.ordered.gplot[,1]),nchar))
heatmap_metabolic_pathways = ggplot(data.heat.ordered.gplot , aes(time, term)) + 
  geom_tile(aes(fill = pval), colour = "white") + 
  scale_fill_gradientn(colours = rev(brewer.pal(9,"RdBu")), values = scales::rescale(range(data.heat.ordered.gplot$pval))) + 
  scale_x_discrete("", expand = c(0, 0)) + 
  scale_y_discrete("", expand = c(0, 0)) + 
  geom_text(aes(label = ifelse(abs(pval)>0.88,"*","")), size=10) + # To put asterisk in significant
  theme_grey(base_size = 28) +
  theme(legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text.x = element_text(angle = 330, hjust = 0,size = 16),
        axis.text.y = element_text(size = 15)) 
print(heatmap_metabolic_pathways)


##### Signalling Pathways
data.heat = as.matrix(final_path_clean.data [type.path == "s", ])
row.order <- d3heatmap(data.heat, colors = rev(brewer.pal(9,"RdBu")),Colv = FALSE,
          dendrogram = 'none',Rowv = TRUE,cexRow = 0.9, cexCol = 0.9, 
          yaxis_width=300)$x$matrix$rows
new.order <- row.order[c(1:7,15:23,24:28,8:12,14,13)]

## Using d3heatmap
d3heatmap(data.heat[new.order,], colors = rev(brewer.pal(9,"RdBu")),Colv = FALSE,
                       dendrogram = 'none',cexRow = 0.9, cexCol = 0.9, 
                       yaxis_width=300)

# Using ggplot
data.heat.ordered.gplot <- format2gplot(data.heat[new.order,])
data.heat.ordered.gplot[data.heat.ordered.gplot$term == "Chemokine gene expression in HMC1",]
levels(data.heat.ordered.gplot$term)[levels(data.heat.ordered.gplot$term)=="Chemokine gene expression in HMC1"] <- "Chemokine gene expression in MHC1"
data.heat.ordered.gplot[data.heat.ordered.gplot$term == "Chemokine gene expression in MHC1",]
heatmap_signaling_pathways = ggplot(data.heat.ordered.gplot , aes(time, term)) + 
  geom_tile(aes(fill = pval), colour = "white") + 
  scale_fill_gradientn(colours = rev(brewer.pal(9,"RdBu")), 
                       limits = c(-1,1),
                       labels=c(-1,-0.5,0,0.5,1),
                       values = scales::rescale(range(data.heat.ordered.gplot$pval))#,
                       #guide = guide_legend(label.vjust = 40
                       ) + 
  scale_x_discrete("", expand = c(0, 0)) + 
  scale_y_discrete("", expand = c(0, 0)) + 
  geom_text(aes(label = ifelse(abs(pval)>0.9,"*","")), size=10) + # To put asterisk in significant
  theme_grey(base_size = 28) +
  theme(legend.position = "bottom",
        legend.title=element_text(size=0),
        legend.text=element_text(size=18),
        axis.ticks = element_blank(), 
        axis.text.x = element_text(angle = 330, hjust = 0,size = 16),
        axis.text.y = element_text(size = 16)) +
        guides(fill = guide_colourbar(barwidth = 15, barheight = 3,nbin=50,
                                      #border=element_line(color='black'),
                                      label.vjust = 20
                                      ))
print(heatmap_signaling_pathways)

# Heatmaps for metabolites, just to have everything in the same format
#######################################################################

pianoDFALLselectionSorted2 <- read.delim("Data/pianoDFALLselectionSorted2.txt", as.is = TRUE, sep = " ")
pianoDFALLselectionSorted2[,1] <- sapply(pianoDFALLselectionSorted2[,1], Caps)

direction_classes = c("Distinct-directional (up)","Mixed-directional (up)","Non-directional","Mixed-directional (dn)","Distinct-directional (dn)")
directioncolors = c("red","red", "white","blue1","blue1")
LEGEND_IX = which(direction_classes %in% unique(pianoDFALLselectionSorted2$DIRECTION))

direction_classes = direction_classes[LEGEND_IX]
pianoDFALLselectionSorted2$DIRECTION = factor(pianoDFALLselectionSorted2$DIRECTION, levels = direction_classes)

#Assign the same colors with the neww palette
brewer.pal(9,"RdBu")
directioncolors = c("#B2182B","#B2182B", "#F7F7F7","#2166AC","#2166AC")
directioncolors = directioncolors[LEGEND_IX]
colorends = 1:(length(direction_classes)*2)
colorends[(1:(length(direction_classes)) * 2) - 1] = "white"
colorends[(1:(length(direction_classes)) * 2) ] = directioncolors

# Reescaling values so we can plot different categories
Nclasses = sort((unique(100 * (as.numeric(pianoDFALLselectionSorted2$DIRECTION) ))))
scalerange <- range(pianoDFALLselectionSorted2$PVALcomp)
gradientends <- scalerange + rep(Nclasses, each=2)

# This sorting I liked
pianoDFALLselectionSorted2 <- pianoDFALLselectionSorted2[order(pianoDFALLselectionSorted2$directionnumbers, pianoDFALLselectionSorted2$TIME,pianoDFALLselectionSorted2$PVAL),]
pianoDFALLselectionSorted2$MetSET <- factor(pianoDFALLselectionSorted2$MetSET, levels = rev(unique(as.character(pianoDFALLselectionSorted2$MetSET))))
vectitoseleccion<-as.vector(unique(pianoDFALLselectionSorted2$MetSET))
vectitoseleccion<-rev(c("Fatty Acid","Phosphatidylethanolamine","Lysophosphatidylethanolamine",
                    "Cholesterol","Aminoacid Precursor","Glucosylceramide","Carbohydrate",              
                    "Monoglyceride","Sphingomyelin","Ceramide","Phosphatidylcholine",
                    "Triglyceride","Lysophosphatidylcholine") )

pianoDFALLselectionSorted2<-pianoDFALLselectionSorted2[pianoDFALLselectionSorted2$MetSET %in% vectitoseleccion,]
pianoDFALLselectionSorted2$MetSET <- factor(pianoDFALLselectionSorted2$MetSET,levels=vectitoseleccion)
pianoDFALLselectionSorted2$MetSET

#set vector of levels
mylevels <-  c("m12metabs","m9metabs","m6metabs","m3metabs","m0metabs")
#reorder factors
pianoDFALLselectionSorted2$TIME <- factor(pianoDFALLselectionSorted2$TIME,levels=mylevels)
pianoDFALLselectionSorted2$PVALcomp
pianoDFALLselectionSorted2$PVALresc
pianoDFALLselectionSorted2$PVAL
pianoDFALLselectionSorted2[pianoDFALLselectionSorted2$MetSET == "Glucosylceramide",]

pselectionSortedselected = ggplot(pianoDFALLselectionSorted2, aes(TIME, MetSET)) + 
  geom_tile(aes(fill = PVALresc), colour = "white") + 
  #scale_fill_gradientn(colours = brewer.pal(9,"RdBu"), values = scales::rescale(gradientends)) + 
  scale_fill_gradientn(colours = colorends, values = scales::rescale(gradientends)) + 
  scale_x_discrete("", expand = c(0, 0), labels=c("m12metabs" = "M.12", "m9metabs" = "M.9","m6metabs" = "M.6",                                                "m3metabs" = "M.3", "m0metabs" = "M.0")) + 
  scale_y_discrete("", expand = c(0, 0)) + 
  geom_text(aes(label = ifelse(PVAL<0.2,"*","")), size=10) + # To put asterisk in significant
  theme_grey(base_size = 28) +
  theme(legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text.x = element_text(angle = 330, hjust = 0,size = 16),
        axis.text.y = element_text(size = 16)) 
print(pselectionSortedselected)

    
#### All plots in one figure for paper
######################################

library(grid)
pdf(file="Heatmaps_Figure3.pdf",height = 13,width = 20)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 40, ncol = 40)))
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(heatmap_metabolic_pathways, vp = define_region(row = 1:22, col = 1:19))   # Span over two columns
print(pselectionSortedselected, vp = define_region(row = 23:40, col = 2:19))
print(heatmap_signaling_pathways, vp = define_region(row = 1:39, col = 21:40))
dev.off()

######### Supplementary table SD13

rownames(final_path_clean)
final_path_clean[,1]

KEGG_result <- read.delim("KEGG_result.txt", row.names = 1, header = T, as.is = TRUE)
BP_result <- read.delim("BP_result.txt", row.names = 1, header = T, as.is = TRUE)
MF_result <- read.delim("MF_result.txt", row.names = 1, header = T, as.is = TRUE)
CART_result2 <- read.delim("CART_result2.txt", row.names = 1, header = T, as.is = TRUE)
HMARKS_result <- read.delim("HMARKS_result.txt", row.names = 1, header = T, as.is = TRUE)
REACTOME_result2 <- read.delim("REACTOME_result2.txt", row.names = 1, header = T, as.is = TRUE)


table <- rbind(KEGG_result[is.element(rownames(KEGG_result), final_path_clean[,1]),c(1:5,8,7,9)],
               BP_result[is.element(rownames(BP_result), final_path_clean[,1]),c(1:5,8,7,9)],
               MF_result[is.element(rownames(MF_result), final_path_clean[,1]),c(1:5,8,7,9)],
               CART_result2[is.element(rownames(CART_result2), final_path_clean[,1]),c(1:5,8,7,9)],
               HMARKS_result[is.element(rownames(HMARKS_result), final_path_clean[,1]),c(1:5,8,7,9)],
               REACTOME_result2[is.element(rownames(REACTOME_result2), final_path_clean[,1]),c(1:5,8,7,9)]) 
final_path_clean2 = data.frame(Pathway = rownames(final_path_clean), final_path_clean)
table2 <- merge(table, final_path_clean2, by.x = 0, by.y = 2 )
table2 <- table2[,c(1,10,2:9)]
write.table(table2, file = "SD13_Pathway_Enrichment.txt", row.names = T, col.names = T, sep = "\t", quote = F )


grep(names(REACTOME), pattern = "glucagon", value = TRUE)
grep(rownames(REACTOME_result[[1]]), pattern = "glucagon", value = TRUE)
grep(rownames(REACTOME_result2), pattern = "glucagon", value = TRUE)
head(REACTOME_result[[1]])


#dev.off()