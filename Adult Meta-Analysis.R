## Meta-analysis on Adult HR/LR Data

#This script contains code running the meta-analysis of 5 adult HR/LR datasets
#To run the meta-analysis, t-test of phenotype effect on gene expression is read-in for each dataset. Cohen's d is then extracted from t-test values and joined together for the meta-analysis.
#Note: one dataset, the RNA-seq F29 dataset, is the only dataset not re-analyzed prior to the meta-analysis as original data files were unavailable. Instead, Cohen's d is extracted from the p-values obtained through a prior analysis.

### Loading necessary libraries

library(plyr)
library(metafor)
library(compute.es)
library(bmp)

#This package is necessary to run the tes function

source("https://bioconductor.org/biocLite.R")
biocLite("multtest")
library(multtest)


## Reading in data

# Reading in test statistic data

RNAseq_F29 <- read.csv("data/RNAseq F29/LRHR_rgsc34_th2out_g1tx_diff_AdLvsAdH_out_gene_exp.csv")
RNAseq_F37 <- read.csv("output/Blandino RNAseq F37/BlandinoTT.csv")
RNAseq_F43 <- read.csv("output/Aydin RNAseq F43/F43_AdultTT.csv")
Affy_F4    <- read.csv("output/Stead Affy F4/Stead_AdultF4TT.csv")
NimbleGen_F34 <- read.csv("output/Alabama New Colony Development/newdevelopmentdataAdultTT.csv")

# Reading in cell-type data

cellRNAseq_F37 <- read.csv("output/Blandino RNAseq F37/BlandinoCellTypeTT_forMeta.csv")
cellRNAseq_F43 <- read.csv("output/Aydin RNAseq F43/VehCellTypeTT_forMeta.csv")
cellAffy_F4 <- read.csv("output/Stead Affy F4/Stead_AdultF4CellTypeTT_forMeta.csv")
cellNimbleGen_F34 <- read.csv("output/Alabama New Colony Development/NewColonyDevelopmentCellTypeTT_forMeta.csv")


#setting directories

outDir <- "output/MetaAnalysis Adult/"
outForestPlots <- "output/MetaAnalysis Adult/ForestPlots_topGenes_Adult/"
outGOIPlots <- "output/MetaAnalysis Adult/Forest Plots for GOI Adult/"
outCellPlots <- "output/MetaAnalysis Adult/ForestPlots_CellType_Adult/"


## Extracting Effect Size
### The following code converts t-test values into effect sizes d (Cohen's D) and g (Hedge's g) and extracts the corresponding variances for each of the five adult datasets. 


# RNAseq_F29

colnames(RNAseq_F29)

#need to omit data with no value (p-value=1)
rnaseq<-subset(RNAseq_F29, RNAseq_F29$p_value!=1)
dim(rnaseq)

sum(duplicated(rnaseq$gene))
#lots of entries have associated gene IDs but no gene symbol

#removing genes with no gene symbol

rnaseq <- subset(rnaseq, !(rnaseq$gene=="-"))

###Removing extreme values (genes where transcript info was only present for one group HR or LR)

rnaseq_nobad<-subset(rnaseq, rnaseq$test_stat > -1.79769e+308 & rnaseq$test_stat < 1.79769e+308)


sum(duplicated(rnaseq_nobad$gene))
#still lots of duplicated genes. Many of these seem to have very different expression values.

#Averaging across duplicate genes
colnames(rnaseq_nobad)

rnaseq_avg<-tapply(rnaseq_nobad[,11], as.character(rnaseq_nobad[,3]), function(y) mean(y))


#extracting the test statistic averaged by gene
rnaseqTT<-cbind.data.frame(rnaseq_avg)

CohDRNAseq<-matrix(0, length(rnaseqTT[,1]), 4)

for(i in 1:length(rnaseqTT[,1])){
  CohDRNAseq[i, 1]<-tes(rnaseqTT[i,], 2,2)$d
  CohDRNAseq[i, 2]<-tes(rnaseqTT[i,], 2,2)$'var.d'
  CohDRNAseq[i, 3]<-tes(rnaseqTT[i,], 2,2)$g
  CohDRNAseq[i, 4]<-tes(rnaseqTT[i,], 2,2)$'var.g'
}

colnames(CohDRNAseq)[1:4]<-c("d RNAseq F29", "var d RNAseq F29", "g RNAseq F29", "var g RNAseq F29")


#adding in gene symbol information
CohDRNAseq<-cbind.data.frame(CohDRNAseq, Symbol=row.names(rnaseqTT))


#Flipping direction of effect so it matches other datasets (HR positive LR negative)
CohDRNAseq[,1] <- (-1*CohDRNAseq[[1]])



# RNAseq_F37

F37TT <- cbind.data.frame(RNAseq_F37$t.test)

CohDF37<-matrix(0, length(RNAseq_F37[,1]), 4)

for(i in 1:length(F37TT[,1])){
  CohDF37[i, 1]<-tes(F37TT[i,], 6,6)$d
  CohDF37[i, 2]<-tes(F37TT[i,], 6,6)$'var.d'
  CohDF37[i, 3]<-tes(F37TT[i,], 6,6)$g
  CohDF37[i, 4]<-tes(F37TT[i,], 6,6)$'var.g'
}


colnames(CohDF37)<-c("d F37", "var d F37", "g F37", "var g F37")
CohDF37 <- cbind.data.frame(Symbol=RNAseq_F37$X, CohDF37)



#RNAseq_F43

#Getting Cohen's d from Cigdem Aydin's dataset
F43TT<-cbind.data.frame(RNAseq_F43$t.test)

CohDF43<-matrix(0, length(F43TT[,1]), 4)

for(i in 1:length(F43TT[,1])){
  CohDF43[i,1]<-tes(F43TT[i,], 5,5)$d
  CohDF43[i,2]<-tes(F43TT[i,], 5,5)$'var.d'
  CohDF43[i, 3]<-tes(F43TT[i,], 5,5)$g
  CohDF43[i, 4]<-tes(F43TT[i,], 5,5)$'var.g'
}


colnames(CohDF43)<-c("d F43", "var d F43", "g F43", "var g F43")
CohDF43<-cbind.data.frame(Symbol=RNAseq_F43$X, CohDF43)


#Affy_F4

#extracting TTstat for tes function
affyF4TT<-cbind.data.frame(Affy_F4$t.test)


CohDAffyF4<-matrix(0, length(Affy_F4[,1]), 4)

for(i in 1:length(affyF4TT[,1])){
  CohDAffyF4[i, 1]<-tes(affyF4TT[i,], 5,5)$d
  CohDAffyF4[i, 2]<-tes(affyF4TT[i,], 5,5)$'var.d'
  CohDAffyF4[i, 3]<-tes(affyF4TT[i,], 5,5)$g
  CohDAffyF4[i, 4]<-tes(affyF4TT[i,], 5,5)$'var.g'
}

colnames(CohDAffyF4)<-c("d Affy F4", "var d Affy F4", "g Affy F4", "var g Affy F4")
CohDAffyF4<-cbind.data.frame(Symbol=(Affy_F4$'X'), CohDAffyF4)


#NimbleGen_F34

ncF34TT<-cbind.data.frame(NimbleGen_F34[,2])

CohDncF34<-matrix(0, length(ncF34TT[,1]), 4)

for(i in 1:length(ncF34TT[,1])){
  CohDncF34[i, 1]<-tes(ncF34TT[i,], 5,5)$d
  CohDncF34[i, 2]<-tes(ncF34TT[i,], 5,5)$'var.d'
  CohDncF34[i, 3]<-tes(ncF34TT[i,], 5,5)$g
  CohDncF34[i, 4]<-tes(ncF34TT[i,], 5,5)$'var.g'
}

colnames(CohDncF34)<-c("d NC F34", "var d NC F34", "g NC F34", "var g NC F34")
CohDncF34<-cbind.data.frame(Symbol=NimbleGen_F34$X, CohDncF34)


## Meta-analysis
#The following code joins together the effect size dataframes for each of the five adult datasets and runs the meta-analysis

dfs<-list(
  CohDF43,
  CohDF37,
  CohDncF34,
  CohDAffyF4,
  CohDRNAseq
)

meta<-join_all(dfs, by="Symbol", type="full", match="all")
dim(meta)

sum(duplicated(meta$Symbol))


#Checking number of overlapping genes
allFive<-join_all(dfs, by="Symbol", type="inner")
length(!duplicated(allFive$Symbol))
#Number of genes found in ALL FIVE datasets.

colnames(meta)

#####Looking at correlation between cohen's d of different datasets

colnames(allFive) 
#Have to use one with only genes found in all datasets because cor can't handle NAs

#Grabbing only cohen's d columns
temp<-cbind(allFive$`d F43`, allFive$`d F37`, allFive$`d NC F34`, 
            allFive$`d RNAseq`, allFive$`d Affy F4`)

temp1<-as.matrix(temp)
cormatrix<-cor(temp1)

#output cor matrix
write.csv(cormatrix, paste0(outDir, "CohensdCorMatrix_AdultMeta.csv"))


#### Running Meta-Analysis

metaOutput<-matrix(0, length(meta$Symbol), 6)

for(i in 1:length(meta$Symbol)){
  effect<-c(meta$`d F43`[i], meta$`d F37`[i], meta$`d NC F34`[i], meta$`d RNAseq F29`[i], meta$`d Affy F4`[i])
  var<-c(meta$`var d F43`[i], meta$`var d F37`[i], meta$`var d NC F34`[i], meta$`var d RNAseq F29`[i], meta$`var d Affy F4`[i])
  if(sum(is.na(effect))>3){}
  else{
    metaOutput[i, 1]<-rma.mv(effect, var)$b #gives estimate
    metaOutput[i, 2]<-rma.mv(effect, var)$se #gives standard error
    metaOutput[i, 3]<-rma.mv(effect, var)$pval #gives pval
    metaOutput[i, 4]<-rma.mv(effect, var)$ci.lb #gives confidence interval lower bound
    metaOutput[i, 5]<-rma.mv(effect, var)$ci.ub #gives confidence interval upper bound
  }
  metaOutput[i, 6] <- sum(5 - sum(is.na(effect)))
}

colnames(metaOutput)<-c("estimate", "SE", "pval", "CI_lb", "CI_ub", "datasets")

metaOutputSymbol<-cbind.data.frame(GeneSymbol=meta$Symbol, metaOutput)


sum(metaOutputSymbol$pval==0) 
#Some genes are found in only one dataset and 
#therefor have no values associated with them from the meta-analysis

metaOutputSymbol <- subset(metaOutputSymbol, !(metaOutputSymbol$pval==0))
dim(metaOutputSymbol)
#total genes present in at least two datasets

write.csv(metaOutputSymbol, paste0(outDir, "AdultMetaAnalysisOutput.csv"))



## Multiple Comparison corrections
#The following code applies two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli) 
#Meta-analysis output with adjusted p-values is then output along with effect size information.


colnames(metaOutputSymbol)

tempPvalAdjMeta<-mt.rawp2adjp(metaOutputSymbol$pval, proc=c("BH", "BY"))
metaPvalAdj<-tempPvalAdjMeta$adjp[order(tempPvalAdjMeta$index),]

#adjusted pvalue object is in same orientation as metaoutput so can simply be binded together


metaAnalysisOutputFDR<-cbind.data.frame(metaPvalAdj, metaOutputSymbol)
dim(metaAnalysisOutputFDR)


write.csv(metaAnalysisOutputFDR, paste0(outDir, "AdultMetaAnalysisOutputFDR.csv"))


#Extracting sig genes (FDR < 0.05)

adultMetaSigGeneList<-subset(metaAnalysisOutputFDR, metaAnalysisOutputFDR$'BH' < 0.05)


write.csv(adultMetaSigGeneList, paste0(outDir, "AdultMetaSigGeneList.csv"))


#putting coh d and adjusted pval info together

#Renaming Symbol column from meta dataframe to join with meta output dataframe
temp<-meta
colnames(temp)[1] <- "GeneSymbol"

cohDandPval<-join(metaAnalysisOutputFDR, temp, by="GeneSymbol", type="inner")

write.csv(cohDandPval, paste0(outDir, "AdultMeta_CohDandPval.csv"))



## Cell Type MetaAnalysis
#The following code performs the meta-analysis on the cell type test statistic values, very similar to previous code.

# Extracting Cohen's D from all

#RNAseq_F43          
#sample size is 5
cellF43TT <- cbind.data.frame(cellRNAseq_F43$t.test)

CohD_CellType_F43<-matrix(0, length(cellF43TT[,1]), 4)

for(i in 1:length(cellF43TT[,1])){
  CohD_CellType_F43[i, 1]<-tes(cellF43TT[i,], 5,5)$d
  CohD_CellType_F43[i, 2]<-tes(cellF43TT[i,], 5,5)$'var.d'
  CohD_CellType_F43[i, 3]<-tes(cellF43TT[i,], 5,5)$g
  CohD_CellType_F43[i, 4]<-tes(cellF43TT[i,], 5,5)$'var.g'
}

row.names(CohD_CellType_F43)<-(cellRNAseq_F43$X)
colnames(CohD_CellType_F43)<-c("d F43", "var d F43", "g F43", "var g F43")



#RNAseq_F37         
#sample size is 6
cellF37TT <- cbind.data.frame(cellRNAseq_F37$t.test)

CohD_CellType_F37<-matrix(0, length(cellF37TT[,1]), 4)

for(i in 1:length(cellF37TT[,1])){
  CohD_CellType_F37[i, 1]<-tes(cellF37TT[i,], 6,6)$d
  CohD_CellType_F37[i, 2]<-tes(cellF37TT[i,], 6,6)$'var.d'
  CohD_CellType_F37[i, 3]<-tes(cellF37TT[i,], 6,6)$g
  CohD_CellType_F37[i, 4]<-tes(cellF37TT[i,], 6,6)$'var.g'
}

row.names(CohD_CellType_F37)<-cellRNAseq_F37$X
colnames(CohD_CellType_F37)<-c("d F37", "var d F37", "g F37", "var g F37")

#New Colony         
#sample size is 5
cellF34TT <- cbind.data.frame(cellNimbleGen_F34$Adult.TT)

CohD_CellType_NCF34<-matrix(0, length(cellF34TT[,1]), 4)

for(i in 1:length(cellF34TT[,1])){
  CohD_CellType_NCF34[i, 1]<-tes(cellF34TT[i,], 5,5)$d
  CohD_CellType_NCF34[i, 2]<-tes(cellF34TT[i,], 5,5)$'var.d'
  CohD_CellType_NCF34[i, 3]<-tes(cellF34TT[i,], 5,5)$g
  CohD_CellType_NCF34[i, 4]<-tes(cellF34TT[i,], 5,5)$'var.g'
}

row.names(CohD_CellType_NCF34)<-cellNimbleGen_F34$X
colnames(CohD_CellType_NCF34)<-c("d NC F34", "var d NC F34", "g NC F34", "var g NC F34")


#John Stead          
#sample size is 5
cellF4TT <- cbind.data.frame(cellAffy_F4$Stead_TT)

CohD_CellType_F4<-matrix(0, length(SteadTT[,1]), 4)

for(i in 1:length(cellF4TT[,1])){
  CohD_CellType_F4[i, 1]<-tes(cellF4TT[i,], 5,5)$d
  CohD_CellType_F4[i, 2]<-tes(cellF4TT[i,], 5,5)$'var.d'
  CohD_CellType_F4[i, 3]<-tes(cellF4TT[i,], 5,5)$g
  CohD_CellType_F4[i, 4]<-tes(cellF4TT[i,], 5,5)$'var.g'
}


row.names(CohD_CellType_F4)<-cellAffy_F4$X
colnames(CohD_CellType_F4)<-c("d Affy F4", "var d Affy F4", "g Affy F4", "var g Affy F4")



#### Cbinding Cohen's D for all datasets
#since each cohen's d dataset is in the same order and for the same 10 cell types, they can simply be binded together

CohD_Adult_CellType<-cbind.data.frame(CohD_CellType_F43, CohD_CellType_F37, CohD_CellType_NCF34, CohD_CellType_F4)

CellTypeMetaAnalysisOutput<-matrix(0, 10, 5)

for(i in 1:10){
  effect<-c(CohD_Adult_CellType$`d F43`[i], CohD_Adult_CellType$`d F37`[i], 
            CohD_Adult_CellType$`d NC F34`[i], CohD_Adult_CellType$`d Affy F4`[i])
  var<-c(CohD_Adult_CellType$`var d F43`[i], CohD_Adult_CellType$`var d F37`[i],
         CohD_Adult_CellType$`var d NC F34`[i], CohD_Adult_CellType$`var d Affy F4`[i])
  CellTypeMetaAnalysisOutput[i, 1]<-rma.mv(effect, var)$b #gives estimate
  CellTypeMetaAnalysisOutput[i, 2]<-rma.mv(effect, var)$se #gives standard error
  CellTypeMetaAnalysisOutput[i, 3]<-rma.mv(effect, var)$pval #gives pval
  CellTypeMetaAnalysisOutput[i, 4]<-rma.mv(effect, var)$ci.lb #gives confidence interval lower bound
  CellTypeMetaAnalysisOutput[i, 5]<-rma.mv(effect, var)$ci.ub #gives confidence interval upper bound
}

colnames(CellTypeMetaAnalysisOutput)<-c("estimate", "SE", "pval", "CI_lb", "CI_ub")
row.names(CellTypeMetaAnalysisOutput)<-row.names(CohD_Adult_CellType)

#joining with Cohen's D info

colnames(CohD_Adult_CellType)
temp<-cbind.data.frame(CellTypeMetaAnalysisOutput, CohD_Adult_CellType)

write.csv(temp, paste0(outDir, "Adult_CellTypeMeta_CohDandPval.csv"))


###############

## Forest Plots

#First exporting all significant genes as .png files and then exporting genes of interests as TIFF files

#Subsetting significant genes
sigMeta<-subset(cohDandPval, cohDandPval$BH < 0.05)
length(sigMeta$GeneSymbol)
#[1] 74


##Forest Plots with official thesis names for datasets


#Forloop over all sig genes and create forest plots for them

for(i in 1:length(sigMeta$GeneSymbol)){
  
  png(paste0(outForestPlots, "Forest Plot ", paste((sigMeta$GeneSymbol)[i]), ".png"))
  
  effect<-c(sigMeta$`d Affy F4`[i], sigMeta$`d RNAseq F29`[i], sigMeta$`d NC F34`[i], 
            sigMeta$`d F37`[i], sigMeta$`d F43`[i])
  
  var<-c(sigMeta$`var d Affy F4`[i], sigMeta$`var d RNAseq F29`[i], sigMeta$`var d NC F34`[i], 
         sigMeta$`var d F37`[i], sigMeta$`var d F43`[i])
  
  forest.rma(rma.mv(effect, var), slab=c("MBNI_AffymetrixRae230_F4", "MBNI_RNASeq_F29",  
                                         "Alabama_NimbleGen_F34", "MBNI_RNASeq_F37 ", 
                                         "MBNI_RNASeq_F43"), xlim=c(-18, 17))
  
  text(-12.7, 5.7, "Investigator & Study", cex=1)
  text(-13, 6.5, "LR", pos=4, cex=1.5)
  mtext(paste(sigMeta$GeneSymbol[i]), line=-1.5, cex=2)
  text(11.7, 5.7, "Cohen's D [95% CI] ", cex=1)
  text(12, 6.5, "HR ", pos=2, cex=1.5)
  
  dev.off()
  
}


### Export high res forest plots for genes of interest

#subsetting genes of interest

temp<-subset(meta, meta$Symbol=="Bmp4" | meta$Symbol=="Ncan" | meta$Symbol=="Trhr" |
               meta$Symbol=="C1qa" | meta$Symbol=="Fos" | meta$Symbol=="C1qb" | 
               meta$Symbol=="C1qc" | meta$Symbol=="Ucp2" |
               meta$Symbol=="Mfge8" | meta$Symbol=="Fxyd7")


#Export forest plots as TIFF files
for(i in 1:length(temp$Symbol)){
  
  tiff(paste0(outGOIPlots, "Forest Plot ", 
              paste((temp$Symbol)[i]), ".tiff"), res=300, compression="lzw",
       width = 8.5, height = 8, units = 'in')
  
  effect<-c(temp$`d Affy F4`[i], temp$`d RNAseq F29`[i], temp$`d NC F34`[i], 
            temp$`d F37`[i], temp$`d F43`[i])
  var<-c(temp$`var d Affy F4`[i], temp$`var d RNAseq F29`[i], temp$`var d NC F34`[i], 
         temp$`var d F37`[i], temp$`var d F43`[i])
  
  forest.rma(rma.mv(effect, var), slab=c("MBNI_AffymetrixRae230_F4", "MBNI_RNASeq_F29", 
                                         "Alabama_NimbleGen_F34", "MBNI_RNASeq_F37",
                                         "MBNI_RNASeq_F43"), 
             xlim=c(-18, 15), cex=1.6)
  text(-11, 5.6, "Investigator & Study", cex=1.6)
  text(-14, 6.6, "LR", cex=1.8)
  mtext(paste("Adult\n", temp$Symbol[i]), line=-1.5, cex=2.5)
  text(8, 5.6, "Cohen's D [95% CI]", cex=1.6)
  text(11, 6.6, "HR", cex=1.8)
  
  dev.off()
}



### Cell Type Forest Plots


#Exporting cell type forest plots as TIFF files

for(i in 1:10){
  
  tiff(paste0(outCellPlots, "Cell Type Forest Plot ", 
              paste(row.names(CohD_Adult_CellType)[i]), ".tiff"), res=300, 
       compression="lzw", width = 7.5, height = 7.5, units = 'in')
  
  effect<-c(CohD_Adult_CellType$`d Affy F4`[i], CohD_Adult_CellType$`d NC F34`[i], 
            CohD_Adult_CellType$`d F37`[i], CohD_Adult_CellType$`d F43`[i])
  var<-c(CohD_Adult_CellType$`var d Affy F4`[i], CohD_Adult_CellType$`var d NC F34`[i],
         CohD_Adult_CellType$`var d F37`[i], CohD_Adult_CellType$`var d F43`[i])
  
  forest.rma(rma.mv(effect, var), slab=c("MBNI_AffymetrixRae230_F4", "Alabama_NimbleGen_F34",
                                         "MBNI_RNASeq_F37", "MBNI_RNASeq_F43"), 
             xlim=c(-19, 15), cex=1.6)
  text(-12, 4.7, "Investigator & Study", cex=1.6)
  text(-15, 5.5, "LR", cex=1.8)
  mtext(paste("Adult\n", row.names(CohD_Adult_CellType)[i]), line=-1.5, cex=2.4)
  text(8, 4.7, "Cohen's D [95% CI]", cex=1.6)
  text(11, 5.5, "HR", cex=1.8)
  
  dev.off()
  
}




######### Volcano Plots


### Generating Volcano Plots to Observe Overall Spread of Expression Values (using estimate)

colnames(metaAnalysisOutputFDR)


############Adult expression colored by FDR and estimate

tiff(paste0(outDir, "VolcanoPlotAdult.tiff"), width = 5, height = 5, 
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))

par(mai=c(1.02, 1,0.9,0.40))
with(metaAnalysisOutputFDR, plot(estimate, -log10(pval), pch=19, main="Overall Expression", 
                                 xlim=c(-3,3), cex.lab=1.8, cex.main=2, cex=0.6))
# Add colored points: red if padj<0.05, blue of estimate>1)
with(subset(metaAnalysisOutputFDR, abs(estimate)>1), points(estimate, -log10(pval), 
                                                            pch=19, col="red", cex=0.6))
with(subset(metaAnalysisOutputFDR, BH< .05 ), points(estimate, -log10(pval), 
                                                     pch=19, col="blue", cex=0.6))
legend(-1.5, 7.6, legend=c("estimate > 1", "FDR < 0.05"), col=c("red", "blue"), pch=19, cex=1.2)

dev.off()



#############Overall Adult expression colored by dataset

tiff(paste0(outDir, "VolcanoPlotAdult_byDatasets.tiff"), width = 5, height = 5, 
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))

with(metaAnalysisOutputFDR, plot(estimate, -log10(pval), pch=19, main="Overall Expression by\n Dataset Number", 
                                 xlim=c(-3,3), cex.lab=1.8, cex.main=2, cex=0.6))

#color by dataset
with(subset(metaAnalysisOutputFDR, datasets > 4), points(estimate, -log10(pval), pch=19, col="blue", cex=0.6))
with(subset(metaAnalysisOutputFDR, datasets < 5 & datasets > 3), points(estimate, -log10(pval), pch=19, col="red", cex=0.6))
with(subset(metaAnalysisOutputFDR, datasets < 4 & datasets > 2), points(estimate, -log10(pval), pch=19, col="green3", cex=0.6))
with(subset(metaAnalysisOutputFDR, datasets < 3), points(estimate, -log10(pval), pch=19, col="purple", cex=0.6))

legend(-.5, 7.5, legend=c("2", "3", "4", "5"), 
       col=c("purple", "green3", "red", "blue"), pch=19, cex=1.2)


dev.off()




############Overall Adult Expression of cell-type specific genes
## Adding meta genes to cell type specific gene matrix in order to view overall expression of only the genes used in the cell type meta-analysis


#Read in list of cell type specific genes
cellInfo <- read.csv("data/CellTypeSpecificGenes_Master3.csv")
colnames(cellInfo)[5]<-"GeneSymbol" #renamae gene column to join

#Join data with cell type info to get only cell type specific genes
temp<-join(metaAnalysisOutputFDR, cellInfo, by="GeneSymbol", type="inner")

#Remove duplicates
temp<-subset(temp, !duplicated(temp$GeneSymbol))



#Output volcano plot
tiff(paste0(outDir, "VolcanoPlot CellTypeAdult.tiff"), width = 5, height = 5, units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.50))
with(temp, plot(estimate, -log10(pval), pch=19, main="Adult Overall Expression\nof Cell Type Specific Genes", 
                xlim=c(-3,3), cex.lab=1.8, cex.main=2, cex=0.6))

# Add colored points: red if padj<0.05, blue of estimate>1)
with(subset(temp, abs(estimate)>1), points(estimate, -log10(pval), 
                                           pch=19, col="red", cex=0.6))
with(subset(temp, BH< .05 ), points(estimate, -log10(pval), 
                                    pch=19, col="blue", cex=0.6))
legend(-1.5, 7.6, legend=c("estimate > 1", "FDR < 0.05"), col=c("red", "blue"), pch=19, cex=1.2)

dev.off()
