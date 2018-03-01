## Meta-analysis on Developmental HR/LR Data

#This script contains code running the meta-analysis of 5 P14 HR/LR datasets, and 2 HR/LR datasets for P7 and P21
#To run the meta-analysis, t-test of phenotype effect on gene expression is read-in for each dataset. Cohen's d is then extracted from t-test values and joined together for the meta-analysis.
#Note: one dataset, the RNA-seq F29 dataset, is the only dataset not re-analyzed prior to the meta-analysis as original data files were unavailable. Instead, Cohen's d is extracted from the p-values obtained through a prior analysis.

### Loading necessary libraries

library(plyr)
library(metafor)
library(compute.es)
library(bmp)
library(fields)
library(gplots)


source("https://bioconductor.org/biocLite.R")
biocLite("multtest")
library(multtest)

### Reading in data

#will read in all data and then extract Cohen's d

#Read in all the development t.test files
setwd("~/Phenotype Project/ibirt/output/MetaAnalysis Development/Read-in Files")

ttestData<-list.files()#must remove everything except for data files for the code to run this step
for (i in 1:length(ttestData)) assign(paste("",ttestData[i], sep = ""), read.csv(ttestData[i], row.names=1, header=TRUE))
rm(ttestData)


#Set directory back to main directory
setwd("~/Phenotype Project/ibirt")


#Set output directories

outDir <- "output/MetaAnalysis Development/"
outForestPlotsP7 <- "output/MetaAnalysis Development/ForestPlots_topGenesforP7/"
outForestPlotsP14 <- "output/MetaAnalysis Development/ForestPlots_topGenesforP14/"
outForestPlotsP21 <- "output/MetaAnalysis Development/ForestPlots_topGenesforP21/"
outCellPlotsP7 <- "output/MetaAnalysis Development/ForestPlots_CellType_P7/"
outCellPlotsP14 <- "output/MetaAnalysis Development/ForestPlots_CellType_P14/"
outCellPlotsP21 <- "output/MetaAnalysis Development/ForestPlots_CellType_P21/"
outGOIPlots<-"output/MetaAnalysis Development/Forest Plots for GOI in Development/"


###
## Obtaining Coh d for all using compute.es

#affy230 P14 *** sample size is 6 HR and 5 LR (outlier was thrown out)
CohDAffy230<-matrix(0, length(affy230P14TT.csv[,1]), 4)

for(i in 1:length(affy230P14TT.csv[,1])){
  CohDAffy230[i,1]<-tes(affy230P14TT.csv[i,], 6,5)$d
  CohDAffy230[i,2]<-tes(affy230P14TT.csv[i,], 6,5)$'var.d'
  CohDAffy230[i, 3]<-tes(affy230P14TT.csv[i,], 6,5)$g
  CohDAffy230[i, 4]<-tes(affy230P14TT.csv[i,], 6,5)$'var.g'
}

#assigning colnames and adding gene column
colnames(CohDAffy230)<-c("d affy", "var.d affy", "g affy", "var.g affy")
CohDAffy230 <- cbind.data.frame(gene=row.names(affy230P14TT.csv), CohDAffy230)

sum(!duplicated(CohDAffy230$gene))
#9958


##Coh d for development P7 *** sample size is 6
CohDdevelopmentP7<-matrix(0, length(DevelopmentP7TT.csv[,1]), 4)

for(i in 1:length(DevelopmentP7TT.csv[,1])){
  CohDdevelopmentP7[i,1]<-tes(DevelopmentP7TT.csv[i,], 6,6)$d
  CohDdevelopmentP7[i,2]<-tes(DevelopmentP7TT.csv[i,], 6,6)$'var.d'
  CohDdevelopmentP7[i, 3]<-tes(DevelopmentP7TT.csv[i,], 6,6)$g
  CohDdevelopmentP7[i, 4]<-tes(DevelopmentP7TT.csv[i,], 6,6)$'var.g'
}

#assigning colnames and adding gene column
colnames(CohDdevelopmentP7)<-c("d P7", "var.d P7", "g P7", "var.g P7")
CohDdevelopmentP7 <- cbind.data.frame(gene=row.names(DevelopmentP7TT.csv), CohDdevelopmentP7)


##Coh d for development P14 *** sample size is 6
CohDdevelopmentP14<-matrix(0, length(DevelopmentP14TT.csv[,1]), 4)

for(i in 1:length(DevelopmentP14TT.csv[,1])){
  CohDdevelopmentP14[i,1]<-tes(DevelopmentP14TT.csv[i,], 6,6)$d
  CohDdevelopmentP14[i,2]<-tes(DevelopmentP14TT.csv[i,], 6,6)$'var.d'
  CohDdevelopmentP14[i, 3]<-tes(DevelopmentP14TT.csv[i,], 6,6)$g
  CohDdevelopmentP14[i, 4]<-tes(DevelopmentP14TT.csv[i,], 6,6)$'var.g'
}


#assigning colnames and adding gene column
colnames(CohDdevelopmentP14)<-c("d P14", "var.d P14", "g P14", "var.g P14")
CohDdevelopmentP14 <- cbind.data.frame(gene=row.names(DevelopmentP14TT.csv), CohDdevelopmentP14)

sum(!duplicated(CohDdevelopmentP14$gene))
#[1] 4588


##Coh d for development P21 *** sample size is 6
CohDdevelopmentP21<-matrix(0, length(DevelopmentP21TT.csv[,1]), 4)

for(i in 1:length(DevelopmentP21TT.csv[,1])){
  CohDdevelopmentP21[i,1]<-tes(DevelopmentP21TT.csv[i,], 6,6)$d
  CohDdevelopmentP21[i,2]<-tes(DevelopmentP21TT.csv[i,], 6,6)$'var.d'
  CohDdevelopmentP21[i, 3]<-tes(DevelopmentP21TT.csv[i,], 6,6)$g
  CohDdevelopmentP21[i, 4]<-tes(DevelopmentP21TT.csv[i,], 6,6)$'var.g'
}

#assigning colnames and adding gene column
colnames(CohDdevelopmentP21)<-c("d P21", "var.d P21", "g P21", "var.g P21")
CohDdevelopmentP21 <- cbind.data.frame(gene=row.names(DevelopmentP21TT.csv), CohDdevelopmentP21)


##Coh d for illumina P14 *** sample size is 6
CohDilluminaP14<-matrix(0, length(illuminaP14TT.csv[,1]), 4)

for(i in 1:length(illuminaP14TT.csv[,1])){
  CohDilluminaP14[i,1]<-tes(illuminaP14TT.csv[i,], 6,6)$d
  CohDilluminaP14[i,2]<-tes(illuminaP14TT.csv[i,], 6,6)$'var.d'
  CohDilluminaP14[i, 3]<-tes(illuminaP14TT.csv[i,], 6,6)$g
  CohDilluminaP14[i, 4]<-tes(illuminaP14TT.csv[i,], 6,6)$'var.g'
}

#assigning colnames and adding gene column
colnames(CohDilluminaP14)<-c("d illumina", "var.d illumina", "g illumina", "var.g illumina")
CohDilluminaP14 <- cbind.data.frame(gene=row.names(illuminaP14TT.csv), CohDilluminaP14)

sum(!duplicated(CohDilluminaP14$gene))
#[1] 21568


##Coh d for New Colony development P7 *** sample size is 5
CohDnewcolonydevelopmentP7<-matrix(0, length(newdevelopmentdataP7TT.csv[,1]), 4)

for(i in 1:length(newdevelopmentdataP7TT.csv[,1])){
  CohDnewcolonydevelopmentP7[i,1]<-tes(newdevelopmentdataP7TT.csv[i,], 5,5)$d
  CohDnewcolonydevelopmentP7[i,2]<-tes(newdevelopmentdataP7TT.csv[i,], 5,5)$'var.d'
  CohDnewcolonydevelopmentP7[i, 3]<-tes(newdevelopmentdataP7TT.csv[i,], 5,5)$g
  CohDnewcolonydevelopmentP7[i, 4]<-tes(newdevelopmentdataP7TT.csv[i,], 5,5)$'var.g'
}

#assigning colnames and adding gene column
colnames(CohDnewcolonydevelopmentP7)<-c("d new colony P7", "var.d new colony P7", "g new colony P7", "var.g new colony P7")
CohDnewcolonydevelopmentP7 <- cbind.data.frame(gene=row.names(newdevelopmentdataP7TT.csv), CohDnewcolonydevelopmentP7)


##Coh d for New Colony development P14 *** sample size is 5
CohDnewcolonydevelopmentP14<-matrix(0, length(newdevelopmentdataP14TT.csv[,1]), 4)

for(i in 1:length(newdevelopmentdataP14TT.csv[,1])){
  CohDnewcolonydevelopmentP14[i,1]<-tes(newdevelopmentdataP14TT.csv[i,], 5,5)$d
  CohDnewcolonydevelopmentP14[i,2]<-tes(newdevelopmentdataP14TT.csv[i,], 5,5)$'var.d'
  CohDnewcolonydevelopmentP14[i, 3]<-tes(newdevelopmentdataP14TT.csv[i,], 5,5)$g
  CohDnewcolonydevelopmentP14[i, 4]<-tes(newdevelopmentdataP14TT.csv[i,], 5,5)$'var.g'
}


#assigning colnames and adding gene column
colnames(CohDnewcolonydevelopmentP14)<-c("d new colony P14", "var.d new colony P14", "g new colony P14", "var.g new colony P14")
CohDnewcolonydevelopmentP14 <- cbind.data.frame(gene=row.names(newdevelopmentdataP14TT.csv), CohDnewcolonydevelopmentP14)


sum(!duplicated(CohDnewcolonydevelopmentP14$gene))
#[1] 10674


##Coh d for new colony development P21 *** sample size is 5 HR and 4 LR (one removed as outlier)
CohDnewcolonydevelopmentP21<-matrix(0, length(newdevelopmentdataP21TT.csv[,1]), 4)

for(i in 1:length(newdevelopmentdataP21TT.csv[,1])){
  CohDnewcolonydevelopmentP21[i,1]<-tes(newdevelopmentdataP21TT.csv[i,], 5,4)$d
  CohDnewcolonydevelopmentP21[i,2]<-tes(newdevelopmentdataP21TT.csv[i,], 5,4)$'var.d'
  CohDnewcolonydevelopmentP21[i, 3]<-tes(newdevelopmentdataP21TT.csv[i,], 5,4)$g
  CohDnewcolonydevelopmentP21[i, 4]<-tes(newdevelopmentdataP21TT.csv[i,], 5,4)$'var.g'
}


#assigning colnames and adding gene column
colnames(CohDnewcolonydevelopmentP21)<-c("d new colony P21", "var.d new colony P21", "g new colony P21", "var.g new colony P21")
CohDnewcolonydevelopmentP21 <- cbind.data.frame(gene=row.names(newdevelopmentdataP21TT.csv), CohDnewcolonydevelopmentP21)



### RNAseq_F29 P14 data

##removing data that has inaccurate pvals and no gene symbol
RNAseqdata<-subset(LRHR_rgsc34_th2out_g1tx_diff_P14LvsP14H_RNAseqnof2.csv, LRHR_rgsc34_th2out_g1tx_diff_P14LvsP14H_RNAseqnof2.csv$p_value != 1)
RNAseqdata<-subset(RNAseqdata, RNAseqdata$gene!="-")


sum(duplicated(RNAseqdata$gene)) #still some duplicated genes after removal of data without a gene symbol
#[1] 310

rnaseq_nobad<-subset(RNAseqdata, RNAseqdata$test_stat > -1.79769e+308 & RNAseqdata$test_stat < 1.79769e+308)

dim(rnaseq_nobad)

sum(duplicated(rnaseq_nobad$gene))
#[1] 304    still lots of duplicated genes. Many of these seem to have very different expression values.

#Averaging across duplicate genes

colnames(rnaseq_nobad)
#[1] "gene_id"           "gene"              "locus"             "sample_1"          "sample_2"         
#[6] "status"            "value_1"           "value_2"           "log2.fold_change." "test_stat"        
#[11] "p_value"           "q_value"           "significant"  

rnaseq_avg<-tapply(rnaseq_nobad[,10], as.character(rnaseq_nobad[,2]), function(y) mean(y))


#extracting the test statistic averaged by gene
RNAseqdataTT<-cbind.data.frame(rnaseq_avg)


#extracting coh D

CohDRNAseqTT<-matrix(0, length(RNAseqdataTT[,1]), 4)

for(i in 1:length(RNAseqdataTT[,1])){
  CohDRNAseqTT[i, 1]<-tes(RNAseqdataTT[i,], 2,2)$d
  CohDRNAseqTT[i, 2]<-tes(RNAseqdataTT[i,], 2,2)$'var.d'
  CohDRNAseqTT[i, 3]<-tes(RNAseqdataTT[i,], 2,2)$g
  CohDRNAseqTT[i, 4]<-tes(RNAseqdataTT[i,], 2,2)$'var.g'
}

colnames(CohDRNAseqTT)<-c("d RNAseq", "var.d RNAseq", "g RNAseq", "var.g RNAseq")
CohDRNAseqTT<-cbind.data.frame(gene=row.names(RNAseqdataTT), CohDRNAseqTT)

#Need to multiply effect sizes by -1 because it is in the opposite direction of the t.tests from the other datasets

CohDRNAseq_samedoe<-(-1*CohDRNAseqTT[[2]]) ###multiplying all Coh D by -1 so it matches DOE from other ttest output

#adding gene, variance of d and corrected cohen's d into one dataframe
CohDRNAseq_samedoe<-cbind.data.frame(gene=CohDRNAseqTT$gene, 'd RNAseq'=CohDRNAseq_samedoe, 'var.d RNAseq'=CohDRNAseqTT$`var.d RNAseq`)

length(CohDRNAseq_samedoe$gene)
#[1] 16342

sum(!duplicated(CohDRNAseq_samedoe$gene))
#[1] 16342



## Meta-analysis by Age Group
#Join cohen's d dataframes for each dataset together for each age group then run meta-analysis on each age group using metafor


#### P7 meta-anaysis


#Join by inner since there is only two datasets, therefor genes must be present in both to be run through the meta-analysis
P7<-join(CohDdevelopmentP7, CohDnewcolonydevelopmentP7, by="gene", type="inner") 
dim(P7)
#[1] 3527    9 ### 3,527 genes found in both datasets


MetaAnalysisOutputP7<-matrix(0, length(P7$gene), 5)

######  
for(i in 1:length(P7$gene)){
  effect<-c(P7$`d P7`[i], P7$`d new colony P7`[i])
  var<-c(P7$`var.d P7`[i], P7$`var.d new colony P7`[i])
  MetaAnalysisOutputP7[i, 1]<-rma.mv(effect, var)$b #gives estimate
  MetaAnalysisOutputP7[i, 2]<-rma.mv(effect, var)$se #gives standard error
  MetaAnalysisOutputP7[i, 3]<-rma.mv(effect, var)$pval #gives pval
  MetaAnalysisOutputP7[i, 4]<-rma.mv(effect, var)$ci.lb #gives confidence interval lower bound
  MetaAnalysisOutputP7[i, 5]<-rma.mv(effect, var)$ci.ub #gives confidence interval upper bound
}


colnames(MetaAnalysisOutputP7)<-c("estimate", "SE", "pval", "CI_lb", "CI_ub")
row.names(MetaAnalysisOutputP7)<-(P7$gene)

write.csv(MetaAnalysisOutputP7, paste0(outDir, "MetaAnalysisOutputP7.csv"))


### Multiple Comparison corrections for P7 Meta-Analysis

#Applying two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli):

TempPvalAdjP7Meta<-mt.rawp2adjp(MetaAnalysisOutputP7[,3], proc=c("BH", "BY"))
MetaPvalAdjP7<-TempPvalAdjP7Meta$adjp[order(TempPvalAdjP7Meta$index),]
dim(MetaPvalAdjP7)
#[1] 3527    3

row.names(MetaPvalAdjP7)<-row.names(MetaAnalysisOutputP7)

temp<-cbind.data.frame(gene=row.names(MetaPvalAdjP7), MetaPvalAdjP7)
P7metagene<-cbind.data.frame(gene=row.names(MetaAnalysisOutputP7), MetaAnalysisOutputP7)
MetaAnalysisOutputP7FDR<-join(temp, P7metagene, by="gene", type="inner")

write.csv(MetaAnalysisOutputP7FDR, paste0(outDir, "MetaAnalysisOutputP7FDR.csv"))


##putting all P7 info together
temp<-cbind.data.frame(gene=MetaAnalysisOutputP7FDR[1], MetaAnalysisOutputP7FDR[3:9])
P7cohDandPval<-join(P7, temp, by="gene", type="inner")

write.csv(P7cohDandPval, paste0(outDir, "P7cohDandPval.csv"))

#looking at number of significant genes
sum(P7cohDandPval$pval < 0.05)
#[1] 392

sum(P7cohDandPval$BH < 0.1)
#[1] 1

###### P14 Meta-Analysis
#Join all P14 datasets together, run meta-analysis, and perform multiple comparisons correction


dfs<-list(
  
  CohDdevelopmentP14,
  CohDnewcolonydevelopmentP14,
  CohDAffy230,
  CohDilluminaP14,
  CohDRNAseq_samedoe
  
)

P14<-join_all(dfs, by="gene", type="full")
dim(P14)
#[1] 30556   19
sum(duplicated(P14$gene))
#[1] 0


#Checking number of genes found in all 5 datasets
temp<-join_all(dfs, by="gene", type="inner")
length(temp$gene)
#[1] 2353      genes found in all datasets, only roughly 1/13              



## Meta-analysis

MetaAnalysisOutputP14<-matrix(0, length(P14$gene), 6)

for(i in 1:length(P14$gene)){
  effect<-c(P14$`d P14`[i], P14$`d new colony P14`[i], P14$`d affy`[i], P14$`d illumina`[i], P14$`d RNAseq`[i])
  var<-c(P14$`var.d P14`[i], P14$`var.d new colony P14`[i], P14$`var.d affy`[i], P14$`var.d illumina`[i], P14$`var.d RNAseq`[i])
  if(sum(is.na(effect))>3){}
  else{
    MetaAnalysisOutputP14[i, 1]<-rma.mv(effect, var)$b #gives estimate
    MetaAnalysisOutputP14[i, 2]<-rma.mv(effect, var)$se #gives standard error
    MetaAnalysisOutputP14[i, 3]<-rma.mv(effect, var)$pval #gives pval
    MetaAnalysisOutputP14[i, 4]<-rma.mv(effect, var)$ci.lb #gives confidence interval lower bound
    MetaAnalysisOutputP14[i, 5]<-rma.mv(effect, var)$ci.ub #gives confidence interval upper bound
  }
  MetaAnalysisOutputP14[i, 6] <- sum(5 - sum(is.na(effect)))
}

colnames(MetaAnalysisOutputP14)<-c("estimate", "SE", "pval", "CI_lb", "CI_ub", "Datasets")
MetaAnalysisOutputP14<-cbind.data.frame(gene=P14$gene, MetaAnalysisOutputP14)
dim(MetaAnalysisOutputP14)
#[1] 30556     7



###Multiple Comparison corrections###
#Applying two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli):


#Removing genes with estimates & pvalues=0
nozeros <- subset(MetaAnalysisOutputP14, !(MetaAnalysisOutputP14$estimate==0))

length(nozeros$gene)
#[1] 15682
colnames(nozeros)
#[1] "gene"     "estimate" "SE"       "pval"     "CI_lb"    "CI_ub"    "Datasets" 

TempPvalAdjP14Meta<-mt.rawp2adjp(nozeros[,4], proc=c("BH", "BY"))
MetaPvalAdjP14<-TempPvalAdjP14Meta$adjp[order(TempPvalAdjP14Meta$index),]
dim(MetaPvalAdjP14)
#[1] 15682     3                  

row.names(MetaPvalAdjP14) <- nozeros$gene


MetaAnalysisOutputP14FDRnozeros<-cbind.data.frame(nozeros, MetaPvalAdjP14)

#looking at number of significant genes
sum(MetaAnalysisOutputP14FDRnozeros$pval < 0.05)
#[1] 921

sum(MetaAnalysisOutputP14FDRnozeros$BH < 0.1)
#[1] 0


write.csv(MetaAnalysisOutputP14FDRnozeros, paste0(outDir, "MetaAnalysisOutputP14FDR.csv"))


P14cohDandPval <- join(MetaAnalysisOutputP14FDRnozeros, P14, by="gene", type="inner")


write.csv(P14cohDandPval, paste0(outDir, "p14Meta_cohDandpval.csv"))



### P21 Meta-Analysis
#very similar to the P7 meta-analysis, run in the same way

DevelopmentP21<-cbind.data.frame(gene=row.names(CohDdevelopmentP21), CohDdevelopmentP21)
NewcolonyP21<-cbind.data.frame(gene=row.names(CohDnewcolonydevelopmentP21), CohDnewcolonydevelopmentP21)


P21<-join(CohDdevelopmentP21, CohDnewcolonydevelopmentP21, by="gene", type="inner")
dim(P21)
#[1] 3527    9

MetaAnalysisOutputP21<-matrix(0, length(P21$gene), 5)

for(i in 1:length(P21$gene)){
  effect<-c(P21$`d P21`[i], P21$`d new colony P21`[i])
  var<-c(P21$`var.d P21`[i], P21$`var.d new colony P21`[i])
  MetaAnalysisOutputP21[i, 1]<-rma.mv(effect, var)$b #gives estimate
  MetaAnalysisOutputP21[i, 2]<-rma.mv(effect, var)$se #gives standard error
  MetaAnalysisOutputP21[i, 3]<-rma.mv(effect, var)$pval #gives pval
  MetaAnalysisOutputP21[i, 4]<-rma.mv(effect, var)$ci.lb #gives confidence interval lower bound
  MetaAnalysisOutputP21[i, 5]<-rma.mv(effect, var)$ci.ub #gives confidence interval upper bound
}

colnames(MetaAnalysisOutputP21)<-c("estimate", "SE", "pval", "CI_lb", "CI_ub")
row.names(MetaAnalysisOutputP21)<-(P21$gene)

write.csv(MetaAnalysisOutputP21, paste0(outDir, "MetaAnalysisOutputP21.csv"))



##FDR on P21 Meta pvals

###Multiple Comparison corrections###

#Applying two different types of multiple-comparison corrections to the raw p-values (Benjamini-Hochberg and Benjamini-Yekutieli):

TempPvalAdjP21Meta<-mt.rawp2adjp(MetaAnalysisOutputP21[,3], proc=c("BH", "BY"))
MetaPvalAdjP21<-TempPvalAdjP21Meta$adjp[order(TempPvalAdjP21Meta$index),]
dim(MetaPvalAdjP21)

row.names(MetaPvalAdjP21)<-row.names(MetaAnalysisOutputP21)

temp<-cbind.data.frame(gene=row.names(MetaPvalAdjP21), MetaPvalAdjP21)
P21metagene<-cbind.data.frame(gene=row.names(MetaAnalysisOutputP21), MetaAnalysisOutputP21)
MetaAnalysisOutputP21FDR<-join(temp, P21metagene, by="gene", type="inner")


write.csv(MetaAnalysisOutputP21FDR, paste0(outDir, "MetaAnalysisOutputP21FDR.csv"))

##putting all P21 info together
temp<-cbind.data.frame(gene=MetaAnalysisOutputP21FDR[1], MetaAnalysisOutputP21FDR[3:9])
P21cohDandPval<-join(temp, P21, by="gene", type="inner")

write.csv(P21cohDandPval, paste0(outDir, "P21cohDandPval.csv"))


## Looking at number of significant genes

sum(P21cohDandPval$pval < 0.05)
#[1] 263

sum(P21cohDandPval$BH < 0.1)
#[1] 0

####### 


################ Forest Plots
#Outputting forest plots for a subset of significant genes for each age group


#P7
sigP7<-subset(P7cohDandPval, P7cohDandPval$pval < 0.001)
length(sigP7$gene)
#[1] 22


for(i in 1:22){
  
  png(paste0(outForestPlotsP7, "Forest Plot ", paste((sigP7$gene)[i]), ".png"))
  
  effect<-c(sigP7$`d P7`[i], sigP7$`d new colony P7`[i])
  var<-c(sigP7$`var.d P7`[i], sigP7$`var.d new colony P7`[i])
  
  forest.rma(rma.mv(effect, var), slab=c("MBNI_AffymetrixRgU34A_F6", "Alabama_NimbleGen_F34"), 
             xlim=c(-24, 17))
  text(-21.5, 2.6, "Investigator & Study", pos=4, cex=1)
  text(-23, 3.5, "LR", pos=4, cex=1.5)
  mtext(paste(sigP7$gene[i]), line=-2.2, cex=2.5)
  text(14.5, 2.6, "Cohen's D [95% CI]", pos=2, cex=1)
  text(15.5, 3.5, "HR ", pos=2, cex=1.5)
  dev.off()
  
}




########   P14 Forest Plots

sigP14<-subset(P14cohDandPval, P14cohDandPval$pval < 0.001)
length(sigP14$gene)
#[1] 29



for(i in 1:29){

  png(paste0(outForestPlotsP14, "Forest Plot ", 
             paste((sigP14$gene)[i]), ".png"))
  
  effect<-c(sigP14$`d P14`[i], sigP14$`d affy`[i], sigP14$`d illumina`[i], sigP14$`d RNAseq`[i], 
            sigP14$`d new colony P14`[i])
  var<-c(sigP14$`var.d P14`[i], sigP14$`var.d affy`[i], sigP14$`var.d illumina`[i], 
         sigP14$`var.d RNAseq`[i], sigP14$`var.d new colony P14`[i])
  
  forest.rma(rma.mv(effect, var), slab=c("MBNI_AffymetrixRgU34A_F6", "MBNI_AffymetrixRae230_F15", 
                                         "MBNI_IlluminaRatRef12v1_F15", "MBNI_RNASeq_F29", "Alabama_NimbleGen_F34"), 
             xlim=c(-19.5, 17.2))
  text(-19, 5.6, "Investigator & Study", pos=4, cex=.75)
  text(-16, 6.8, "LR")
  mtext(paste(sigP14$gene[i]), line=-1.5, cex=2)
  text(16.6, 5.6, "Cohen's D [95% CI]", pos=2, cex=.75)
  text(14, 6.8, "HR")
  dev.off()
}



#########   P21 Forest Plots

sigP21<-subset(P21cohDandPval, P21cohDandPval$pval < 0.001)
length(sigP21$gene)
#[1] 9


for(i in 1:9){

  png(paste0(outForestPlotsP21, "Forest Plot ", 
            paste((sigP21$gene)[i]), ".png"))
  
  effect<-c(sigP21$`d P21`[i], sigP21$`d new colony P21`[i])
  var<-c(sigP21$`var.d P21`[i], sigP21$`var.d new colony P21`[i])
  
  forest.rma(rma.mv(effect, var), slab=c("MBNI_AffymetrixRgU34A_F6", "Alabama_NimbleGen_F34"), 
             xlim=c(-15, 17))
  text(-12.5, 2.6, "Investigator & Study", pos=4, cex=1)
  text(-14, 3.5, "LR", pos=4, cex=1.5)
  mtext(paste(sigP21$gene[i]), line=-2.2, cex=2.5)
  text(14.5, 2.6, "Cohen's D [95% CI]", pos=2, cex=1)
  text(15.5, 3.5, "HR ", pos=2, cex=1.5)
  dev.off()
}




###########################################################################


#####   Cell Type Meta-analysis
#Perform cell type meta-analysis for each age group


colnames(DevelopmentCellTypeTT_forMeta.csv)
#[1] "P7.TT"  "P14.TT" "P21.TT"


#splitting development cell type t-test dataset into 3 different age groups

DevelopmentP7TT_celltype<-cbind.data.frame(DevelopmentCellTypeTT_forMeta.csv$P7.TT)
row.names(DevelopmentP7TT_celltype)<-row.names(DevelopmentCellTypeTT_forMeta.csv)

DevelopmentP14TT_celltype<-cbind.data.frame(DevelopmentCellTypeTT_forMeta.csv$P14.TT)
row.names(DevelopmentP14TT_celltype)<-row.names(DevelopmentCellTypeTT_forMeta.csv)

DevelopmentP21TT_celltype<-cbind.data.frame(DevelopmentCellTypeTT_forMeta.csv$P21.TT)
row.names(DevelopmentP21TT_celltype)<-row.names(DevelopmentCellTypeTT_forMeta.csv)

###

#Splitting New Colony TT dataset into 3 different age groups (don't need the adult group for this analysis)

colnames(NewColonyDevelopmentCellTypeTT_forMeta.csv)
#[1] "P7.TT"    "P14.TT"   "P21.TT"   "Adult.TT"

NewColonyDevelopmentP7TT_celltype<-cbind.data.frame(NewColonyDevelopmentCellTypeTT_forMeta.csv$P7.TT)
row.names(NewColonyDevelopmentP7TT_celltype)<-row.names(NewColonyDevelopmentCellTypeTT_forMeta.csv)

NewColonyDevelopmentP14TT_celltype<-cbind.data.frame(NewColonyDevelopmentCellTypeTT_forMeta.csv$P14.TT)
row.names(NewColonyDevelopmentP14TT_celltype)<-row.names(NewColonyDevelopmentCellTypeTT_forMeta.csv)

NewColonyDevelopmentP21TT_celltype<-cbind.data.frame(NewColonyDevelopmentCellTypeTT_forMeta.csv$P21.TT)
row.names(NewColonyDevelopmentP21TT_celltype)<-row.names(NewColonyDevelopmentCellTypeTT_forMeta.csv)



### Obtaining Coh D for P7 Groups

##Coh d for cell type development P7 *** sample size is 6
CellType_CohDdevelopmentP7<-matrix(0, length(DevelopmentP7TT_celltype[,1]), 4)

for(i in 1:length(DevelopmentP7TT_celltype[,1])){
  CellType_CohDdevelopmentP7[i,1]<-tes(DevelopmentP7TT_celltype[i,], 6,6)$d
  CellType_CohDdevelopmentP7[i,2]<-tes(DevelopmentP7TT_celltype[i,], 6,6)$'var.d'
  CellType_CohDdevelopmentP7[i, 3]<-tes(DevelopmentP7TT_celltype[i,], 6,6)$g
  CellType_CohDdevelopmentP7[i, 4]<-tes(DevelopmentP7TT_celltype[i,], 6,6)$'var.g'
}

row.names(CellType_CohDdevelopmentP7)<-row.names(DevelopmentP7TT_celltype)
colnames(CellType_CohDdevelopmentP7)<-c("d P7", "var.d P7", "g P7", "var.g P7")


##Coh d for New Colony cell type P7 *** sample size is 5

CellType_CohDNewColonyP7<-matrix(0, length(NewColonyDevelopmentP7TT_celltype[,1]), 4)

for(i in 1:length(NewColonyDevelopmentP7TT_celltype[,1])){
  CellType_CohDNewColonyP7[i,1]<-tes(NewColonyDevelopmentP7TT_celltype[i,], 5,5)$d
  CellType_CohDNewColonyP7[i,2]<-tes(NewColonyDevelopmentP7TT_celltype[i,], 5,5)$'var.d'
  CellType_CohDNewColonyP7[i, 3]<-tes(NewColonyDevelopmentP7TT_celltype[i,], 5,5)$g
  CellType_CohDNewColonyP7[i, 4]<-tes(NewColonyDevelopmentP7TT_celltype[i,], 5,5)$'var.g'
}

row.names(CellType_CohDNewColonyP7)<-row.names(NewColonyDevelopmentCellTypeTT_forMeta.csv)
colnames(CellType_CohDNewColonyP7)<-c("d NC P7", "var.d NC P7", "g NC P7", "var.g NC P7")


########## Meta-Analysis for P7 cell type


##Cell types (row.names) are all in same order so I can simply cbind the different datasets together
CellType_P7<-cbind.data.frame(cell.type=row.names(CellType_CohDdevelopmentP7), CellType_CohDdevelopmentP7, CellType_CohDNewColonyP7)

colnames(CellType_P7)
#[1] "cell.type"   "d P7"        "var.d P7"    "g P7"        "var.g P7"    "d NC P7"     "var.d NC P7" "g NC P7"    
#[9] "var.g NC P7"

CellType_MetaAnalysisOutputP7<-matrix(0, 10, 5)

for(i in 1:10){
  effect<-c(CellType_P7$`d P7`[i], CellType_P7$`d NC P7`[i])
  var<-c(CellType_P7$`var.d P7`[i], CellType_P7$`var.d NC P7`[i])
  CellType_MetaAnalysisOutputP7[i, 1]<-rma.mv(effect, var)$b #gives estimate
  CellType_MetaAnalysisOutputP7[i, 2]<-rma.mv(effect, var)$se #gives standard error
  CellType_MetaAnalysisOutputP7[i, 3]<-rma.mv(effect, var)$pval #gives pval
  CellType_MetaAnalysisOutputP7[i, 4]<-rma.mv(effect, var)$ci.lb #gives confidence interval lower bound
  CellType_MetaAnalysisOutputP7[i, 5]<-rma.mv(effect, var)$ci.ub #gives confidence interval upper bound
}

colnames(CellType_MetaAnalysisOutputP7)<-c("estimate", "SE", "pval", "CI_lb", "CI_ub")
row.names(CellType_MetaAnalysisOutputP7)<-row.names(CellType_P7)


write.csv(CellType_MetaAnalysisOutputP7, paste0(outDir, "CellType_MetaAnalysisOutputP7.csv"))



##### Forest Plots for Cell Analysis of P7


for(i in 1:10){
  
  png(paste0(outCellPlotsP7, "Forest Plot ", 
             paste((row.names(CellType_P7))[i]), ".png"))
  
  effect<-c(CellType_P7$`d P7`[i], CellType_P7$`d NC P7`[i])
  var<-c(CellType_P7$`var.d P7`[i], CellType_P7$`var.d NC P7`[i])
  
  forest.rma(rma.mv(effect, var), slab=c("MBNI_AffymetrixRgU34A_F6", "Alabama_NimbleGen_F34"), 
             xlim=c(-15, 17))
  text(-15, 2.6, "Investigator & Study", pos=4, cex=.75)
  mtext(paste(row.names(CellType_P7)[i]), line=-1.5, cex=2)
  text(16.5, 2.6, "Cohen's D [95% CI]", pos=2, cex=.75)
  dev.off()
  
}



####  P14 cell type meta-analysis

colnames(P14Affy230CellTypeTT_forMeta.csv)
#[1] "Affy.P14.t.test"

colnames(P14IlluminaCellTypeTT_forMeta.csv)
#[1] "Illumina.t.test"


### Obtaining Coh d for all P14 cell type groups


#Coh d for New Colony cell type P14 *** sample size is 5

CellType_CohDNewColonyP14<-matrix(0, length(NewColonyDevelopmentP14TT_celltype[,1]), 4)

for(i in 1:length(NewColonyDevelopmentP14TT_celltype[,1])){
  CellType_CohDNewColonyP14[i,1]<-tes(NewColonyDevelopmentP14TT_celltype[i,], 5,5)$d
  CellType_CohDNewColonyP14[i,2]<-tes(NewColonyDevelopmentP14TT_celltype[i,], 5,5)$'var.d'
  CellType_CohDNewColonyP14[i, 3]<-tes(NewColonyDevelopmentP14TT_celltype[i,], 5,5)$g
  CellType_CohDNewColonyP14[i, 4]<-tes(NewColonyDevelopmentP14TT_celltype[i,], 5,5)$'var.g'
}

row.names(CellType_CohDNewColonyP14)<-row.names(NewColonyDevelopmentCellTypeTT_forMeta.csv)
colnames(CellType_CohDNewColonyP14)<-c("d NC P14", "var.d NC P14", "g NC P14", "var.g NC P14")

##Coh d for cell type development P14 *** sample size is 6
CellType_CohDdevelopmentP14<-matrix(0, length(DevelopmentP14TT_celltype[,1]), 4)

for(i in 1:length(DevelopmentP14TT_celltype[,1])){
  CellType_CohDdevelopmentP14[i,1]<-tes(DevelopmentP14TT_celltype[i,], 6,6)$d
  CellType_CohDdevelopmentP14[i,2]<-tes(DevelopmentP14TT_celltype[i,], 6,6)$'var.d'
  CellType_CohDdevelopmentP14[i, 3]<-tes(DevelopmentP14TT_celltype[i,], 6,6)$g
  CellType_CohDdevelopmentP14[i, 4]<-tes(DevelopmentP14TT_celltype[i,], 6,6)$'var.g'
}

row.names(CellType_CohDdevelopmentP14)<-row.names(DevelopmentP14TT_celltype)
colnames(CellType_CohDdevelopmentP14)<-c("d P14", "var.d P14", "g P14", "var.g P14")

##Coh d for Cell Type of illumina P14 *** sample size is 6
CellType_CohDilluminaP14<-matrix(0, length(P14IlluminaCellTypeTT_forMeta.csv[,1]), 4)

for(i in 1:length(P14IlluminaCellTypeTT_forMeta.csv[,1])){
  CellType_CohDilluminaP14[i,1]<-tes(P14IlluminaCellTypeTT_forMeta.csv[i,], 6,6)$d
  CellType_CohDilluminaP14[i,2]<-tes(P14IlluminaCellTypeTT_forMeta.csv[i,], 6,6)$'var.d'
  CellType_CohDilluminaP14[i, 3]<-tes(P14IlluminaCellTypeTT_forMeta.csv[i,], 6,6)$g
  CellType_CohDilluminaP14[i, 4]<-tes(P14IlluminaCellTypeTT_forMeta.csv[i,], 6,6)$'var.g'
}

row.names(CellType_CohDilluminaP14)<-row.names(P14IlluminaCellTypeTT_forMeta.csv)
colnames(CellType_CohDilluminaP14)<-c("d illumina", "var.d illumina", "g illumina", "var.g illumina")


##Coh d Cell Type for affy230 P14 *** sample size is 6 HR and 5 LR (outlier was thrown out)
CellType_CohDAffy230<-matrix(0, length(P14Affy230CellTypeTT_forMeta.csv[,1]), 4)

for(i in 1:length(P14Affy230CellTypeTT_forMeta.csv[,1])){
  CellType_CohDAffy230[i,1]<-tes(P14Affy230CellTypeTT_forMeta.csv[i,], 6,5)$d
  CellType_CohDAffy230[i,2]<-tes(P14Affy230CellTypeTT_forMeta.csv[i,], 6,5)$'var.d'
  CellType_CohDAffy230[i, 3]<-tes(P14Affy230CellTypeTT_forMeta.csv[i,], 6,5)$g
  CellType_CohDAffy230[i, 4]<-tes(P14Affy230CellTypeTT_forMeta.csv[i,], 6,5)$'var.g'
}
row.names(CellType_CohDAffy230)<-row.names(P14Affy230CellTypeTT_forMeta.csv)
colnames(CellType_CohDAffy230)<-c("d affy", "var.d affy", "g affy", "var.g affy")


########### Meta-Analysis for P14 cell type

##Cell types (row.names) are all in same order so I can simply cbind the different datasets together

CellType_P14<-cbind.data.frame(CellType_CohDdevelopmentP14, CellType_CohDNewColonyP14, CellType_CohDAffy230, CellType_CohDilluminaP14)
colnames(CellType_P14)
#[1] "d P14"          "var.d P14"      "g P14"          "var.g P14"      "d NC P14"      
#[6] "var.d NC P14"   "g NC P14"       "var.g NC P14"   "d affy"         "var.d affy"    
#[11] "g affy"         "var.g affy"     "d illumina"     "var.d illumina" "g illumina"    
#[16] "var.g illumina"


CellType_MetaAnalysisOutputP14<-matrix(0, 10, 5)


for(i in 1:10){
  effect<-c(CellType_P14$`d P14`[i], CellType_P14$`d NC P14`[i], CellType_P14$`d affy`[i], CellType_P14$`d illumina`[i])
  var<-c(CellType_P14$`var.d P14`[i], CellType_P14$`var.d NC P14`[i], CellType_P14$`var.d affy`[i], CellType_P14$`var.d illumina`[i])
  CellType_MetaAnalysisOutputP14[i, 1]<-rma.mv(effect, var)$b #gives estimate
  CellType_MetaAnalysisOutputP14[i, 2]<-rma.mv(effect, var)$se #gives standard error
  CellType_MetaAnalysisOutputP14[i, 3]<-rma.mv(effect, var)$pval #gives pval
  CellType_MetaAnalysisOutputP14[i, 4]<-rma.mv(effect, var)$ci.lb #gives confidence interval lower bound
  CellType_MetaAnalysisOutputP14[i, 5]<-rma.mv(effect, var)$ci.ub #gives confidence interval upper bound
}

colnames(CellType_MetaAnalysisOutputP14)<-c("estimate", "SE", "pval", "CI_lb", "CI_ub")
row.names(CellType_MetaAnalysisOutputP14)<-row.names(CellType_P14)



#joining meta output with cohens d info

CellType_P14Output_Pval_CohD<-cbind.data.frame(CellType_MetaAnalysisOutputP14, CellType_P14)
write.csv(CellType_P14Output_Pval_CohD, "output/MetaAnalysis Development/CellType_P14Output_Pval_CohD.csv")

### Comparing overall P14 expression for negative and positive direction of effect
# focusing on closest to the mean/least extreme values (between 0 and 1 estimate value)

higher <- subset(MetaAnalysisOutputP14FDRnozeros, MetaAnalysisOutputP14FDRnozeros$estimate < 1 & MetaAnalysisOutputP14FDRnozeros$estimate > 0)
mean(higher$estimate)
# [1] 0.3103888

lower <- subset(MetaAnalysisOutputP14FDRnozeros, (MetaAnalysisOutputP14FDRnozeros$estimate > -1 & MetaAnalysisOutputP14FDRnozeros$estimate < 0))
mean(lower$estimate)
#[1] -0.3200704

#Pretty well balanced


#Checking extreme values
higher <- subset(MetaAnalysisOutputP14FDRnozeros, MetaAnalysisOutputP14FDRnozeros$estimate > 1)
mean(higher$estimate)
# [1] 1.24222

lower <- subset(MetaAnalysisOutputP14FDRnozeros, MetaAnalysisOutputP14FDRnozeros$estimate < -1)
mean(lower$estimate)
#[1] -1.233892

#Also balanced

######## Forest Plots for Cell Analysis of P14


for(i in 1:10){

  png(paste0(outCellPlotsP14, "Forest Plot ", paste((row.names(CellType_P14))[i]), ".png"))
  
  effect<-c(CellType_P14$`d P14`[i], CellType_P14$`d affy`[i], CellType_P14$`d illumina`[i], 
            CellType_P14$`d NC P14`[i])
  var<-c(CellType_P14$`var.d P14`[i], CellType_P14$`var.d affy`[i], CellType_P14$`var.d illumina`[i], 
         CellType_P14$`var.d NC P14`[i])
  
  forest.rma(rma.mv(effect, var), slab=c("MBNI_AffymetrixRgU34A_F6", "MBNI_AffymetrixRae230_ F15", 
                                         "MBNI_IlluminaRatRef12v1_F15", "Alabama_NimbleGen_F34"), 
             xlim=c(-17, 17))
  text(-17, 4.6, "Investigator & Study", pos=4, cex=.75)
  mtext(paste(row.names(CellType_P14)[i]), line=-1.5, cex=2)
  text(16.5, 4.6, "Cohen's D [95% CI]", pos=2, cex=.75)
  dev.off()
}




#### Meta-Analysis for P21 cell type


##Coh d for cell type development P21 *** sample size is 6
CellType_CohDdevelopmentP21<-matrix(0, length(DevelopmentP21TT_celltype[,1]), 4)

for(i in 1:length(DevelopmentP21TT_celltype[,1])){
  CellType_CohDdevelopmentP21[i,1]<-tes(DevelopmentP21TT_celltype[i,], 6,6)$d
  CellType_CohDdevelopmentP21[i,2]<-tes(DevelopmentP21TT_celltype[i,], 6,6)$'var.d'
  CellType_CohDdevelopmentP21[i, 3]<-tes(DevelopmentP21TT_celltype[i,], 6,6)$g
  CellType_CohDdevelopmentP21[i, 4]<-tes(DevelopmentP21TT_celltype[i,], 6,6)$'var.g'
}

row.names(CellType_CohDdevelopmentP21)<-row.names(DevelopmentP21TT_celltype)
colnames(CellType_CohDdevelopmentP21)<-c("d P21", "var.d P21", "g P21", "var.g P21")


##Coh d for New Colony cell type P21 *** sample size is 5

CellType_CohDNewColonyP21<-matrix(0, length(NewColonyDevelopmentP21TT_celltype[,1]), 4)

for(i in 1:length(NewColonyDevelopmentP21TT_celltype[,1])){
  CellType_CohDNewColonyP21[i,1]<-tes(NewColonyDevelopmentP21TT_celltype[i,], 5,5)$d
  CellType_CohDNewColonyP21[i,2]<-tes(NewColonyDevelopmentP21TT_celltype[i,], 5,5)$'var.d'
  CellType_CohDNewColonyP21[i, 3]<-tes(NewColonyDevelopmentP21TT_celltype[i,], 5,5)$g
  CellType_CohDNewColonyP21[i, 4]<-tes(NewColonyDevelopmentP21TT_celltype[i,], 5,5)$'var.g'
}

row.names(CellType_CohDNewColonyP21)<-row.names(NewColonyDevelopmentCellTypeTT_forMeta.csv)
colnames(CellType_CohDNewColonyP21)<-c("d NC P21", "var.d NC P21", "g NC P21", "var.g NC P21")


### P21 cell type meta-analysis
##Cell types (row.names) are all in same order so I can simply cbind the different datasets together

CellType_P21<-cbind.data.frame(cell.type=row.names(CellType_CohDdevelopmentP21), CellType_CohDdevelopmentP21, CellType_CohDNewColonyP21)
colnames(CellType_P21)
#[1] "cell.type"    "d P21"        "var.d P21"    "g P21"        "var.g P21"    "d NC P21"     "var.d NC P21"
#[8] "g NC P21"     "var.g NC P21"


CellType_MetaAnalysisOutputP21<-matrix(0, 10, 5)

for(i in 1:10){
  effect<-c(CellType_P21$`d P21`[i], CellType_P21$`d NC P21`[i])
  var<-c(CellType_P21$`var.d P21`[i], CellType_P21$`var.d NC P21`[i])
  CellType_MetaAnalysisOutputP21[i, 1]<-rma.mv(effect, var)$b #gives estimate
  CellType_MetaAnalysisOutputP21[i, 2]<-rma.mv(effect, var)$se #gives standard error
  CellType_MetaAnalysisOutputP21[i, 3]<-rma.mv(effect, var)$pval #gives pval
  CellType_MetaAnalysisOutputP21[i, 4]<-rma.mv(effect, var)$ci.lb #gives confidence interval lower bound
  CellType_MetaAnalysisOutputP21[i, 5]<-rma.mv(effect, var)$ci.ub #gives confidence interval upper bound
}

colnames(CellType_MetaAnalysisOutputP21)<-c("estimate", "SE", "pval", "CI_lb", "CI_ub")
row.names(CellType_MetaAnalysisOutputP21)<-row.names(CellType_P21)


write.csv(CellType_MetaAnalysisOutputP21, paste0(outDir, "CellType_MetaAnalysisOutputP21.csv"))


###  Forest Plots for Cell Analysis of P21


for(i in 1:10){
  
  png(paste0(outCellPlotsP21, "Forest Plot ", 
            paste((row.names(CellType_P21))[i]), ".png"))
  
  effect<-c(CellType_P21$`d P21`[i], CellType_P21$`d NC P21`[i])
  var<-c(CellType_P21$`var.d P21`[i], CellType_P21$`var.d NC P21`[i])
  
  forest.rma(rma.mv(effect, var), slab=c("MBNI_AffymetrixRgU34A_F6", 
                                         "Alabama_NimbleGen_F34"), xlim=c(-17, 17))
  text(-17, 2.6, "Investigator & Study", pos=4, cex=.75)
  mtext(paste(row.names(CellType_P21)[i]), line=-1.5, cex=2)
  text(16.5, 2.6, "Cohen's D [95% CI]", pos=2, cex=.75)
  dev.off()
}



#######


########## Generating Forest Plots for Genes of Interest

#This code produces forest plot output for various genes of interest within the developmental datasets

#P7 genes

temp1<-subset(P7, P7$gene=="C1qc" | P7$gene=="C1qb" | P7$gene=="Ncan" | P7$gene =="Bmp4")

for(i in 1:length(temp1$gene)){

  tiff(paste0(outGOIPlots, "P7 Forest Plot ", 
              paste((temp1$gene)[i]), ".tiff"), res=300, compression="lzw",
       width = 8.5, height = 8, units = 'in')
  
  effect<-c(temp1$`d P7`[i], temp1$`d new colony P7`[i])
  var<-c(temp1$`var.d P7`[i], temp1$`var.d new colony P7`[i])
  
  forest.rma(rma.mv(effect, var), slab=c("MBNI\nAffymetrixRgU34A_F6", "Alabama\nNimbleGen_F34"), 
             xlim=c(-21, 14), cex=1.6)
  text(-14, 2.7, "Investigator & Study", cex=1.6)
  text(-17, 3.5, "LR", cex=1.8)
  mtext(paste(temp1$gene[i]), line=-1.5, cex=2.5)
  text(7, 2.7, "Cohen's D [95% CI]", cex=1.6)
  text(10, 3.5, "HR", cex=1.8)
  dev.off()
}


#P14
temp1<-subset(P14cohDandPval, P14cohDandPval$gene== "Ncan"|P14cohDandPval$gene=="Bmp4"| 
                P14cohDandPval$gene== "C1qa" | P14cohDandPval$gene=="C1qb" | P14cohDandPval$gene=="C1qc")


for(i in 1:length(temp1$gene)){
  
  tiff(paste0(outGOIPlots, "P14 Forest Plot ", 
              paste((temp1$gene)[i]), ".tiff"), res=300, compression="lzw",
       width = 8.5, height = 8, units = 'in')
  
  effect<-c(temp1$`d P14`[i], temp1$`d affy`[i], temp1$`d illumina`[i], temp1$`d RNAseq`[i], temp1$`d new colony P14`[i])
  var<-c(temp1$`var.d P14`[i], temp1$`var.d affy`[i], temp1$`var.d illumina`[i], temp1$`var.d RNAseq`[i], temp1$`var.d new colony P14`[i])
  forest.rma(rma.mv(effect, var), slab=c("MBNI_AffymetrixRgU34A_F6", "MBNI_AffymetrixRae230_F15", 
                                         "MBNI_IlluminaRatRef12v1_F15","MBNI_RNASeq_F29", "Alabama_NimbleGen_F34"), 
             xlim=c(-20, 13), cex=1.6)
  
  text(-13, 5.6, "Investigator & Study", cex=1.6)
  text(-16, 6.6, "LR", cex=1.8)
  mtext(paste("P14\n", temp1$gene[i]), line=-1.5, cex=2.5)
  text(6, 5.6, "Cohen's D [95% CI]", cex=1.6)
  text(9, 6.6, "HR", cex=1.8)
  dev.off()
}


# P21
temp1<-subset(P21, P21$gene== "Ncan"|P21$gene=="Bmp4"|P21$gene== "C1qb" |P21$gene=="C1qc")

for(i in 1:length(temp1$gene)){
  tiff(paste0(outGOIPlots, "P21 Forest Plot ", 
              paste((temp1$gene)[i]), ".tiff"), res=300, compression="lzw",
       width = 8.5, height = 8, units = 'in')
  
  effect<-c(temp1$`d P21`[i], temp1$`d new colony P21`[i])
  var<-c(temp1$`var.d P21`[i], temp1$`var.d new colony P21`[i])
  
  forest.rma(rma.mv(effect, var), slab=c("MBNI\nAffymetrixRgU34A_F6", "Alabama\nNimbleGen_F34"), 
             xlim=c(-19, 17), cex=1.6)
  
  text(-12, 2.7, "Investigator & Study", cex=1.6)
  text(-15, 3.5, "LR", cex=1.8)
  mtext(paste(temp1$gene[i]), line=-1.5, cex=2.5)
  text(10, 2.7, "Cohen's D [95% CI]", cex=1.6)
  text(13, 3.5, "HR", cex=1.8)
  dev.off()
}



### Cell Types of Interest

#output for significant cell types in P14 - none will be outputted for P7 or P21
sigCellP14 <- subset(CellType_P14Output_Pval_CohD, CellType_P14Output_Pval_CohD$pval < 0.05)


for(i in 1:length(row.names(sigCellP14))){
  
  tiff(paste0(outCellPlotsP14, "Forest Plot ", (row.names(sigCellP14)[i]),
              ".tiff"), res=300, compression="lzw",
       width = 8, height = 8, units = 'in')
  
  effect<-c(sigCellP14$`d P14`[i], sigCellP14$`d affy`[i], sigCellP14$`d illumina`[i], 
            sigCellP14$`d NC P14`[i])
  var<-c(sigCellP14$`var.d P14`[i], sigCellP14$`var.d affy`[i], sigCellP14$`var.d illumina`[i], 
         sigCellP14$`var.d NC P14`[i])
  forest.rma(rma.mv(effect, var), slab=c("MBNI_AffymetrixRgU34A_F6", "MBNI_AffymetrixRae230_ F15", 
                                         "MBNI_IlluminaRatRef12v1_F15", "Alabama_NimbleGen_F34"), 
             xlim=c(-19, 13), cex=1.6)
  
  text(-13, 4.7, "Investigator & Study", cex=1.6)
  text(-15, 5.5, "LR", cex=1.8)
  mtext(paste("P14\n", row.names(sigCellP14)[i]), line=-1.5, cex=2.4)
  text(7, 4.7, "Cohen's D [95% CI]", cex=1.6)
  text(9, 5.5, "HR", cex=1.8)
  dev.off()
}


##########################################

### Heatmap of Cell Type Results

#Need dataframe of all betas for cell types - rows cell types, columns are age groups

#Reading in cell type meta output from adult analysis
cellTypeAdult <- read.csv("output/MetaAnalysis Adult/Adult_CellTypeMeta_CohDandPval.csv")
colnames(cellTypeAdult)


#Creating data frame of cell type estimated effects (Betas)

colnames(CellType_MetaAnalysisOutputP7) #same order for all
#[1] "estimate" "SE"       "pval"     "CI_lb"    "CI_ub" 

cellTypeBetas <- cbind.data.frame(CellType_MetaAnalysisOutputP7[,1], CellType_MetaAnalysisOutputP14[,1],
                                  CellType_MetaAnalysisOutputP21[,1], cellTypeAdult$estimate)

colnames(cellTypeBetas) <- c("P7", "P14", "P21", "Adult")

grep("Neuron", row.names(cellTypeBetas)) #rows to remove - do not want neuron cell types

cellTypeBetas <- rbind.data.frame(cellTypeBetas[c(1:4, 8:10),])

dim(cellTypeBetas)

#Outputting high resolution image of cell type heatmap (neuron types not included)

tiff(paste0(outDir, "CellTypeHeatmap.tiff"), 
     width = 16, height = 6, units = 'in', res = 300, compression = "lzw")

par(mar=c(5,16,4,5))
image(t(as.matrix(cellTypeBetas)), axes=F, col=bluered(100))
axis(2, at=c(0, .166, .33, .5, .67, .84, 1), labels=row.names(cellTypeBetas), las=1, cex.axis=1.4)
axis(1, at=c(0, .33, .67, 1), labels=colnames(cellTypeBetas), las=2, cex.axis=1.5)

# Add * symbol for significant cell typse
points(1, .33, pch=8, cex=2) #adult microglia
points(.33, .33, pch=8, cex=2) #P14 microglia
points(.33, .17, pch=8, cex=2) #P14 endothelial
mtext("Cell Type Estimate by Age", cex=2.2, line=1)

#Adding color gradient legend
image.plot(t(as.matrix(cellTypeBetas)), axes=F, col=bluered(100), 
           legend.only=T, axis.args = list(cex.axis = 1.5))

dev.off()



#########################################   Volcano Plots


## Generating Volcano plots
# Observing HR/LR overall expression from meta-analysis results

###P7 

colnames(MetaAnalysisOutputP7FDR)
#[1] "gene"     "rawp"     "BH"       "BY"       "estimate" "SE"       "pval"     "CI_lb"    "CI_ub" 


#Volcano Plot overall P7 expression

tiff(paste0(outDir, "VolcanoPlotP7.tiff"), width = 5, height = 5, units = 'in', 
            res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(MetaAnalysisOutputP7FDR, plot(estimate, -log10(pval), pch=19, main="Overall Expression", 
                                   xlim=c(-4,3), cex.lab=1.8, cex.main=2, cex=0.6))

# Add colored points: red if padj<0.05, blue of estimate>1)
with(subset(MetaAnalysisOutputP7FDR, abs(estimate)>1), points(estimate, -log10(pval), 
                                                                      pch=19, col="red", cex=0.6))
with(subset(MetaAnalysisOutputP7FDR, BH< .05 ), points(estimate, -log10(pval), 
                                                               pch=19, col="blue", cex=0.6))
legend(-1.5, 7, legend=c("estimate > 1", "FDR < 0.05"), col=c("red", "blue"), pch=19, cex=1.2)

dev.off()

##Overall P7 expression of cell type specific genes

cellInfo <- read.csv("data/CellTypeSpecificGenes_Master3.csv")
colnames(cellInfo)[5]<-"gene"

temp<-join(MetaAnalysisOutputP7FDR, cellInfo, by="gene", type="inner")

temp<-subset(temp, !duplicated(temp$gene))

colnames(temp)

tiff(paste0(outDir, "VolcanoPlot CellTypeP7.tiff"), height=5, 
     width = 5, units="in", res=300, compression="lzw")
par(mai=c(1.02, 1,0.9,0.50))

with(temp, plot(estimate, -log10(pval), pch=19, main="P7 Overall Expression\nof Cell Type Specific Genes", 
                xlim=c(-4,3), cex.lab=1.8, cex.main=2, cex=0.6))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(temp, abs(estimate)>1), points(estimate, -log10(pval), pch=19, col="red", cex=0.6))
with(subset(temp, BH<.05 ), points(estimate, -log10(pval), pch=19, col="blue", cex=0.6))
legend(-1.5, 7, legend=c("estimate > 1", "FDR < 0.05"), col=c("red", "blue"), pch=19, cex=1.2)

dev.off()



### P14

colnames(MetaAnalysisOutputP14FDRnozeros)
#[1] "gene"     "estimate" "SE"       "pval"     "CI_lb"    "CI_ub"    "Datasets" "rawp"     "BH"      
#[10] "BY" 


#Volcano Plot overall P14 expression

tiff(paste0(outDir, "VolcanoPlotP14.tiff"), width = 5, height = 5, 
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(MetaAnalysisOutputP14FDRnozeros, plot(estimate, -log10(pval), pch=19, main="Overall Expression", 
                                 xlim=c(-3,3), ylim=c(0, 5.5), cex.lab=1.8, cex.main=2, cex=0.6))
# Add colored points: red if padj<0.05, blue of estimate>1)
with(subset(MetaAnalysisOutputP14FDRnozeros, abs(estimate)>1), points(estimate, -log10(pval), 
                                                            pch=19, col="red", cex=0.6))
with(subset(MetaAnalysisOutputP14FDRnozeros, BH< .05 ), points(estimate, -log10(pval), 
                                                     pch=19, col="blue", cex=0.6))
legend(-1.2, 5.499, legend=c("estimate > 1", "FDR < 0.05"), col=c("red", "blue"), pch=19, cex=1.2)

dev.off()


#Overall P14 expression colored by dataset

tiff(paste0(outDir, "VolcanoPlotP14_byDatasets.tiff"), width = 5, height = 5, 
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))

with(MetaAnalysisOutputP14FDRnozeros, plot(estimate, -log10(pval), pch=19, main="Overall Expression by 
     Dataset Number", xlim=c(-3,3), ylim=c(0, 5.5), cex.lab=1.8, cex.main=2, cex=0.6))

#color by dataset
with(subset(MetaAnalysisOutputP14FDRnozeros, Datasets > 4), points(estimate, -log10(pval), pch=19, col="blue", cex=0.6))
with(subset(MetaAnalysisOutputP14FDRnozeros, Datasets < 5 & Datasets > 3), points(estimate, -log10(pval), pch=19, col="red", cex=0.6))
with(subset(MetaAnalysisOutputP14FDRnozeros, Datasets < 4 & Datasets > 2), points(estimate, -log10(pval), pch=19, col="green3", cex=0.6))
with(subset(MetaAnalysisOutputP14FDRnozeros, Datasets < 3), points(estimate, -log10(pval), pch=19, col="purple", cex=0.6))

legend(-.5, 5.499, legend=c("2", "3", "4", "5"), 
       col=c("purple", "green3", "red", "blue"), pch=19, cex=1.2)


dev.off()


### checking overall expression of P14 cell type specific genes

temp<-join(MetaAnalysisOutputP14FDRnozeros, cellInfo, by="gene", type="inner")

temp<-subset(temp, !duplicated(temp$gene))

colnames(temp)

tiff(paste0(outDir, "VolcanoPlot CellTypeP14.tiff"), height=5, 
     width = 5, units="in", res=300, compression="lzw")
par(mai=c(1.02, 1,0.9,0.50))

with(temp, plot(estimate, -log10(pval), pch=19, main="P14 Overall Expression\nof Cell Type Specific Genes",
                xlim=c(-3,3), cex.lab=1.8, cex.main=2, cex=0.6))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(temp, abs(estimate)>1), points(estimate, -log10(pval), pch=19, col="red", cex=0.6))
with(subset(temp, BH<.05 ), points(estimate, -log10(pval), pch=19, col="blue", cex=0.6))
legend(-1.2, 5, legend=c("estimate > 1", "FDR < 0.05"), col=c("red", "blue"), pch=19, cex=1.2)

dev.off()



###P21

colnames(MetaAnalysisOutputP21FDR)
#[1] "gene"     "rawp"     "BH"       "BY"       "estimate" "SE"       "pval"     "CI_lb"    "CI_ub" 

#Volcano Plot overall P21 expression
tiff(paste0(outDir, "VolcanoPlotP21.tiff"), width = 5, height = 5, units = 'in', 
     res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))

with(MetaAnalysisOutputP21FDR, plot(estimate, -log10(pval), pch=19, main="Overall Expression", 
                                   xlim=c(-3,3), ylim=c(0, 5.5), cex.lab=1.8, cex.main=2, cex=0.6))
# Add colored points: red if padj<0.05, blue of estimate>1)
with(subset(MetaAnalysisOutputP21FDR, abs(estimate)>1), points(estimate, -log10(pval), 
                                                              pch=19, col="red", cex=0.6))
with(subset(MetaAnalysisOutputP21FDR, BH< .05 ), points(estimate, -log10(pval), 
                                                       pch=19, col="blue", cex=0.6))
legend(-1.5, 5.499, legend=c("estimate > 1", "FDR < 0.05"), col=c("red", "blue"), pch=19, cex=1.2)

dev.off()



##Overall P21 expression of cell type specific genes

temp<-join(MetaAnalysisOutputP21FDR, cellInfo, by="gene", type="inner")

temp<-subset(temp, !duplicated(temp$gene))

colnames(temp)

tiff(paste0(outDir, "VolcanoPlot CellTypeP21.tiff"), height=5, width = 5, 
     units="in", res=300, compression="lzw")
par(mai=c(1.02, 1,0.9,0.50))

with(temp, plot(estimate, -log10(pval), pch=19, main="P21 Overall Expression\nof Cell Type Specific Genes",
                xlim=c(-3,3), cex.lab=1.8, cex.main=2, cex=0.6))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(temp, BH<.05 ), points(estimate, -log10(pval), pch=19, col="blue", cex=0.6))
with(subset(temp, abs(estimate)>1), points(estimate, -log10(pval), pch=19, col="red", cex=0.6))
legend(-1.5, 4.2, legend=c("estimate > 1", "FDR < 0.05"), col=c("red", "blue"), pch=19, cex=1.2)

dev.off()





