#Installing fgsea package
source("https://bioconductor.org/biocLite.R")
biocLite("fgsea")
library("fgsea")

packageVersion("fgsea")


#Reading in meta-analysis results to make ranked gene list
setwd("~/Phenotype Project/ibirt")

adultMeta <- read.csv("output/MetaAnalysis Adult/AdultMetaAnalysisOutputFDR.csv")
colnames(adultMeta)
#[1] "X"          "rawp"       "BH"         "BY"         "GeneSymbol" "estimate"   "SE"        
#[8] "pval"       "CI_lb"      "CI_ub"      "datasets"  


## Reading in .gmt file 

#Using biological processes gene sets
rat_Pathways<-gmtPathways("Read-in files/gsea/Rattus_norvegicus_GSEA_GO_sets_bp_symbols_highquality_April_2015.gmt")
names(rat_Pathways)

rat_Pathways[[1]]
# [1] "AADAT" "DLD"   "DLST"  "GOT1"  "GOT2"  "IDH1"  "IDH2"  "IDH3A" "IDH3B" "IDH3G" "OGDH" 

#They're all caps (despite being rat gene symbols) so will have to change gene symbols in the ranked list to all caps


#creating ranked list based on estimate
adultrnklist <- sort(adultMeta$estimate, decreasing=TRUE)
adultrnknames <- order(adultMeta$estimate, decreasing=TRUE)

#naming ranked gene list with uppercase gene symbols to be used in fgsea
names(adultrnklist) <- toupper(adultMeta$GeneSymbol[adultrnknames])

#making data frame of ranked genes for output
adultRnk<-cbind.data.frame(GeneSymbol=toupper(adultMeta$GeneSymbol[adultrnknames]), estimate=adultrnklist)

colnames(adultRnk)

#exporting ranked list
write.table(adultRnk, file="output/FGSEA/adultmetaUPPERCASE_rankedlist.rnk", 
            sep="\t", row.names=FALSE, col.names = FALSE, quote = FALSE)


### FGSEA

#Running fgsea with biological processes GO sets .gmt

gseaMetaGenes<-fgsea(rat_Pathways, adultrnklist, nperm=1000, minSize = 1, maxSize = 500)

#condensing leadingEdge list so that it can be exported
gseaMetaGenes$leadingEdge<-vapply(gseaMetaGenes$leadingEdge, paste, collapse = ",", character(1L))
write.csv(gseaMetaGenes, "output/FGSEA/adultmetagene_FGSEA_1000perm.csv")

#Running with 100,000 perm 

gseaMetaGeneshun<-fgsea(rat_Pathways, adultrnklist, nperm=100000, minSize = 1, maxSize = 500)

gseaMetaGeneshun$leadingEdge<-vapply(gseaMetaGeneshun$leadingEdge, paste, collapse = ",", character(1L))
write.csv(gseaMetaGeneshun, "output/FGSEA/adultmetagene_FGSEA_100,000perm.csv")



### P14 FGSEA

P14Meta <- read.csv("output/MetaAnalysis Development/MetaAnalysisOutputP14FDR.csv")
colnames(P14Meta)
#[1] "X"        "gene"     "estimate" "SE"       "pval"     "CI_lb"    "CI_ub"    "Datasets" "rawp"    
#[10] "BH"       "BY"   

#creating ranked list based on estimate
P14rnklist <- sort(P14Meta$estimate, decreasing=TRUE)
P14rnknames <- order(P14Meta$estimate, decreasing=TRUE)

#naming ranked gene list with uppercase gene symbols to be used in fgsea
names(P14rnklist) <- toupper(P14Meta$gene[P14rnknames])


#Running with 100,000 perm 

gseaP14MetaGenes<-fgsea(rat_Pathways, P14rnklist, nperm=100000, minSize = 1, maxSize = 500)

gseaP14MetaGenes$leadingEdge<-vapply(gseaP14MetaGenes$leadingEdge, paste, collapse = ",", character(1L))
write.csv(gseaP14MetaGenes, "output/FGSEA/P14metagene_FGSEA_100,000perm.csv")





##
### Visualizing gsea output


##Creating table of top pathways for adult

topPathways <- gseaMetaGeneshun[head(order(pval), n=15)][order(pval), pathway]
## 

png("output/FGSEA/MetaAdult_TopPathways.png", width=1000)
plotGseaTable(rat_Pathways[topPathways], adultrnklist,
              gseaMetaGeneshun, gseaParam=0.5)
dev.off()

#looking at top pathway
(topPathways[1])
#[1] "protein_folding(6)"


#plotting top adult fgsea pathway

png("output/FGSEA/Protein Folding Pathway Enrichment for Adult.png")
plotEnrichment(rat_Pathways[["protein_folding(6)"]], adultrnklist)
dev.off()

#Visualizing gsea output


##Creating table of top pathways for P14

topP14Pathways <- gseaP14MetaGenes[head(order(pval), n=15)][order(pval), pathway]
## 

png("output/FGSEA/MetaP14_TopPathways.png", width=1000)
plotGseaTable(rat_Pathways[topP14Pathways], P14rnklist,
              gseaP14MetaGenes, gseaParam=0.5)
dev.off()


#Getting top pathway
(topP14Pathways[1])
#[1] "neuron_differentiation(6)"


#plotting top P14 pathway

png("output/FGSEA/Neuron Differentiation Pathway Enrichment for P14.png")
plotEnrichment(rat_Pathways[["neuron_differentiation(6)"]])
dev.off()






########################### ***********************************************


## Creating Custom gene set files


#Reading in Module information .csv files

humanCoexpressionModules<-read.csv("Read-in files/gsea/HumanCoexpressionModules.csv")
colnames(humanCoexpressionModules)
#[1] "X"                 "Illumina.Probe.ID" "Ensembl.Gene.ID"   "GeneSymbol_Human"  "Gene.Description" 
#[6] "Chromosome"        "Gene.Start..bp."   "Gene.End..bp."     "Module"  

mouseCoexpressionModules<-read.csv("Read-in files/gsea/MouseCoexpressionModules.csv")
colnames(mouseCoexpressionModules)
#[1] "X"                   "Search_Key"          "Transcript"          "GeneSymbol_Human"    "Source_Reference_ID"
#[6] "RefSeq_ID"           "Probe_Id"            "cis_minP_hippo"      "cis_minP_striatum"   "chromosome"         
#[11] "chrStart"            "moduleHippo"         "moduleStriatum"      "kTotalHippo"         "kWithinHippo"       
#[16] "kOutHippo"           "kDiffHippo"          "kTotalStriatum"      "kWithinStriatum"     "kOutStriatum"       
#[21] "kDiffStriatum" 


###
#Making list of genes according to module

humanModules<-split(humanCoexpressionModules$GeneSymbol_Human, humanCoexpressionModules$Module)


#splitting according to hippocampal modules
mouseModules<-split(mouseCoexpressionModules$GeneSymbol_Human, mouseCoexpressionModules$moduleHippo)


#Making gmt files for both

#human gmt
str(humanModules)
names(humanModules[1])
humanModules[1]


humanModulesGMT<-matrix(0, 25, 3)
for(i in 1:25){
  humanModulesGMT[i,1]<-names(humanModules[i]) #module names
  humanModulesGMT[,2]<-"NA" #dummy row
  humanModulesGMT[i,3]<-vapply(humanModules[i], paste, collapse="\t", character(1L)) #gene list per module
}

colnames(humanModulesGMT)<-c("Modules", "Dummy", "Genes")
write.table(humanModulesGMT, "output/FGSEA/humanCoexpressionModulesGMT.gmt", col.names = FALSE, row.names = FALSE, sep="\t")


#mouse gmt
length(mouseModules)
#[1] 30

mouseModulesGMT<-matrix(0, 30, 3)
for(i in 1:30){
  mouseModulesGMT[i,1]<-names(mouseModules[i]) #module names
  mouseModulesGMT[,2]<-"NA" #dummy row
  mouseModulesGMT[i,3]<-vapply(mouseModules[i], paste, collapse="\t", character(1L)) #gene list per module
}

colnames(mouseModulesGMT)<-c("Mouse Modules", "Dummy", "Genes")
write.table(mouseModulesGMT, "output/FGSEA/mouseCoexpressionModulesGMT.gmt", col.names = FALSE, row.names = FALSE, sep="\t")




###################################################################

# Marker Genes

dataDir<-("Read-in files/gsea/")

dendrogramMarkers<-read.csv(paste0(dataDir, "HippoSeq_DendrogramMarkerGenes2.csv"))
colnames(dendrogramMarkers)
#[1] "X"                    "enriched"             "depleted"             "gene_id"             
#[5] "GeneSymbol_Mouse"     "uniprot_protein_name" "coronalIshAvailable"  "coronalIshConsistent"
#[9] "GeneSymbol_Human"     "EnrichedVSDepleted"

dgGranuleMarkers<-read.csv(paste0(dataDir, "HippoSeq_DGGranuleCellMarkerGenes.csv"))
colnames(dgGranuleMarkers)
#[1] "X"                    "gene_id"              "GeneSymbol_Mouse"     "uniprot_protein_name" "dg_d_fpkm"           
#[6] "dg_v_fpkm"            "ratio"                "coronalIshExists"     "coronalIshCorrect"  

dgMossyMarkers<-read.csv(paste0(dataDir, "HippoSeq_DGMossyCellMarkerGenes.csv"))
colnames(dgMossyMarkers)
#[1] "X"                    "enriched"             "gene_id"              "GeneSymbol_Mouse"     "uniprot_protein_name"
#[6] "novel"                "coronalIshAvailable"  "coronalIshConsistent"


###################################
#Extracting enriched vs depleted for dendrogram makers with gene lists for each value

colnames(dendrogramMarkers)
#[1] "X"                    "enriched"             "depleted"            
#[4] "gene_id"              "GeneSymbol_Mouse"     "uniprot_protein_name"
#[7] "coronalIshAvailable"  "coronalIshConsistent" "GeneSymbol_Human"    
#[10] "EnrichedVSDepleted" 


#Splitting genes according to which module they fall in (enriched vs depleted category)
dendrogramEvsD_genes<-split(dendrogramMarkers$GeneSymbol_Human, dendrogramMarkers$EnrichedVSDepleted)
length(dendrogramEvsD_genes)
#[1] 12

#Splitting genes using mouse gene symbols (first letter uppercase rest lowercase)
dendrogramEvsD_mousegenes<-split(dendrogramMarkers$GeneSymbol_Mouse, dendrogramMarkers$EnrichedVSDepleted)
length(dendrogramEvsD_mousegenes)
#[1] 12

#Condensing genes to list for both dg granule and dg mossy files

mossy<-split(toupper(dgMossyMarkers$GeneSymbol_Mouse), "Mossy")
(mossy[[1]])[1:5]
#[1] "6330406I15RIK" "AJAP1"         "AP2S1"         "ASS1"          "B230216N24RIK"

granule<-split(toupper(dgGranuleMarkers$GeneSymbol_Mouse), "Granule")
(granule[[1]])[1:5]
#[1] "AKAP5"  "DES"    "TNNT2"  "TAGLN2" "IFITM2"

#######################

#  making dendrogram gmt with granule and mossy cell markers

dendrogramandCellTypeMarkers<-matrix(0, 14, 3)
for(i in 1:12){
  dendrogramandCellTypeMarkers[i,1]<-names(dendrogramEvsD_genes[i]) #enriched versus depleted names
  dendrogramandCellTypeMarkers[13,1]<-names(mossy)
  dendrogramandCellTypeMarkers[14,1]<-names(granule)
  dendrogramandCellTypeMarkers[(1:14),2]<-"NA"  #dummy column
  dendrogramandCellTypeMarkers[i,3]<-vapply(dendrogramEvsD_genes[i], paste, collapse="\t", character(1L)) #gene list per module
  dendrogramandCellTypeMarkers[13,3]<-vapply(mossy, paste, collapse="\t", character(1L))
  dendrogramandCellTypeMarkers[14,3]<-vapply(granule, paste, collapse="\t", character(1L))
}

colnames(dendrogramandCellTypeMarkers)<-c("EnrichedVSDepleted", "dummy", "Genes")
#
write.table(dendrogramandCellTypeMarkers, "output/FGSEA/dendrogram_CellTypeMarkers.gmt", col.names = FALSE, row.names = FALSE, sep="\t")



##################################

#              Running fgsea with all homemade gmts
#

outDir <- ("output/FGSEA/")

#Mouse Coexpression
homemadeMouseCoexpression<-gmtPathways("output/FGSEA/mouseCoexpressionModulesGMT.gmt")

adultMetaGenes_MouseCoexpressiongsea<-fgsea(homemadeMouseCoexpression, adultrnklist, nperm=100000, minSize = 1, maxSize = 1000)

adultMetaGenes_MouseCoexpressiongsea$leadingEdge<-vapply(adultMetaGenes_MouseCoexpressiongsea$leadingEdge, paste, collapse= ",", character(1L))
write.csv(adultMetaGenes_MouseCoexpressiongsea, paste0(outDir, "fgsea_Adult_MouseCoexpression.csv"))


#Human Coexpression
homemadeHumanCoexpression<-gmtPathways("output/FGSEA/humanCoexpressionModulesGMT.gmt")

adultMetaGenes_HumanCoexpressiongsea<-fgsea(homemadeHumanCoexpression, adultrnklist, nperm=100000, minSize = 1, maxSize = 1000)


sigadultHuman<-subset(adultMetaGenes_HumanCoexpressiongsea, adultMetaGenes_HumanCoexpressiongsea$pval < 0.05)

adultMetaGenes_HumanCoexpressiongsea$leadingEdge<-vapply(adultMetaGenes_HumanCoexpressiongsea$leadingEdge, paste, collapse= ",", character(1L))
write.csv(adultMetaGenes_HumanCoexpressiongsea, paste0(outDir, "fgsea_Adult_HumanCoexpression.csv"))



#Dendrogram Cell Markers Encrichment vs Depletion
homemadeDendrogram_CellMarkers<-gmtPathways("output/FGSEA/dendrogram_CellTypeMarkers.gmt")

adultMetaGenes_Dendrogram_CellMarkersgsea<-fgsea(homemadeDendrogram_CellMarkers, adultrnklist, nperm=100000, minSize = 1, maxSize = 1000)

adultMetaGenes_Dendrogram_CellMarkersgsea$leadingEdge<-vapply(adultMetaGenes_Dendrogram_CellMarkersgsea$leadingEdge, paste, collapse= ",", character(1L))
write.csv(adultMetaGenes_Dendrogram_CellMarkersgsea, paste0(outDir, "FGSEA_AdultMetaGenes_DendrogramCellMarkers.csv"))




########## Making tables of top coexpressin modules/cell enrichment for each fgsea output


topMarkers <- adultMetaGenes_Dendrogram_CellMarkersgsea[head(order(pval), n=15)][order(pval), pathway]
## 

png(paste0(outDir, "MetaAdult_TopCellMarkerEnrichmentvsDepletion.png"), width=800)
plotGseaTable(homemadeDendrogram_CellMarkers[topMarkers], adultrnklist,
              adultMetaGenes_Dendrogram_CellMarkersgsea, gseaParam=0.5)
dev.off()



######## Mouse coexpression

topMouseCoexpression <- adultMetaGenes_MouseCoexpressiongsea[head(order(pval), n=15)][order(pval), pathway]
## 

png(paste0(outDir, "MetaAdult_TopMouseCoexpressionModules.png"))
plotGseaTable(homemadeMouseCoexpression[topMouseCoexpression], adultrnklist,
              adultMetaGenes_MouseCoexpressiongsea, gseaParam=0.5)
dev.off()


######## Human coexpression

topHumanCoexpression <- adultMetaGenes_HumanCoexpressiongsea[head(order(pval), n=15)][order(pval), pathway]
## 

png(paste0(outDir, "MetaAdult_TopHumanCoexpressionModules.png"))
plotGseaTable(homemadeHumanCoexpression[topHumanCoexpression], adultrnklist,
              adultMetaGenes_HumanCoexpressiongsea, gseaParam=0.5)
dev.off()





###################################################################################

######## Running homemade .gmts with P14 meta genes

str(P14rnklist)


#Mouse Coexpression

P14MetaGenes_MouseCoexpressiongsea<-fgsea(homemadeMouseCoexpression, P14rnklist, nperm=100000, minSize = 1, maxSize = 1000)

P14MetaGenes_MouseCoexpressiongsea$leadingEdge<-vapply(P14MetaGenes_MouseCoexpressiongsea$leadingEdge, paste, collapse= ",", character(1L))
write.csv(P14MetaGenes_MouseCoexpressiongsea, paste0(outDir, "fgsea_P14_MouseCoexpression.csv"))


#Human Coexpression
P14MetaGenes_HumanCoexpressiongsea<-fgsea(homemadeHumanCoexpression, P14rnklist, nperm=100000, minSize = 1, maxSize = 1000)

P14MetaGenes_HumanCoexpressiongsea$leadingEdge<-vapply(P14MetaGenes_HumanCoexpressiongsea$leadingEdge, paste, collapse= ",", character(1L))
write.csv(P14MetaGenes_HumanCoexpressiongsea, paste0(outDir, "fgsea_P14_HumanCoexpression.csv"))


#Dendrogram Cell Markers/region Encrichment vs Depletion
P14MetaGenes_Dendrogram_CellMarkersgsea<-fgsea(homemadeDendrogram_CellMarkers, P14rnklist, nperm=100000, minSize = 1, maxSize = 1000)

P14MetaGenes_Dendrogram_CellMarkersgsea$leadingEdge<-vapply(P14MetaGenes_Dendrogram_CellMarkersgsea$leadingEdge, paste, collapse= ",", character(1L))
write.csv(P14MetaGenes_Dendrogram_CellMarkersgsea, paste0(outDir, "fgsea_P14MetaGenes_DendrogramCellMarkers.csv"))

######

############# Creating Tables for P14 homemade gmt fgsea output

#Subsetting most significant pathways and ordering by NES
topMarkersP14 <- P14MetaGenes_Dendrogram_CellMarkersgsea[head(order(pval), n=15)][order(pval), pathway]
## 

png(paste0(outDir, "MetaP14_TopCellMarkerEnrichmentvsDepletion.png"), width=800)
plotGseaTable(homemadeDendrogram_CellMarkers[topMarkersP14], P14rnklist,
              P14MetaGenes_Dendrogram_CellMarkersgsea, gseaParam=0.5)
dev.off()



######## Mouse coexpression

topMouseCoexpressionP14 <- P14MetaGenes_MouseCoexpressiongsea[head(order(pval), n=15)][order(pval), pathway]
## 

png(paste0(outDir, "MetaP14_TopMouseCoexpressionModules.png"))
plotGseaTable(homemadeMouseCoexpression[topMouseCoexpressionP14], P14rnklist,
              P14MetaGenes_MouseCoexpressiongsea, gseaParam=0.5)
dev.off()


######## Human coexpression
topHumanCoexpression <- P14MetaGenes_HumanCoexpressiongsea[head(order(pval), n=15)][order(pval), pathway]
 

png(paste0(outDir, "MetaP14_TopHumanCoexpressionModules.png"))
plotGseaTable(homemadeHumanCoexpression[topHumanCoexpression], P14rnklist,
              P14MetaGenes_HumanCoexpressiongsea, gseaParam=0.5)
dev.off()





##################### Observing overlap between gene ontology and module/expression gene sets



#        ()_()
#        =\"/=   Mouse Coexpression/Rat pathway overlap
#          O


temp<-vector("list", length(rat_Pathways))
intersectionMouseCoexp_RatPathway<-matrix(0, length(rat_Pathways), 30)

for(i in 1:length(rat_Pathways)){
  for(l in 1:length(homemadeMouseCoexpression)){
    temp[[i]]<-(intersect(homemadeMouseCoexpression[[l]], rat_Pathways[[i]]))
    intersectionMouseCoexp_RatPathway[i,l]<-vapply(temp[i], paste, collapse=",", character(1L))
  }  
}

colnames(intersectionMouseCoexp_RatPathway)<-names(homemadeMouseCoexpression)
row.names(intersectionMouseCoexp_RatPathway)<-names(rat_Pathways)

write.csv(intersectionMouseCoexp_RatPathway, paste0(outDir, "intersectionMouseCoexp_RatGOPathways.csv"))

##### Making dataframe of sums of overlapping gene

intersectionMouseCoexp_RatPathway_SUMS<-matrix(0, length(rat_Pathways), length(homemadeMouseCoexpression))

for(i in 1:length(rat_Pathways)){
  for(l in 1:length(homemadeMouseCoexpression)){
    intersectionMouseCoexp_RatPathway_SUMS[i,l]<-sum(rat_Pathways[[i]]%in%homemadeMouseCoexpression[[l]])
  }  
}

colnames(intersectionMouseCoexp_RatPathway_SUMS)<-names(homemadeMouseCoexpression)
row.names(intersectionMouseCoexp_RatPathway_SUMS)<-names(rat_Pathways)

write.csv(intersectionMouseCoexp_RatPathway_SUMS, paste0(outDir, "SUMSintersectionMouseCoexp_RatGOPathways.csv"))


#################################################################################
#  Now running for Human co-expression modules


temp<-vector("list", length(rat_Pathways))
intersectionHumanCoexp_RatPathway<-matrix(0, length(rat_Pathways), length(homemadeHumanCoexpression))

for(i in 1:length(rat_Pathways)){
  for(l in 1:length(homemadeHumanCoexpression)){
    temp[[i]]<-(intersect(homemadeHumanCoexpression[[l]], rat_Pathways[[i]]))
    intersectionHumanCoexp_RatPathway[i,l]<-vapply(temp[i], paste, collapse=",", character(1L))
  }  
}

colnames(intersectionHumanCoexp_RatPathway)<-names(homemadeHumanCoexpression)
row.names(intersectionHumanCoexp_RatPathway)<-names(rat_Pathways)

write.csv(intersectionHumanCoexp_RatPathway, paste0(outDir, "intersectionHumanCoexp_RatGOPathways.csv"))

######### Making dataframe of sums of overlapping gene

intersectionHumanCoexp_RatPathway_SUMS<-matrix(0, length(rat_Pathways), length(homemadeHumanCoexpression))

for(i in 1:length(rat_Pathways)){
  for(l in 1:length(homemadeHumanCoexpression)){
    intersectionHumanCoexp_RatPathway_SUMS[i,l]<-sum(rat_Pathways[[i]]%in%homemadeHumanCoexpression[[l]])
  }  
}

colnames(intersectionHumanCoexp_RatPathway_SUMS)<-names(homemadeHumanCoexpression)
row.names(intersectionHumanCoexp_RatPathway_SUMS)<-names(rat_Pathways)

write.csv(intersectionHumanCoexp_RatPathway_SUMS, paste0(outDir, "SUMSintersectionHumanCoexp_RatGOPathways.csv"))

                                           

###############################
#number of genes present overlapping between cell marker/region enrichment vs rat pathway sum

intersectionHippoSeqModules_RatPathway_SUMS<-matrix(0, length(rat_Pathways), length(homemadeDendrogram_CellMarkers))

for(i in 1:length(rat_Pathways)){
  for(l in 1:length(homemadeDendrogram_CellMarkers)){
    intersectionHippoSeqModules_RatPathway_SUMS[i,l]<-sum(rat_Pathways[[i]]%in%homemadeDendrogram_CellMarkers[[l]])
  }  
}

colnames(intersectionHippoSeqModules_RatPathway_SUMS)<-names(homemadeDendrogram_CellMarkers)
row.names(intersectionHippoSeqModules_RatPathway_SUMS)<-names(rat_Pathways)
write.csv(intersectionHippoSeqModules_RatPathway_SUMS, paste0(outDir, "SUMSintersectionHippoSeqModules_RatGOPathways.csv"))




###################################################


############ Creating fgsea figure for manuscript

#GO fgsea sig pathways tables
FDRP14Pathways <- gseaP14MetaGenes[which(gseaP14MetaGenes$padj < 0.1)][order(pval), pathway]
FDRAdultPathways <- gseaMetaGeneshun[which(gseaMetaGeneshun$padj < 0.1)][order(pval), pathway]


tiff("output/FGSEA/P14topGSEApathways.tiff", width=11.8, height=6.5, units = 'in', res = 300, compression = "lzw")
plotGseaTable(rat_Pathways[FDRP14Pathways], P14rnklist,
              gseaP14MetaGenes, gseaParam=0.5, colwidths = c(5.2, 3, 0.6, .8, .8))
dev.off()


tiff("output/FGSEA/AdulttopGSEApathways.tiff", width=9, height=3.5, units = 'in', res = 300, compression = "lzw")
plotGseaTable(rat_Pathways[FDRAdultPathways], adultrnklist,
              gseaMetaGeneshun, gseaParam=0.5)

dev.off()



###Custom GMT fgsea sig pathways tables for manuscript


## Mouse coexpression Adult

topMarkers <- adultMetaGenes_MouseCoexpressiongsea[which(adultMetaGenes_MouseCoexpressiongsea$padj < 0.1)][order(pval), pathway]

tiff("output/FGSEA/AdulttopMouseCoexpressionmodules.tiff", width=6, 
     height=4, units = 'in', res = 300, compression = "lzw")
plotGseaTable(homemadeMouseCoexpression[topMarkers], adultrnklist,
              adultMetaGenes_MouseCoexpressiongsea, gseaParam=0.5, colwidths = c(2.5, 3, 0.8, 1, 1))
dev.off()


##  Mouse coexpression P14

topMarkers <- P14MetaGenes_MouseCoexpressiongsea[which(P14MetaGenes_MouseCoexpressiongsea$padj < 0.1)][order(pval), pathway]

tiff("output/FGSEA/P14topMouseCoexpressionmodules.tiff", width=6, 
     height=4, units = 'in', res = 300, compression = "lzw")
plotGseaTable(homemadeMouseCoexpression[topMarkers], P14rnklist,
              P14MetaGenes_MouseCoexpressiongsea, gseaParam=0.5, colwidths = c(2.5, 3, 0.8, 1, 1))
dev.off()



## Human coexpression Adult

topMarkers <- adultMetaGenes_HumanCoexpressiongsea[which(adultMetaGenes_HumanCoexpressiongsea$padj < 0.1)][order(pval), pathway]

tiff("output/FGSEA/AdulttopHumanCoexpressionmodules.tiff", width=6, 
     height=2, units = 'in', res = 300, compression = "lzw")
plotGseaTable(homemadeHumanCoexpression[topMarkers], adultrnklist,
              adultMetaGenes_HumanCoexpressiongsea, gseaParam=0.5, colwidths = c(2, 3, 0.8, 1, 1))
dev.off()


## Human coexpression P14

topMarkers <- P14MetaGenes_HumanCoexpressiongsea[which(P14MetaGenes_HumanCoexpressiongsea$padj < 0.1)][order(pval), pathway]

tiff("output/FGSEA/P14topHumanCoexpressionmodules.tiff", width=6, 
     height=2, units = 'in', res = 300, compression = "lzw")
plotGseaTable(homemadeHumanCoexpression[topMarkers], P14rnklist,
              P14MetaGenes_HumanCoexpressiongsea, gseaParam=0.5, colwidths = c(2, 3, 0.8, 1, 1))
dev.off()



## Region/cell type enrichment Adult

topMarkers <- adultMetaGenes_Dendrogram_CellMarkersgsea[which(adultMetaGenes_Dendrogram_CellMarkersgsea$padj < 0.1)][order(pval), pathway]

tiff("output/FGSEA/AdulttopCell_RegionEnrichmentpathways.tiff", width=8.5, height=2, units = 'in', res = 300, compression = "lzw")
plotGseaTable(homemadeDendrogram_CellMarkers[topMarkers], adultrnklist,
              adultMetaGenes_Dendrogram_CellMarkersgsea, gseaParam=0.5, colwidths = c(4, 3, 0.8, 1.3, 1.3))

dev.off()


## Region/cell type enrichment P14

topMarkers <- P14MetaGenes_Dendrogram_CellMarkersgsea[which(P14MetaGenes_Dendrogram_CellMarkersgsea$padj < 0.1)][order(pval), pathway]

tiff("output/FGSEA/P14topCell_RegionEnrichmentpathways.tiff", width=8.5, height=2, units = 'in', res = 300, compression = "lzw")
plotGseaTable(homemadeDendrogram_CellMarkers[topMarkers], P14rnklist,
              P14MetaGenes_Dendrogram_CellMarkersgsea, gseaParam=0.5, colwidths = c(6.1, 3, 0.7, 1, 1))

dev.off()



