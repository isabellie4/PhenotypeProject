# PhenotypeProject

### All the code used to perform the meta-analysis and subsequent analyses on basal differences between the HR and LR rats.

Script containing the updated adult meta-analysis.
  - Adult Meta-Analysis.R

Development meta-analysis code for P7, P14, and P21.
  - MetaAnalysis_Development.R

FGSEA script containing the code used to run fgsea on adult and P14 meta-analysis genes. 
Also contains code for generating custom GMT files from several co-expression studies and running fgsea using these gene set files.
  - fgsea.R

##### Individual dataset re-analysis (7 total – the old RNA-Seq study was not re-analyzed)

Dr. Sarah Clinton New Colony data obtained from Gene Expression Omnibus Website. Contains data for P7, P14, P21, and adult basal rats. 
Only dataset not from our lab. Because the data was obtained online, it is raw data and requires additional processing and annotation. 
An R package from Bioconductor is needed to annotate the data. Data files exist as separate .txt files of expression data for each rat.
  - Alabama_NimbleGen_F34.Rmd

Dr. Sarah Clinton Development. Contains data for P7, P14, and P21 basal rats. 
Data requires annotation so separate annotation file is needed.
  - MBNI_AffymetrixRgU34A_F6.Rmd

Dr. John Stead adult basal rat data. Requires re-annotated and is completely re-analyzed (using original data files) 
because it’s old as sin. Original data files used exist as separate .CEL files for each rat. 
Specific Affymetrix R package from Bioconductor is needed to read in files and annotate the data.
  - MBNI_AffymetrixRae230_F4.Rmd

Dr. Cigdem Aydin vehicle HR/LR adult rat data.
  - MBNI_RNASeq_F43.Rmd

Dr. Peter Blandino adult basal rat data.
  - MBNI_RNASeq_F37.Rmd

Dr. Sarah Clinton P14 Affy. Sarah Clinton’s basal P14 data performed on an Affymetrix platform. 
Data requires annotation so multiple input files are needed.
  - MBNI_AﬀymetrixRae230_ F15.Rmd

Sarah Clinton P14 Illumina. Sarah Clinton’s basal P14 data performed on an Illumina platform.
  - MBNI_IlluminaRatRef12v1_F15.Rmd
