# Predicting the DNA binding specificity of transcription factors variants using family-level biophysically interpretable machine learning
Shaoxun Liu, Pilar Gomez-Alcala, Christ Leemans, William J. Glassford, Lucas A.N. Melo, Richard S. Mann, Harmen J. Bussemaker

BioRxiv link: https://www.biorxiv.org/content/10.1101/2024.01.24.577115v1
Link to the Bussemaker lab: https://bussemakerlab.org/site/

## DOWNLOAD DATA

rawData\HD\Zscores.txt and rawData\bHLH\Zscores.txt should be downloaded from https://cisbp.ccbr.utoronto.ca/bulk.php. Select bHLH/Homeodomain in selection field and select Z-score to download. The Zscore.txt files shall be unziped and placed to the corresponding rawData directories. 

## FIGURE GENERATION

Running FamilyCodeFigures.Rmd alone will generate all figures using precomputed intermediate data. To generate intermediate data locally, run the data processing Rmds. 

## DATA PROCESSING

FamilyCodeOnSELEX_bHLH.Rmd: Data processing using SELEX-seq data as input to perform Family Code prediction. Generates supplimental data 1-4.

FamilyCodeOnPBM_bHLH.Rmd: Data processing using PBM data as input to perform Family Code prediction.

FamilyCodeOnSELEX+PBM_bHLH.Rmd: Data processing using SELEX-seq and PBM data as input to perform Family Code prediction.

FamilyCodeOnSELEX_HD.Rmd: Data processing using SELEX-seq data as input to perform Family Code prediction. 

FamilyCodeOnPBM_HD.Rmd: Data processing using PBM data as input to perform Family Code prediction.

FamilyCodeOnSELEX+PBM_HD.Rmd: Data processing using SELEX-seq and PBM data as input to perform Family Code prediction.

## SUPPLEMENTAL DATA

Supplemental Data S1: DNA recognition models for the 52 bHLH factors analyzed in this study.

Supplemental Data S2: Interactive 3D representations of tetrahedrons. Open HTML files in browser to view and manipulate. Related to Figure 2A-C.

Supplemental Data S3: Empirical cumulative distribution of tetrahedral position along each principal component direction for various amino-acid positions. Related to Figure 4B.

Supplemental Data S4: Statistical significance of ANOVA test of PCs at DNA position –2/+2 and –3/+3. Related to Figure 4C.

Supplemental Data S5: JSON configuration file used for all ProBound analyses




