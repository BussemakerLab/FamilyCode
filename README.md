# Predicting the DNA binding specificity of transcription factors variants using family-level biophysically interpretable machine learning
Shaoxun Liu, Pilar Gomez-Alcala, Christ Leemans, William J. Glassford, Lucas A.N. Melo, Richard S. Mann, Harmen J. Bussemaker

BioRxiv link: https://www.biorxiv.org/content/10.1101/2024.01.24.577115v1
Link to the Bussemaker lab: https://bussemakerlab.org

## CISBP DATA DOWNLOAD

rawData/HD/Zscores.txt and rawData/bHLH/Zscores.txt should be downloaded from https://cisbp.ccbr.utoronto.ca/bulk.php (version 3.00). Choose "bHLH" or "Homeodomain" (one at a time) when selecting by family, check boxes "Z-score" and "TF info", and click "Download Family Archive" button to download. The two downloaded ZIP files need to be uncompressed and placed to the corresponding rawData directories.

## FIGURE GENERATION SCRIPT

Running FamilyCodeFigures.Rmd will use the pre-computed files in folder "intermediateData" to generate all figures using precomputed intermediate data. To generate these intermediate data files from the raw data, run the data processing Rmds in the next section.

## DATA PROCESSING SCRIPTS

FamilyCodeOnSELEX_bHLH.Rmd: Data processing using SELEX-seq data as input to perform FamilyCode prediction. Generates Supplemental Data S2-S5.

FamilyCodeOnPBM_bHLH.Rmd: Data processing using PBM data as input to perform FamilyCode prediction.

FamilyCodeOnSELEX+PBM_bHLH.Rmd: Data processing using SELEX-seq and PBM data as input to perform FamilyCode prediction.

FamilyCodeOnSELEX_HD.Rmd: Data processing using SELEX-seq data as input to perform FamilyCode prediction. 

FamilyCodeOnPBM_HD.Rmd: Data processing using PBM data as input to perform FamilyCode prediction.

FamilyCodeOnSELEX+PBM_HD.Rmd: Data processing using SELEX-seq and PBM data as input to perform FamilyCode prediction. Generates Supplemental Data S6-S7.

## SUPPLEMENTAL DATA FILES

Supplemental Data S1: JSON configuration file used for all ProBound analyses

Supplemental Data S2: DNA recognition models for the 52 bHLH factors analyzed in this study.

Supplemental Data S3: Interactive 3D representations of tetrahedrons. Open HTML files in browser to view and manipulate. Related to Figure 1A-C.

Supplemental Data S4: Empirical cumulative distribution of tetrahedral position along each principal component direction for various amino-acid positions. Related to Figure 3B.

Supplemental Data S5: Statistical significance of ANOVA test of PCs at DNA position –2/+2 and –3/+3. Related to Figure 4C.

Supplemental Data S6: DNA recognition models for the 414 HD factors analyzed in this study.

Supplemental Data S7: Interactive 3D representations of tetrahedrons for HD examples. Open HTML files in browser to view and manipulate. Related to Figure 5B.






