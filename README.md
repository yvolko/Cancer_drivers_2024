# Identification of new cancer driver genes and associated cellular pathways in high-level amplification regions of genome
The project is done by Iurii Slepov and Yulia Volkova.
Project supervisor: Andrey Yurchenko, PhD.
## Background
Extrachromosomal DNA (ecDNA) is circular DNA outside of chromosomes. It has been often found in different types of cancers. Many known oncogenes reside within ecDNA. Moreover, ecDNA is often present in high copy number in tumor cells and increased expression of genes located in the ecDNA is observed. Maintaining many copies shows that this feature states under positive selection. We use high-amplification regions from different cancer types samples to estimate how many regions have known oncogenes and search for new oncogenes in the regions where no known oncogenes are present.
## Goal
Identify novel cancer drivers and progression genes located in the regions of genome with high-level amplifications.
## Methods and Pipeline
Two datasets were used for the study: ASCAT copy number profiles for WES TCGA Pan-Cancer dataset containing data of 9699 patients [1] and copy number profiles for targeted DNA sequencing using gene panel IMPACT-468 containing data of 17602 patients with primary and metastatic cancer [2]. Genes located in the high-level amplification regions were annotated as presented in the pipeline.
![Pipeline](./pipline_split_in_two_small.png)
Three approaches have been used for the identification of putative cancer driver genes:
1) *System of ranks:* all features characterising genes ranked for each region and genes with total rank above the ninth decile were selected as potential oncogenes. After that, all ranks were normalized within each region and divided into two datasets: a dataset containing regions with known oncogenes, and a dataset containing regions without known oncogenes. The first dataset has been used to train machine learning models such as RandomForestClassifier and Ensemble consisting of RandomForestClassifier, ExtraTreesClassifier, CatBoostClassifier, XGBRFClassifier and GaussianNB. Predictions were made on the second dataset. Genes with a probability of being oncogenes greater than 0.9 were selected as potential cancer driver genes. 
