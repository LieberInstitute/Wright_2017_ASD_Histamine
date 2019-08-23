# Wright_2017_ASD_Histamine
[![DOI](https://zenodo.org/badge/78884248.svg)](https://zenodo.org/badge/latestdoi/78884248)

Code for differential expression analysis and gene set enrichment analysis using RNA-seq data

This code executes the analyses presented in the "Altered expression of histamine signaling genes in Autism Spectrum Disorder" manuscript. 

Wright, C et al. “Altered expression of histamine signaling genes in autism spectrum disorder.” Translational psychiatry vol. 7,5 e1126. 9 May. 2017, doi:10.1038/tp.2017.87

Two independent data sets are presented in the manuscript. The gene count files for one of the datasets, from the Lieber Institute for Brain Development (LIBD) is available on Gene Omnibus (Accession: GSE102741). The replication data was obtained from the Geschwind lab at UCLA. The code included here is for the LIBD data only.

FeatureCounts (v1.4.4) was used to count the total number of reads overlapping each gene using the default settings, with paired-end and reverse-stranded counting specified using the Ensembl Build GRCh37.67 gtf file. These counts were merged from both reads of the paired-end sequencing. This was run using a shell script in linux. See `featureCounts_script` for an example.

1. Expression Data Normalization:
Code in the `GeneRpkm2.R` file normalizes the aligned/assigned read counts from the LIBD into Reads Per Kilobase of Gene per Million mapped reads (RPKM) values. This normalizes the expression values by library size and coding gene length. These RPKM values were transformed using log2 after applying an offset of 1 to each count to stabilize the variance among lowly expressed genes and to prevent negative normalized counts, resulting in normalized gene expression levels: `log2(RPKM+1)`

2. Differential Expression Analysis
Code `Analysis.R` first filters the gene expression data to only include genes with a mean normalized value > 0.5.
Then the number of principal components to correct for known and unknown latent confounds is determined using the num.sv "be" method. Next, the influence of diagnosis status on all individual normalized and filtered gene abundance estimates is evaluated using a linear regression analysis covarying for known RNA-seq confounder variables and the principal components.
Finally, gene set analysis is performed to determine if genes of interest are differentially expressed collectively as a group.


This code depends on R and the Bioconductor sva and limma R packages.
To download these packages in R, run:

```R
source("http://bioconductor.org/biocLite.R")
biocLite("sva")
biocLite("limma")
```

