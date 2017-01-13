# Wright_2017_ASD_Histamine

Code for differential expression analysis and gene set enrichment anaysis using RNAseq data

This code executes the analyses presented in the "Increased expression of histamine signaling genes in Autism Spectrum Disorder" manuscript.

Two independent data sets are presented in the manuscript. The raw fastq files for one of the datasets, from the Lieber Institute for Brain Development (LIBD) is available on Gene Ominbus....  The replication data was obtained from the Geschwind lab at UCLA. The data and code included here is for the LIBD data only. The Geschwind lab recently published an analysis with their data but they have not yet made their data publicly available. Therfore, the code for the replication analysis will not been shown here untill thier data is made publicly available.

FeatureCounts (v1.4.4) was used to count the total number of reads overlapping each gene using the default settings, with paired-end and reverse-stranded counting specified using the Ensembl Build GRCh37.67 gtf file. These counts were merged from both reads of the paired-end sequencing. This was run using a shell script in linux.

1. Expression Data Normalization:
Code in the RPKM folder normalizes the aligned/assigned read counts from the LIBD into Reads Per Kilobase of Gene per Million mapped reads (RPKM) values.  This normalizes the expression values by library size and coding gene length. These RPKM values were transformed using log2 after applying an offset of 1 to each count to stabilize the variance among lowly expressed genes and to prevent negative normalized counts, resulting in normalized gene expression levels: `log2(RPKM+1)`

2. Differential Expression Analysis
Code in the Analysis Folder first filters the gene expression data to only include genes with a mean normalized value > 0.1.
Then the number of principal components to correct for known and unknown latent confounds using the num.sv "be" method.
Then the influence of diagnosis status on all normalized and filtered gene abundance estimates is evaluated using a linear regression analysis covarying for known RNA-seqconfounder variables and the principal components.

This code depends on R and the Rscript command line utility and the Bioconductor sva and Limma R packages.
To download these packages in R, run:

```R
source("http://bioconductor.org/biocLite.R")
biocLite("sva")
biocLite("limma")
```

