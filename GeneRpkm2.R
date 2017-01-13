load("RiboHistGeneCounts.rda")
Pheno<-read.table("RiboPheno.txt", header = T)
Phenodup<-read.table("RiboPhenonobad.txt", header = T)
Phenodup$NumMappedcurrentFlowcell<-as.numeric(gsub(",", "", Phenodup$NumMappedcurrentFlowcell)) 
Phenodup$rRNA<-as.numeric(gsub(",", "", Phenodup$rRNA)) 
Phenodup$mRNA<-as.numeric(gsub(",", "", Phenodup$mRNA)) 
Phenodup$rperc<- round(((100)*(Phenodup$rRNA/Phenodup$NumMappedcurrentFlowcell)), digits = 2)
Phenodup$chrmperc<- round(((100)*(Phenodup$mRNA/Phenodup$NumMappedcurrentFlowcell)), digits = 2)
save(Phenodup, file = "RiboPhenodupN=60")
Phenofull<-read.table("Ribo52fullPheno.txt", header = T)
Phenooverlap <- read.table("Ribooverlap.txt", header = T)#####
PhenoFC <- read.table("RibooverlapFC.txt", header = T)
Geneanno = GeneCounts[,1:6]
gCounts = GeneCounts[, 7:length(colnames(GeneCounts))]
GoodqualityIndex = which(Pheno$Dx !="BADquality")
GeneCountsGood = gCounts[GoodqualityIndex]
GeneCounts<-data.frame(GeneCountsGood)####N=60
dim(GeneCounts)###N=60
head(colnames(GeneCounts))####12325, 12325, 12327...
Pheno<-Phenodup
test = (nrow(Pheno) == ncol(GeneCounts))
if (test == TRUE) {print("good")
} else {
  print("bad - Check on things!!!")
}
# Test<-as.matrix(Pheno)
# Pheno<-as.data.frame(Test)

#IF you aren't doing multiple specific flowcells per subject than use TotalNumMapped...otherwise use Pheno$NumMappedcurrentFlowcell:


Pheno$NumMappedcurrentFlowcell<-as.numeric(gsub(",", "", Pheno$NumMappedcurrentFlowcell)) 
genebgE= matrix(rep(Pheno$NumMappedcurrentFlowcell), nc = nrow(Pheno), nr = nrow(GeneCounts), byrow = T)
genewidE = matrix(rep(Geneanno$Length), nr = nrow(GeneCounts), 
                  nc = nrow(Pheno),  byrow=F)
geneRpkm<- GeneCounts/(genewidE/1000)/(genebgE/1e6)### we only want to calculate Rpkm on the counts! Not the annotation of location!!!
geneRpkm_withanno <- data.frame(Geneanno, geneRpkm)
geneRpkm2 <- data.frame(Geneanno, log2(geneRpkm + 1))

###Now to check that it worked:
###Correct dimensions????
test = ((length(geneRpkm_withanno) == length(geneRpkm2)) &  (nrow(geneRpkm_withanno) == nrow(geneRpkm2)))
if (test == TRUE) {print("good")
} else {
  print("bad - Check on things!!!")
}
###Nonzero values????
geneRpkm2Index = which(as.numeric(geneRpkm2$X.media.DATA.GroupStudy.DLPFC_RiboZeroGold.Autism_RiboZero.Sample_R12325_C4UYCACXX_out.accepted_hits.bam)>0)##### to check if we have nonzero counts
test = (length(geneRpkm2Index) > 0)
if (test == TRUE) {print("good")
} else {
  print("bad - Check on things!!!")
}

save(Phenodup, geneRpkm2, geneRpkm, geneRpkm_withanno, file = "RiboHistgeneRpkm2.rda")