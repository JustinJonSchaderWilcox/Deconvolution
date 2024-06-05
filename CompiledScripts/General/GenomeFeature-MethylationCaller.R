#Make Methylation Calls across Genomic Features for future analyses
library(data.table)

system("subst x: \"C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation\"")
setwd("x:/")

DeconBrainIslandCpG=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/CpG-Island/mAbel-BrainMix.Blood-AllIslandShared-deduplicated.bismark.deconvoluted-brain-cov.CpG-Island_iCpG.bed", header=TRUE)
DeconBrainIslandNonCpG=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/CpG-Island/mAbel-BrainMix.Blood-AllIslandShared-deduplicated.bismark.deconvoluted-brain-cov.CpG-Island_nCpG.bed", header=TRUE)
DeconBrainIslandMethCpG=read.delim("x:/CpG-Island/mAbel-BrainMix.Blood-AllIslandShared-deduplicated.bismark.deconvoluted-brain-cov.Meth-CpG-Island_iCpG.bed", header=TRUE)
DeconBrainIslandMethNonCpG=read.delim("x:/CpG-Island/mAbel-BrainMix.Blood-AllIslandShared-deduplicated.bismark.deconvoluted-brain-cov.Meth-CpG-Island_nCpG.bed", header=TRUE)
system("subst x: /D")

DeconBrainGeneCpG=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/GeneBody/mAbel-BrainMix.Blood-AllGeneShared-deduplicated.bismark.deconvoluted-brain-cov.CpG-GeneBody.bed", header=TRUE)
DeconBrainGeneNonCpG=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/GeneBody/mAbel-BrainMix.Blood-AllGeneShared-deduplicated.bismark.deconvoluted-brain-cov.nCpG-GeneBody.bed", header=TRUE)
DeconBrainGeneMethCpG=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/GeneBody/mAbel-BrainMix.Blood-AllGeneShared-deduplicated.bismark.deconvoluted-brain-cov.Meth-CpG-GeneBody.bed", header=TRUE)
DeconBrainGeneMethNonCpG=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/GeneBody/mAbel-BrainMix.Blood-AllGeneShared-deduplicated.bismark.deconvoluted-brain-cov.Meth-nCpG-GeneBody.bed", header=TRUE)


DeconBrainIslandCpG$DECONvBRAINdMP=1
for (i in 1:dim(DeconBrainIslandCpG)[1]){
DeconBrainIslandCpG$DECONvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(DeconBrainIslandCpG$DeconMethReads[i], DeconBrainIslandCpG$DeconUnmethReads[i]), c(DeconBrainIslandCpG$BrainMethReads[i], DeconBrainIslandCpG$BrainUnmethReads[i]))))$p.value
}
DeconBrainIslandCpG$tpDECONvBRAINdMP=1
for (i in 1:dim(DeconBrainIslandCpG)[1]){
DeconBrainIslandCpG$tpDECONvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(DeconBrainIslandCpG$tpMethReads[i], DeconBrainIslandCpG$tpUnmethReads[i]), c(DeconBrainIslandCpG$BrainMethReads[i], DeconBrainIslandCpG$BrainUnmethReads[i]))))$p.value
}
DeconBrainIslandCpG$MIXvBRAINdMP=1
for (i in 1:dim(DeconBrainIslandCpG)[1]){
DeconBrainIslandCpG$MIXvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(DeconBrainIslandCpG$MixMethReads[i], DeconBrainIslandCpG$MixUnmethReads[i]), c(DeconBrainIslandCpG$BrainMethReads[i], DeconBrainIslandCpG$BrainUnmethReads[i]))))$p.value
}

DeconBrainIslandCpG$TrueDiffMethP=1
for (i in 1:dim(DeconBrainIslandCpG)[1]){
DeconBrainIslandCpG$TrueDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainIslandCpG$BrainMethReads[i], DeconBrainIslandCpG$BrainUnmethReads[i]), c(DeconBrainIslandCpG$BloodMethReads[i], DeconBrainIslandCpG$BloodUnmethReads[i]))))$p.value
}
DeconBrainIslandCpG$PredictDiffMethP=1
for (i in 1:dim(DeconBrainIslandCpG)[1]){
DeconBrainIslandCpG$PredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainIslandCpG$DeconMethReads[i], DeconBrainIslandCpG$DeconUnmethReads[i]), c(DeconBrainIslandCpG$BloodMethReads[i], DeconBrainIslandCpG$BloodUnmethReads[i]))))$p.value
}
DeconBrainIslandCpG$tpPredictDiffMethP=1
for (i in 1:dim(DeconBrainIslandCpG)[1]){
DeconBrainIslandCpG$tpPredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainIslandCpG$tpMethReads[i], DeconBrainIslandCpG$tpUnmethReads[i]), c(DeconBrainIslandCpG$BloodMethReads[i], DeconBrainIslandCpG$BloodUnmethReads[i]))))$p.value
}
DeconBrainIslandCpG$MixDiffMethP=1
for (i in 1:dim(DeconBrainIslandCpG)[1]){
DeconBrainIslandCpG$MixDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainIslandCpG$MixMethReads[i], DeconBrainIslandCpG$MixUnmethReads[i]), c(DeconBrainIslandCpG$BloodMethReads[i], DeconBrainIslandCpG$BloodUnmethReads[i]))))$p.value
}

write.table(DeconBrainIslandCpG, file="C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/DeconBrainIslandCpG.TRUE-P.txt", sep="\t")

DeconBrainIslandNonCpG$TrueDiffMethP=1
for (i in 1:dim(DeconBrainIslandNonCpG)[1]){
DeconBrainIslandNonCpG$TrueDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainIslandNonCpG$BrainMethReads[i], DeconBrainIslandNonCpG$BrainUnmethReads[i]), c(DeconBrainIslandNonCpG$BloodMethReads[i], DeconBrainIslandNonCpG$BloodUnmethReads[i]))))$p.value
}
DeconBrainIslandNonCpG$PredictDiffMethP=1
for (i in 1:dim(DeconBrainIslandNonCpG)[1]){
DeconBrainIslandNonCpG$PredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainIslandNonCpG$DeconMethReads[i], DeconBrainIslandNonCpG$DeconUnmethReads[i]), c(DeconBrainIslandNonCpG$BloodMethReads[i], DeconBrainIslandNonCpG$BloodUnmethReads[i]))))$p.value
}
DeconBrainIslandNonCpG$tpPredictDiffMethP=1
for (i in 1:dim(DeconBrainIslandNonCpG)[1]){
DeconBrainIslandNonCpG$tpPredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainIslandNonCpG$tpMethReads[i], DeconBrainIslandNonCpG$tpUnmethReads[i]), c(DeconBrainIslandNonCpG$BloodMethReads[i], DeconBrainIslandNonCpG$BloodUnmethReads[i]))))$p.value
}
DeconBrainIslandNonCpG$MixDiffMethP=1
for (i in 1:dim(DeconBrainIslandNonCpG)[1]){
DeconBrainIslandNonCpG$MixDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainIslandNonCpG$MixMethReads[i], DeconBrainIslandNonCpG$MixUnmethReads[i]), c(DeconBrainIslandNonCpG$BloodMethReads[i], DeconBrainIslandNonCpG$BloodUnmethReads[i]))))$p.value
}

write.table(DeconBrainIslandNonCpG, file="C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/DeconBrainIslandNonCpG.TRUE-P.txt", sep="\t")


DeconBrainGeneCpG$DECONvBRAINdMP=1
for (i in 1:dim(DeconBrainGeneCpG)[1]){
DeconBrainGeneCpG$DECONvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneCpG$DeconMethReads[i], DeconBrainGeneCpG$DeconUnmethReads[i]), c(DeconBrainGeneCpG$BrainMethReads[i], DeconBrainGeneCpG$BrainUnmethReads[i]))))$p.value
}
DeconBrainGeneNonCpG$DECONvBRAINdMP=1
for (i in 1:dim(DeconBrainGeneNonCpG)[1]){
DeconBrainGeneNonCpG$DECONvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneNonCpG$DeconMethReads[i], DeconBrainGeneNonCpG$DeconUnmethReads[i]), c(DeconBrainGeneNonCpG$BrainMethReads[i], DeconBrainGeneNonCpG$BrainUnmethReads[i]))))$p.value
}

DeconBrainGeneCpG$tpDECONvBRAINdMP=1
for (i in 1:dim(DeconBrainGeneCpG)[1]){
DeconBrainGeneCpG$tpDECONvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneCpG$tpMethReads[i], DeconBrainGeneCpG$tpUnmethReads[i]), c(DeconBrainGeneCpG$BrainMethReads[i], DeconBrainGeneCpG$BrainUnmethReads[i]))))$p.value
}
DeconBrainGeneNonCpG$tpDECONvBRAINdMP=1
for (i in 1:dim(DeconBrainGeneNonCpG)[1]){
DeconBrainGeneNonCpG$tpDECONvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneNonCpG$tpMethReads[i], DeconBrainGeneNonCpG$tpUnmethReads[i]), c(DeconBrainGeneNonCpG$BrainMethReads[i], DeconBrainGeneNonCpG$BrainUnmethReads[i]))))$p.value
}

DeconBrainGeneCpG$MIXvBRAINdMP=1
for (i in 1:dim(DeconBrainGeneCpG)[1]){
DeconBrainGeneCpG$MIXvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneCpG$MixMethReads[i], DeconBrainGeneCpG$MixUnmethReads[i]), c(DeconBrainGeneCpG$BrainMethReads[i], DeconBrainGeneCpG$BrainUnmethReads[i]))))$p.value
}
DeconBrainGeneNonCpG$MIXvBRAINdMP=1
for (i in 1:dim(DeconBrainGeneNonCpG)[1]){
DeconBrainGeneNonCpG$MIXvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneNonCpG$MixMethReads[i], DeconBrainGeneNonCpG$MixUnmethReads[i]), c(DeconBrainGeneNonCpG$BrainMethReads[i], DeconBrainGeneNonCpG$BrainUnmethReads[i]))))$p.value
}

DeconBrainGeneCpG$TrueDiffMethP=1
for (i in 1:dim(DeconBrainGeneCpG)[1]){
DeconBrainGeneCpG$TrueDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneCpG$BrainMethReads[i], DeconBrainGeneCpG$BrainUnmethReads[i]), c(DeconBrainGeneCpG$BloodMethReads[i], DeconBrainGeneCpG$BloodUnmethReads[i]))))$p.value
}
DeconBrainGeneNonCpG$TrueDiffMethP=1
for (i in 1:dim(DeconBrainGeneNonCpG)[1]){
DeconBrainGeneNonCpG$TrueDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneNonCpG$BrainMethReads[i], DeconBrainGeneNonCpG$BrainUnmethReads[i]), c(DeconBrainGeneNonCpG$BloodMethReads[i], DeconBrainGeneNonCpG$BloodUnmethReads[i]))))$p.value
}

DeconBrainGeneCpG$PredictDiffMethP=1
for (i in 1:dim(DeconBrainGeneCpG)[1]){
DeconBrainGeneCpG$PredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneCpG$DeconMethReads[i], DeconBrainGeneCpG$DeconUnmethReads[i]), c(DeconBrainGeneCpG$BloodMethReads[i], DeconBrainGeneCpG$BloodUnmethReads[i]))))$p.value
}
DeconBrainGeneNonCpG$PredictDiffMethP=1
for (i in 1:dim(DeconBrainGeneNonCpG)[1]){
DeconBrainGeneNonCpG$PredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneNonCpG$DeconMethReads[i], DeconBrainGeneNonCpG$DeconUnmethReads[i]), c(DeconBrainGeneNonCpG$BloodMethReads[i], DeconBrainGeneNonCpG$BloodUnmethReads[i]))))$p.value
}

DeconBrainGeneCpG$tpPredictDiffMethP=1
for (i in 1:dim(DeconBrainGeneCpG)[1]){
DeconBrainGeneCpG$tpPredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneCpG$tpMethReads[i], DeconBrainGeneCpG$tpUnmethReads[i]), c(DeconBrainGeneCpG$BloodMethReads[i], DeconBrainGeneCpG$BloodUnmethReads[i]))))$p.value
}
DeconBrainGeneNonCpG$tpPredictDiffMethP=1
for (i in 1:dim(DeconBrainGeneNonCpG)[1]){
DeconBrainGeneNonCpG$tpPredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneNonCpG$tpMethReads[i], DeconBrainGeneNonCpG$tpUnmethReads[i]), c(DeconBrainGeneNonCpG$BloodMethReads[i], DeconBrainGeneNonCpG$BloodUnmethReads[i]))))$p.value
}

DeconBrainGeneCpG$MixDiffMethP=1
for (i in 1:dim(DeconBrainGeneCpG)[1]){
DeconBrainGeneCpG$MixDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneCpG$MixMethReads[i], DeconBrainGeneCpG$MixUnmethReads[i]), c(DeconBrainGeneCpG$BloodMethReads[i], DeconBrainGeneCpG$BloodUnmethReads[i]))))$p.value
}
DeconBrainGeneNonCpG$MixDiffMethP=1
for (i in 1:dim(DeconBrainGeneNonCpG)[1]){
DeconBrainGeneNonCpG$MixDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneNonCpG$MixMethReads[i], DeconBrainGeneNonCpG$MixUnmethReads[i]), c(DeconBrainGeneNonCpG$BloodMethReads[i], DeconBrainGeneNonCpG$BloodUnmethReads[i]))))$p.value
}

write.table(DeconBrainGeneCpG, file="C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/DeconBrainGeneCpG.TRUE-P.txt", sep="\t")
write.table(DeconBrainGeneNonCpG, file="C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/DeconBrainGeneNonCpG.TRUE-P.txt", sep="\t")

DeconBrainIslandMethCpG$DECONvBRAINdMP=1
for (i in 1:dim(DeconBrainIslandMethCpG)[1]){
DeconBrainIslandMethCpG$DECONvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(DeconBrainIslandMethCpG$DeconMethReads[i], DeconBrainIslandMethCpG$DeconUnmethReads[i]), c(DeconBrainIslandMethCpG$BrainMethReads[i], DeconBrainIslandMethCpG$BrainUnmethReads[i]))))$p.value
}
DeconBrainIslandMethCpG$tpDECONvBRAINdMP=1
for (i in 1:dim(DeconBrainIslandMethCpG)[1]){
DeconBrainIslandMethCpG$tpDECONvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(DeconBrainIslandMethCpG$tpMethReads[i], DeconBrainIslandMethCpG$tpUnmethReads[i]), c(DeconBrainIslandMethCpG$BrainMethReads[i], DeconBrainIslandMethCpG$BrainUnmethReads[i]))))$p.value
}
DeconBrainIslandMethCpG$MIXvBRAINdMP=1
for (i in 1:dim(DeconBrainIslandMethCpG)[1]){
DeconBrainIslandMethCpG$MIXvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(DeconBrainIslandMethCpG$MixMethReads[i], DeconBrainIslandMethCpG$MixUnmethReads[i]), c(DeconBrainIslandMethCpG$BrainMethReads[i], DeconBrainIslandMethCpG$BrainUnmethReads[i]))))$p.value
}

DeconBrainIslandMethCpG$TrueDiffMethP=1
for (i in 1:dim(DeconBrainIslandMethCpG)[1]){
DeconBrainIslandMethCpG$TrueDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainIslandMethCpG$BrainMethReads[i], DeconBrainIslandMethCpG$BrainUnmethReads[i]), c(DeconBrainIslandMethCpG$BloodMethReads[i], DeconBrainIslandMethCpG$BloodUnmethReads[i]))))$p.value
}
DeconBrainIslandMethCpG$PredictDiffMethP=1
for (i in 1:dim(DeconBrainIslandMethCpG)[1]){
DeconBrainIslandMethCpG$PredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainIslandMethCpG$DeconMethReads[i], DeconBrainIslandMethCpG$DeconUnmethReads[i]), c(DeconBrainIslandMethCpG$BloodMethReads[i], DeconBrainIslandMethCpG$BloodUnmethReads[i]))))$p.value
}
DeconBrainIslandMethCpG$tpPredictDiffMethP=1
for (i in 1:dim(DeconBrainIslandMethCpG)[1]){
DeconBrainIslandMethCpG$tpPredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainIslandMethCpG$tpMethReads[i], DeconBrainIslandMethCpG$tpUnmethReads[i]), c(DeconBrainIslandMethCpG$BloodMethReads[i], DeconBrainIslandMethCpG$BloodUnmethReads[i]))))$p.value
}
DeconBrainIslandMethCpG$MixDiffMethP=1
for (i in 1:dim(DeconBrainIslandMethCpG)[1]){
DeconBrainIslandMethCpG$MixDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainIslandMethCpG$MixMethReads[i], DeconBrainIslandMethCpG$MixUnmethReads[i]), c(DeconBrainIslandMethCpG$BloodMethReads[i], DeconBrainIslandMethCpG$BloodUnmethReads[i]))))$p.value
}

write.table(DeconBrainIslandMethCpG, file="C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/DeconBrainIslandMethCpG.TRUE-P.txt", sep="\t")

DeconBrainIslandMethNonCpG$TrueDiffMethP=1
for (i in 1:dim(DeconBrainIslandMethNonCpG)[1]){
DeconBrainIslandMethNonCpG$TrueDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainIslandMethNonCpG$BrainMethReads[i], DeconBrainIslandMethNonCpG$BrainUnmethReads[i]), c(DeconBrainIslandMethNonCpG$BloodMethReads[i], DeconBrainIslandMethNonCpG$BloodUnmethReads[i]))))$p.value
}
DeconBrainIslandMethNonCpG$PredictDiffMethP=1
for (i in 1:dim(DeconBrainIslandMethNonCpG)[1]){
DeconBrainIslandMethNonCpG$PredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainIslandMethNonCpG$DeconMethReads[i], DeconBrainIslandMethNonCpG$DeconUnmethReads[i]), c(DeconBrainIslandMethNonCpG$BloodMethReads[i], DeconBrainIslandMethNonCpG$BloodUnmethReads[i]))))$p.value
}
DeconBrainIslandMethNonCpG$tpPredictDiffMethP=1
for (i in 1:dim(DeconBrainIslandMethNonCpG)[1]){
DeconBrainIslandMethNonCpG$tpPredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainIslandMethNonCpG$tpMethReads[i], DeconBrainIslandMethNonCpG$tpUnmethReads[i]), c(DeconBrainIslandMethNonCpG$BloodMethReads[i], DeconBrainIslandMethNonCpG$BloodUnmethReads[i]))))$p.value
}
DeconBrainIslandMethNonCpG$MixDiffMethP=1
for (i in 1:dim(DeconBrainIslandMethNonCpG)[1]){
DeconBrainIslandMethNonCpG$MixDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainIslandMethNonCpG$MixMethReads[i], DeconBrainIslandMethNonCpG$MixUnmethReads[i]), c(DeconBrainIslandMethNonCpG$BloodMethReads[i], DeconBrainIslandMethNonCpG$BloodUnmethReads[i]))))$p.value
}

write.table(DeconBrainIslandMethNonCpG, file="C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/DeconBrainIslandMethNonCpG.TRUE-P.txt", sep="\t")


DeconBrainGeneMethCpG$DECONvBRAINdMP=1
for (i in 1:dim(DeconBrainGeneMethCpG)[1]){
DeconBrainGeneMethCpG$DECONvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneMethCpG$DeconMethReads[i], DeconBrainGeneMethCpG$DeconUnmethReads[i]), c(DeconBrainGeneMethCpG$BrainMethReads[i], DeconBrainGeneMethCpG$BrainUnmethReads[i]))))$p.value
}
DeconBrainGeneMethNonCpG$DECONvBRAINdMP=1
for (i in 1:dim(DeconBrainGeneMethNonCpG)[1]){
DeconBrainGeneMethNonCpG$DECONvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneMethNonCpG$DeconMethReads[i], DeconBrainGeneMethNonCpG$DeconUnmethReads[i]), c(DeconBrainGeneMethNonCpG$BrainMethReads[i], DeconBrainGeneMethNonCpG$BrainUnmethReads[i]))))$p.value
}

DeconBrainGeneMethCpG$tpDECONvBRAINdMP=1
for (i in 1:dim(DeconBrainGeneMethCpG)[1]){
DeconBrainGeneMethCpG$tpDECONvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneMethCpG$tpMethReads[i], DeconBrainGeneMethCpG$tpUnmethReads[i]), c(DeconBrainGeneMethCpG$BrainMethReads[i], DeconBrainGeneMethCpG$BrainUnmethReads[i]))))$p.value
}
DeconBrainGeneMethNonCpG$tpDECONvBRAINdMP=1
for (i in 1:dim(DeconBrainGeneMethNonCpG)[1]){
DeconBrainGeneMethNonCpG$tpDECONvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneMethNonCpG$tpMethReads[i], DeconBrainGeneMethNonCpG$tpUnmethReads[i]), c(DeconBrainGeneMethNonCpG$BrainMethReads[i], DeconBrainGeneMethNonCpG$BrainUnmethReads[i]))))$p.value
}

DeconBrainGeneMethCpG$MIXvBRAINdMP=1
for (i in 1:dim(DeconBrainGeneMethCpG)[1]){
DeconBrainGeneMethCpG$MIXvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneMethCpG$MixMethReads[i], DeconBrainGeneMethCpG$MixUnmethReads[i]), c(DeconBrainGeneMethCpG$BrainMethReads[i], DeconBrainGeneMethCpG$BrainUnmethReads[i]))))$p.value
}
DeconBrainGeneMethNonCpG$MIXvBRAINdMP=1
for (i in 1:dim(DeconBrainGeneMethNonCpG)[1]){
DeconBrainGeneMethNonCpG$MIXvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneMethNonCpG$MixMethReads[i], DeconBrainGeneMethNonCpG$MixUnmethReads[i]), c(DeconBrainGeneMethNonCpG$BrainMethReads[i], DeconBrainGeneMethNonCpG$BrainUnmethReads[i]))))$p.value
}

DeconBrainGeneMethCpG$TrueDiffMethP=1
for (i in 1:dim(DeconBrainGeneMethCpG)[1]){
DeconBrainGeneMethCpG$TrueDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneMethCpG$BrainMethReads[i], DeconBrainGeneMethCpG$BrainUnmethReads[i]), c(DeconBrainGeneMethCpG$BloodMethReads[i], DeconBrainGeneMethCpG$BloodUnmethReads[i]))))$p.value
}
DeconBrainGeneMethNonCpG$TrueDiffMethP=1
for (i in 1:dim(DeconBrainGeneMethNonCpG)[1]){
DeconBrainGeneMethNonCpG$TrueDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneMethNonCpG$BrainMethReads[i], DeconBrainGeneMethNonCpG$BrainUnmethReads[i]), c(DeconBrainGeneMethNonCpG$BloodMethReads[i], DeconBrainGeneMethNonCpG$BloodUnmethReads[i]))))$p.value
}

DeconBrainGeneMethCpG$PredictDiffMethP=1
for (i in 1:dim(DeconBrainGeneMethCpG)[1]){
DeconBrainGeneMethCpG$PredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneMethCpG$DeconMethReads[i], DeconBrainGeneMethCpG$DeconUnmethReads[i]), c(DeconBrainGeneMethCpG$BloodMethReads[i], DeconBrainGeneMethCpG$BloodUnmethReads[i]))))$p.value
}
DeconBrainGeneMethNonCpG$PredictDiffMethP=1
for (i in 1:dim(DeconBrainGeneMethNonCpG)[1]){
DeconBrainGeneMethNonCpG$PredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneMethNonCpG$DeconMethReads[i], DeconBrainGeneMethNonCpG$DeconUnmethReads[i]), c(DeconBrainGeneMethNonCpG$BloodMethReads[i], DeconBrainGeneMethNonCpG$BloodUnmethReads[i]))))$p.value
}

DeconBrainGeneMethCpG$tpPredictDiffMethP=1
for (i in 1:dim(DeconBrainGeneMethCpG)[1]){
DeconBrainGeneMethCpG$tpPredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneMethCpG$tpMethReads[i], DeconBrainGeneMethCpG$tpUnmethReads[i]), c(DeconBrainGeneMethCpG$BloodMethReads[i], DeconBrainGeneMethCpG$BloodUnmethReads[i]))))$p.value
}
DeconBrainGeneMethNonCpG$tpPredictDiffMethP=1
for (i in 1:dim(DeconBrainGeneMethNonCpG)[1]){
DeconBrainGeneMethNonCpG$tpPredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneMethNonCpG$tpMethReads[i], DeconBrainGeneMethNonCpG$tpUnmethReads[i]), c(DeconBrainGeneMethNonCpG$BloodMethReads[i], DeconBrainGeneMethNonCpG$BloodUnmethReads[i]))))$p.value
}

DeconBrainGeneMethCpG$MixDiffMethP=1
for (i in 1:dim(DeconBrainGeneMethCpG)[1]){
DeconBrainGeneMethCpG$MixDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneMethCpG$MixMethReads[i], DeconBrainGeneMethCpG$MixUnmethReads[i]), c(DeconBrainGeneMethCpG$BloodMethReads[i], DeconBrainGeneMethCpG$BloodUnmethReads[i]))))$p.value
}
DeconBrainGeneMethNonCpG$MixDiffMethP=1
for (i in 1:dim(DeconBrainGeneMethNonCpG)[1]){
DeconBrainGeneMethNonCpG$MixDiffMethP[i]=fisher.test(as.matrix(cbind(c(DeconBrainGeneMethNonCpG$MixMethReads[i], DeconBrainGeneMethNonCpG$MixUnmethReads[i]), c(DeconBrainGeneMethNonCpG$BloodMethReads[i], DeconBrainGeneMethNonCpG$BloodUnmethReads[i]))))$p.value
}

write.table(DeconBrainGeneMethCpG, file="C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/DeconBrainGeneMethCpG.TRUE-P.txt", sep="\t")
write.table(DeconBrainGeneMethNonCpG, file="C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/DeconBrainGeneMethNonCpG.TRUE-P.txt", sep="\t")
