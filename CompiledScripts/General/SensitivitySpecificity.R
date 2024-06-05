#Methylated Sites Sensitivity and Specificity Analysis
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(DescTools)
library(GenBinomApps)

BloodDeconP=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/BasicStatistics/blood-deconvolution-validation-table.txt", header=TRUE, sep="\t")

MethCallDeconBrainMethWGoCpG=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/MergedDeconBrainMethWGoCpG.TRUE-P.txt", header=TRUE)
MethCallDeconBrainMethWGnCpG=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/MergedDeconBrainMethWGnCpG.TRUE-P.txt", header=TRUE)
MethCallDeconBrainMethWGoCpG=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/MergedDeconBrainMethWGoCpG.TRUE-P.txt", header=TRUE)
MethCallDeconBrainMethWGnCpG=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/MergedDeconBrainMethWGnCpG.TRUE-P.txt", header=TRUE)


MethCallDeconBrainIslandMethCpG=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/DeconBrainIslandMethCpG.TRUE-P.txt", header=TRUE)
MethCallDeconBrainIslandMethNonCpG=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/DeconBrainIslandNonCpG.TRUE-P.txt", header=TRUE)

MethCallDeconBrainGeneMethCpG=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/DeconBrainGeneMethCpG.TRUE-P.txt", header=TRUE)
MethCallDeconBrainGeneMethNonCpG=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/DeconBrainGeneNonCpG.TRUE-P.txt", header=TRUE)

MethCallDeconBrainIslandMethCpG$DeconPROP=MethCallDeconBrainIslandMethCpG$DeconMethReads/(MethCallDeconBrainIslandMethCpG$DeconMethReads+MethCallDeconBrainIslandMethCpG$DeconUnmethReads)
MethCallDeconBrainIslandMethCpG$BrainPROP=MethCallDeconBrainIslandMethCpG$BrainMethReads/(MethCallDeconBrainIslandMethCpG$BrainMethReads+MethCallDeconBrainIslandMethCpG$BrainUnmethReads)
MethCallDeconBrainIslandMethCpG$BloodPROP=MethCallDeconBrainIslandMethCpG$BloodMethReads/(MethCallDeconBrainIslandMethCpG$BloodMethReads+MethCallDeconBrainIslandMethCpG$BloodUnmethReads)
MethCallDeconBrainIslandMethNonCpG$BrainPROP=MethCallDeconBrainIslandMethNonCpG$BrainMethReads/(MethCallDeconBrainIslandMethNonCpG$BrainMethReads+MethCallDeconBrainIslandMethNonCpG$BrainUnmethReads)
MethCallDeconBrainIslandMethNonCpG$BloodPROP=MethCallDeconBrainIslandMethNonCpG$BloodMethReads/(MethCallDeconBrainIslandMethNonCpG$BloodMethReads+MethCallDeconBrainIslandMethNonCpG$BloodUnmethReads)
MethCallDeconBrainGeneMethCpG$DeconPROP=MethCallDeconBrainGeneMethCpG$DeconMethReads/(MethCallDeconBrainGeneMethCpG$DeconMethReads+MethCallDeconBrainGeneMethCpG$DeconUnmethReads)
MethCallDeconBrainGeneMethCpG$BrainPROP=MethCallDeconBrainGeneMethCpG$BrainMethReads/(MethCallDeconBrainGeneMethCpG$BrainMethReads+MethCallDeconBrainGeneMethCpG$BrainUnmethReads)
MethCallDeconBrainGeneMethCpG$BloodPROP=MethCallDeconBrainGeneMethCpG$BloodMethReads/(MethCallDeconBrainGeneMethCpG$BloodMethReads+MethCallDeconBrainGeneMethCpG$BloodUnmethReads)
MethCallDeconBrainGeneMethNonCpG$DeconPROP=MethCallDeconBrainGeneMethNonCpG$DeconMethReads/(MethCallDeconBrainGeneMethNonCpG$DeconMethReads+MethCallDeconBrainGeneMethNonCpG$DeconUnmethReads)
MethCallDeconBrainGeneMethNonCpG$BrainPROP=MethCallDeconBrainGeneMethNonCpG$BrainMethReads/(MethCallDeconBrainGeneMethNonCpG$BrainMethReads+MethCallDeconBrainGeneMethNonCpG$BrainUnmethReads)
MethCallDeconBrainGeneMethNonCpG$BloodPROP=MethCallDeconBrainGeneMethNonCpG$BloodMethReads/(MethCallDeconBrainGeneMethNonCpG$BloodMethReads+MethCallDeconBrainGeneMethNonCpG$BloodUnmethReads)

MethCallDeconBrainMethWGoCpG.diffmeth=MethCallDeconBrainMethWGoCpG[0,]
MethCallDeconBrainMethWGnCpG.diffmeth=MethCallDeconBrainMethWGnCpG[0,]
MethCallDeconBrainIslandMethCpG.diffmeth=MethCallDeconBrainIslandMethCpG[0,]
MethCallDeconBrainIslandMethNonCpG.diffmeth=MethCallDeconBrainIslandMethNonCpG[0,]
MethCallDeconBrainGeneMethCpG.diffmeth=MethCallDeconBrainGeneMethCpG[0,]
MethCallDeconBrainGeneMethNonCpG.diffmeth=MethCallDeconBrainGeneMethNonCpG[0,]

MethCallDeconBrainMethWGoCpG.predict=MethCallDeconBrainMethWGoCpG[0,]
MethCallDeconBrainMethWGnCpG.predict=MethCallDeconBrainMethWGnCpG[0,]
MethCallDeconBrainIslandMethCpG.predict=MethCallDeconBrainIslandMethCpG[0,]
MethCallDeconBrainIslandMethNonCpG.predict=MethCallDeconBrainIslandMethNonCpG[0,]
MethCallDeconBrainGeneMethCpG.predict=MethCallDeconBrainGeneMethCpG[0,]
MethCallDeconBrainGeneMethNonCpG.predict=MethCallDeconBrainGeneMethNonCpG[0,]

MethCallDeconBrainMethWGoCpG.tp_predict=MethCallDeconBrainMethWGoCpG[0,]
MethCallDeconBrainMethWGnCpG.tp_predict=MethCallDeconBrainMethWGnCpG[0,]
MethCallDeconBrainIslandMethCpG.tp_predict=MethCallDeconBrainIslandMethCpG[0,]
MethCallDeconBrainIslandMethNonCpG.tp_predict=MethCallDeconBrainIslandMethNonCpG[0,]
MethCallDeconBrainGeneMethCpG.tp_predict=MethCallDeconBrainGeneMethCpG[0,]
MethCallDeconBrainGeneMethNonCpG.tp_predict=MethCallDeconBrainGeneMethNonCpG[0,]

for (p in c(0,10,20,30,40,50,60,70,80,90,100)){
MethCallDeconBrainMethWGoCpG.diffmeth.temp=MethCallDeconBrainMethWGoCpG[which(MethCallDeconBrainMethWGoCpG$P==p),][which(p.adjust(MethCallDeconBrainMethWGoCpG[which(MethCallDeconBrainMethWGoCpG$P==p),]$TrueDiffMethP, method="none")<0.01),]
MethCallDeconBrainMethWGnCpG.diffmeth.temp=MethCallDeconBrainMethWGnCpG[which(MethCallDeconBrainMethWGnCpG$P==p),][which(p.adjust(MethCallDeconBrainMethWGnCpG[which(MethCallDeconBrainMethWGnCpG$P==p),]$TrueDiffMethP, method="none")<0.01),]
MethCallDeconBrainIslandMethCpG.diffmeth.temp=MethCallDeconBrainIslandMethCpG[which(MethCallDeconBrainIslandMethCpG$P==p),][which(p.adjust(MethCallDeconBrainIslandMethCpG[which(MethCallDeconBrainIslandMethCpG$P==p),]$TrueDiffMethP, method="none")<0.01),]
MethCallDeconBrainIslandMethNonCpG.diffmeth.temp=MethCallDeconBrainIslandMethNonCpG[which(MethCallDeconBrainIslandMethNonCpG$P==p),][which(p.adjust(MethCallDeconBrainIslandMethNonCpG[which(MethCallDeconBrainIslandMethNonCpG$P==p),]$TrueDiffMethP, method="none")<0.01),]
MethCallDeconBrainGeneMethCpG.diffmeth.temp=MethCallDeconBrainGeneMethCpG[which(MethCallDeconBrainGeneMethCpG$P==p),][which(p.adjust(MethCallDeconBrainGeneMethCpG[which(MethCallDeconBrainGeneMethCpG$P==p),]$TrueDiffMethP, method="none")<0.01),]
MethCallDeconBrainGeneMethNonCpG.diffmeth.temp=MethCallDeconBrainGeneMethNonCpG[which(MethCallDeconBrainGeneMethNonCpG$P==p),][which(p.adjust(MethCallDeconBrainGeneMethNonCpG[which(MethCallDeconBrainGeneMethNonCpG$P==p),]$TrueDiffMethP, method="none")<0.01),]
MethCallDeconBrainMethWGoCpG.diffmeth=rbind(MethCallDeconBrainMethWGoCpG.diffmeth, MethCallDeconBrainMethWGoCpG.diffmeth.temp)
MethCallDeconBrainMethWGnCpG.diffmeth=rbind(MethCallDeconBrainMethWGnCpG.diffmeth, MethCallDeconBrainMethWGnCpG.diffmeth.temp)
MethCallDeconBrainIslandMethCpG.diffmeth=rbind(MethCallDeconBrainIslandMethCpG.diffmeth, MethCallDeconBrainIslandMethCpG.diffmeth.temp)
MethCallDeconBrainIslandMethNonCpG.diffmeth=rbind(MethCallDeconBrainIslandMethNonCpG.diffmeth, MethCallDeconBrainIslandMethNonCpG.diffmeth.temp)
MethCallDeconBrainGeneMethCpG.diffmeth=rbind(MethCallDeconBrainGeneMethCpG.diffmeth, MethCallDeconBrainGeneMethCpG.diffmeth.temp)
MethCallDeconBrainGeneMethNonCpG.diffmeth=rbind(MethCallDeconBrainGeneMethNonCpG.diffmeth,MethCallDeconBrainGeneMethNonCpG.diffmeth.temp)

MethCallDeconBrainMethWGoCpG.predict.temp=MethCallDeconBrainMethWGoCpG[which(MethCallDeconBrainMethWGoCpG$P==p),][which(p.adjust(MethCallDeconBrainMethWGoCpG[which(MethCallDeconBrainMethWGoCpG$P==p),]$PredictDiffMethP, method="none")<0.01),]
MethCallDeconBrainMethWGnCpG.predict.temp=MethCallDeconBrainMethWGnCpG[which(MethCallDeconBrainMethWGnCpG$P==p),][which(p.adjust(MethCallDeconBrainMethWGnCpG[which(MethCallDeconBrainMethWGnCpG$P==p),]$PredictDiffMethP, method="none")<0.01),]
MethCallDeconBrainIslandMethCpG.predict.temp=MethCallDeconBrainIslandMethCpG[which(MethCallDeconBrainIslandMethCpG$P==p),][which(p.adjust(MethCallDeconBrainIslandMethCpG[which(MethCallDeconBrainIslandMethCpG$P==p),]$PredictDiffMethP, method="none")<0.01),]
MethCallDeconBrainIslandMethNonCpG.predict.temp=MethCallDeconBrainIslandMethNonCpG[which(MethCallDeconBrainIslandMethNonCpG$P==p),][which(p.adjust(MethCallDeconBrainIslandMethNonCpG$PredictDiffMethP, method="none")<0.01),]
MethCallDeconBrainGeneMethCpG.predict.temp=MethCallDeconBrainGeneMethCpG[which(MethCallDeconBrainGeneMethCpG$P==p),][which(p.adjust(MethCallDeconBrainGeneMethCpG[which(MethCallDeconBrainGeneMethCpG$P==p),]$PredictDiffMethP, method="none")<0.01),]
MethCallDeconBrainGeneMethNonCpG.predict.temp=MethCallDeconBrainGeneMethNonCpG[which(MethCallDeconBrainGeneMethNonCpG$P==p),][which(p.adjust(MethCallDeconBrainGeneMethNonCpG[which(MethCallDeconBrainGeneMethNonCpG$P==p),]$PredictDiffMethP, method="none")<0.01),]
MethCallDeconBrainMethWGoCpG.predict=rbind(MethCallDeconBrainMethWGoCpG.predict, MethCallDeconBrainMethWGoCpG.predict.temp)
MethCallDeconBrainMethWGnCpG.predict=rbind(MethCallDeconBrainMethWGnCpG.predict, MethCallDeconBrainMethWGnCpG.predict.temp)
MethCallDeconBrainIslandMethCpG.predict=rbind(MethCallDeconBrainIslandMethCpG.predict, MethCallDeconBrainIslandMethCpG.predict.temp)
MethCallDeconBrainIslandMethNonCpG.predict=rbind(MethCallDeconBrainIslandMethNonCpG.predict, MethCallDeconBrainIslandMethNonCpG.predict.temp)
MethCallDeconBrainGeneMethCpG.predict=rbind(MethCallDeconBrainGeneMethCpG.predict, MethCallDeconBrainGeneMethCpG.predict.temp)
MethCallDeconBrainGeneMethNonCpG.predict=rbind(MethCallDeconBrainGeneMethNonCpG.predict,MethCallDeconBrainGeneMethNonCpG.predict.temp)

MethCallDeconBrainMethWGoCpG.tp_predict.temp=MethCallDeconBrainMethWGoCpG[which(MethCallDeconBrainMethWGoCpG$P==p),][which(p.adjust(MethCallDeconBrainMethWGoCpG[which(MethCallDeconBrainMethWGoCpG$P==p),]$tpPredictDiffMethP, method="none")<0.01),]
MethCallDeconBrainMethWGnCpG.tp_predict.temp=MethCallDeconBrainMethWGnCpG[which(MethCallDeconBrainMethWGnCpG$P==p),][which(p.adjust(MethCallDeconBrainMethWGnCpG[which(MethCallDeconBrainMethWGnCpG$P==p),]$tpPredictDiffMethP, method="none")<0.01),]
MethCallDeconBrainIslandMethCpG.tp_predict.temp=MethCallDeconBrainIslandMethCpG[which(MethCallDeconBrainIslandMethCpG$P==p),][which(p.adjust(MethCallDeconBrainIslandMethCpG[which(MethCallDeconBrainIslandMethCpG$P==p),]$tpPredictDiffMethP, method="none")<0.01),]
MethCallDeconBrainIslandMethNonCpG.tp_predict.temp=MethCallDeconBrainIslandMethNonCpG[which(MethCallDeconBrainIslandMethNonCpG$P==p),][which(p.adjust(MethCallDeconBrainIslandMethNonCpG[which(MethCallDeconBrainIslandMethNonCpG$P==p),]$tpPredictDiffMethP, method="none")<0.01),]
MethCallDeconBrainGeneMethCpG.tp_predict.temp=MethCallDeconBrainGeneMethCpG[which(MethCallDeconBrainGeneMethCpG$P==p),][which(p.adjust(MethCallDeconBrainGeneMethCpG[which(MethCallDeconBrainGeneMethCpG$P==p),]$tpPredictDiffMethP, method="none")<0.01),]
MethCallDeconBrainGeneMethNonCpG.tp_predict.temp=MethCallDeconBrainGeneMethNonCpG[which(MethCallDeconBrainGeneMethNonCpG$P==p),][which(p.adjust(MethCallDeconBrainGeneMethNonCpG[which(MethCallDeconBrainGeneMethNonCpG$P==p),]$tpPredictDiffMethP, method="none")<0.01),]
MethCallDeconBrainMethWGoCpG.tp_predict=rbind(MethCallDeconBrainMethWGoCpG.tp_predict, MethCallDeconBrainMethWGoCpG.tp_predict.temp)
MethCallDeconBrainMethWGnCpG.tp_predict=rbind(MethCallDeconBrainMethWGnCpG.tp_predict, MethCallDeconBrainMethWGnCpG.tp_predict.temp)
MethCallDeconBrainIslandMethCpG.tp_predict=rbind(MethCallDeconBrainIslandMethCpG.tp_predict, MethCallDeconBrainIslandMethCpG.tp_predict.temp)
MethCallDeconBrainIslandMethNonCpG.tp_predict=rbind(MethCallDeconBrainIslandMethNonCpG.tp_predict, MethCallDeconBrainIslandMethNonCpG.tp_predict.temp)
MethCallDeconBrainGeneMethCpG.tp_predict=rbind(MethCallDeconBrainGeneMethCpG.tp_predict, MethCallDeconBrainGeneMethCpG.tp_predict.temp)
MethCallDeconBrainGeneMethNonCpG.tp_predict=rbind(MethCallDeconBrainGeneMethNonCpG.tp_predict,MethCallDeconBrainGeneMethNonCpG.tp_predict.temp)
}

#Site Sensitivity and Specificity: CpG
SpecificitySensitivityCpG=as.data.frame(matrix(nrow=11, ncol=29, dimnames=list(c(0,10,20,30,40,50,60,70,80,90,100), c("P", "PredictPositives", "WrongSign", "PositivesTP", "WrongSignTP", "PositivesMix", "WrongSignMix", "TruePositives", "PredictNegatives", "NegativesTP", "NegativesMix", "TrueNegatives", "Sensitivity", "Specificity", "SensitivityTP", "SpecificityTP", "SensitivityMix", "SpecificityMix", "SensitivityTP_LC", "SensitivityMix_LC", "Sensitivity_UC", "SensitivityTP_UC", "SensitivityMix_UC", "Specificity_LC", "SpecificityTP_LC", "SpecificityMix_LC", "Specificity_UC", "SpecificityTP_UC", "SpecificityMix_UC"))))
count=0
for (p in c(0,10,20,30,40,50,60,70,80,90,100)){
count=count+1
P=BloodDeconP$EstimatedP[count]
MethCallDeconBrainMethWGoCpG.P=MethCallDeconBrainMethWGoCpG[which(MethCallDeconBrainMethWGoCpG$P==p),]
MethCallDeconBrainMethWGoCpG.diffmethP=MethCallDeconBrainMethWGoCpG.diffmeth[which(MethCallDeconBrainMethWGoCpG.diffmeth$P==p),]
SpecificitySensitivityCpG$P[count]=p
SpecificitySensitivityCpG$PredictPositives[count]=dim(MethCallDeconBrainMethWGoCpG.diffmethP[which(p.adjust(MethCallDeconBrainMethWGoCpG.diffmethP$PredictDiffMethP, method="none", n=dim(MethCallDeconBrainMethWGoCpG.P)[1])<0.01),])[1]
SpecificitySensitivityCpG$TruePositives[count]=dim(MethCallDeconBrainMethWGoCpG.diffmethP)[1]
SpecificitySensitivityCpG$WrongSign[count]=dim(MethCallDeconBrainMethWGoCpG.diffmethP[which(p.adjust(MethCallDeconBrainMethWGoCpG.diffmethP$PredictDiffMeth, method="none", n=dim(MethCallDeconBrainMethWGoCpG.P)[1])<0.01 & ((MethCallDeconBrainMethWGoCpG.diffmethP$DeconMethPercent < MethCallDeconBrainMethWGoCpG.diffmethP$BloodMethPercent & MethCallDeconBrainMethWGoCpG.diffmethP$BrainMethPercent > MethCallDeconBrainMethWGoCpG.diffmethP$BloodMethPercent) | (MethCallDeconBrainMethWGoCpG.diffmethP$DeconMethPercent > MethCallDeconBrainMethWGoCpG.diffmethP$BloodMethPercent & MethCallDeconBrainMethWGoCpG.diffmethP$BrainMethPercent < MethCallDeconBrainMethWGoCpG.diffmethP$BloodMethPercent))),])[1]
SpecificitySensitivityCpG$Sensitivity[count]=(SpecificitySensitivityCpG$PredictPositives[count]-SpecificitySensitivityCpG$WrongSign[count])/SpecificitySensitivityCpG$TruePositives[count]
SpecificitySensitivityCpG$PredictNegatives[count]=dim(MethCallDeconBrainMethWGoCpG.P[which(p.adjust(MethCallDeconBrainMethWGoCpG.P$PredictDiffMethP,  method="none")>0.01 & p.adjust(MethCallDeconBrainMethWGoCpG.P$TrueDiffMethP, method="none")>0.01),])[1]
SpecificitySensitivityCpG$TrueNegatives[count]=dim(MethCallDeconBrainMethWGoCpG.P[which(p.adjust(MethCallDeconBrainMethWGoCpG.P$TrueDiffMethP, method="none")>0.01), ])[1]
SpecificitySensitivityCpG$Specificity[count]=SpecificitySensitivityCpG$PredictNegatives[count]/SpecificitySensitivityCpG$TrueNegatives[count]

SpecificitySensitivityCpG$PositivesTP[count]=dim(MethCallDeconBrainMethWGoCpG.diffmethP[which(p.adjust(MethCallDeconBrainMethWGoCpG.diffmethP$tpPredictDiffMethP, method="none", n=dim(MethCallDeconBrainMethWGoCpG.P)[1])<0.01),])[1]
SpecificitySensitivityCpG$PositivesMix[count]=dim(MethCallDeconBrainMethWGoCpG.diffmethP[which(p.adjust(MethCallDeconBrainMethWGoCpG.diffmethP$MixDiffMethP, method="none", n=dim(MethCallDeconBrainMethWGoCpG.P)[1])<0.01),])[1]
SpecificitySensitivityCpG$WrongSignTP[count]=dim(MethCallDeconBrainMethWGoCpG.diffmethP[which(p.adjust(MethCallDeconBrainMethWGoCpG.diffmethP$tpPredictDiffMeth, method="none", n=dim(MethCallDeconBrainMethWGoCpG.P)[1])<0.01 & ((MethCallDeconBrainMethWGoCpG.diffmethP$tpMethPercent < MethCallDeconBrainMethWGoCpG.diffmethP$BloodMethPercent & MethCallDeconBrainMethWGoCpG.diffmethP$BrainMethPercent > MethCallDeconBrainMethWGoCpG.diffmethP$BloodMethPercent) | (MethCallDeconBrainMethWGoCpG.diffmethP$tpMethPercent > MethCallDeconBrainMethWGoCpG.diffmethP$BloodMethPercent & MethCallDeconBrainMethWGoCpG.diffmethP$BrainMethPercent < MethCallDeconBrainMethWGoCpG.diffmethP$BloodMethPercent))),])[1]
SpecificitySensitivityCpG$WrongSignMix[count]=dim(MethCallDeconBrainMethWGoCpG.diffmethP[which(p.adjust(MethCallDeconBrainMethWGoCpG.diffmethP$MixDiffMeth, method="none", n=dim(MethCallDeconBrainMethWGoCpG.P)[1])<0.01 & ((MethCallDeconBrainMethWGoCpG.diffmethP$MixMethPercent < MethCallDeconBrainMethWGoCpG.diffmethP$BloodMethPercent & MethCallDeconBrainMethWGoCpG.diffmethP$BrainMethPercent > MethCallDeconBrainMethWGoCpG.diffmethP$BloodMethPercent) | (MethCallDeconBrainMethWGoCpG.diffmethP$MixMethPercent > MethCallDeconBrainMethWGoCpG.diffmethP$BloodMethPercent & MethCallDeconBrainMethWGoCpG.diffmethP$BrainMethPercent < MethCallDeconBrainMethWGoCpG.diffmethP$BloodMethPercent))),])[1]
SpecificitySensitivityCpG$NegativesTP[count]=dim(MethCallDeconBrainMethWGoCpG.P[which(p.adjust(MethCallDeconBrainMethWGoCpG.P$tpPredictDiffMethP, method="none")>0.01 & p.adjust(MethCallDeconBrainMethWGoCpG.P$TrueDiffMethP, method="none")>0.01),])[1]
SpecificitySensitivityCpG$NegativesMix[count]=dim(MethCallDeconBrainMethWGoCpG.P[which(p.adjust(MethCallDeconBrainMethWGoCpG.P$MixDiffMethP, method="none")>0.01 & p.adjust(MethCallDeconBrainMethWGoCpG.P$TrueDiffMethP, method="none")>0.01),])[1]
SpecificitySensitivityCpG$SensitivityTP[count]=(SpecificitySensitivityCpG$PositivesTP[count]-SpecificitySensitivityCpG$WrongSignTP[count])/SpecificitySensitivityCpG$TruePositives[count]
SpecificitySensitivityCpG$SensitivityMix[count]=(SpecificitySensitivityCpG$PositivesMix[count]-SpecificitySensitivityCpG$WrongSignMix[count])/SpecificitySensitivityCpG$TruePositives[count]
SpecificitySensitivityCpG$SpecificityTP[count]=SpecificitySensitivityCpG$NegativesTP[count]/SpecificitySensitivityCpG$TrueNegatives[count]
SpecificitySensitivityCpG$SpecificityMix[count]=SpecificitySensitivityCpG$NegativesMix[count]/SpecificitySensitivityCpG$TrueNegatives[count]

SpecificitySensitivityCpG$Sensitivity_LC[count]=clopper.pearson.ci(SpecificitySensitivityCpG$PredictPositives[count], SpecificitySensitivityCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityCpG$SensitivityTP_LC[count]=clopper.pearson.ci(SpecificitySensitivityCpG$PositivesTP[count], SpecificitySensitivityCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityCpG$SensitivityMix_LC[count]=clopper.pearson.ci(SpecificitySensitivityCpG$PositivesMix[count], SpecificitySensitivityCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityCpG$Sensitivity_UC[count]=clopper.pearson.ci(SpecificitySensitivityCpG$PredictPositives[count], SpecificitySensitivityCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityCpG$SensitivityTP_UC[count]=clopper.pearson.ci(SpecificitySensitivityCpG$PositivesTP[count], SpecificitySensitivityCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityCpG$SensitivityMix_UC[count]=clopper.pearson.ci(SpecificitySensitivityCpG$PositivesMix[count], SpecificitySensitivityCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Upper.limit

SpecificitySensitivityCpG$Specificity_LC[count]=clopper.pearson.ci(SpecificitySensitivityCpG$PredictNegatives[count], SpecificitySensitivityCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityCpG$SpecificityTP_LC[count]=clopper.pearson.ci(SpecificitySensitivityCpG$NegativesTP[count], SpecificitySensitivityCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityCpG$SpecificityMix_LC[count]=clopper.pearson.ci(SpecificitySensitivityCpG$NegativesMix[count], SpecificitySensitivityCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityCpG$Specificity_UC[count]=clopper.pearson.ci(SpecificitySensitivityCpG$PredictNegatives[count], SpecificitySensitivityCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityCpG$SpecificityTP_UC[count]=clopper.pearson.ci(SpecificitySensitivityCpG$NegativesTP[count], SpecificitySensitivityCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityCpG$SpecificityMix_UC[count]=clopper.pearson.ci(SpecificitySensitivityCpG$NegativesMix[count], SpecificitySensitivityCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Upper.limit
}

#Site Sensitivity and Specificity: Non-CpG
SpecificitySensitivityNonCpG=as.data.frame(matrix(nrow=11, ncol=29, dimnames=list(c(0,10,20,30,40,50,60,70,80,90,100), c("P", "PredictPositives", "WrongSign", "PositivesTP", "WrongSignTP", "PositivesMix", "WrongSignMix", "TruePositives", "PredictNegatives", "NegativesTP", "NegativesMix", "TrueNegatives", "Sensitivity", "Specificity", "SensitivityTP", "SpecificityTP", "SensitivityMix", "SpecificityMix", "SensitivityTP_LC", "SensitivityMix_LC", "Sensitivity_UC", "SensitivityTP_UC", "SensitivityMix_UC", "Specificity_LC", "SpecificityTP_LC", "SpecificityMix_LC", "Specificity_UC", "SpecificityTP_UC", "SpecificityMix_UC"))))
count=0
for (p in c(0,10,20,30,40,50,60,70,80,90,100)){
count=count+1
P=BloodDeconP$EstimatedP[count]
MethCallDeconBrainMethWGnCpG.P=MethCallDeconBrainMethWGnCpG[which(MethCallDeconBrainMethWGnCpG$P==p),]
MethCallDeconBrainMethWGnCpG.diffmethP=MethCallDeconBrainMethWGnCpG.diffmeth[which(MethCallDeconBrainMethWGnCpG.diffmeth$P==p),]
SpecificitySensitivityNonCpG$P[count]=p
SpecificitySensitivityNonCpG$PredictPositives[count]=dim(MethCallDeconBrainMethWGnCpG.diffmethP[which(p.adjust(MethCallDeconBrainMethWGnCpG.diffmethP$PredictDiffMethP, method="none", n=dim(MethCallDeconBrainMethWGnCpG.P)[1])<0.01),])[1]
SpecificitySensitivityNonCpG$TruePositives[count]=dim(MethCallDeconBrainMethWGnCpG.diffmethP)[1]
SpecificitySensitivityNonCpG$WrongSign[count]=dim(MethCallDeconBrainMethWGnCpG.diffmethP[which(p.adjust(MethCallDeconBrainMethWGnCpG.diffmethP$PredictDiffMeth, method="none", n=dim(MethCallDeconBrainMethWGnCpG.P)[1])<0.01 & ((MethCallDeconBrainMethWGnCpG.diffmethP$DeconMethPercent < MethCallDeconBrainMethWGnCpG.diffmethP$BloodMethPercent & MethCallDeconBrainMethWGnCpG.diffmethP$BrainMethPercent > MethCallDeconBrainMethWGnCpG.diffmethP$BloodMethPercent) | (MethCallDeconBrainMethWGnCpG.diffmethP$DeconMethPercent > MethCallDeconBrainMethWGnCpG.diffmethP$BloodMethPercent & MethCallDeconBrainMethWGnCpG.diffmethP$BrainMethPercent < MethCallDeconBrainMethWGnCpG.diffmethP$BloodMethPercent))),])[1]
SpecificitySensitivityNonCpG$Sensitivity[count]=(SpecificitySensitivityNonCpG$PredictPositives[count]-SpecificitySensitivityNonCpG$WrongSign[count])/SpecificitySensitivityNonCpG$TruePositives[count]
SpecificitySensitivityNonCpG$PredictNegatives[count]=dim(MethCallDeconBrainMethWGnCpG.P[which(p.adjust(MethCallDeconBrainMethWGnCpG.P$PredictDiffMethP,  method="none")>0.01 & p.adjust(MethCallDeconBrainMethWGnCpG.P$TrueDiffMethP, method="none")>0.01),])[1]
SpecificitySensitivityNonCpG$TrueNegatives[count]=dim(MethCallDeconBrainMethWGnCpG.P[which(p.adjust(MethCallDeconBrainMethWGnCpG.P$TrueDiffMethP, method="none")>0.01), ])[1]
SpecificitySensitivityNonCpG$Specificity[count]=SpecificitySensitivityNonCpG$PredictNegatives[count]/SpecificitySensitivityNonCpG$TrueNegatives[count]

SpecificitySensitivityNonCpG$PositivesTP[count]=dim(MethCallDeconBrainMethWGnCpG.diffmethP[which(p.adjust(MethCallDeconBrainMethWGnCpG.diffmethP$tpPredictDiffMethP, method="none", n=dim(MethCallDeconBrainMethWGnCpG.P)[1])<0.01),])[1]
SpecificitySensitivityNonCpG$PositivesMix[count]=dim(MethCallDeconBrainMethWGnCpG.diffmethP[which(p.adjust(MethCallDeconBrainMethWGnCpG.diffmethP$MixDiffMethP, method="none", n=dim(MethCallDeconBrainMethWGnCpG.P)[1])<0.01),])[1]
SpecificitySensitivityNonCpG$WrongSignTP[count]=dim(MethCallDeconBrainMethWGnCpG.diffmethP[which(p.adjust(MethCallDeconBrainMethWGnCpG.diffmethP$tpPredictDiffMeth, method="none", n=dim(MethCallDeconBrainMethWGnCpG.P)[1])<0.01 & ((MethCallDeconBrainMethWGnCpG.diffmethP$tpMethPercent < MethCallDeconBrainMethWGnCpG.diffmethP$BloodMethPercent & MethCallDeconBrainMethWGnCpG.diffmethP$BrainMethPercent > MethCallDeconBrainMethWGnCpG.diffmethP$BloodMethPercent) | (MethCallDeconBrainMethWGnCpG.diffmethP$tpMethPercent > MethCallDeconBrainMethWGnCpG.diffmethP$BloodMethPercent & MethCallDeconBrainMethWGnCpG.diffmethP$BrainMethPercent < MethCallDeconBrainMethWGnCpG.diffmethP$BloodMethPercent))),])[1]
SpecificitySensitivityNonCpG$WrongSignMix[count]=dim(MethCallDeconBrainMethWGnCpG.diffmethP[which(p.adjust(MethCallDeconBrainMethWGnCpG.diffmethP$MixDiffMeth, method="none", n=dim(MethCallDeconBrainMethWGnCpG.P)[1])<0.01 & ((MethCallDeconBrainMethWGnCpG.diffmethP$MixMethPercent < MethCallDeconBrainMethWGnCpG.diffmethP$BloodMethPercent & MethCallDeconBrainMethWGnCpG.diffmethP$BrainMethPercent > MethCallDeconBrainMethWGnCpG.diffmethP$BloodMethPercent) | (MethCallDeconBrainMethWGnCpG.diffmethP$MixMethPercent > MethCallDeconBrainMethWGnCpG.diffmethP$BloodMethPercent & MethCallDeconBrainMethWGnCpG.diffmethP$BrainMethPercent < MethCallDeconBrainMethWGnCpG.diffmethP$BloodMethPercent))),])[1]
SpecificitySensitivityNonCpG$NegativesTP[count]=dim(MethCallDeconBrainMethWGnCpG.P[which(p.adjust(MethCallDeconBrainMethWGnCpG.P$tpPredictDiffMethP, method="none")>0.01 & p.adjust(MethCallDeconBrainMethWGnCpG.P$TrueDiffMethP, method="none")>0.01),])[1]
SpecificitySensitivityNonCpG$NegativesMix[count]=dim(MethCallDeconBrainMethWGnCpG.P[which(p.adjust(MethCallDeconBrainMethWGnCpG.P$MixDiffMethP, method="none")>0.01 & p.adjust(MethCallDeconBrainMethWGnCpG.P$TrueDiffMethP, method="none")>0.01),])[1]
SpecificitySensitivityNonCpG$SensitivityTP[count]=(SpecificitySensitivityNonCpG$PositivesTP[count]-SpecificitySensitivityNonCpG$WrongSignTP[count])/SpecificitySensitivityNonCpG$TruePositives[count]
SpecificitySensitivityNonCpG$SensitivityMix[count]=(SpecificitySensitivityNonCpG$PositivesMix[count]-SpecificitySensitivityNonCpG$WrongSignMix[count])/SpecificitySensitivityNonCpG$TruePositives[count]
SpecificitySensitivityNonCpG$SpecificityTP[count]=SpecificitySensitivityNonCpG$NegativesTP[count]/SpecificitySensitivityNonCpG$TrueNegatives[count]
SpecificitySensitivityNonCpG$SpecificityMix[count]=SpecificitySensitivityNonCpG$NegativesMix[count]/SpecificitySensitivityNonCpG$TrueNegatives[count]

SpecificitySensitivityNonCpG$Sensitivity_LC[count]=clopper.pearson.ci(SpecificitySensitivityNonCpG$PredictPositives[count], SpecificitySensitivityNonCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityNonCpG$SensitivityTP_LC[count]=clopper.pearson.ci(SpecificitySensitivityNonCpG$PositivesTP[count], SpecificitySensitivityNonCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityNonCpG$SensitivityMix_LC[count]=clopper.pearson.ci(SpecificitySensitivityNonCpG$PositivesMix[count], SpecificitySensitivityNonCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityNonCpG$Sensitivity_UC[count]=clopper.pearson.ci(SpecificitySensitivityNonCpG$PredictPositives[count], SpecificitySensitivityNonCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityNonCpG$SensitivityTP_UC[count]=clopper.pearson.ci(SpecificitySensitivityNonCpG$PositivesTP[count], SpecificitySensitivityNonCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityNonCpG$SensitivityMix_UC[count]=clopper.pearson.ci(SpecificitySensitivityNonCpG$PositivesMix[count], SpecificitySensitivityNonCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Upper.limit

SpecificitySensitivityNonCpG$Specificity_LC[count]=clopper.pearson.ci(SpecificitySensitivityNonCpG$PredictNegatives[count], SpecificitySensitivityNonCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityNonCpG$SpecificityTP_LC[count]=clopper.pearson.ci(SpecificitySensitivityNonCpG$NegativesTP[count], SpecificitySensitivityNonCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityNonCpG$SpecificityMix_LC[count]=clopper.pearson.ci(SpecificitySensitivityNonCpG$NegativesMix[count], SpecificitySensitivityNonCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityNonCpG$Specificity_UC[count]=clopper.pearson.ci(SpecificitySensitivityNonCpG$PredictNegatives[count], SpecificitySensitivityNonCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityNonCpG$SpecificityTP_UC[count]=clopper.pearson.ci(SpecificitySensitivityNonCpG$NegativesTP[count], SpecificitySensitivityNonCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityNonCpG$SpecificityMix_UC[count]=clopper.pearson.ci(SpecificitySensitivityNonCpG$NegativesMix[count], SpecificitySensitivityNonCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Upper.limit
}

#CpG Island Sensitivity and Specificity: CpG
SpecificitySensitivityIslandCpG=as.data.frame(matrix(nrow=11, ncol=29, dimnames=list(c(0,10,20,30,40,50,60,70,80,90,100), c("P", "PredictPositives", "WrongSign", "PositivesTP", "WrongSignTP", "PositivesMix", "WrongSignMix", "TruePositives", "PredictNegatives", "NegativesTP", "NegativesMix", "TrueNegatives", "Sensitivity", "Specificity", "SensitivityTP", "SpecificityTP", "SensitivityMix", "SpecificityMix", "SensitivityTP_LC", "SensitivityMix_LC", "Sensitivity_UC", "SensitivityTP_UC", "SensitivityMix_UC", "Specificity_LC", "SpecificityTP_LC", "SpecificityMix_LC", "Specificity_UC", "SpecificityTP_UC", "SpecificityMix_UC"))))
count=0
for (p in c(0,10,20,30,40,50,60,70,80,90,100)){
count=count+1
P=BloodDeconP$EstimatedP[count]
MethCallDeconBrainIslandMethCpG.P=MethCallDeconBrainIslandMethCpG[which(MethCallDeconBrainIslandMethCpG$P==p),]
MethCallDeconBrainIslandMethCpG.diffmethP=MethCallDeconBrainIslandMethCpG.diffmeth[which(MethCallDeconBrainIslandMethCpG.diffmeth$P==p),]
SpecificitySensitivityIslandCpG$P[count]=p
SpecificitySensitivityIslandCpG$PredictPositives[count]=dim(MethCallDeconBrainIslandMethCpG.diffmethP[which(p.adjust(MethCallDeconBrainIslandMethCpG.diffmethP$PredictDiffMethP, method="none", n=dim(MethCallDeconBrainIslandMethCpG.P)[1])<0.01),])[1]
SpecificitySensitivityIslandCpG$TruePositives[count]=dim(MethCallDeconBrainIslandMethCpG.diffmethP)[1]
SpecificitySensitivityIslandCpG$WrongSign[count]=dim(MethCallDeconBrainIslandMethCpG.diffmethP[which(p.adjust(MethCallDeconBrainIslandMethCpG.diffmethP$PredictDiffMeth, method="none", n=dim(MethCallDeconBrainIslandMethCpG.P)[1])<0.01 & ((MethCallDeconBrainIslandMethCpG.diffmethP$DeconMethPercent < MethCallDeconBrainIslandMethCpG.diffmethP$BloodMethPercent & MethCallDeconBrainIslandMethCpG.diffmethP$BrainMethPercent > MethCallDeconBrainIslandMethCpG.diffmethP$BloodMethPercent) | (MethCallDeconBrainIslandMethCpG.diffmethP$DeconMethPercent > MethCallDeconBrainIslandMethCpG.diffmethP$BloodMethPercent & MethCallDeconBrainIslandMethCpG.diffmethP$BrainMethPercent < MethCallDeconBrainIslandMethCpG.diffmethP$BloodMethPercent))),])[1]
SpecificitySensitivityIslandCpG$Sensitivity[count]=(SpecificitySensitivityIslandCpG$PredictPositives[count]-SpecificitySensitivityIslandCpG$WrongSign[count])/SpecificitySensitivityIslandCpG$TruePositives[count]
SpecificitySensitivityIslandCpG$PredictNegatives[count]=dim(MethCallDeconBrainIslandMethCpG.P[which(p.adjust(MethCallDeconBrainIslandMethCpG.P$PredictDiffMethP,  method="none")>0.01 & p.adjust(MethCallDeconBrainIslandMethCpG.P$TrueDiffMethP, method="none")>0.01),])[1]
SpecificitySensitivityIslandCpG$TrueNegatives[count]=dim(MethCallDeconBrainIslandMethCpG.P[which(p.adjust(MethCallDeconBrainIslandMethCpG.P$TrueDiffMethP, method="none")>0.01), ])[1]
SpecificitySensitivityIslandCpG$Specificity[count]=SpecificitySensitivityIslandCpG$PredictNegatives[count]/SpecificitySensitivityIslandCpG$TrueNegatives[count]

SpecificitySensitivityIslandCpG$PositivesTP[count]=dim(MethCallDeconBrainIslandMethCpG.diffmethP[which(p.adjust(MethCallDeconBrainIslandMethCpG.diffmethP$tpPredictDiffMethP, method="none", n=dim(MethCallDeconBrainIslandMethCpG.P)[1])<0.01),])[1]
SpecificitySensitivityIslandCpG$PositivesMix[count]=dim(MethCallDeconBrainIslandMethCpG.diffmethP[which(p.adjust(MethCallDeconBrainIslandMethCpG.diffmethP$MixDiffMethP, method="none", n=dim(MethCallDeconBrainIslandMethCpG.P)[1])<0.01),])[1]
SpecificitySensitivityIslandCpG$WrongSignTP[count]=dim(MethCallDeconBrainIslandMethCpG.diffmethP[which(p.adjust(MethCallDeconBrainIslandMethCpG.diffmethP$tpPredictDiffMeth, method="none", n=dim(MethCallDeconBrainIslandMethCpG.P)[1])<0.01 & ((MethCallDeconBrainIslandMethCpG.diffmethP$tpMethPercent < MethCallDeconBrainIslandMethCpG.diffmethP$BloodMethPercent & MethCallDeconBrainIslandMethCpG.diffmethP$BrainMethPercent > MethCallDeconBrainIslandMethCpG.diffmethP$BloodMethPercent) | (MethCallDeconBrainIslandMethCpG.diffmethP$tpMethPercent > MethCallDeconBrainIslandMethCpG.diffmethP$BloodMethPercent & MethCallDeconBrainIslandMethCpG.diffmethP$BrainMethPercent < MethCallDeconBrainIslandMethCpG.diffmethP$BloodMethPercent))),])[1]
SpecificitySensitivityIslandCpG$WrongSignMix[count]=dim(MethCallDeconBrainIslandMethCpG.diffmethP[which(p.adjust(MethCallDeconBrainIslandMethCpG.diffmethP$MixDiffMeth, method="none", n=dim(MethCallDeconBrainIslandMethCpG.P)[1])<0.01 & ((MethCallDeconBrainIslandMethCpG.diffmethP$MixMethPercent < MethCallDeconBrainIslandMethCpG.diffmethP$BloodMethPercent & MethCallDeconBrainIslandMethCpG.diffmethP$BrainMethPercent > MethCallDeconBrainIslandMethCpG.diffmethP$BloodMethPercent) | (MethCallDeconBrainIslandMethCpG.diffmethP$MixMethPercent > MethCallDeconBrainIslandMethCpG.diffmethP$BloodMethPercent & MethCallDeconBrainIslandMethCpG.diffmethP$BrainMethPercent < MethCallDeconBrainIslandMethCpG.diffmethP$BloodMethPercent))),])[1]
SpecificitySensitivityIslandCpG$NegativesTP[count]=dim(MethCallDeconBrainIslandMethCpG.P[which(p.adjust(MethCallDeconBrainIslandMethCpG.P$tpPredictDiffMethP, method="none")>0.01 & p.adjust(MethCallDeconBrainIslandMethCpG.P$TrueDiffMethP, method="none")>0.01),])[1]
SpecificitySensitivityIslandCpG$NegativesMix[count]=dim(MethCallDeconBrainIslandMethCpG.P[which(p.adjust(MethCallDeconBrainIslandMethCpG.P$MixDiffMethP, method="none")>0.01 & p.adjust(MethCallDeconBrainIslandMethCpG.P$TrueDiffMethP, method="none")>0.01),])[1]
SpecificitySensitivityIslandCpG$SensitivityTP[count]=(SpecificitySensitivityIslandCpG$PositivesTP[count]-SpecificitySensitivityIslandCpG$WrongSignTP[count])/SpecificitySensitivityIslandCpG$TruePositives[count]
SpecificitySensitivityIslandCpG$SensitivityMix[count]=(SpecificitySensitivityIslandCpG$PositivesMix[count]-SpecificitySensitivityIslandCpG$WrongSignMix[count])/SpecificitySensitivityIslandCpG$TruePositives[count]
SpecificitySensitivityIslandCpG$SpecificityTP[count]=SpecificitySensitivityIslandCpG$NegativesTP[count]/SpecificitySensitivityIslandCpG$TrueNegatives[count]
SpecificitySensitivityIslandCpG$SpecificityMix[count]=SpecificitySensitivityIslandCpG$NegativesMix[count]/SpecificitySensitivityIslandCpG$TrueNegatives[count]

SpecificitySensitivityIslandCpG$Sensitivity_LC[count]=clopper.pearson.ci(SpecificitySensitivityIslandCpG$PredictPositives[count], SpecificitySensitivityIslandCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityIslandCpG$SensitivityTP_LC[count]=clopper.pearson.ci(SpecificitySensitivityIslandCpG$PositivesTP[count], SpecificitySensitivityIslandCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityIslandCpG$SensitivityMix_LC[count]=clopper.pearson.ci(SpecificitySensitivityIslandCpG$PositivesMix[count], SpecificitySensitivityIslandCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityIslandCpG$Sensitivity_UC[count]=clopper.pearson.ci(SpecificitySensitivityIslandCpG$PredictPositives[count], SpecificitySensitivityIslandCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityIslandCpG$SensitivityTP_UC[count]=clopper.pearson.ci(SpecificitySensitivityIslandCpG$PositivesTP[count], SpecificitySensitivityIslandCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityIslandCpG$SensitivityMix_UC[count]=clopper.pearson.ci(SpecificitySensitivityIslandCpG$PositivesMix[count], SpecificitySensitivityIslandCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Upper.limit

SpecificitySensitivityIslandCpG$Specificity_LC[count]=clopper.pearson.ci(SpecificitySensitivityIslandCpG$PredictNegatives[count], SpecificitySensitivityIslandCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityIslandCpG$SpecificityTP_LC[count]=clopper.pearson.ci(SpecificitySensitivityIslandCpG$NegativesTP[count], SpecificitySensitivityIslandCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityIslandCpG$SpecificityMix_LC[count]=clopper.pearson.ci(SpecificitySensitivityIslandCpG$NegativesMix[count], SpecificitySensitivityIslandCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityIslandCpG$Specificity_UC[count]=clopper.pearson.ci(SpecificitySensitivityIslandCpG$PredictNegatives[count], SpecificitySensitivityIslandCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityIslandCpG$SpecificityTP_UC[count]=clopper.pearson.ci(SpecificitySensitivityIslandCpG$NegativesTP[count], SpecificitySensitivityIslandCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityIslandCpG$SpecificityMix_UC[count]=clopper.pearson.ci(SpecificitySensitivityIslandCpG$NegativesMix[count], SpecificitySensitivityIslandCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Upper.limit
}

#CpG Island Sensitivity and Specificity: Non-CpG
SpecificitySensitivityIslandNonCpG=as.data.frame(matrix(nrow=11, ncol=29, dimnames=list(c(0,10,20,30,40,50,60,70,80,90,100), c("P", "PredictPositives", "WrongSign", "PositivesTP", "WrongSignTP", "PositivesMix", "WrongSignMix", "TruePositives", "PredictNegatives", "NegativesTP", "NegativesMix", "TrueNegatives", "Sensitivity", "Specificity", "SensitivityTP", "SpecificityTP", "SensitivityMix", "SpecificityMix", "SensitivityTP_LC", "SensitivityMix_LC", "Sensitivity_UC", "SensitivityTP_UC", "SensitivityMix_UC", "Specificity_LC", "SpecificityTP_LC", "SpecificityMix_LC", "Specificity_UC", "SpecificityTP_UC", "SpecificityMix_UC"))))
count=0
for (p in c(0,10,20,30,40,50,60,70,80,90,100)){
count=count+1
P=BloodDeconP$EstimatedP[count]
MethCallDeconBrainIslandMethNonCpG.P=MethCallDeconBrainIslandMethNonCpG[which(MethCallDeconBrainIslandMethNonCpG$P==p),]
MethCallDeconBrainIslandMethNonCpG.diffmethP=MethCallDeconBrainIslandMethNonCpG.diffmeth[which(MethCallDeconBrainIslandMethNonCpG.diffmeth$P==p),]
SpecificitySensitivityIslandNonCpG$P[count]=p
SpecificitySensitivityIslandNonCpG$PredictPositives[count]=dim(MethCallDeconBrainIslandMethNonCpG.diffmethP[which(p.adjust(MethCallDeconBrainIslandMethNonCpG.diffmethP$PredictDiffMethP, method="none", n=dim(MethCallDeconBrainIslandMethNonCpG.P)[1])<0.01),])[1]
SpecificitySensitivityIslandNonCpG$TruePositives[count]=dim(MethCallDeconBrainIslandMethNonCpG.diffmethP)[1]
SpecificitySensitivityIslandNonCpG$WrongSign[count]=dim(MethCallDeconBrainIslandMethNonCpG.diffmethP[which(p.adjust(MethCallDeconBrainIslandMethNonCpG.diffmethP$PredictDiffMeth, method="none", n=dim(MethCallDeconBrainIslandMethNonCpG.P)[1])<0.01 & ((MethCallDeconBrainIslandMethNonCpG.diffmethP$DeconMethPercent < MethCallDeconBrainIslandMethNonCpG.diffmethP$BloodMethPercent & MethCallDeconBrainIslandMethNonCpG.diffmethP$BrainMethPercent > MethCallDeconBrainIslandMethNonCpG.diffmethP$BloodMethPercent) | (MethCallDeconBrainIslandMethNonCpG.diffmethP$DeconMethPercent > MethCallDeconBrainIslandMethNonCpG.diffmethP$BloodMethPercent & MethCallDeconBrainIslandMethNonCpG.diffmethP$BrainMethPercent < MethCallDeconBrainIslandMethNonCpG.diffmethP$BloodMethPercent))),])[1]
SpecificitySensitivityIslandNonCpG$Sensitivity[count]=(SpecificitySensitivityIslandNonCpG$PredictPositives[count]-SpecificitySensitivityIslandNonCpG$WrongSign[count])/SpecificitySensitivityIslandNonCpG$TruePositives[count]
SpecificitySensitivityIslandNonCpG$PredictNegatives[count]=dim(MethCallDeconBrainIslandMethNonCpG.P[which(p.adjust(MethCallDeconBrainIslandMethNonCpG.P$PredictDiffMethP,  method="none")>0.01 & p.adjust(MethCallDeconBrainIslandMethNonCpG.P$TrueDiffMethP, method="none")>0.01),])[1]
SpecificitySensitivityIslandNonCpG$TrueNegatives[count]=dim(MethCallDeconBrainIslandMethNonCpG.P[which(p.adjust(MethCallDeconBrainIslandMethNonCpG.P$TrueDiffMethP, method="none")>0.01), ])[1]
SpecificitySensitivityIslandNonCpG$Specificity[count]=SpecificitySensitivityIslandNonCpG$PredictNegatives[count]/SpecificitySensitivityIslandNonCpG$TrueNegatives[count]

SpecificitySensitivityIslandNonCpG$PositivesTP[count]=dim(MethCallDeconBrainIslandMethNonCpG.diffmethP[which(p.adjust(MethCallDeconBrainIslandMethNonCpG.diffmethP$tpPredictDiffMethP, method="none", n=dim(MethCallDeconBrainIslandMethNonCpG.P)[1])<0.01),])[1]
SpecificitySensitivityIslandNonCpG$PositivesMix[count]=dim(MethCallDeconBrainIslandMethNonCpG.diffmethP[which(p.adjust(MethCallDeconBrainIslandMethNonCpG.diffmethP$MixDiffMethP, method="none", n=dim(MethCallDeconBrainIslandMethNonCpG.P)[1])<0.01),])[1]
SpecificitySensitivityIslandNonCpG$WrongSignTP[count]=dim(MethCallDeconBrainIslandMethNonCpG.diffmethP[which(p.adjust(MethCallDeconBrainIslandMethNonCpG.diffmethP$tpPredictDiffMeth, method="none", n=dim(MethCallDeconBrainIslandMethNonCpG.P)[1])<0.01 & ((MethCallDeconBrainIslandMethNonCpG.diffmethP$tpMethPercent < MethCallDeconBrainIslandMethNonCpG.diffmethP$BloodMethPercent & MethCallDeconBrainIslandMethNonCpG.diffmethP$BrainMethPercent > MethCallDeconBrainIslandMethNonCpG.diffmethP$BloodMethPercent) | (MethCallDeconBrainIslandMethNonCpG.diffmethP$tpMethPercent > MethCallDeconBrainIslandMethNonCpG.diffmethP$BloodMethPercent & MethCallDeconBrainIslandMethNonCpG.diffmethP$BrainMethPercent < MethCallDeconBrainIslandMethNonCpG.diffmethP$BloodMethPercent))),])[1]
SpecificitySensitivityIslandNonCpG$WrongSignMix[count]=dim(MethCallDeconBrainIslandMethNonCpG.diffmethP[which(p.adjust(MethCallDeconBrainIslandMethNonCpG.diffmethP$MixDiffMeth, method="none", n=dim(MethCallDeconBrainIslandMethNonCpG.P)[1])<0.01 & ((MethCallDeconBrainIslandMethNonCpG.diffmethP$MixMethPercent < MethCallDeconBrainIslandMethNonCpG.diffmethP$BloodMethPercent & MethCallDeconBrainIslandMethNonCpG.diffmethP$BrainMethPercent > MethCallDeconBrainIslandMethNonCpG.diffmethP$BloodMethPercent) | (MethCallDeconBrainIslandMethNonCpG.diffmethP$MixMethPercent > MethCallDeconBrainIslandMethNonCpG.diffmethP$BloodMethPercent & MethCallDeconBrainIslandMethNonCpG.diffmethP$BrainMethPercent < MethCallDeconBrainIslandMethNonCpG.diffmethP$BloodMethPercent))),])[1]
SpecificitySensitivityIslandNonCpG$NegativesTP[count]=dim(MethCallDeconBrainIslandMethNonCpG.P[which(p.adjust(MethCallDeconBrainIslandMethNonCpG.P$tpPredictDiffMethP, method="none")>0.01 & p.adjust(MethCallDeconBrainIslandMethNonCpG.P$TrueDiffMethP, method="none")>0.01),])[1]
SpecificitySensitivityIslandNonCpG$NegativesMix[count]=dim(MethCallDeconBrainIslandMethNonCpG.P[which(p.adjust(MethCallDeconBrainIslandMethNonCpG.P$MixDiffMethP, method="none")>0.01 & p.adjust(MethCallDeconBrainIslandMethNonCpG.P$TrueDiffMethP, method="none")>0.01),])[1]
SpecificitySensitivityIslandNonCpG$SensitivityTP[count]=(SpecificitySensitivityIslandNonCpG$PositivesTP[count]-SpecificitySensitivityIslandNonCpG$WrongSignTP[count])/SpecificitySensitivityIslandNonCpG$TruePositives[count]
SpecificitySensitivityIslandNonCpG$SensitivityMix[count]=(SpecificitySensitivityIslandNonCpG$PositivesMix[count]-SpecificitySensitivityIslandNonCpG$WrongSignMix[count])/SpecificitySensitivityIslandNonCpG$TruePositives[count]
SpecificitySensitivityIslandNonCpG$SpecificityTP[count]=SpecificitySensitivityIslandNonCpG$NegativesTP[count]/SpecificitySensitivityIslandNonCpG$TrueNegatives[count]
SpecificitySensitivityIslandNonCpG$SpecificityMix[count]=SpecificitySensitivityIslandNonCpG$NegativesMix[count]/SpecificitySensitivityIslandNonCpG$TrueNegatives[count]

SpecificitySensitivityIslandNonCpG$Sensitivity_LC[count]=clopper.pearson.ci(SpecificitySensitivityIslandNonCpG$PredictPositives[count], SpecificitySensitivityIslandNonCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityIslandNonCpG$SensitivityTP_LC[count]=clopper.pearson.ci(SpecificitySensitivityIslandNonCpG$PositivesTP[count], SpecificitySensitivityIslandNonCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityIslandNonCpG$SensitivityMix_LC[count]=clopper.pearson.ci(SpecificitySensitivityIslandNonCpG$PositivesMix[count], SpecificitySensitivityIslandNonCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityIslandNonCpG$Sensitivity_UC[count]=clopper.pearson.ci(SpecificitySensitivityIslandNonCpG$PredictPositives[count], SpecificitySensitivityIslandNonCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityIslandNonCpG$SensitivityTP_UC[count]=clopper.pearson.ci(SpecificitySensitivityIslandNonCpG$PositivesTP[count], SpecificitySensitivityIslandNonCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityIslandNonCpG$SensitivityMix_UC[count]=clopper.pearson.ci(SpecificitySensitivityIslandNonCpG$PositivesMix[count], SpecificitySensitivityIslandNonCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Upper.limit

SpecificitySensitivityIslandNonCpG$Specificity_LC[count]=clopper.pearson.ci(SpecificitySensitivityIslandNonCpG$PredictNegatives[count], SpecificitySensitivityIslandNonCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityIslandNonCpG$SpecificityTP_LC[count]=clopper.pearson.ci(SpecificitySensitivityIslandNonCpG$NegativesTP[count], SpecificitySensitivityIslandNonCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityIslandNonCpG$SpecificityMix_LC[count]=clopper.pearson.ci(SpecificitySensitivityIslandNonCpG$NegativesMix[count], SpecificitySensitivityIslandNonCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityIslandNonCpG$Specificity_UC[count]=clopper.pearson.ci(SpecificitySensitivityIslandNonCpG$PredictNegatives[count], SpecificitySensitivityIslandNonCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityIslandNonCpG$SpecificityTP_UC[count]=clopper.pearson.ci(SpecificitySensitivityIslandNonCpG$NegativesTP[count], SpecificitySensitivityIslandNonCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityIslandNonCpG$SpecificityMix_UC[count]=clopper.pearson.ci(SpecificitySensitivityIslandNonCpG$NegativesMix[count], SpecificitySensitivityIslandNonCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Upper.limit
}

#Gene Sensitivity and Specificity: CpG
SpecificitySensitivityGeneCpG=as.data.frame(matrix(nrow=11, ncol=29, dimnames=list(c(0,10,20,30,40,50,60,70,80,90,100), c("P", "PredictPositives", "WrongSign", "PositivesTP", "WrongSignTP", "PositivesMix", "WrongSignMix", "TruePositives", "PredictNegatives", "NegativesTP", "NegativesMix", "TrueNegatives", "Sensitivity", "Specificity", "SensitivityTP", "SpecificityTP", "SensitivityMix", "SpecificityMix", "SensitivityTP_LC", "SensitivityMix_LC", "Sensitivity_UC", "SensitivityTP_UC", "SensitivityMix_UC", "Specificity_LC", "SpecificityTP_LC", "SpecificityMix_LC", "Specificity_UC", "SpecificityTP_UC", "SpecificityMix_UC"))))
count=0
for (p in c(0,10,20,30,40,50,60,70,80,90,100)){
count=count+1
P=BloodDeconP$EstimatedP[count]
MethCallDeconBrainGeneMethCpG.P=MethCallDeconBrainGeneMethCpG[which(MethCallDeconBrainGeneMethCpG$P==p),]
MethCallDeconBrainGeneMethCpG.diffmethP=MethCallDeconBrainGeneMethCpG.diffmeth[which(MethCallDeconBrainGeneMethCpG.diffmeth$P==p),]
SpecificitySensitivityGeneCpG$P[count]=p
SpecificitySensitivityGeneCpG$PredictPositives[count]=dim(MethCallDeconBrainGeneMethCpG.diffmethP[which(p.adjust(MethCallDeconBrainGeneMethCpG.diffmethP$PredictDiffMethP, method="none", n=dim(MethCallDeconBrainGeneMethCpG.P)[1])<0.01),])[1]
SpecificitySensitivityGeneCpG$TruePositives[count]=dim(MethCallDeconBrainGeneMethCpG.diffmethP)[1]
SpecificitySensitivityGeneCpG$WrongSign[count]=dim(MethCallDeconBrainGeneMethCpG.diffmethP[which(p.adjust(MethCallDeconBrainGeneMethCpG.diffmethP$PredictDiffMeth, method="none", n=dim(MethCallDeconBrainGeneMethCpG.P)[1])<0.01 & ((MethCallDeconBrainGeneMethCpG.diffmethP$DeconMethPercent < MethCallDeconBrainGeneMethCpG.diffmethP$BloodMethPercent & MethCallDeconBrainGeneMethCpG.diffmethP$BrainMethPercent > MethCallDeconBrainGeneMethCpG.diffmethP$BloodMethPercent) | (MethCallDeconBrainGeneMethCpG.diffmethP$DeconMethPercent > MethCallDeconBrainGeneMethCpG.diffmethP$BloodMethPercent & MethCallDeconBrainGeneMethCpG.diffmethP$BrainMethPercent < MethCallDeconBrainGeneMethCpG.diffmethP$BloodMethPercent))),])[1]
SpecificitySensitivityGeneCpG$Sensitivity[count]=(SpecificitySensitivityGeneCpG$PredictPositives[count]-SpecificitySensitivityGeneCpG$WrongSign[count])/SpecificitySensitivityGeneCpG$TruePositives[count]
SpecificitySensitivityGeneCpG$PredictNegatives[count]=dim(MethCallDeconBrainGeneMethCpG.P[which(p.adjust(MethCallDeconBrainGeneMethCpG.P$PredictDiffMethP,  method="none")>0.01 & p.adjust(MethCallDeconBrainGeneMethCpG.P$TrueDiffMethP, method="none")>0.01),])[1]
SpecificitySensitivityGeneCpG$TrueNegatives[count]=dim(MethCallDeconBrainGeneMethCpG.P[which(p.adjust(MethCallDeconBrainGeneMethCpG.P$TrueDiffMethP, method="none")>0.01), ])[1]
SpecificitySensitivityGeneCpG$Specificity[count]=SpecificitySensitivityGeneCpG$PredictNegatives[count]/SpecificitySensitivityGeneCpG$TrueNegatives[count]

SpecificitySensitivityGeneCpG$PositivesTP[count]=dim(MethCallDeconBrainGeneMethCpG.diffmethP[which(p.adjust(MethCallDeconBrainGeneMethCpG.diffmethP$tpPredictDiffMethP, method="none", n=dim(MethCallDeconBrainGeneMethCpG.P)[1])<0.01),])[1]
SpecificitySensitivityGeneCpG$PositivesMix[count]=dim(MethCallDeconBrainGeneMethCpG.diffmethP[which(p.adjust(MethCallDeconBrainGeneMethCpG.diffmethP$MixDiffMethP, method="none", n=dim(MethCallDeconBrainGeneMethCpG.P)[1])<0.01),])[1]
SpecificitySensitivityGeneCpG$WrongSignTP[count]=dim(MethCallDeconBrainGeneMethCpG.diffmethP[which(p.adjust(MethCallDeconBrainGeneMethCpG.diffmethP$tpPredictDiffMeth, method="none", n=dim(MethCallDeconBrainGeneMethCpG.P)[1])<0.01 & ((MethCallDeconBrainGeneMethCpG.diffmethP$tpMethPercent < MethCallDeconBrainGeneMethCpG.diffmethP$BloodMethPercent & MethCallDeconBrainGeneMethCpG.diffmethP$BrainMethPercent > MethCallDeconBrainGeneMethCpG.diffmethP$BloodMethPercent) | (MethCallDeconBrainGeneMethCpG.diffmethP$tpMethPercent > MethCallDeconBrainGeneMethCpG.diffmethP$BloodMethPercent & MethCallDeconBrainGeneMethCpG.diffmethP$BrainMethPercent < MethCallDeconBrainGeneMethCpG.diffmethP$BloodMethPercent))),])[1]
SpecificitySensitivityGeneCpG$WrongSignMix[count]=dim(MethCallDeconBrainGeneMethCpG.diffmethP[which(p.adjust(MethCallDeconBrainGeneMethCpG.diffmethP$MixDiffMeth, method="none", n=dim(MethCallDeconBrainGeneMethCpG.P)[1])<0.01 & ((MethCallDeconBrainGeneMethCpG.diffmethP$MixMethPercent < MethCallDeconBrainGeneMethCpG.diffmethP$BloodMethPercent & MethCallDeconBrainGeneMethCpG.diffmethP$BrainMethPercent > MethCallDeconBrainGeneMethCpG.diffmethP$BloodMethPercent) | (MethCallDeconBrainGeneMethCpG.diffmethP$MixMethPercent > MethCallDeconBrainGeneMethCpG.diffmethP$BloodMethPercent & MethCallDeconBrainGeneMethCpG.diffmethP$BrainMethPercent < MethCallDeconBrainGeneMethCpG.diffmethP$BloodMethPercent))),])[1]
SpecificitySensitivityGeneCpG$NegativesTP[count]=dim(MethCallDeconBrainGeneMethCpG.P[which(p.adjust(MethCallDeconBrainGeneMethCpG.P$tpPredictDiffMethP, method="none")>0.01 & p.adjust(MethCallDeconBrainGeneMethCpG.P$TrueDiffMethP, method="none")>0.01),])[1]
SpecificitySensitivityGeneCpG$NegativesMix[count]=dim(MethCallDeconBrainGeneMethCpG.P[which(p.adjust(MethCallDeconBrainGeneMethCpG.P$MixDiffMethP, method="none")>0.01 & p.adjust(MethCallDeconBrainGeneMethCpG.P$TrueDiffMethP, method="none")>0.01),])[1]
SpecificitySensitivityGeneCpG$SensitivityTP[count]=(SpecificitySensitivityGeneCpG$PositivesTP[count]-SpecificitySensitivityGeneCpG$WrongSignTP[count])/SpecificitySensitivityGeneCpG$TruePositives[count]
SpecificitySensitivityGeneCpG$SensitivityMix[count]=(SpecificitySensitivityGeneCpG$PositivesMix[count]-SpecificitySensitivityGeneCpG$WrongSignMix[count])/SpecificitySensitivityGeneCpG$TruePositives[count]
SpecificitySensitivityGeneCpG$SpecificityTP[count]=SpecificitySensitivityGeneCpG$NegativesTP[count]/SpecificitySensitivityGeneCpG$TrueNegatives[count]
SpecificitySensitivityGeneCpG$SpecificityMix[count]=SpecificitySensitivityGeneCpG$NegativesMix[count]/SpecificitySensitivityGeneCpG$TrueNegatives[count]

SpecificitySensitivityGeneCpG$Sensitivity_LC[count]=clopper.pearson.ci(SpecificitySensitivityGeneCpG$PredictPositives[count], SpecificitySensitivityGeneCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityGeneCpG$SensitivityTP_LC[count]=clopper.pearson.ci(SpecificitySensitivityGeneCpG$PositivesTP[count], SpecificitySensitivityGeneCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityGeneCpG$SensitivityMix_LC[count]=clopper.pearson.ci(SpecificitySensitivityGeneCpG$PositivesMix[count], SpecificitySensitivityGeneCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityGeneCpG$Sensitivity_UC[count]=clopper.pearson.ci(SpecificitySensitivityGeneCpG$PredictPositives[count], SpecificitySensitivityGeneCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityGeneCpG$SensitivityTP_UC[count]=clopper.pearson.ci(SpecificitySensitivityGeneCpG$PositivesTP[count], SpecificitySensitivityGeneCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityGeneCpG$SensitivityMix_UC[count]=clopper.pearson.ci(SpecificitySensitivityGeneCpG$PositivesMix[count], SpecificitySensitivityGeneCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Upper.limit

SpecificitySensitivityGeneCpG$Specificity_LC[count]=clopper.pearson.ci(SpecificitySensitivityGeneCpG$PredictNegatives[count], SpecificitySensitivityGeneCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityGeneCpG$SpecificityTP_LC[count]=clopper.pearson.ci(SpecificitySensitivityGeneCpG$NegativesTP[count], SpecificitySensitivityGeneCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityGeneCpG$SpecificityMix_LC[count]=clopper.pearson.ci(SpecificitySensitivityGeneCpG$NegativesMix[count], SpecificitySensitivityGeneCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityGeneCpG$Specificity_UC[count]=clopper.pearson.ci(SpecificitySensitivityGeneCpG$PredictNegatives[count], SpecificitySensitivityGeneCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityGeneCpG$SpecificityTP_UC[count]=clopper.pearson.ci(SpecificitySensitivityGeneCpG$NegativesTP[count], SpecificitySensitivityGeneCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityGeneCpG$SpecificityMix_UC[count]=clopper.pearson.ci(SpecificitySensitivityGeneCpG$NegativesMix[count], SpecificitySensitivityGeneCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Upper.limit
}

#Gene Sensitivity and Specificity: Non-CpG
SpecificitySensitivityGeneNonCpG=as.data.frame(matrix(nrow=11, ncol=29, dimnames=list(c(0,10,20,30,40,50,60,70,80,90,100), c("P", "PredictPositives", "WrongSign", "PositivesTP", "WrongSignTP", "PositivesMix", "WrongSignMix", "TruePositives", "PredictNegatives", "NegativesTP", "NegativesMix", "TrueNegatives", "Sensitivity", "Specificity", "SensitivityTP", "SpecificityTP", "SensitivityMix", "SpecificityMix", "SensitivityTP_LC", "SensitivityMix_LC", "Sensitivity_UC", "SensitivityTP_UC", "SensitivityMix_UC", "Specificity_LC", "SpecificityTP_LC", "SpecificityMix_LC", "Specificity_UC", "SpecificityTP_UC", "SpecificityMix_UC"))))
count=0
for (p in c(0,10,20,30,40,50,60,70,80,90,100)){
count=count+1
P=BloodDeconP$EstimatedP[count]
MethCallDeconBrainGeneMethNonCpG.P=MethCallDeconBrainGeneMethNonCpG[which(MethCallDeconBrainGeneMethNonCpG$P==p),]
MethCallDeconBrainGeneMethNonCpG.diffmethP=MethCallDeconBrainGeneMethNonCpG.diffmeth[which(MethCallDeconBrainGeneMethNonCpG.diffmeth$P==p),]
SpecificitySensitivityGeneNonCpG$P[count]=p
SpecificitySensitivityGeneNonCpG$PredictPositives[count]=dim(MethCallDeconBrainGeneMethNonCpG.diffmethP[which(p.adjust(MethCallDeconBrainGeneMethNonCpG.diffmethP$PredictDiffMethP, method="none", n=dim(MethCallDeconBrainGeneMethNonCpG.P)[1])<0.01),])[1]
SpecificitySensitivityGeneNonCpG$TruePositives[count]=dim(MethCallDeconBrainGeneMethNonCpG.diffmethP)[1]
SpecificitySensitivityGeneNonCpG$WrongSign[count]=dim(MethCallDeconBrainGeneMethNonCpG.diffmethP[which(p.adjust(MethCallDeconBrainGeneMethNonCpG.diffmethP$PredictDiffMeth, method="none", n=dim(MethCallDeconBrainGeneMethNonCpG.P)[1])<0.01 & ((MethCallDeconBrainGeneMethNonCpG.diffmethP$DeconMethPercent < MethCallDeconBrainGeneMethNonCpG.diffmethP$BloodMethPercent & MethCallDeconBrainGeneMethNonCpG.diffmethP$BrainMethPercent > MethCallDeconBrainGeneMethNonCpG.diffmethP$BloodMethPercent) | (MethCallDeconBrainGeneMethNonCpG.diffmethP$DeconMethPercent > MethCallDeconBrainGeneMethNonCpG.diffmethP$BloodMethPercent & MethCallDeconBrainGeneMethNonCpG.diffmethP$BrainMethPercent < MethCallDeconBrainGeneMethNonCpG.diffmethP$BloodMethPercent))),])[1]
SpecificitySensitivityGeneNonCpG$Sensitivity[count]=(SpecificitySensitivityGeneNonCpG$PredictPositives[count]-SpecificitySensitivityGeneNonCpG$WrongSign[count])/SpecificitySensitivityGeneNonCpG$TruePositives[count]
SpecificitySensitivityGeneNonCpG$PredictNegatives[count]=dim(MethCallDeconBrainGeneMethNonCpG.P[which(p.adjust(MethCallDeconBrainGeneMethNonCpG.P$PredictDiffMethP,  method="none")>0.01 & p.adjust(MethCallDeconBrainGeneMethNonCpG.P$TrueDiffMethP, method="none")>0.01),])[1]
SpecificitySensitivityGeneNonCpG$TrueNegatives[count]=dim(MethCallDeconBrainGeneMethNonCpG.P[which(p.adjust(MethCallDeconBrainGeneMethNonCpG.P$TrueDiffMethP, method="none")>0.01), ])[1]
SpecificitySensitivityGeneNonCpG$Specificity[count]=SpecificitySensitivityGeneNonCpG$PredictNegatives[count]/SpecificitySensitivityGeneNonCpG$TrueNegatives[count]

SpecificitySensitivityGeneNonCpG$PositivesTP[count]=dim(MethCallDeconBrainGeneMethNonCpG.diffmethP[which(p.adjust(MethCallDeconBrainGeneMethNonCpG.diffmethP$tpPredictDiffMethP, method="none", n=dim(MethCallDeconBrainGeneMethNonCpG.P)[1])<0.01),])[1]
SpecificitySensitivityGeneNonCpG$PositivesMix[count]=dim(MethCallDeconBrainGeneMethNonCpG.diffmethP[which(p.adjust(MethCallDeconBrainGeneMethNonCpG.diffmethP$MixDiffMethP, method="none", n=dim(MethCallDeconBrainGeneMethNonCpG.P)[1])<0.01),])[1]
SpecificitySensitivityGeneNonCpG$WrongSignTP[count]=dim(MethCallDeconBrainGeneMethNonCpG.diffmethP[which(p.adjust(MethCallDeconBrainGeneMethNonCpG.diffmethP$tpPredictDiffMeth, method="none", n=dim(MethCallDeconBrainGeneMethNonCpG.P)[1])<0.01 & ((MethCallDeconBrainGeneMethNonCpG.diffmethP$tpMethPercent < MethCallDeconBrainGeneMethNonCpG.diffmethP$BloodMethPercent & MethCallDeconBrainGeneMethNonCpG.diffmethP$BrainMethPercent > MethCallDeconBrainGeneMethNonCpG.diffmethP$BloodMethPercent) | (MethCallDeconBrainGeneMethNonCpG.diffmethP$tpMethPercent > MethCallDeconBrainGeneMethNonCpG.diffmethP$BloodMethPercent & MethCallDeconBrainGeneMethNonCpG.diffmethP$BrainMethPercent < MethCallDeconBrainGeneMethNonCpG.diffmethP$BloodMethPercent))),])[1]
SpecificitySensitivityGeneNonCpG$WrongSignMix[count]=dim(MethCallDeconBrainGeneMethNonCpG.diffmethP[which(p.adjust(MethCallDeconBrainGeneMethNonCpG.diffmethP$MixDiffMeth, method="none", n=dim(MethCallDeconBrainGeneMethNonCpG.P)[1])<0.01 & ((MethCallDeconBrainGeneMethNonCpG.diffmethP$MixMethPercent < MethCallDeconBrainGeneMethNonCpG.diffmethP$BloodMethPercent & MethCallDeconBrainGeneMethNonCpG.diffmethP$BrainMethPercent > MethCallDeconBrainGeneMethNonCpG.diffmethP$BloodMethPercent) | (MethCallDeconBrainGeneMethNonCpG.diffmethP$MixMethPercent > MethCallDeconBrainGeneMethNonCpG.diffmethP$BloodMethPercent & MethCallDeconBrainGeneMethNonCpG.diffmethP$BrainMethPercent < MethCallDeconBrainGeneMethNonCpG.diffmethP$BloodMethPercent))),])[1]
SpecificitySensitivityGeneNonCpG$NegativesTP[count]=dim(MethCallDeconBrainGeneMethNonCpG.P[which(p.adjust(MethCallDeconBrainGeneMethNonCpG.P$tpPredictDiffMethP, method="none")>0.01 & p.adjust(MethCallDeconBrainGeneMethNonCpG.P$TrueDiffMethP, method="none")>0.01),])[1]
SpecificitySensitivityGeneNonCpG$NegativesMix[count]=dim(MethCallDeconBrainGeneMethNonCpG.P[which(p.adjust(MethCallDeconBrainGeneMethNonCpG.P$MixDiffMethP, method="none")>0.01 & p.adjust(MethCallDeconBrainGeneMethNonCpG.P$TrueDiffMethP, method="none")>0.01),])[1]
SpecificitySensitivityGeneNonCpG$SensitivityTP[count]=(SpecificitySensitivityGeneNonCpG$PositivesTP[count]-SpecificitySensitivityGeneNonCpG$WrongSignTP[count])/SpecificitySensitivityGeneNonCpG$TruePositives[count]
SpecificitySensitivityGeneNonCpG$SensitivityMix[count]=(SpecificitySensitivityGeneNonCpG$PositivesMix[count]-SpecificitySensitivityGeneNonCpG$WrongSignMix[count])/SpecificitySensitivityGeneNonCpG$TruePositives[count]
SpecificitySensitivityGeneNonCpG$SpecificityTP[count]=SpecificitySensitivityGeneNonCpG$NegativesTP[count]/SpecificitySensitivityGeneNonCpG$TrueNegatives[count]
SpecificitySensitivityGeneNonCpG$SpecificityMix[count]=SpecificitySensitivityGeneNonCpG$NegativesMix[count]/SpecificitySensitivityGeneNonCpG$TrueNegatives[count]

SpecificitySensitivityGeneNonCpG$Sensitivity_LC[count]=clopper.pearson.ci(SpecificitySensitivityGeneNonCpG$PredictPositives[count], SpecificitySensitivityGeneNonCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityGeneNonCpG$SensitivityTP_LC[count]=clopper.pearson.ci(SpecificitySensitivityGeneNonCpG$PositivesTP[count], SpecificitySensitivityGeneNonCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityGeneNonCpG$SensitivityMix_LC[count]=clopper.pearson.ci(SpecificitySensitivityGeneNonCpG$PositivesMix[count], SpecificitySensitivityGeneNonCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityGeneNonCpG$Sensitivity_UC[count]=clopper.pearson.ci(SpecificitySensitivityGeneNonCpG$PredictPositives[count], SpecificitySensitivityGeneNonCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityGeneNonCpG$SensitivityTP_UC[count]=clopper.pearson.ci(SpecificitySensitivityGeneNonCpG$PositivesTP[count], SpecificitySensitivityGeneNonCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityGeneNonCpG$SensitivityMix_UC[count]=clopper.pearson.ci(SpecificitySensitivityGeneNonCpG$PositivesMix[count], SpecificitySensitivityGeneNonCpG$TruePositives[count], alpha=0.01, CI="two.sided")$Upper.limit

SpecificitySensitivityGeneNonCpG$Specificity_LC[count]=clopper.pearson.ci(SpecificitySensitivityGeneNonCpG$PredictNegatives[count], SpecificitySensitivityGeneNonCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityGeneNonCpG$SpecificityTP_LC[count]=clopper.pearson.ci(SpecificitySensitivityGeneNonCpG$NegativesTP[count], SpecificitySensitivityGeneNonCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityGeneNonCpG$SpecificityMix_LC[count]=clopper.pearson.ci(SpecificitySensitivityGeneNonCpG$NegativesMix[count], SpecificitySensitivityGeneNonCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Lower.limit
SpecificitySensitivityGeneNonCpG$Specificity_UC[count]=clopper.pearson.ci(SpecificitySensitivityGeneNonCpG$PredictNegatives[count], SpecificitySensitivityGeneNonCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityGeneNonCpG$SpecificityTP_UC[count]=clopper.pearson.ci(SpecificitySensitivityGeneNonCpG$NegativesTP[count], SpecificitySensitivityGeneNonCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Upper.limit
SpecificitySensitivityGeneNonCpG$SpecificityMix_UC[count]=clopper.pearson.ci(SpecificitySensitivityGeneNonCpG$NegativesMix[count], SpecificitySensitivityGeneNonCpG$TrueNegatives[count], alpha=0.01, CI="two.sided")$Upper.limit
}


SpecificitySensitivityCpG$Context="CpG"
SpecificitySensitivityNonCpG$Context="Non-CpG"
SpecificitySensitivityIslandCpG$Context="CpG Island: CpG"
SpecificitySensitivityIslandNonCpG$Context="CpG Island: Non-CpG"
SpecificitySensitivityGeneCpG$Context="Gene Body: CpG"
SpecificitySensitivityGeneNonCpG$Context="Gene Body: Non-CpG"

SpecificitySensitivityCpG$SensitivityRatio=SpecificitySensitivityCpG$Sensitivity/SpecificitySensitivityCpG$SensitivityMix
SpecificitySensitivityNonCpG$SensitivityRatio=SpecificitySensitivityNonCpG$Sensitivity/SpecificitySensitivityNonCpG$SensitivityMix
SpecificitySensitivityIslandCpG$SensitivityRatio=SpecificitySensitivityIslandCpG$Sensitivity/SpecificitySensitivityIslandCpG$SensitivityMix
SpecificitySensitivityIslandNonCpG$SensitivityRatio=SpecificitySensitivityIslandNonCpG$Sensitivity/SpecificitySensitivityIslandNonCpG$SensitivityMix
SpecificitySensitivityGeneCpG$SensitivityRatio=SpecificitySensitivityGeneCpG$Sensitivity/SpecificitySensitivityGeneCpG$SensitivityMix
SpecificitySensitivityGeneNonCpG$SensitivityRatio=SpecificitySensitivityGeneNonCpG$Sensitivity/SpecificitySensitivityGeneNonCpG$SensitivityMix

SpecificitySensitivityCpG$SpecificityRatio=SpecificitySensitivityCpG$Specificity/SpecificitySensitivityCpG$SpecificityMix
SpecificitySensitivityNonCpG$SpecificityRatio=SpecificitySensitivityNonCpG$Specificity/SpecificitySensitivityNonCpG$SpecificityMix
SpecificitySensitivityIslandCpG$SpecificityRatio=SpecificitySensitivityIslandCpG$Specificity/SpecificitySensitivityIslandCpG$SpecificityMix
SpecificitySensitivityIslandNonCpG$SpecificityRatio=SpecificitySensitivityIslandNonCpG$Specificity/SpecificitySensitivityIslandNonCpG$SpecificityMix
SpecificitySensitivityGeneCpG$SpecificityRatio=SpecificitySensitivityGeneCpG$Specificity/SpecificitySensitivityGeneCpG$SpecificityMix
SpecificitySensitivityGeneNonCpG$SpecificityRatio=SpecificitySensitivityGeneNonCpG$Specificity/SpecificitySensitivityGeneNonCpG$SpecificityMix

SpecificitySensitivityCpG$DetectionRatio=(SpecificitySensitivityCpG$Sensitivity*SpecificitySensitivityCpG$Specificity)/(SpecificitySensitivityCpG$SensitivityMix*SpecificitySensitivityCpG$SpecificityMix)
SpecificitySensitivityNonCpG$DetectionRatio=(SpecificitySensitivityNonCpG$Sensitivity*SpecificitySensitivityNonCpG$Specificity)/(SpecificitySensitivityNonCpG$SensitivityMix*SpecificitySensitivityNonCpG$SpecificityMix)
SpecificitySensitivityIslandCpG$DetectionRatio=(SpecificitySensitivityIslandCpG$Sensitivity*SpecificitySensitivityIslandCpG$Specificity)/(SpecificitySensitivityIslandCpG$SensitivityMix*SpecificitySensitivityIslandCpG$SpecificityMix)
SpecificitySensitivityIslandNonCpG$DetectionRatio=(SpecificitySensitivityIslandNonCpG$Sensitivity*SpecificitySensitivityIslandNonCpG$Specificity)/(SpecificitySensitivityIslandNonCpG$SensitivityMix*SpecificitySensitivityIslandNonCpG$SpecificityMix)
SpecificitySensitivityGeneCpG$DetectionRatio=(SpecificitySensitivityGeneCpG$Sensitivity*SpecificitySensitivityGeneCpG$Specificity)/(SpecificitySensitivityGeneCpG$SensitivityMix*SpecificitySensitivityGeneCpG$SpecificityMix)
SpecificitySensitivityGeneNonCpG$DetectionRatio=(SpecificitySensitivityGeneNonCpG$Sensitivity*SpecificitySensitivityGeneNonCpG$Specificity)/(SpecificitySensitivityGeneNonCpG$SensitivityMix*SpecificitySensitivityGeneNonCpG$SpecificityMix)

RatioSensitivitySpecificity=as.data.frame(matrix(nrow=dim(rbind(SpecificitySensitivityCpG, SpecificitySensitivityNonCpG, SpecificitySensitivityIslandCpG, SpecificitySensitivityIslandNonCpG, SpecificitySensitivityGeneCpG, SpecificitySensitivityGeneNonCpG))[1]*2, ncol=1, dimnames=list(1:(dim(rbind(SpecificitySensitivityCpG, SpecificitySensitivityNonCpG, SpecificitySensitivityIslandCpG, SpecificitySensitivityIslandNonCpG, SpecificitySensitivityGeneCpG, SpecificitySensitivityGeneNonCpG))[1]*2), "P")))
RatioSensitivitySpecificity$Calculation=rep(c(rep("No Correction", 11), rep("Correction", 11)), 6)
RatioSensitivitySpecificity$Category=paste(RatioSensitivitySpecificity$Context, RatioSensitivitySpecificity$Calculation, sep=":")
RatioSensitivitySpecificity$SensitivityRatio=c(SpecificitySensitivityCpG$SensitivityRatio, SpecificitySensitivityNonCpG$SensitivityRatio, SpecificitySensitivityIslandCpG$SensitivityRatio, SpecificitySensitivityIslandNonCpG$SensitivityRatio, SpecificitySensitivityGeneCpG$SensitivityRatio, SpecificitySensitivityGeneNonCpG$SensitivityRatio)
RatioSensitivitySpecificity$SpecificityRatio=c(SpecificitySensitivityCpG$SpecificityRatio, SpecificitySensitivityNonCpG$SpecificityRatio, SpecificitySensitivityIslandCpG$SpecificityRatio, SpecificitySensitivityIslandNonCpG$SpecificityRatio, SpecificitySensitivityGeneCpG$SpecificityRatio, SpecificitySensitivityGeneNonCpG$SpecificityRatio)
RatioSensitivitySpecificity$DetectionRatio=c(SpecificitySensitivityCpG$DetectionRatio, SpecificitySensitivityNonCpG$DetectionRatio, SpecificitySensitivityIslandCpG$DetectionRatio, SpecificitySensitivityIslandNonCpG$DetectionRatio, SpecificitySensitivityGeneCpG$DetectionRatio, SpecificitySensitivityGeneNonCpG$DetectionRatio)

ultSensitivitySpecificity=as.data.frame(matrix(nrow=dim(rbind(SpecificitySensitivityCpG, SpecificitySensitivityNonCpG, SpecificitySensitivityIslandCpG, SpecificitySensitivityIslandNonCpG, SpecificitySensitivityGeneCpG, SpecificitySensitivityGeneNonCpG))[1]*3, ncol=1, dimnames=list(1:(dim(rbind(SpecificitySensitivityCpG, SpecificitySensitivityNonCpG, SpecificitySensitivityIslandCpG, SpecificitySensitivityIslandNonCpG, SpecificitySensitivityGeneCpG, SpecificitySensitivityGeneNonCpG))[1]*3), "P")))
ultSensitivitySpecificity$P=c(rep(SpecificitySensitivityCpG$P, 3), rep(SpecificitySensitivityNonCpG$P, 3), rep(SpecificitySensitivityIslandCpG$P, 3), rep(SpecificitySensitivityIslandNonCpG$P, 3), rep(SpecificitySensitivityGeneCpG$P, 3), rep(SpecificitySensitivityGeneNonCpG$P, 3))
ultSensitivitySpecificity$Context=c(rep(SpecificitySensitivityCpG$Context, 3), rep(SpecificitySensitivityNonCpG$Context, 3), rep(SpecificitySensitivityIslandCpG$Context, 3), rep(SpecificitySensitivityIslandNonCpG$Context, 3), rep(SpecificitySensitivityGeneCpG$Context, 3), rep(SpecificitySensitivityGeneNonCpG$Context, 3))
ultSensitivitySpecificity$Estimate=rep(c(rep("Est. P", 11), rep("True P", 11), rep("Mix", 11)), 6)

ultSensitivitySpecificity$Sensitivity=c(c(SpecificitySensitivityCpG$Sensitivity, SpecificitySensitivityCpG$SensitivityTP, SpecificitySensitivityCpG$SensitivityMix), c(SpecificitySensitivityNonCpG$Sensitivity, SpecificitySensitivityNonCpG$SensitivityTP, SpecificitySensitivityNonCpG$SensitivityMix), c(SpecificitySensitivityIslandCpG$Sensitivity, SpecificitySensitivityIslandCpG$SensitivityTP, SpecificitySensitivityIslandCpG$SensitivityMix), c(SpecificitySensitivityIslandNonCpG$Sensitivity, SpecificitySensitivityIslandNonCpG$SensitivityTP, SpecificitySensitivityIslandNonCpG$SensitivityMix), c(SpecificitySensitivityGeneCpG$Sensitivity, SpecificitySensitivityGeneCpG$SensitivityTP, SpecificitySensitivityGeneCpG$SensitivityMix), c(SpecificitySensitivityGeneNonCpG$Sensitivity, SpecificitySensitivityGeneNonCpG$SensitivityTP, SpecificitySensitivityGeneNonCpG$SensitivityMix))
ultSensitivitySpecificity$Sensitivity_LC=c(c(SpecificitySensitivityCpG$Sensitivity_LC, SpecificitySensitivityCpG$SensitivityTP_LC, SpecificitySensitivityCpG$SensitivityMix_LC), c(SpecificitySensitivityNonCpG$Sensitivity_LC, SpecificitySensitivityNonCpG$SensitivityTP_LC, SpecificitySensitivityNonCpG$SensitivityMix_LC), c(SpecificitySensitivityIslandCpG$Sensitivity_LC, SpecificitySensitivityIslandCpG$SensitivityTP_LC, SpecificitySensitivityIslandCpG$SensitivityMix_LC), c(SpecificitySensitivityIslandNonCpG$Sensitivity_LC, SpecificitySensitivityIslandNonCpG$SensitivityTP_LC, SpecificitySensitivityIslandNonCpG$SensitivityMix_LC), c(SpecificitySensitivityGeneCpG$Sensitivity_LC, SpecificitySensitivityGeneCpG$SensitivityTP_LC, SpecificitySensitivityGeneCpG$SensitivityMix_LC), c(SpecificitySensitivityGeneNonCpG$Sensitivity_LC, SpecificitySensitivityGeneNonCpG$SensitivityTP_LC, SpecificitySensitivityGeneNonCpG$SensitivityMix_LC))
ultSensitivitySpecificity$Sensitivity_UC=c(c(SpecificitySensitivityCpG$Sensitivity_UC, SpecificitySensitivityCpG$SensitivityTP_UC, SpecificitySensitivityCpG$SensitivityMix_UC), c(SpecificitySensitivityNonCpG$Sensitivity_UC, SpecificitySensitivityNonCpG$SensitivityTP_UC, SpecificitySensitivityNonCpG$SensitivityMix_UC), c(SpecificitySensitivityIslandCpG$Sensitivity_UC, SpecificitySensitivityIslandCpG$SensitivityTP_UC, SpecificitySensitivityIslandCpG$SensitivityMix_UC), c(SpecificitySensitivityIslandNonCpG$Sensitivity_UC, SpecificitySensitivityIslandNonCpG$SensitivityTP_UC, SpecificitySensitivityIslandNonCpG$SensitivityMix_UC), c(SpecificitySensitivityGeneCpG$Sensitivity_UC, SpecificitySensitivityGeneCpG$SensitivityTP_UC, SpecificitySensitivityGeneCpG$SensitivityMix_UC), c(SpecificitySensitivityGeneNonCpG$Sensitivity_UC, SpecificitySensitivityGeneNonCpG$SensitivityTP_UC, SpecificitySensitivityGeneNonCpG$SensitivityMix_UC))

ultSensitivitySpecificity$Specificity=c(c(SpecificitySensitivityCpG$Specificity, SpecificitySensitivityCpG$SpecificityTP, SpecificitySensitivityCpG$SpecificityMix), c(SpecificitySensitivityNonCpG$Specificity, SpecificitySensitivityNonCpG$SpecificityTP, SpecificitySensitivityNonCpG$SpecificityMix), c(SpecificitySensitivityIslandCpG$Specificity, SpecificitySensitivityIslandCpG$SpecificityTP, SpecificitySensitivityIslandCpG$SpecificityMix), c(SpecificitySensitivityIslandNonCpG$Specificity, SpecificitySensitivityIslandNonCpG$SpecificityTP, SpecificitySensitivityIslandNonCpG$SpecificityMix), c(SpecificitySensitivityGeneCpG$Specificity, SpecificitySensitivityGeneCpG$SpecificityTP, SpecificitySensitivityGeneCpG$SpecificityMix), c(SpecificitySensitivityGeneNonCpG$Specificity, SpecificitySensitivityGeneNonCpG$SpecificityTP, SpecificitySensitivityGeneNonCpG$SpecificityMix))
ultSensitivitySpecificity$Specificity_LC=c(c(SpecificitySensitivityCpG$Specificity_LC, SpecificitySensitivityCpG$SpecificityTP_LC, SpecificitySensitivityCpG$SpecificityMix_LC), c(SpecificitySensitivityNonCpG$Specificity_LC, SpecificitySensitivityNonCpG$SpecificityTP_LC, SpecificitySensitivityNonCpG$SpecificityMix_LC), c(SpecificitySensitivityIslandCpG$Specificity_LC, SpecificitySensitivityIslandCpG$SpecificityTP_LC, SpecificitySensitivityIslandCpG$SpecificityMix_LC), c(SpecificitySensitivityIslandNonCpG$Specificity_LC, SpecificitySensitivityIslandNonCpG$SpecificityTP_LC, SpecificitySensitivityIslandNonCpG$SpecificityMix_LC), c(SpecificitySensitivityGeneCpG$Specificity_LC, SpecificitySensitivityGeneCpG$SpecificityTP_LC, SpecificitySensitivityGeneCpG$SpecificityMix_LC), c(SpecificitySensitivityGeneNonCpG$Specificity_LC, SpecificitySensitivityGeneNonCpG$SpecificityTP_LC, SpecificitySensitivityGeneNonCpG$SpecificityMix_LC))
ultSensitivitySpecificity$Specificity_UC=c(c(SpecificitySensitivityCpG$Specificity_UC, SpecificitySensitivityCpG$SpecificityTP_UC, SpecificitySensitivityCpG$SpecificityMix_UC), c(SpecificitySensitivityNonCpG$Specificity_UC, SpecificitySensitivityNonCpG$SpecificityTP_UC, SpecificitySensitivityNonCpG$SpecificityMix_UC), c(SpecificitySensitivityIslandCpG$Specificity_UC, SpecificitySensitivityIslandCpG$SpecificityTP_UC, SpecificitySensitivityIslandCpG$SpecificityMix_UC), c(SpecificitySensitivityIslandNonCpG$Specificity_UC, SpecificitySensitivityIslandNonCpG$SpecificityTP_UC, SpecificitySensitivityIslandNonCpG$SpecificityMix_UC), c(SpecificitySensitivityGeneCpG$Specificity_UC, SpecificitySensitivityGeneCpG$SpecificityTP_UC, SpecificitySensitivityGeneCpG$SpecificityMix_UC), c(SpecificitySensitivityGeneNonCpG$Specificity_UC, SpecificitySensitivityGeneNonCpG$SpecificityTP_UC, SpecificitySensitivityGeneNonCpG$SpecificityMix_UC))
ultSensitivitySpecificity$WrongSign=c(c(SpecificitySensitivityCpG$WrongSign, SpecificitySensitivityCpG$WrongSignTP, SpecificitySensitivityCpG$WrongSignMix), c(SpecificitySensitivityNonCpG$WrongSign, SpecificitySensitivityNonCpG$WrongSignTP, SpecificitySensitivityNonCpG$WrongSignMix), c(SpecificitySensitivityIslandCpG$WrongSign, SpecificitySensitivityIslandCpG$WrongSignTP, SpecificitySensitivityIslandCpG$WrongSignMix), c(SpecificitySensitivityIslandNonCpG$WrongSign, SpecificitySensitivityIslandNonCpG$WrongSignTP, SpecificitySensitivityIslandNonCpG$WrongSignMix), c(SpecificitySensitivityGeneCpG$WrongSign, SpecificitySensitivityGeneCpG$WrongSignTP, SpecificitySensitivityGeneCpG$WrongSignMix), c(SpecificitySensitivityGeneNonCpG$WrongSign, SpecificitySensitivityGeneNonCpG$WrongSignTP, SpecificitySensitivityGeneNonCpG$WrongSignMix))


SensitivityOverAWARENESS=aov(Sensitivity~P*Context*Estimate, data=ultSensitivitySpecificity)
SensitivityANOVA=anova(SensitivityOverAWARENESS)
summary(SensitivityOverAWARENESS)
summary.lm(lm(Sensitivity~P*Context*Estimate, data=ultSensitivitySpecificity))
#Partial R-square for P; Context; Deconvolultion
(anova(aov(Sensitivity~P, data=ultSensitivitySpecificity))$`Sum Sq`[2]-SensitivityANOVA$`Sum Sq`[8])/anova(aov(Sensitivity~P, data=ultSensitivitySpecificity))$`Sum Sq`[2]
(anova(aov(Sensitivity~Context, data=ultSensitivitySpecificity))$`Sum Sq`[2]-SensitivityANOVA$`Sum Sq`[8])/anova(aov(Sensitivity~Context, data=ultSensitivitySpecificity))$`Sum Sq`[2]
(anova(aov(Sensitivity~Estimate, data=ultSensitivitySpecificity))$`Sum Sq`[2]-SensitivityANOVA$`Sum Sq`[8])/anova(aov(Sensitivity~Estimate, data=ultSensitivitySpecificity))$`Sum Sq`[2]
TukeyHSD(SensitivityOverAWARENESS, "Context")
TukeyHSD(SensitivityOverAWARENESS, "Estimate")


SensitivityAWARENESS=aov(Sensitivity~P+Context+Estimate, data=ultSensitivitySpecificity)
SensitivityANOVA=anova(SensitivityAWARENESS)
summary(SensitivityAWARENESS)
summary.lm(lm(Sensitivity~P+Context+Estimate, data=ultSensitivitySpecificity))
#Partial R-square for P; Context; Deconvolultion
(anova(aov(Sensitivity~P, data=ultSensitivitySpecificity))$`Sum Sq`[2]-SensitivityANOVA$`Sum Sq`[4])/anova(aov(Sensitivity~P, data=ultSensitivitySpecificity))$`Sum Sq`[2]
(anova(aov(Sensitivity~Context, data=ultSensitivitySpecificity))$`Sum Sq`[2]-SensitivityANOVA$`Sum Sq`[4])/anova(aov(Sensitivity~Context, data=ultSensitivitySpecificity))$`Sum Sq`[2]
(anova(aov(Sensitivity~Estimate, data=ultSensitivitySpecificity))$`Sum Sq`[2]-SensitivityANOVA$`Sum Sq`[4])/anova(aov(Sensitivity~Estimate, data=ultSensitivitySpecificity))$`Sum Sq`[2]
TukeyHSD(SensitivityAWARENESS, "Context")
TukeyHSD(SensitivityAWARENESS, "Estimate")


SpecificityOverAWARENESS=aov(Specificity~P*Context*Estimate, data=ultSensitivitySpecificity)
SpecificityANOVA=anova(SpecificityOverAWARENESS)
summary(SpecificityOverAWARENESS)
summary.lm(lm(Specificity~P*Context*Estimate, data=ultSensitivitySpecificity))
#Partial R-square for P; Context; Deconvolultion
(anova(aov(Specificity~P, data=ultSensitivitySpecificity))$`Sum Sq`[2]-SpecificityANOVA$`Sum Sq`[8])/anova(aov(Specificity~P, data=ultSensitivitySpecificity))$`Sum Sq`[2]
(anova(aov(Specificity~Context, data=ultSensitivitySpecificity))$`Sum Sq`[2]-SpecificityANOVA$`Sum Sq`[8])/anova(aov(Specificity~Context, data=ultSensitivitySpecificity))$`Sum Sq`[2]
(anova(aov(Specificity~Estimate, data=ultSensitivitySpecificity))$`Sum Sq`[2]-SpecificityANOVA$`Sum Sq`[8])/anova(aov(Specificity~Estimate, data=ultSensitivitySpecificity))$`Sum Sq`[2]
TukeyHSD(SpecificityOverAWARENESS, "Context")
TukeyHSD(SpecificityOverAWARENESS, "Estimate")

SpecificityAWARENESS=aov(Specificity~P+Context+Estimate, data=ultSensitivitySpecificity)
SpecificityANOVA=anova(SpecificityAWARENESS)
summary(SpecificityAWARENESS)
summary.lm(lm(Specificity~P+Context+Estimate, data=ultSensitivitySpecificity))
#Partial R-square for P; Context; Deconvolultion
(anova(aov(Specificity~P, data=ultSensitivitySpecificity))$`Sum Sq`[2]-SpecificityANOVA$`Sum Sq`[4])/anova(aov(Specificity~P, data=ultSensitivitySpecificity))$`Sum Sq`[2]
(anova(aov(Specificity~Context, data=ultSensitivitySpecificity))$`Sum Sq`[2]-SpecificityANOVA$`Sum Sq`[4])/anova(aov(Specificity~Context, data=ultSensitivitySpecificity))$`Sum Sq`[2]
(anova(aov(Specificity~Estimate, data=ultSensitivitySpecificity))$`Sum Sq`[2]-SpecificityANOVA$`Sum Sq`[4])/anova(aov(Specificity~Estimate, data=ultSensitivitySpecificity))$`Sum Sq`[2]
TukeyHSD(SpecificityAWARENESS, "Context")
TukeyHSD(SpecificityAWARENESS, "Estimate")

ultSensitivitySpecificity$Detection=ultSensitivitySpecificity$Sensitivity*ultSensitivitySpecificity$Specificity

ultSensitivitySpecificity$Estimate=factor(ultSensitivitySpecificity$Estimate, levels=c("Mix", "True P", "Est. P"))
ultSensitivitySpecificity$Context=factor(ultSensitivitySpecificity$Context, levels=c("CpG", "Non-CpG", "CpG Island: CpG", "CpG Island: Non-CpG", "Gene Body: CpG", "Gene Body: Non-CpG"))
DeconEstimate=c(brewer.pal(3, "Dark2")[3], rev(brewer.pal(6, "Paired")))
DeconEstimate=c(rev(brewer.pal(10, "Paired")[10]), rev(brewer.pal(4, "Paired")[c(4,2)]))

ggplot(ultSensitivitySpecificity, aes(as.factor(P), WrongSign))+geom_bar(aes(fill=Estimate), stat="identity", position="dodge", width=0.8)+labs(x="Contamination Blood (%)")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")+scale_fill_manual(values=DeconEstimate)+facet_wrap(nrow=1, "Context")
ggsave("allAllWrongSign-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/SensitivitySpecificity", scale = 1, width = 3.5, height = 4, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(ultSensitivitySpecificity, aes(x=as.factor(P), y=Sensitivity, ymin=Sensitivity_LC, ymax=Sensitivity_UC, group=Estimate))+geom_bar(aes(fill=Estimate), color="black", stat="identity", position=position_dodge(width=0.8))+geom_errorbar(color="black", stat="identity", position="dodge", width=0.8)+labs(x="Contamination Blood (%)")+scale_y_continuous(breaks=c(seq(0, 1, 0.25)))+theme(text=element_text(size=14), axis.text.x=element_text(size=8), axis.text.y=element_text(size=5), strip.text.x=element_text(face="bold"), legend.position="none")+scale_fill_manual(values=DeconEstimate)+facet_wrap(nrow=6, "Context")
ggsave("allSensitivity-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/SensitivitySpecificity", scale = 1, width = 3.5, height = 4, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(ultSensitivitySpecificity, aes(x=as.factor(P), y=Specificity, ymin=Specificity_LC, ymax=Specificity_UC, group=Estimate))+geom_bar(aes(fill=Estimate), color="black", stat="identity", position=position_dodge(width=0.8))+geom_errorbar(color="black", stat="identity", position="dodge", width=0.8)+labs(x="Contamination Blood (%)")+scale_y_continuous(breaks=c(seq(0, 1, 0.25)))+theme(text=element_text(size=14), axis.text.x=element_text(size=8), axis.text.y=element_text(size=5), strip.text.x=element_text(face="bold"), legend.position="none")+scale_fill_manual(values=DeconEstimate)+facet_wrap(nrow=6, "Context")
ggsave("allSpecificity-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/SensitivitySpecificity", scale = 1, width = 3.5, height = 4, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(ultSensitivitySpecificity, aes(x=as.factor(P), y=Specificity, ymin=Specificity_LC, ymax=Specificity_UC, group=Estimate))+geom_bar(aes(fill=Estimate), color="black", stat="identity", position=position_dodge(width=0.8))+geom_errorbar(color="black", stat="identity", position="dodge", width=0.8)+labs(x="Contamination Blood (%)")+scale_y_continuous(breaks=c(seq(0, 1, 0.25)))+theme(text=element_text(size=14), axis.text.x=element_text(size=8), axis.text.y=element_text(size=5), strip.text.x=element_text(face="bold"), legend.position="bottom")+scale_fill_manual(values=DeconEstimate)+facet_wrap(nrow=6, "Context")
ggsave("allSpecificityLegend-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/SensitivitySpecificity", scale = 1, width = 7, height = 4, units = c("in"), dpi=300, limitsize = TRUE)
