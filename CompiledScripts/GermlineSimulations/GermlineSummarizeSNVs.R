#AlleleFrequency Spectra and Summary Statistics in R

library(DescTools)
library(ggplot2)

AllAlleleFreq=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Data/AlleleFreq/Summary/sbS100X-SNV-AlleleFreq.txt", header=TRUE)
TrueHetAlleleFreq=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Data/AlleleFreq/Summary/sbS100X-SNV-TrueHetAlleleFreq.txt", header=TRUE)
ErrorAlleleFreq=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Data/AlleleFreq/Summary/sbS100X-SNV-ErrorAlleleFreq.txt", header=TRUE)

AllAlleleFreq$NullAltAlleleFreq=rbinom(length(AllAlleleFreq$AlleleFreq), size=AllAlleleFreq$Depth, prob=0.5)/AllAlleleFreq$Depth
TrueHetAlleleFreq$NullAltAlleleFreq=rbinom(length(TrueHetAlleleFreq$AlleleFreq), size=TrueHetAlleleFreq$Depth, prob=0.5)/TrueHetAlleleFreq$Depth

BinomAllAlleleFreq=data.frame(matrix(nrow=10000, ncol=6, dimnames=list(1:10000, c("Iteration", "Mean", "Median", "StDev", "AbsDev", "Prop25"))))
BinomAllAlleleFreq$sMean=mean(AllAlleleFreq$AlleleFreq)
BinomAllAlleleFreq$sMedian=median(AllAlleleFreq$AlleleFreq)
BinomAllAlleleFreq$sStDev=sd(AllAlleleFreq$AlleleFreq)
BinomAllAlleleFreq$sAbsDev=mad(AllAlleleFreq$AlleleFreq)
BinomAllAlleleFreq$sProp25=dim(AllAlleleFreq[which(AllAlleleFreq$AlleleFreq<0.25 | AllAlleleFreq$AlleleFreq>0.75),])[1]/dim(AllAlleleFreq)[1]
for(i in 1:10000){
SimBinomAltAlleleFreq=rbinom(length(AllAlleleFreq$AlleleFreq), size=AllAlleleFreq$Depth, prob=0.5)/AllAlleleFreq$Depth
BinomAllAlleleFreq[i,]$Iteration=i
BinomAllAlleleFreq[i,]$Mean=mean(SimBinomAltAlleleFreq)
BinomAllAlleleFreq[i,]$Median=median(SimBinomAltAlleleFreq)
BinomAllAlleleFreq[i,]$StDev=sd(SimBinomAltAlleleFreq)
BinomAllAlleleFreq[i,]$AbsDev=mad(SimBinomAltAlleleFreq)
BinomAllAlleleFreq[i,]$Prop25=length(SimBinomAltAlleleFreq[which(SimBinomAltAlleleFreq<0.25 | SimBinomAltAlleleFreq>0.75)])/dim(AllAlleleFreq)[1]
}

BinomTrueHetAlleleFreq=data.frame(matrix(nrow=10000, ncol=6, dimnames=list(1:10000, c("Iteration", "Mean", "Median", "StDev", "AbsDev", "Prop25"))))
BinomTrueHetAlleleFreq$sMean=mean(TrueHetAlleleFreq$AlleleFreq)
BinomTrueHetAlleleFreq$sMedian=median(TrueHetAlleleFreq$AlleleFreq)
BinomTrueHetAlleleFreq$sStDev=sd(TrueHetAlleleFreq$AlleleFreq)
BinomTrueHetAlleleFreq$sAbsDev=mad(TrueHetAlleleFreq$AlleleFreq)
BinomTrueHetAlleleFreq$sProp25=dim(TrueHetAlleleFreq[which(TrueHetAlleleFreq$AlleleFreq<0.25 | TrueHetAlleleFreq$AlleleFreq>0.75),])[1]/dim(TrueHetAlleleFreq)[1]
for(i in 1:10000){
SimBinomAltAlleleFreq=rbinom(length(TrueHetAlleleFreq$AlleleFreq), size=TrueHetAlleleFreq$Depth, prob=0.5)/TrueHetAlleleFreq$Depth
BinomTrueHetAlleleFreq[i,]$Iteration=i
BinomTrueHetAlleleFreq[i,]$Mean=mean(SimBinomAltAlleleFreq)
BinomTrueHetAlleleFreq[i,]$Median=median(SimBinomAltAlleleFreq)
BinomTrueHetAlleleFreq[i,]$StDev=sd(SimBinomAltAlleleFreq)
BinomTrueHetAlleleFreq[i,]$AbsDev=mad(SimBinomAltAlleleFreq)
BinomTrueHetAlleleFreq[i,]$Prop25=length(SimBinomAltAlleleFreq[which(SimBinomAltAlleleFreq<0.25 | SimBinomAltAlleleFreq>0.75)])/dim(TrueHetAlleleFreq)[1]
}

filtAllAlleleFreq=AllAlleleFreq[which(AllAlleleFreq$Mutation=="CA" | AllAlleleFreq$Mutation=="AC" | AllAlleleFreq$Mutation=="AT" | AllAlleleFreq$Mutation=="TA" | AllAlleleFreq$Mutation=="GT" | AllAlleleFreq$Mutation=="TG"),]
filtTrueHetAlleleFreq=TrueHetAlleleFreq[which(TrueHetAlleleFreq$Mutation=="CA" | TrueHetAlleleFreq$Mutation=="AC" | TrueHetAlleleFreq$Mutation=="AT" | TrueHetAlleleFreq$Mutation=="TA" | TrueHetAlleleFreq$Mutation=="GT" | TrueHetAlleleFreq$Mutation=="TG"),]
filtErrorAlleleFreq=ErrorAlleleFreq[which(ErrorAlleleFreq$Mutation=="CA" | ErrorAlleleFreq$Mutation=="AC" | ErrorAlleleFreq$Mutation=="AT" | ErrorAlleleFreq$Mutation=="TA" | ErrorAlleleFreq$Mutation=="GT" | ErrorAlleleFreq$Mutation=="TG"),]

filtAllAlleleFreq$NullAltAlleleFreq=rbinom(length(filtAllAlleleFreq$AlleleFreq), size=filtAllAlleleFreq$Depth, prob=0.5)/filtAllAlleleFreq$Depth
filtTrueHetAlleleFreq$NullAltAlleleFreq=rbinom(length(filtTrueHetAlleleFreq$AlleleFreq), size=filtTrueHetAlleleFreq$Depth, prob=0.5)/filtTrueHetAlleleFreq$Depth

BinomfiltAllAlleleFreq=data.frame(matrix(nrow=10000, ncol=6, dimnames=list(1:10000, c("Iteration", "Mean", "Median", "StDev", "AbsDev", "Prop25"))))
BinomfiltAllAlleleFreq$sMean=mean(filtAllAlleleFreq$AlleleFreq)
BinomfiltAllAlleleFreq$sMedian=median(filtAllAlleleFreq$AlleleFreq)
BinomfiltAllAlleleFreq$sStDev=sd(filtAllAlleleFreq$AlleleFreq)
BinomfiltAllAlleleFreq$sAbsDev=mad(filtAllAlleleFreq$AlleleFreq)
BinomfiltAllAlleleFreq$sProp25=dim(filtAllAlleleFreq[which(filtAllAlleleFreq$AlleleFreq<0.25 | filtAllAlleleFreq$AlleleFreq>0.75),])[1]/dim(filtAllAlleleFreq)[1]
for(i in 1:10000){
SimBinomAltAlleleFreq=rbinom(length(filtAllAlleleFreq$AlleleFreq), size=filtAllAlleleFreq$Depth, prob=0.5)/filtAllAlleleFreq$Depth
BinomfiltAllAlleleFreq[i,]$Iteration=i
BinomfiltAllAlleleFreq[i,]$Mean=mean(SimBinomAltAlleleFreq)
BinomfiltAllAlleleFreq[i,]$Median=median(SimBinomAltAlleleFreq)
BinomfiltAllAlleleFreq[i,]$StDev=sd(SimBinomAltAlleleFreq)
BinomfiltAllAlleleFreq[i,]$AbsDev=mad(SimBinomAltAlleleFreq)
BinomfiltAllAlleleFreq[i,]$Prop25=length(SimBinomAltAlleleFreq[which(SimBinomAltAlleleFreq<0.25 | SimBinomAltAlleleFreq>0.75)])/dim(filtAllAlleleFreq)[1]
}

BinomfiltTrueHetAlleleFreq=data.frame(matrix(nrow=10000, ncol=6, dimnames=list(1:10000, c("Iteration", "Mean", "Median", "StDev", "AbsDev", "Prop25"))))
BinomfiltTrueHetAlleleFreq$sMean=mean(filtTrueHetAlleleFreq$AlleleFreq)
BinomfiltTrueHetAlleleFreq$sMedian=median(filtTrueHetAlleleFreq$AlleleFreq)
BinomfiltTrueHetAlleleFreq$sStDev=sd(filtTrueHetAlleleFreq$AlleleFreq)
BinomfiltTrueHetAlleleFreq$sAbsDev=mad(filtTrueHetAlleleFreq$AlleleFreq)
BinomfiltTrueHetAlleleFreq$sProp25=dim(filtTrueHetAlleleFreq[which(filtTrueHetAlleleFreq$AlleleFreq<0.25 | filtTrueHetAlleleFreq$AlleleFreq>0.75),])[1]/dim(filtTrueHetAlleleFreq)[1]
for(i in 1:10000){
SimBinomAltAlleleFreq=rbinom(length(filtTrueHetAlleleFreq$AlleleFreq), size=filtTrueHetAlleleFreq$Depth, prob=0.5)/filtTrueHetAlleleFreq$Depth
BinomfiltTrueHetAlleleFreq[i,]$Iteration=i
BinomfiltTrueHetAlleleFreq[i,]$Mean=mean(SimBinomAltAlleleFreq)
BinomfiltTrueHetAlleleFreq[i,]$Median=median(SimBinomAltAlleleFreq)
BinomfiltTrueHetAlleleFreq[i,]$StDev=sd(SimBinomAltAlleleFreq)
BinomfiltTrueHetAlleleFreq[i,]$AbsDev=mad(SimBinomAltAlleleFreq)
BinomfiltTrueHetAlleleFreq[i,]$Prop25=length(SimBinomAltAlleleFreq[which(SimBinomAltAlleleFreq<0.25 | SimBinomAltAlleleFreq>0.75)])/dim(filtTrueHetAlleleFreq)[1]
}

write.table(BinomAllAlleleFreq, file="C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Tables/BinomAllAlleleFreq.txt", sep="\t")
write.table(BinomTrueHetAlleleFreq, file="C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Tables/BinomTrueHetAlleleFreq.txt", sep="\t")
write.table(BinomfiltAllAlleleFreq, file="C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Tables/BinomfiltAllAlleleFreq.txt", sep="\t")
write.table(BinomfiltTrueHetAlleleFreq, file="C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Tables/BinomfiltTrueHetAlleleFreq.txt", sep="\t")

ggplot(AllAlleleFreq, aes(x=AlleleFreq))+geom_histogram(binwidth=0.01, aes(y=after_stat(density)))+geom_histogram(binwidth=0.01, aes(x=NullAltAlleleFreq, y=after_stat(density)), alpha=0.0, color="red")+labs(x="Alternate Allele Frequency", y="SNVs (%)")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("AllAlleleFreqHetSomHistogram.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(TrueHetAlleleFreq, aes(x=AlleleFreq))+geom_histogram(binwidth=0.01, aes(y=after_stat(density)))+geom_histogram(binwidth=0.01, aes(x=NullAltAlleleFreq, y=after_stat(density)), alpha=0.0, color="red")+labs(x="Alternate Allele Frequency", y="SNVs (%)")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("TrueHetAlleleFreqHetSomHistogram.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(ErrorAlleleFreq, aes(x=AlleleFreq))+geom_histogram(binwidth=0.01)+labs(x="Alternate Allele Frequency", y="SNVs (%)")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("ErrorAlleleFreqHetSomHistogram.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(ErrorAlleleFreq[which(ErrorAlleleFreq$Depth>(0.5*mean(AllAlleleFreq$Depth)) & ErrorAlleleFreq$Depth<(2*mean(AllAlleleFreq$Depth))),], aes(x=AlleleFreq))+geom_histogram(binwidth=0.01)+labs(x="Alternate Allele Frequency", y="SNVs (%)")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("dfErrorAlleleFreqHetSomHistogram.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(filtAllAlleleFreq, aes(x=AlleleFreq))+geom_histogram(binwidth=0.01, aes(y=after_stat(density)))+geom_histogram(binwidth=0.01, aes(x=NullAltAlleleFreq, y=after_stat(density)), alpha=0.0, color="red")+labs(x="Alternate Allele Frequency", y="SNVs (%)")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("filtAllAlleleFreqHetSomHistogram.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(filtTrueHetAlleleFreq, aes(x=AlleleFreq))+geom_histogram(binwidth=0.01, aes(y=after_stat(density)))+geom_histogram(binwidth=0.01, aes(x=NullAltAlleleFreq, y=after_stat(density)), alpha=0.0, color="red")+labs(x="Alternate Allele Frequency", y="SNVs (%)")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("filtTrueHetAlleleFreqHetSomHistogram.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(filtErrorAlleleFreq, aes(x=AlleleFreq))+geom_histogram(binwidth=0.01)+labs(x="Alternate Allele Frequency", y="SNVs (%)")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("filtErrorAlleleFreqHetSomHistogram.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(filtErrorAlleleFreq[which(ErrorAlleleFreq$Depth>(0.5*mean(filtAllAlleleFreq$Depth)) & ErrorAlleleFreq$Depth<(2*mean(filtAllAlleleFreq$Depth))),], aes(x=AlleleFreq))+geom_histogram(binwidth=0.01)+labs(x="Alternate Allele Frequency", y="SNVs (%)")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("filtdfErrorAlleleFreqHetSomHistogram.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
dim(filtErrorAlleleFreq[which(ErrorAlleleFreq$Depth>(0.5*mean(filtAllAlleleFreq$Depth)) & ErrorAlleleFreq$Depth<(2*mean(filtAllAlleleFreq$Depth))),])[1]
