#FiltAlleleFrequency Spectra and Summary Statistics in R

library(DescTools)
library(ggplot2)

AllFiltAlleleFreq=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Data/FiltAlleleFreq/Summary/sbS100X-SNV-AlleleFreq.txt", header=TRUE)
TrueHetFiltAlleleFreq=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Data/FiltAlleleFreq/Summary/sbS100X-SNV-TrueHetAlleleFreq.txt", header=TRUE)
ErrorFiltAlleleFreq=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Data/FiltAlleleFreq/Summary/sbS100X-SNV-ErrorAlleleFreq.txt", header=TRUE)

AllFiltAlleleFreq$NullAltFiltAlleleFreq=rbinom(length(AllFiltAlleleFreq$AlleleFreq), size=AllFiltAlleleFreq$Depth, prob=0.5)/AllFiltAlleleFreq$Depth
TrueHetFiltAlleleFreq$NullAltFiltAlleleFreq=rbinom(length(TrueHetFiltAlleleFreq$AlleleFreq), size=TrueHetFiltAlleleFreq$Depth, prob=0.5)/TrueHetFiltAlleleFreq$Depth

BinomAllFiltAlleleFreq=data.frame(matrix(nrow=10000, ncol=6, dimnames=list(1:10000, c("Iteration", "Mean", "Median", "StDev", "AbsDev", "Prop25"))))
BinomAllFiltAlleleFreq$sMean=mean(AllFiltAlleleFreq$AlleleFreq)
BinomAllFiltAlleleFreq$sMedian=median(AllFiltAlleleFreq$AlleleFreq)
BinomAllFiltAlleleFreq$sStDev=sd(AllFiltAlleleFreq$AlleleFreq)
BinomAllFiltAlleleFreq$sAbsDev=mad(AllFiltAlleleFreq$AlleleFreq)
BinomAllFiltAlleleFreq$sProp25=dim(AllFiltAlleleFreq[which(AllFiltAlleleFreq$AlleleFreq<0.25 | AllFiltAlleleFreq$AlleleFreq>0.75),])[1]/dim(AllFiltAlleleFreq)[1]
for(i in 1:10000){
SimBinomAltFiltAlleleFreq=rbinom(length(AllFiltAlleleFreq$AlleleFreq), size=AllFiltAlleleFreq$Depth, prob=0.5)/AllFiltAlleleFreq$Depth
BinomAllFiltAlleleFreq[i,]$Iteration=i
BinomAllFiltAlleleFreq[i,]$Mean=mean(SimBinomAltFiltAlleleFreq)
BinomAllFiltAlleleFreq[i,]$Median=median(SimBinomAltFiltAlleleFreq)
BinomAllFiltAlleleFreq[i,]$StDev=sd(SimBinomAltFiltAlleleFreq)
BinomAllFiltAlleleFreq[i,]$AbsDev=mad(SimBinomAltFiltAlleleFreq)
BinomAllFiltAlleleFreq[i,]$Prop25=length(SimBinomAltFiltAlleleFreq[which(SimBinomAltFiltAlleleFreq<0.25 | SimBinomAltFiltAlleleFreq>0.75)])/dim(AllFiltAlleleFreq)[1]
}

BinomTrueHetFiltAlleleFreq=data.frame(matrix(nrow=10000, ncol=6, dimnames=list(1:10000, c("Iteration", "Mean", "Median", "StDev", "AbsDev", "Prop25"))))
BinomTrueHetFiltAlleleFreq$sMean=mean(TrueHetFiltAlleleFreq$AlleleFreq)
BinomTrueHetFiltAlleleFreq$sMedian=median(TrueHetFiltAlleleFreq$AlleleFreq)
BinomTrueHetFiltAlleleFreq$sStDev=sd(TrueHetFiltAlleleFreq$AlleleFreq)
BinomTrueHetFiltAlleleFreq$sAbsDev=mad(TrueHetFiltAlleleFreq$AlleleFreq)
BinomTrueHetFiltAlleleFreq$sProp25=dim(TrueHetFiltAlleleFreq[which(TrueHetFiltAlleleFreq$AlleleFreq<0.25 | TrueHetFiltAlleleFreq$AlleleFreq>0.75),])[1]/dim(TrueHetFiltAlleleFreq)[1]
for(i in 1:10000){
SimBinomAltFiltAlleleFreq=rbinom(length(TrueHetFiltAlleleFreq$AlleleFreq), size=TrueHetFiltAlleleFreq$Depth, prob=0.5)/TrueHetFiltAlleleFreq$Depth
BinomTrueHetFiltAlleleFreq[i,]$Iteration=i
BinomTrueHetFiltAlleleFreq[i,]$Mean=mean(SimBinomAltFiltAlleleFreq)
BinomTrueHetFiltAlleleFreq[i,]$Median=median(SimBinomAltFiltAlleleFreq)
BinomTrueHetFiltAlleleFreq[i,]$StDev=sd(SimBinomAltFiltAlleleFreq)
BinomTrueHetFiltAlleleFreq[i,]$AbsDev=mad(SimBinomAltFiltAlleleFreq)
BinomTrueHetFiltAlleleFreq[i,]$Prop25=length(SimBinomAltFiltAlleleFreq[which(SimBinomAltFiltAlleleFreq<0.25 | SimBinomAltFiltAlleleFreq>0.75)])/dim(TrueHetFiltAlleleFreq)[1]
}

filtAllFiltAlleleFreq=AllFiltAlleleFreq[which(AllFiltAlleleFreq$Mutation=="CA" | AllFiltAlleleFreq$Mutation=="AC" | AllFiltAlleleFreq$Mutation=="AT" | AllFiltAlleleFreq$Mutation=="TA" | AllFiltAlleleFreq$Mutation=="GT" | AllFiltAlleleFreq$Mutation=="TG"),]
filtTrueHetFiltAlleleFreq=TrueHetFiltAlleleFreq[which(TrueHetFiltAlleleFreq$Mutation=="CA" | TrueHetFiltAlleleFreq$Mutation=="AC" | TrueHetFiltAlleleFreq$Mutation=="AT" | TrueHetFiltAlleleFreq$Mutation=="TA" | TrueHetFiltAlleleFreq$Mutation=="GT" | TrueHetFiltAlleleFreq$Mutation=="TG"),]
filtErrorFiltAlleleFreq=ErrorFiltAlleleFreq[which(ErrorFiltAlleleFreq$Mutation=="CA" | ErrorFiltAlleleFreq$Mutation=="AC" | ErrorFiltAlleleFreq$Mutation=="AT" | ErrorFiltAlleleFreq$Mutation=="TA" | ErrorFiltAlleleFreq$Mutation=="GT" | ErrorFiltAlleleFreq$Mutation=="TG"),]

filtAllFiltAlleleFreq$NullAltFiltAlleleFreq=rbinom(length(filtAllFiltAlleleFreq$AlleleFreq), size=filtAllFiltAlleleFreq$Depth, prob=0.5)/filtAllFiltAlleleFreq$Depth
filtTrueHetFiltAlleleFreq$NullAltFiltAlleleFreq=rbinom(length(filtTrueHetFiltAlleleFreq$AlleleFreq), size=filtTrueHetFiltAlleleFreq$Depth, prob=0.5)/filtTrueHetFiltAlleleFreq$Depth

BinomfiltAllFiltAlleleFreq=data.frame(matrix(nrow=10000, ncol=6, dimnames=list(1:10000, c("Iteration", "Mean", "Median", "StDev", "AbsDev", "Prop25"))))
BinomfiltAllFiltAlleleFreq$sMean=mean(filtAllFiltAlleleFreq$AlleleFreq)
BinomfiltAllFiltAlleleFreq$sMedian=median(filtAllFiltAlleleFreq$AlleleFreq)
BinomfiltAllFiltAlleleFreq$sStDev=sd(filtAllFiltAlleleFreq$AlleleFreq)
BinomfiltAllFiltAlleleFreq$sAbsDev=mad(filtAllFiltAlleleFreq$AlleleFreq)
BinomfiltAllFiltAlleleFreq$sProp25=dim(filtAllFiltAlleleFreq[which(filtAllFiltAlleleFreq$AlleleFreq<0.25 | filtAllFiltAlleleFreq$AlleleFreq>0.75),])[1]/dim(filtAllFiltAlleleFreq)[1]
for(i in 1:10000){
SimBinomAltFiltAlleleFreq=rbinom(length(filtAllFiltAlleleFreq$AlleleFreq), size=filtAllFiltAlleleFreq$Depth, prob=0.5)/filtAllFiltAlleleFreq$Depth
BinomfiltAllFiltAlleleFreq[i,]$Iteration=i
BinomfiltAllFiltAlleleFreq[i,]$Mean=mean(SimBinomAltFiltAlleleFreq)
BinomfiltAllFiltAlleleFreq[i,]$Median=median(SimBinomAltFiltAlleleFreq)
BinomfiltAllFiltAlleleFreq[i,]$StDev=sd(SimBinomAltFiltAlleleFreq)
BinomfiltAllFiltAlleleFreq[i,]$AbsDev=mad(SimBinomAltFiltAlleleFreq)
BinomfiltAllFiltAlleleFreq[i,]$Prop25=length(SimBinomAltFiltAlleleFreq[which(SimBinomAltFiltAlleleFreq<0.25 | SimBinomAltFiltAlleleFreq>0.75)])/dim(filtAllFiltAlleleFreq)[1]
}

BinomfiltTrueHetFiltAlleleFreq=data.frame(matrix(nrow=10000, ncol=6, dimnames=list(1:10000, c("Iteration", "Mean", "Median", "StDev", "AbsDev", "Prop25"))))
BinomfiltTrueHetFiltAlleleFreq$sMean=mean(filtTrueHetFiltAlleleFreq$AlleleFreq)
BinomfiltTrueHetFiltAlleleFreq$sMedian=median(filtTrueHetFiltAlleleFreq$AlleleFreq)
BinomfiltTrueHetFiltAlleleFreq$sStDev=sd(filtTrueHetFiltAlleleFreq$AlleleFreq)
BinomfiltTrueHetFiltAlleleFreq$sAbsDev=mad(filtTrueHetFiltAlleleFreq$AlleleFreq)
BinomfiltTrueHetFiltAlleleFreq$sProp25=dim(filtTrueHetFiltAlleleFreq[which(filtTrueHetFiltAlleleFreq$AlleleFreq<0.25 | filtTrueHetFiltAlleleFreq$AlleleFreq>0.75),])[1]/dim(filtTrueHetFiltAlleleFreq)[1]
for(i in 1:10000){
SimBinomAltFiltAlleleFreq=rbinom(length(filtTrueHetFiltAlleleFreq$AlleleFreq), size=filtTrueHetFiltAlleleFreq$Depth, prob=0.5)/filtTrueHetFiltAlleleFreq$Depth
BinomfiltTrueHetFiltAlleleFreq[i,]$Iteration=i
BinomfiltTrueHetFiltAlleleFreq[i,]$Mean=mean(SimBinomAltFiltAlleleFreq)
BinomfiltTrueHetFiltAlleleFreq[i,]$Median=median(SimBinomAltFiltAlleleFreq)
BinomfiltTrueHetFiltAlleleFreq[i,]$StDev=sd(SimBinomAltFiltAlleleFreq)
BinomfiltTrueHetFiltAlleleFreq[i,]$AbsDev=mad(SimBinomAltFiltAlleleFreq)
BinomfiltTrueHetFiltAlleleFreq[i,]$Prop25=length(SimBinomAltFiltAlleleFreq[which(SimBinomAltFiltAlleleFreq<0.25 | SimBinomAltFiltAlleleFreq>0.75)])/dim(filtTrueHetFiltAlleleFreq)[1]
}

write.table(BinomAllFiltAlleleFreq, file="C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Tables/BinomAllFiltAlleleFreq.txt", sep="\t")
write.table(BinomTrueHetFiltAlleleFreq, file="C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Tables/BinomTrueHetFiltAlleleFreq.txt", sep="\t")
write.table(BinomfiltAllFiltAlleleFreq, file="C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Tables/BinomfiltAllFiltAlleleFreq.txt", sep="\t")
write.table(BinomfiltTrueHetFiltAlleleFreq, file="C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Tables/BinomfiltTrueHetFiltAlleleFreq.txt", sep="\t")

ggplot(AllFiltAlleleFreq, aes(x=AlleleFreq))+geom_histogram(binwidth=0.01, aes(y=after_stat(density)))+geom_histogram(binwidth=0.01, aes(x=NullAltFiltAlleleFreq, y=after_stat(density)), alpha=0.0, color="red")+labs(x="Alternate Allele Frequency", y="SNVs (%)")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("AllFiltAlleleFreqHetSomHistogram.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(TrueHetFiltAlleleFreq, aes(x=AlleleFreq))+geom_histogram(binwidth=0.01, aes(y=after_stat(density)))+geom_histogram(binwidth=0.01, aes(x=NullAltFiltAlleleFreq, y=after_stat(density)), alpha=0.0, color="red")+labs(x="Alternate Allele Frequency", y="SNVs (%)")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("TrueHetFiltAlleleFreqHetSomHistogram.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(ErrorFiltAlleleFreq, aes(x=AlleleFreq))+geom_histogram(binwidth=0.01)+labs(x="Alternate Allele Frequency", y="SNVs (%)")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("ErrorFiltAlleleFreqHetSomHistogram.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(ErrorFiltAlleleFreq[which(ErrorFiltAlleleFreq$Depth>(0.5*mean(AllFiltAlleleFreq$Depth)) & ErrorFiltAlleleFreq$Depth<(2*mean(AllFiltAlleleFreq$Depth))),], aes(x=AlleleFreq))+geom_histogram(binwidth=0.01)+labs(x="Alternate Allele Frequency", y="SNVs (%)")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("dfErrorFiltAlleleFreqHetSomHistogram.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(filtAllFiltAlleleFreq, aes(x=AlleleFreq))+geom_histogram(binwidth=0.01, aes(y=after_stat(density)))+geom_histogram(binwidth=0.01, aes(x=NullAltFiltAlleleFreq, y=after_stat(density)), alpha=0.0, color="red")+labs(x="Alternate Allele Frequency", y="SNVs (%)")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("filtAllFiltAlleleFreqHetSomHistogram.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(filtTrueHetFiltAlleleFreq, aes(x=AlleleFreq))+geom_histogram(binwidth=0.01, aes(y=after_stat(density)))+geom_histogram(binwidth=0.01, aes(x=NullAltFiltAlleleFreq, y=after_stat(density)), alpha=0.0, color="red")+labs(x="Alternate Allele Frequency", y="SNVs (%)")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("filtTrueHetFiltAlleleFreqHetSomHistogram.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(filtErrorFiltAlleleFreq, aes(x=AlleleFreq))+geom_histogram(binwidth=0.01)+labs(x="Alternate Allele Frequency", y="SNVs (%)")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("filtErrorFiltAlleleFreqHetSomHistogram.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(filtErrorFiltAlleleFreq[which(ErrorFiltAlleleFreq$Depth>(0.5*mean(filtAllFiltAlleleFreq$Depth)) & ErrorFiltAlleleFreq$Depth<(2*mean(filtAllFiltAlleleFreq$Depth))),], aes(x=AlleleFreq))+geom_histogram(binwidth=0.01)+labs(x="Alternate Allele Frequency", y="SNVs (%)")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("filtdfErrorFiltAlleleFreqHetSomHistogram.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
dim(filtErrorFiltAlleleFreq[which(ErrorFiltAlleleFreq$Depth>(0.5*mean(filtAllFiltAlleleFreq$Depth)) & ErrorFiltAlleleFreq$Depth<(2*mean(filtAllFiltAlleleFreq$Depth))),])[1]
