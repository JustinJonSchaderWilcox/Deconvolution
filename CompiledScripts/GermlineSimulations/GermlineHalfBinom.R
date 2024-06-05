#Showing Half-Binomial Distribution Works as approximate
library(DescTools)
library(ggplot2)

AllFiltAlleleFreq=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Data/FiltAlleleFreq/Summary/sbS100X-SNV-AlleleFreq.txt", header=TRUE)
TrueHetFiltAlleleFreq=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Data/FiltAlleleFreq/Summary/sbS100X-SNV-TrueHetAlleleFreq.txt", header=TRUE)
ErrorFiltAlleleFreq=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Data/FiltAlleleFreq/Summary/sbS100X-SNV-ErrorAlleleFreq.txt", header=TRUE)

AllFiltAlleleFreq$NullAltFiltAlleleFreq=rbinom(length(AllFiltAlleleFreq$AlleleFreq), size=round((AllFiltAlleleFreq$Depth/2)+0.5), prob=mean(AllFiltAlleleFreq$AlleleFreq))/round((AllFiltAlleleFreq$Depth/2)+0.5)
TrueHetFiltAlleleFreq$NullAltFiltAlleleFreq=rbinom(length(TrueHetFiltAlleleFreq$AlleleFreq), size=round((TrueHetFiltAlleleFreq$Depth/2)+0.5), prob=mean(TrueHetFiltAlleleFreq$AlleleFreq))/round((TrueHetFiltAlleleFreq$Depth/2)+0.5)

HalfBinomAllFiltAlleleFreq=data.frame(matrix(nrow=10000, ncol=6, dimnames=list(1:10000, c("Iteration", "Mean", "Median", "StDev", "AbsDev", "Prop25"))))
HalfBinomAllFiltAlleleFreq$sMean=mean(AllFiltAlleleFreq$AlleleFreq)
HalfBinomAllFiltAlleleFreq$sMedian=median(AllFiltAlleleFreq$AlleleFreq)
HalfBinomAllFiltAlleleFreq$sStDev=sd(AllFiltAlleleFreq$AlleleFreq)
HalfBinomAllFiltAlleleFreq$sAbsDev=mad(AllFiltAlleleFreq$AlleleFreq)
HalfBinomAllFiltAlleleFreq$sProp25=dim(AllFiltAlleleFreq[which(AllFiltAlleleFreq$AlleleFreq<0.25 | AllFiltAlleleFreq$AlleleFreq>0.75),])[1]/dim(AllFiltAlleleFreq)[1]
for(i in 1:10000){
SimBinomAltFiltAlleleFreq=rbinom(length(AllFiltAlleleFreq$AlleleFreq), size=round((AllFiltAlleleFreq$Depth/2)+0.5), prob=0.5)/round((AllFiltAlleleFreq$Depth/2)+0.5)
HalfBinomAllFiltAlleleFreq[i,]$Iteration=i
HalfBinomAllFiltAlleleFreq[i,]$Mean=mean(SimBinomAltFiltAlleleFreq)
HalfBinomAllFiltAlleleFreq[i,]$Median=median(SimBinomAltFiltAlleleFreq)
HalfBinomAllFiltAlleleFreq[i,]$StDev=sd(SimBinomAltFiltAlleleFreq)
HalfBinomAllFiltAlleleFreq[i,]$AbsDev=mad(SimBinomAltFiltAlleleFreq)
HalfBinomAllFiltAlleleFreq[i,]$Prop25=length(SimBinomAltFiltAlleleFreq[which(SimBinomAltFiltAlleleFreq<0.25 | SimBinomAltFiltAlleleFreq>0.75)])/dim(AllFiltAlleleFreq)[1]
}

HalfBinomTrueHetFiltAlleleFreq=data.frame(matrix(nrow=10000, ncol=6, dimnames=list(1:10000, c("Iteration", "Mean", "Median", "StDev", "AbsDev", "Prop25"))))
HalfBinomTrueHetFiltAlleleFreq$sMean=mean(TrueHetFiltAlleleFreq$AlleleFreq)
HalfBinomTrueHetFiltAlleleFreq$sMedian=median(TrueHetFiltAlleleFreq$AlleleFreq)
HalfBinomTrueHetFiltAlleleFreq$sStDev=sd(TrueHetFiltAlleleFreq$AlleleFreq)
HalfBinomTrueHetFiltAlleleFreq$sAbsDev=mad(TrueHetFiltAlleleFreq$AlleleFreq)
HalfBinomTrueHetFiltAlleleFreq$sProp25=dim(TrueHetFiltAlleleFreq[which(TrueHetFiltAlleleFreq$AlleleFreq<0.25 | TrueHetFiltAlleleFreq$AlleleFreq>0.75),])[1]/dim(TrueHetFiltAlleleFreq)[1]
for(i in 1:10000){
SimBinomAltFiltAlleleFreq=rbinom(length(TrueHetFiltAlleleFreq$AlleleFreq), size=round((TrueHetFiltAlleleFreq$Depth/2)+0.5), prob=0.5)/round((TrueHetFiltAlleleFreq$Depth/2)+0.5)
HalfBinomTrueHetFiltAlleleFreq[i,]$Iteration=i
HalfBinomTrueHetFiltAlleleFreq[i,]$Mean=mean(SimBinomAltFiltAlleleFreq)
HalfBinomTrueHetFiltAlleleFreq[i,]$Median=median(SimBinomAltFiltAlleleFreq)
HalfBinomTrueHetFiltAlleleFreq[i,]$StDev=sd(SimBinomAltFiltAlleleFreq)
HalfBinomTrueHetFiltAlleleFreq[i,]$AbsDev=mad(SimBinomAltFiltAlleleFreq)
HalfBinomTrueHetFiltAlleleFreq[i,]$Prop25=length(SimBinomAltFiltAlleleFreq[which(SimBinomAltFiltAlleleFreq<0.25 | SimBinomAltFiltAlleleFreq>0.75)])/dim(TrueHetFiltAlleleFreq)[1]
}

dim(HalfBinomTrueHetFiltAlleleFreq[which(HalfBinomTrueHetFiltAlleleFreq$Prop25>HalfBinomTrueHetFiltAlleleFreq$sProp25),])[1]
wilcox.test(TrueHetFiltAlleleFreq$AlleleFreq, TrueHetFiltAlleleFreq$NullAltFiltAlleleFreq, paired=TRUE)
dbinom(round(dim(TrueHetFiltAlleleFreq[which(TrueHetFiltAlleleFreq$AlleleFreq<0.25 | TrueHetFiltAlleleFreq$AlleleFreq>0.75),])[1]/dim(TrueHetFiltAlleleFreq)[1]*length(TrueHetFiltAlleleFreq$AlleleFreq)+0.5), size=length(TrueHetFiltAlleleFreq$AlleleFreq), prob=dbinom(round((mean(TrueHetFiltAlleleFreq$Depth)/8)+0.5), size=round((mean(TrueHetFiltAlleleFreq$Depth)/2)+0.5), prob=mean(TrueHetFiltAlleleFreq$AlleleFreq))*2)

write.table(HalfBinomAllFiltAlleleFreq, file="C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Tables/HalfBinomAllFiltAlleleFreq.txt", sep="\t")
write.table(HalfBinomTrueHetFiltAlleleFreq, file="C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Tables/HalfBinomTrueHetFiltAlleleFreq.txt", sep="\t")

ggplot(AllFiltAlleleFreq, aes(x=AlleleFreq))+geom_histogram(binwidth=0.01, aes(y=after_stat(density)))+geom_histogram(binwidth=0.01, aes(x=NullAltFiltAlleleFreq, y=after_stat(density)), alpha=0.0, color="red")+labs(x="Alternate Allele Frequency", y="SNVs (%)")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("HalfBinomAllFiltAlleleFreqHetSomHistogram.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(TrueHetFiltAlleleFreq, aes(x=AlleleFreq))+geom_histogram(binwidth=0.01, aes(y=after_stat(density)))+geom_histogram(binwidth=0.01, aes(x=NullAltFiltAlleleFreq, y=after_stat(density)), alpha=0.0, color="red")+labs(x="Alternate Allele Frequency", y="SNVs (%)")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("HalfBinomTrueHetFiltAlleleFreqHetSomHistogram.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/HeterozygoteSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
