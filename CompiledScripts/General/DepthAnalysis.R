#Ultimate Depth Analysis
library(ggplot2)
library(RColorBrewer)

Deduplications=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Data/DeduplicationExtractor/AllP-DeduplicationSEE.txt")

ggplot(Deduplications, aes(x=as.factor(P), y=ALIGNMENTS))+geom_col()+labs(x="Contamination P (%)", y="Sequence Alignments")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("BismarkAlignments-x-P(Count).png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Images", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(Deduplications, aes(x=as.factor(P), y=ALIGNMENTS/Deduplications[which(Deduplications$P==100),]$ALIGNMENTS))+geom_col()+labs(x="Contamination P (%)", y="Sequence Alignments")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("BismarkAlignments-x-P(Prop).png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Images", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)

ggplot(Deduplications, aes(x=as.factor(P), y=DUPLICATES))+geom_col()+labs(x="Contamination P (%)", y="Duplicated Alignments")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("BismarkDuplicates-x-P(Count).png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Images", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(Deduplications, aes(x=as.factor(P), y=DUPLICATES/Deduplications[which(Deduplications$P==100),]$DUPLICATES))+geom_col()+labs(x="Contamination P (%)", y="Duplicated Alignments")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("BismarkDuplicates-x-P(Prop).png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Images", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)

ggplot(Deduplications, aes(x=as.factor(P), y=POSITIONS))+geom_col()+labs(x="Contamination P (%)", y="Deduplicated Positions")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("BismarkDuplicatePositions-x-P(Count).png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Images", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(Deduplications, aes(x=as.factor(P), y=POSITIONS/Deduplications[which(Deduplications$P==100),]$POSITIONS))+geom_col()+labs(x="Contamination P (%)", y="Deduplicated Positions")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("BismarkDuplicatePositions-x-P(Prop).png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Images", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)

ggplot(Deduplications, aes(x=as.factor(P), y=100*DUPLICATES/ALIGNMENTS))+geom_col()+labs(x="Contamination P (%)", y="Duplicated Alignments (%)")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("BismarkDuplicatesPercent-x-P(Prop).png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Images", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)

ggplot(Deduplications, aes(x=DUPLICATES, y=POSITIONS, groups=as.factor(P)))+geom_point(aes(color=(as.factor(P))), size=6, shape=18)+labs(x="Duplicates", y="Duplicated Positions")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"))+scale_color_manual(values=rev(brewer.pal("RdYlGn", n=11)))+guides(color=guide_legend(title="Blood %"))+geom_abline(intercept=0, slope=1, linewidth=1.2, linetype="dotted")
ggsave("BismarkPositionsDuplicated-x-Duplicates.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Images", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(Deduplications, aes(x=as.factor(P), y=POSITIONS/DUPLICATES, groups=as.factor(P)))+geom_col()+labs(x="Contamination P (%)", y="Duplicated: Positions per Alignment")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"))
ggsave("BismarkPositionsDuplicated-per-Duplicates.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Images", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)


ggplot(Deduplications, aes(x=as.factor(P), y=SEQUENCES))+geom_col()+labs(x="Contamination P (%)", y="Remaining Alignments")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("BismarkRemainingSequences-x-P(Count).png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Images", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(Deduplications, aes(x=as.factor(P), y=SEQUENCES/Deduplications[which(Deduplications$P==100),]$SEQUENCES))+geom_col()+labs(x="Contamination P (%)", y="Remaining Alignments")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"), legend.position="none")
ggsave("BismarkRemainingSequences-x-P(Prop).png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Images", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)

ggplot(Deduplications, aes(x=ALIGNMENTS, y=SEQUENCES, groups=as.factor(P)))+geom_point(aes(color=(as.factor(P))), size=6, shape=18)+labs(x="Alignments", y="Remaining Sequences")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"))+scale_color_manual(values=rev(brewer.pal("RdYlGn", n=11)))+guides(color=guide_legend(title="Blood %"))+geom_abline(intercept=0, slope=1, linewidth=1.2, linetype="dotted")
ggsave("BismarkRemainingSequences-x-Alignments.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Images", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)


cor.test(Deduplications$ALIGNMENTS, Deduplications$P)
cor.test(Deduplications$ALIGNMENTS, Deduplications$P)$estimate^2
cor.test(Deduplications$DUPLICATES, -1*abs((Deduplications$P)*(100-Deduplications$P)))
cor.test(Deduplications$DUPLICATES, -1*abs((Deduplications$P)*(100-Deduplications$P)))$estimate^2
cor.test(Deduplications$POSITIONS, -1*abs((Deduplications$P)*(100-Deduplications$P)))
cor.test(Deduplications$POSITIONS, -1*abs((Deduplications$P)*(100-Deduplications$P)))$estimate^2
cor.test(Deduplications$POSITIONS, Deduplications$DUPLICATES)
cor.test(Deduplications$POSITIONS, Deduplications$DUPLICATES)$estimate^2
cor.test(Deduplications$SEQUENCES, Deduplications$P)
cor.test(Deduplications$SEQUENCES, Deduplications$P)$estimate^2
cor.test(Deduplications$SEQUENCES, Deduplications$ALIGNMENTS)
cor.test(Deduplications$SEQUENCES, Deduplications$ALIGNMENTS)$estimate^2


Depths0=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Data/DepthFinder/gtSub_0.chrom-depth.txt", header=TRUE)
Depths10=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Data/DepthFinder/gtSub_10.chrom-depth.txt", header=TRUE)
Depths20=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Data/DepthFinder/gtSub_20.chrom-depth.txt", header=TRUE)
Depths30=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Data/DepthFinder/gtSub_30.chrom-depth.txt", header=TRUE)
Depths40=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Data/DepthFinder/gtSub_40.chrom-depth.txt", header=TRUE)
Depths50=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Data/DepthFinder/gtSub_50.chrom-depth.txt", header=TRUE)
Depths60=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Data/DepthFinder/gtSub_60.chrom-depth.txt", header=TRUE)
Depths70=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Data/DepthFinder/gtSub_70.chrom-depth.txt", header=TRUE)
Depths80=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Data/DepthFinder/gtSub_80.chrom-depth.txt", header=TRUE)
Depths90=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Data/DepthFinder/gtSub_90.chrom-depth.txt", header=TRUE)
Depths100=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Data/DepthFinder/gtSub_100.chrom-depth.txt", header=TRUE)


Depths0$P=0
Depths10$P=10
Depths20$P=20
Depths30$P=30
Depths40$P=40
Depths50$P=50
Depths60$P=60
Depths70$P=70
Depths80$P=80
Depths90$P=90
Depths100$P=100

Depths0$Chrom=c("1", "2", "3", "4", "4A", "1A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "17", "18", "19", "20", "21", "22", "23", "24",
"25LG1", "25LG2", "26", "27", "28", "LGE22", "Z")
Depths10$Chrom=c("1", "2", "3", "4", "4A", "1A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "17", "18", "19", "20", "21", "22", "23", "24",
"25LG1", "25LG2", "26", "27", "28", "LGE22", "Z")
Depths20$Chrom=c("1", "2", "3", "4", "4A", "1A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "17", "18", "19", "20", "21", "22", "23", "24",
"25LG1", "25LG2", "26", "27", "28", "LGE22", "Z")
Depths30$Chrom=c("1", "2", "3", "4", "4A", "1A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "17", "18", "19", "20", "21", "22", "23", "24",
"25LG1", "25LG2", "26", "27", "28", "LGE22", "Z")
Depths40$Chrom=c("1", "2", "3", "4", "4A", "1A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "17", "18", "19", "20", "21", "22", "23", "24",
"25LG1", "25LG2", "26", "27", "28", "LGE22", "Z")
Depths50$Chrom=c("1", "2", "3", "4", "4A", "1A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "17", "18", "19", "20", "21", "22", "23", "24",
"25LG1", "25LG2", "26", "27", "28", "LGE22", "Z")
Depths60$Chrom=c("1", "2", "3", "4", "4A", "1A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "17", "18", "19", "20", "21", "22", "23", "24",
"25LG1", "25LG2", "26", "27", "28", "LGE22", "Z")
Depths70$Chrom=c("1", "2", "3", "4", "4A", "1A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "17", "18", "19", "20", "21", "22", "23", "24",
"25LG1", "25LG2", "26", "27", "28", "LGE22", "Z")
Depths80$Chrom=c("1", "2", "3", "4", "4A", "1A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "17", "18", "19", "20", "21", "22", "23", "24",
"25LG1", "25LG2", "26", "27", "28", "LGE22", "Z")
Depths90$Chrom=c("1", "2", "3", "4", "4A", "1A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "17", "18", "19", "20", "21", "22", "23", "24",
"25LG1", "25LG2", "26", "27", "28", "LGE22", "Z")
Depths100$Chrom=c("1", "2", "3", "4", "4A", "1A", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "17", "18", "19", "20", "21", "22", "23", "24",
"25LG1", "25LG2", "26", "27", "28", "LGE22", "Z")

theDEPTHs=rbind(Depths0, Depths10, Depths20, Depths30, Depths40, Depths50, Depths60, Depths70, Depths80, Depths90, Depths100)
theDEPTHs$Chrom=factor(theDEPTHs$Chrom, levels=unique(theDEPTHs[order(theDEPTHs$Basepairs, decreasing=TRUE),]$Chrom))
theDEPTHs$D100P=rep(Depths100$Depth, 11)
theDEPTHs$RD100P=theDEPTHs$Depth/theDEPTHs$D100P

ggplot(theDEPTHs, aes(x=log10(Basepairs), y=Depth))+geom_point()+facet_wrap("P")+labs(x="log(Chromosome Size)", y="Depth (Reads)")
ggsave("depth-x-ChromSize.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Images", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)

ggplot(theDEPTHs[which(theDEPTHs$P==0 | theDEPTHs$P==100),], aes(x=log10(Basepairs), y=Depth))+geom_point()+facet_wrap("P")+labs(x="log(Chromosome Size)", y="Depth (Reads)")
ggsave("depth-x-ChromSizeExtremeP.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Images", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)

ggplot(theDEPTHs, aes(x=Chrom, y=Depth, groups=P))+geom_point(aes(color=as.factor(P)), size=5)+labs(x="Chromosome (Ordered by Size)", y="Depth (Reads)")+scale_color_manual(values=rev(brewer.pal("RdYlGn", n=11)))+guides(color=guide_legend(title="Blood %"))+theme(text=element_text(size=14), axis.text.x=element_text(size=4, face="bold"))
ggsave("depth-x-ChromAllP.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Images", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)

ggplot(theDEPTHs[which(theDEPTHs$Chrom=="25LG1"),], aes(x=as.factor(P), y=Depth, groups=P))+geom_point(aes(color=as.factor(P)), size=5)+labs(x="Chromosome", y="Depth (Reads)")+scale_color_manual(values=rev(brewer.pal("RdYlGn", n=11)))+guides(color=guide_legend(title="Blood %"))
ggsave("depth-x-NC_031793AllP.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Images", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)

ggplot(theDEPTHs, aes(x=Chrom, y=RD100P, groups=P))+geom_point(aes(color=as.factor(P)), size=5)+labs(x="Chromosome (Ordered by Size)", y="Depth (Reads)")+scale_color_manual(values=rev(brewer.pal("RdYlGn", n=11)))+guides(color=guide_legend(title="Blood %"))+theme(text=element_text(size=14), axis.text.x=element_text(size=4, face="bold"))
ggsave("RD100P-x-ChromAllP.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Depths-of-the-Soul/Images", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)


summary(lm(Depth~P+log10(Basepairs), data=theDEPTHs))
pRESID=resid(lm(Depth~P, data=theDEPTHs))
summary(lm(pRESID~log10(Basepairs), data=theDEPTHs))
