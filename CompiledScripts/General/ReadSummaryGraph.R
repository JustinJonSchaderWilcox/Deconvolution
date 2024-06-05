#Making Read Summary Graph
library(ggplot2)
library(RColorBrewer)

ReadEstimates=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ManuscriptMisc/Outputs/ReadCounts/Summary/DeconAbelSummaryReadCountSampled.txt", header=TRUE)
ReadNumbers=data.frame(P=rep(ReadEstimates$P, 4), Stage=c(rep("Subsampled", 11), rep("Trimmed", 11), rep("Aligned", 11), rep("Deduplicated", 11)), Reads=c(ReadEstimates$Subsampled, ReadEstimates$TrimmedPE, ReadEstimates$Aligned, ReadEstimates$Deduplicated))
ReadNumbers$Stage=factor(ReadNumbers$Stage, levels=c("Subsampled", "Trimmed", "Aligned", "Deduplicated"))
BloodGradient=rev(brewer.pal("Spectral", n=11))
ggplot(ReadNumbers, aes(x=Stage, y=Reads, fill=as.factor(P)))+geom_bar(position="stack", stat="identity")+theme(text=element_text(size=14))+labs(x="Stage of Processing", y="Paired-End Reads")+guides(fill=guide_legend(title="Blood %"))+scale_fill_manual(values=BloodGradient)
ggsave("peREAD-x-Stage.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ManuscriptMisc/Images/Reads", scale = 1, width = 5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)

ggplot(ReadEstimates, aes(x=PropBloodReads*100, y=BloodProportion*100))+geom_point(size=5, shape=18)+geom_abline(intercept=0, slope=1, linewidth=1.2, linetype="dotted")+labs(x="Blood Reads (%)", y="Blood Alignments (%)")+theme(text=element_text(size=14))
ggsave("corPbloodRead.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ManuscriptMisc/Images/Reads", scale = 1, width = 5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)

cor.test(ReadEstimates$PropBloodReads, ReadEstimates$BloodProportion)
(cor.test(ReadEstimates$PropBloodReads, ReadEstimates$BloodProportion)$estimate)^2
