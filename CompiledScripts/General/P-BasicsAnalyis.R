#Making Basic Statistic Validation Graphs
library(ggplot2)
library(ggnewscale)
library(RColorBrewer)

setwd("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/BasicStatistics")
SummarySNVsP=read.delim("deconSNV-and-P-Summary.txt", header=TRUE, sep="\t")
SomaticMarkers=data.frame(P=c(SummarySNVsP$P, SummarySNVsP$P), mSNVs=c(SummarySNVsP$SomaticSNVsMarkers-SummarySNVsP$SomaticSNVs, SummarySNVsP$SomaticSNVs), Type=c(rep("Reference", dim(SummarySNVsP)[1]), rep("Non-Ref", dim(SummarySNVsP)[1])))

allSNV_AFREQ=read.delim("decon-allSNV-AFREQ.txt", header=TRUE, sep="\t")
ultFiltSNV_AFREQ=read.delim("decon-ultFiltSNV-AFREQ.txt", header=TRUE, sep="\t")
mSomaticSNV_AFREQ=read.delim("decon-MarkersSomaticSNV-AFREQ.txt", header=TRUE, sep="\t")

allSNV_Depth=read.delim("decon-allSNV-Depth.txt", header=TRUE, sep="\t")
ultFiltSNV_Depth=read.delim("decon-ultFiltSNV-Depth.txt", header=TRUE, sep="\t")
mSomaticSNV_Depth=read.delim("decon-MarkersSomaticSNV-Depth.txt", header=TRUE, sep="\t")

allSNV_Phred=read.delim("decon-allSNV-Phred.txt", header=TRUE, sep="\t")
ultFiltSNV_Phred=read.delim("decon-ultFiltSNV-Phred.txt", header=TRUE, sep="\t")
mSomaticSNV_Phred=read.delim("decon-MarkersSomaticSNV-Phred.txt", header=TRUE, sep="\t")

aFreqMarkers=data.frame(P=c(allSNV_AFREQ$P, ultFiltSNV_AFREQ$P, mSomaticSNV_AFREQ$P), AFREQ=c(allSNV_AFREQ$AFREQ, ultFiltSNV_AFREQ$AFREQ, mSomaticSNV_AFREQ$AFREQ), SNV=c(rep("All", dim(allSNV_AFREQ)[1]), rep("Select", dim(ultFiltSNV_AFREQ)[1]), rep("Somatic", dim(mSomaticSNV_AFREQ)[1])))
DepthMarkers=data.frame(P=c(allSNV_Depth$P, ultFiltSNV_Depth$P, mSomaticSNV_Depth$P), Depth=c(allSNV_Depth$Depth, ultFiltSNV_Depth$Depth, mSomaticSNV_Depth$Depth), SNV=c(rep("All", dim(allSNV_Depth)[1]), rep("Select", dim(ultFiltSNV_Depth)[1]), rep("Somatic", dim(mSomaticSNV_Depth)[1])))
PhredMarkers=data.frame(P=c(allSNV_Phred$P, ultFiltSNV_Phred$P, mSomaticSNV_Phred$P), Phred=c(allSNV_Phred$Phred, ultFiltSNV_Phred$Phred, mSomaticSNV_Phred$Phred), SNV=c(rep("All", dim(allSNV_Phred)[1]), rep("Select", dim(ultFiltSNV_Phred)[1]), rep("Somatic", dim(mSomaticSNV_Phred)[1])))

allSNV_AFREQ$fAFREQ=allSNV_AFREQ$AFREQ
ultFiltSNV_AFREQ$fAFREQ=ultFiltSNV_AFREQ$AFREQ
mSomaticSNV_AFREQ$fAFREQ=mSomaticSNV_AFREQ$AFREQ

allSNV_AFREQ[allSNV_AFREQ$fAFREQ>0.5,]$fAFREQ=1-allSNV_AFREQ[allSNV_AFREQ$fAFREQ>0.5,]$fAFREQ
ultFiltSNV_AFREQ[ultFiltSNV_AFREQ$fAFREQ>0.5,]$fAFREQ=1-ultFiltSNV_AFREQ[ultFiltSNV_AFREQ$fAFREQ>0.5,]$fAFREQ
mSomaticSNV_AFREQ[mSomaticSNV_AFREQ$fAFREQ>0.5,]$fAFREQ=1-mSomaticSNV_AFREQ[mSomaticSNV_AFREQ$fAFREQ>0.5,]$fAFREQ

##SNVs and Markers
CountSNVs=data.frame(P=rep(SummarySNVsP$P, 3), SNVs=c(SummarySNVsP$SNVs, SummarySNVsP$ultFiltSNVs, SummarySNVsP$SomaticSNVsMarkers), Type=c(rep("All", 11), rep("Select", 11), rep("Somatic", 11)))
TypesSNV=c(brewer.pal(3, "Dark2")[c(1,2)], "#993333")
ggplot(CountSNVs, aes(x=as.factor(P), y=log10(SNVs), groups=Type))+geom_col(aes(fill=Type), position="dodge")+labs(x="Contamination: P (%)" , y=expression(Log[10]("SNVs")))+scale_fill_manual(values=TypesSNV)+scale_y_continuous(breaks=seq(0, 6, 1))+theme(text=element_text(size=18))
ggsave("markersCountSNVs.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/BasicStatistics", scale = 1, width = 7.5, height = 4, units = c("in"), dpi=300, limitsize = TRUE)

TypesSNV=c(brewer.pal(3, "Dark2")[c(1,2)], "#993333")
ggplot(aFreqMarkers[aFreqMarkers$AFREQ<=1,], aes(x=as.factor(P), y=AFREQ, groups=SNV))+geom_violin(aes(fill=SNV), draw_quantiles=c(0.25, 0.5, 0.75))+labs(x="Contamination: P (%)" , y="Allele Frequency")+scale_fill_manual(values=TypesSNV)+theme(text=element_text(size=10), legend.position="none")
ggsave("markersAfreqSNVs.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/BasicStatistics", scale = 1, width = 7.5, height = 2.5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(DepthMarkers, aes(x=as.factor(P), y=Depth, groups=SNV))+geom_violin(aes(fill=SNV), draw_quantiles=c(0.25, 0.5, 0.75))+scale_y_continuous(limits=quantile(DepthMarkers$Depth, c(0.01, 0.99)))+labs(x="Contamination: P (%)" , y="Depth (Reads)")+scale_fill_manual(values=TypesSNV)+theme(text=element_text(size=10), legend.position="none")
ggsave("markersDepthSNVs.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/BasicStatistics", scale = 1, width = 7.5, height = 2.5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(PhredMarkers, aes(x=as.factor(P), y=Phred, groups=SNV))+geom_violin(aes(fill=SNV), draw_quantiles=c(0.25, 0.5, 0.75))+labs(x="Contamination: P (%)" , y="Quality (Phred)")+scale_fill_manual(values=TypesSNV)+theme(text=element_text(size=10), legend.position="none")
ggsave("markersPhredSNVs.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/BasicStatistics", scale = 1, width = 7.5, height = 2.5, units = c("in"), dpi=300, limitsize = TRUE)

REFnonREF=brewer.pal(6, "Paired")[c(5,6)]
ggplot(SomaticMarkers, aes(x=as.factor(P), y=mSNVs, groups=Type))+geom_col(aes(fill=Type))+labs(x="Contamination: P (%)" , y="Somatic Marker SNVs")+scale_fill_manual(values=REFnonREF)+theme(text=element_text(size=18), legend.position="bottom")
ggsave("markersRefMatchSNVs.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/BasicStatistics", scale = 1, width = 7.5, height = 3, units = c("in"), dpi=300, limitsize = TRUE)

##Blood Somatic Metrics
bloodSomaticSNV_AFREQ=read.delim("decon-bloodSomaticSNV-AFREQ.txt", header=TRUE, sep="\t")
bloodSomaticSNV_Depth=read.delim("decon-bloodSomaticSNV-Depth.txt", header=TRUE, sep="\t")
bloodSomaticSNV_Phred=read.delim("decon-bloodSomaticSNV-Phred.txt", header=TRUE, sep="\t")

aFreqSomatic=data.frame(P=c(mSomaticSNV_AFREQ$P, bloodSomaticSNV_AFREQ$P), AFREQ=c(mSomaticSNV_AFREQ$AFREQ, bloodSomaticSNV_AFREQ$AFREQ), Tissue=c(rep("Blood", dim(mSomaticSNV_AFREQ)[1]), rep("Mix", dim(bloodSomaticSNV_AFREQ)[1])))
DepthSomatic=data.frame(P=c(mSomaticSNV_Depth$P, bloodSomaticSNV_Depth$P), Depth=c(mSomaticSNV_Depth$Depth, bloodSomaticSNV_Depth$Depth), Tissue=c(rep("Blood", dim(mSomaticSNV_Depth)[1]), rep("Mix", dim(bloodSomaticSNV_Depth)[1])))
PhredSomatic=data.frame(P=c(mSomaticSNV_Phred$P, bloodSomaticSNV_Phred$P), Phred=c(mSomaticSNV_Phred$Phred, bloodSomaticSNV_Phred$Phred), Tissue=c(rep("Blood", dim(mSomaticSNV_Phred)[1]), rep("Mix", dim(bloodSomaticSNV_Phred)[1])))

bloodSomaticSNV_AFREQ$fAFREQ=bloodSomaticSNV_AFREQ$AFREQ
bloodSomaticSNV_AFREQ[bloodSomaticSNV_AFREQ$fAFREQ>0.5,]$fAFREQ=1-bloodSomaticSNV_AFREQ[bloodSomaticSNV_AFREQ$fAFREQ>0.5,]$fAFREQ

cor.test(SummarySNVsP$SomaticSNVsMarkersAFREQ, SummarySNVsP$P)
cor.test(SummarySNVsP$SomaticSNVsAFREQ, SummarySNVsP$P)
summary(glm(SummarySNVsP$SomaticSNVsMarkersAFREQ~SummarySNVsP$P))
summary(glm(SummarySNVsP$SomaticSNVsAFREQ~SummarySNVsP$P))

cor.test(SummarySNVsP$P, SummarySNVsP$SomaticSNVsMarkersDEPTH)
cor.test(SummarySNVsP$P, SummarySNVsP$SomaticSNVsDEPTH)
summary(glm(SummarySNVsP$SomaticSNVsMarkersDEPTH~SummarySNVsP$P))
summary(glm(SummarySNVsP$SomaticSNVsDEPTH~SummarySNVsP$P))


ggplot(SummarySNVsP, aes(x=P))+geom_line(aes(y=SomaticSNVsMarkersAFREQ, color="Blood"), linewidth=3)+geom_line(aes(y=SomaticSNVsAFREQ, color="Mix"), linewidth=3)+labs(x="Contamination: P (%)" , y="Minor Allele Frequency")+scale_x_continuous(breaks=seq(0, 100, 10))+theme(text=element_text(size=18), legend.position="bottom")+scale_color_manual(name='Tissue', breaks=c("Blood", "Mix"), values=c("Blood"="#993333", "Mix"=brewer.pal(3, "Dark2")[3]))
ggsave("markersSummaryAfreqSNVs.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/BasicStatistics", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(SummarySNVsP, aes(x=P))+geom_line(aes(y=SomaticSNVsMarkersDEPTH, color="Blood"), linewidth=3)+geom_line(aes(y=SomaticSNVsDEPTH, color="Mix"), linewidth=3)+labs(x="Contamination: P (%)" , y="Allelic Depth")+scale_x_continuous(breaks=seq(0, 100, 10))+theme(text=element_text(size=18), legend.position="bottom")+scale_color_manual(name='Tissue', breaks=c("Blood", "Mix"), values=c("Blood"="#993333", "Mix"=brewer.pal(3, "Dark2")[3]))
ggsave("markersSummaryDepthSNVs.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/BasicStatistics", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)

pairwise.wilcox.test(mSomaticSNV_AFREQ$AFREQ, mSomaticSNV_AFREQ$P, p.adujst="holm")
pairwise.wilcox.test(mSomaticSNV_AFREQ$fAFREQ, mSomaticSNV_AFREQ$P, p.adujst="holm")
pairwise.wilcox.test(bloodSomaticSNV_AFREQ$AFREQ, bloodSomaticSNV_AFREQ$P, p.adjust="holm")
pairwise.wilcox.test(bloodSomaticSNV_AFREQ$fAFREQ, bloodSomaticSNV_AFREQ$P, p.adjust="holm")

PureMixSNV=c("#993333", brewer.pal(3, "Dark2")[3])
ggplot(aFreqSomatic, aes(x=as.factor(P), y=AFREQ, groups=Tissue))+geom_violin(aes(fill=Tissue), draw_quantiles=c(0.25, 0.5, 0.75))+labs(x="Contamination: P (%)" , y="Allele Frequency")+scale_fill_manual(values=PureMixSNV)+theme(text=element_text(size=18), legend.position="bottom")
ggsave("somaticAfreqSNVs.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/BasicStatistics", scale = 1, width = 7, height = 2.5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(DepthSomatic, aes(x=as.factor(P), y=Depth, groups=Tissue))+geom_violin(aes(fill=Tissue), draw_quantiles=c(0.25, 0.5, 0.75))+labs(x="Contamination: P (%)" , y="Allelic Depth")+scale_fill_manual(values=PureMixSNV)+theme(text=element_text(size=18), legend.position="bottom")
ggsave("somaticDepthSNVs.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/BasicStatistics", scale = 1, width = 7, height = 2.5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(PhredSomatic, aes(x=as.factor(P), y=Phred, groups=Tissue))+geom_violin(aes(fill=Tissue), draw_quantiles=c(0.25, 0.5, 0.75))+labs(x="Contamination: P (%)" , y="Quality (Phred)")+scale_fill_manual(values=PureMixSNV)+theme(text=element_text(size=18), legend.position="bottom")
ggsave("somaticPhredSNVs.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/BasicStatistics", scale = 1, width = 7, height = 2.5, units = c("in"), dpi=300, limitsize = TRUE)



BloodDeconP=read.delim("blood-deconvolution-validation-table.txt", header=TRUE, sep="\t")
ReadEstimates=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ManuscriptMisc/Outputs/ReadCounts/Summary/DeconAbelSummaryReadCountSampled.txt", header=TRUE)

ggplot(ReadEstimates, aes(x=PropBloodReads*100, y=BloodProportion*100))+geom_point(size=5, shape=18)+geom_abline(intercept=0, slope=1, linewidth=1.2, linetype="dotted")+labs(x="Blood Reads (%)", y="Blood Alignments (%)")
cor.test(ReadEstimates$PropBloodReads, ReadEstimates$BloodProportion)
cor.test(ReadEstimates$PropBloodReads, BloodDeconP$EstimatedP)
(cor.test(ReadEstimates$PropBloodReads, ReadEstimates$BloodProportion)$estimate)^2
(cor.test(ReadEstimates$PropBloodReads, BloodDeconP$EstimatedP)$estimate)^2


setwd("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation//BootStrap-P")
bootstrapBloodDeconP=read.delim("blood-deconvolution-bootstrap100-table.txt", header=TRUE, sep="\t")

BloodDeconP$Tissue="Blood"

corrBloodDeconP=cor.test(BloodDeconP$P,BloodDeconP$EstimatedP)
corrBloodDeconP
(corrBloodDeconP$estimate)^2

bootstrapBloodDeconP$BloodProportion=NA
bootstrapBloodDeconP[which(bootstrapBloodDeconP$P==0),]$BloodProportion=ReadEstimates[which(ReadEstimates$P==0),]$BloodProportion
bootstrapBloodDeconP[which(bootstrapBloodDeconP$P==10),]$BloodProportion=ReadEstimates[which(ReadEstimates$P==10),]$BloodProportion
bootstrapBloodDeconP[which(bootstrapBloodDeconP$P==20),]$BloodProportion=ReadEstimates[which(ReadEstimates$P==20),]$BloodProportion
bootstrapBloodDeconP[which(bootstrapBloodDeconP$P==30),]$BloodProportion=ReadEstimates[which(ReadEstimates$P==30),]$BloodProportion
bootstrapBloodDeconP[which(bootstrapBloodDeconP$P==40),]$BloodProportion=ReadEstimates[which(ReadEstimates$P==40),]$BloodProportion
bootstrapBloodDeconP[which(bootstrapBloodDeconP$P==50),]$BloodProportion=ReadEstimates[which(ReadEstimates$P==50),]$BloodProportion
bootstrapBloodDeconP[which(bootstrapBloodDeconP$P==60),]$BloodProportion=ReadEstimates[which(ReadEstimates$P==60),]$BloodProportion
bootstrapBloodDeconP[which(bootstrapBloodDeconP$P==70),]$BloodProportion=ReadEstimates[which(ReadEstimates$P==70),]$BloodProportion
bootstrapBloodDeconP[which(bootstrapBloodDeconP$P==80),]$BloodProportion=ReadEstimates[which(ReadEstimates$P==80),]$BloodProportion
bootstrapBloodDeconP[which(bootstrapBloodDeconP$P==90),]$BloodProportion=ReadEstimates[which(ReadEstimates$P==90),]$BloodProportion
bootstrapBloodDeconP[which(bootstrapBloodDeconP$P==100),]$BloodProportion=ReadEstimates[which(ReadEstimates$P==100),]$BloodProportion



colors=c(Blood="#FF0000",Brain="#CC33FF")
ggplot(data=BloodDeconP, aes(x=P, y=EstimatedP*100))+geom_line(linewidth=3, color="#FF0000")+scale_color_manual(values = colors)+labs(x="Contamination (%)" , y="Estimated (%)", color="Contaminant")+scale_x_continuous(breaks=seq(0, 100, 20))+scale_y_continuous(breaks=seq(0, 100, 20))+theme(text=element_text(size=14))
ggsave("decon-bloodEstP.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/ValidateP", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
Miscellanea=c("Estimate"="red", "Range"="Black", "Reads"="black", "Aligned"="darkgreen")
BOX=c("Est. Range"="firebrick")
ESTIMATE=c("Estimate"="firebrick")
READs=c("Reads"="black")
ALIGNED=c("Aligned"="darkgreen")
ggplot(data=bootstrapBloodDeconP, aes(x=as.factor(P), y=EstimatedP*100))+geom_boxplot(linewidth=1.4, aes(color="Est. Range"))+scale_color_manual(values=BOX)+guides(color=guide_legend(order=4))+labs(x="Contamination: P (%)" , y="Estimated P (%)")+new_scale_color()+scale_y_continuous(breaks=seq(0, 100, 10))+theme(text=element_text(size=14), legend.title=element_blank(), legend.position="bottom")+stat_summary(fun="mean", geom="line", aes(group=1, color="Estimate"), linewidth=01.8)+scale_color_manual(values=ESTIMATE)+guides(color=guide_legend(order=3))+new_scale_color()+geom_abline(aes(intercept=-10, slope=10, color="Reads"), linewidth=1.2, linetype="dotted")+scale_color_manual(values=READs)+guides(color=guide_legend(order=1))+new_scale_color()+geom_point(aes(y=BloodProportion*100, color="Aligned"), size=7, shape=8, stroke=1.4)+scale_color_manual(values=ALIGNED)+guides(color=guide_legend(order=2))+new_scale_color()
ggsave("decon-ValidateP.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/ValidateP", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
