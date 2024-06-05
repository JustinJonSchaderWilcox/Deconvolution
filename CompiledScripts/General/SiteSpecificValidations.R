#Methylated Site-Specific Validations
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(DescTools)
library(GenBinomApps)

MethCallDeconBrainMethSCoCpG=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/MergedDeconBrainMethSCoCpG.TRUE-P.txt", header=TRUE)
MethCallDeconBrainMethSCnCpG=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/MergedDeconBrainMethSCnCpG.TRUE-P.txt", header=TRUE)
MethCallDeconBrainMethWGoCpG=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/MergedDeconBrainMethWGoCpG.TRUE-P.txt", header=TRUE)
MethCallDeconBrainMethWGnCpG=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/MergedDeconBrainMethWGnCpG.TRUE-P.txt", header=TRUE)


#Depth Metrics
MethCallDeconBrainMethSCoCpG$MixDepth=MethCallDeconBrainMethSCoCpG$MixMethReads+MethCallDeconBrainMethSCoCpG$MixUnmethReads
MethCallDeconBrainMethSCnCpG$MixDepth=MethCallDeconBrainMethSCnCpG$MixMethReads+MethCallDeconBrainMethSCnCpG$MixUnmethReads
MethCallDeconBrainMethWGoCpG$MixDepth=MethCallDeconBrainMethWGoCpG$MixMethReads+MethCallDeconBrainMethWGoCpG$MixUnmethReads
MethCallDeconBrainMethWGnCpG$MixDepth=MethCallDeconBrainMethWGnCpG$MixMethReads+MethCallDeconBrainMethWGnCpG$MixUnmethReads


MethCallDeconBrainMethSCoCpG$DeconDepth=MethCallDeconBrainMethSCoCpG$DeconMethReads+MethCallDeconBrainMethSCoCpG$DeconUnmethReads
MethCallDeconBrainMethSCnCpG$DeconDepth=MethCallDeconBrainMethSCnCpG$DeconMethReads+MethCallDeconBrainMethSCnCpG$DeconUnmethReads
MethCallDeconBrainMethWGoCpG$DeconDepth=MethCallDeconBrainMethWGoCpG$DeconMethReads+MethCallDeconBrainMethWGoCpG$DeconUnmethReads
MethCallDeconBrainMethWGnCpG$DeconDepth=MethCallDeconBrainMethWGnCpG$DeconMethReads+MethCallDeconBrainMethWGnCpG$DeconUnmethReads

MethCallDeconBrainMethSCoCpG$tpDepth=MethCallDeconBrainMethSCoCpG$tpMethReads+MethCallDeconBrainMethSCoCpG$tpUnmethReads
MethCallDeconBrainMethSCnCpG$tpDepth=MethCallDeconBrainMethSCnCpG$tpMethReads+MethCallDeconBrainMethSCnCpG$tpUnmethReads
MethCallDeconBrainMethWGoCpG$tpDepth=MethCallDeconBrainMethWGoCpG$tpMethReads+MethCallDeconBrainMethWGoCpG$tpUnmethReads
MethCallDeconBrainMethWGnCpG$tpDepth=MethCallDeconBrainMethWGnCpG$tpMethReads+MethCallDeconBrainMethWGnCpG$tpUnmethReads

AllDepthSC=data.frame(P=c(rep(MethCallDeconBrainMethSCoCpG$P, 3), rep(MethCallDeconBrainMethSCnCpG$P, 3)), Context=c(rep("CpG", dim(MethCallDeconBrainMethSCoCpG)[1]*3), rep("Non-CpG", dim(MethCallDeconBrainMethSCnCpG)[1]*3)), Estimate=c(rep("Mix", dim(MethCallDeconBrainMethSCoCpG)[1]), rep("Est. P", dim(MethCallDeconBrainMethSCoCpG)[1]), rep("True P", dim(MethCallDeconBrainMethSCoCpG)[1]), rep("Mix", dim(MethCallDeconBrainMethSCnCpG)[1]), rep("Est. P", dim(MethCallDeconBrainMethSCnCpG)[1]), rep("True P", dim(MethCallDeconBrainMethSCnCpG)[1])), Depth=c(c(MethCallDeconBrainMethSCoCpG$MixDepth, MethCallDeconBrainMethSCoCpG$DeconDepth, MethCallDeconBrainMethSCoCpG$tpDepth), c(MethCallDeconBrainMethSCnCpG$MixDepth, MethCallDeconBrainMethSCnCpG$DeconDepth, MethCallDeconBrainMethSCnCpG$tpDepth)))
AllDepthSC$Estimate=factor(AllDepthSC$Estimate, levels=c("Mix", "True P", "Est. P"))
AllDepthWG=data.frame(P=c(rep(MethCallDeconBrainMethWGoCpG$P, 3), rep(MethCallDeconBrainMethWGnCpG$P, 3)), Context=c(rep("CpG", dim(MethCallDeconBrainMethWGoCpG)[1]*3), rep("Non-CpG", dim(MethCallDeconBrainMethWGnCpG)[1]*3)), Estimate=c(rep("Mix", dim(MethCallDeconBrainMethWGoCpG)[1]), rep("Est. P", dim(MethCallDeconBrainMethWGoCpG)[1]), rep("True P", dim(MethCallDeconBrainMethWGoCpG)[1]), rep("Mix", dim(MethCallDeconBrainMethWGnCpG)[1]), rep("Est. P", dim(MethCallDeconBrainMethWGnCpG)[1]), rep("True P", dim(MethCallDeconBrainMethWGnCpG)[1])), Depth=c(c(MethCallDeconBrainMethWGoCpG$MixDepth, MethCallDeconBrainMethWGoCpG$DeconDepth, MethCallDeconBrainMethWGoCpG$tpDepth), c(MethCallDeconBrainMethWGnCpG$MixDepth, MethCallDeconBrainMethWGnCpG$DeconDepth, MethCallDeconBrainMethWGnCpG$tpDepth)))
AllDepthWG$Estimate=factor(AllDepthWG$Estimate, levels=c("Mix", "True P", "Est. P"))

DeconEstimate=c(brewer.pal(3, "Dark2")[3], rev(brewer.pal(4, "Paired")))
ggplot(AllDepthSC, aes(x=as.factor(P), y=Depth, groups=Estimate))+geom_boxplot(aes(color=Estimate))+labs(x="Contamination Blood (%)", y="Depth (Reads)")+scale_color_manual(values=DeconEstimate)+theme(text=element_text(size=18))+facet_wrap("Context", nrow=2)
ggsave("AllDepthDeconSC-BoxPlot-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/Sites", scale = 1, width = 7.5, height = 3, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(AllDepthWG, aes(x=as.factor(P), y=Depth, groups=Estimate))+geom_boxplot(aes(color=Estimate))+labs(x="Contamination Blood (%)", y="Depth (Reads)")+scale_color_manual(values=DeconEstimate)+theme(text=element_text(size=18))+facet_wrap("Context", nrow=2)
ggsave("AllDepthDeconWG-BoxPlot-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/Sites", scale = 1, width = 7.5, height = 3, units = c("in"), dpi=300, limitsize = TRUE)

#Testing Binary Concordance of Deconvoluted P based on Significance
MethCallDeconBrainMethWGoCpG$Context="CpG"
MethCallDeconBrainMethWGnCpG$Context="Non-CpG"
MethCallDeconBrainMethWGmergedC=data.frame(rbind(MethCallDeconBrainMethWGoCpG, MethCallDeconBrainMethWGnCpG))

MethCallDeconBrainMethWGmergedC$DeconACC=1
MethCallDeconBrainMethWGmergedC$tpDeconACC=1
MethCallDeconBrainMethWGmergedC$MixACC=1

MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$DECONvBRAINdMP<0.01),]$DeconACC=0
MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$tpDECONvBRAINdMP<0.01),]$tpDeconACC=0
MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$MIXvBRAINdMP<0.01),]$MixACC=0

MethCallDeconBrainMethWGmergedC$DeconDirectionACC=0
MethCallDeconBrainMethWGmergedC$tpDeconDirectionACC=0
MethCallDeconBrainMethWGmergedC$MixDirectionACC=0

if(length(MethCallDeconBrainMethWGmergedC$DeconMethProportion>MethCallDeconBrainMethWGmergedC$BrainMethProportion)>=1){MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$DeconACC==0 & MethCallDeconBrainMethWGmergedC$DeconMethProportion>MethCallDeconBrainMethWGmergedC$BrainMethProportion),]$DeconDirectionACC=1}
if(length(MethCallDeconBrainMethWGmergedC$tpDeconMethProportion>MethCallDeconBrainMethWGmergedC$BrainMethProportion)>=1){MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$tpDeconACC==0 & MethCallDeconBrainMethWGmergedC$tpDeconMethProportion>MethCallDeconBrainMethWGmergedC$BrainMethProportion),]$tpDeconDirectionACC=1}
if(length(MethCallDeconBrainMethWGmergedC$MixMethProportion>MethCallDeconBrainMethWGmergedC$BrainMethProportion)>=1){MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$MixACC==0 & MethCallDeconBrainMethWGmergedC$MixMethProportion>MethCallDeconBrainMethWGmergedC$BrainMethProportion),]$MixDirectionACC=1}

if(length(MethCallDeconBrainMethWGmergedC$DeconMethProportion<MethCallDeconBrainMethWGmergedC$BrainMethProportion)>=1){MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$DeconACC==0 & MethCallDeconBrainMethWGmergedC$DeconMethProportion<MethCallDeconBrainMethWGmergedC$BrainMethProportion),]$DeconDirectionACC=-1}
if(length(MethCallDeconBrainMethWGmergedC$tpDeconMethProportion<MethCallDeconBrainMethWGmergedC$BrainMethProportion)>=1){MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$tpDeconACC==0 & MethCallDeconBrainMethWGmergedC$tpDeconMethProportion<MethCallDeconBrainMethWGmergedC$BrainMethProportion),]$tpDeconDirectionACC=-1}
if(length(MethCallDeconBrainMethWGmergedC$MixMethProportion<MethCallDeconBrainMethWGmergedC$BrainMethProportion)>=1){MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$MixACC==0 & MethCallDeconBrainMethWGmergedC$MixMethProportion<MethCallDeconBrainMethWGmergedC$BrainMethProportion),]$MixDirectionACC=-1}

AllCallACC=setDT(MethCallDeconBrainMethWGmergedC)[,.(DeconDepth=mean(DeconDepth), tpDepth=mean(tpDepth), MixDepth=mean(MixDepth), DeconACC=mean(DeconACC), tpDeconACC=mean(tpDeconACC), MixACC=mean(MixACC), DeconDirectionACC=mean(DeconDirectionACC), tpDeconDirectionACC=mean(tpDeconDirectionACC), MixDirectionACC=mean(MixDirectionACC)), by=c("P", "Context")]
AllCallACC$DeconConfLowACC=NA
AllCallACC$tpDeconConfLowACC=NA
AllCallACC$MixConfLowACC=NA
AllCallACC$DeconConfUpACC=NA
AllCallACC$tpDeconConfUpACC=NA
AllCallACC$MixConfUpACC=NA
for (i in 1:dim(AllCallACC)[1]){
AllCallACC$DeconConfLowACC[i]=clopper.pearson.ci(AllCallACC$DeconACC[i]*30000, 30000, alpha=0.01, CI="two.sided")$Lower.limit
AllCallACC$tpDeconConfLowACC[i]=clopper.pearson.ci(AllCallACC$tpDeconACC[i]*30000, 30000, alpha=0.01, CI="two.sided")$Lower.limit
AllCallACC$MixConfLowACC[i]=clopper.pearson.ci(AllCallACC$MixACC[i]*30000, 30000, alpha=0.01, CI="two.sided")$Lower.limit

AllCallACC$DeconConfUpACC[i]=clopper.pearson.ci(AllCallACC$DeconACC[i]*30000, 30000, alpha=0.01, CI="two.sided")$Upper.limit
AllCallACC$tpDeconConfUpACC[i]=clopper.pearson.ci(AllCallACC$tpDeconACC[i]*30000, 30000, alpha=0.01, CI="two.sided")$Upper.limit
AllCallACC$MixConfUpACC[i]=clopper.pearson.ci(AllCallACC$MixACC[i]*30000, 30000, alpha=0.01, CI="two.sided")$Upper.limit
}

AllCallACC$DeconBrainCCCP=NA
AllCallACC$DeconBrainLowerCCCP=NA
AllCallACC$DeconBrainUpperCCCP=NA
AllCallACC$DeconBrainCCCP.CB=NA
AllCallACC$tpDeconBrainCCCP=NA
AllCallACC$tpDeconBrainLowerCCCP=NA
AllCallACC$tpDeconBrainUpperCCCP=NA
AllCallACC$tpDeconBrainCCCP.CB=NA
AllCallACC$MixBrainCCCP=NA
AllCallACC$MixBrainLowerCCCP=NA
AllCallACC$MixBrainUpperCCCP=NA
AllCallACC$MixBrainCCCP.CB=NA
AllCallACC$WilcoxDEVp[i]=NA
AllCallACC$WilcoxDEVtV=NA
AllCallACC$DeconDevMean=NA
AllCallACC$MixDevMean=NA
AllCallACC$DECONvMIXbsCCCp=NA
AllCallACC$tpDECONvMIXbsCCCp=NA
AllCallACC$DECONvtpDECONbsCCCp[i]=NA
#Concordance Correlations
for (i in 1:dim(AllCallACC)[1]){
p=AllCallACC$P[i]
c=AllCallACC$Context[i]
AllCallACC$DeconBrainCCCP[i]=as.numeric(CCC(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$DeconMethProportion, MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$BrainMethProportion, ci="z-transform")$rho.c[1])
AllCallACC$DeconBrainLowerCCCP[i]=as.numeric(CCC(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$DeconMethProportion, MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$BrainMethProportion, ci="z-transform")$rho.c[2])
AllCallACC$DeconBrainUpperCCCP[i]=as.numeric(CCC(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$DeconMethProportion, MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$BrainMethProportion, ci="z-transform")$rho.c[3])
AllCallACC$DeconBrainCCCP.CB[i]=as.numeric(CCC(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$DeconMethProportion, MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$BrainMethProportion, ci="z-transform")$C.b)

AllCallACC$tpDeconBrainCCCP[i]=as.numeric(CCC(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$tpMethProportion, MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$BrainMethProportion, ci="z-transform")$rho.c[1])
AllCallACC$tpDeconBrainLowerCCCP[i]=as.numeric(CCC(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$tpMethProportion, MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$BrainMethProportion, ci="z-transform")$rho.c[2])
AllCallACC$tpDeconBrainUpperCCCP[i]=as.numeric(CCC(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$tpMethProportion, MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$BrainMethProportion, ci="z-transform")$rho.c[3])
AllCallACC$tpDeconBrainCCCP.CB[i]=as.numeric(CCC(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$tpMethProportion, MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$BrainMethProportion, ci="z-transform")$C.b)

AllCallACC$MixBrainCCCP[i]=as.numeric(CCC(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$MixMethProportion, MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$BrainMethProportion, ci="z-transform")$rho.c[1])
AllCallACC$MixBrainLowerCCCP[i]=as.numeric(CCC(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$MixMethProportion, MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$BrainMethProportion, ci="z-transform")$rho.c[2])
AllCallACC$MixBrainUpperCCCP[i]=as.numeric(CCC(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$MixMethProportion, MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$BrainMethProportion, ci="z-transform")$rho.c[3])
AllCallACC$MixBrainCCCP.CB[i]=as.numeric(CCC(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$MixMethProportion, MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$BrainMethProportion, ci="z-transform")$C.b)
WILCOX=wilcox.test(unlist(CCC(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$DeconMethProportion, MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$BrainMethProportion, ci="z-transform")$blalt[2]), unlist(CCC(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$MixMethProportion, MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$BrainMethProportion, ci="z-transform")$blalt[2]), paired=TRUE)
AllCallACC$WilcoxDEVp[i]=WILCOX$p.value
AllCallACC$WilcoxDEVtV[i]=WILCOX$statistic
AllCallACC$DeconDevMean[i]=mean(unlist(CCC(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$DeconMethProportion, MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$BrainMethProportion, ci="z-transform")$blalt[2]))
AllCallACC$MixDevMean[i]=mean(unlist(CCC(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$MixMethProportion, MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),]$BrainMethProportion, ci="z-transform")$blalt[2]))

DECON=0
tpDECON=0
DECONvTP=0
n=dim(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c),])[1]
if(p!=100){
for (b in 1:10000){
SITES=sort(sample(which(MethCallDeconBrainMethWGmergedC$P==p & MethCallDeconBrainMethWGmergedC$Context==c), size=n, replace=TRUE))
BOOTSTRAP=MethCallDeconBrainMethWGmergedC[SITES,]
bsDeconCCC=as.numeric(CCC(BOOTSTRAP$DeconMethProportion, BOOTSTRAP$BrainMethProportion, ci="z-transform")$rho.c[1])
BStpDeconCCC=as.numeric(CCC(BOOTSTRAP$tpMethProportion, BOOTSTRAP$BrainMethProportion, ci="z-transform")$rho.c[1])
bsMixCCC=as.numeric(CCC(BOOTSTRAP$MixMethProportion, BOOTSTRAP$BrainMethProportion, ci="z-transform")$rho.c[1])
if(bsDeconCCC>bsMixCCC){
DECON=DECON+1
}
if(BStpDeconCCC>bsMixCCC){
tpDECON=tpDECON+1
}
if(BStpDeconCCC>bsDeconCCC){
DECONvTP=DECONvTP+1
}
}
AllCallACC$DECONvMIXbsCCCp[i]=(1-DECON/b)
AllCallACC$tpDECONvMIXbsCCCp[i]=(1-tpDECON/b)
AllCallACC$DECONvtpDECONbsCCCp[i]=(1-DECONvTP/b)
}
}

AllCallACC$VerdictBSp="NA"
AllCallACC$VerdictP="NA"
AllCallACC[which(AllCallACC$DECONvMIXbsCCCp<0.0001),]$VerdictBSp="Decon"
AllCallACC[which(AllCallACC$DECONvMIXbsCCCp>0.0001),]$VerdictBSp="Mix"
AllCallACC[which(AllCallACC$DECONvtpDECONbsCCCp<0.0001),]$VerdictP="P"
AllCallACC[which(AllCallACC$DECONvtpDECONbsCCCp>0.0001),]$VerdictP="TrueP"

AllCallACC[which(AllCallACC$Context=="CpG" & AllCallACC$VerdictBSp=="Decon"),]$P
AllCallACC[which(AllCallACC$Context=="Non-CpG" & AllCallACC$VerdictBSp=="Decon"),]$P
AllCallACC[which(AllCallACC$Context=="CpG" & AllCallACC$VerdictP=="P"),]$P
AllCallACC[which(AllCallACC$Context=="Non-CpG" & AllCallACC$VerdictP=="P"),]$P

AllCallACC$Verdict="NA"
AllCallACC[which(abs(AllCallACC$DeconDevMean)<abs(AllCallACC$MixDevMean)),]$Verdict="Decon"
AllCallACC[which(abs(AllCallACC$DeconDevMean)>abs(AllCallACC$MixDevMean)),]$Verdict="Mix"

AllCallACC[which(AllCallACC$Context=="CpG" & AllCallACC$Verdict=="Decon"),]$P
AllCallACC[which(AllCallACC$Context=="Non-CpG" & AllCallACC$Verdict=="Decon"),]$P
max(AllCallACC[which(AllCallACC$P!=0 & AllCallACC$P!=100),]$WilcoxDEVp)
min(AllCallACC[which(AllCallACC$P!=0 & AllCallACC$P!=100),]$WilcoxDEVtV)

AllCallACC.long=data.frame(P=rep(AllCallACC$P, 3), Context=rep(AllCallACC$Context, 3), Depth=c(AllCallACC$DeconDepth, AllCallACC$tpDepth, AllCallACC$MixDepth), Estimate=c(rep("Est. P", 22), rep("True P", 22), rep("Mix", 22)), Accuracy=c(AllCallACC$DeconACC, AllCallACC$tpDeconACC, AllCallACC$MixACC), LowerConf=c(AllCallACC$DeconConfLowACC, AllCallACC$tpDeconConfLowACC, AllCallACC$MixConfLowACC), UpperConf=c(AllCallACC$DeconConfUpACC, AllCallACC$tpDeconConfUpACC, AllCallACC$MixConfUpACC), CCCP=as.numeric(c(AllCallACC$DeconBrainCCCP, AllCallACC$tpDeconBrainCCCP, AllCallACC$MixBrainCCCP)), CB=c(AllCallACC$DeconBrainCCCP.CB, AllCallACC$tpDeconBrainCCCP.CB, AllCallACC$MixBrainCCCP.CB))

MethCallDeconBrainMethWGmergedC.long=data.frame(P=rep(MethCallDeconBrainMethWGmergedC$P,3), Estimate=c(rep("Est. P", dim(MethCallDeconBrainMethWGmergedC)[1]), rep("True P", dim(MethCallDeconBrainMethWGmergedC)[1]), rep("Mix", dim(MethCallDeconBrainMethWGmergedC)[1])), Context=rep(MethCallDeconBrainMethWGmergedC$Context, 3), Position=rep(paste(MethCallDeconBrainMethWGmergedC$Chromosome, MethCallDeconBrainMethWGmergedC$Position, sep=":"), 3), Depth=c(MethCallDeconBrainMethWGmergedC$DeconDepth, MethCallDeconBrainMethWGmergedC$tpDepth, MethCallDeconBrainMethWGmergedC$MixDepth), MethProp=c(MethCallDeconBrainMethWGmergedC$DeconMethProportion, MethCallDeconBrainMethWGmergedC$tpMethProportion, MethCallDeconBrainMethWGmergedC$MixMethProportion), BrainMethProp=rep(MethCallDeconBrainMethWGmergedC$BrainMethProportion, 3))
MethCallDeconBrainMethWGmergedC.long$Deviation=abs(MethCallDeconBrainMethWGmergedC.long$MethProp-MethCallDeconBrainMethWGmergedC.long$BrainMethProp)
MethCallDeconBrainMethWGmergedC.long$DevDepthResid=resid(aov(Deviation~Estimate/P+Depth,data=MethCallDeconBrainMethWGmergedC.long))
t.test(MethCallDeconBrainMethWGmergedC.long[which(MethCallDeconBrainMethWGmergedC.long$Context=="CpG" & MethCallDeconBrainMethWGmergedC.long$Estimate=="Est. P"),]$DevDepthResid, MethCallDeconBrainMethWGmergedC.long[which(MethCallDeconBrainMethWGmergedC.long$Context=="CpG" & MethCallDeconBrainMethWGmergedC.long$Estimate=="Mix"),]$DevDepthResid, pooled=TRUE)
mean(MethCallDeconBrainMethWGmergedC.long[which(MethCallDeconBrainMethWGmergedC.long$Context=="CpG" & MethCallDeconBrainMethWGmergedC.long$Estimate=="Est. P"),]$DevDepthResid)
mean(MethCallDeconBrainMethWGmergedC.long[which(MethCallDeconBrainMethWGmergedC.long$Context=="CpG" & MethCallDeconBrainMethWGmergedC.long$Estimate=="Mix"),]$DevDepthResid)

MethCallDeconBrainMethWGmergedC.long=data.frame(P=rep(MethCallDeconBrainMethWGmergedC$P,3), Estimate=c(rep("Est. P", dim(MethCallDeconBrainMethWGmergedC)[1]), rep("True P", dim(MethCallDeconBrainMethWGmergedC)[1]), rep("Mix", dim(MethCallDeconBrainMethWGmergedC)[1])), Context=rep(MethCallDeconBrainMethWGmergedC$Context, 3), Position=rep(paste(MethCallDeconBrainMethWGmergedC$Chromosome, MethCallDeconBrainMethWGmergedC$Position, sep=":"), 3), Depth=c(MethCallDeconBrainMethWGmergedC$DeconDepth, MethCallDeconBrainMethWGmergedC$tpDepth, MethCallDeconBrainMethWGmergedC$MixDepth), MethProp=c(MethCallDeconBrainMethWGmergedC$DeconMethProportion, MethCallDeconBrainMethWGmergedC$tpMethProportion, MethCallDeconBrainMethWGmergedC$MixMethProportion), BrainMethProp=rep(MethCallDeconBrainMethWGmergedC$BrainMethProportion, 3))
MethCallDeconBrainMethWGmergedC.long$Deviation=abs(MethCallDeconBrainMethWGmergedC.long$MethProp-MethCallDeconBrainMethWGmergedC.long$BrainMethProp)
MethCallDeconBrainMethWGmergedC.long$DevDepthResid=resid(aov(Deviation~Estimate/P+Depth,data=MethCallDeconBrainMethWGmergedC.long))
t.test(MethCallDeconBrainMethWGmergedC.long[which(MethCallDeconBrainMethWGmergedC.long$Context=="Non-CpG" & MethCallDeconBrainMethWGmergedC.long$Estimate=="Est. P"),]$DevDepthResid, MethCallDeconBrainMethWGmergedC.long[which(MethCallDeconBrainMethWGmergedC.long$Context=="Non-CpG" & MethCallDeconBrainMethWGmergedC.long$Estimate=="Mix"),]$DevDepthResid, pooled=TRUE)
mean(MethCallDeconBrainMethWGmergedC.long[which(MethCallDeconBrainMethWGmergedC.long$Context=="Non-CpG" & MethCallDeconBrainMethWGmergedC.long$Estimate=="Est. P"),]$DevDepthResid)
mean(MethCallDeconBrainMethWGmergedC.long[which(MethCallDeconBrainMethWGmergedC.long$Context=="Non-CpG" & MethCallDeconBrainMethWGmergedC.long$Estimate=="Mix"),]$DevDepthResid)


t.test(abs(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context=="CpG"),]$DeconMethProportion-MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context=="CpG"),]$BrainMethProportion), abs(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context=="CpG"),]$MixMethProportion-MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context=="CpG"),]$BrainMethProportion), paired=TRUE) 
mean(abs(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context=="CpG"),]$DeconMethProportion-MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context=="CpG"),]$BrainMethProportion))
mean(abs(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context=="CpG"),]$MixMethProportion-MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context=="CpG"),]$BrainMethProportion))


t.test(abs(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context=="Non-CpG"),]$DeconMethProportion-MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context=="Non-CpG"),]$BrainMethProportion), abs(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context=="Non-CpG"),]$MixMethProportion-MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context=="Non-CpG"),]$BrainMethProportion), paired=TRUE) 
mean(abs(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context=="Non-CpG"),]$DeconMethProportion-MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context=="Non-CpG"),]$BrainMethProportion))
mean(abs(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context=="Non-CpG"),]$MixMethProportion-MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context=="Non-CpG"),]$BrainMethProportion))


t.test(AllCallACC[which(AllCallACC$Context=="CpG"),]$DeconBrainCCCP, AllCallACC[which(AllCallACC$Context=="CpG"),]$MixBrainCCCP, paired=TRUE)
mean(na.omit(AllCallACC[which(AllCallACC$Context=="CpG"),]$DeconBrainCCCP))
mean(na.omit(AllCallACC[which(AllCallACC$Context=="CpG"),]$MixBrainCCCP))

t.test(AllCallACC[which(AllCallACC$Context=="Non-CpG"),]$DeconBrainCCCP, AllCallACC[which(AllCallACC$Context=="Non-CpG"),]$MixBrainCCCP, paired=TRUE)
mean(na.omit(AllCallACC[which(AllCallACC$Context=="Non-CpG"),]$DeconBrainCCCP))
mean(na.omit(AllCallACC[which(AllCallACC$Context=="Non-CpG"),]$MixBrainCCCP))

CCC(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context=="CpG"),]$DeconMethProportion, MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context==c),]$BrainMethProportion, ci="z-transform")$rho.c
CCC(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context=="CpG"),]$MixMethProportion, MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context==c),]$BrainMethProportion, ci="z-transform")$rho.c

CCC(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context=="Non-CpG"),]$DeconMethProportion, MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context==c),]$BrainMethProportion, ci="z-transform")$rho.c
CCC(MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context=="Non-CpG"),]$MixMethProportion, MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$Context==c),]$BrainMethProportion, ci="z-transform")$rho.c


MethCallDeconBrainMethWGoCpG$Context="CpG"
MethCallDeconBrainMethWGnCpG$Context="Non-CpG"
MethCallDeconBrainMethWGmergedC=data.frame(rbind(MethCallDeconBrainMethWGoCpG, MethCallDeconBrainMethWGnCpG))

MethCallDeconBrainMethWGmergedC$TrueSCC=0
MethCallDeconBrainMethWGmergedC$DeconSCC=0
MethCallDeconBrainMethWGmergedC$tpDeconSCC=0
MethCallDeconBrainMethWGmergedC$MixSCC=0

MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$TrueDiffMethP<0.01),]$TrueSCC=1
MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$PredictDiffMethP<0.01 & MethCallDeconBrainMethWGmergedC$TrueDiffMethP<0.01 & ((MethCallDeconBrainMethWGmergedC$DeconMethProportion>MethCallDeconBrainMethWGmergedC$BloodMethProportion & MethCallDeconBrainMethWGmergedC$BrainMethProportion>MethCallDeconBrainMethWGmergedC$BloodMethProportion) | (MethCallDeconBrainMethWGmergedC$DeconMethProportion<MethCallDeconBrainMethWGmergedC$BloodMethProportion & MethCallDeconBrainMethWGmergedC$BrainMethProportion<MethCallDeconBrainMethWGmergedC$BloodMethProportion))),]$DeconSCC=1
MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$tpPredictDiffMethP<0.01 & MethCallDeconBrainMethWGmergedC$TrueDiffMethP<0.01 & ((MethCallDeconBrainMethWGmergedC$tpMethProportion>MethCallDeconBrainMethWGmergedC$BloodMethProportion & MethCallDeconBrainMethWGmergedC$BrainMethProportion>MethCallDeconBrainMethWGmergedC$BloodMethProportion) | (MethCallDeconBrainMethWGmergedC$tpMethProportion<MethCallDeconBrainMethWGmergedC$BloodMethProportion & MethCallDeconBrainMethWGmergedC$BrainMethProportion<MethCallDeconBrainMethWGmergedC$BloodMethProportion))),]$tpDeconSCC=1
MethCallDeconBrainMethWGmergedC[which(MethCallDeconBrainMethWGmergedC$MixDiffMeth<0.01 & MethCallDeconBrainMethWGmergedC$TrueDiffMethP<0.01 & ((MethCallDeconBrainMethWGmergedC$MixMethProportion>MethCallDeconBrainMethWGmergedC$BloodMethProportion & MethCallDeconBrainMethWGmergedC$BrainMethProportion>MethCallDeconBrainMethWGmergedC$BloodMethProportion) | (MethCallDeconBrainMethWGmergedC$MixMethProportion<MethCallDeconBrainMethWGmergedC$BloodMethProportion & MethCallDeconBrainMethWGmergedC$BrainMethProportion<MethCallDeconBrainMethWGmergedC$BloodMethProportion))),]$MixSCC=1

AllCallACC.long=data.frame(P=rep(AllCallACC$P, 3), Context=rep(AllCallACC$Context, 3), Depth=c(AllCallACC$DeconDepth, AllCallACC$tpDepth, AllCallACC$MixDepth), Estimate=c(rep("Est. P", 22), rep("True P", 22), rep("Mix", 22)), Accuracy=c(AllCallACC$DeconACC, AllCallACC$tpDeconACC, AllCallACC$MixACC), LowerConf=c(AllCallACC$DeconConfLowACC, AllCallACC$tpDeconConfLowACC, AllCallACC$MixConfLowACC), UpperConf=c(AllCallACC$DeconConfUpACC, AllCallACC$tpDeconConfUpACC, AllCallACC$MixConfUpACC), CCCP=as.numeric(c(AllCallACC$DeconBrainCCCP, AllCallACC$tpDeconBrainCCCP, AllCallACC$MixBrainCCCP)), CB=c(AllCallACC$DeconBrainCCCP.CB, AllCallACC$tpDeconBrainCCCP.CB, AllCallACC$MixBrainCCCP.CB))

AllCallACCcorDEPTHxRHO=data.frame(Context=rep(c("CpG", "Non-CpG"), 3), Estimate=rep(c("Est. P", "True P", "Mix"), 2))
AllCallACCcorDEPTHxRHO$CorrCoeficient=NA
for (i in 1:dim(AllCallACCcorDEPTHxRHO)[1]){
AllCallACCcorDEPTHxRHO$CorrCoeficient[i]=cor.test(as.numeric(na.omit(AllCallACC.long[which(AllCallACC.long$Context==AllCallACCcorDEPTHxRHO$Context[i] & AllCallACC.long$Estimate==AllCallACCcorDEPTHxRHO$Estimate[i]),]$CCCP)), as.numeric(na.omit(AllCallACC.long[which(AllCallACC.long$Context==AllCallACCcorDEPTHxRHO$Context[i] & AllCallACC.long$Estimate==AllCallACCcorDEPTHxRHO$Estimate[i] & AllCallACC.long$Depth!=0),]$Depth)))$estimate
}


AllCallSCC=setDT(MethCallDeconBrainMethWGmergedC)[,.(DeconSCC=mean(DeconSCC)/mean(TrueSCC), tpDeconSCC=mean(tpDeconSCC)/mean(TrueSCC), MixSCC=mean(MixSCC)/mean(TrueSCC)), , by=c("P", "Context")]

AllCallSCC$DeconConfLowSCC=NA
AllCallSCC$tpDeconConfLowSCC=NA
AllCallSCC$MixConfLowSCC=NA
AllCallSCC$DeconConfUpSCC=NA
AllCallSCC$tpDeconConfUpSCC=NA
AllCallSCC$MixConfUpSCC=NA
for (i in 1:dim(AllCallSCC)[1]){
AllCallSCC$DeconConfLowSCC[i]=clopper.pearson.ci(round(AllCallSCC$DeconSCC[i]*30000), 30000, alpha=0.01, CI="two.sided")$Lower.limit
AllCallSCC$tpDeconConfLowSCC[i]=clopper.pearson.ci(round(AllCallSCC$tpDeconSCC[i]*30000), 30000, alpha=0.01, CI="two.sided")$Lower.limit
AllCallSCC$MixConfLowSCC[i]=clopper.pearson.ci(round(AllCallSCC$MixSCC[i]*30000), 30000, alpha=0.01, CI="two.sided")$Lower.limit

AllCallSCC$DeconConfUpSCC[i]=clopper.pearson.ci(round(AllCallSCC$DeconSCC[i]*30000), 30000, alpha=0.01, CI="two.sided")$Upper.limit
AllCallSCC$tpDeconConfUpSCC[i]=clopper.pearson.ci(round(AllCallSCC$tpDeconSCC[i]*30000), 30000, alpha=0.01, CI="two.sided")$Upper.limit
AllCallSCC$MixConfUpSCC[i]=clopper.pearson.ci(round(AllCallSCC$MixSCC[i]*30000), 30000, alpha=0.01, CI="two.sided")$Upper.limit
}

AllCallSCC.long=data.frame(P=rep(AllCallSCC$P, 3), Context=rep(AllCallSCC$Context, 3), Estimate=c(rep("Est. P", 22), rep("True P", 22), rep("Mix", 22)), Accuracy=c(AllCallSCC$DeconSCC, AllCallSCC$tpDeconSCC, AllCallSCC$MixSCC), LowerConf=c(AllCallSCC$DeconConfLowSCC, AllCallSCC$tpDeconConfLowSCC, AllCallSCC$MixConfLowSCC), UpperConf=c(AllCallSCC$DeconConfUpSCC, AllCallSCC$tpDeconConfUpSCC, AllCallSCC$MixConfUpSCC))

AllCallACC.long$Estimate=factor(AllCallACC.long$Estimate, levels=c("Mix", "True P", "Est. P"))
AllCallSCC.long$Estimate=factor(AllCallSCC.long$Estimate, levels=c("Mix", "True P", "Est. P"))
AllCallACCcorDEPTHxRHO$Estimate=factor(AllCallACCcorDEPTHxRHO$Estimate, levels=c("Mix", "True P", "Est. P"))
DeconEstimate=c(brewer.pal(3, "Dark2")[3], rev(brewer.pal(4, "Paired")))

ggplot(AllCallACC.long, aes(x=P, y=Estimate, fill=CCCP))+geom_tile()+facet_wrap(nrow=1, "Context")+scale_fill_gradient2(low = "red", high = "green", mid = "yellow", midpoint = 0.5, limit = c(0,1), space = "Lab", name=expression(ρ[c]))+labs(x="Contamination Blood: P (%)", y="")+theme(text=element_text(size=18))
ggsave("allcCCCP-Tile-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/Sites", scale = 1, width = 7.5, height = 3, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(AllCallACC.long[which(AllCallACC.long$Context=="CpG"),], aes(x=P, y=Estimate, fill=CCCP))+geom_tile()+scale_fill_gradient2(low = "red", high = "green", mid = "yellow", midpoint = 0.5, limit = c(0,1), space = "Lab", name=expression(ρ[c]))+facet_wrap(nrow=1, "Context")+labs(x="Contamination Blood: P (%)", y="")+theme(text=element_text(size=18))
ggsave("cpgCCCP-Tile-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/Sites", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(AllCallACC.long[which(AllCallACC.long$Context=="Non-CpG"),], aes(x=P, y=Estimate, fill=CCCP))+geom_tile()+scale_fill_gradient2(low = "red", high = "green", mid = "yellow", midpoint = 0.5, limit = c(0,1), space = "Lab", name=expression(ρ[c]))+facet_wrap(nrow=1, "Context")+labs(x="Contamination Blood: P (%)", y="")+theme(text=element_text(size=18))
ggsave("ncpgCCCP-Tile-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/Sites", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(AllCallACC.long[which(AllCallACC.long$Context=="CpG"),], aes(x=P, y=Estimate, fill=CCCP))+geom_tile()+scale_fill_gradient2(low = "red", high = "green", mid = "yellow", midpoint = (min(na.omit(AllCallACC.long[which(AllCallACC.long$Context=="CpG"),]$CCCP))+1)/2, limit = c(min(na.omit(AllCallACC.long[which(AllCallACC.long$Context=="CpG"),]$CCCP)),1), space = "Lab", name="Rho")+facet_wrap(nrow=1, "Context")+labs(x="Contamination Blood: P (%)", y="")+theme(text=element_text(size=18))
ggsave("cpgCCCPv-Tile-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/Sites", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(AllCallACC.long[which(AllCallACC.long$Context=="Non-CpG"),], aes(x=P, y=Estimate, fill=CCCP))+geom_tile()+scale_fill_gradient2(low = "red", high = "green", mid = "yellow", midpoint = (min(na.omit(AllCallACC.long[which(AllCallACC.long$Context=="Non-CpG"),]$CCCP))+1)/2, limit = c(min(na.omit(AllCallACC.long[which(AllCallACC.long$Context=="Non-CpG"),]$CCCP)),1), space = "Lab", name="Rho")+facet_wrap(nrow=1, "Context")+labs(x="Contamination Blood: P (%)", y="")+theme(text=element_text(size=18))
ggsave("ncpgCCCPv-Tile-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/Sites", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)


ggplot(AllCallACC.long, aes(x=as.factor(P), y=Estimate, fill=Depth))+geom_tile()+scale_fill_gradient2(low="red", high="green", mid="yellow", midpoint = max(AllCallACC.long$Depth)/2, limit=c(0, max(AllCallACC.long$Depth)), space="Lab", name="Depth")+labs(x="Contamination Blood: P (%)", y="")+scale_color_manual(values=DeconEstimate)+theme(text=element_text(size=18))+facet_wrap(nrow=1, "Context")
ggsave("allcDEPTH-Tile-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/Sites", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(AllCallACC.long[which(AllCallACC.long$Context=="CpG"),], aes(x=as.factor(P), y=Estimate, fill=Depth))+geom_tile()+scale_fill_gradient2(low="red", high="green", mid="yellow", midpoint = max(AllCallACC.long[which(AllCallACC.long$Context=="CpG"),]$Depth)/2, limit=c(0, max(AllCallACC.long[which(AllCallACC.long$Context=="CpG"),]$Depth)), space="Lab", name="Depth")+labs(x="Contamination Blood: P (%)", y="")+scale_color_manual(values=DeconEstimate)+theme(text=element_text(size=18))+facet_wrap(nrow=1, "Context")
ggsave("cpgDEPTH-Tile-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/Sites", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(AllCallACC.long[which(AllCallACC.long$Context=="Non-CpG"),], aes(x=as.factor(P), y=Estimate, fill=Depth))+geom_tile()+scale_fill_gradient2(low="red", high="green", mid="yellow", midpoint = max(AllCallACC.long[which(AllCallACC.long$Context=="Non-CpG"),]$Depth)/2, limit=c(0, max(AllCallACC.long[which(AllCallACC.long$Context=="Non-CpG"),]$Depth)), space="Lab", name="Depth")+labs(x="Contamination Blood: P (%)", y="")+scale_color_manual(values=DeconEstimate)+theme(text=element_text(size=18))+facet_wrap(nrow=1, "Context")
ggsave("ncpgDEPTH-Tile-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/Sites", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)

ggplot(AllCallACC.long, aes(x=Depth, y=CCCP, groups=Estimate))+geom_point(aes(color=P, shape=Estimate), size=6)+scale_color_gradient2(low="green", high="red", mid="yellow", midpoint=50, limit=c(0,100))+scale_shape_manual(values=c(15,18,19))+labs(x="Depth", y="Pearson's Correlation Coefficient (r)")+facet_wrap("Context")
ggsave("mMetricAllDEPTH-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/Sites", scale = 1, width = 7.5, height = 3.5, units = c("in"), dpi=300, limitsize = TRUE)

#Make below a graph of direction of Bias
ggplot(AllCallACC.long, aes(x=P, y=CB, group=Estimate))+geom_line(aes(color=Estimate), linewidth=3, alpha=0.6)+labs(x="Contamination Blood: P (%)", y=expression("C"[b]))+scale_color_manual(values=DeconEstimate)+theme(text=element_text(size=18))+facet_wrap(nrow=1, "Context")
ggsave("allcCB-Line-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/Sites", scale = 1, width = 7.5, height = 3, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(AllCallACC.long, aes(x=P, y=Estimate, fill=CB))+geom_tile()+facet_wrap(nrow=1, "Context")+scale_fill_gradient2(low = "red", high = "green", mid = "yellow", midpoint = 0.5, limit = c(0,1), space = "Lab", name=expression("C"[b]))+labs(x="Contamination Blood: P (%)", y="")+theme(text=element_text(size=18))
ggsave("allcCB-Tile-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/Sites", scale = 1, width = 7.5, height = 3, units = c("in"), dpi=300, limitsize = TRUE)+theme(text=element_text(size=18))
ggplot(AllCallACC.long, aes(x=as.factor(P), y=Accuracy, group=Estimate, ymin=LowerConf, ymax=UpperConf))+geom_point(aes(color=Estimate), position=position_dodge(width=0.8), size=2.3, shape=18)+geom_errorbar(aes(color=Estimate), position=position_dodge(width=0.8))+labs(x="Contamination Blood: P (%)", y="Concordant Brain")+scale_color_manual(values=DeconEstimate)+theme(text=element_text(size=18))+facet_wrap(nrow=2, "Context", scales="free_y")
ggsave("allcCONCORD-Point-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/Sites", scale = 1, width = 7.5, height = 3, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(AllCallACC.long, aes(x=P, y=Estimate, fill=Accuracy))+facet_wrap("Context")+geom_tile()+scale_fill_gradient2(low = "red", high = "green", mid="yellow", midpoint = 0.9, limit = c(0.80,1), space = "Lab", name="∪")+labs(x="Contamination Blood: P (%)", y="")+theme(text=element_text(size=18))
ggsave("allCCONCORD-Tile-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/Sites", scale = 1, width = 7.5, height = 3, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(AllCallACC.long[which(AllCallACC.long$Context=="CpG"),], aes(x=as.factor(P), y=Estimate, fill=Accuracy))+facet_wrap("Context")+geom_tile()+scale_fill_gradient2(low = "red", high = "green", mid="yellow", midpoint = 0.94, limit = c(0.88,1), space = "Lab", name="Concordant (%)")+labs(x="Contamination Blood: P (%)", y="Concordant")
ggsave("cpgCONCORD-Tile-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/Sites", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(AllCallACC.long[which(AllCallACC.long$Context=="Non-CpG"),], aes(x=as.factor(P), y=Estimate, fill=Accuracy))+geom_tile()+facet_wrap("Context")+scale_fill_gradient2(low = "red", high = "green", mid="yellow", midpoint = 0.9985, limit = c(0.997,1), space = "Lab", name="Concordant (%)")+labs(x="Contamination Blood: P (%)", y="Concordant")
ggsave("ncpgCONCORD-Tile-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/Sites", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(AllCallSCC.long, aes(x=as.factor(P), y=Accuracy, group=Estimate, ymin=LowerConf, ymax=UpperConf))+geom_point(aes(color=Estimate), position=position_dodge(width=0.8), size=2.3, shape=18)+geom_errorbar(aes(color=Estimate), position=position_dodge(width=0.8))+labs(x="Contamination Blood: P (%)", y="Discordant Blood")+scale_color_manual(values=DeconEstimate)+theme(text=element_text(size=18))+facet_wrap(nrow=2, "Context", scales="free_y")
ggsave("cpgDISCORDv-Point-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/Sites", scale = 1, width = 7.5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(AllCallSCC.long, aes(x=as.factor(P), y=Estimate, fill=Accuracy))+geom_tile()+facet_wrap(nrow=1, "Context")+scale_fill_gradient2(low = "red", high = "green", mid = "yellow", midpoint = 0.5, limit = c(0,1), space = "Lab", name="Discordant (%)")+labs(x="Contamination Blood: P (%)", y="Discordant Blood")
ggsave("allcDISCORDv-Point-Meth.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/Images/Sites", scale = 1, width = 7.5, height = 3, units = c("in"), dpi=300, limitsize = TRUE)
