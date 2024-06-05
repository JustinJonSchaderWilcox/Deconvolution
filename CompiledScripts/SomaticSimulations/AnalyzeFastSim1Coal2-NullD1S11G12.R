#Anaylze and Create Summary Graphs in R
library(ggplot2)

FastSimSummary=read.delim('C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/SomaticSimulation/Data/FastSimCoal2/Summary/FastSimCoalSomaticAF-D1S100G12MUall.summary.txt')

#note that all results are negative even without a correction
summary.glm(glm(Mutations~Mu, family=poisson, data=FastSimSummary))
with(summary.glm(glm(Mutations~Mu, family=poisson, data=FastSimSummary)), 1-deviance/null.deviance)

summary(lm(afMEAN~Mu, data=FastSimSummary))
summary(lm(afMEDIAN~Mu, data=FastSimSummary))
summary(lm(afMODE~Mu, data=FastSimSummary))

pairwise.wilcox.test(FastSimSummary$Mutations, as.factor(FastSimSummary$Mu), p.adjust.method="bonferroni")

pairwise.wilcox.test(FastSimSummary$afMEAN, as.factor(FastSimSummary$Mu), p.adjust.method="bonferroni")
pairwise.wilcox.test(FastSimSummary$afMEDIAN, as.factor(FastSimSummary$Mu), p.adjust.method="bonferroni")
pairwise.wilcox.test(FastSimSummary$afMODE, as.factor(FastSimSummary$Mu), p.adjust.method="bonferroni")
pairwise.wilcox.test(FastSimSummary$afMODEprop, as.factor(FastSimSummary$Mu), p.adjust.method="bonferroni")
pairwise.wilcox.test(FastSimSummary$afMEDIANAbsDEV, as.factor(FastSimSummary$Mu), p.adjust.method="bonferroni")

pairwise.wilcox.test(FastSimSummary$FOLDEDafMEAN, as.factor(FastSimSummary$Mu), p.adjust.method="bonferroni")
pairwise.wilcox.test(FastSimSummary$FOLDEDafMEDIAN, as.factor(FastSimSummary$Mu), p.adjust.method="bonferroni")
pairwise.wilcox.test(FastSimSummary$FOLDEDafMODE, as.factor(FastSimSummary$Mu), p.adjust.method="bonferroni")
pairwise.wilcox.test(FastSimSummary$FOLDEDafMODEprop, as.factor(FastSimSummary$Mu), p.adjust.method="bonferroni")
pairwise.wilcox.test(FastSimSummary$FOLDEDafMEDIANAbsDEV, as.factor(FastSimSummary$Mu), p.adjust.method="bonferroni")

pairwise.wilcox.test(FastSimSummary$cut25PROP, as.factor(FastSimSummary$Mu), p.adjust.method="bonferroni")

ggplot(FastSimSummary, aes(x=as.factor(Mu), y=Mutations))+geom_boxplot()+labs(x="μ", y="Mutations")+theme(text=element_text(size=14), axis.text.x=element_text(size=6), strip.text.x=element_text(face="bold"))
ggsave("somMutations-x-MuBxPl.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/SomaticSimulation/Analyses/Images", scale = 1, width = 5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(FastSimSummary, aes(x=as.factor(Mu), y=log10(Mutations)))+geom_boxplot()+labs(x="μ", y=expression(Log[10](Mutations)))+theme(text=element_text(size=14), axis.text.x=element_text(size=6), strip.text.x=element_text(face="bold"))
ggsave("somMutationsLog10-x-MuBxPl.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/SomaticSimulation/Analyses/Images", scale = 1, width = 5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(FastSimSummary, aes(x=as.factor(Mu), y=afMEAN))+geom_boxplot()+labs(x="μ", y="Mean Allele Frequency")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"))
ggsave("afMEAN-x-MuBxPl.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/SomaticSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(FastSimSummary, aes(x=as.factor(Mu), y=afMEDIAN))+geom_boxplot()+labs(x="μ", y="Median Allele Frequency")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"))
ggsave("afMEDIAN-x-MuBxPl.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/SomaticSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(FastSimSummary, aes(x=as.factor(Mu), y=afMODE))+geom_boxplot()+labs(x="μ", y="Modal Allele Frequency")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"))
ggsave("afMODE-x-MuBxPl.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/SomaticSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(FastSimSummary, aes(x=as.factor(Mu), y=afMEDIANAbsDEV))+geom_boxplot()+labs(x="μ", y="AAD of Allele Freq")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"))
ggsave("afAbsMedDEV-x-MuBxPl.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/SomaticSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(FastSimSummary, aes(x=as.factor(Mu), y=FOLDEDafMEAN))+geom_boxplot()+labs(x="μ", y="Mean MAF Allele Frequency")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"))
ggsave("fafMEAN-x-MuBxPl.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/SomaticSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(FastSimSummary, aes(x=as.factor(Mu), y=FOLDEDafMEDIAN))+geom_boxplot()+labs(x="μ", y="Median MAF Allele Frequency")+theme(text=element_text(size=14), axis.text.x=element_text(size=6), strip.text.x=element_text(face="bold"))
ggsave("fafMEDIAN-x-MuBxPl.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/SomaticSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(FastSimSummary, aes(x=as.factor(Mu), y=FOLDEDafMODE))+geom_boxplot()+labs(x="μ", y="Mode MAF Allele Freq")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"))
ggsave("fafMODE-x-MuBxPl.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/SomaticSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(FastSimSummary, aes(x=as.factor(Mu), y=FOLDEDafMEDIANAbsDEV))+geom_boxplot()+labs(x="μ", y="AAD MAF Allele Freq")+theme(text=element_text(size=14), axis.text.x=element_text(size=8), strip.text.x=element_text(face="bold"))
ggsave("fafAbsMedDEV-x-MuBxPl.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/SomaticSimulation/Analyses/Images", scale = 1, width = 7, height = 5, units = c("in"), dpi=300, limitsize = TRUE)
ggplot(FastSimSummary, aes(x=as.factor(Mu), y=cut25PROP))+geom_boxplot()+labs(x="μ", y="MAF Allele Freq<0.25")+theme(text=element_text(size=14), axis.text.x=element_text(size=6), strip.text.x=element_text(face="bold"))
ggsave("cutPROP25-x-MuBxPl.png", plot = last_plot(), device = png(), path = "C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/SomaticSimulation/Analyses/Images", scale = 1, width = 5, height = 5, units = c("in"), dpi=300, limitsize = TRUE)

range(FastSimSummary$mean)

min(FastSimSummary[which(FastSimSummary$Mu==2.3e-9),]$Mutations)
max(FastSimSummary[which(FastSimSummary$Mu==4.6e-8),]$Mutations)
mean(FastSimSummary$afMEAN)
range(FastSimSummary$afMEAN)
mean(FastSimSummary$afMEDIAN)

mean(FastSimSummary$afMODE)
range(FastSimSummary$afMODE)

mean(FastSimSummary$cut25PROP)
range(FastSimSummary$cut25PROP)
