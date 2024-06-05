#Make Methylation Calls across whole genome for future analyses
library(data.table)

system("subst x: \"C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation\"")
setwd("x:/Subset/WholeGenome/DeconvolutedBrain")
DeconBrainWGoCpG0=read.delim("./CpG/mAbel-BrainMix.Blood-0.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.CpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainWGoCpG10=read.delim("./CpG/mAbel-BrainMix.Blood-10.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.CpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainWGoCpG20=read.delim("./CpG/mAbel-BrainMix.Blood-20.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.CpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainWGoCpG30=read.delim("./CpG/mAbel-BrainMix.Blood-30.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.CpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainWGoCpG40=read.delim("./CpG/mAbel-BrainMix.Blood-40.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.CpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainWGoCpG50=read.delim("./CpG/mAbel-BrainMix.Blood-50.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.CpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainWGoCpG60=read.delim("./CpG/mAbel-BrainMix.Blood-60.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.CpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainWGoCpG70=read.delim("./CpG/mAbel-BrainMix.Blood-70.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.CpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainWGoCpG80=read.delim("./CpG/mAbel-BrainMix.Blood-80.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.CpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainWGoCpG90=read.delim("./CpG/mAbel-BrainMix.Blood-90.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.CpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainWGoCpG100=read.delim("./CpG/mAbel-BrainMix.Blood-100.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.CpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))

DeconBrainWGnCpG0=read.delim("./nCpG/mAbel-BrainMix.Blood-0.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.nCpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainWGnCpG10=read.delim("./nCpG/mAbel-BrainMix.Blood-10.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.nCpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainWGnCpG20=read.delim("./nCpG/mAbel-BrainMix.Blood-20.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.nCpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainWGnCpG30=read.delim("./nCpG/mAbel-BrainMix.Blood-30.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.nCpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainWGnCpG40=read.delim("./nCpG/mAbel-BrainMix.Blood-40.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.nCpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainWGnCpG50=read.delim("./nCpG/mAbel-BrainMix.Blood-50.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.nCpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainWGnCpG60=read.delim("./nCpG/mAbel-BrainMix.Blood-60.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.nCpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainWGnCpG70=read.delim("./nCpG/mAbel-BrainMix.Blood-70.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.nCpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainWGnCpG80=read.delim("./nCpG/mAbel-BrainMix.Blood-80.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.nCpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainWGnCpG90=read.delim("./nCpG/mAbel-BrainMix.Blood-90.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.nCpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainWGnCpG100=read.delim("./nCpG/mAbel-BrainMix.Blood-100.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.nCpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))

DeconBrainMethWGoCpG0=read.delim("./CpG/mAbel-BrainMix.Blood-0.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-CpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainMethWGoCpG10=read.delim("./CpG/mAbel-BrainMix.Blood-10.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-CpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainMethWGoCpG20=read.delim("./CpG/mAbel-BrainMix.Blood-20.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-CpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainMethWGoCpG30=read.delim("./CpG/mAbel-BrainMix.Blood-30.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-CpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainMethWGoCpG40=read.delim("./CpG/mAbel-BrainMix.Blood-40.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-CpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainMethWGoCpG50=read.delim("./CpG/mAbel-BrainMix.Blood-50.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-CpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainMethWGoCpG60=read.delim("./CpG/mAbel-BrainMix.Blood-60.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-CpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainMethWGoCpG70=read.delim("./CpG/mAbel-BrainMix.Blood-70.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-CpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainMethWGoCpG80=read.delim("./CpG/mAbel-BrainMix.Blood-80.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-CpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainMethWGoCpG90=read.delim("./CpG/mAbel-BrainMix.Blood-90.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-CpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainMethWGoCpG100=read.delim("./CpG/mAbel-BrainMix.Blood-100.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-CpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))

DeconBrainMethWGnCpG0=read.delim("./nCpG/mAbel-BrainMix.Blood-0.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-nCpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainMethWGnCpG10=read.delim("./nCpG/mAbel-BrainMix.Blood-10.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-nCpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainMethWGnCpG20=read.delim("./nCpG/mAbel-BrainMix.Blood-20.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-nCpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainMethWGnCpG30=read.delim("./nCpG/mAbel-BrainMix.Blood-30.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-nCpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainMethWGnCpG40=read.delim("./nCpG/mAbel-BrainMix.Blood-40.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-nCpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainMethWGnCpG50=read.delim("./nCpG/mAbel-BrainMix.Blood-50.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-nCpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainMethWGnCpG60=read.delim("./nCpG/mAbel-BrainMix.Blood-60.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-nCpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainMethWGnCpG70=read.delim("./nCpG/mAbel-BrainMix.Blood-70.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-nCpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainMethWGnCpG80=read.delim("./nCpG/mAbel-BrainMix.Blood-80.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-nCpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainMethWGnCpG90=read.delim("./nCpG/mAbel-BrainMix.Blood-90.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-nCpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))
DeconBrainMethWGnCpG100=read.delim("./nCpG/mAbel-BrainMix.Blood-100.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-nCpG.shared-site30K.bed", header=FALSE, sep="\t", col.names=c("Chromosome", "Position","Position2", "DeconMethProportion", "DeconMethReads", "DeconUnmethReads", "BrainMethProportion", "BrainMethReads", "BrainUnmethReads", "BloodMethProportion", "BloodMethReads", "BloodUnmethReads", "MixMethProportion", "MixMethReads", "MixUnmethReads", "tpMethProportion", "tpMethReads", "tpUnmethReads", "P"))

DeconBrainIslandCpG=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/CpG-Island/mAbel-BrainMix.Blood-AllIslandShared-deduplicated.bismark.deconvoluted-brain-cov.CpG-Island_iCpG.bed", header=TRUE)
DeconBrainIslandNonCpG=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/CpG-Island/mAbel-BrainMix.Blood-AllIslandShared-deduplicated.bismark.deconvoluted-brain-cov.CpG-Island_nCpG.bed", header=TRUE)

DeconBrainGeneCpG=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/GeneBody/mAbel-BrainMix.Blood-AllGeneShared-deduplicated.bismark.deconvoluted-brain-cov.CpG-GeneBody.bed", header=TRUE)
DeconBrainGeneNonCpG=read.delim("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/GeneBody/mAbel-BrainMix.Blood-AllGeneShared-deduplicated.bismark.deconvoluted-brain-cov.nCpG-GeneBody.bed", header=TRUE)

setwd("C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation")
system("subst x: /D")

MergedDeconBrainWGoCpG=rbind(DeconBrainWGoCpG0, DeconBrainWGoCpG10, DeconBrainWGoCpG20, DeconBrainWGoCpG30, DeconBrainWGoCpG40, DeconBrainWGoCpG50, DeconBrainWGoCpG60, DeconBrainWGoCpG70, DeconBrainWGoCpG80, DeconBrainWGoCpG90, DeconBrainWGoCpG100)
MergedDeconBrainWGnCpG=rbind(DeconBrainWGnCpG0, DeconBrainWGnCpG10, DeconBrainWGnCpG20, DeconBrainWGnCpG30, DeconBrainWGnCpG40, DeconBrainWGnCpG50, DeconBrainWGnCpG60, DeconBrainWGnCpG70, DeconBrainWGnCpG80, DeconBrainWGnCpG90, DeconBrainWGnCpG100)
MergedPurgedDeconBrainWGoCpG=setDT(MergedDeconBrainWGoCpG)[,.(Chromosome=unique(Chromosome), DeconMethProportion=mean(DeconMethProportion), BrainMethProportion=mean(BrainMethProportion), BloodMethProportion=mean(BloodMethProportion), MixMethProportion=mean(MixMethProportion), tpMethProportion=mean(tpMethProportion)), by=Position]
MergedPurgedDeconBrainWGnCpG=setDT(MergedDeconBrainWGnCpG)[,.(Chromosome=unique(Chromosome), DeconMethProportion=mean(DeconMethProportion), BrainMethProportion=mean(BrainMethProportion), BloodMethProportion=mean(BloodMethProportion), MixMethProportion=mean(MixMethProportion), tpMethProportion=mean(tpMethProportion)), by=Position]

MergedDeconBrainWGoCpG$DECONvBRAINdMP=1
for (i in 1:dim(MergedDeconBrainWGoCpG)[1]){
MergedDeconBrainWGoCpG$DECONvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainWGoCpG$DeconMethReads[i], MergedDeconBrainWGoCpG$DeconUnmethReads[i]), c(MergedDeconBrainWGoCpG$BrainMethReads[i], MergedDeconBrainWGoCpG$BrainUnmethReads[i]))))$p.value
}
MergedDeconBrainWGnCpG$DECONvBRAINdMP=1
for (i in 1:dim(MergedDeconBrainWGnCpG)[1]){
MergedDeconBrainWGnCpG$DECONvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainWGnCpG$DeconMethReads[i], MergedDeconBrainWGnCpG$DeconUnmethReads[i]), c(MergedDeconBrainWGnCpG$BrainMethReads[i], MergedDeconBrainWGnCpG$BrainUnmethReads[i]))))$p.value
}

MergedDeconBrainWGoCpG$tpDECONvBRAINdMP=1
for (i in 1:dim(MergedDeconBrainWGoCpG)[1]){
MergedDeconBrainWGoCpG$tpDECONvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainWGoCpG$tpMethReads[i], MergedDeconBrainWGoCpG$tpUnmethReads[i]), c(MergedDeconBrainWGoCpG$BrainMethReads[i], MergedDeconBrainWGoCpG$BrainUnmethReads[i]))))$p.value
}
MergedDeconBrainWGnCpG$tpDECONvBRAINdMP=1
for (i in 1:dim(MergedDeconBrainWGnCpG)[1]){
MergedDeconBrainWGnCpG$tpDECONvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainWGnCpG$tpMethReads[i], MergedDeconBrainWGnCpG$tpUnmethReads[i]), c(MergedDeconBrainWGnCpG$BrainMethReads[i], MergedDeconBrainWGnCpG$BrainUnmethReads[i]))))$p.value
}

MergedDeconBrainWGoCpG$MIXvBRAINdMP=1
for (i in 1:dim(MergedDeconBrainWGoCpG)[1]){
MergedDeconBrainWGoCpG$MIXvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainWGoCpG$MixMethReads[i], MergedDeconBrainWGoCpG$MixUnmethReads[i]), c(MergedDeconBrainWGoCpG$BrainMethReads[i], MergedDeconBrainWGoCpG$BrainUnmethReads[i]))))$p.value
}
MergedDeconBrainWGnCpG$MIXvBRAINdMP=1
for (i in 1:dim(MergedDeconBrainWGnCpG)[1]){
MergedDeconBrainWGnCpG$MIXvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainWGnCpG$MixMethReads[i], MergedDeconBrainWGnCpG$MixUnmethReads[i]), c(MergedDeconBrainWGnCpG$BrainMethReads[i], MergedDeconBrainWGnCpG$BrainUnmethReads[i]))))$p.value
}


MergedDeconBrainWGoCpG$TrueDiffMethP=1
for (i in 1:dim(MergedDeconBrainWGoCpG)[1]){
MergedDeconBrainWGoCpG$TrueDiffMethP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainWGoCpG$BrainMethReads[i], MergedDeconBrainWGoCpG$BrainUnmethReads[i]), c(MergedDeconBrainWGoCpG$BloodMethReads[i], MergedDeconBrainWGoCpG$BloodUnmethReads[i]))))$p.value
}
MergedDeconBrainWGnCpG$TrueDiffMethP=1
for (i in 1:dim(MergedDeconBrainWGnCpG)[1]){
MergedDeconBrainWGnCpG$TrueDiffMethP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainWGnCpG$BrainMethReads[i], MergedDeconBrainWGnCpG$BrainUnmethReads[i]), c(MergedDeconBrainWGnCpG$BloodMethReads[i], MergedDeconBrainWGnCpG$BloodUnmethReads[i]))))$p.value
}

MergedDeconBrainWGoCpG$PredictDiffMethP=1
for (i in 1:dim(MergedDeconBrainWGoCpG)[1]){
MergedDeconBrainWGoCpG$PredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainWGoCpG$DeconMethReads[i], MergedDeconBrainWGoCpG$DeconUnmethReads[i]), c(MergedDeconBrainWGoCpG$BloodMethReads[i], MergedDeconBrainWGoCpG$BloodUnmethReads[i]))))$p.value
}
MergedDeconBrainWGnCpG$PredictDiffMethP=1
for (i in 1:dim(MergedDeconBrainWGnCpG)[1]){
MergedDeconBrainWGnCpG$PredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainWGnCpG$DeconMethReads[i], MergedDeconBrainWGnCpG$DeconUnmethReads[i]), c(MergedDeconBrainWGnCpG$BloodMethReads[i], MergedDeconBrainWGnCpG$BloodUnmethReads[i]))))$p.value
}

MergedDeconBrainWGoCpG$tpPredictDiffMethP=1
for (i in 1:dim(MergedDeconBrainWGoCpG)[1]){
MergedDeconBrainWGoCpG$tpPredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainWGoCpG$tpMethReads[i], MergedDeconBrainWGoCpG$tpUnmethReads[i]), c(MergedDeconBrainWGoCpG$BloodMethReads[i], MergedDeconBrainWGoCpG$BloodUnmethReads[i]))))$p.value
}
MergedDeconBrainWGnCpG$tpPredictDiffMethP=1
for (i in 1:dim(MergedDeconBrainWGnCpG)[1]){
MergedDeconBrainWGnCpG$tpPredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainWGnCpG$tpMethReads[i], MergedDeconBrainWGnCpG$tpUnmethReads[i]), c(MergedDeconBrainWGnCpG$BloodMethReads[i], MergedDeconBrainWGnCpG$BloodUnmethReads[i]))))$p.value
}

MergedDeconBrainWGoCpG$MixDiffMethP=1
for (i in 1:dim(MergedDeconBrainWGoCpG)[1]){
MergedDeconBrainWGoCpG$MixDiffMethP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainWGoCpG$MixMethReads[i], MergedDeconBrainWGoCpG$MixUnmethReads[i]), c(MergedDeconBrainWGoCpG$BloodMethReads[i], MergedDeconBrainWGoCpG$BloodUnmethReads[i]))))$p.value
}
MergedDeconBrainWGnCpG$MixDiffMethP=1
for (i in 1:dim(MergedDeconBrainWGnCpG)[1]){
MergedDeconBrainWGnCpG$MixDiffMethP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainWGnCpG$MixMethReads[i], MergedDeconBrainWGnCpG$MixUnmethReads[i]), c(MergedDeconBrainWGnCpG$BloodMethReads[i], MergedDeconBrainWGnCpG$BloodUnmethReads[i]))))$p.value
}

write.table(MergedDeconBrainWGoCpG, file="C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/MergedDeconBrainWGoCpG.TRUE-P.txt", sep="\t")
write.table(MergedDeconBrainWGnCpG, file="C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/MergedDeconBrainWGnCpG.TRUE-P.txt", sep="\t")

MergedDeconBrainMethWGoCpG=rbind(DeconBrainMethWGoCpG0, DeconBrainMethWGoCpG10, DeconBrainMethWGoCpG20, DeconBrainMethWGoCpG30, DeconBrainMethWGoCpG40, DeconBrainMethWGoCpG50, DeconBrainMethWGoCpG60, DeconBrainMethWGoCpG70, DeconBrainMethWGoCpG80, DeconBrainMethWGoCpG90, DeconBrainMethWGoCpG100)
MergedDeconBrainMethWGnCpG=rbind(DeconBrainMethWGnCpG0, DeconBrainMethWGnCpG10, DeconBrainMethWGnCpG20, DeconBrainMethWGnCpG30, DeconBrainMethWGnCpG40, DeconBrainMethWGnCpG50, DeconBrainMethWGnCpG60, DeconBrainMethWGnCpG70, DeconBrainMethWGnCpG80, DeconBrainMethWGnCpG90, DeconBrainMethWGnCpG100)
MergedPurgedDeconBrainMethWGoCpG=setDT(MergedDeconBrainMethWGoCpG)[,.(Chromosome=unique(Chromosome), DeconMethProportion=mean(DeconMethProportion), BrainMethProportion=mean(BrainMethProportion), BloodMethProportion=mean(BloodMethProportion), MixMethProportion=mean(MixMethProportion), tpMethProportion=mean(tpMethProportion)), by=Position]
MergedPurgedDeconBrainMethWGnCpG=setDT(MergedDeconBrainMethWGnCpG)[,.(Chromosome=unique(Chromosome), DeconMethProportion=mean(DeconMethProportion), BrainMethProportion=mean(BrainMethProportion), BloodMethProportion=mean(BloodMethProportion), MixMethProportion=mean(MixMethProportion), tpMethProportion=mean(tpMethProportion)), by=Position]


MergedDeconBrainMethWGoCpG$DECONvBRAINdMP=1
for (i in 1:dim(MergedDeconBrainMethWGoCpG)[1]){
MergedDeconBrainMethWGoCpG$DECONvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainMethWGoCpG$DeconMethReads[i], MergedDeconBrainMethWGoCpG$DeconUnmethReads[i]), c(MergedDeconBrainMethWGoCpG$BrainMethReads[i], MergedDeconBrainMethWGoCpG$BrainUnmethReads[i]))))$p.value
}
MergedDeconBrainMethWGnCpG$DECONvBRAINdMP=1
for (i in 1:dim(MergedDeconBrainMethWGnCpG)[1]){
MergedDeconBrainMethWGnCpG$DECONvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainMethWGnCpG$DeconMethReads[i], MergedDeconBrainMethWGnCpG$DeconUnmethReads[i]), c(MergedDeconBrainMethWGnCpG$BrainMethReads[i], MergedDeconBrainMethWGnCpG$BrainUnmethReads[i]))))$p.value
}

MergedDeconBrainMethWGoCpG$tpDECONvBRAINdMP=1
for (i in 1:dim(MergedDeconBrainMethWGoCpG)[1]){
MergedDeconBrainMethWGoCpG$tpDECONvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainMethWGoCpG$tpMethReads[i], MergedDeconBrainMethWGoCpG$tpUnmethReads[i]), c(MergedDeconBrainMethWGoCpG$BrainMethReads[i], MergedDeconBrainMethWGoCpG$BrainUnmethReads[i]))))$p.value
}
MergedDeconBrainMethWGnCpG$tpDECONvBRAINdMP=1
for (i in 1:dim(MergedDeconBrainMethWGnCpG)[1]){
MergedDeconBrainMethWGnCpG$tpDECONvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainMethWGnCpG$tpMethReads[i], MergedDeconBrainMethWGnCpG$tpUnmethReads[i]), c(MergedDeconBrainMethWGnCpG$BrainMethReads[i], MergedDeconBrainMethWGnCpG$BrainUnmethReads[i]))))$p.value
}

MergedDeconBrainMethWGoCpG$MIXvBRAINdMP=1
for (i in 1:dim(MergedDeconBrainMethWGoCpG)[1]){
MergedDeconBrainMethWGoCpG$MIXvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainMethWGoCpG$MixMethReads[i], MergedDeconBrainMethWGoCpG$MixUnmethReads[i]), c(MergedDeconBrainMethWGoCpG$BrainMethReads[i], MergedDeconBrainMethWGoCpG$BrainUnmethReads[i]))))$p.value
}
MergedDeconBrainMethWGnCpG$MIXvBRAINdMP=1
for (i in 1:dim(MergedDeconBrainMethWGnCpG)[1]){
MergedDeconBrainMethWGnCpG$MIXvBRAINdMP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainMethWGnCpG$MixMethReads[i], MergedDeconBrainMethWGnCpG$MixUnmethReads[i]), c(MergedDeconBrainMethWGnCpG$BrainMethReads[i], MergedDeconBrainMethWGnCpG$BrainUnmethReads[i]))))$p.value
}


MergedDeconBrainMethWGoCpG$TrueDiffMethP=1
for (i in 1:dim(MergedDeconBrainMethWGoCpG)[1]){
MergedDeconBrainMethWGoCpG$TrueDiffMethP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainMethWGoCpG$BrainMethReads[i], MergedDeconBrainMethWGoCpG$BrainUnmethReads[i]), c(MergedDeconBrainMethWGoCpG$BloodMethReads[i], MergedDeconBrainMethWGoCpG$BloodUnmethReads[i]))))$p.value
}
MergedDeconBrainMethWGnCpG$TrueDiffMethP=1
for (i in 1:dim(MergedDeconBrainMethWGnCpG)[1]){
MergedDeconBrainMethWGnCpG$TrueDiffMethP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainMethWGnCpG$BrainMethReads[i], MergedDeconBrainMethWGnCpG$BrainUnmethReads[i]), c(MergedDeconBrainMethWGnCpG$BloodMethReads[i], MergedDeconBrainMethWGnCpG$BloodUnmethReads[i]))))$p.value
}

MergedDeconBrainMethWGoCpG$PredictDiffMethP=1
for (i in 1:dim(MergedDeconBrainMethWGoCpG)[1]){
MergedDeconBrainMethWGoCpG$PredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainMethWGoCpG$DeconMethReads[i], MergedDeconBrainMethWGoCpG$DeconUnmethReads[i]), c(MergedDeconBrainMethWGoCpG$BloodMethReads[i], MergedDeconBrainMethWGoCpG$BloodUnmethReads[i]))))$p.value
}
MergedDeconBrainMethWGnCpG$PredictDiffMethP=1
for (i in 1:dim(MergedDeconBrainMethWGnCpG)[1]){
MergedDeconBrainMethWGnCpG$PredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainMethWGnCpG$DeconMethReads[i], MergedDeconBrainMethWGnCpG$DeconUnmethReads[i]), c(MergedDeconBrainMethWGnCpG$BloodMethReads[i], MergedDeconBrainMethWGnCpG$BloodUnmethReads[i]))))$p.value
}

MergedDeconBrainMethWGoCpG$tpPredictDiffMethP=1
for (i in 1:dim(MergedDeconBrainMethWGoCpG)[1]){
MergedDeconBrainMethWGoCpG$tpPredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainMethWGoCpG$tpMethReads[i], MergedDeconBrainMethWGoCpG$tpUnmethReads[i]), c(MergedDeconBrainMethWGoCpG$BloodMethReads[i], MergedDeconBrainMethWGoCpG$BloodUnmethReads[i]))))$p.value
}
MergedDeconBrainMethWGnCpG$tpPredictDiffMethP=1
for (i in 1:dim(MergedDeconBrainMethWGnCpG)[1]){
MergedDeconBrainMethWGnCpG$tpPredictDiffMethP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainMethWGnCpG$tpMethReads[i], MergedDeconBrainMethWGnCpG$tpUnmethReads[i]), c(MergedDeconBrainMethWGnCpG$BloodMethReads[i], MergedDeconBrainMethWGnCpG$BloodUnmethReads[i]))))$p.value
}

MergedDeconBrainMethWGoCpG$MixDiffMethP=1
for (i in 1:dim(MergedDeconBrainMethWGoCpG)[1]){
MergedDeconBrainMethWGoCpG$MixDiffMethP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainMethWGoCpG$MixMethReads[i], MergedDeconBrainMethWGoCpG$MixUnmethReads[i]), c(MergedDeconBrainMethWGoCpG$BloodMethReads[i], MergedDeconBrainMethWGoCpG$BloodUnmethReads[i]))))$p.value
}
MergedDeconBrainMethWGnCpG$MixDiffMethP=1
for (i in 1:dim(MergedDeconBrainMethWGnCpG)[1]){
MergedDeconBrainMethWGnCpG$MixDiffMethP[i]=fisher.test(as.matrix(cbind(c(MergedDeconBrainMethWGnCpG$MixMethReads[i], MergedDeconBrainMethWGnCpG$MixUnmethReads[i]), c(MergedDeconBrainMethWGnCpG$BloodMethReads[i], MergedDeconBrainMethWGnCpG$BloodUnmethReads[i]))))$p.value
}

write.table(MergedDeconBrainMethWGoCpG, file="C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/MergedDeconBrainMethWGoCpG.TRUE-P.txt", sep="\t")
write.table(MergedDeconBrainMethWGnCpG, file="C:/Users/mjuswilc/OneDrive - tu-dortmund.de/Documents/TU-Dortmund/GreatTits/Projects/Deconvolution/UltimateDeconvolution/ultAnalysis/Validation/MethylationCalls/MergedDeconBrainMethWGnCpG.TRUE-P.txt", sep="\t")
