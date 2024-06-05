#!/bin/bash
#SBATCH --time=7:58:00
#SBATCH --mem=16GB
#SBATCH -p med
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=24
# Output and error files
#SBATCH -o job.%J.out

source activate Genomics

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SharedSites
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SharedSites/CpG
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SharedSites/nCpG
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SmallChromosome
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SmallChromosome/Random
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SmallChromosome/DeconvolutedBrain
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SmallChromosome/DeconvolutedBrain/CpG
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SmallChromosome/DeconvolutedBrain/nCpG
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/Random

CHROMOSOMES=`awk '{print $1}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/MethylationBed/Genome/Abel-GCF_001522545.3_Parus_major1.1_genomic.chrom.genome`

for c in $CHROMOSOMES
do
BRAINhDEPTH=`awk -v chrom=$c '$1==chrom{print int(($3*2)+0.5); exit}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_0.chrom-depth.txt`
BRAINlDEPTH=`awk -v chrom=$c '$1==chrom{print int(($3*0.5)+0.5); exit}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_0.chrom-depth.txt`
BLOODhDEPTH=`awk -v chrom=$c '$1==chrom{print int(($3*2)+0.5); exit}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_100.chrom-depth.txt`
BLOODlDEPTH=`awk -v chrom=$c '$1==chrom{print int(($3*0.5)+0.5); exit}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_100.chrom-depth.txt`
bedtools intersect -sorted -f 1 -F 1 -a \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/0/mAbel-BrainMix.Blood-0.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz -b \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/10/mAbel-BrainMix.Blood-10.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/20/mAbel-BrainMix.Blood-20.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/30/mAbel-BrainMix.Blood-30.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/40/mAbel-BrainMix.Blood-40.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/50/mAbel-BrainMix.Blood-50.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/60/mAbel-BrainMix.Blood-60.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/70/mAbel-BrainMix.Blood-70.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/80/mAbel-BrainMix.Blood-80.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/90/mAbel-BrainMix.Blood-90.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/100/mAbel-BrainMix.Blood-100.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz | \
awk -v hBRAIN=$BRAINhDEPTH -v lBRAIN=$BRAINlDEPTH  -v hBLOOD=$BLOODhDEPTH -v lBLOOD=$BLOODlDEPTH 'BEGIN{OFS=FS="\t"}($8+$9)>hBRAIN || ($8+$9)<lBRAIN{next}($11+$12)>hBLOOD || ($11+$12)<lBLOOD{next}($8/($8+$9))>=0.1 || ($11/($11+$12))>=0.1{print $1 OFS $2 OFS $3}' | uniq | gzip -9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SharedSites/CpG/mAbel-BrainMix.$c.Meth-CpG.shared-sites.bed.gz &
JOBS=`jobs -l | grep -c '^\['`
while [ $JOBS -ge 2 ]
do
sleep 3s
JOBS=`jobs -l | grep -c '^\['`
done
done
for c in $CHROMOSOMES
do
BRAINhDEPTH=`awk -v chrom=$c '$1==chrom{print int(($3*2)+0.5); exit}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_0.chrom-depth.txt`
BRAINlDEPTH=`awk -v chrom=$c '$1==chrom{print int(($3*0.5)+0.5); exit}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_0.chrom-depth.txt`
BLOODhDEPTH=`awk -v chrom=$c '$1==chrom{print int(($3*2)+0.5); exit}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_100.chrom-depth.txt`
BLOODlDEPTH=`awk -v chrom=$c '$1==chrom{print int(($3*0.5)+0.5); exit}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_100.chrom-depth.txt`
bedtools intersect -sorted -f 1 -F 1 -a \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/0/mAbel-BrainMix.Blood-0.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz -b \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/10/mAbel-BrainMix.Blood-10.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/20/mAbel-BrainMix.Blood-20.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/30/mAbel-BrainMix.Blood-30.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/40/mAbel-BrainMix.Blood-40.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/50/mAbel-BrainMix.Blood-50.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/60/mAbel-BrainMix.Blood-60.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/70/mAbel-BrainMix.Blood-70.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/80/mAbel-BrainMix.Blood-80.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/90/mAbel-BrainMix.Blood-90.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/100/mAbel-BrainMix.Blood-100.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz | \
awk -v hBRAIN=$BRAINhDEPTH -v lBRAIN=$BRAINlDEPTH  -v hBLOOD=$BLOODhDEPTH -v lBLOOD=$BLOODlDEPTH 'BEGIN{OFS=FS="\t"}($8+$9)>hBRAIN || ($8+$9)<lBRAIN{next}($11+$12)>hBLOOD || ($11+$12)<lBLOOD{next}($8/($8+$9))>=0.1 || ($11/($11+$12))>=0.1{print $1 OFS $2 OFS $3}' | uniq | gzip -9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SharedSites/nCpG/mAbel-BrainMix.$c.Meth-nCpG.shared-sites.bed.gz &
JOBS=`jobs -l | grep -c '^\['`
while [ $JOBS -ge 2 ]
do
sleep 3s
JOBS=`jobs -l | grep -c '^\['`
done
done



for c in $CHROMOSOMES
do
BRAINhDEPTH=`awk -v chrom=$c '$1==chrom{print int(($3*2)+0.5); exit}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_0.chrom-depth.txt`
BRAINlDEPTH=`awk -v chrom=$c '$1==chrom{print int(($3*0.5)+0.5); exit}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_0.chrom-depth.txt`
BLOODhDEPTH=`awk -v chrom=$c '$1==chrom{print int(($3*2)+0.5); exit}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_100.chrom-depth.txt`
BLOODlDEPTH=`awk -v chrom=$c '$1==chrom{print int(($3*0.5)+0.5); exit}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_100.chrom-depth.txt`
bedtools intersect -sorted -f 1 -F 1 -a \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/0/mAbel-BrainMix.Blood-0.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz -b \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/10/mAbel-BrainMix.Blood-10.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/20/mAbel-BrainMix.Blood-20.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/30/mAbel-BrainMix.Blood-30.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/40/mAbel-BrainMix.Blood-40.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/50/mAbel-BrainMix.Blood-50.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/60/mAbel-BrainMix.Blood-60.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/70/mAbel-BrainMix.Blood-70.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/80/mAbel-BrainMix.Blood-80.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/90/mAbel-BrainMix.Blood-90.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/100/mAbel-BrainMix.Blood-100.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz | \
awk -v hBRAIN=$BRAINhDEPTH -v lBRAIN=$BRAINlDEPTH  -v hBLOOD=$BLOODhDEPTH -v lBLOOD=$BLOODlDEPTH 'BEGIN{OFS=FS="\t"}($8+$9)>hBRAIN || ($8+$9)<lBRAIN{next}($11+$12)>hBLOOD || ($11+$12)<lBLOOD{next}{print $1 OFS $2 OFS $3}' | uniq | gzip -9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SharedSites/CpG/mAbel-BrainMix.$c.CpG.shared-sites.bed.gz &
JOBS=`jobs -l | grep -c '^\['`
while [ $JOBS -ge 2 ]
do
sleep 3s
JOBS=`jobs -l | grep -c '^\['`
done
done
for c in $CHROMOSOMES
do
BRAINhDEPTH=`awk -v chrom=$c '$1==chrom{print int(($3*2)+0.5); exit}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_0.chrom-depth.txt`
BRAINlDEPTH=`awk -v chrom=$c '$1==chrom{print int(($3*0.5)+0.5); exit}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_0.chrom-depth.txt`
BLOODhDEPTH=`awk -v chrom=$c '$1==chrom{print int(($3*2)+0.5); exit}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_100.chrom-depth.txt`
BLOODlDEPTH=`awk -v chrom=$c '$1==chrom{print int(($3*0.5)+0.5); exit}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_100.chrom-depth.txt`
bedtools intersect -sorted -f 1 -F 1 -a \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/0/mAbel-BrainMix.Blood-0.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz -b \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/10/mAbel-BrainMix.Blood-10.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/20/mAbel-BrainMix.Blood-20.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/30/mAbel-BrainMix.Blood-30.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/40/mAbel-BrainMix.Blood-40.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/50/mAbel-BrainMix.Blood-50.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/60/mAbel-BrainMix.Blood-60.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/70/mAbel-BrainMix.Blood-70.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/80/mAbel-BrainMix.Blood-80.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/90/mAbel-BrainMix.Blood-90.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/100/mAbel-BrainMix.Blood-100.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz | \
awk -v hBRAIN=$BRAINhDEPTH -v lBRAIN=$BRAINlDEPTH  -v hBLOOD=$BLOODhDEPTH -v lBLOOD=$BLOODlDEPTH 'BEGIN{OFS=FS="\t"}($8+$9)>hBRAIN || ($8+$9)<lBRAIN{next}($11+$12)>hBLOOD || ($11+$12)<lBLOOD{next}{print $1 OFS $2 OFS $3}' | uniq | gzip -9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SharedSites/nCpG/mAbel-BrainMix.$c.nCpG.shared-sites.bed.gz &
JOBS=`jobs -l | grep -c '^\['`
while [ $JOBS -ge 2 ]
do
sleep 3s
JOBS=`jobs -l | grep -c '^\['`
done
done
wait

c=NC_031797.1
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SharedSites/CpG/mAbel-BrainMix.$c.Meth-CpG.shared-sites.bed.gz | shuf -n 30000 | sort -T /work/mjuswilc/tmp -n -k2 | gzip -9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SmallChromosome/Random/mAbel-BrainMix.$c.Meth-CpG.shared-sites30K.bed.gz &
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SharedSites/nCpG/mAbel-BrainMix.$c.Meth-nCpG.shared-sites.bed.gz | shuf -n 30000 | sort -T /work/mjuswilc/tmp -n -k2 | gzip -9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SmallChromosome/Random/mAbel-BrainMix.$c.Meth-nCpG.shared-sites30K.bed.gz &
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SharedSites/CpG/mAbel-BrainMix.$c.CpG.shared-sites.bed.gz | shuf -n 30000 | sort -T /work/mjuswilc/tmp -n -k2 | gzip -9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SmallChromosome/Random/mAbel-BrainMix.$c.CpG.shared-sites30K.bed.gz &
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SharedSites/nCpG/mAbel-BrainMix.$c.nCpG.shared-sites.bed.gz | shuf -n 30000 | sort -T /work/mjuswilc/tmp -n -k2 | gzip -9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SmallChromosome/Random/mAbel-BrainMix.$c.nCpG.shared-sites30K.bed.gz &
wait


for P in 0 10 20 30 40 50 60 70 80 90 100
do
bedtools intersect -sorted -f 1 -F 1 -a \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz -b \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SmallChromosome/Random/mAbel-BrainMix.$c.Meth-CpG.shared-sites30K.bed.gz \
> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SmallChromosome/DeconvolutedBrain/CpG/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.Meth-CpG.shared-site30K.bed &
bedtools intersect -sorted -f 1 -F 1 -a \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz -b \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SmallChromosome/Random/mAbel-BrainMix.$c.Meth-nCpG.shared-sites30K.bed.gz \
> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SmallChromosome/DeconvolutedBrain/nCpG/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.Meth-nCpG.shared-site30K.bed &
bedtools intersect -sorted -f 1 -F 1 -a \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz -b \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SmallChromosome/Random/mAbel-BrainMix.$c.CpG.shared-sites30K.bed.gz \
> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SmallChromosome/DeconvolutedBrain/CpG/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.shared-site30K.bed &
bedtools intersect -sorted -f 1 -F 1 -a \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz -b \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SmallChromosome/Random/mAbel-BrainMix.$c.nCpG.shared-sites30K.bed.gz \
> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SmallChromosome/DeconvolutedBrain/nCpG/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.shared-site30K.bed &
JOBS=`jobs -l | grep -c '^\['`
while [ $JOBS -ge 21 ]
do
sleep 3s
JOBS=`jobs -l | grep -c '^\['`
done
done
wait

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/TEMP
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/TEMP
for c in $CHROMOSOMES
do
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SharedSites/CpG/mAbel-BrainMix.$c.Meth-CpG.shared-sites.bed.gz >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/TEMP/mAbel-BrainMix.AllChromosome.Meth-CpG.shared-sites.bed &
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SharedSites/nCpG/mAbel-BrainMix.$c.Meth-nCpG.shared-sites.bed.gz >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/TEMP/mAbel-BrainMix.AllChromosome.Meth-nCpG.shared-sites.bed &
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SharedSites/CpG/mAbel-BrainMix.$c.CpG.shared-sites.bed.gz >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/TEMP/mAbel-BrainMix.AllChromosome.CpG.shared-sites.bed &
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/SharedSites/nCpG/mAbel-BrainMix.$c.nCpG.shared-sites.bed.gz >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/TEMP/mAbel-BrainMix.AllChromosome.nCpG.shared-sites.bed &
wait
done

shuf -n 30000 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/TEMP/mAbel-BrainMix.AllChromosome.Meth-CpG.shared-sites.bed  | sort -T /work/mjuswilc/tmp -k1,1V -k2,2n | gzip -9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/Random/mAbel-BrainMix.AllChromosome.Meth-CpG.shared-sites30K.bed.gz &
shuf -n 30000 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/TEMP/mAbel-BrainMix.AllChromosome.Meth-nCpG.shared-sites.bed  | sort -T /work/mjuswilc/tmp -k1,1V -k2,2n | gzip -9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/Random/mAbel-BrainMix.AllChromosome.Meth-nCpG.shared-sites30K.bed.gz &
shuf -n 30000 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/TEMP/mAbel-BrainMix.AllChromosome.CpG.shared-sites.bed  | sort -T /work/mjuswilc/tmp -k1,1V -k2,2n | gzip -9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/Random/mAbel-BrainMix.AllChromosome.CpG.shared-sites30K.bed.gz &
shuf -n 30000 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/TEMP/mAbel-BrainMix.AllChromosome.nCpG.shared-sites.bed  | sort -T /work/mjuswilc/tmp -k1,1V -k2,2n | gzip -9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/Random/mAbel-BrainMix.AllChromosome.nCpG.shared-sites30K.bed.gz &
wait

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/DeconvolutedBrain
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/DeconvolutedBrain
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/DeconvolutedBrain/CpG
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/DeconvolutedBrain/nCpG
for c in $CHROMOSOMES
do
for P in 0 10 20 30 40 50 60 70 80 90 100
do
bedtools intersect -sorted -f 1 -F 1 -a \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz -b \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/Random/mAbel-BrainMix.AllChromosome.Meth-CpG.shared-sites30K.bed.gz \
>> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/DeconvolutedBrain/CpG/mAbel-BrainMix.Blood-$P.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-CpG.shared-site30K.bed &
bedtools intersect -sorted -f 1 -F 1 -a \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz -b \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/Random/mAbel-BrainMix.AllChromosome.Meth-nCpG.shared-sites30K.bed.gz \
>> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/DeconvolutedBrain/nCpG/mAbel-BrainMix.Blood-$P.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.Meth-nCpG.shared-site30K.bed &
bedtools intersect -sorted -f 1 -F 1 -a \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz -b \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/Random/mAbel-BrainMix.AllChromosome.CpG.shared-sites30K.bed.gz \
>> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/DeconvolutedBrain/CpG/mAbel-BrainMix.Blood-$P.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.CpG.shared-site30K.bed &
bedtools intersect -sorted -f 1 -F 1 -a \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz -b \
/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/Random/mAbel-BrainMix.AllChromosome.nCpG.shared-sites30K.bed.gz \
>> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/DeconvolutedBrain/nCpG/mAbel-BrainMix.Blood-$P.AllChromosome.deduplicated.bismark.deconvoluted-brain-cov.nCpG.shared-site30K.bed &
wait
done
done
wait

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Subset/WholeGenome/TEMP

echo 'Subsetting Job Complete!!!'
