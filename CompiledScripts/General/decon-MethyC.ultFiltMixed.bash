#!/bin/bash
#SBATCH --time=1:48:00
#SBATCH --mem=6GB
#SBATCH -p short
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=18
# Output and error files
#SBATCH -o job.%J.out

P=$1

source activate Genomics

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylDecon/Chromosomes/$P
rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylDecon/DeconBrain/$P
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylDecon/Chromosomes/$P
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylDecon/DeconBrain/$P


CHROMOSOMES=`awk '{print $1}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/MethylationBed/Genome/Abel-GCF_001522545.3_Parus_major1.1_genomic.chrom.genome`
for c in $CHROMOSOMES
do
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/MethylationBed/Bed/mAbel-BrainMix.Blood-$P.deduplicated.bismark.cov-sorted.bed.gz | grep "^$c" | gzip -9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylDecon/Chromosomes/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.cov-sorted.bed.gz &
JOBS=`jobs -l | grep -c '^\['`
while [ $JOBS -ge 6 ]
do
sleep 3s
JOBS=`jobs -l | grep -c '^\['`
done
done
wait


BLOODBEDS=`ls -laht /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylDecon/Chromosomes/100/mAbel-BrainMix.Blood-100.*.deduplicated.bismark.cov-sorted.bed.gz | wc -l`
if [ $BLOODBEDS -lt 32 ]
then
echo "Waiting for Blood Bed Files to Complete."
fi
while [ $BLOODBEDS -lt 32 ]
do
sleep 3s
BLOODBEDS=`ls -laht /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylDecon/Chromosomes/100/mAbel-BrainMix.Blood-100.*.deduplicated.bismark.cov-sorted.bed.gz | wc -l`
if [ $BLOODBEDS -eq 32 ]
then
echo "Blood Bed Files have Completed."
fi
done


#Global P
for c in $CHROMOSOMES
do
bloodP=`awk -v p=$P 'BEGIN{OFS=FS="\t"}NR==1{next}$1==p{print $5}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/Output/blood-deconvolution-validation-table.txt`
tbloodP=`echo $P | awk '{print $0/100}'`
bedtools intersect -sorted -f 1 -F 1 -wo -a /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylDecon/Chromosomes/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.cov-sorted.bed.gz -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylDecon/Chromosomes/100/mAbel-BrainMix.Blood-100.$c.deduplicated.bismark.cov-sorted.bed.gz | awk -v PREC=256 -v p=$bloodP 'BEGIN{OFS=FS="\t"}{mDEPTH=($5+$6); pDEPTH=($11+$12); $5=$5-mDEPTH*($11/pDEPTH)*p; $6=$6-mDEPTH*($12/pDEPTH)*p}$5<0{$5=0}$6<0{$6=0}{$4=0}($5+$6)>=0.5{$4=$5/($5+$6)}{print $1 OFS $2 OFS $3 OFS $4 OFS int($5+0.5) OFS int($6+0.5)}' | gzip -9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylDecon/DeconBrain/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.bed.gz &
bedtools intersect -sorted -f 1 -F 1 -wo -a /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylDecon/Chromosomes/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.cov-sorted.bed.gz -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylDecon/Chromosomes/100/mAbel-BrainMix.Blood-100.$c.deduplicated.bismark.cov-sorted.bed.gz | awk -v PREC=256 -v p=$tbloodP 'BEGIN{OFS=FS="\t"}{mDEPTH=($5+$6); pDEPTH=($11+$12); $5=$5-mDEPTH*($11/pDEPTH)*p; $6=$6-mDEPTH*($12/pDEPTH)*p}$5<0{$5=0}$6<0{$6=0}{$4=0}($5+$6)>=0.5{$4=$5/($5+$6)}{print $1 OFS $2 OFS $3 OFS $4 OFS int($5+0.5) OFS int($6+0.5)}' | gzip -9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylDecon/DeconBrain/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.true-P.bed.gz &
JOBS=`jobs -l | grep -c '^\['`
while [ $JOBS -ge 6 ]
do
sleep 3s
JOBS=`jobs -l | grep -c '^\['`
done
done
wait

echo 'Deconvolution Job Complete!!!'
