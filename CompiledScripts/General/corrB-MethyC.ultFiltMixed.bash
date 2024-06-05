#!/bin/bash
#SBATCH --time=7:46:00
#SBATCH --mem=58GB
#SBATCH -p med
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=20
#SBATCH --constraint=cstd01
#SBATCH --exclusive
# Output and error files
#SBATCH -o job.%J.out

P=$1
t=$2

sleep "$t"m

#Paths
source activate Genomics

deconJOBS=`squeue -u mjuswilc | grep -P -c "\s+decon-Me"`
if [ $deconJOBS -gt 0 ]
then
echo "Waiting for Deconvolution Jobs to Complete."
fi
while [ $deconJOBS -gt 0 ]
do
sleep 3s
deconJOBS=`squeue -u mjuswilc | grep -P -c "\s+decon-Me"`
if [ $deconJOBS -eq 0 ]
then
echo "Deconvolution Jobs have Completed."
fi
done

rm -r /scratch/mjuswilc
mkdir /scratch/mjuswilc
mkdir /scratch/mjuswilc/DeconBrain
cp -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylDecon/Chromosomes /scratch/mjuswilc/Chromosomes
cp -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylDecon/DeconBrain/$P /scratch/mjuswilc/DeconBrain/$P

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Inputs/Beds/$P
rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/All-C/$P
rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/$P
rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/$P

mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Inputs/Beds/$P
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Inputs/Beds/$P/CpG
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Inputs/Beds/$P/nCpG
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/All-C/$P
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/$P
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/$P

function clear_scratch {
rm -r /scratch/mjuswilc
}

function ChromosomalCorrelation {
#Creating Input Context Information
CHROMOSOMES=`awk '{print $1}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/MethylationBed/Genome/Abel-GCF_001522545.3_Parus_major1.1_genomic.chrom.genome`
for c in $CHROMOSOMES
do
grep "^$c" /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.CG  | awk 'BEGIN{OFS=FS="\t"}{print $1 OFS $2 OFS $2 OFS "CpG"}' | gzip -9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Inputs/Beds/$P/CpG/BloodBrainMix.$P-$c.CpG.bed.gz &
cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.CHG /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.CHH | grep "^$c" | awk 'BEGIN{OFS=FS="\t"}{print $1 OFS $2 OFS $2 OFS "nCpG"}' | sort -T /scratch/mjuswilc -n -k 2 | gzip -9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Inputs/Beds/$P/nCpG/BloodBrainMix.$P-$c.nCpG.bed.gz &
JOBS=`jobs -l | grep -c '^\['`
while [ $JOBS -ge 4 ]
do
sleep 3s
JOBS=`jobs -l | grep -c '^\['`
done
done
wait

#Order of deconvolution bed files: deconvolution, brain, blood,convoluted mix, deconvolution with true P, global-P deconvolution, global-P true-P; proportion methylated, methylated reads unmethylated reads
for c in $CHROMOSOMES
do
bedtools intersect -sorted -f 1 -F 1 -wo -a /scratch/mjuswilc/DeconBrain/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.bed.gz -b /scratch/mjuswilc/Chromosomes/0/mAbel-BrainMix.Blood-0.$c.deduplicated.bismark.cov-sorted.bed.gz | awk 'BEGIN{OFS=FS="\t"}{print $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $10/100 OFS $11 OFS $12}' | \
bedtools intersect -sorted -f 1 -F 1 -wo -a - -b /scratch/mjuswilc/Chromosomes/100/mAbel-BrainMix.Blood-100.$c.deduplicated.bismark.cov-sorted.bed.gz | awk 'BEGIN{OFS=FS="\t"}{print $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $13/100 OFS $14 OFS $15}' | \
bedtools intersect -sorted -f 1 -F 1 -wo -a - -b /scratch/mjuswilc/Chromosomes/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.cov-sorted.bed.gz | awk 'BEGIN{OFS=FS="\t"}{print $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11 OFS $12 OFS $16/100 OFS $17 OFS $18}' | \
bedtools intersect -sorted -f 1 -F 1 -wo -a - -b /scratch/mjuswilc/DeconBrain/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.true-P.bed.gz | awk -v p=$P 'BEGIN{OFS=FS="\t"}{print $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11 OFS $12 OFS $13 OFS $14 OFS $15 OFS $19 OFS $20 OFS $21 OFS p}' | gzip -9 \
> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/All-C/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.all-c.bed.gz &
JOBS=`jobs -l | grep -c '^\['`
while [ $JOBS -ge 2 ]
do
sleep 3s
JOBS=`jobs -l | grep -c '^\['`
done
done

for c in $CHROMOSOMES
do
bedtools intersect -sorted -f 1 -F 1 -wo -a /scratch/mjuswilc/DeconBrain/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.bed.gz -b /scratch/mjuswilc/Chromosomes/0/mAbel-BrainMix.Blood-0.$c.deduplicated.bismark.cov-sorted.bed.gz | awk 'BEGIN{OFS=FS="\t"}{print $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $10/100 OFS $11 OFS $12}' | \
bedtools intersect -sorted -f 1 -F 1 -wo -a - -b /scratch/mjuswilc/Chromosomes/100/mAbel-BrainMix.Blood-100.$c.deduplicated.bismark.cov-sorted.bed.gz | awk 'BEGIN{OFS=FS="\t"}{print $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $13/100 OFS $14 OFS $15}' | \
bedtools intersect -sorted -f 1 -F 1 -wo -a - -b /scratch/mjuswilc/Chromosomes/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.cov-sorted.bed.gz | awk 'BEGIN{OFS=FS="\t"}{print $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11 OFS $12 OFS $16/100 OFS $17 OFS $18}' | \
bedtools intersect -sorted -f 1 -F 1 -wo -a - -b /scratch/mjuswilc/DeconBrain/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.true-P.bed.gz | awk -v p=$P 'BEGIN{OFS=FS="\t"}{print $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11 OFS $12 OFS $13 OFS $14 OFS $15 OFS $19 OFS $20 OFS $21 OFS p}' | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Inputs/Beds/$P/CpG/BloodBrainMix.$P-$c.CpG.bed.gz | gzip -9 \
> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/CpG/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.CpG.bed.gz &
JOBS=`jobs -l | grep -c '^\['`
while [ $JOBS -ge 2 ]
do
sleep 3s
JOBS=`jobs -l | grep -c '^\['`
done
done

for c in $CHROMOSOMES
do
bedtools intersect -sorted -f 1 -F 1 -wo -a /scratch/mjuswilc/DeconBrain/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.bed.gz -b /scratch/mjuswilc/Chromosomes/0/mAbel-BrainMix.Blood-0.$c.deduplicated.bismark.cov-sorted.bed.gz | awk 'BEGIN{OFS=FS="\t"}{print $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $10/100 OFS $11 OFS $12}' | \
bedtools intersect -sorted -f 1 -F 1 -wo -a - -b /scratch/mjuswilc/Chromosomes/100/mAbel-BrainMix.Blood-100.$c.deduplicated.bismark.cov-sorted.bed.gz | awk 'BEGIN{OFS=FS="\t"}{print $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $13/100 OFS $14 OFS $15}' | \
bedtools intersect -sorted -f 1 -F 1 -wo -a - -b /scratch/mjuswilc/Chromosomes/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.cov-sorted.bed.gz | awk 'BEGIN{OFS=FS="\t"}{print $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11 OFS $12 OFS $16/100 OFS $17 OFS $18}' | \
bedtools intersect -sorted -f 1 -F 1 -wo -a - -b /scratch/mjuswilc/DeconBrain/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.true-P.bed.gz | awk -v p=$P 'BEGIN{OFS=FS="\t"}{print $1 OFS $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9 OFS $10 OFS $11 OFS $12 OFS $13 OFS $14 OFS $15 OFS $19 OFS $20 OFS $21 OFS p}' | \
bedtools intersect -sorted -f 1 -F 1 -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Inputs/Beds/$P/nCpG/BloodBrainMix.$P-$c.nCpG.bed.gz | gzip -9 \
> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Outputs/DeconvolutedBrain/nCpG/$P/mAbel-BrainMix.Blood-$P.$c.deduplicated.bismark.deconvoluted-brain-cov.nCpG.bed.gz &
JOBS=`jobs -l | grep -c '^\['`
while [ $JOBS -ge 2 ]
do
sleep 3s
JOBS=`jobs -l | grep -c '^\['`
done
done
wait
}



trap -- 'clear_scratch' USR1 SIGQUIT SIGINT SIGABRT SIGKILL SIGHUP SIGTERM EXIT TERM
ChromosomalCorrelation &>  /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/MethylCorreleations/Logs/"$P".chromosomal-correlation.log &
wait
rm -r /scratch/mjuswilc
trap - USR1 SIGQUIT SIGINT SIGABRT SIGKILL SIGHUP SIGTERM EXIT TERM

echo 'Correleation Configuration Job Complete!!!'
