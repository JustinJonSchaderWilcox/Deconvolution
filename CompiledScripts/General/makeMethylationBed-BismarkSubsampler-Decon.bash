#!/bin/bash
#SBATCH --time=5:46:00
#SBATCH --mem=3GB
#SBATCH -p med
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=4
# Output and error files
#SBATCH -o job.%J.out

P=$1

rm -r /scratch/mjuswilc/$P
mkdir -p /scratch/mjuswilc/$P

function clear_scratch {
rm -r /scratch/mjuswilc/$P
}


function ChromosomalConfiguration {
for c in $(awk 'BEGIN{OFS=FS="\t"}{print $1}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/MethylationBed/Genome/Abel-GCF_001522545.3_Parus_major1.1_genomic.chrom.genome)
do
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Methylation/$P/mAbel-BrainMix.Blood-$P.deduplicated.bismark.cov.gz | grep $c | sort -T /scratch/mjuswilc/$P  -n -k 2 | gzip -9 >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/MethylationBed/Bed/mAbel-BrainMix.Blood-$P.deduplicated.bismark.cov-sorted.bed.gz
done
}

trap -- 'clear_scratch' USR1 SIGQUIT SIGINT SIGABRT SIGKILL SIGHUP SIGTERM EXIT TERM
ChromosomalConfiguration &> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/MethylationBed/Log/MakeMethylationBed.$P.log &
wait
rm -r /scratch/mjuswilc/$P
trap - USR1 SIGQUIT SIGINT SIGABRT SIGKILL SIGHUP SIGTERM EXIT TERM

echo 'Methylation Bed Job Complete!!!'
