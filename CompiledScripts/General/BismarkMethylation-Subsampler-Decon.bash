#!/bin/bash
#SBATCH --time=23:40:00
#SBATCH --mem=58GB
#SBATCH -p long
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=20
#SBATCH --exclusive
#SBATCH --constraint=cstd01
# Output and error files
#SBATCH -o job.%J.out

source activate Genomics

P=$1
t=$2

sleep $t

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Methylation/$P
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Methylation/$P

function bismark_methylation {
rm -r /scratch/mjuswilc
mkdir /scratch/mjuswilc
mkdir /scratch/mjuswilc/BismarkMethylation
mkdir /scratch/mjuswilc/BismarkMethylation/$P
cp /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Dedup/$P/mAbel-BrainMix.Blood-$P.deduplicated.bam /scratch/mjuswilc/BismarkMethylation/$P/mAbel-BrainMix.Blood-$P.deduplicated.bam
/home/mjuswilc/Programs/Bismark-0.22.3/bismark_methylation_extractor -p --comprehensive --include_overlap --ignore 2 --ignore_r2 2 --ignore_3prime 2 --ignore_3prime_r2 2 --cutoff 1 --parallel 6 --bedGraph --CX --ample_memory --gzip -o /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Methylation/$P /scratch/mjuswilc/BismarkMethylation/$P/mAbel-BrainMix.Blood-$P.deduplicated.bam
rm -r /scratch/mjuswilc
}

function clear_scratch {
rm -r /scratch/mjuswilc
}

trap -- 'clear_scratch' USR1 SIGQUIT SIGINT SIGABRT SIGKILL SIGHUP SIGTERM EXIT TERM
bismark_methylation &> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Methylation/$P/bismark-methylation.$P.log &
wait
trap - USR1 SIGQUIT SIGINT SIGABRT SIGKILL SIGHUP SIGTERM EXIT TERM

echo 'Bismark Methylation Extraction Job Complete!!!'
