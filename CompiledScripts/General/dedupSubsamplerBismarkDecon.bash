#!/bin/bash
#SBATCH --time=7:40:00
#SBATCH --mem=34GB
#SBATCH -p med
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

P=$1

source activate Genomics

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Dedup/$P
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Dedup/$P

#Run Bismark Deduplicate
/home/mjuswilc/Programs/Bismark-0.22.3/deduplicate_bismark -p /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/$P/Bismark-Bowtie.mAbel-BrainMix.Blood-"$P".read1.trimmomatic.paired_bismark_bt2_pe.bam -o mAbel-BrainMix.Blood-$P --output_dir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Dedup/$P

echo 'Job Complete!!!'
