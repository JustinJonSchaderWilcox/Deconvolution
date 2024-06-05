#!/bin/bash
#SBATCH --time=4:58:00
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

rm /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Sort/mAbel-BrainMix.Blood-$P.deduplicated.sorted.bam
samtools sort /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Dedup/$P/mAbel-BrainMix.Blood-$P.deduplicated.bam -m 28G -o /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Sort/mAbel-BrainMix.Blood-$P.deduplicated.sorted.bam
echo 'Job Complete!!!'
