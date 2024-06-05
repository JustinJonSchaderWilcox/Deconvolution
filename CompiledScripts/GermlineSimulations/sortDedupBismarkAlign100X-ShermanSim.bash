#!/bin/bash
#SBATCH --time=23:40:00
#SBATCH --mem=74GB
#SBATCH -p long
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

source activate Genomics

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/BAMs
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/BAMs
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/BAMs/Dedup
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/BAMs/Sort

#Run Bismark Deduplicate
/home/mjuswilc/Programs/Bismark-0.22.3/deduplicate_bismark -p /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/Bismark/Alignment/Bismark-Bowtie.sbS100X-ShermanSimulation-Diploid_1.paired_bismark_bt2_pe.bam  -o Bismark-Bowtie.sbS100X --output_dir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/BAMs/Dedup

samtools sort -m 28G /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/BAMs/Dedup/Bismark-Bowtie.sbS100X.deduplicated.bam -o /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/BAMs/Sort/Bismark-Bowtie.sbS100X.deduplicated-sorted.bam
echo 'Heterozygote Simulation Deduplication and Sorting Job Complete!!!'
