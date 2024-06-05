#!/bin/bash
#SBATCH --time=4-23:56:00
#SBATCH --mem=4GB
#SBATCH -p ultralong
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=3
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

rm /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/MergedReads/*.fastq.gz
gzip -9 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/MergedReads/merged-ShermanSimulation-Diploid_1.fastq &
gzip -9 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/MergedReads/merged-ShermanSimulation-Diploid_2.fastq &
wait

echo 'Sherman Simulation Merged Read GZIP Complete!!!'
