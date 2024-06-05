#!/bin/bash
#SBATCH --time=1-23:00:00
#SBATCH --mem=18GB
#SBATCH -p long
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=3
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out


rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SubsampledReads
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SubsampledReads

/home/mjuswilc/Programs/seqtk-master/seqtk sample -s 23 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/MergedReads/merged-ShermanSimulation-Diploid_1.fastq.gz 0.25 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SubsampledReads/sbS100X-ShermanSimulation-Diploid_1.fastq &
/home/mjuswilc/Programs/seqtk-master/seqtk sample -s 23 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/MergedReads/merged-ShermanSimulation-Diploid_2.fastq.gz 0.25 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SubsampledReads/sbS100X-ShermanSimulation-Diploid_2.fastq &
wait

gzip -9 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SubsampledReads/sbS100X-ShermanSimulation-Diploid_1.fastq &
gzip -9 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SubsampledReads/sbS100X-ShermanSimulation-Diploid_2.fastq &
wait

echo 'Sherman Simulation Subsample Job Complete!!!'
