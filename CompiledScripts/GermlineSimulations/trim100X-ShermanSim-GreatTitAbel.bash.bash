#!/bin/bash
#SBATCH --time=6:48:00
#SBATCH --mem=26GB
#SBATCH -p med
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=16
# Output and error files
#SBATCH -o job.%J.out

module purge
module load java/1.8.0u201

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/TrimmedReads
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/TrimmedReads
cat /home/mjuswilc/Programs/Trimmomatic-0.39/adapters/*.fa > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/TrimmedReads/IlluminaAdapters.fasta

java -jar /home/mjuswilc/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 8 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SubsampledReads/sbS100X-ShermanSimulation-Diploid_1.fastq.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SubsampledReads/sbS100X-ShermanSimulation-Diploid_2.fastq.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/TrimmedReads/sbS100X-ShermanSimulation-Diploid_1.paired.fastq.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/TrimmedReads/sbS100X-ShermanSimulation-Diploid_1.fastq.unpaired.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/TrimmedReads/sbS100X-ShermanSimulation-Diploid_2.paired.fastq.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/TrimmedReads/sbS100X-ShermanSimulation-Diploid_2.unpaired.fastq.gz ILLUMINACLIP:/work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/TrimmedReads/IlluminaAdapters.fasta:2:26:8:4:TRUE LEADING:3 TRAILING:3

echo 'Trimmomatic Job Complete!!!'
