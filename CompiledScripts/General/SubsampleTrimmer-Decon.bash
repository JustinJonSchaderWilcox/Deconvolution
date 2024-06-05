#!/bin/bash
#SBATCH --time=7:40:00
#SBATCH --mem=12GB
#SBATCH -p med
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=4
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

P=$1

module purge
module load java/1.8.0u201

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/FastQs/Trimmed/$P
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/FastQs/Trimmed/$P
java -jar /home/mjuswilc/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 4 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/FastQs/$P/mAbel-BrainMix.Blood-$P.read1.fastq.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/FastQs/$P/mAbel-BrainMix.Blood-$P.read2.fastq.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/FastQs/Trimmed/$P/mAbel-BrainMix.Blood-$P.read1.trimmomatic.paired.fastq.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/FastQs/Trimmed/$P/mAbel-BrainMix.Blood-$P.read1.trimmomatic.unpaired.fastq.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/FastQs/Trimmed/$P/mAbel-BrainMix.Blood-$P.read2.trimmomatic.paired.fastq.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/FastQs/Trimmed/$P/mAbel-BrainMix.Blood-$P.read2.trimmomatic.unpaired.fastq.gz ILLUMINACLIP:/work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/FastQs/Trimmed/IlluminaAdapters.fasta:2:26:8:4:TRUE LEADING:3 TRAILING:3
echo 'Job Complete!!!'
