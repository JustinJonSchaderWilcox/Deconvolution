#!/bin/bash
#SBATCH --time=23:00:00
#SBATCH --mem=12GB
#SBATCH -p long
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=5
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

#Download Genome w/annotations and the SRA Fastq Files
rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/RawData/Sequences
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/RawData/Sequences
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/RawData/Sequences/Bisulphite
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/RawData/Sequences/Bisulphite/Brain
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/RawData/Sequences/Bisulphite/Blood
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/RawData/Sequences/Bisulphite/Brain/Fastq
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/RawData/Sequences/Bisulphite/Blood/Fastq

rm -r /scratch/mjuswilc
mkdir /scratch/mjuswilc
mkdir /scratch/mjuswilc/SRA
mkdir /scratch/mjuswilc/SRA_Cache
mkdir /scratch/mjuswilc/BrainFastq
mkdir /scratch/mjuswilc/BloodFastq


/home/mjuswilc/Programs/sratoolkit.3.0.5-ubuntu64/bin/prefetch -X 100g https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR2070790/SRR2070790 --output-directory /scratch/mjuswilc/SRA &
/home/mjuswilc/Programs/sratoolkit.3.0.5-ubuntu64/bin/prefetch -X 100g https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR2070791/SRR2070791 --output-directory /scratch/mjuswilc/SRA &
wait

/home/mjuswilc/Programs/sratoolkit.3.0.5-ubuntu64/bin/fasterq-dump /scratch/mjuswilc/SRA/SRR2070790/SRR2070790.sra -S -O /scratch/mjuswilc/BrainFastq -t /scratch/mjuswilc/SRA_Cache --threads 2 &
/home/mjuswilc/Programs/sratoolkit.3.0.5-ubuntu64/bin/fasterq-dump /scratch/mjuswilc/SRA/SRR2070791/SRR2070791.sra -S -O /scratch/mjuswilc/BloodFastq -t /scratch/mjuswilc/SRA_Cache --threads 2 &
wait

gzip -9 /scratch/mjuswilc/BrainFastq/SRR2070790_1.fastq &
gzip -9 /scratch/mjuswilc/BrainFastq/SRR2070790_2.fastq &
gzip -9 /scratch/mjuswilc/BloodFastq/SRR2070791_1.fastq &
gzip -9 /scratch/mjuswilc/BloodFastq/SRR2070791_2.fastq &
wait
cp /scratch/mjuswilc/BrainFastq/SRR2070790_1.fastq.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/RawData/Sequences/Bisulphite/Brain/Fastq/SRR2070790_1.fastq.gz &
cp /scratch/mjuswilc/BrainFastq/SRR2070790_2.fastq.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/RawData/Sequences/Bisulphite/Brain/Fastq/SRR2070790_2.fastq.gz &
cp /scratch/mjuswilc/BloodFastq/SRR2070791_1.fastq.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/RawData/Sequences/Bisulphite/Blood/Fastq/SRR2070791_1.fastq.gz &
cp /scratch/mjuswilc/BloodFastq/SRR2070791_2.fastq.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/RawData/Sequences/Bisulphite/Blood/Fastq/SRR2070791_2.fastq.gz &
wait
echo 'Job Complete!!!'
