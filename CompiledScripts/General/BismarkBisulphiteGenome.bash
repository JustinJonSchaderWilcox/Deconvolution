#!/bin/bash
#SBATCH --time=1:40:00
#SBATCH --mem=64GB
#SBATCH -p short
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=20
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

source activate Genomics
export PATH="$PATH:/home/mjuswilc/Programs/bowtie2-2.5.1"

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/BisulphiteGenome
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/BisulphiteGenome
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/BisulphiteGenome/Bismark

zcat /work/mjuswilc/Domain/GreatTits/General/RawData/Abel-Reference/GCF_001522545.3_Parus_major1.1_genomic.fna.gz | awk 'BEGIN{count=0}/^>/{count=count+1}count>32{exit}{print $0}' > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/BisulphiteGenome/Bismark/Abel-GCF_001522545.3_Parus_major1.1_genomic.chrom.fa
gzip -9 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/BisulphiteGenome/Bismark/Abel-GCF_001522545.3_Parus_major1.1_genomic.chrom.fa

/home/mjuswilc/Programs/Bismark-0.22.3/bismark_genome_preparation --bowtie2 --parallel 10 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/BisulphiteGenome/Bismark
echo 'Job Complete!!!'
