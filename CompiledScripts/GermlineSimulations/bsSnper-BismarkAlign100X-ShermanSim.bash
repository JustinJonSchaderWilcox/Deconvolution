#!/bin/bash
#SBATCH --time=7:48:00
#SBATCH --mem=23GB
#SBATCH -p med
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=2
# Output and error files
#SBATCH -o job.%J.out

source activate Genomics

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SNPs
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SNPs
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SNPs/Reference

zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/BisulphiteGenome/Bismark/Abel-GCF_001522545.3_Parus_major1.1_genomic.chrom.fa.gz > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SNPs/Reference/Abel-GCF_001522545.3_Parus_major1.1_genomic.chrom.fasta

perl /home/mjuswilc/Programs/BS-Snper-master/BS-Snper.pl --fa /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SNPs/Reference/Abel-GCF_001522545.3_Parus_major1.1_genomic.chrom.fasta /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/BAMs/Sort/Bismark-Bowtie.sbS100X.deduplicated-sorted.bam --output /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SNPs/sbS100X.snp --methcg /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SNPs/sbS100X.CG --methchg /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SNPs/sbS100X.CHG --methchh /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SNPs/sbS100X.CHH --minhetfreq  0.01 --minhomfreq  0.85 --minquali 10 --mincover 1 --maxcover 100 --minread2 1 --errorate 0.02 --mapvalue 20 >  /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SNPs/sbS100X-SNV.vcf

bgzip --threads 1 -l 9 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SNPs/sbS100X-SNV.vcf
bcftools index /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SNPs/sbS100X-SNV.vcf.gz

echo 'Sherman Sim SNV Calling Job Complete!!!'
