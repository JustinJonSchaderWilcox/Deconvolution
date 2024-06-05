#!/bin/bash
#SBATCH --time=12:00
#SBATCH --mem=4GB
#SBATCH -p short
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=3
# Output and error files
#SBATCH -o job.%J.out

source activate Genomics

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/AlleleFreq
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/AlleleFreq
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/AlleleFreq/Summary

echo -e "Chromosome\tPosition\tRefAllele\tAltAllele\tMutation\tPhred\tDepth\tAlleleFreq" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/AlleleFreq/Summary/sbS100X-SNV-AlleleFreq.txt
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SNPs/sbS100X-SNV.vcf.gz | grep -v "#" | perl -pe 's#(\S+)\s+(\S+)\s+\.\s+(\S+)\s+(\S+)\s+(\S+).*DP=(\d+).*(\d\.\d+)$#$1\t$2\t$3\t$4\t$3$4\t$5\t$6\t$7#' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/AlleleFreq/Summary/sbS100X-SNV-AlleleFreq.txt

echo -e "Chromosome\tPosition\tRefAllele\tAltAllele\tMutation\tPhred\tDepth\tAlleleFreq" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/AlleleFreq/Summary/sbS100X-SNV-TrueHetAlleleFreq.txt
bedtools intersect -f 1 -F 1 -wa -a /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SNPs/sbS100X-SNV.vcf.gz -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/SimulatedGenome/Haplotype2/great-tit-sim.Abel-GCF_001522545v3-Hap2_ms.vcf | grep -v "#" | perl -pe 's#(\S+)\s+(\S+)\s+\.\s+(\S+)\s+(\S+)\s+(\S+).*DP=(\d+).*(\d\.\d+)$#$1\t$2\t$3\t$4\t$3$4\t$5\t$6\t$7#' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/AlleleFreq/Summary/sbS100X-SNV-TrueHetAlleleFreq.txt

echo -e "Chromosome\tPosition\tRefAllele\tAltAllele\tMutation\tPhred\tDepth\tAlleleFreq" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/AlleleFreq/Summary/sbS100X-SNV-ErrorAlleleFreq.txt
bedtools subtract -A -a /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/SNPs/sbS100X-SNV.vcf.gz -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/SimulatedGenome/Haplotype2/great-tit-sim.Abel-GCF_001522545v3-Hap2_ms.vcf | grep -v "#" | perl -pe 's#(\S+)\s+(\S+)\s+\.\s+(\S+)\s+(\S+)\s+(\S+).*DP=(\d+).*(\d\.\d+)$#$1\t$2\t$3\t$4\t$3$4\t$5\t$6\t$7#' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/AlleleFreq/Summary/sbS100X-SNV-ErrorAlleleFreq.txt

echo 'Allele Frequecny Spectra Extracted!!!'
