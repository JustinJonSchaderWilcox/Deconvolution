#!/bin/bash
#SBATCH --time=5:40:00
#SBATCH --mem=168GB
#SBATCH -p med
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=12
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

source activate Genomics

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/Reference

zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/BisulphiteGenome/Bismark/Abel-GCF_001522545.3_Parus_major1.1_genomic.chrom.fa.gz > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/Reference/Abel-GCF_001522545.3_Parus_major1.1_genomic.chrom.fasta

for P in 0 10 20 30 40 50 60 70 80 90 100
do
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/"$P"
perl /home/mjuswilc/Programs/BS-Snper-master/BS-Snper.pl --fa /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/Reference/Abel-GCF_001522545.3_Parus_major1.1_genomic.chrom.fasta /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Sort/mAbel-BrainMix.Blood-$P.deduplicated.sorted.bam --output /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp --methcg /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.CG --methchg /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.CHG --methchh /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.CHH --minhetfreq  0.01 --minhomfreq  0.85 --minquali 10 --mincover 1 --maxcover 100 --minread2 1 --errorate 0.02 --mapvalue 20 >  /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf &
done
wait

function bgzipit {
bgzip --threads 1 -l 9 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf
bcftools index /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz
}
for P in 0 10 20 30 40 50 60 70 80 90 100
do
bgzipit &
done
wait

echo 'Job Complete!!!'
