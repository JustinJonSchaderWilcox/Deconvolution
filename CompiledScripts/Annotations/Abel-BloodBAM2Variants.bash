#!/bin/bash
#SBATCH --time=1-23:46:00
#SBATCH --mem=18GB
#SBATCH -p long
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=12
# Output and error files
#SBATCH -o job.%J.out

source activate Genomics

rm -r /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP
mkdir /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP
mkdir /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/Inputs
mkdir /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs
mkdir /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/AllSite
mkdir /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/SNVs
mkdir /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/Indels

gzip -d -c /work/mjuswilc/Domain/GreatTits/General/RawData/Abel-Reference/GCF_001522545.3_Parus_major1.1_genomic.fna.gz > /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/Inputs/GCF_001522545.3_Parus_major1.1_genomic.fna
samtools faidx /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/Inputs/GCF_001522545.3_Parus_major1.1_genomic.fna

bcftools mpileup -f /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/Inputs/GCF_001522545.3_Parus_major1.1_genomic.fna -d 2300 /work/mjuswilc/Domain/GreatTits/General/RawData/BAM/Abel/Blood-SRR3085443_recal_reads.bam | bcftools call -c - | bcftools norm -f /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/Inputs/GCF_001522545.3_Parus_major1.1_genomic.fna  | bgzip -l 9 --threads 8 > /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/AllSite/Abel-GCF_001522545.3_Parus_major1.1_genomic.all-site.vcf.gz

#Extract Heterozygous Sites from all site VCF: Q30 candidates, Q90 certain
bgzip -d -c /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/AllSite/Abel-GCF_001522545.3_Parus_major1.1_genomic.all-site.vcf.gz | awk 'BEGIN{OFS=FS="\t"}/^#/{print $0; next}length($4)>1 || length($5)>1{next}/\s0\/1/{print $0}' | bgzip -l 9 --threads 2 > /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/SNVs/Abel-GCF_001522545.3_Parus_major1.1_genomic.snvs-any.vcf.gz &
bgzip -d -c /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/AllSite/Abel-GCF_001522545.3_Parus_major1.1_genomic.all-site.vcf.gz | awk 'BEGIN{OFS=FS="\t"}/^#/{print $0; next}length($4)>1 || length($5)>1{next}$6>=20 && /\s0\/1/{print $0}' | bgzip -l 9 --threads 2 > /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/SNVs/Abel-GCF_001522545.3_Parus_major1.1_genomic.snvs-Q20.vcf.gz &
bgzip -d -c /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/AllSite/Abel-GCF_001522545.3_Parus_major1.1_genomic.all-site.vcf.gz | awk 'BEGIN{OFS=FS="\t"}/^#/{print $0; next}length($4)>1 || length($5)>1{next}$6>=30 && /\s0\/1/{print $0}' | bgzip -l 9 --threads 2 > /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/SNVs/Abel-GCF_001522545.3_Parus_major1.1_genomic.snvs-Q30.vcf.gz &
bgzip -d -c /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/AllSite/Abel-GCF_001522545.3_Parus_major1.1_genomic.all-site.vcf.gz | awk 'BEGIN{OFS=FS="\t"}/^#/{print $0; next}length($4)>1 || length($5)>1{next}$6>=90 && /\s0\/1/{print $0}' | bgzip -l 9 --threads 2 > /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/SNVs/Abel-GCF_001522545.3_Parus_major1.1_genomic.snvs-Q90.vcf.gz &
wait

#Extract Indels Sites from all site VCF: Q30 candidates, Q90 high confidence
bgzip -d -c /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/AllSite/Abel-GCF_001522545.3_Parus_major1.1_genomic.all-site.vcf.gz | awk 'BEGIN{OFS=FS="\t"}/^#/{print $0; next}length($4)==1 || length($5)==1{next}/\s0\/1/{print $0}' | bgzip -l 9 --threads 2 > /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/Indels/Abel-GCF_001522545.3_Parus_major1.1_genomic.Indels-any.vcf.gz &
bgzip -d -c /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/AllSite/Abel-GCF_001522545.3_Parus_major1.1_genomic.all-site.vcf.gz | awk 'BEGIN{OFS=FS="\t"}/^#/{print $0; next}length($4)==1 || length($5)==1{next}$6>=20 && /\s0\/1/{print $0}' | bgzip -l 9 --threads 2 > /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/Indels/Abel-GCF_001522545.3_Parus_major1.1_genomic.Indels-Q20.vcf.gz &
bgzip -d -c /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/AllSite/Abel-GCF_001522545.3_Parus_major1.1_genomic.all-site.vcf.gz | awk 'BEGIN{OFS=FS="\t"}/^#/{print $0; next}length($4)==1 || length($5)==1{next}$6>=30 && /\s0\/1/{print $0}' | bgzip -l 9 --threads 2 > /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/Indels/Abel-GCF_001522545.3_Parus_major1.1_genomic.Indels-Q30.vcf.gz &
bgzip -d -c /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/AllSite/Abel-GCF_001522545.3_Parus_major1.1_genomic.all-site.vcf.gz | awk 'BEGIN{OFS=FS="\t"}/^#/{print $0; next}length($4)==1 || length($5)==1{next}$6>=90 && /\s0\/1/{print $0}' | bgzip -l 9 --threads 2 > /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/Indels/Abel-GCF_001522545.3_Parus_major1.1_genomic.Indels-Q90.vcf.gz &
wait

bcftools index /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/AllSite/Abel-GCF_001522545.3_Parus_major1.1_genomic.all-site.vcf.gz &
bcftools index /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/SNVs/Abel-GCF_001522545.3_Parus_major1.1_genomic.snvs-any.vcf.gz &
bcftools index /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/SNVs/Abel-GCF_001522545.3_Parus_major1.1_genomic.snvs-Q20.vcf.gz &
bcftools index /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/SNVs/Abel-GCF_001522545.3_Parus_major1.1_genomic.snvs-Q30.vcf.gz &
bcftools index /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/SNVs/Abel-GCF_001522545.3_Parus_major1.1_genomic.snvs-Q90.vcf.gz &
bcftools index /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/Indels/Abel-GCF_001522545.3_Parus_major1.1_genomic.Indels-any.vcf.gz &
bcftools index /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/Indels/Abel-GCF_001522545.3_Parus_major1.1_genomic.Indels-Q20.vcf.gz &
bcftools index /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/Indels/Abel-GCF_001522545.3_Parus_major1.1_genomic.Indels-Q30.vcf.gz &
bcftools index /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/Indels/Abel-GCF_001522545.3_Parus_major1.1_genomic.Indels-Q90.vcf.gz &
wait

echo 'Abel Blood BAM MPILEUP Complete!!!'
