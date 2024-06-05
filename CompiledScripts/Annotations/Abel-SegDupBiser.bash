#!/bin/bash
#SBATCH --time=34:00
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

rm -r /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/SegmentalDuplications/BISER
mkdir /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/SegmentalDuplications/BISER
mkdir /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/SegmentalDuplications/BISER/Input
mkdir /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/SegmentalDuplications/BISER/Output

gzip -d -c /work/mjuswilc/Domain/GreatTits/General/RawData/Abel-Reference/GCF_001522545.3_Parus_major1.1_genomic.fna.gz > /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/SegmentalDuplications/BISER/Input/GCF_001522545.3_Parus_major1.1_genomic.fna
samtools faidx /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/SegmentalDuplications/BISER/Input/GCF_001522545.3_Parus_major1.1_genomic.fna

biser --threads 20 --keep-contigs --keep-temp --gc-heap 1G -o /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/SegmentalDuplications/BISER/Output/Abel-GCF_001522545.3_Parus_major1.1_genomic.fna.seg_dup.bedpe /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/SegmentalDuplications/BISER/Input/GCF_001522545.3_Parus_major1.1_genomic.fna
echo 'Job Complete'
