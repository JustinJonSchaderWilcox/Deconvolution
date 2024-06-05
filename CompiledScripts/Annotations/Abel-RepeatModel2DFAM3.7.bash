#!/bin/bash
#SBATCH --time=1-23:00:00
#SBATCH --mem=256GB
#SBATCH -p long
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=48
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

module purge
source activate Genomics

rm -r /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatModeler2
mkdir /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatModeler2
mkdir /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatModeler2

#Make Database
cd /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatModeler2
zcat /work/mjuswilc/Domain/GreatTits/General/RawData/Abel-Reference/GCF_001522545.3_Parus_major1.1_genomic.fna.gz > /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatModeler2/Abel-GCF_001522545.3_Parus_major1.1_genomic.fna
#Build Database
/home/mjuswilc/Programs/RepeatModeler-2.0.4/BuildDatabase -engine ncbi -name Abel-GCF_001522545.3_Parus_major1.1_genomic /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatModeler2/Abel-GCF_001522545.3_Parus_major1.1_genomic.fna

##RepeatModeler
/home/mjuswilc/Programs/RepeatModeler-2.0.4/RepeatModeler -threads 48 -database Abel-GCF_001522545.3_Parus_major1.1_genomic -engine ncbi -genomeSampleSizeMax 243000000 -LTRStruct

rm -r /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatModeler2/RM_*
echo 'Job Complete'
