#!/bin/bash
#SBATCH --time=1:46:00
#SBATCH --mem=8GB
#SBATCH -p short
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

source activate EMBOSS

#Annotate CpG Islands
rm -r /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/CpG-Islands
mkdir /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/CpG-Islands

gzip -d -c /work/mjuswilc/Domain/GreatTits/General/RawData/Abel-Reference/GCF_001522545.3_Parus_major1.1_genomic.fna.gz > /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/CpG-Islands/GCF_001522545.3_Parus_major1.1_genomic.fna
cpgplot /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/CpG-Islands/GCF_001522545.3_Parus_major1.1_genomic.fna -window 100 -minlen 200 -minoe 0.6 -minpc 50 -outfile /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/CpG-Islands/GCF_001522545.3_Parus_major1.1_genomic-CpG.svg -outfeat /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/CpG-Islands/GCF_001522545.3_Parus_major1.1_genomic-CpG.gff -graph svg
'echo Job Complete'
