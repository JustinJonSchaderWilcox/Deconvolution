#!/bin/bash
#SBATCH --time=7:00
#SBATCH --mem=4GB
#SBATCH -p short
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

rm -r /work/mjuswilc/Domain/GreatTits/General/RawData/Abel-Reference
mkdir /work/mjuswilc/Domain/GreatTits/General/RawData/Abel-Reference

rsync --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/522/545/GCF_001522545.3_Parus_major1.1/ /work/mjuswilc/Domain/GreatTits/General/RawData/Abel-Reference
echo 'Job Complete'
