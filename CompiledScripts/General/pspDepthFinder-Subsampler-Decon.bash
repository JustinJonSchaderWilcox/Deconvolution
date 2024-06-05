#!/bin/bash
#SBATCH --time=3:46:00
#SBATCH --mem=8GB
#SBATCH -p med
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=24
# Output and error files
#SBATCH -o job.%J.out

source activate Genomics

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder


for P in 0 10 20 30 40 50 60 70 80 90 100
do
samtools depth -a /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Sort/mAbel-BrainMix.Blood-$P.deduplicated.sorted.bam | awk 'BEGIN{OFS=FS="\t"; bases=0; print "Chromosome" OFS "Basepairs" OFS "Depth"}NR==1{chrom=$1; depth=$3; bases=1; next}$1!=chrom{print chrom OFS bases OFS depth/bases; chrom=$1; depth=$3; bases=1; next}{depth=depth+$3; bases++}END{print chrom OFS bases OFS depth/bases}' >  /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_"$P".chrom-depth.txt &
done
wait

for P in 0 10 20 30 40 50 60 70 80 90 100
do
samtools depth -a /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Sort/mAbel-BrainMix.Blood-$P.deduplicated.sorted.bam | awk 'BEGIN{OFS=FS="\t"; depth=0}{depth=depth+$3}END{print depth/NR}' > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_"$P".depth.txt &
done
wait

echo 'Depth Finder Complete!!!'
