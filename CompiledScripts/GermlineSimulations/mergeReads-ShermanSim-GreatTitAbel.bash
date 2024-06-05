#!/bin/bash
#SBATCH --time=7:56:00
#SBATCH --mem=3GB
#SBATCH -p med
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=4
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/MergedReads
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/MergedReads

for r in $(seq 1 10)
do
HAPLOTYPE=1H
perl -pe "s#^@(\d+_)#\@$r\Batch$HAPLOTYPE.\1#" /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/Haplotype1/$r/simulated_1.fastq >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/MergedReads/merged-ShermanSimulation-Haplotype1_1.fastq &
perl -pe "s#^@(\d+_)#\@$r\Batch$HAPLOTYPE.\1#" /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/Haplotype1/$r/simulated_2.fastq >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/MergedReads/merged-ShermanSimulation-Haplotype1_2.fastq &
HAPLOTYPE=2H
perl -pe "s#^@(\d+_)#\@$r\Batch$HAPLOTYPE.\1#" /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/Haplotype2/$r/simulated_1.fastq >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/MergedReads/merged-ShermanSimulation-Haplotype2_1.fastq &
perl -pe "s#^@(\d+_)#\@$r\Batch$HAPLOTYPE.\1#" /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/Haplotype2/$r/simulated_2.fastq >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/MergedReads/merged-ShermanSimulation-Haplotype2_2.fastq &
wait
done

cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/MergedReads/merged-ShermanSimulation-Haplotype1_1.fastq /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/MergedReads/merged-ShermanSimulation-Haplotype2_1.fastq > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/MergedReads/merged-ShermanSimulation-Diploid_1.fastq &
cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/MergedReads/merged-ShermanSimulation-Haplotype1_2.fastq /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/MergedReads/merged-ShermanSimulation-Haplotype2_2.fastq > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/MergedReads/merged-ShermanSimulation-Diploid_2.fastq &
wait

rm /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/MergedReads/merged-ShermanSimulation-Haplotype1_1.fastq
rm /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/MergedReads/merged-ShermanSimulation-Haplotype1_2.fastq
rm /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/MergedReads/merged-ShermanSimulation-Haplotype2_1.fastq
rm /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/MergedReads/merged-ShermanSimulation-Haplotype2_2.fastq

echo 'Sherman Simulation Read Merge Job Complete!!!'
