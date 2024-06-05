#!/bin/bash
#SBATCH --time=1-23:56:00
#SBATCH --mem=1GB
#SBATCH -p long
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=16
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out


rm /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/Haplotype1/*/*fastq.gz
rm /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/Haplotype2/*/*fastq.gz
for r in $(seq 1 10)
do
gzip -9 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/Haplotype1/$r/simulated_1.fastq &
gzip -9 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/Haplotype1/$r/simulated_2.fastq &
gzip -9 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/Haplotype2/$r/simulated_1.fastq &
gzip -9 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/Haplotype2/$r/simulated_2.fastq &
JOBS=`jobs | wc | awk '{print $1}'`
while [ $JOBS -ge 12 ]
do
sleep 3s
JOBS=`jobs | wc | awk '{print $1}'`
done
done
wait

echo 'Sherman Simulation Raw Read GZIP Complete!!!'
