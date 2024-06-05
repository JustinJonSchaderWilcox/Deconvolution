#!/bin/bash
#SBATCH --time=1-23:40:00
#SBATCH --mem=8GB
#SBATCH -p long
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

r=$1

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/Haplotype1/$r
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/Haplotype1/$r

#Haplotype 1 Simulation
cd /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/Haplotype1/$r
/home/mjuswilc/Programs/Sherman/Sherman --genome_folder /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/SimulatedGenome/Haplotype1 -l 100 -n 102031077 -e 0.1 -pe -I 200 -X 600 --CG 54 -CH 99
echo 'Read Simulation Job Complete!!!'
