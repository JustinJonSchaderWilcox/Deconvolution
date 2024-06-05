#!/bin/bash
#SBATCH --time=1:58:00
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

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/SimulatedGenome
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/SimulatedGenome
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/SimulatedGenome/Input
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/SimulatedGenome/Haplotype1
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/SimulatedGenome/Haplotype2

zcat /work/mjuswilc/Domain/GreatTits/General/RawData/Abel-Reference/GCF_001522545.3_Parus_major1.1_genomic.fna.gz > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/SimulatedGenome/Input/GCF_001522545.3_Parus_major1.1_genomic.fa

cp /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/SimulatedGenome/Input/GCF_001522545.3_Parus_major1.1_genomic.fa /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/SimulatedGenome/Haplotype1/great-tit-sim.Abel-GCF_001522545v3-Hap1.fasta
mutation-simulator -o /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/SimulatedGenome/Haplotype2/great-tit-sim.Abel-GCF_001522545v3-Hap2.fasta /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/SimulatedGenome/Input/GCF_001522545.3_Parus_major1.1_genomic.fa args -sn 0.001 -titv 3 -in 0.00005 -inmin 1 -inmax 8 -de 0.00005 -demin 1 -demax 8

echo 'Heterozygote Simulation Job Complete!!!'
