#!/bin/bash
#SBATCH --time=1-23:40:00
#SBATCH --mem=168GB
#SBATCH -p long
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=24
# Output and error files
#SBATCH -o job.%J.out

source activate Genomics
export PATH="$PATH:/home/mjuswilc/Programs/bowtie2-2.5.1"

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/Bismark
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/Bismark
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/Bismark/Alignment

function clear_scratch {
rm -r /scratch/mjuswilc
}


function bismark_align {
rm -r /scratch/mjuswilc/BismarkAlign-HeterozygoteSim
mkdir -p /scratch/mjuswilc/BismarkAlign-HeterozygoteSim

cp /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/TrimmedReads/sbS100X-ShermanSimulation-Diploid_1.paired.fastq.gz /scratch/mjuswilc/BismarkAlign-HeterozygoteSim/sbS100X-ShermanSimulation-Diploid_1.paired.fastq.gz &
cp /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/TrimmedReads/sbS100X-ShermanSimulation-Diploid_2.paired.fastq.gz /scratch/mjuswilc/BismarkAlign-HeterozygoteSim/sbS100X-ShermanSimulation-Diploid_2.paired.fastq.gz &
wait

#Run Bismark Alignment
/home/mjuswilc/Programs/Bismark-0.22.3/bismark --bowtie2 --genome_folder /work/mjuswilc/Domain/GreatTits/ultDeconvolution/BisulphiteGenome/Bismark --parallel 2 -p 4 -L 32 -o /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/Bismark/Alignment --temp_dir /scratch/mjuswilc/BismarkAlign-HeterozygoteSim --prefix Bismark-Bowtie -1 /scratch/mjuswilc/BismarkAlign-HeterozygoteSim/sbS100X-ShermanSimulation-Diploid_1.paired.fastq.gz -2 /scratch/mjuswilc/BismarkAlign-HeterozygoteSim/sbS100X-ShermanSimulation-Diploid_2.paired.fastq.gz
rm -r /scratch/mjuswilc
}

trap -- 'clear_scratch' USR1 SIGQUIT SIGINT SIGABRT SIGKILL SIGHUP SIGTERM EXIT TERM
bismark_align &> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/HeterozygoteSimulation/GreatTit/Sherman/Bismark/Alignment/bismark-align.ShermanSimulation.log &
wait
trap - USR1 SIGQUIT SIGINT SIGABRT SIGKILL SIGHUP SIGTERM EXIT TERM

echo 'Heterozygote Simulation Bismark Aligment Job Complete!!!'
