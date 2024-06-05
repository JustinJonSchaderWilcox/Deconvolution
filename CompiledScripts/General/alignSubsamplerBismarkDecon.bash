#!/bin/bash
#SBATCH --time=1-23:40:00
#SBATCH --mem=128GB
#SBATCH -p long
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=12
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

P=$1
t=$2

source activate Genomics
export PATH="$PATH:/home/mjuswilc/Programs/bowtie2-2.5.1"

sleep "$t"m

rm -r /scratch/Wilcox/Bismark/$P
mkdir -p /scratch/Wilcox/Bismark/$P

function clear_scratch {
rm -r /scratch/Wilcox
}

function bismark_align {
cp /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/FastQs/Trimmed/$P/mAbel-BrainMix.Blood-$P.read1.trimmomatic.paired.fastq.gz /scratch/Wilcox/Bismark/$P/mAbel-BrainMix.Blood-$P.read1.trimmomatic.paired.fastq.gz
cp /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/FastQs/Trimmed/$P/mAbel-BrainMix.Blood-$P.read2.trimmomatic.paired.fastq.gz /scratch/Wilcox/Bismark/$P/mAbel-BrainMix.Blood-$P.read2.trimmomatic.paired.fastq.gz

#Run Bismark Alignment
rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/$P
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/$P
/home/mjuswilc/Programs/Bismark-0.22.3/bismark --bowtie2 --genome_folder /work/mjuswilc/Domain/GreatTits/ultDeconvolution/BisulphiteGenome/Bismark --parallel 1 -p 4 -L 32 -o /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/$P --temp_dir /scratch/Wilcox/Bismark/$P --prefix Bismark-Bowtie -1 /scratch/Wilcox/Bismark/$P/mAbel-BrainMix.Blood-$P.read1.trimmomatic.paired.fastq.gz -2 /scratch/Wilcox/Bismark/$P/mAbel-BrainMix.Blood-$P.read2.trimmomatic.paired.fastq.gz
mv /scratch/Wilcox/Bismark/$P/bismark-align.$P.log /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/$P
rm -r /scratch/Wilcox
}

trap -- 'clear_scratch' USR1 SIGQUIT SIGINT SIGABRT SIGKILL SIGHUP SIGTERM EXIT TERM
bismark_align &> /scratch/Wilcox/Bismark/$P/bismark-align.$P.log &
wait
trap - USR1 SIGQUIT SIGINT SIGABRT SIGKILL SIGHUP SIGTERM EXIT TERM

echo 'Job Complete!!!'
