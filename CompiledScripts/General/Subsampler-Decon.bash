#!/bin/bash
#SBATCH --time=1-23:00:00
#SBATCH --mem=64GB
#SBATCH -p long
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=20
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Records
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/FastQs

rm -r /scratch/mjuswilc
mkdir /scratch/mjuswilc
mkdir /scratch/mjuswilc/Subsampled

#Subsample Blood and Brain reads for mixtures of equal reads at P percent blood--Blood Seqs=351593992; Brain Seqs=437735158
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/RawData/Sequences/Bisulphite/Blood/Fastq/SRR2070791_1.fastq.gz | grep -c "@" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Records/fBLOODSeqs.txt &
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/RawData/Sequences/Bisulphite/Brain/Fastq/SRR2070790_1.fastq.gz | grep -c "@" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Records/fBRAINSeqs.txt &
wait
BloodSeqs=`cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Records/fBLOODSeqs.txt`
BrainSeqs=`cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Records/fBRAINSeqs.txt`
rADJUST=`echo $BloodSeqs | awk -v brain=$BrainSeqs '{print $0/brain}'`
SEED=0
for P in 0 10 20 30 40 50 60 70 80 90 100
do
SEED=`expr $SEED + 23`
rBLOOD=`echo $P | awk -v reads=$BloodSeqs -v adjust=$rADJUST '$0==100{printf "%.9f", 0.999999999; next}{print ($0/100)}'`
rBRAIN=`echo $P | awk -v reads=$BrainSeqs -v adjust=$rADJUST '{print (1-($0/100))*adjust}'`
/home/mjuswilc/Programs/seqtk-master/seqtk sample -s$SEED  /work/mjuswilc/Domain/GreatTits/ultDeconvolution/RawData/Sequences/Bisulphite/Blood/Fastq/SRR2070791_1.fastq.gz $rBLOOD | gzip -9 > /scratch/mjuswilc/Subsampled/mAbel-Blood.Blood-$P.read1.fastq.gz &
/home/mjuswilc/Programs/seqtk-master/seqtk sample -s$SEED  /work/mjuswilc/Domain/GreatTits/ultDeconvolution/RawData/Sequences/Bisulphite/Blood/Fastq/SRR2070791_2.fastq.gz $rBLOOD | gzip -9 > /scratch/mjuswilc/Subsampled/mAbel-Blood.Blood-$P.read2.fastq.gz &
/home/mjuswilc/Programs/seqtk-master/seqtk sample -s$SEED  /work/mjuswilc/Domain/GreatTits/ultDeconvolution/RawData/Sequences/Bisulphite/Brain/Fastq/SRR2070790_1.fastq.gz $rBRAIN | gzip -9 > /scratch/mjuswilc/Subsampled/mAbel-Brain.Blood-$P.read1.fastq.gz &
/home/mjuswilc/Programs/seqtk-master/seqtk sample -s$SEED  /work/mjuswilc/Domain/GreatTits/ultDeconvolution/RawData/Sequences/Bisulphite/Brain/Fastq/SRR2070790_2.fastq.gz $rBRAIN | gzip -9 > /scratch/mjuswilc/Subsampled/mAbel-Brain.Blood-$P.read2.fastq.gz &
JOBS=`jobs | wc | awk '{print $1}'`
while [ $JOBS -ge 16 ]
do
sleep 3s
JOBS=`jobs | wc | awk '{print $1}'`
done
done
wait

#Mix blood and brain reads
for P in 0 10 20 30 40 50 60 70 80 90 100
do
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/FastQs/$P
zcat /scratch/mjuswilc/Subsampled/mAbel-Blood.Blood-$P.read1.fastq.gz /scratch/mjuswilc/Subsampled/mAbel-Brain.Blood-$P.read1.fastq.gz | gzip -9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/FastQs/$P/mAbel-BrainMix.Blood-$P.read1.fastq.gz &
zcat /scratch/mjuswilc/Subsampled/mAbel-Blood.Blood-$P.read2.fastq.gz /scratch/mjuswilc/Subsampled/mAbel-Brain.Blood-$P.read2.fastq.gz | gzip -9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/FastQs/$P/mAbel-BrainMix.Blood-$P.read2.fastq.gz &
JOBS=`jobs | wc | awk '{print $1}'`
while [ $JOBS -ge 12 ]
do
sleep 3s
JOBS=`jobs | wc | awk '{print $1}'`
done
done
wait
rm -r /scratch/mjuswilc
echo 'Job Complete!!!'
