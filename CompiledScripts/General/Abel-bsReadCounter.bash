#!/bin/bash
#SBATCH --time=5:48:00
#SBATCH --mem=4GB
#SBATCH -p med
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=12
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out


source activate Genomics

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Data
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Summary

zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/RawData/Sequences/Bisulphite/Blood/Fastq/SRR2070791_1.fastq.gz | grep -P -c "^@SRR\d+\.\d+\s+" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Data/BloodReadCount_1.txt &
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/RawData/Sequences/Bisulphite/Blood/Fastq/SRR2070791_2.fastq.gz | grep -P -c "^@SRR\d+\.\d+\s+" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Data/BloodReadCount_2.txt &
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/RawData/Sequences/Bisulphite/Brain/Fastq/SRR2070790_1.fastq.gz | grep -P -c "^@SRR\d+\.\d+\s+" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Data/BrainReadCount_1.txt &
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/RawData/Sequences/Bisulphite/Brain/Fastq/SRR2070790_2.fastq.gz | grep -P -c "^@SRR\d+\.\d+\s+" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Data/BrainReadCount_2.txt &
wait

BLOOD1=`cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Data/BloodReadCount_1.txt`
BLOOD2=`cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Data/BloodReadCount_2.txt`
BRAIN1=`cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Data/BrainReadCount_1.txt`
BRAIN2=`cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Data/BrainReadCount_2.txt`
BLOODSUM=`expr $BLOOD1 + $BLOOD2`
BRAINSUM=`expr $BRAIN1 + $BRAIN2`

echo -e "Name\tFile\tTissue\tPair\tCount" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Summary/DeconAbelSummaryReadCountFASTQ.txt
echo -e "Blood Forward Reads:\tSRR2070791_1.fastq.gz\tBlood\tForward\t$BLOOD1" >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Summary/DeconAbelSummaryReadCountFASTQ.txt
echo -e "Blood Reverse Reads:\tSRR2070791_2.fastq.gz\tBlood\tReverse\t$BLOOD2" >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Summary/DeconAbelSummaryReadCountFASTQ.txt
echo -e "Brain Forward Reads:\tSRR2070790_1.fastq.gz\tBrain\tForward\t$BRAIN1" >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Summary/DeconAbelSummaryReadCountFASTQ.txt
echo -e "Brain Reverse Reads:\tSRR2070790_2.fastq.gz\tBrain\tReverse\t$BRAIN2" >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Summary/DeconAbelSummaryReadCountFASTQ.txt

echo 'Raw FastQ Read Summary Job Complete!!!'

function samplefaction {
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/FastQs/$P/mAbel-BrainMix.Blood-$P.read1.fastq.gz | grep -c "^@SRR2070791" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Data/BloodSampReads1-$P.txt
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/FastQs/$P/mAbel-BrainMix.Blood-$P.read2.fastq.gz | grep -c "^@SRR2070791" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Data/BloodSampReads2-$P.txt
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/FastQs/$P/mAbel-BrainMix.Blood-$P.read1.fastq.gz | grep -c "^@SRR2070790" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Data/BrainSampReads1-$P.txt
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/FastQs/$P/mAbel-BrainMix.Blood-$P.read2.fastq.gz | grep -c "^@SRR2070790" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Data/BrainSampReads2-$P.txt
}

for P in 0 10 20 30 40 50 60 70 80 90 100
do
samplefaction &
done
wait

function read_through {
samtools view /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Sort/mAbel-BrainMix.Blood-$P.deduplicated.sorted.bam | grep -c "^SRR2070791" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Data/AlignedDedupBloodCount-$P.txt
samtools view /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Sort/mAbel-BrainMix.Blood-$P.deduplicated.sorted.bam | grep -c "^SRR2070790" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Data/AlignedDedupBrainCount-$P.txt
}
for P in 0 10 20 30 40 50 60 70 80 90 100
do
read_through &
done
wait

P=0
echo -e "P\tSubsampled\tBloodReads\tBrainReads\tPropBloodReads\tTrimmedPE\tAligned\tDeduplicated\tBlood\tBrain\tBloodProportion" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Summary/DeconAbelSummaryReadCountSampled.txt
for l in $(ls /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Scripts/Submissions/job.318073*.out)
do
SUBSAMPLED=`grep '^Input Read Pairs:' $l | perl -pe 's#Input Read Pairs:\s+(\d+).*#$1#'`
blREADS1=`cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Data/BloodSampReads1-$P.txt`
blREADS2=`cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Data/BloodSampReads2-$P.txt`
brREADS1=`cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Data/BrainSampReads1-$P.txt`
brREADS2=`cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Data/BrainSampReads2-$P.txt`
blALL=`expr $blREADS1 + $blREADS2`
brALL=`expr $brREADS1 + $brREADS2`
blReadsPROP=`echo $blALL | awk -v brain=$brALL '{print ($0/($0+brain))}'`
TRIMMED=`grep '^Input Read Pairs:' $l | perl -pe 's#.*Both Surviving:\s+(\d+).*#$1#'`
ALIGN=`grep "Total number of alignments analysed in" /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Dedup/$P/Bismark-Bowtie.mAbel-BrainMix.Blood-$P.read1.trimmomatic.paired_bismark_bt2_pe.deduplication_report.txt | perl -pe 's#.*\s+(\d+)$#$1#'`
DEDUP=`grep "Total count of deduplicated leftover sequences:" /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Dedup/$P/Bismark-Bowtie.mAbel-BrainMix.Blood-$P.read1.trimmomatic.paired_bismark_bt2_pe.deduplication_report.txt | perl -pe 's#Total count of deduplicated leftover sequences:\s+(\d+).*#$1#'`
BLOOD=`cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Data/AlignedDedupBloodCount-$P.txt`
BRAIN=`cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Data/AlignedDedupBrainCount-$P.txt`
blPROP=`echo $BLOOD | awk -v brain=$BRAIN '{print $0/($0+brain)}'`
echo -e "$P\t$SUBSAMPLED\t$blALL\t$brALL\t$blReadsPROP\t$TRIMMED\t$ALIGN\t$DEDUP\t$BLOOD\t$BRAIN\t$blPROP" >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/ReadCounts/Summary/DeconAbelSummaryReadCountSampled.txt
P=`expr $P + 10`
done

echo 'Subsampled Read Summary Job Complete!!!'
