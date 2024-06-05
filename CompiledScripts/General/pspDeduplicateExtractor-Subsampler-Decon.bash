#!/bin/bash
#SBATCH --time=1:23
#SBATCH --mem=1GB
#SBATCH -p short
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DeduplicationExtractor
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DeduplicationExtractor

echo -e "P\tALIGNMENTS\tDUPLICATES\tPOSITIONS\tSEQUENCES" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DeduplicationExtractor/AllP-DeduplicationSEE.txt
for P in 0 10 20 30 40 50 60 70 80 90 100
do
ALIGNMENTS=`grep "^Total number of alignments analysed" /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Dedup/$P/Bismark-Bowtie.mAbel-BrainMix.Blood-$P.read1.trimmomatic.paired_bismark_bt2_pe.deduplication_report.txt | perl -pe 's#.*\s+(\d+)$#$1#'`
DUPLICATES=`grep "^Total number duplicated alignments removed" /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Dedup/$P/Bismark-Bowtie.mAbel-BrainMix.Blood-$P.read1.trimmomatic.paired_bismark_bt2_pe.deduplication_report.txt | perl -pe 's#.*\s+(\d+)\s+\(.*\)$#$1#'`
POSITIONS=`grep "^Duplicated alignments were found at" /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Dedup/$P/Bismark-Bowtie.mAbel-BrainMix.Blood-$P.read1.trimmomatic.paired_bismark_bt2_pe.deduplication_report.txt | perl -pe 's#.*\s+(\d+)\s+.*$#$1#'`
SEQUENCES=`grep "^Total count of deduplicated leftover sequences" /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Bismark/Dedup/$P/Bismark-Bowtie.mAbel-BrainMix.Blood-$P.read1.trimmomatic.paired_bismark_bt2_pe.deduplication_report.txt | perl -pe 's#.*\s+(\d+)\s+.*$#$1#'`

echo -e "$P\t$ALIGNMENTS\t$DUPLICATES\t$POSITIONS\t$SEQUENCES" >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DeduplicationExtractor/AllP-DeduplicationSEE.txt
done

echo 'Deduplication Information Extractor Complete!!!'
