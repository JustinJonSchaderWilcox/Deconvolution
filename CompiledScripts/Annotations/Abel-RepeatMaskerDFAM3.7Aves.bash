#!/bin/bash
#SBATCH --time=23:00:00
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

source activate Genomics

rm -r /work/mjuswilc/Programs/RepeatMasker/Libraries/HMM-Dfam_3.7
rm -r /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker
mkdir /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker
mkdir /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/Genome
mkdir /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/RepeatMasker-DFAM3.7

gzip -d -c /work/mjuswilc/Domain/GreatTits/General/RawData/Abel-Reference/GCF_001522545.3_Parus_major1.1_genomic.fna.gz > /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/Genome/Abel-GCF_001522545.3_Parus_major1.1_genomic.fna
samtools faidx /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/Genome/Abel-GCF_001522545.3_Parus_major1.1_genomic.fna

conda deactivate

cd /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/RepeatMasker-DFAM3.7
/work/mjuswilc/Programs/RepeatMasker/RepeatMasker -e hmmer -gccalc -s -a -pa 10 -species aves  -dir /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/RepeatMasker-DFAM3.7 /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/Genome/Abel-GCF_001522545.3_Parus_major1.1_genomic.fna

GENOMESIZE=`grep -v ">" /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/Genome/Abel-GCF_001522545.3_Parus_major1.1_genomic.fna | perl -pe -chomp | wc -m`

perl /work/mjuswilc/Programs/RepeatMasker/util/calcDivergenceFromAlign.pl -s /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/RepeatMasker-DFAM3.7/Abel-GCF_001522545.3_Parus_major1.1_genomic.fna.divsum -a  /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/RepeatMasker-DFAM3.7/Abel-GCF_001522545.3_Parus_major1.1_genomic.fna.GC-Adjusted.align /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/RepeatMasker-DFAM3.7/Abel-GCF_001522545.3_Parus_major1.1_genomic.fna.align
perl /work/mjuswilc/Programs/RepeatMasker/util/createRepeatLandscape.pl -div /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/RepeatMasker-DFAM3.7/Abel-GCF_001522545.3_Parus_major1.1_genomic.fna.divsum -t "Great Tit (Abel) Repeat Landscape" -g $GENOMESIZE > /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/RepeatMasker-DFAM3.7/Abel-GCF_001522545.3_Parus_major1.1.repeat_landscape.html
echo 'Job Complete'
