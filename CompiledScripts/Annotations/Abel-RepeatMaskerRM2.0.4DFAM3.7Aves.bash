#!/bin/bash
#SBATCH --time=7:40:00
#SBATCH --mem=64GB
#SBATCH -p med
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=20
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

rm -r /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7
mkdir /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7

cd /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7
cp /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatModeler2/Abel-GCF_001522545.3_Parus_major1.1_genomic-families.fa /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/Abel-GCF_001522545.3_Parus_major1.1_genomic-families.fa
zcat /work/mjuswilc/Domain/GreatTits/General/RawData/Abel-Reference/GCF_001522545.3_Parus_major1.1_genomic.fna.gz > /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/GCF_001522545.3_Parus_major1.1_genomic.fna
/work/mjuswilc/Programs/RepeatMasker/RepeatMasker -e rmblast -gccalc -s -a -pa 20 -lib  /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/Abel-GCF_001522545.3_Parus_major1.1_genomic-families.fa /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatModeler2/Abel-GCF_001522545.3_Parus_major1.1_genomic.fna /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/GCF_001522545.3_Parus_major1.1_genomic.fna

GENOMESIZE=`grep -v ">" /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/Genome/Abel-GCF_001522545.3_Parus_major1.1_genomic.fna | perl -pe -chomp | wc -m`

perl /work/mjuswilc/Programs/RepeatMasker/util/calcDivergenceFromAlign.pl -s /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/Abel-GCF_001522545.3_Parus_major1.1_genomic.fna.divsum -a  /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/Abel-GCF_001522545.3_Parus_major1.1_genomic.fna.GC-Adjusted.align /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/GCF_001522545.3_Parus_major1.1_genomic.fna.align
perl /work/mjuswilc/Programs/RepeatMasker/util/createRepeatLandscape.pl -div /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/Abel-GCF_001522545.3_Parus_major1.1_genomic.fna.divsum -t "Great Tit (Abel) Repeat Landscape" -g $GENOMESIZE > /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/Abel-GCF_001522545.3_Parus_major1.1.repeat_landscape.html
echo 'Job Complete'
