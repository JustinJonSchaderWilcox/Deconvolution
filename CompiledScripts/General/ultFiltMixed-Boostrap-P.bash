#!/bin/bash
#SBATCH --time=13:46
#SBATCH --mem=1GB
#SBATCH -p short
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=8
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ValidateP
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ValidateP
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ValidateP/TEMP
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ValidateP/BootStrap-P

echo -e "P\tBootStrap\tPureAverage\tMixAverage\tAverageDiffAlleleFreq\tEstimatedP\tMarkerSNVs" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ValidateP/BootStrap-P/blood-deconvolution-bootstrap100-table.txt
for P in 0 10 20 30 40 50 60 70 80 90 100
do
hDEPTH=`cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_"$P".depth.txt | awk '{print int(($0*2)+0.5)}'`
lDEPTH=`cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_"$P".depth.txt | awk '{print int(($0*0.5)+0.5)}'`

for c in $(seq 1 100)
do
grep -v "#" /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodMarkersMerged/BloodBrainMix.100-$P.snp.ultfFilt-sym.blood-mix.vcf | awk 'BEGIN{OFS=FS="\t"}NR==1{scaffold=$1; site=$2; prev=$0; next}$1==scaffold && ($2-site)<100000{scaffold=$1; site=$2; prev=$0; next}{scaffold=$1; site=$2; prev=$0; print prev}END{print $prev}' > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ValidateP/TEMP/BloodBrainMix.100-$P.snp.ultfFilt-sym.blood-mix.pruned.pseudo.vcf

cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ValidateP/TEMP/BloodBrainMix.100-$P.snp.ultfFilt-sym.blood-mix.pruned.pseudo.vcf | shuf -n 100 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ValidateP/TEMP/blood-bootsrap"$c".pseudo_vcf

cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ValidateP/TEMP/blood-bootsrap"$c".pseudo_vcf | awk 'BEGIN{OFS=FS="\t"}NR==1{scaffold=$1; site=$2; prev=$0; next}$1==scaffold && ($2-site)<100000{scaffold=$1; site=$2; prev=$0; next}{scaffold=$1; site=$2; prev=$0; print prev}END{print $prev}' | perl -pe 's#:.$#,0#' | perl -pe 's#(\d\/\d:\d+:)(\d+),(\d+):(\d+),(\d+)([\d,:\.]+)$#$1$2,$3:$4,$5$6\t$2\t$3\t$4\t$5#' | awk -v low=$lDEPTH -v high=$hDEPTH 'BEGIN{OFS=FS="\t"}2*($12+$13)<($14+$15){next}2*($14+$15)<($12+$13){next}{depth=($12+$13+$14+$15)}(depth>0 && depth<low) || depth>high{next}{print $0}' | perl -pe 's#\s+[\d\.]+\s+[\d\.]+\s+[\d\.]+\s+[\d\.]+$##' | perl -pe 's#.*,([0-9\.]+)\s+.*,([0-9\.]+)$#$1\t$2#' | awk 'BEGIN{OFS=FS="\t"; diff=0; count=0; afreq1=0; afreq2=0}$1>0.5{$1=1-$1}$2>0.5{$2=1-$2}$1==0{next}{print $0}' | awk -v p=$P -v bootstrap=$c 'BEGIN{OFS=FS="\t"; diff=0; count=0; afreq1=0; afreq2=0}{afreq1=afreq1+$1; afreq2=afreq2+$2; diff=diff+$1-$2; count=count+1}count==0{count=1}END{print p OFS bootstrap OFS afreq1/count OFS afreq2/count OFS diff/count OFS 1-(diff/count)/(afreq1/count) OFS count}' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ValidateP/BootStrap-P/blood-deconvolution-bootstrap100-table.txt
done
done

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ValidateP/TEMP

echo 'Bootstrap of P Estimate Complete!!!!!'
