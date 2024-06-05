#!/bin/bash
#SBATCH --time=14:00
#SBATCH --mem=58GB
#SBATCH -p short
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=20
# Output and error files
#SBATCH -o job.%J.out


source activate Genomics

#Directory Structure

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter
rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/BEDS
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/BEDS/Genomic
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/filtVCFs
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/ultVCFs
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/Mod-Filt99
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodSomatic
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BrainSomatic
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodMixSomatic
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BrainMixSomatic
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodMarkers
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BrainMarkers
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodMarkersMerged
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BrainMarkersMerged
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/Output


(awk 'BEGIN{OFS=FS="\t"}NR<3{next}{print $0}' /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/RepeatMasker-DFAM3.7/Abel-GCF_001522545.3_Parus_major1.1_genomic.fna.out; awk 'BEGIN{OFS=FS="\t"}NR<3{next}{print $0}' /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Repeats/RepeatMasker/RepeatMasker-RM2.0.4DFAM3.7/GCF_001522545.3_Parus_major1.1_genomic.fna.out) | awk 'BEGIN{OFS="\t"}NR<=3{next}{print $5 OFS $6 OFS $7}' | sort -k1,1V -k2,2n | uniq | awk 'BEGIN{OFS=FS="\t"}NR==1{scaff=$1; start=$2; end=$3; next}$1==scaff && ($2<end && $3<end){next}$1==scaff && ($2<end && $3>end){end=$3;next}$1!=scaff || $2>end{print scaff OFS start OFS end; scaff=$1; start=$2; end=$3}END{print scaff OFS start OFS end}' > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/BEDS/Genomic/GCF_001522545.3_Parus_major1.1.repeats.bed &
awk 'BEGIN{OFS="\t"}{print $1 OFS $2 OFS $3; print $4 OFS $5 OFS $6}' /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/SegmentalDuplications/BISER/Output/Abel-GCF_001522545.3_Parus_major1.1_genomic.fna.seg_dup.bedpe > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/BEDS/Genomic/GCF_001522545.3_Parus_major1.1.seg_dup.bed &
for P in 0 10 20 30 40 50 60 70 80 90 100
do
hDEPTH=`cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_"$P".depth.txt | awk '{print int(($0*2)+0.5)}'`
lDEPTH=`cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_"$P".depth.txt | awk '{print int(($0*0.5)+0.5)}'`
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/BEDS/$P
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz | grep -v "#" | perl -pe 's#(\S+)\s+(\d+).*ADF=(\d+),(\d+);ADR=(\d+),(\d+).*#$1\t$2\t$3\t$4\t$5\t$6#' | awk 'BEGIN{OFS=FS="\t"}2*($3+$4)<($5+$6){next}2*($5+$6)<($3+$4){next}{print $1 OFS $2 OFS $2}' > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/BEDS/$P/BloodBrainMix.$P.snp.vcf.read-symetrical.bed &
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz | awk 'BEGIN{OFS=FS="\t"}/^#/{next}$6<30 || $7=="Low"{print $1 OFS $2 OFS $2}' > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/BEDS/$P/BloodBrainMix.$P.snp.vcf.low-quality.bed &
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz | grep -v "#" | perl -pe 's#(\S+)\s+(\d+).*ADF=(\d+),(\d+);ADR=(\d+),(\d+).*#$1\t$2\t$3\t$4\t$5\t$6#' | awk -v low=$lDEPTH -v high=$hDEPTH 'BEGIN{OFS=FS="\t"}{depth=$3+$4+$5+$6}depth<low{print $1 OFS $2 OFS $2; next}depth>high{print $1 OFS $2 OFS $2; next}' > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/BEDS/$P/BloodBrainMix.$P.snp.vcf.low-depth.bed &
awk 'BEGIN{OFS=FS="\t"}NR==1{next}{print $1 OFS $2 OFS $2+1}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.CG > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/BEDS/$P/BloodBrainMix.$P.snp.vcf.CG.bed &
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz | awk 'BEGIN{OFS=FS="\t"}/^#/{next}length($5)>1{print $1 OFS $2 OFS $2}' > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/BEDS/$P/BloodBrainMix.$P.snp.vcf.multi-allelic.bed &
wait
done


function VCFerate {
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/filtVCFs/$P
(zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz | awk '!/^#/{exit}{print}'; bedtools subtract -a /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/BEDS/Genomic/GCF_001522545.3_Parus_major1.1.repeats.bed) | bgzip -l 9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/filtVCFs/$P/BloodBrainMix.$P.snp.no-repeats.vcf.gz
(zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz | awk '!/^#/{exit}{print}'; bedtools subtract -a /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/BEDS/Genomic/GCF_001522545.3_Parus_major1.1.seg_dup.bed) | bgzip -l 9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/filtVCFs/$P/BloodBrainMix.$P.snp.no-seg_dups.vcf.gz
(zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz | awk '!/^#/{exit}{print}'; bedtools intersect -a /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/BEDS/$P/BloodBrainMix.$P.snp.vcf.read-symetrical.bed) | bgzip -l 9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/filtVCFs/$P/BloodBrainMix.$P.snp.read-symetrical.vcf.gz
(zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz | awk '!/^#/{exit}{print}'; bedtools subtract -a /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/BEDS/$P/BloodBrainMix.$P.snp.vcf.low-quality.bed) | bgzip -l 9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/filtVCFs/$P/BloodBrainMix.$P.snp.Q30-pass.vcf.gz
(zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz | awk '!/^#/{exit}{print}'; bedtools subtract -a /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/BEDS/$P/BloodBrainMix.$P.snp.vcf.low-depth.bed) | bgzip -l 9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/filtVCFs/$P/BloodBrainMix.$P.snp.high-depth28X.vcf.gz
(zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz | awk '!/^#/{exit}{print}'; bedtools subtract -a /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/BEDS/$P/BloodBrainMix.$P.snp.vcf.CG.bed) | bgzip -l 9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/filtVCFs/$P/BloodBrainMix.$P.snp.nCpG.vcf.gz
(zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz | awk '!/^#/{exit}{print}'; bedtools subtract -a /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/BEDS/$P/BloodBrainMix.$P.snp.vcf.multi-allelic.bed) | bgzip -l 9 > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/filtVCFs/$P/BloodBrainMix.$P.snp.no-multi-allelic.vcf.gz
for v in $(ls /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/filtVCFs/$P/BloodBrainMix.$P.snp.*.vcf.gz)
do
bcftools index $v
done
bcftools isec -o /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/filtVCFs/$P/BloodBrainMix.$P.snp.fullFilt.vcf -n 7 -w 1 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/filtVCFs/$P/BloodBrainMix.$P.snp.no-repeats.vcf.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/filtVCFs/$P/BloodBrainMix.$P.snp.no-seg_dups.vcf.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/filtVCFs/$P/BloodBrainMix.$P.snp.read-symetrical.vcf.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/filtVCFs/$P/BloodBrainMix.$P.snp.Q30-pass.vcf.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/filtVCFs/$P/BloodBrainMix.$P.snp.high-depth28X.vcf.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/filtVCFs/$P/BloodBrainMix.$P.snp.nCpG.vcf.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/filtVCFs/$P/BloodBrainMix.$P.snp.no-multi-allelic.vcf.gz
(awk '!/^#/{exit}{print}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/filtVCFs/$P/BloodBrainMix.$P.snp.fullFilt.vcf; awk 'BEGIN{OFS=FS="\t"}/^#/{print $0; next}$4=="A" && $5=="T"{print $0; next}$4=="T" && $5=="A"{print $0; next}$4=="C" && $5=="A"{print $0; next}$4=="G" && $5=="T"{print $0; next}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/filtVCFs/$P/BloodBrainMix.$P.snp.fullFilt.vcf) > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/ultVCFs/BloodBrainMix.$P.snp.ultfFilt-sym.vcf
bgzip -l 9 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/ultVCFs/BloodBrainMix.$P.snp.ultfFilt-sym.vcf
bcftools index /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/ultVCFs/BloodBrainMix.$P.snp.ultfFilt-sym.vcf.gz
}

for P in 0 10 20 30 40 50 60 70 80 90 100
do
VCFerate &
done
wait
echo 'Filtration Files Created!!!'


function bgzipit {
bgzip -l 9 $v
bcftools index "$v".gz
}
for P in 0 10 20 30 40 50 60 70 80 90 100
do
bgzip -d  -c /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz | perl -pe 's#(.*DP=)(\d+)(.*,)(\d\.*\d*)$#$1$2$3$4\t$2\t$4#' | awk 'BEGIN{OFS=FS="\t"}/^#/{print $0; next}$12>0.25 && $12<0.75{print $0}' | perl -pe 's#\b\s+\d+\s+[\d\.]*\d$##' > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/Mod-Filt99/BloodBrainMix.$P.snp.vcf.mod-filt99.vcf &
done
wait
for v in $(ls /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/Mod-Filt99/BloodBrainMix.*.snp.vcf.mod-filt99.vcf)
do
bgzipit &
done
wait

for P in 0 10 20 30 40 50 60 70 80 90 100
do
bcftools isec -o /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodSomatic/BloodBrainMix.$P.snp.ultfFilt-sym.blood-somatic.vcf -n 1 -w 1 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/ultVCFs/BloodBrainMix.100.snp.ultfFilt-sym.vcf.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/Mod-Filt99/BloodBrainMix.$P.snp.vcf.mod-filt99.vcf.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/Mod-Filt99/BloodBrainMix.100.snp.vcf.mod-filt99.vcf.gz &
JOBS=`jobs -l | grep -c '^\['`
done
wait
for v in $(ls /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/*Somatic/BloodBrainMix.*.snp.ultfFilt-sym.*-somatic.vcf)
do
bgzipit &
done
wait

for P in 0 10 20 30 40 50 60 70 80 90 100
do
bcftools isec -o /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodMarkers/BloodBrainMix.$P.snp.ultfFilt-sym.blood-markers.vcf -n 2 -w 1 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodSomatic/BloodBrainMix.$P.snp.ultfFilt-sym.blood-somatic.vcf.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz &
bcftools isec -o /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodMarkers/BloodBrainMix.$P.snp.ultfFilt-sym.blood-mix-markers.vcf -n 2 -w 1 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodSomatic/BloodBrainMix.$P.snp.ultfFilt-sym.blood-somatic.vcf.gz &
JOBS=`jobs -l | grep -c '^\['`
while [ $JOBS -ge 18 ]
do
sleep 3s
JOBS=`jobs -l | grep -c '^\['`
done
done
wait
for v in $(ls /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/*Markers/BloodBrainMix.*.snp.ultfFilt-sym.*-markers.vcf)
do
bgzipit &
JOBS=`jobs -l | grep -c '^\['`
while [ $JOBS -ge 20 ]
do
sleep 3s
JOBS=`jobs -l | grep -c '^\['`
done
done
wait

for P in 0 10 20 30 40 50 60 70 80 90 100
do
bcftools merge --force-samples /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodSomatic/BloodBrainMix.$P.snp.ultfFilt-sym.blood-somatic.vcf.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodMarkers/BloodBrainMix.$P.snp.ultfFilt-sym.blood-mix-markers.vcf.gz | perl -pe 's#(.AC=*)(\d+)(.*)#$1$2$3\t$2#' | awk 'BEGIN{OFS=FS="\t"}/^#/{print $0; next}$12>2{next}{print $0}' | perl -pe 's#\s+\d+$##' > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodMarkersMerged/BloodBrainMix.100-$P.snp.ultfFilt-sym.blood-mix.vcf &
done
wait

echo -e "P\tPureAverage\tMixAverage\tAverageDiffAlleleFreq\tEstimatedP\tMarkerSNVs" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/Output/blood-deconvolution-validation-table.txt
for P in 0 10 20 30 40 50 60 70 80 90 100
do
hDEPTH=`cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_"$P".depth.txt | awk '{print int(($0*2)+0.5)}'`
lDEPTH=`cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_"$P".depth.txt | awk '{print int(($0*0.5)+0.5)}'`
grep -v "#" /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodMarkersMerged/BloodBrainMix.100-$P.snp.ultfFilt-sym.blood-mix.vcf | awk 'BEGIN{OFS=FS="\t"}NR==1{scaffold=$1; site=$2; prev=$0; next}$1==scaffold && ($2-site)<100000{scaffold=$1; site=$2; prev=$0; next}{scaffold=$1; site=$2; prev=$0; print prev}END{print $prev}' | perl -pe 's#:.$#,0#' | perl -pe 's#(\d\/\d:\d+:)(\d+),(\d+):(\d+),(\d+)([\d,:\.]+)$#$1$2,$3:$4,$5$6\t$2\t$3\t$4\t$5#' | awk -v low=$lDEPTH -v high=$hDEPTH 'BEGIN{OFS=FS="\t"}2*($12+$13)<($14+$15){next}2*($14+$15)<($12+$13){next}{depth=($12+$13+$14+$15)}(depth>0 && depth<low) || depth>high{next}{print $0}' | perl -pe 's#\s+[\d\.]+\s+[\d\.]+\s+[\d\.]+\s+[\d\.]+$##' | perl -pe 's#.*,([0-9\.]+)\s+.*,([0-9\.]+)$#$1\t$2#' | awk 'BEGIN{OFS=FS="\t"; diff=0; count=0; afreq1=0; afreq2=0}$1>0.5{$1=1-$1}$2>0.5{$2=1-$2}$1==0{next}{print $0}' | awk -v p=$P 'BEGIN{OFS=FS="\t"; diff=0; count=0; afreq1=0; afreq2=0}{afreq1=afreq1+$1; afreq2=afreq2+$2; diff=diff+$1-$2; count=count+1}count==0{count=1}END{print p OFS afreq1/count OFS afreq2/count OFS diff/count OFS 1-(diff/count)/(afreq1/count) OFS count}' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/Output/blood-deconvolution-validation-table.txt
done

echo 'Estimate of P Complete!!!!!'
