#!/bin/bash
#SBATCH --time=8:40
#SBATCH --mem=1GB
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

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/SomaticSiteIntercept

cp /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/Output/blood-deconvolution-validation-table.txt /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/blood-deconvolution-validation-table.txt

function bgzipit {
bgzip -l 9 $v
bcftools index "$v".gz
}
for P in 0 10 20 30 40 50 60 70 80 90 100
do
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/SomaticSiteIntercept/$P
bcftools isec -o /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/SomaticSiteIntercept/$P/BloodBrainMix.$P.SNV-MarkerIntercpt.vcf -n2 -w1 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodMarkers/BloodBrainMix.$P.snp.ultfFilt-sym.blood-mix-markers.vcf.gz &
done
wait
for v in $(ls /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/SomaticSiteIntercept/*/BloodBrainMix.*.SNV-MarkerIntercpt.vcf)
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

function SNV {
echo -e "P\tAFREQ" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-allSNV-AFREQ.txt
echo -e "P\tDepth" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-allSNV-Depth.txt
echo -e "P\tPhred" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-allSNV-Phred.txt
for P in 0 10 20 30 40 50 60 70 80 90 100
do
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz | grep -v "#" | perl -pe 's#.*,([\d\.]+)$#$1#' | awk -v p=$P 'BEGIN{OFS=FS="\t"}{print p OFS $0}' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-allSNV-AFREQ.txt &
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz | grep -v "#" | perl -pe 's#.*DP=([\d\.]+).*#$1#' | awk -v p=$P 'BEGIN{OFS=FS="\t"}{print p OFS $0}' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-allSNV-Depth.txt &
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz | grep -v "#" | awk -v p=$P 'BEGIN{OFS=FS="\t"}/^#/{next}{print p OFS $6}' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-allSNV-Phred.txt &
wait
done
}
function ultFiltSNV {
echo -e "P\tAFREQ" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-ultFiltSNV-AFREQ.txt
echo -e "P\tDepth" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-ultFiltSNV-Depth.txt
echo -e "P\tPhred" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-ultFiltSNV-Phred.txt
for P in 0 10 20 30 40 50 60 70 80 90 100
do
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/ultVCFs/BloodBrainMix.$P.snp.ultfFilt-sym.vcf.gz | grep -v "#" | perl -pe 's#.*,([\d\.]+)$#$1#' | awk -v p=$P 'BEGIN{OFS=FS="\t"}{print p OFS $0}' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-ultFiltSNV-AFREQ.txt &
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/ultVCFs/BloodBrainMix.$P.snp.ultfFilt-sym.vcf.gz | grep -v "#" | perl -pe 's#.*DP=([\d\.]+).*#$1#' | awk -v p=$P 'BEGIN{OFS=FS="\t"}{print p OFS $0}' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-ultFiltSNV-Depth.txt &
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/ultVCFs/BloodBrainMix.$P.snp.ultfFilt-sym.vcf.gz | grep -v "#" | awk -v p=$P 'BEGIN{OFS=FS="\t"}/^#/{next}{print p OFS $6}' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-ultFiltSNV-Phred.txt &
wait
done
}
function markersSomaticSNV {
echo -e "P\tAFREQ" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-MarkersSomaticSNV-AFREQ.txt
echo -e "P\tDepth" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-MarkersSomaticSNV-Depth.txt
echo -e "P\tPhred" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-MarkersSomaticSNV-Phred.txt
for P in 0 10 20 30 40 50 60 70 80 90 100
do
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodSomatic/BloodBrainMix.$P.snp.ultfFilt-sym.blood-somatic.vcf.gz | grep -v "#" | perl -pe 's#.*,([\d\.]+)$#$1#' | awk -v p=$P 'BEGIN{OFS=FS="\t"}{print p OFS $0}' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-MarkersSomaticSNV-AFREQ.txt &
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodSomatic/BloodBrainMix.$P.snp.ultfFilt-sym.blood-somatic.vcf.gz | grep -v "#" | perl -pe 's#.*DP=([\d\.]+).*#$1#' | awk -v p=$P 'BEGIN{OFS=FS="\t"}{print p OFS $0}' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-MarkersSomaticSNV-Depth.txt &
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodSomatic/BloodBrainMix.$P.snp.ultfFilt-sym.blood-somatic.vcf.gz | grep -v "#" | awk -v p=$P 'BEGIN{OFS=FS="\t"}/^#/{next}{print p OFS $6}' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-MarkersSomaticSNV-Phred.txt &
wait
done
}
function SomaticSNV {
echo -e "P\tAFREQ" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-bloodSomaticSNV-AFREQ.txt
echo -e "P\tDepth" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-bloodSomaticSNV-Depth.txt
echo -e "P\tPhred" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-bloodSomaticSNV-Phred.txt
for P in 0 10 20 30 40 50 60 70 80 90 100
do
hDEPTH=`cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_"$P".depth.txt | awk '{print int(($0*2)+0.5)}'`
lDEPTH=`cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_"$P".depth.txt | awk '{print int(($0*0.5)+0.5)}'`

grep -v "#" /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodMarkersMerged/BloodBrainMix.100-$P.snp.ultfFilt-sym.blood-mix.vcf | awk 'BEGIN{OFS=FS="\t"}NR==1{scaffold=$1; site=$2; prev=$0; next}$1==scaffold && ($2-site)<100000{scaffold=$1; site=$2; prev=$0; next}{scaffold=$1; site=$2; prev=$0; print prev}END{print $prev}' | perl -pe 's#:.$#,0#' | perl -pe 's#(\d\/\d:\d+:)(\d+),(\d+):(\d+),(\d+)([\d,:\.]+)$#$1$2,$3:$4,$5$6\t$2\t$3\t$4\t$5#' | awk -v low=$lDEPTH -v high=$hDEPTH 'BEGIN{OFS=FS="\t"}2*($12+$13)<($14+$15){next}2*($14+$15)<($12+$13){next}{depth=($12+$13+$14+$15)}(depth>0 && depth<low) || depth>high{next}{print $0}' | perl -pe 's#\s+[\d\.]+\s+[\d\.]+\s+[\d\.]+\s+[\d\.]+$##' | perl -pe 's#.*,([0-9\.]+)\s+.*,([0-9\.]+)$#$1\t$2#' | awk -v p=$P 'BEGIN{OFS=FS="\t"; diff=0; count=0; afreq1=0; afreq2=0}$1>0.5{fold1=1-$1}$2>0.5{fold2=1-$2}fold1==0{next}{print p OFS $2}' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-bloodSomaticSNV-AFREQ.txt &
grep -v "#" /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodMarkersMerged/BloodBrainMix.100-$P.snp.ultfFilt-sym.blood-mix.vcf | awk 'BEGIN{OFS=FS="\t"}NR==1{scaffold=$1; site=$2; prev=$0; next}$1==scaffold && ($2-site)<100000{scaffold=$1; site=$2; prev=$0; next}{scaffold=$1; site=$2; prev=$0; print prev}END{print $prev}' | grep -v './.:.:.:.:.:.:.' | perl -pe 's#:.$#,0#' | perl -pe 's#(\d\/\d:\d+:)(\d+),(\d+):(\d+),(\d+)([\d,:\.]+)$#$1$2,$3:$4,$5$6\t$2\t$3\t$4\t$5#' | awk -v low=$lDEPTH -v high=$hDEPTH 'BEGIN{OFS=FS="\t"}2*($12+$13)<($14+$15){next}2*($14+$15)<($12+$13){next}{depth=($12+$13+$14+$15)}(depth>0 && depth<low) || depth>high{next}{print $0}' | perl -pe 's#\s+[\d\.]+\s+[\d\.]+\s+[\d\.]+\s+[\d\.]+$##' | perl -pe 's#.*\s+\d/\d:\S+\s+\d/\d:(\d+):.*([0-9\.]+),([0-9\.]+)$#$1\t$2\t$3#' | awk -v p=$P 'BEGIN{OFS=FS="\t"}$2>0.5{fold1=1-$2}$3>0.5{fold2=1-$3}fold1==0{next}{print p OFS $1}' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-bloodSomaticSNV-Depth.txt &
grep -v "#" /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodMarkersMerged/BloodBrainMix.100-$P.snp.ultfFilt-sym.blood-mix.vcf | awk 'BEGIN{OFS=FS="\t"}NR==1{scaffold=$1; site=$2; prev=$0; next}$1==scaffold && ($2-site)<100000{scaffold=$1; site=$2; prev=$0; next}{scaffold=$1; site=$2; prev=$0; print prev}END{print $prev}' | perl -pe 's#:.$#,0#' | perl -pe 's#(\d\/\d:\d+:)(\d+),(\d+):(\d+),(\d+)([\d,:\.]+)$#$1$2,$3:$4,$5$6\t$2\t$3\t$4\t$5#' | awk -v low=$lDEPTH -v high=$hDEPTH 'BEGIN{OFS=FS="\t"}2*($12+$13)<($14+$15){next}2*($14+$15)<($12+$13){next}{depth=($12+$13+$14+$15)}(depth>0 && depth<low) || depth>high{next}{print $0}' | perl -pe 's#\s+[\d\.]+\s+[\d\.]+\s+[\d\.]+\s+[\d\.]+$##' | perl -pe 's#(\S+)\s+(\d+).*,([0-9\.]+)\s+.*,([0-9\.]+)$#$1\t$2\t$3\t$4#' | awk -v p=$P 'BEGIN{OFS=FS="\t"}$3>0.5{fold1=1-$3}$4>0.5{fold2=1-$4}fold1==0{next}(1-(fold1-fold2)/fold1){print $1 OFS $2 OFS $2}' | bedtools intersect -wb -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/SomaticSiteIntercept/$P/BloodBrainMix.$P.SNV-MarkerIntercpt.vcf.gz | awk -v p=$P 'BEGIN{OFS=FS="\t"; phred=0}{print p OFS $9}' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-bloodSomaticSNV-Phred.txt &
wait
done
}

SNV &
ultFiltSNV &
markersSomaticSNV &
SomaticSNV &
wait


#BS-SNper called  SNVs; Filter-Site SNVs; Putative-Somatic SNVS; Depth SNV; Freq SNVs; Quality Score SNVs
echo -e "P\tSNVs\tultFiltSNVs\tSomaticSNVsMarkers\tSomaticSNVs\tSNVsDEPTH\tultFiltSNVsDEPTH\tSomaticSNVsMarkersDEPTH\tSomaticSNVsDEPTH\tSNVsAFREQ\tultFiltSNVsAFREQ\tSomaticSNVsMarkersAFREQ\tSomaticSNVsAFREQ\tSNVsPFRED\tultFiltSNVsPFRED\tSomaticSNVsMarkersPFRED\tSomaticSNVsPFRED" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/deconSNV-and-P-Summary.txt
for P in 0 10 20 30 40 50 60 70 80 90 100
do
hDEPTH=`cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_"$P".depth.txt | awk '{print int(($0*2)+0.5)}'`
lDEPTH=`cat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/DepthFinder/gtSub_"$P".depth.txt | awk '{print int(($0*0.5)+0.5)}'`

SNVs=`zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz | grep -v -c "#"`
ultFiltSNVs=`zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/ultVCFs/BloodBrainMix.$P.snp.ultfFilt-sym.vcf.gz | grep -v -c "#"`
SomaticSNVsMarkers=`zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodSomatic/BloodBrainMix.$P.snp.ultfFilt-sym.blood-somatic.vcf.gz | grep -v -c "#"`
SomaticSNVs=`awk -v p=$P 'BEGIN{OFS=FS="\t"; count=0}NR==1{next}$1==p{count++}END{print count}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-bloodSomaticSNV-Depth.txt`
SNVsDEPTH=`zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz |  perl -pe 's#.*DP=([\d\.]+).*#$1#' | awk 'BEGIN{count=0}/^#/{next}{count=count+$0}END{print count/NR}'`
ultFiltSNVsDEPTH=`zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/ultVCFs/BloodBrainMix.$P.snp.ultfFilt-sym.vcf.gz | perl -pe 's#.*DP=([\d\.]+).*#$1#' | awk 'BEGIN{count=0}/^#/{next}{count=count+$0}END{print count/NR}'`
SomaticSNVsMarkersDEPTH=`zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodSomatic/BloodBrainMix.$P.snp.ultfFilt-sym.blood-somatic.vcf.gz | perl -pe 's#.*DP=([\d\.]+).*#$1#' | awk 'BEGIN{count=0}/^#/{next}{count=count+$0}END{print count/NR}'`
SomaticSNVsDEPTH=`awk -v p=$P 'BEGIN{OFS=FS="\t"; depth=0; count=0}NR==1{next}$1==p{depth=depth+$2; count++}count==0{count=1}END{print depth/count}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-bloodSomaticSNV-Depth.txt`
SNVsAFREQ=`zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz | perl -pe 's#.*,([\d\.]+)$#$1#' | awk 'BEGIN{count=0}/^#/{next}$0>0.5{$0=1-$0}{count=count+$0}END{print count/NR}'`
ultFiltSNVsAFREQ=`zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/ultVCFs/BloodBrainMix.$P.snp.ultfFilt-sym.vcf.gz | perl -pe 's#.*,([\d\.]+)$#$1#' | awk 'BEGIN{count=0}/^#/{next}$0>0.5{$0=1-$0}{count=count+$0}END{print count/NR}'`
SomaticSNVsMarkersAFREQ=`zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodSomatic/BloodBrainMix.$P.snp.ultfFilt-sym.blood-somatic.vcf.gz | perl -pe 's#.*,([\d\.]+)$#$1#' | awk 'BEGIN{count=0}/^#/{next}$0>0.5{$0=1-$0}{count=count+$0}END{print count/NR}'`
SomaticSNVsAFREQ=`awk -v p=$P 'BEGIN{OFS=FS="\t"; afreq=0; count=0}NR==1{next}$1==p{afreq=afreq+$2; count++}count==0{count=1}END{print afreq/count}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/decon-bloodSomaticSNV-AFREQ.txt`


SNVsPFRED=`zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/SNPs/$P/BloodBrainMix.$P.snp.vcf.gz | awk 'BEGIN{count=0}/^#/{next}{count=count+$6}END{print count/NR}'`
ultFiltSNVsPFRED=`zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/ultFilter/ultVCFs/BloodBrainMix.$P.snp.ultfFilt-sym.vcf.gz | awk 'BEGIN{count=0}/^#/{next}{count=count+$6}END{print count/NR}'`
SomaticSNVsMarkersPFRED=`zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodSomatic/BloodBrainMix.$P.snp.ultfFilt-sym.blood-somatic.vcf.gz | awk 'BEGIN{count=0}/^#/{next}{count=count+$6}END{print count/NR}'`
SomaticSNVsPFRED=`grep -v "#" /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodMarkersMerged/BloodBrainMix.100-$P.snp.ultfFilt-sym.blood-mix.vcf | awk 'BEGIN{OFS=FS="\t"}NR==1{scaffold=$1; site=$2; prev=$0; next}$1==scaffold && ($2-site)<100000{scaffold=$1; site=$2; prev=$0; next}{scaffold=$1; site=$2; prev=$0; print prev}END{print $prev}' | perl -pe 's#:.$#,0#' | perl -pe 's#(\d\/\d:\d+:)(\d+),(\d+):(\d+),(\d+)([\d,:\.]+)$#$1$2,$3:$4,$5$6\t$2\t$3\t$4\t$5#' | awk -v low=$lDEPTH -v high=$hDEPTH 'BEGIN{OFS=FS="\t"}2*($12+$13)<($14+$15){next}2*($14+$15)<($12+$13){next}{depth=($12+$13+$14+$15)}(depth>0 && depth<low) || depth>high{next}{print $0}' | perl -pe 's#\s+[\d\.]+\s+[\d\.]+\s+[\d\.]+\s+[\d\.]+$##' | perl -pe 's#(\S+)\s+(\d+).*,([0-9\.]+)\s+.*,([0-9\.]+)$#$1\t$2\t$3\t$4#' | awk -v p=$P 'BEGIN{OFS=FS="\t"}$3>0.5{fold1=1-$3}$4>0.5{fold2=1-$4}fold1==0{next}{print $1 OFS $2 OFS $2}' | bedtools intersect -wb -a - -b /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/SomaticSiteIntercept/$P/BloodBrainMix.$P.SNV-MarkerIntercpt.vcf.gz | awk 'BEGIN{phred=0}{phred=phred+$9}END{print phred/NR}'`
echo -e "$P\t$SNVs\t$ultFiltSNVs\t$SomaticSNVsMarkers\t$SomaticSNVs\t$SNVsDEPTH\t$ultFiltSNVsDEPTH\t$SomaticSNVsMarkersDEPTH\t$SomaticSNVsDEPTH\t$SNVsAFREQ\t$ultFiltSNVsAFREQ\t$SomaticSNVsMarkersAFREQ\t$SomaticSNVsAFREQ\t$SNVsPFRED\t$ultFiltSNVsPFRED\t$SomaticSNVsMarkersPFRED\t$SomaticSNVsPFRED" >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/Summary/BasicStatistics/deconSNV-and-P-Summary.txt
done
wait

echo 'Basic Statistics for Estimate of P Compiled!!!'
