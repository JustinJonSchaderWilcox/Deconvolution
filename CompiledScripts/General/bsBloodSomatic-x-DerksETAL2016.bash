srun --pty -n 1 -c 12 --mem 18GB -t 1:23:00 /bin/bash
source activate Genomics

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/SomaticValidate
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/SomaticValidate
for P in 0 10 20 30 40 50 60 70 80 90 100
do
bcftools isec -o /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/SomaticValidate/$P-bsBlood-wo-iBlood.vcf -n 1 -w 1 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Mixture/EstimateP/BloodSomatic/BloodBrainMix.$P.snp.ultfFilt-sym.blood-somatic.vcf.gz /work/mjuswilc/Domain/GreatTits/General/Annotations/Abel/Variants/MPILEUP/VCFs/SNVs/Abel-GCF_001522545.3_Parus_major1.1_genomic.snvs-any.vcf.gz &
done
wait
for P in 0 10 20 30 40 50 60 70 80 90 100
do
fpsSNVs=`grep -v -c "#" /work/mjuswilc/Domain/GreatTits/ultDeconvolution/ManuscriptMisc/SomaticValidate/$P-bsBlood-wo-iBlood.vcf`
echo -e "$P\t$fpsSNVs"
done
