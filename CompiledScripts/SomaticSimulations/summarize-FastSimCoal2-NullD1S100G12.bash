#!/bin/bash
#SBATCH --time=1-23:54:00
#SBATCH --mem=8GB
#SBATCH -p long
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=1
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=34
# Output and error files
#SBATCH -o job.%J.out

rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP

function summary_execution {
n=D1S100G12MU"$mu"
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Output/FastSimCoalSimAF-$n.txt.gz > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimAF-$n.txt &
zcat /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Output/FastSimCoalSimAF-$n.txt.gz | awk 'BEGIN{OFS=FS="\t"}$6>0.5{$6=1-$6}{print $0}' > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimFoldedAF-$n.txt &
wait
g=1
grep -P "^$g\s" /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimAF-$n.txt > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimAF-$n.$g.txt &
grep -P "^$g\s" /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimFoldedAF-$n.txt > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimFoldedAF-$n.$g.txt &
wait

echo -e "Simulation\tMu\tIteration\tMutations\tafMEAN\tafMEDIAN\tafMODE\tafMODEprop\tafMEDIANAbsDEV\tFOLDEDafMEAN\tFOLDEDafMEDIAN\tFOLDEDafMODE\tFOLDEDafMODEprop\tFOLDEDafMEDIANAbsDEV\tcut25PROP" > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/FastSimCoalSomaticAF-$n.summary.txt
for g in $(seq 1 100)
do

if [ $g -lt 100 ]
then
NEXT=`expr $g + 1`
grep -P "^$NEXT\s" /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimAF-$n.txt > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimAF-$n.$NEXT.txt &
grep -P "^$NEXT\s" /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimFoldedAF-$n.txt > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimFoldedAF-$n.$NEXT.txt &
fi

MUTATIONS=`wc -l /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimAF-$n.$g.txt | awk '{print $1}'`
afMEAN=`awk 'BEGIN{OFS=FS="\t"; afreq=0}{afreq+=$6}END{print afreq/NR}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimAF-$n.$g.txt`
afMEDIAN=`awk 'BEGIN{OFS=FS="\t"}{print $6}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimAF-$n.$g.txt | sort -n | awk -v mutations=$MUTATIONS 'BEGIN{median=int(mutations/2)}NR==median{print $0; exit}'`
afMODE=`awk 'BEGIN{OFS=FS="\t"}{print $6}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimAF-$n.$g.txt | sort -n | uniq -c | sort -k1,1rn | awk 'BEGIN{OFS="\t"}{print $2; exit}'`
afMODEprop=`awk 'BEGIN{OFS=FS="\t"}{print $6}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimAF-$n.$g.txt | sort -n | uniq -c | sort -k1,1rn | awk -v mutations=$MUTATIONS 'BEGIN{OFS="\t"}{print $1/mutations; exit}'`
afMEDIANAbsDEV=`awk -v median=$afMEDIAN 'BEGIN{OFS=FS="\t"; dev=0}{dev+=sqrt(($6-median)^2)}END{print dev/NR}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimAF-$n.$g.txt`

FOLDEDafMEAN=`awk 'BEGIN{OFS=FS="\t"; afreq=0}{afreq+=$6}END{print afreq/NR}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimFoldedAF-$n.$g.txt`
FOLDEDafMEDIAN=`awk 'BEGIN{OFS=FS="\t"}{print $6}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimFoldedAF-$n.$g.txt | sort -n | awk -v mutations=$MUTATIONS 'BEGIN{median=int(mutations/2)}NR==median{print $0; exit}'`
FOLDEDafMODE=`awk 'BEGIN{OFS=FS="\t"}{print $6}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimFoldedAF-$n.$g.txt | sort -n | uniq -c | sort -k1,1rn | awk 'BEGIN{OFS="\t"}{print $2; exit}'`
FOLDEDafMODEprop=`awk 'BEGIN{OFS=FS="\t"}{print $6}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimFoldedAF-$n.$g.txt | sort -n | uniq -c | sort -k1,1rn | awk -v mutations=$MUTATIONS 'BEGIN{OFS="\t"}{print $1/mutations; exit}'`
FOLDEDafMEDIANAbsDEV=`awk -v median=$afMEDIAN 'BEGIN{OFS=FS="\t"; dev=0}{dev+=sqrt(($6-median)^2)}END{print dev/NR}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimFoldedAF-$n.$g.txt`

cut25PROP=`awk 'BEGIN{OFS=FS="\t"; count}$6<0.25{count++}END{print count/NR}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimFoldedAF-$n.$g.txt`

echo -e "$n\t$mu\t$g\t$MUTATIONS\t$afMEAN\t$afMEDIAN\t$afMODE\t$afMODEprop\t$afMEDIANAbsDEV\t$FOLDEDafMEAN\t$FOLDEDafMEDIAN\t$FOLDEDafMODE\t$FOLDEDafMODEprop\t$FOLDEDafMEDIANAbsDEV\t$cut25PROP" >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/FastSimCoalSomaticAF-$n.summary.txt
rm /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimAF-$n.$g.txt
rm /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP/FastSimCoalSimFoldedAF-$n.$g.txt
wait
done
}

for mu in 0.0000000023 0.0000000046 0.0000000092 0.0000000138 0.0000000184 0.0000000230 0.0000000276 0.0000000322 0.0000000368 0.0000000414 0.0000000460
do
summary_execution &
done
wait

echo -e "Simulation\tMu\tIteration\tMutations\tafMEAN\tafMEDIAN\tafMODE\tafMODEprop\tafMEDIANAbsDEV\tFOLDEDafMEAN\tFOLDEDafMEDIAN\tFOLDEDafMODE\tFOLDEDafMODEprop\tFOLDEDafMEDIANAbsDEV\tcut25PROP"  > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/FastSimCoalSomaticAF-D1S100G12MUall.summary.txt
for mu in 0.0000000023 0.0000000046 0.0000000092 0.0000000138 0.0000000184 0.0000000230 0.0000000276 0.0000000322 0.0000000368 0.0000000414 0.0000000460
do
n=D1S100G12MU"$mu"
awk 'BEGIN{OFS=FS="\t"}NR==1{next}{print $0}' /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/FastSimCoalSomaticAF-$n.summary.txt >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/FastSimCoalSomaticAF-D1S100G12MUall.summary.txt
done
rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/Summary/TEMP

echo 'Simulation Summary Job Complete!!!'
