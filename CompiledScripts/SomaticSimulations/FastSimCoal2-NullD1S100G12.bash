#!/bin/bash
#SBATCH --time=26:40:00
#SBATCH --mem=58GB
#SBATCH -p long
# Set number of nodes to run
#SBATCH --nodes=1
# Set number of tasks to run
#SBATCH --ntasks=20
# Set number of cores per task (default is 1)
#SBATCH --cpus-per-task=1
# Output and error files
#SBATCH -o job.%J.out


#Directory Structure
rm -r /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Input
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Output
mkdir /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Records

function fastball {
#Parameters
DEMES=1
cPOPSIZE=819200
sHAPLOTYPES=100
cGROWTH=-0.69
cMIGRATION=0
hNUM=1
BASEPAIRS=1000000000
RECRATE=0
dMU=$mu
RECRATE=0
TRANSITION=0.6757415

#Parameter File
n=D1S100G12MU"$mu"
echo '//Number of population samples (demes)' > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Input/SomaticMutation$n.fsc27.par
echo $DEMES >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Input/SomaticMutation$n.fsc27.par
echo '//Population effective sizes (number of genes)' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Input/SomaticMutation$n.fsc27.par
echo $cPOPSIZE >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Input/SomaticMutation$n.fsc27.par
echo '//Sample sizes' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Input/SomaticMutation$n.fsc27.par
echo $sHAPLOTYPES >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Input/SomaticMutation$n.fsc27.par
echo '//Growth rates  : negative growth implies population expansion' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Input/SomaticMutation$n.fsc27.par
echo $cGROWTH >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Input/SomaticMutation$n.fsc27.par
echo '//Number of migration matrices : 0 implies no migration between demes' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Input/SomaticMutation$n.fsc27.par
echo $cMIGRATION >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Input/SomaticMutation$n.fsc27.par
echo '//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Input/SomaticMutation$n.fsc27.par
echo "$hNUM historical event" >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Input/SomaticMutation$n.fsc27.par
#echo "$hTIME1 $hSOURCE1 $hSINK1 $hMIGRANTS1 $hSIZE1 $hGROWTH1 $hMIGRATION1"
hTIME=12
hSIZE=200
echo "$hTIME 0 0 0 $hSIZE 0 0" >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Input/SomaticMutation$n.fsc27.par
echo '//Number of independent loci [chromosome]' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Input/SomaticMutation$n.fsc27.par
echo '1 0' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Input/SomaticMutation$n.fsc27.par
echo '//Per chromosome: Number of linkage blocks' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Input/SomaticMutation$n.fsc27.par
echo 1 >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Input/SomaticMutation$n.fsc27.par
echo '//per Block: data type, num loci, rec. rate and mut rate + optional parameters' >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Input/SomaticMutation$n.fsc27.par
echo "DNA $BASEPAIRS $RECRATE $dMU $TRANSITION" >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Input/SomaticMutation$n.fsc27.par
cd /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/
/home/mjuswilc/Programs/fsc27_linux64/fsc27093 --seed 23 -i /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Input/SomaticMutation$n.fsc27.par -n 100 -k 1000000000 -G -g -p -x -c 1
}

function frequency {
n=D1S100G12MU"$mu"
echo -e 'Generation\tChrom\tPositon\tOriginal\tDerived\tAlleFreq' > /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Output/FastSimCoalSimDepth-$n.txt
COUNT=0
for g in $(ls /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/SomaticMutation$n.fsc27/SomaticMutation$n.fsc27_1_*.gen | sort -V)
do
COUNT=`expr $COUNT + 1`
awk -v gen=$COUNT 'BEGIN{OFS=FS="\t"}NR==1{next}{AC=$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22+$23+$24+$25+$26+$27+$28+$29+$30+$31+$32+$33+$34+$35+$36+$37+$38+$39+$40+$41+$42+$43+$44+$45+$46+$47+$48+$49+$50+$51+$52+$53+$54}AC>0{print gen OFS $1 OFS $2 OFS $3 OFS $4 OFS AC/(2*50)}' $g >> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Output/FastSimCoalSimAF-$n.txt
done
}

for mu in 0.0000000023 0.0000000046 0.0000000092 0.0000000138 0.0000000184 0.0000000230 0.0000000276 0.0000000322 0.0000000368 0.0000000414 0.0000000460
do
fastball &> /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Records/run_output.$mu.txt &
JOBS=`jobs | wc | awk '{print $1}'`
while [ $JOBS -ge 4 ]
do
sleep 3s
JOBS=`jobs | wc | awk '{print $1}'`
done
done
wait

for mu in 0.0000000023 0.0000000046 0.0000000092 0.0000000138 0.0000000184 0.0000000230 0.0000000276 0.0000000322 0.0000000368 0.0000000414 0.0000000460
do
frequency &
done
wait

for mu in 0.0000000023 0.0000000046 0.0000000092 0.0000000138 0.0000000184 0.0000000230 0.0000000276 0.0000000322 0.0000000368 0.0000000414 0.0000000460
do
n=D1S100G12MU"$mu"
gzip -9 /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/Output/FastSimCoalSimAF-$n.txt &
done
wait

for mu in 0.0000000023 0.0000000046 0.0000000092 0.0000000138 0.0000000184 0.0000000230 0.0000000276 0.0000000322 0.0000000368 0.0000000414 0.0000000460
do
n=D1S100G12MU"$mu"
for g in $(ls /work/mjuswilc/Domain/GreatTits/ultDeconvolution/Simulations/FastSimCoal2/D1S100G12/SomaticMutation$n.fsc27/SomaticMutation$n.fsc27_1_*.gen)
do
gzip -9 $g &
JOBS=`jobs | wc | awk '{print $1}'`
while [ $JOBS -ge 18 ]
do
sleep 3s
JOBS=`jobs | wc | awk '{print $1}'`
done
done
done
wait

echo 'Simulation Job Complete!!!'
