#!/bin/bash
#$ -cwd

#...............................................................................#

#  System      | Mating             | Contributions    | Couples  | Limitation of K  | Selfing  |
# -------------|--------------------|------------------|----------|------------------|----------|
#  (0) RC      | Random             | Random           | Monogamy | Yes              | No       |
#  (1) EC      | Random or avoid FS | Equalization (2) | Monogamy | Yes              | No       |
#  (2) CM      | Circular           | Equalization (2) | Monogamy | Yes              | No       |
#  (3) PFS     | Random (% FS*)     | Random           | Monogamy | Yes              | No       |
#  (4) RCpol   | Random             | Random           | Polygamy | No               | Optional |

# *FS = full-sib mating
#...............................................................................#

rm script_extgwas_PURGd2.1.sh.*

################################# ARGUMENTS ###############################

#script_extgwas_PURGd2.1.sh <d> <NMAX> <TYPE> <PFS> <SEED> <FV> <MFEC> 
#To be used only in a machine with /state/partition1 directory

#Check number of arguments
if [ $# -ne 7 ]  
then
	echo "Usage: $0 <d> <NMAX> <TYPE> <PFS> <SEED> <FV> <MFEC>" 
	exit 1
fi

#Set arguments
d=$1
NMAX=$2		#population size
TYPE=$3		#mating system (R=0; EC=1; CM=2; PFS=3; RCpol=4)
PFS=$4		#proportion of full-sib matings
SEED=$5
FV=$6		#fec_via (0:equal, 1:1/3-2/3, 2:only viability; 3:only fecundity)
MFEC=$7         #k

################################ VARIABLES ################################

#Number of lines
NLINES=10000

#Parameters
LAMB=0.2
NCRO=300
AVEs=0.2
AVEh=0.283
BETA=0.33
LAMBL=0

#PURGd
Edelta=1	#0:estimation of delta from ind.Fa=0, 1:at the same time as any other parameter (default)

############################## DIRECTORIES ################################

#Working directory
WDIR=$PWD 

#Output directory
if [[ $TYPE == 0 ]]; then case="RC" 
elif [[ $TYPE == 1 ]]; then case="EC"
elif [[ $TYPE == 2 ]]; then case="CM"
elif [[ $TYPE == 3 ]]; then case="PFS.$PFS"
else case="RCpol"
fi

if [[ $FV == 0 ]]; then fit="fec&via"
elif [[ $FV == 1 ]]; then fit="fec&via_13_23"
elif [[ $FV == 2 ]]; then fit="via"
else fit="fec"
fi

mkdir -p $WDIR/extgwas_PURGd_Results/ch_L$LAMB.K$NCRO.s$AVEs.h$AVEh.N$NMAX/$case.$fit.mfec$MFEC
DIR="extgwas_PURGd_Results/ch_L$LAMB.K$NCRO.s$AVEs.h$AVEh.N$NMAX/$case.$fit.mfec$MFEC"

mkdir $WDIR/$DIR/GENFILES
mkdir $WDIR/$DIR/PURGd
mkdir $WDIR/$DIR/PEDFILES
if [[ $FV == 0 ]] || [[ $FV == 1 ]]; then 
	mkdir $WDIR/$DIR/PEDFILES/fecundity
	mkdir $WDIR/$DIR/PEDFILES/viability
fi

#Scratch directory
mkdir -p /state/partition1/noeliaEXT$d/$SLURM_JOBID/
SCDIR="/state/partition1/noeliaEXT$d" 

######################## TRANSFER OF FILES TO SCRATCH #####################

#Copy all files in scratch directory
cp seedfile $SCDIR/$SLURM_JOBID/
cp extinctiongwas_ms $SCDIR/$SLURM_JOBID/
cp POPFILE_L$LAMB.K$NCRO.s$AVEs.h$AVEh $SCDIR/$SLURM_JOBID/popfile
cp DATAFILE_L$LAMB.K$NCRO.s$AVEs.h$AVEh $SCDIR/$SLURM_JOBID/datafile
cp $WDIR/PURGd-2.1/PURGd $SCDIR/$SLURM_JOBID/
cp $WDIR/PURGd-2.1/settings.txt $SCDIR/$SLURM_JOBID/

#File with information of node and directory
touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`

#Move to scratch directory
cd $SCDIR/$SLURM_JOBID

###########################################################################

echo "Pedigree file,Analysis,d coefficient,RSS,AICc,RL,Chi2,p-value (Chi2),p-value (bootstrap),W0,SE(W0),b(g)" >> pedfile_ip.csv
if [[ "$FV" != 2 ]] ; then
echo "Pedigree file,Analysis,d coefficient,RSS,AICc,RL,Chi2,p-value (Chi2),p-value (bootstrap),W0,SE(W0),b(g)" >> fec_pedfile_ip.csv
echo "Pedigree file,Analysis,d coefficient,RSS,AICc,RL,Chi2,p-value (Chi2),p-value (bootstrap),W0,SE(W0),b(g)" >> via_pedfile_ip.csv
fi

START=$(date +%s)
for ((i=1; i<= $NLINES ; i++)); do
seed=$(( $SEED + $i - 1 ))

############# EXTINCTIONGWAS #############
 
time ./extinctiongwas_ms>>out<<@
0
$seed
$NMAX	NMAX
$MFEC	MAXFEC
$TYPE	type (R=0; EC=1; CM=2; PFS=3; RCpol=4)
1	Avoidance of full-sib matings (EC)(0 no, 1 yes)
0	Selfing in case of polygamy (0 no, 1 yes)
$PFS	PFS (For PFS, proportion of partial full-sib mating; 0=random)
99	L (morgans; 99=FreeRec)
$NCRO	NCRO (max 2000)(Neu=Ncro)
30	NLOCI (2-30)
$LAMB	Lambda_s
$LAMBL	Lambda_L
$BETA	beta_s
$AVEs	ave |s|
1	dom (0 constant, 1 CK94)
$AVEh	ave h_s
0.0	VE
0	relaxation factor (0: no, 1:yes)
0	neutral (0: no, 1:yes)
0	scaling (0: no, 1:yes)
$FV	fec_via (0: equal, 1: 1/3-2/3, 2: viability, 3: fecundity)
1	Fecundity as a maternal component (0: no; 1: yes)
1	Pedfile information (0: survivors; 1: zygotes)
$NMAX	generations
99	Change type to: RC(0), EC(1), CM(2), PFS(3), without change(99)
99	Generation of change
1	replicates
@

mv pedfile.dat pedfile${i}.csv
rm pedfile.dat
if [[ "$FV" != 2 ]] ; then
mv fec_pedfile.dat fec_pedfile${i}.csv
mv via_pedfile.dat via_pedfile${i}.csv
rm fec_pedfile.dat
rm via_pedfile.dat
fi

cp -r $SCDIR/$SLURM_JOBID/genfile.dat $WDIR/$DIR/GENFILES/GENFILE${i}
rm genfile.dat

################## PURGd #################

if [[ "$Edelta" == 1 ]] ; then	./PURGd --seed=1234 pedfile${i}.csv
else	./PURGd --seed=1234 --delta=s pedfile${i}.csv
fi

cp -r $SCDIR/$SLURM_JOBID/pedfile${i}.csv $WDIR/$DIR/PEDFILES/
rm pedfile${i}.csv
tail -n +2 pedfile${i}_ip.csv >> pedfile_ip.csv
rm pedfile${i}_ip.csv
cp -r $SCDIR/$SLURM_JOBID/pedfile_ip.csv $WDIR/$DIR/PURGd/

if [[ "$FV" < 2 ]] ; then
	if [[ "$Edelta" == 1 ]] ; then	./PURGd --seed=1234 fec_pedfile${i}.csv
	else	./PURGd --seed=1234 --delta=s fec_pedfile${i}.csv
	fi
	if [[ "$Edelta" == 1 ]] ; then	./PURGd --seed=1234 via_pedfile${i}.csv
	else	./PURGd --seed=1234 --delta=s via_pedfile${i}.csv
	fi

	cp -r $SCDIR/$SLURM_JOBID/fec_pedfile${i}.csv $WDIR/$DIR/PEDFILES/fecundity/
	cp -r $SCDIR/$SLURM_JOBID/via_pedfile${i}.csv $WDIR/$DIR/PEDFILES/viability/
	rm fec_pedfile${i}.csv
	rm via_pedfile${i}.csv

	tail -n +2 fec_pedfile${i}_ip.csv >> fec_pedfile_ip.csv
	tail -n +2 via_pedfile${i}_ip.csv >> via_pedfile_ip.csv
	rm fec_pedfile${i}_ip.csv
	rm via_pedfile${i}_ip.csv
	cp -r $SCDIR/$SLURM_JOBID/fec_pedfile_ip.csv $WDIR/$DIR/PURGd/
	cp -r $SCDIR/$SLURM_JOBID/via_pedfile_ip.csv $WDIR/$DIR/PURGd/
fi

echo "Line${i}" >> timefile

done

#########################################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "extinctiongwas and PURGd (E[delta] $Edelta) took $DIFF seconds\n" >> timefile

################# TRANSFER OF OTHER FILES TO DIRECTORY ####################

cp -r $SCDIR/$SLURM_JOBID/seedfile $WDIR
cp -r $SCDIR/$SLURM_JOBID/timefile $WDIR/$DIR/
cp -r $SCDIR/$SLURM_JOBID/dfilename*.dat $WDIR/$DIR/
cp -r $SCDIR/$SLURM_JOBID/out $WDIR/$DIR/

###########################################################################
######################### CLEANING OF SCRATCH #############################

rm -r $SCDIR/$SLURM_JOBID/
rm $WDIR/$SLURM_JOBID.*
