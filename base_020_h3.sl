#!/bin/bash
#SBATCH --job-name="20.h3"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=2G
#SBATCH -p high
#SBATCH -o /homedtic/malenya/BrainGrowthModel/outs/%x-%j.out
#SBATCH -e /homedtic/malenya/BrainGrowthModel/outs/%x-%j.err
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
module load GCC/6.3.0-2.27
module load GCCcore/6.3.0
cd /homedtic/malenya/BrainGrowthModel/Batch

export H=3
export IDE=20
export GROWTH_RELATIVE=1.829
export MESH_FILE="../../Brains_LargerStudy/Meshes/F-s020-w26-820k-gt.mesh"
#export LABEL_FILE="../../Meshes/fet-020-820k-labelfile.txt"
export PATH_DIR="../Results/Growth_Slope"
export NAME="GS_"$IDE

export gaFile="/homedtic/malenya/MeshesBrains/gestational_weeks_complete.csv"
declare -A FGA NGA BGA GENDER ONSET
INPUT=$gaFile
OLDIFS=$IFS
IFS=,
[ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
while read id_data fetal_ga neo_ga birth gender onset
do
FGA[$id_data]=$fetal_ga
NGA[$id_data]=$neo_ga
BIRTH[$id_data]=$birth
GENDER[$id_data]=$gender
ONSET[$id_data]=$onset
done < $INPUT
IFS=$OLDIFS
GA=${FGA[$IDE]}
nga=${NGA[$IDE]}
birth=${BIRTH[$IDE]}
gender=${GENDER[$IDE]}
onset=${ONSET[$IDE]}

g++ -Ofast -fopenmp base_0$IDE_h$H.cpp vema.cpp eig3.cpp -o $IDE.h$H
./$IDE.h$H $H $GROWTH_RELATIVE $MESH_FILE $IDE $PATH_DIR $NAME $GA $nga $birth $gender $onset
