#!/bin/bash
#SBATCH --job-name="growth.maps"
#SBATCH -p high
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH -o /homedtic/malenya/Simulations/Atlas/GIT-outs/%x-%j.out
#SBATCH -e /homedtic/malenya/Simulations/Atlas/GIT-outs/%x-%j.err
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
module load GCC/10.2.0
module load GCCcore/10.2.0
cd /homedtic/malenya/Simulations/Atlas/Regional_batch
export TEST=1
export IDE=125 					#100 = Gholipour Atlas. 122, 125, 128, 131.. are meshes generated from STA22, STA25, STA28, STA31 MRI Images.

export MESH_FILE="/homedtic/malenya/Simulations/Atlas/Mesh/STA25_1.2M.mesh"					# Path to mesh file (neutral format)
export LABEL_FILE="/homedtic/malenya/Simulations/Atlas/Mesh/STA25_1.2M_INDICATOR.txt"		# NO EFFECT NOW. In Brain Multiphysics paper here I added the label files. Now in the cpp file it has no effect.
export NODES_FILE="/homedtic/malenya/Simulations/Atlas/Mesh/STA25_1.2M_INDICATOR.txt"			# Cortical thickness values: 1 = Grey matter; 2 = white matter
export RESULTS_DIR="/homedtic/malenya/Simulations/Atlas/GIT-Results"						# Path to results folder (Global folder of the whole study. Subfolders with "NAME" will be created inside))
export MAP="/homedtic/malenya/Simulations/Atlas/Maps/STA2537_MAP.txt"						# Path to txt file with "distance maps" information. 
export NAME="GITHUB-STA25.MAP2537."$TEST															# Name of folder with results, inside the Results folder

export GA=25
export H=1.48

export GROWTH_RELATIVE=1.829 	# = sqrt(8)-1. Value from Tallinen, Xiaoyu.. 
export GROWTH_FACTOR=0			# Factor defined as multiplicative factor of whole brain (to increase normal growth of all tets). Here has no effect.

export HEMI=2;  				# 2 = Both Hemispheres. Some tests were performed with only 1 hemisphere but folding was not "OK". We can discuss this soon.
export RL="B";  				# B = Both Hemispheres.RL = "R", RL="L" or RL="B" for right, left or both hemis.


g++ -Ofast -fopenmp growthgithub.cpp vema.cpp eig3.cpp -o $TEST".growthmap"
./$TEST".growthmap" $H $MESH_FILE $LABEL_FILE $NODES_FILE $RESULTS_DIR $MAP $IDE $NAME $TEST $GA $GROWTH_RELATIVE $GROWTH_FACTOR $HEMI $RL 
