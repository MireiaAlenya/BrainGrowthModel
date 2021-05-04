#!/bin/bash
#SBATCH --job-name="20.h3"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=2G
#SBATCH -p high
#SBATCH -o /homedtic/malenya/CMMBE/outs.new/%x-%j.out
#SBATCH -e /homedtic/malenya/CMMBE/outs.new/%x-%j.err
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
module load GCC/6.3.0-2.27
module load GCCcore/6.3.0
cd /homedtic/malenya/CMMBE/Batch.new
g++ -Ofast -fopenmp 020.cpp vema.cpp eig3.cpp -o 20.h3
./20.h3 "3"
