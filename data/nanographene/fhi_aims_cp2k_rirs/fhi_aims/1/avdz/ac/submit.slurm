#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH --exclusive
#SBATCH --ntasks-per-node=192
#SBATCH --cpus-per-task=1
#SBATCH -N 1
#SBATCH -J "rirs"
#SBATCH -p normal
#SBATCH --account=hpc-prf-metdyn

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

ulimit -s unlimited

module use /pc2/groups/hpc-prf-metdyn/metdyn07_Ritaj/modulefiles
module load aims/dev

mpirun -np $SLURM_NTASKS /pc2/groups/hpc-prf-metdyn/metdyn07_Ritaj/bin/aims.x > output.out
