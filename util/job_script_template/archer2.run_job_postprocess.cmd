#!/bin/bash
#SBATCH --acount=e769-jennywang
#SBATCH --job-name=heartpp
#SBATCH --time=00:20:00
#SBATCH --nodes=<<computational_nodes>>
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --qos=short
#SBATCH --partition=standard
#SBATCH --export-none

module load PrgEnv-gnu
module load gcc/11.2.0
module load cray-python/3.9.13.1

srun -n <<computational_cores>> python alya2csvensight_mpi.py
