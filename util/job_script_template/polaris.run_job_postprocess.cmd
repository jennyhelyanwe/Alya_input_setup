#!/bin/bash
#PBS -l select=<<computational_nodes>>:ncpus=<<computational_cores>>
#PBS -l walltime=0<<job_time>>:00:00
#PBS -q prod
#PBS -l filesystems=grand
#PBS -A 11568
#PBS -N AlyaHeart

module load PrgEnv-gnu
module load gcc/11.2.0
cd $PBS_O_WORKDIR
mpiexec -n <<computational_cores>> python3.6 alya2csvensight_mpi.py
