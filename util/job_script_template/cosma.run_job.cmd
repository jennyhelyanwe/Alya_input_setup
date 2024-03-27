#!/bin/bash -l
#SBATCH --account="dp287"
#SBATCH --job-name="<<job_name>>"
#SBATCH --time=0<<job_time>>:00:00
#SBATCH --nodes=<<computational_nodes>>
#SBATCH --ntasks-per-node=<<tasks_per_node>>
#SBATCH --cpus-per-task=1
#SBATCH --partition=<<job_type>>

module load intel_comp/2020-update2
module load intel_mpi/2020-update2
mpirun -n <<computational_cores>> python <<python_script_name>>