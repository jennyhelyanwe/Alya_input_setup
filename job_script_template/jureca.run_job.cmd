#!/bin/bash -l
#SBATCH --account=icei-prace-2022-0003
#SBATCH --job-name="<<job_name>>"
#SBATCH --time=0<<job_time>>:00:00
#SBATCH --nodes=<<computational_nodes>>
#SBATCH --ntasks-per-node=<<tasks_per_node>>
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --partition=<<job_type>>

alyaenv
srun -n <<computational_cores>> python main.py