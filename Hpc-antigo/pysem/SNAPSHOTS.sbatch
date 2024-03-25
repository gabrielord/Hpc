#!/bin/bash
#SBATCH --job-name=SNAPSHOTS 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --output=so.txt 
#SBATCH --error=se.txt 
#SBATCH --partition=cpu_prod

mpirun -n $SLURM_NTASKS python3 ./pysem/parse_h5_snapshots.py