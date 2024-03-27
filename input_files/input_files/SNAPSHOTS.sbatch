#!/bin/bash
#SBATCH --job-name=SOLVER_BACKWARD
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=00:30:00
#SBATCH --output=ro.txt 
#SBATCH --error=re.txt 
#SBATCH --partition=cpu_prod

output = $(python3 ./pysem/parse_h5_snapshots.py)
echo output > g_lamb_mu.txt 
