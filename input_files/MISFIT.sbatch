#!/bin/bash
#SBATCH --job-name=SOLVER 
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH --time=00:30:00
#SBATCH --output=ro.txt 
#SBATCH --error=re.txt 
#SBATCH --partition=cpu_prod

python ./pysem/compute_misfit.py 