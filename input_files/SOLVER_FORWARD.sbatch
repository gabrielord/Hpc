#!/bin/bash
#SBATCH --job-name=SOLVER 
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=00:30:00
#SBATCH --output=ro.txt 
#SBATCH --error=re.txt 
#SBATCH --partition=cpu_prod

module load sem3d/2024.02.22/intel-2021.9.0-intel-oneapi-mpi

mpirun -n $SLURM_NTASKS sem3d.exe > output.solver

# Rename traces and res to traces_backward and res_backward
mv traces traces_forward
mv res res_forward