#!/bin/bash
#SBATCH --job-name=MESHER
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=mo.txt 
#SBATCH --error=me.txt
#SBATCH --partition=cpu_prod

module load sem3d/2024.02.22/intel-2021.9.0-intel-oneapi-mpi

## On copie le fichier input_forward.spec dans input.spec
cp input_forward.spec input.spec

## On lance le mesher
mesher < mesh.input > output.mesher
mkdir -p sem; mv mesh4spec.* ./sem/
