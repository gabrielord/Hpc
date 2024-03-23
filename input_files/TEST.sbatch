#!/bin/bash
#SBATCH --job-name=SOLVER_BACKWARD
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --output=ro.txt 
#SBATCH --error=re.txt 
#SBATCH --partition=cpu_prod


sbatch  "./input_files/test1.sbatch"
sbatch --dependency=afterok:$J "./input_files/test2.sbatch
valeur=$(cat test.txt)
((valeur += 1))
echo valeur test.txt