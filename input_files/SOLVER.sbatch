#!/bin/bash
#SBATCH --job-name=SOLVER 
#SBATCH --nodes=4
#SBATCH --ntasks=72
#SBATCH --time=00:30:00
#SBATCH --output=ro.txt 
#SBATCH --error=re.txt 
#SBATCH --partition=cpu_med

# redistribute allocated resources
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
np=$[${numnodes}*${mpi_tasks_per_node}]


# loading modules
source $HOME/load_sem3d_ruche_modules.sh

export buildfolder=/gpfs/users/gattif/SEMq/build

SOLVER="$buildfolder/SEM3D/sem3d.exe"

srun -n 72 $SOLVER > output.run
