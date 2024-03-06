#!/bin/bash
#SBATCH --job-name=MESHER
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50GB
#SBATCH --time=00:30:00
#SBATCH --output=mo.txt 
#SBATCH --error=me.txt
#SBATCH --partition=cpu_short

# loading modules
source $HOME/load_sem3d_ruche_modules.sh

export buildfolder=/gpfs/users/gattif/SEMq/build
MESHER="$buildfolder/MESH/mesher"

$MESHER < mesh.input > output.mesh
mkdir -p sem; mv mesh4spec.* ./sem/
