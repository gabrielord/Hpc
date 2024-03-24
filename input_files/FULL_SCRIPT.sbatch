#!/bin/bash
#SBATCH --job-name=FULL_PROJECT_SEISM
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=00:30:00
#SBATCH --output=output.txt
#SBATCH --error=error.txt 
#SBATCH --partition=cpu_prod

echo "*******************************************************************************************"
echo "*******************************************************************************************"
echo "********************* Project of determination of the Lamé parameters *********************"
echo "*********************      of a material using the adjoint method     *********************"
echo "*******************************************************************************************"
echo "*******************************************************************************************"
echo -e "\n\t- Course: High Performance Computing"
echo -e "\t- Authors: Alexandre FAURE, Paul HUGUET, Gabriel LIMA, Vitor OPSFELDER."
echo -e "\n-----------------------------------"

# Import des modules
module load py-mpi4py/3.1.4/gcc-12.3.0-openmpi
module load py-scipy/1.11.2/gcc-13.1.0
module load py-h5py/3.8.0/gcc-12.3.0-openmpi
module load py-numpy/1.24.3/gcc-13.1.0
module load py-matplotlib/3.7.1/gcc-12.3.0
module load sem3d/2024.02.22/intel-2021.9.0-intel-oneapi-mpi

echo "Modules loaded successfully!"
echo -e "-----------------------------------"

# Launch the python script
python3 full_script.py