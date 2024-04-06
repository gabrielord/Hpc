"""
Main script to determine the Lamé parameters of the geological layer.
"""
import os
import sys
import numpy as np
import subprocess
from math import floor
from datetime import datetime
import json

### Constants
N_iter = 10
alpha_lamb = 1
alpha_mu = 1
c1 = 1e-4
xi = 0.5
tol = 1e-4

def write_output(txt):
    """
    Function to write in the output file.
    """
    with open("output.txt", "a") as output_file:
        output_file.write(txt+"\n")

def get_current_time():
    """
    Return the current date with format DD/MM/YYYY à HHhMMminSS
    """
    current_datetime = datetime.now()
    return current_datetime.strftime("%d/%m/%Y à %Hh%Mmin%S")


def create_useful_folders():
    """
    Create the folders that are missing.
    """
    folders = ["output_files", "monitors_misfit", "sem", "stats"]
    for folder in folders:
        os.makedirs(folder, exist_ok=True)
    
    sem3d_folders = ["traces", "res"]
    for folder in sem3d_folders:
        os.makedirs(f"forward/{folder}", exist_ok=True)
        os.makedirs(f"backward/{folder}", exist_ok=True)

def copy_material_files():
    """
    Copy the material files in the material folder to keep a version
    of the initial state of the layer.
    """
    os.makedirs("material", exist_ok=True)
    # On copie tous les fichiers de la forme example_*.h5 et example_*.xml
    # au format *.h5 et *.xml dans le dossier material
    os.system("""for file in material/initial_*; do cp "$file" "$(echo $file | sed 's/initial_//')"; done""")

def main():
    """Main function of the script."""
    ### Initialization
    create_useful_folders()
    copy_material_files()

    ### Constantes
    dict_sizePb_nbProc = {
        1:[8,16,32],
        4:[4,8,16,32],
        8:[1,2,4,8,16,32]
    }

    for size_pb in dict_sizePb_nbProc:
        for nb_proc in dict_sizePb_nbProc[size_pb]:
            write_output(f"=============> Résolution d'un problème 1/{size_pb} sur {nb_proc} process")
            ### Redéfinition des fichiers
            os.system(f"cp material_perf_{size_pb}.spec material.spec")
            os.system(f"cp mat_perf_{size_pb}.dat mat.dat")
            with open("mesh.input", "w") as f:
                f.write(f"{nb_proc}\n1")

            ### We solve the forward problem
            write_output(f"\n\n1. [{get_current_time()}] Solving the forward problem...")

            # Launch the mesher
            write_output(f"\n==> [{get_current_time()}] Launching the mesher...")
            write_output("\t- Update the input.spec file...")
            os.system("cp input_forward.spec input.spec")
            write_output("\t- Launch the mesher...")
            os.system(f"mesher < mesh.input > output_files/forward_test.mesher")
            write_output("\t- Move the mesh files...")
            os.system("mv mesh4spec.* ./sem/")

            # Launch the solver
            write_output(f"\n==> [{get_current_time()}] Launching the solver...")
            write_output("\t- Launch the solver...")
            cmd = f"mpirun -np {nb_proc} " + \
                f"-map-by ppr:1:core:PE=1 " + \
                f"sem3d.exe"
            with open(f"output_files/forward_test.solver", "w") as f:
                subprocess.run(cmd, shell=True, stdout=f)
            
            write_output("\t- Move the output files...")
            os.system(f"mv stat.log output_files/stat_forward_test_{nb_proc}_{size_pb}.log")


if __name__ == "__main__":
    ### Call the main function
    main()