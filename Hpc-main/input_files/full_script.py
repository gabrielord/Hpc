"""
Main script to determine the Lamé parameters of the geological layer.
"""
import os
import sys
import numpy as np
import subprocess
from math import floor
from datetime import datetime

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

def main():
    """Main function of the script."""
    ### Initialization
    create_useful_folders()

    ### Variables
    iter = 0
    J = None

    while iter < N_iter and (J is None or J > tol):
        write_output("\n" + 2*f"***********************{'*'*floor(np.log10(max(1, iter)))}\n"+ \
                              f"***** ITERATION {iter} *****\n" + \
                            2*f"***********************{'*'*floor(np.log10(max(1, iter)))}\n")
        if J is not None:
            write_output(f"Current cost function value: {J}")
        

        ### We solve the forward problem
        write_output(f"\n\n1. [{get_current_time()}] Solving the forward problem...")

        # Launch the mesher
        write_output(f"\n==> [{get_current_time()}] Launching the mesher...")
        write_output("\t- Update the input.spec file...")
        os.system("cp input_forward.spec input.spec")
        write_output("\t- Launch the mesher...")
        os.system(f"mesher < mesh.input > output_files/forward_{iter}.mesher")
        write_output("\t- Move the mesh files...")
        os.system("mv mesh4spec.* ./sem/")

        # Launch the solver
        write_output(f"\n==> [{get_current_time()}] Launching the solver...")
        write_output("\t- Launch the solver...")
        cmd = f"mpirun -np 32 " + \
              f"-map-by ppr:1:core:PE=1 " + \
              f"sem3d.exe"
        f = open(f"output_files/forward_{iter}.solver", "w")
        subprocess.run(cmd, shell=True, stdout=f)
        f.close()
        write_output("\t- Move the output files...")
        os.system("mv traces forward/traces")
        os.system("mv res forward/res")
        os.system(f"mv stat.log output_files/stat_forward_{iter}.log")



        ### We determine the misfit
        write_output(f"\n\n2. [{get_current_time()}] Determining the misfit...")
        # Generate the misfit files
        write_output(f"\n==> [{get_current_time()}] Generating the misfit files...")
        cmd = f"mpirun -np 1 " + \
              f"-map-by ppr:1:core:PE=1 " + \
              f"python3 ../pysem/compute_misfit.py"
        f = open(f"output_files/misfit_{iter}.output", "w")
        subprocess.run(cmd, shell=True, stdout=f)
        f.close()


        ### We solve the adjoint (backward) problem
        write_output(f"\n\n3. [{get_current_time()}] Solving the adjoint problem...")
        # Launch the mesher
        write_output(f"\n==> [{get_current_time()}] Launching the mesher...")
        write_output("\t- Update the input.spec file...")
        os.system("cp input_backward.spec input.spec")
        write_output("\t- Launch the mesher...")
        os.system(f"mesher < mesh.input > output_files/backward_{iter}.mesher")
        write_output("\t- Move the mesh files...")
        os.system("mv mesh4spec.* ./sem/")

        # Launch the solver
        write_output(f"\n==> [{get_current_time()}] Launching the solver...")
        write_output("\t- Launch the solver...")
        cmd = f"mpirun -np 32 " + \
                f"-map-by ppr:1:core:PE=1 " + \
                f"sem3d.exe"
        f = open(f"output_files/backward_{iter}.solver", "w")
        subprocess.run(cmd, shell=True, stdout=f)
        f.close()
        write_output("\t- Move the output files...")
        os.system("mv traces backward/traces")
        os.system("mv res backward/res")
        os.system(f"mv stat.log output_files/stat_backward_{iter}.log")

        ### We determine the gradient
        write_output(f"\n\n4. [{get_current_time()}] Determining the gradient...")
        # Calculate the gradient
        write_output(f"\n==> [{get_current_time()}] Calculating the gradient...")
        cmd = f"mpirun -np 1 " + \
              f"-map-by ppr:1:core:PE=1 " + \
              f"python3 ../pysem/gradient.py"
        f = open(f"output_files/gradient_{iter}.output", "w")
        subprocess.run(cmd, shell=True, stdout=f)
        f.close()

        # Calculate the search direction
        write_output(f"\n==> [{get_current_time()}] Calculating the search direction...")

        # Calculate the step lengths
        write_output(f"\n==> [{get_current_time()}] Calculating the step lengths...")

        ### We update the material parameters
        write_output(f"\n\n5. [{get_current_time()}] Updating the material parameters...")

        ### We increment the iteration counter
        iter += 10


if __name__ == "__main__":
    ### Call the main function
    main()