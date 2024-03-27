"""
Main script to determine the Lamé parameters of the geological layer.
"""
import os
import sys
import numpy as np
import subprocess
from math import floor
from datetime import datetime
from ..pysem.parse_h5_snapshots import GetSnapshots
from mpi4py import MPI
from step_direction import L_bfgs, step_length, recursion_two_loop
from ..pysem.gradient import gradient_mu, gradient_lamb
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

    #Init varibale
    N_iter = 5
    alpha_lamb = 1
    alpha_mu = 1
    c1= 10**(-4)
    xi = 0.5

    #Update les nodes pour changer le matériaux
    MPI.Init()
    comm = MPI.COMM_WORLD               # Get communicator
    size = MPI.COMM_WORLD.Get_size()    # Get size of communicator
    rank = MPI.COMM_WORLD.Get_rank()    # Get the current rank
    hostname = MPI.Get_processor_name() # Get the hostname

    with open('cout.txt', 'r') as file:
        J_k_1 = file.readlines()[-1].strip()

    ### Variables
    iter = 0
    J = None

    while iter < N_iter and (J is None or J > tol):
        write_output("\n" + 2*f"***********************{'*'*floor(np.log10(max(1, iter)))}\n"+ \
                              f"***** ITERATION {iter} *****\n" + \
                            2*f"***********************{'*'*floor(np.log10(max(1, iter)))}\n")
        if J is not None:
            write_output(f"Current cost function value: {J}")
        

        '''### We solve the forward problem
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
        os.system("mv traces forward")
        os.system("mv res forward")
        os.system(f"mv stat.log output_files/stat_forward_{iter}.log")'''

        # Backtracking
        while iter < N_iter and (J_k is None or J_k_1 > J_k + c1*(alpha_lamb*np.dot(s_lamb, g_lamb)+ alpha_mu*np.dot(s_mu, g_mu)) ):
            write_output("\n" + 2*f"***********************{'*'*floor(np.log10(max(1, iter)))}\n"+ \
                                f"***** ITERATION {iter} *****\n" + \
                                2*f"***********************{'*'*floor(np.log10(max(1, iter)))}\n")
            if J_k_1 is not None:
                write_output(f"Current cost function value, J_k_1 and J_k: {(J_k_1,J_k)}")
            alpha_lamb *= xi
            alpha_mu *= xi
        
            snp_forward, snp_backward = GetSnapshots(comm,size,rank)
            nodes = snp_forward.dset['Nodes']
            MPI.Finalize()

            ### We update the material parameters
            write_output(f"\n\n5. [{get_current_time()}] Updating the material parameters...")

            ###Change parametes in material file 
            write_output(f"\n\n5. [{get_current_time()}] Updating the material parameters...")

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

            ### We determine the cost_function
            write_output(f"\n\n2. [{get_current_time()}] Determining the cost...")
            # Generate the misfit files
            write_output(f"\n==> [{get_current_time()}] Generating the cost files...")
            cmd = f"mpirun -np 1 " + \
                f"-map-by ppr:1:core:PE=1 " + \
                f"python3 ../pysem/costFonction.py"
            f = open(f"output_files/misfit_{iter}.output", "w")
            subprocess.run(cmd, shell=True, stdout=f)
            f.close()

            J_k = J_k_1
            with open('cout.txt', 'r') as file:
                J_k_1 = file.readlines()[-1].strip()
            iter +=10



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
        os.system("mv traces backward")
        os.system("mv res backward")
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

        m = 10
        if iter == 0:
            y_mu_stored = []
            y_lamb_stored = []
            s_mu_stored = []
            s_lamb_stored = []
            nabla_mu = gradient_mu(mu, lamb, g_reg_lamb, g_reg_mu, g_lamb_mis, g_mu_mis) # initial gradient 
            nabla_lamb = gradient_lamb(mu, lamb, g_reg_lamb, g_reg_mu, g_lamb_mis, g_mu_mis)
        
        grad_mu_old = nabla_mu[:]
        grad_lamb_old = nabla_lamb[:]

        s_mu = - recursion_two_loop(nabla_mu, np.array(s_mu_stored), np.array(y_mu_stored), m)
        s_lamb = - recursion_two_loop(nabla_lamb, np.array(s_lamb_stored), np.array(y_lamb_stored), m)

        s_mu_dict = {"s_mu": list(s_mu)}
        s_lamb_dict = {"s_lamb": list(s_lamb)}

        with open("s_mu_dict.json") as json_file:
            json.dump(s_mu_dict, json_file)

        with open("s_lamb_dict.json") as json_file:
            json.dump(s_lamb_dict, json_file)

        # Calculate the step length
        write_output(f"\n==> [{get_current_time()}] Calculating the step length...")

        alpha_lamb, alpha_mu = step_length(s_lamb, g_lamb, s_mu, g_mu, comm, rank, size)

        ### We update the material parameters
        step_mu = - alpha_mu * nabla_mu
        step_lamb = - alpha_lamb * nabla_lamb

        s_mu_stored.append(step_mu)
        s_lamb_stored.append(step_lamb)
        
        mu = mu + step_mu
        lamb = lamb + step_lamb

        nabla_mu = gradient_mu(mu, lamb, g_reg_lamb, g_reg_mu, g_lamb_mis, g_mu_mis) # initial gradient 
        nabla_lamb = gradient_lamb(mu, lamb, g_reg_lamb, g_reg_mu, g_lamb_mis, g_mu_mis)

        y_mu_stored.append(nabla_mu - grad_mu_old)
        y_lamb_stored.append(nabla_lamb - grad_lamb_old)

        write_output(f"\n\n5. [{get_current_time()}] Updating the material parameters...")
        cmd = f"mpirun -np 1 " + \
              f"-map-by ppr:1:core:PE=1 " + \
              f"python3 generate_h5_materials.py @@prop lamb mu @@tag linear_gradient @@pfx example @@dir z @@xlim -100.0 100.0 @@ylim -100.0 100.0 @zlim -100.0 100.0
            @@step 20 20 20"
        
        subprocess.run(cmd, shell=True, stdout=f)
        f.close()

        ### We increment the iteration counter
        iter += 10


if __name__ == "__main__":
    ### Call the main function
    main()